#!/usr/bin/env python3
"""
Correlation between PRS and phenotypes: one GLM fit per phenotype, parallelized as independent
OS processes (via GNU parallel), not Python multiprocessing -- avoids large-shared-object memory
duplication across forked workers. Reads a prebuilt phenotype parquet (build_pheno_parquet.sh),
not the raw gz/tsv.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import argparse
import shlex
import shutil
import subprocess
import sys
import time
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
import pyarrow
import pyarrow.parquet as pq
pyarrow.set_cpu_count(1)
pyarrow.set_io_thread_count(1)

from utils import file_exists, make_sure_path_exists, pretty_string, progressBar

cov = ['SEX_IMPUTED', 'AGE_AT_DEATH_OR_END_OF_FOLLOWUP', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']


def calculate_nested_f_statistic(small_model, big_model):
    addtl_params = big_model.df_model - small_model.df_model
    f_stat = (small_model.deviance - big_model.deviance) / (addtl_params * big_model.scale)
    df_denom = big_model.fittedvalues.shape[0] - big_model.df_model
    return stats.f.sf(f_stat, addtl_params, df_denom)


def _fit_from_data(data, pheno, cov, score_col='SCORE1_AVG'):
    sub = data[[pheno, score_col] + cov].copy()
    # phenotype columns are int8 with NA as a -1 sentinel (see build_base) -- real values are 0/1
    sub[pheno] = sub[pheno].replace(-1, np.nan)
    sub = sub.dropna()
    sub[pheno] = sub[pheno].astype('float64')

    ctrls = (sub[pheno] == 0)
    avg = sub[score_col][ctrls].mean()
    std = sub[score_col][ctrls].std()
    sub['effect'] = (sub[score_col] - avg) / std
    # drop covariates constant within this fit's sample (e.g. SEX_IMPUTED for a sex-specific
    # pheno) -- would make the design matrix singular
    usable_cov = [c for c in cov if sub[c].nunique() > 1]
    c = '+'.join(usable_cov)

    try:
        model = sm.GLM.from_formula(f'{pheno} ~ effect + {c}', data=sub, family=sm.families.Binomial()).fit()
        null = sm.GLM.from_formula(f'{pheno} ~ {c}', data=sub, family=sm.families.Binomial()).fit()
        ss_full = np.sum(model.resid_response ** 2)
        ss_red = np.sum(null.resid_response ** 2)
        results = [model.params.effect, model.pvalues.effect, 1 - ss_full / ss_red, calculate_nested_f_statistic(null, model)]
        return [pheno] + results, model.summary()
    except Exception:
        return [pheno] + ['NA'] * 4, False


def build_base(pheno_file, sfile, cov, score_col, out_file):
    """Small numeric-only base (covariates + score + _orig_idx, no phenotype columns, no id
    strings) that every per-phenotype fit loads instead of the huge parquet. Rebuilt every run
    (a few seconds) rather than cached, so it can never silently go stale."""
    available_cols = pq.ParquetFile(pheno_file).schema.names
    id_col = 'FINNGENID' if 'FINNGENID' in available_cols else 'IID'

    base = pd.read_parquet(pheno_file, columns=[id_col] + cov)
    if id_col != 'FINNGENID':
        base.rename(columns={id_col: 'FINNGENID'}, inplace=True)
    base['_orig_idx'] = np.arange(len(base))

    score_data = pd.read_csv(sfile, sep='\t', usecols=['IID', score_col]).rename(columns={"IID": "FINNGENID"})
    merged = pd.merge(base, score_data, on='FINNGENID').drop(columns='FINNGENID')
    del base, score_data

    merged.to_parquet(out_file, index=False)
    print(f"wrote {out_file}: {len(merged)} matched samples, {len(cov)} covariates", flush=True)


def fit_one(base_file, pheno_file, pheno, cov, score_col, score_pheno, study, tmp_dir):
    base = pd.read_parquet(base_file)
    orig_idx = base['_orig_idx'].to_numpy()
    base = base.drop(columns='_orig_idx')

    col = pq.ParquetFile(pheno_file).read(columns=[pheno]).column(0).to_numpy(zero_copy_only=False)
    data = base.assign(**{pheno: col[orig_idx]})
    result, summary = _fit_from_data(data, pheno, cov, score_col)

    with open(os.path.join(tmp_dir, f"{pheno}.tsv"), 'wt') as o:
        values = ['{:.2e}'.format(v) if isinstance(v, float) else v for v in result]
        values.append(f"{score_pheno}_{study}")
        o.write('\t'.join(values) + '\n')
    with open(os.path.join(tmp_dir, f"{pheno}.log"), 'wt') as o:
        if summary:
            o.write(str(summary) + '\n')


def run_parallel(phenos, cpus, common_argv, tmp_dir, memfree='2G', retries=5):
    """GNU parallel instead of xargs: xargs aborts its *entire remaining queue* the moment any
    child is killed by a signal (verified directly -- not just the killed item, everything not
    yet started), which made an OOM kill far more disruptive than it needed to be. parallel
    continues past a killed job by default, and --memfree proactively avoids most OOM kills in
    the first place (won't start a new job below `memfree` free memory; kills+requeues the
    youngest job, not an OS-chosen one, if free memory drops below half that) with --retries as
    a backstop for whatever --memfree doesn't catch in time."""
    par = shutil.which('parallel')
    if not par:
        raise RuntimeError("GNU parallel not found on PATH")
    worker_argv = [sys.executable, os.path.abspath(__file__)] + common_argv + ['--_worker', '{}']
    cmd = [par, '-j', str(cpus), '--memfree', memfree, '--retries', str(retries)] + worker_argv
    print(f"launching {len(phenos)} phenotype fits, {cpus} at a time (memfree={memfree}, retries={retries}):\n  {shlex.join(cmd)}", flush=True)

    total = len(phenos)
    def n_done():
        return sum(1 for p in phenos if os.path.exists(os.path.join(tmp_dir, f"{p}.tsv")))

    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, text=True)
    proc.stdin.write('\n'.join(phenos))
    proc.stdin.close()
    while proc.poll() is None:
        progressBar(n_done(), total)
        time.sleep(1)
    progressBar(n_done(), total)
    print("\ndone.")

    if proc.returncode != 0:
        print(f"WARNING: parallel exited {proc.returncode} -- some phenotypes may have failed (check {tmp_dir})", file=sys.stderr)


def gather_results(tmp_dir, out_dir, study, score_pheno):
    rows, na_rows, summaries = [], [], []
    for fname in sorted(os.listdir(tmp_dir)):
        if fname.endswith('.tsv'):
            with open(os.path.join(tmp_dir, fname)) as f:
                fields = f.readline().rstrip('\n').split('\t')
            (na_rows if 'NA' in fields else rows).append(fields)
        elif fname.endswith('.log'):
            with open(os.path.join(tmp_dir, fname)) as f:
                text = f.read()
            if text.strip():
                summaries.append(text)

    rows.sort(key=lambda f: (float(f[2]), float(f[1])))

    out_root = os.path.join(out_dir, study) + '_corr'
    with open(out_root + '.txt', 'wt') as o:
        print(f'saving results to {out_root}.txt', flush=True)
        o.write('\t'.join(['PHENO', 'beta', 'pval', 'p_R2', 'pval_F', 'study']) + '\n')
        for fields in rows + na_rows:
            o.write('\t'.join(fields) + '\n')

    with open(out_root + '.log', 'wt') as o:
        print(f'saving logs to {out_root}.log', flush=True)
        l = len(summaries[0].split('\n')[0]) if summaries else 80
        o.write('=' * l + '\n')
        o.write(pretty_string(study, l) + '\n\n')
        for summary in summaries:
            o.write(summary + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Correlation between PRS and other phenos.")
    parser.add_argument('--pheno-file', type=file_exists, required=True, help='Prebuilt phenotype parquet')
    parser.add_argument('--scores', type=file_exists, required=True)
    parser.add_argument('--pheno', type=str, required=True, help="Study's own clinical phenotype label")
    parser.add_argument('--pheno-list', type=file_exists, default=None, help='Required unless --_worker is set')
    parser.add_argument('--cov', type=str, default=','.join(cov))
    parser.add_argument('--scol', type=str, default='SCORE1_AVG')
    parser.add_argument('--cpus', type=int, default=max(1, os.cpu_count() - 1))
    parser.add_argument('--memfree', type=str, default='2G', help="GNU parallel --memfree: min free memory before starting another job")
    parser.add_argument('--retries', type=int, default=5, help="GNU parallel --retries: retry a killed/failed phenotype this many times")
    parser.add_argument('--out', type=str, default='.')
    parser.add_argument('--_worker', type=str, default=None, help=argparse.SUPPRESS)
    args = parser.parse_args()

    covariates = args.cov.split(',')
    study = os.path.basename(args.scores).split('.sscore')[0]
    make_sure_path_exists(args.out)
    base_path = os.path.join(args.out, 'base.parquet')
    tmp_dir = os.path.join(args.out, '_tmp')

    if args._worker:
        fit_one(base_path, args.pheno_file, args._worker, covariates, args.scol, args.pheno, study, tmp_dir)
        sys.exit(0)

    if not args.pheno_list:
        parser.error("--pheno-list is required")

    build_base(args.pheno_file, args.scores, covariates, args.scol, base_path)

    available_cols = pq.ParquetFile(args.pheno_file).schema.names
    with open(args.pheno_list) as f:
        phenos = [l.strip() for l in f if l.strip()]
    shared_phenos = sorted(set(available_cols) & set(phenos))
    print(f"{len(shared_phenos)} phenotypes shared between list and pheno file", flush=True)

    # fresh tmp_dir every run: a shared/reused --out would otherwise let stale .tsv files from a
    # previous (possibly killed, possibly different-phenotype-list) run inflate the progress bar
    # and contaminate the final gathered output
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    make_sure_path_exists(tmp_dir)

    common_argv = ['--pheno-file', args.pheno_file, '--scores', args.scores, '--pheno', args.pheno,
                   '--cov', args.cov, '--scol', args.scol, '--out', args.out]
    run_parallel(shared_phenos, args.cpus, common_argv, tmp_dir, memfree=args.memfree, retries=args.retries)

    # parallel's own --memfree/--retries handle OOM recovery; this is only a backstop for a
    # phenotype that still failed for some other reason after all of parallel's own retries
    missing = [p for p in shared_phenos if not os.path.exists(os.path.join(tmp_dir, f"{p}.tsv"))]
    if missing:
        print(f"ERROR: {len(missing)}/{len(shared_phenos)} phenotypes produced no result even after "
              f"parallel's own retries: {missing[:10]}{' ...' if len(missing) > 10 else ''}", file=sys.stderr)
        sys.exit(1)

    gather_results(tmp_dir, args.out, study, args.pheno)
