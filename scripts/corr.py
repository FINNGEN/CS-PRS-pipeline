import os
# must be set before numpy/scipy/statsmodels import: otherwise their BLAS backend (OpenBLAS/MKL)
# may spawn its own internal worker threads in this (parent) process. multiprocessing.Pool below
# uses the 'fork' start method (see its comment) to share _SHARED_DATA via copy-on-write -- but
# fork() only duplicates the calling thread, so any BLAS-internal mutex held by one of those other
# threads at fork time is inherited already-locked in every child, forever. The first BLAS call in
# a forked worker (sm.GLM.fit()'s IRLS solve) then deadlocks permanently with no error, which is
# indistinguishable from the worker just being slow -- this is what was actually causing
# correlate_pheno shards to stall at a fixed percentage rather than progress. Single-threaded BLAS
# per worker also avoids oversubscribing CPUs we're already parallelizing across at the process level.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import statsmodels.api as sm
from scipy import stats
import pandas as pd
import pyarrow.parquet as pq
import multiprocessing
import numpy as np
import time,argparse
from utils import return_header,progressBar,file_exists,make_sure_path_exists,pretty_string

cov = ['SEX_IMPUTED','AGE_AT_DEATH_OR_END_OF_FOLLOWUP', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']


def _fit_from_data(data,pheno,cov,score_col='SCORE1_AVG'):
    """
    Core GLM-fitting logic, given an already-loaded & merged pheno+scores DataFrame (may contain
    other phenotype columns too -- this subsets to just the ones this fit needs before dropna, so
    unrelated phenotypes' missingness doesn't affect this fit).
    """
    # FINNGENID isn't selected here: it was only needed earlier to join pheno data with scores,
    # nothing below uses it, and not touching it avoids the (small, but free to avoid) per-task
    # refcount overhead of a large object-dtype string column shared across worker processes
    sub = data[[pheno,score_col]+cov].copy()
    # phenotype columns come from the prebuilt parquet matrix (build_pheno_parquet.sh) as a plain
    # int8 with NA encoded as a -1 sentinel, not pandas' nullable Int8 -- real values are 0/1
    # only, so -1 always safely means missing here. Converting it back to NaN before dropna()
    # only touches this one phenotype's small per-fit slice, not the whole ~2800-column matrix.
    sub[pheno] = sub[pheno].replace(-1, np.nan)
    sub = sub.dropna()
    # statsmodels/patsy formulas want a plain float dtype, not int8
    sub[pheno] = sub[pheno].astype('float64')
    # CREATE EFFECT COLUMN
    ctrls = (sub[pheno]==0)
    avg = sub[score_col][ctrls].mean()
    std = sub[score_col][ctrls].std()
    sub['effect'] = (sub[score_col] -avg)/std
    # a covariate that's constant within this phenotype's actual fitted sample (e.g. SEX_IMPUTED,
    # once dropna() above has already restricted a sex-specific phenotype down to one sex, since
    # the other sex is coded NA rather than 0 for such phenotypes) makes the design matrix
    # singular -- drop it rather than let the fit silently fail into the except below
    usable_cov = [col for col in cov if sub[col].nunique() > 1]
    c = '+'.join(usable_cov)

    try:
        #FIT
        model = sm.GLM.from_formula(f'{pheno} ~ effect + {c}', data=sub, family=sm.families.Binomial()).fit()
        null = sm.GLM.from_formula(f'{pheno} ~  {c}', data=sub, family=sm.families.Binomial()).fit()

        # RESULTS
        ss_full = np.sum(model.resid_response**2)
        ss_red  =  np.sum(null.resid_response**2)

        results =  [model.params.effect,model.pvalues.effect,1 - ss_full/ss_red,calculate_nested_f_statistic(null,model)]
        return [pheno] + results,model.summary()

    except:
        return [pheno] + ['NA']*4,False


def fit(pheno_file,sfile,pheno='G6_ALZHEIMER',cov=cov,score_col='SCORE1_AVG',test=True):
    """
    Single-phenotype entrypoint that reads pheno_file/sfile fresh on every call. Kept as-is for
    backward compatibility / standalone use. parallel() below does NOT use this -- re-reading the
    entire pheno_file from scratch once per phenotype is the dominant cost when sweeping
    thousands of phenotypes, so it reads once up front and reuses the data via _fit_from_data.
    """
    # READ PHENO DATA
    # some COV_PHENO exports use FINNGENID as the sample id, others (plink-style FID/IID) use IID
    id_col = 'FINNGENID' if 'FINNGENID' in return_header(pheno_file) else 'IID'
    cols = [id_col,pheno] + cov
    if test:
        data =  pd.read_csv(pheno_file,sep='\t',nrows = 10000,usecols = cols)
    else:
        data = pd.read_csv(pheno_file,sep='\t',usecols = cols)
    if id_col != 'FINNGENID':
        data = data.rename(columns={id_col: 'FINNGENID'})
    #SCORES DATA
    score_data = pd.read_csv(sfile,sep='\t',usecols =['IID',score_col]).rename(columns={"IID": "FINNGENID"})
    # merge data and remove rows with NA
    data = pd.merge(data,score_data,on='FINNGENID')
    return _fit_from_data(data,pheno,cov,score_col)


# Set by parallel() before it creates its multiprocessing.Pool, and read by worker processes via
# fork's copy-on-write (Linux default start method) -- zero-copy, no per-task (re-)pickling of the
# data. Relies on the fork start method: won't auto-populate under 'spawn' (e.g. macOS default).
_SHARED_DATA = None

def wrapper_fit(args):
    pheno,cov,score_col = args
    return _fit_from_data(_SHARED_DATA,pheno,cov,score_col)



def calculate_nested_f_statistic(small_model, big_model):
    """Given two fitted GLMs, the larger of which contains the parameter space of the smaller, return the F Stat and P value corresponding to the larger model adding explanatory power"""
    addtl_params = big_model.df_model - small_model.df_model
    f_stat = (small_model.deviance - big_model.deviance) / (addtl_params * big_model.scale)
    df_numerator = addtl_params
    # use fitted values to obtain n_obs from model object:
    df_denom = (big_model.fittedvalues.shape[0] - big_model.df_model)
    p_value = stats.f.sf(f_stat, df_numerator, df_denom)
    return  p_value


def r2(model,y):
    sst_val = sum(map(lambda x: np.power(x,2),y-np.mean(y)))
    sse_val = sum(map(lambda x: np.power(x,2),model.resid_response))
    return 1.0 - sse_val/sst_val


def parallel(phenos,pheno_file,sfile,cov,score_pheno,out_path,score_col ='SCORE1_AVG',test=False,processes =1):
    """
    Parallelizes fit() across phenotypes. pheno_file is expected to be a prebuilt Parquet matrix
    (see build_pheno_parquet.sh), not the raw gz/tsv -- phenotype columns already downcast to a
    plain int8 with NA encoded as a -1 sentinel, and already restricted to a binary phenotype
    list (this fast path is only correct for binary 0/1/NA phenotypes, same restriction as the
    Binomial GLM fit below).
    """
    global _SHARED_DATA
    available_cols = pq.ParquetFile(pheno_file).schema.names
    id_col = 'FINNGENID' if 'FINNGENID' in available_cols else 'IID'
    shared_phenos = sorted(set(available_cols) & set(phenos))
    print(f"{len(shared_phenos)} phenotypes shared between list and pheno file")

    # read pheno_file + scores ONCE here instead of once per phenotype inside fit() -- reading
    # only [id_col]+shared_phenos+cov columns costs proportional to just those columns' bytes,
    # since parquet is columnar (unlike a gz/tsv, which must decompress/tokenize every column of
    # every row regardless of which ones are actually requested)
    print(f"reading {pheno_file} ({len(shared_phenos)} phenotype columns + {len(cov)} covariates)...")
    data = pd.read_parquet(pheno_file, columns=[id_col]+shared_phenos+cov)
    if test:
        data = data.head(10000)
    if id_col != 'FINNGENID':
        data.rename(columns={id_col: 'FINNGENID'}, inplace=True)
    print(f"loaded {len(data)} rows, {data.memory_usage(deep=True).sum()/1e9:.2f} GB in memory; reading scores and merging...")
    score_data = pd.read_csv(sfile,sep='\t',usecols =['IID',score_col]).rename(columns={"IID": "FINNGENID"})
    # FINNGENID dropped once merged: only needed as the join key above, nothing downstream uses it
    _SHARED_DATA = pd.merge(data,score_data,on='FINNGENID').drop(columns='FINNGENID')
    del data, score_data
    print(f"merged, {len(_SHARED_DATA)} matched samples -- starting {processes} parallel workers")

    params = [[pheno,cov,score_col] for pheno in shared_phenos]
    # explicit 'fork' context, not the bare default: this Python's default start method is
    # 'forkserver' (verified directly, not assumed), under which workers are forked from a
    # separate server process rather than from this live process -- _SHARED_DATA being set here
    # wouldn't actually reach them, and worse, 'forkserver'/'spawn' would each pickle a full
    # independent copy of the (large) shared data per worker instead of sharing pages for free.
    pool = multiprocessing.get_context('fork').Pool(processes=processes)
    results = pool.map_async(wrapper_fit,params,chunksize=1)
    while not results.ready():
        progressBar(len(params) - results._number_left,len(params))
        time.sleep(1)

    results = results.get()
    results,summaries  = zip(*results)
    pool.close()
    progressBar(1,1)
    print("\ndone.")
    # log results
    study = os.path.basename(sfile).split('.sscore')[0]
    out_root =os.path.join(out_path, study) + '_corr'
    out_file = out_root + '.txt'
    with open(out_file,'wt') as o:
        print(f'saving results to f{out_file}')
        o.write("\t".join(["PHENO",'beta','pval','p_R2','pval_F','study']) + '\n')
        sig_res = [elem for elem in results if 'NA' not in elem]
        sorted_sig = sorted(sig_res, key=lambda x: x[2])
        na_res = [elem for elem in results if 'NA'  in elem]
        for res in sorted_sig + na_res:
            values = ['{:.2e}'.format(elem) if isinstance(elem,float) else elem for elem in res]
            values.append(score_pheno +'_'+ study)
            o.write('\t'.join(values) + '\n')

    log_file = out_root + '.log'
    with open(log_file,'wt') as o:
        print(f"saving logs to {log_file}")
        l = len([str(s) for s in summaries if s][0].split('\n')[0])
        o.write('='*l + '\n')
        o.write(pretty_string(study,l) + '\n\n')
        for summary in summaries:
            if summary:
                o.write(str(summary) + '\n')



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description ="Correlation between PRS and other phenos.")

    parser.add_argument('--pheno-file',type = file_exists,help ='File that contains the age of onset',required = True)
    parser.add_argument('--scores',type = file_exists,help = 'score file',required = True)
    parser.add_argument('--pheno',type = str,help ='pheno of scores',required = True)

    parser.add_argument('--out',type = str, help ='output_path',default = '.')
    parser.add_argument('--cov',type = str,help='comma separated list of covariates',default =  ','.join(cov))
    parser.add_argument('--pheno-list',type = file_exists,help='File with list of phenotypes to run correlation for',required = True)
    parser.add_argument('--scol',type=str,help = 'Column name in scores file',default = 'SCORE1_AVG')
    parser.add_argument('--cpus',type=int,help = 'Number of parallel jobs',default = multiprocessing.cpu_count() -1)
    parser.add_argument('--test',action='store_true',default = False)

    args = parser.parse_args()
    phenos = np.loadtxt(args.pheno_list,dtype = str)
    covariates = args.cov.split(',')
    make_sure_path_exists(args.out)
    print(f"{len(phenos)} phenos provided")
    parallel(phenos,args.pheno_file,args.scores,covariates,args.pheno,args.out,args.scol,processes= min(len(phenos),args.cpus),test = args.test)
