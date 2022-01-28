# CS-PRS-pipeline

Pipeline to calculates PRS based on a list of sumstats.
Weights are calculated with PRScs: https://github.com/getian107/PRScs

The sumstats used for R8 are the following:

| filename | pheno | publication |
|---|---|---|
| AD_sumstats_Jansenetal.txt.gz | G6_ALZHEIMER | https://www.nature.com/articles/s41588-018-0311-9 |
| adhd_eur_jun2017.gz | F5_ADHD | https://www.biorxiv.org/content/early/2017/06/03/145581 |
| alsMetaSummaryStats_march21st2018.tab.gz | G6_ALS | https://www.cell.com/neuron/references/S0896-6273(18)30148-X |
| Bipolar_vs_control_PGC_Cell_2018_formatted.txt.gz | F5_BIPO | PGC;Cell;2018;https://doi.org/10.1016/j.cell.2018.05.046 |
| cad.add.160614.website.txt.gz | I9_CHD | https://www.nature.com/articles/ng.3396 |
| ckqny.scz2snpres.gz | F5_SCHZPHR | https://www.nature.com/articles/nature13595 |
| daner_PGC_SCZ52_0513a.resultfiles_PGC_SCZ52_0513.sh2_nofin.gz | F5_SCHZPHR | PGC_schzizophrenia_GWAS_Ripke_et_al._excluding_Finnish |
| dpw_excludingFinnishStudies_4finngen.txt.gz | alcohol_consumption | https://www.nature.com/articles/s41398-019-0676-2 |
| EAGLE_AD_GWAS_results_2015.txt.gz | L12_ATOPIC | https://www.nature.com/articles/ng.3424 |
| Educational_Attainment_excl23andme_Lee_2018_NatGen_formatted.txt.gz | educational_attainment | Lee;Nat_Gen;2018;https://doi.org/10.1038/s41588-018-0147-3 |
| EUR.CD.gwas_info03_filtered.assoc.gz | K11_CROHN | https://www.nature.com/articles/ng.3359 |
| EUR.IBD.gwas_info03_filtered.assoc.gz | K11_IBD | https://www.nature.com/articles/ng.3359 |
| EUR.UC.gwas_info03_filtered.assoc.gz | K11_ULCER | https://www.nature.com/articles/ng.3359 |
| focal_epilepsy_METAL_4finngen.txt.gz | FE | https://www.nature.com/articles/s41467-018-07524-z |
| gabriel_asthma_meta-analysis_36studies_format_repository_NEJM.txt.gz | J10_ASTHMA | https://www.nejm.org/doi/full/10.1056/nejmoa0906312 |
| generalised_epilepsy_METAL_4finngen.txt.gz | GE | https://www.nature.com/articles/s41467-018-07524-z |
| GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz | Height | https://www.nature.com/articles/ng.3097 |
| GWAS_CP_all.txt.gz | Cognitive_performance_meta | https://www.nature.com/articles/s41588-018-0147-3 |
| Hb_gwas_summary_fromNealeLab.tsv.gz | Haemoglobin_concentration | UKBB |
| HbA1c_METAL_European.txt.gz | HbA1c | http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1002383 |
| HDL_Meta_ENGAGE_1000G_non_FINRISK_QCd.txt.gz | HDL-C | https://www.nature.com/articles/ng.3300 |
| IGAP_stage_1.txt.gz | G6_ALZHEIMER | https://www.nature.com/articles/ng.2802 |
| iPSYCH-PGC_ASD_Nov2017.gz | F5_AUTISM | https://www.biorxiv.org/content/early/2017/11/27/224774 |
| kunkle_etal_stage1.txt.gz | G6_ALZHEIMER | https://www.niagads.org/datasets/ng00075 |
| LDL_Meta_ENGAGE_1000G_non_FINRISK_QCd.txt.gz | LDL-C | https://www.nature.com/articles/ng.3300 |
| Mahajan.NatGenet2018b.T2D.European.txt.gz | T2D | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6287706/ |
| MDD2018_ex23andMe.19fields.gz | F5_DEPRESSIO | https://www.nature.com/articles/s41588-018-0090-3 |
| meta_v3_onco_euro_overall_ChrAll_1_release.txt.gz | C3_PROSTATE_EXALLC | https://www.nature.com/articles/s41588-018-0142-8 |
| METAANALYSIS_DIAGRAM_SE1.txt.gz | T2D | http://diabetes.diabetesjournals.org/content/66/11/2888 |
| metastroke.all.chr.bp.gz | I9_STR | https://www.thelancet.com/journals/laneur/article/PIIS1474-4422(12)70234-X/abstract |
| MTAG_CP.to10K.txt.gz | Cognitive_Performance | https://www.nature.com/articles/s41588-018-0147-3 |
| MTAG_EA.to10K.txt.gz | Educational_attainment | https://www.nature.com/articles/s41588-018-0147-3 |
| pgc.bip.full.2012-04.txt.gz | F5_BIPO | https://www.nature.com/articles/ng.943 |
| pgc.ed.freeze1.summarystatistics.July2017.txt.gz | F5_ANOREX | Missing |
| RA_GWASmeta_European_v2.txt.gz | M13_RHEUMA | https://www.nature.com/articles/nature12873 |
| SavageJansen_2018_intelligence_metaanalysis_formatted.txt.gz | intelligence | Savage;Nat_Gen-2018;https://doi.org/10.1038/s41588-018-0152-6 |
| Saxena_fullUKBB_Longsleep_summary_stats_formatted.txt.gz | long_sleep | Dashti;Natcomm-2019;https://doi.org/10.1038/s41467-019-08917-4 |
| Shrine_30804560_FEV1_meta-analysis.txt.gz | FEV1 | https://www.nature.com/articles/s41588-018-0321-7 |
| Shrine_30804560_FEV1_to_FVC_RATIO_meta-analysis.txt.gz | FEV1/FVC | https://www.nature.com/articles/s41588-018-0321-7 |
| Shrine_30804560_FVC_meta-analysis.txt.gz | FVC | https://www.nature.com/articles/s41588-018-0321-7 |
| Shrine_30804560_PEF_meta-analysis.txt.gz | Peak_expiratory_flow | https://www.nature.com/articles/s41588-018-0321-7 |
| SNP_gwas_mc_merge_nogc.tbl.uniq.gz | BMI | https://www.nature.com/articles/nature14177 |
| SORTED_PTSD_EA9_ALL_study_specific_PCs1.txt.gz | F5_PTSD | https://www.nature.com/articles/mp201777 |
| sumstats_neuroticism_ctg_format_formatted.txt.gz | neuroticism | Nagel;NatGen2018-https://doi.org/10.1038/s41588-018-0151-7 |
| TC_Meta_ENGAGE_1000G_non_FINRISK_QCd.txt.gz | Total-cholesterol | https://www.nature.com/articles/ng.3300 |
| TG_Meta_ENGAGE_1000G_non_FINRISK_QCd.txt.gz | Triglycerides | https://www.nature.com/articles/ng.3300 |
| UKBB_gwas_Neale_Chronotype_formatted.txt.gz | chronotype | UKBB |
| UKBB_gwas_Neale_pulserate_formatted.txt.gz | heart_rate | UKBB |
| UKBB_gwas_Neale_SleepDuration_formatted.txt.gz | sleepduration | UKBB |
| UKBB_gwas_Neale_SleeplessnessInsomnia_formatted.txt.gz | insomnia | UKBB |
| UKBB_HYPO.txt.gz | HYPOTHYROIDISM | UKBB |
| UKB-ICBPmeta750k_DBPsummaryResults.txt.gz | I9_HYPTENSESS | UKBB_ICBP |
| UKB-ICBPmeta750k_PPsummaryResults.edited.txt.gz | I9_HYPTENSESS | UKBB_ICBP |
| UKB-ICBPmeta750k_SBPsummaryResults.txt.gz | I9_HYPTENSESS | UKBB_ICBP |
| GSCAN_AgeofInitiation.txt.gz | smoking_age_of_initiation | https://www.nature.com/articles/s41588-018-0307-5 |
| GSCAN_CigarettesPerDay.txt.gz | smoking_cigarettes_per_day | https://www.nature.com/articles/s41588-018-0307-5 |
| GSCAN_DrinksPerWeek.txt.gz | drinks_per_week | https://www.nature.com/articles/s41588-018-0307-5 |
| GSCAN_SmokingInitiation.txt.gz | smoking_initiation | https://www.nature.com/articles/s41588-018-0307-5 |
| PGC_UKB_depression_genome-wide.txt.gz | depression | https://www.nature.com/articles/s41588-018-0090-3 |
| PGC_UKB_23andMe_depression_10000.txt.gz | depression_23andme | https://www.nature.com/articles/s41588-018-0090-3 |
| PGC3_SCZ_wave3_public_without_frequencies.v1.tsv.gz | schizophrenia | https://www.medrxiv.org/content/10.1101/2020.09.12.20192922v1 |
| PGC3_SCZ_wave3_public_without_frequencies.clumped.v1.tsv.gz | schizophrenia_clumped | https://www.medrxiv.org/content/10.1101/2020.09.12.20192922v1 |
| daner_PGC_BIP32b_mds7a_0416a.gz | bipolar_disorder | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6956732/ |
| CIMBA_BRCA1_BCAC_TN_meta_summary_level_statistics.txt.gz | breast_cancer | https://www.nature.com/articles/s41588-020-0609-2 |
| AF_HRC_GWAS_ALLv11.txt.gz | atrial_fibrillation | https://www.nature.com/articles/s41588-018-0133-9 |
| lifegen_phase2_bothpl_alldr_2017_09_18.tsv.gz | lifespan | https://elifesciences.org/articles/39856 |
| RISK_GWAS_MA_UKB+23andMe+replication.txt.gz | risk_tolerance | https://www.nature.com/articles/s41588-018-0309-3 |
| chronic_pain-bgen.stats.gz | multisite_chronic_pain | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6592570/ |
| Saxena_fullUKBB_Insomnia_summary_stats.txt.gz | insomnia | https://www.nature.com/articles/s41588-019-0361-7 |
| continuous-LDLC-both_sexes-medadj_irnt.tsv.gz | ldl_adj_irnt | https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=511623409 |
| biomarkers-30750-both_sexes-irnt.tsv.gz | hba1c_inrt | https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=511623409 |
| biomarkers-30870-both_sexes-irnt.tsv.gz | triglycerides_irnt | https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=511623409 |
| MEGASTROKE.1.AS.TRANS.out.gz | any_stroke | https://www.nature.com/articles/s41588-018-0058-3 |
| UKBB.asthma-2.assoc.gz | asthma | https://www.nature.com/articles/s41467-020-15649-3 |
| ICP_FG_DC_EG_DBDS_202012.gz | ICP_WIDE_EXMORE | FG |
| Meta-analysis_Wood_et_al+UKBiobank_2018_height.txt.gz | height | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6488973/ |
| Meta-analysis_Locke_et_al+UKBiobank_2018_bmi.txt.gz | bmi | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6488973/ |
| GSEM.GWAS.EXTERNALIZING.SHARE.v20191014.txt.gz | externalizing | https://www.biorxiv.org/content/10.1101/2020.10.16.342501v1 |
| MTAG_glaucoma_four_traits_summary_statistics.txt.gz | glaucoma_ | https://www.nature.com/articles/s41588-019-0556-y |
| phecode-696.4-both_sexes.tsv.gz | Psoriasis | UKBB |
| phecode-250.1-both_sexes.tsv.gz | T1D | UKBB |
| BCX2_EOS_Trans_GWAMA.out.gz | EOS_Trans | https://pubmed.ncbi.nlm.nih.gov/32888493/ |
| biomarkers-30600-both_sexes-irnt.tsv.gz | Albumin | https://pan.ukbb.broadinstitute.org/ |
| biomarkers-30610-both_sexes-irnt.tsv.gz | Alkaline_phosphatase | https://pan.ukbb.broadinstitute.org/ |
| biomarkers-30620-both_sexes-irnt.tsv.gz | Alanine_aminotransferase | https://pan.ukbb.broadinstitute.org/ |
| biomarkers-30630-both_sexes-irnt.tsv.gz | Apolipoprotein_A | https://pan.ukbb.broadinstitute.org/ |
| biomarkers-30640-both_sexes-irnt.tsv.gz | Apolipoprotein_B | https://pan.ukbb.broadinstitute.org/ |
| biomarkers-30650-both_sexes-irnt.tsv.gz | Aspartate_aminotransferase | https://pan.ukbb.broadinstitute.org/ |
| biomarkers-30660-both_sexes-irnt.tsv.gz | Direct_bilirubin | https://pan.ukbb.broadinstitute.org/ |
| biomarkers-30670-both_sexes-irnt.tsv.gz | Urea | https://pan.ukbb.broadinstitute.org/ |
| biomarkers-30680-both_sexes-irnt.tsv.gz | Calcium | https://pan.ukbb.broadinstitute.org/ |
| biomarkers-30690-both_sexes-irnt.tsv.gz | Cholesterol | https://pan.ukbb.broadinstitute.org/ |
| biomarkers-30700-both_sexes-irnt.tsv.gz | Creatinine | https://pan.ukbb.broadinstitute.org/ |
| biomarkers-30710-both_sexes-irnt.tsv.gz | C-reactive_protein | https://pan.ukbb.broadinstitute.org/ |
| biomarkers-30720-both_sexes-irnt.tsv.gz | Cystatin_C | https://pan.ukbb.broadinstitute.org/ |
| biomarkers-30730-both_sexes-irnt.tsv.gz | Gamma_glutamyltransferase | https://pan.ukbb.broadinstitute.org/ |
| biomarkers-30740-both_sexes-irnt.tsv.gz | Glucose | https://pan.ukbb.broadinstitute.org/ |
| biomarkers-30760-both_sexes-irnt.tsv.gz | HDL_cholesterol | https://pan.ukbb.broadinstitute.org/ |
| biomarkers-30770-both_sexes-irnt.tsv.gz | IGF-1 | https://pan.ukbb.broadinstitute.org/ |
| biomarkers-30780-both_sexes-irnt.tsv.gz | LDL_direct | https://pan.ukbb.broadinstitute.org/ |
| biomarkers-30790-both_sexes-irnt.tsv.gz | Lipoprotein_A | https://pan.ukbb.broadinstitute.org/ |
| biomarkers-30800-both_sexes-irnt.tsv.gz | Oestradiol | https://pan.ukbb.broadinstitute.org/ |
| biomarkers-30810-both_sexes-irnt.tsv.gz | Phosphate | https://pan.ukbb.broadinstitute.org/ |
| biomarkers-30820-both_sexes-irnt.tsv.gz | Rheumatoid_factor | https://pan.ukbb.broadinstitute.org/ |
| biomarkers-30830-both_sexes-irnt.tsv.gz | SHBG | https://pan.ukbb.broadinstitute.org/ |
| biomarkers-30840-both_sexes-irnt.tsv.gz | Total_bilirubin | https://pan.ukbb.broadinstitute.org/ |
| biomarkers-30850-both_sexes-irnt.tsv.gz | Testosterone | https://pan.ukbb.broadinstitute.org/ |
| biomarkers-30860-both_sexes-irnt.tsv.gz | Total_protein | https://pan.ukbb.broadinstitute.org/ |
| biomarkers-30880-both_sexes-irnt.tsv.gz | Urate | https://pan.ukbb.broadinstitute.org/ |
| biomarkers-30890-both_sexes-irnt.tsv.gz | Vitamin_D | https://pan.ukbb.broadinstitute.org/ |
| ukb.allasthma.upload.final.assoc_v2.txt.gz | asthma | https://erj.ersjournals.com/content/54/6/1901507-0.long |
| michailidoubreastcancerall.txt.gz | Breast_cancer | https://www.nature.com/articles/nature24284 |
| 20171017_MW_eGFR_overall_EA_nstud42.dbgap_v2.txt.gz | estimated_glomerular_filtration_rate | https://www.nature.com/articles/s41588-019-0407-x |
| colorectal_cancer.txt.gz | Colorectal_cancer | https://www.nature.com/articles/s41467-019-09775-w |
| MTAG_glaucoma_four_traits_summary_statistics.txt.gz | Glaucoma_multitrait | https://www.nature.com/articles/s41588-019-0556-y |
| urate_chr1_22_LQ_IQ06_mac10_all_741_rsid_v2.txt.gz | Urate | https://www.nature.com/articles/s41588-019-0504-x |
| Tachmazidou_30664745_HIPOA.txt.gz | hip_osteoarthritis | https://www.nature.com/articles/s41588-018-0327-1 |
| Tachmazidou_30664745_KNEEOA.txt.gz | knee_osteoarthritis | https://www.nature.com/articles/s41588-018-0327-1 |
| BMD_v3_SumStats.txt.gz | Bone_mineral_density | https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0200785 |
| anno.CLEANED.MVP.EUR.VTE.results.EQC.tbl.gz | Venous_thromboembolism | https://www.nature.com/articles/s41588-019-0519-3 |
| GCST90011770_buildGRCh37.tsv.gz | Glaucoma | https://www.nature.com/articles/s41467-020-20851-4 |
| Hysi_et_al_Refractive_error_NatGenet_2020_beta_se.txt.gz | Refractive_error_and_myopia | https://www.nature.com/articles/s41588-020-0599-0 |
| ALS_sumstats_EUR_only.txt.gz | ALS | https://doi.org/10.1101/2021.03.12.21253159 |
|discovery_metav3.0.meta.gz | Multiple_Sclerosis | https://pubmed.ncbi.nlm.nih.gov/31604244/ |


## Rsid map

This step generates a mapping to/from rsid/chrompos based on data available at ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz.

`rsid_map.py` produces:

- finngen.rsid.map.tsv (rsid--> chrompos)
```
rs10	7_92754574
rs1000000	12_126406434
rs1000000219	13_95689463
```

- finngen.variants.tsv (chrompos --> ref/alt)
```
10_100000235	C	T
10_100000979	T	C
10_100001839	CAA	C
```

The first file is used throught the computation to go to/from rsid notation. The second is used to filter out variants that do not have the right alleles.

Also, if a rsid list is provided (e.g. hm3 rsids), it returns the subset of variants in the original bim file that match to those rsids:
- hm3.snplist
```
chr10_100000235_C_T
chr10_100002628_A_C
chr10_100004827_A_C
```

## Munging
PRScs automatically does some allele matching:
- the reference genome (1kg) only contains non ambigous variants :

```
cat snpinfo_1kg_hm3 | sed -E 1d |  cut -f 4,5 | awk '{print $1$2}' | sort | uniq
AC
AG
CA
CT
GA
GT
TC
TG
```

Also, in the parsing phase it checks for the ref/alt order and fixes the beta accordingly.
```
vld_snp = set(zip(vld_dict['SNP'], vld_dict['A1'], vld_dict['A2']))
ref_snp = set(zip(ref_dict['SNP'], ref_dict['A1'], ref_dict['A2'])) | set(zip(ref_dict['SNP'], ref_dict['A2'], ref_dict['A1'])) | \
              set(zip(ref_dict['SNP'], [mapping[aa] for aa in ref_dict['A1']], [mapping[aa] for aa in ref_dict['A2']])) | \
              set(zip(ref_dict['SNP'], [mapping[aa] for aa in ref_dict['A2']], [mapping[aa] for aa in ref_dict['A1']]))
sst_snp = set(zip(sst_dict['SNP'], sst_dict['A1'], sst_dict['A2'])) | set(zip(sst_dict['SNP'], sst_dict['A2'], sst_dict['A1'])) | \
              set(zip(sst_dict['SNP'], [mapping[aa] for aa in sst_dict['A1'] if aa in ATGC], [mapping[aa] for aa in sst_dict['A2'] if aa in ATGC])) | \
              set(zip(sst_dict['SNP'], [mapping[aa] for aa in sst_dict['A2'] if aa in ATGC], [mapping[aa] for aa in sst_dict['A1'] if aa in ATGC]))

comm_snp = ref_snp & vld_snp & sst_snp
with open(sst_file) as ff:
     if (snp, a1, a2) in comm_snp:
     	...   
	beta_std = sp.sign(beta)*abs(norm.ppf(p/2.0))/n_sqrt
     elif (snp, a2, a1) in comm_snp:
     	beta_std = -1*sp.sign(beta)*abs(norm.ppf(p/2.0))/n_sqrt
	 ...

```
The final weights are printed based on the a1/a2 order of the reference panel (i.e. the EUR 1kg panel in this case).

PRScs does check for strand flip.

Our solution is therefore the following.

We build a rsid to chrom pos mapping from `ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz` . This allows to move back and forth from/to rsid/chrompos notations and therefore to merge summary stats with different formats.

Then:
- in the munging phase we split the summary stats entries based on whether variants are identified by rsid or by some of chrom/pos notation
- the rsid file is filtered for rsids present in finngen. chrom and pos information are updated to finngen data
- the chrompos file is updated to have chrom and pos added (if provided, else extracted from variant id) and then lifted to build 38
- the two files are then merged to a FINNGEN chrom_pos_ref_alt notation, making sure that the variant exists in finngen data (checking for strand flip as well).

This produces a file with chrom and pos based on Finngen, but with A1/A2/OR/P based on the original data:

```
CHR     SNP     A1      A2      BP      OR      P
19	chr19_260912	G	A	260912	0.9957872809136792	0.050031874722
19	chr19_261033	A	G	261033	0.9957998422696626	0.0507241244322
19	chr19_266034	C	T	266034	1.0053796666398893	0.150796052266
19	chr19_267039	C	T	267039	0.9961140904000996	0.0691110781996
19	chr19_276245	T	C	276245	0.995420944482037	0.0366293445837
19	chr19_277776	A	G	277776	0.9964618752820534	0.119939685495
19	chr19_280299	C	T	280299	0.9964396546231319	0.120269108388
19	chr19_281360	T	C	281360	0.9968755417956986	0.169347605744
19	chr19_288246	C	T	288246	0.9991125513478095	0.722369442562
19	chr19_288374	C	T	288374	1.0021901113140002	0.299346915861
```


## Weights

Weights are calculated using PRScs.
In order to run PRScs we then convert the file to rsids:
```
SNP	A1	A2	OR	P
rs8100066	G	A	0.9957872809136792	0.050031874722
rs8105536	A	G	0.9957998422696626	0.0507241244322
rs2312724	C	T	1.0053796666398893	0.150796052266
rs1020382	C	T	0.9961140904000996	0.0691110781996
rs12459906	T	C	0.995420944482037	0.0366293445837
rs11084928	A	G	0.9964618752820534	0.119939685495
rs11878315	C	T	0.9964396546231319	0.120269108388
rs7815	T	C	0.9968755417956986	0.169347605744
rs10409452	C	T	0.9991125513478095	0.722369442562
rs12981067	C	T	1.0021901113140002	0.299346915861
```

This guarantees that the beta is still correct,since it's based on the original summary stats. However, this way, we can "recycle" the munged data also for other reference panels, if needed in the future.

Then PRScs is run and weights are calculated only for the subset of variants shared across reference panel, summary stats and validation bim file (finngen).

```
19	rs8100066	260912	G	A	-2.438700e-05
19	rs8105536	261033	A	G	-1.608225e-05
19	rs2312724	266034	C	T	2.586432e-04
19	rs1020382	267039	C	T	8.532887e-06
19	rs12459906	276245	T	C	-3.629306e-05
19	rs11084928	277776	A	G	-3.484712e-05
19	rs11878315	280299	C	T	-6.442467e-06
19	rs7815	281360	T	C	-1.508200e-05
19	rs10409452	288246	C	T	2.060791e-05
19	rs12981067	288374	C	T	4.191780e-05
```


The weight file is converted to chrom_pos again through the finngen rsid/chrom_pos mapping,using the a1/a2 from the weights. However, now there is a double mismatch that needs to be fixed:
1) the weights were calculated based on the a1/a2 order of the reference data set
2) the output positions are based on the reference data set.

```
19      chr19_1208073_C_T       1208072 C       T       7.354914e-06
19      chr19_1218220_T_C       1218219 T       C       5.217805e-06
19      chr19_1220005_G_A       1220004 G       A       9.294323e-05
19      chr19_1221162_T_C       1221161 T       C       1.765576e-05
19      chr19_1226005_A_C       1226004 A       C       8.979659e-05
19      chr19_1232559_C_T       1232558 C       T       4.271571e-05
19      chr19_1238900_C_T       1238899 C       T       2.828170e-05

```

In order to fix this, we replicate each entry, considering all possible permutations of the ref_alt in the variant id. This guarantees that at least one permutation is the matching Finngen variant. Also, the position is updated to the one in the id.

```
19	chr19_1208073_C_T	1208073	C	T	7.354914e-06
19	chr19_1208073_T_C	1208073	C	T	7.354914e-06
19	chr19_1208073_G_A	1208073	C	T	7.354914e-06
19	chr19_1208073_A_G	1208073	C	T	7.354914e-06
19	chr19_1218220_T_C	1218220	T	C	5.217805e-06
19	chr19_1218220_C_T	1218220	T	C	5.217805e-06
19	chr19_1218220_A_G	1218220	T	C	5.217805e-06
19	chr19_1218220_G_A	1218220	T	C	5.217805e-06
19	chr19_1220005_G_A	1220005	G	A	9.294323e-05
19	chr19_1220005_A_G	1220005	G	A	9.294323e-05
19	chr19_1220005_C_T	1220005	G	A	9.294323e-05
19	chr19_1220005_T_C	1220005	G	A	9.294323e-05
19	chr19_1221162_T_C	1221162	T	C	1.765576e-05
19	chr19_1221162_C_T	1221162	T	C	1.765576e-05
```

Now we have all elements in place:
- variants are identified with a Finngen ID
- the position is updated to finngen data
- the effect allele is still the original one
- weights are calcualted accordingly based on the effect allele


## Scores
Finally scores are calculate with `plink2 --sscore` which will only compute if the variant ids match, but still computing the score for the correct allele.


# FINNGEN WEIGHTS

We've added in the wdl section the wdl `finngen_weights.wdl` that allows to calculate weights based on the FG sumstats. In order to do so we've build a custom LD panel based on the Finnish panel used for imputation.

In [the scripts folder  link](https://github.com/FINNGEN/CS-PRS-pipeline/tree/master/scripts/panel) there are the scripts used to generate it.

Here's a breakdown of how it works in each step and how to edit the wdl for one's needs

## Global parameters
```
"finngen_weights.test": False,
"finngen_weights.run_scores": False,
"finngen_weights.pheno_list": "gs://path/to/list/of/phenos.txt",
"finngen_weights.plink_root": "gs://finngen-production-library-red/finngen_R8/genotype_plink_1.0/data/finngen_R8_hm3",
"finngen_weights.prefix": "finngen_R?",

```
`Test` mode cuts the input sumstats to only 10k variants and performs the weights calculation in test mode (very few iterations). The output will be useless, but it will run in a very short time (mins vs hours).
`run_scores` determines whether scores are calculated or not.
`pheno_list` is the path to the list of phenos to analyze (see munging step).
`plink_root` is the base of the .bim file used as validation by cs-prs. If scores are passed, it also needs to have .bed and .fam files in the same directory
`prefix'`is the prefix prepended to the output of all files


## Munging

The task `munge_sumstats` munges the input data to a CHROM_POS_REF_ALT format needed with FG.

Relevant inputs:
```
"finngen_weights.munge_sumstats.file_root": "gs://finngen-production-library-green/finngen_R5/finngen_R5_analysis_data/summary_stats/release/finngen_R5_PHENO.gz",
"finngen_weights.munge_sumstats.columns":"#chrom,pos,ref,alt,beta,pval",
```
The first is the path to the input sumsats. The string `PHENO` is replaced with the phenotypes defined in the global parameter described before.

The columns are the core of the munging step. These column configurations need to be in the order that PRCSs expects them i.e chrom pos ref alt beta pvalue.

## Weights

```
"finngen_weights.weights.ref_list": "gs://finngen-production-library-green/prs/ref_list_fin.txt",
"finngen_weights.weights.rsid_map": "gs://finngen-production-library-green/finngen_R8/finngen_R8_analysis_data/variant_mapping/finngen_R8_hm3.rsid.map.tsv",
 ```

`ref_list` is the list of files containing the custom LD panel built with thes scripts described in the introduction.

`rsid_map` is a filed structures as follows:
```
rs11596870 10_100000235
rs11190363 10_100002628
rs7902856 10_100004827
rs7078766 10_100005136
[...]
```
That is used to output weights in rsid format.


This task outputs three files:
- weights in CHR_POS_REF_ALT format
- weights in rsid format
- the log output of CS-PRS for checking all the relevant metadata


## Scores
If interested in calculating scores, the only edit one can make is to increase/decrease the number of cpus in the `.json`. The .bed file is determined by the `plink_root` parameter passed in the global parameters.