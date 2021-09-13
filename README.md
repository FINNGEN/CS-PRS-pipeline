# CS-PRS-pipeline

Pipeline to calculates PRS based on a list of sumstats.
Weights are calculated with PRScs: https://github.com/getian107/PRScs

The sumstats used for R7 are the following:
```|filename                                                            |pheno                     |n_total|n_cases|n_ctrls|publication                                                                                           |ancestry      |finngen_phenocode|
|--------------------------------------------------------------------|--------------------------|-------|-------|-------|------------------------------------------------------------------------------------------------------|--------------|-----------------|
|AD_sumstats_Jansenetal.txt.gz                                       |G6_ALZHEIMER              |455258 |71880  |383378 |https://www.nature.com/articles/s41588-018-0311-9                                                     |European      |G6_ALZHEIMER     |
|adhd_eur_jun2017.gz                                                 |F5_ADHD                   |53293  |19099  |34194  |https://www.biorxiv.org/content/early/2017/06/03/145581                                               |European      |F5_ADHD          |
|alsMetaSummaryStats_march21st2018.tab.gz                            |G6_ALS                    |80610  |20806  |59804  |https://www.cell.com/neuron/references/S0896-6273(18)30148-X                                          |European      |G6_ALS           |
|Bipolar_vs_control_PGC_Cell_2018_formatted.txt.gz                   |F5_BIPO                   |74194  |20129  |54065  |PGC;Cell;2018;https://doi.org/10.1016/j.cell.2018.05.046                                              |European      |F5_BIPO          |
|cad.add.160614.website.txt.gz                                       |I9_CHD                    |184305 |60801  |123504 |https://www.nature.com/articles/ng.3396                                                               |Trans-ancestry|I9_CHD           |
|ckqny.scz2snpres.gz                                                 |F5_SCHZPHR                |150064 |36989  |113075 |https://www.nature.com/articles/nature13595                                                           |Trans-ancestry|F5_SCHZPHR       |
|daner_PGC_SCZ52_0513a.resultfiles_PGC_SCZ52_0513.sh2_nofin.gz       |F5_SCHZPHR                |150064 |36989  |113075 |PGC_schzizophrenia_GWAS_Ripke_et_al._excluding_Finnish                                                |European      |F5_SCHZPHR       |
|dpw_excludingFinnishStudies_4finngen.txt.gz                         |alcohol_consumption       |527282 |NA     |NA     |https://www.nature.com/articles/s41398-019-0676-2                                                     |European      |NA               |
|EAGLE_AD_GWAS_results_2015.txt.gz                                   |L12_ATOPIC                |116863 |21399  |95464  |https://www.nature.com/articles/ng.3424                                                               |Trans-ancestry|L12_ATOPIC       |
|Educational_Attainment_excl23andme_Lee_2018_NatGen_formatted.txt.gz |educational_attainment    |766343 |NA     |NA     |Lee;Nat_Gen;2018;https://doi.org/10.1038/s41588-018-0147-3                                            |European      |NA               |
|EUR.CD.gwas_info03_filtered.assoc.gz                                |K11_CROHN                 |20883  |5956   |14927  |https://www.nature.com/articles/ng.3359                                                               |European      |K11_CROHN        |
|EUR.IBD.gwas_info03_filtered.assoc.gz                               |K11_IBD                   |34652  |12882  |21770  |https://www.nature.com/articles/ng.3359                                                               |European      |K11_IBD          |
|EUR.UC.gwas_info03_filtered.assoc.gz                                |K11_ULCER                 |27432  |6968   |20464  |https://www.nature.com/articles/ng.3359                                                               |European      |K11_ULCER        |
|focal_epilepsy_METAL_4finngen.txt.gz                                |FE                        |39348  |9671   |29677  |https://www.nature.com/articles/s41467-018-07524-z                                                    |European      |FE               |
|gabriel_asthma_meta-analysis_36studies_format_repository_NEJM.txt.gz|J10_ASTHMA                |26475  |10365  |16110  |https://www.nejm.org/doi/full/10.1056/nejmoa0906312                                                   |European      |J10_ASTHMA       |
|generalised_epilepsy_METAL_4finngen.txt.gz                          |GE                        |33446  |3769   |29677  |https://www.nature.com/articles/s41467-018-07524-z                                                    |European      |GE               |
|GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz     |Height                    |253288 |NA     |NA     |https://www.nature.com/articles/ng.3097                                                               |European      |NA               |
|GWAS_CP_all.txt.gz                                                  |Cognitive_performance_meta|257828 |NA     |NA     |https://www.nature.com/articles/s41588-018-0147-3                                                     |European      |NA               |
|Hb_gwas_summary_fromNealeLab.tsv.gz                                 |Haemoglobin_concentration |350474 |NA     |NA     |UKBB                                                                                                  |European      |NA               |
|HbA1c_METAL_European.txt.gz                                         |HbA1c                     |123665 |NA     |NA     |http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1002383                         |European      |NA               |
|HDL_Meta_ENGAGE_1000G_non_FINRISK_QCd.txt.gz                        |HDL-C                     |57045  |NA     |NA     |https://www.nature.com/articles/ng.3300                                                               |European      |NA               |
|IGAP_stage_1.txt.gz                                                 |G6_ALZHEIMER              |54162  |17008  |37154  |https://www.nature.com/articles/ng.2802                                                               |European      |G6_ALZHEIMER     |
|iPSYCH-PGC_ASD_Nov2017.gz                                           |F5_AUTISM                 |46350  |18381  |27969  |https://www.biorxiv.org/content/early/2017/11/27/224774                                               |European      |F5_AUTISM        |
|kunkle_etal_stage1.txt.gz                                           |G6_ALZHEIMER              |63926  |21982  |41944  |https://www.niagads.org/datasets/ng00075                                                              |European      |G6_ALZHEIMER     |
|LDL_Meta_ENGAGE_1000G_non_FINRISK_QCd.txt.gz                        |LDL-C                     |57045  |NA     |NA     |https://www.nature.com/articles/ng.3300                                                               |European      |NA               |
|Mahajan.NatGenet2018b.T2D.European.txt.gz                           |E4_DM2                    |898130 |74124  |824006 |https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6287706/                                                 |European      |E4_DM2           |
|MDD2018_ex23andMe.19fields.gz                                       |F5_DEPRESSIO              |173005 |59851  |113154 |https://www.nature.com/articles/s41588-018-0090-3                                                     |European      |F5_DEPRESSIO     |
|meta_v3_onco_euro_overall_ChrAll_1_release.txt.gz                   |C3_PROSTATE               |140254 |79148  |61106  |https://www.nature.com/articles/s41588-018-0142-8                                                     |European      |C3_PROSTATE      |
|METAANALYSIS_DIAGRAM_SE1.txt.gz                                     |E4_DM2                    |159208 |26676  |132532 |http://diabetes.diabetesjournals.org/content/66/11/2888                                               |European      |E4_DM2           |
|metastroke.all.chr.bp.gz                                            |I9_STR                    |74393  |12389  |62004  |https://www.thelancet.com/journals/laneur/article/PIIS1474-4422(12)70234-X/abstract                   |European      |I9_STR           |
|MTAG_CP.to10K.txt.gz                                                |Cognitive_Performance     |257828 |NA     |NA     |https://www.nature.com/articles/s41588-018-0147-3                                                     |European      |NA               |
|MTAG_EA.to10K.txt.gz                                                |Educational_attainment    |1131881|NA     |NA     |https://www.nature.com/articles/s41588-018-0147-3                                                     |European      |NA               |
|pgc.bip.full.2012-04.txt.gz                                         |F5_BIPO                   |16731  |7481   |9250   |https://www.nature.com/articles/ng.943                                                                |European      |F5_BIPO          |
|pgc.ed.freeze1.summarystatistics.July2017.txt.gz                    |F5_ANOREX                 |14477  |3495   |10982  |Missing                                                                                               |European      |F5_ANOREX        |
|RA_GWASmeta_European_v2.txt.gz                                      |M13_RHEUMA                |58284  |14361  |43923  |https://www.nature.com/articles/nature12873                                                           |European      |M13_RHEUMA       |
|SavageJansen_2018_intelligence_metaanalysis_formatted.txt.gz        |intelligence              |269867 |NA     |NA     |Savage;Nat_Gen-2018;https://doi.org/10.1038/s41588-018-0152-6                                         |European      |NA               |
|Saxena_fullUKBB_Longsleep_summary_stats_formatted.txt.gz            |long_sleep                |446118 |NA     |NA     |Dashti;Natcomm-2019;https://doi.org/10.1038/s41467-019-08917-4                                        |European      |NA               |
|Shrine_30804560_FEV1_meta-analysis.txt.gz                           |FEV1                      |321047 |NA     |NA     |https://www.nature.com/articles/s41588-018-0321-7                                                     |European      |NA               |
|Shrine_30804560_FEV1_to_FVC_RATIO_meta-analysis.txt.gz              |FEV1/FVC                  |321047 |NA     |NA     |https://www.nature.com/articles/s41588-018-0321-7                                                     |European      |NA               |
|Shrine_30804560_FVC_meta-analysis.txt.gz                            |FVC                       |321047 |NA     |NA     |https://www.nature.com/articles/s41588-018-0321-7                                                     |European      |NA               |
|Shrine_30804560_PEF_meta-analysis.txt.gz                            |Peak_expiratory_flow      |321047 |NA     |NA     |https://www.nature.com/articles/s41588-018-0321-7                                                     |European      |NA               |
|SNP_gwas_mc_merge_nogc.tbl.uniq.gz                                  |BMI                       |322154 |NA     |NA     |https://www.nature.com/articles/nature14177                                                           |European      |NA               |
|SORTED_PTSD_EA9_ALL_study_specific_PCs1.txt.gz                      |F5_PTSD                   |9537   |2424   |7113   |https://www.nature.com/articles/mp201777                                                              |European      |F5_PTSD          |
|sumstats_neuroticism_ctg_format_formatted.txt.gz                    |neuroticism               |449484 |NA     |NA     |Nagel;NatGen2018-https://doi.org/10.1038/s41588-018-0151-7                                            |European      |NA               |
|TC_Meta_ENGAGE_1000G_non_FINRISK_QCd.txt.gz                         |Total-cholesterol         |57045  |NA     |NA     |https://www.nature.com/articles/ng.3300                                                               |European      |NA               |
|TG_Meta_ENGAGE_1000G_non_FINRISK_QCd.txt.gz                         |Triglycerides             |57045  |NA     |NA     |https://www.nature.com/articles/ng.3300                                                               |European      |NA               |
|UKBB_gwas_Neale_Chronotype_formatted.txt.gz                         |chronotype                |301143 |NA     |NA     |UKBB                                                                                                  |European      |NA               |
|UKBB_gwas_Neale_pulserate_formatted.txt.gz                          |heart_rate                |317756 |NA     |NA     |UKBB                                                                                                  |European      |NA               |
|UKBB_gwas_Neale_SleepDuration_formatted.txt.gz                      |sleepduration             |335410 |NA     |NA     |UKBB                                                                                                  |European      |NA               |
|UKBB_gwas_Neale_SleeplessnessInsomnia_formatted.txt.gz              |insomnia                  |336965 |NA     |NA     |UKBB                                                                                                  |European      |NA               |
|UKBB_HYPO.txt.gz                                                    |HYPOTHYROIDISM            |406300 |14871  |391429 |UKBB                                                                                                  |European      |HYPOTHYROIDISM   |
|UKB-ICBPmeta750k_DBPsummaryResults.txt.gz                           |I9_HYPTENSESS             |757601 |NA     |NA     |UKBB                                                                                                  |European      |I9_HYPTENSESS    |
|UKB-ICBPmeta750k_PPsummaryResults.edited.txt.gz                     |I9_HYPTENSESS             |757601 |NA     |NA     |UKBB                                                                                                  |European      |I9_HYPTENSESS    |
|UKB-ICBPmeta750k_SBPsummaryResults.txt.gz                           |I9_HYPTENSESS             |757601 |NA     |NA     |UKBB                                                                                                  |European      |I9_HYPTENSESS    |
|GSCAN_AgeofInitiation.txt.gz                                        |smoking_age_of_initiation |341427 |NA     |NA     |https://www.nature.com/articles/s41588-018-0307-5                                                     |European      |NA               |
|GSCAN_CigarettesPerDay.txt.gz                                       |smoking_cigarettes_per_day|337334 |NA     |NA     |https://www.nature.com/articles/s41588-018-0307-5                                                     |European      |NA               |
|GSCAN_DrinksPerWeek.txt.gz                                          |drinks_per_week           |941280 |NA     |NA     |https://www.nature.com/articles/s41588-018-0307-5                                                     |European      |NA               |
|GSCAN_SmokingInitiation.txt.gz                                      |smoking_initiation        |1232091|NA     |NA     |https://www.nature.com/articles/s41588-018-0307-5                                                     |European      |NA               |
|PGC_UKB_depression_genome-wide.txt.gz                               |depression                |500199 |170756 |329443 |https://www.nature.com/articles/s41588-018-0090-3                                                     |European      |F5_DEPRESSIO     |
|PGC_UKB_23andMe_depression_10000.txt.gz                             |depression_23andme        |807553 |246363 |561190 |https://www.nature.com/articles/s41588-018-0090-3                                                     |European      |F5_DEPRESSIO     |
|PGC3_SCZ_wave3_public_without_frequencies.v1.tsv.gz                 |schizophrenia             |306011 |69369  |236642 |https://www.medrxiv.org/content/10.1101/2020.09.12.20192922v1                                         |European      |F5_SCHZPHR       |
|PGC3_SCZ_wave3_public_without_frequencies.clumped.v1.tsv.gz         |schizophrenia_clumped     |306011 |69369  |236642 |https://www.medrxiv.org/content/10.1101/2020.09.12.20192922v1                                         |European      |F5_SCHZPHR       |
|daner_PGC_BIP32b_mds7a_0416a.gz                                     |bipolar_disorder          |51710  |20352  |31358  |https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6956732/                                                 |European      |F5_BIPO          |
|CIMBA_BRCA1_BCAC_TN_meta_summary_level_statistics.txt.gz            |breast_cancer             |247173 |133384 |113789 |https://www.nature.com/articles/s41588-020-0609-2                                                     |European      |C3_BREAST        |
|AF_HRC_GWAS_ALLv11.txt.gz                                           |atrial_fibrillation       |588190 |65446  |522744 |https://www.nature.com/articles/s41588-018-0133-9                                                     |Trans-ancestry|I9_AF            |
|lifegen_phase2_bothpl_alldr_2017_09_18.tsv.gz                       |lifespan                  |1012240|NA     |NA     |https://elifesciences.org/articles/39856                                                              |European      |NA               |
|RISK_GWAS_MA_UKB+23andMe+replication.txt.gz                         |risk_tolerance            |975353 |NA     |NA     |https://www.nature.com/articles/s41588-018-0309-3                                                     |European      |NA               |
|chronic_pain-bgen.stats.gz                                          |multisite_chronic_pain    |380000 |NA     |NA     |https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6592570/                                                 |European      |NA               |
|Saxena_fullUKBB_Insomnia_summary_stats.txt.gz                       |insomnia                  |453379 |NA     |NA     |https://www.nature.com/articles/s41588-019-0361-7                                                     |European      |F5_INSOMNIA      |
|continuous-LDLC-both_sexes-medadj_irnt.tsv.gz                       |ldl_adj_irnt              |248559 |NA     |NA     |https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=511623409|European      |NA               |
|biomarkers-30750-both_sexes-irnt.tsv.gz                             |hba1c_inrt                |194273 |NA     |NA     |https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=511623409|European      |NA               |
|biomarkers-30870-both_sexes-irnt.tsv.gz                             |triglycerides_irnt        |194702 |NA     |NA     |https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=511623409|European      |NA               |
|MEGASTROKE.1.AS.TRANS.out.gz                                        |any_stroke                |524354 |67162  |454450 |https://www.nature.com/articles/s41588-018-0058-3                                                     |Trans-ancestry|I9_STR           |
|UKBB.asthma-2.assoc.gz                                              |asthma                    |536345 |88486  |447859 |https://www.nature.com/articles/s41467-020-15649-3                                                    |European      |J10_ASTHMA       |
|ICP_FG_DC_EG_DBDS_202012.gz                                         |ICP_WIDE_EXMORE           |2138   |301733 |303871 |FG                                                                                                    |Scandinavian  |ICP_WIDE_EXMORE  |

```

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
