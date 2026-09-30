import numpy as np
import pandas as pd
from scipy import stats
import lifelines,argparse,pylab
from matplotlib import pyplot as plt
import matplotlib as mpl
from lifelines import NelsonAalenFitter,CoxPHFitter,KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.plotting import add_at_risk_counts
from utils import file_exists,make_sure_path_exists,pretty_print,mapcount
import os.path
import seaborn as sns
from sklearn import metrics
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()

sns.set(palette='Set2', font_scale=1.0)

def calculate_AUC(data,scores,out_path,phenocode,tag):
   
 
    pretty_print("PLOTTING")
    pylab.ioff()
    tag = f"{phenocode}_{tag}"
    finalFigPath = os.path.join(out_path, f'{tag}_AUC.pdf')
    
    print(finalFigPath)
    pheno_data = pd.merge(scores,data,on='FINNGENID') #calcuate duration for cox
    # SET UP FIG
    fig = plt.figure()
    gs = mpl.gridspec.GridSpec(2, 1)
    ax = fig.add_subplot(gs[0, 0] )

    y_data = np.array(pheno_data.PHENO)
    y_pred = np.array(pheno_data.SCORE)
    fpr, tpr, _ = metrics.roc_curve(y_data,y_pred)
    auc = round(metrics.roc_auc_score(y_data,y_pred),3)

    #create ROC curve
    ax.plot(fpr,tpr,label="AUC="+str(auc))
    ax.set_ylabel('True Positive Rate')
    ax.set_xlabel('False Positive Rate')
    ax.legend(loc=4)
    ax.set_title(tag, fontsize = 10)

    # now do r2
    ax2 = fig.add_subplot(gs[1, 0] )
    
    scaled_data = scaler.fit_transform(y_pred.reshape(-1, 1))
    model = sm.Logit(y_data, sm.add_constant(scaled_data))
    results = model.fit()
    print(results.summary())

    # for plotting only use 10k samples else it crashes
    idx = np.random.randint(len(y_data),size =10000)
    sample_prs = scaled_data[idx]
    sample_pheno = y_data[idx]
    sample_pred = results.predict(sm.add_constant(sample_prs))

    r2 = round(results.prsquared,3)
    ax2.scatter(sample_prs,sample_pred,s=1,label =f"pseudo r2:{r2}" ) 
    ax2.scatter(sample_prs,sample_pheno,s=1)
    ax2.legend(loc='center left')
    ax2.set_ylabel('Pheno')
    ax2.set_xlabel('Renormed prs')
    fig.savefig(finalFigPath)
    plt.close(fig)

    log_path = finalFigPath.replace('.pdf','.log')
    with open(log_path,'wt') as o: o.write('\t'.join([phenocode,tag,str(auc),str(r2)]) +'\n')
 

def get_prs(score_file):
    """
    Returns pandas data frame with z score of PRS.
    """
    pretty_print('returning PRS...')
    scores = pd.read_csv(score_file,usecols =['IID','SCORE1_AVG'],sep = '\t')
    scores.columns = ['FINNGENID','SCORE']
    scores['z'] = stats.zscore(scores['SCORE'])
    scores['percentile'] = pd.qcut(scores.SCORE,100,labels=False)
    return scores

def get_age_data(out_path,age_file,pheno = "G6_ALZHEIMER",onset_suffix="_FU_AGE",test=False):
    """
    Returns the pandas data frame with the following columns:
    PHENO: --> boolean of case/controls
    PHENO_ONSET : --> age of onset for cases, last age for controls
    DEATH : boolean of death
    DEATH_AGE: age at death (or last age if not dead)
    """
    pretty_print('returning AGE...')
    save_path = os.path.join(out_path,'age_data')
    make_sure_path_exists(save_path)
    out_file = os.path.join(save_path,pheno +'.age')

    pheno_column,age_column,death_column,death_age_column  = pheno,f"{pheno}{onset_suffix}",'DEATH',f'DEATH{onset_suffix}'
    nrows = 1000 if args.test else None
    if os.path.isfile(out_file):
        print('loading...')
        data = pd.read_csv(out_file)
        return data

    else:
        print('generating...')
        pheno_data = pd.read_csv(age_file,usecols = ['FINNGENID',age_column,pheno_column],nrows=nrows,sep='\t').rename(columns = {pheno:"PHENO",age_column:"PHENO_ONSET"}).dropna()
        death_data = pd.read_csv(age_file,usecols = ['FINNGENID',death_column,death_age_column],nrows=nrows,sep='\t').rename(columns = {death_column:"DEATH",death_age_column:"DEATH_AGE"}).dropna()
        # remove inconsistent values
        data = pd.merge(pheno_data,death_data,on='FINNGENID')

        data = data[(data['DEATH_AGE'] - data['PHENO_ONSET']) >= 0]
        print('saving...')
        data.to_csv(out_file,index=False)
        return data

def plot_onset(data,scores,out_path,phenocode,tag):
    """
    Survival analysis where we focus on time to death from diagnosis only for cases.
    """
    pretty_print("PLOTTING")
    pylab.ioff()
    tag = f"{phenocode}_{tag}"
    finalFigPath = os.path.join(out_path, f'{tag}_age_onset.pdf')
    print(finalFigPath)

    # merge scores data with pheno data
    pheno_data =   pd.merge(scores,data,on='FINNGENID') #calcuate duration for cox
    pheno_data['DELTA_AGE'] = pheno_data.DEATH_AGE -    pheno_data.PHENO_ONSET
    pheno_data =  pheno_data[pheno_data.PHENO==1][['DEATH','DELTA_AGE','z','PHENO_ONSET']]
    cph = CoxPHFitter() #DELTA_AGE is duration, DEATH is event, z and    pheno_onset covariates
    cph.fit(pheno_data, duration_col='DELTA_AGE', event_col='DEATH', show_progress=True)
    cph.summary.to_csv(finalFigPath.replace('pdf','csv'),sep='\t')

    # SET UP FIG
    fig = plt.figure()
    gs = mpl.gridspec.GridSpec(2, 1)

    ax1 = fig.add_subplot(gs[0, 0] )
    ax1.set_xlim([0,pheno_data.DELTA_AGE.max()])

    quantiles = [scores.z.quantile(q) for q in [0.1,0.3,0.5,0.7,0.9]]

    cph.plot_covariate_groups('z',quantiles, cmap='coolwarm',ax=ax1)

    #pheno_data['DECILE'] = 10-
    #np.digitize(pheno_data.z,np.percentile(pheno_data.z,np.linspace(10,100,10)))
    #cph.plot_covariate_groups('DECILE',list(set(pheno_data.DECILE)),
    #cmap='coolwarm',ax=ax)

    ax2 = fig.add_subplot(gs[1, 0] )
    cph.plot(ax=ax2)

    plt.suptitle(tag, fontsize = 10)
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)
    fig.savefig(finalFigPath)
    plt.close(fig)


def plot_percentile_risk(data,scores,out_path,phenocode,tag):

    '''
    Plots the risk of being diagnosed the disease for each percentile of PRS
    '''
    pretty_print("PLOTTING")
    pylab.ioff()
    tag = f"{phenocode}_{tag}"
    finalFigPath = os.path.join(out_path, f'{tag}_risk.pdf')
    print(finalFigPath)

    data = pd.merge(scores,data,on='FINNGENID')

    # get average age of diagnosis:
    avg = data[data.PHENO==1].PHENO_ONSET.mean()
    std = data[data.PHENO==1].PHENO_ONSET.std()
    age_bins = [avg-std,avg,avg+std]
    percentiles = np.linspace(0,99,100,dtype = int)
    plot_data = np.empty((len(age_bins),len(percentiles),2))
    for i,p in enumerate(percentiles) :
        q_data = data[(data.percentile == p)]
        #kmf = KaplanMeierFitter()
        naf = NelsonAalenFitter()
        naf.fit(q_data.PHENO_ONSET, event_observed = q_data.PHENO, label = p)
        #kmf.fit(q_data.PHENO_ONSET, event_observed = q_data.PHENO, label = p)
        for j,age in enumerate(age_bins):
            age_index = np.abs(naf.timeline-age).argmin()
            lower,upper = naf.confidence_interval_.to_numpy()[age_index]
            plot_data[j,i] = lower,upper

    # SET UP FIG
    print(plot_data.shape)
    fig = plt.figure()
    gs = mpl.gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0, 0] )
    for i,entry in enumerate(plot_data):
        upper,lower = entry.T
        ax.fill_between(percentiles,lower,upper, alpha=0.35)
        ax.plot(percentiles,(upper+lower)/2,label = round(age_bins[i],2))

    ax.set_title(tag, fontsize = 10)
    ax.set_xlabel('Percentile')
    ax.set_ylabel('Risk at AGE')
    ax.legend(loc = 'upper left')
    fig.savefig(finalFigPath)
    plt.close(fig)

def plot_survival(data,scores,out_path,phenocode,tag,quantiles = [(0,0.1,'low PRS'),(0.9,1,'high PRS')]):
    '''
    Plots the KM and NA estimates
    '''
    pretty_print("PLOTTING")
    pylab.ioff()
    tag = f"{phenocode}_{tag}"
    finalFigPath = os.path.join(out_path, f'{tag}_survival.pdf')
    print(finalFigPath)

    # SET UP FIG
    # height_ratios: ax3 (Cox forest plot, a single 'z' covariate) needs far less vertical
    # space than the two survival-curve panels -- equal thirds left most of it blank
    fig = plt.figure(figsize=(8, 9))
    gs = mpl.gridspec.GridSpec(3, 1, height_ratios=[4, 4, 1])

    ax1 = fig.add_subplot(gs[0, 0] )
    ax2 = fig.add_subplot(gs[1, 0] )

    colors = ['blue','red']

    ax1.set_ylabel(r"$\hat{S}(t)$", fontsize = 8)
    ax2.set_ylabel(r"$\hat{\Lambda} (t)$", fontsize = 8)
    ax2.set_xlabel('Age')

    # SURVIVAL CURVES

    # lists where I store data needed for the logRank test and the printing of at_risk_counts
    logrankData = []
    kmfList = []

    for i,q in enumerate(quantiles) :
        lower,upper,q_label = q
        # filter samples to only top and bottom quantiles
        q_scores= scores[(scores.z > scores.z.quantile(lower)) & (scores.z < scores.z.quantile(upper))]
        q_data = pd.merge(q_scores,data,on='FINNGENID')
        durations,observed = q_data.PHENO_ONSET,q_data.PHENO

        # append data
        logrankData.extend((durations,observed))

        #fit
        #import fitters
        kmf,naf = KaplanMeierFitter(),NelsonAalenFitter()

        kmf.fit(durations, event_observed = observed, label = 'Q ' +str(q_label),ci_labels = ['lower bound','upper bound'])
        naf.fit(durations, event_observed = observed, label = 'Q ' +str(q_label))

        #plot
        naf.plot(ax = ax2, color = colors[i])
        kmf.plot(ax = ax1, color = colors[i])
        kmfList.append(kmf)


    add_at_risk_counts(kmfList[0],kmfList[1], ax=ax1)

    results = logrank_test(logrankData[0],logrankData[2],logrankData[1],logrankData[3], alpha=.99)
    results.print_summary()

    ax1.set_xlabel('')
    ax1.legend(loc = 'lower left')
    ax2.set_xlabel('Age')
    ax2.legend(loc = 'upper left')

    for ax in [ax1,ax2]:
        ax.set_ylim([0,1])
        ax.set_yticks(np.linspace(0, 1,6))

    ax3 = fig.add_subplot(gs[2, 0] )

    # PLOT COX MODEL ESTIMATIONS
    cph = CoxPHFitter()
    cox_data = pd.merge(scores,data,on='FINNGENID')[['z',"PHENO_ONSET","PHENO"]]
    # PHENO_ONSET is duration, PHENO is event, other columns are "covariates"
    cph.fit(cox_data,'PHENO_ONSET','PHENO', show_progress=True)
    cph.summary.to_csv(finalFigPath.replace('pdf','csv'),sep='\t')
    cph.plot(ax=ax3)

    plt.suptitle(tag, fontsize = 10)
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)

    fig.savefig(finalFigPath)
    plt.close(fig)

def main(args):
    print(args)

    make_sure_path_exists(args.out)

    data = get_age_data(args.out,args.age_file,args.pheno,args.onset_suffix,args.test)
    scores = get_prs(args.scores)

    root_name = os.path.basename(args.scores).split('.sscore')[0]
    tag = args.tag if args.tag else root_name

    plot_survival(data,scores,args.out,args.pheno,tag)
    plot_onset(data,scores,args.out,args.pheno,tag)
    plot_percentile_risk(data,scores,args.out,args.pheno,tag)
    calculate_AUC(data,scores,args.out,args.pheno,tag)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description ="Calculation of PRS for summary stats.")

    parser.add_argument('--age_file',type = file_exists,help ='File that contains the age of onset',default = '/mnt/disks/r4/Data/pheno/R4_COV_PHENO_V1_AGE.txt.gz')
    parser.add_argument('--out',type = str, help ='output_path',default = '/mnt/disks/r4/PRS/survival/')
    parser.add_argument('--pheno',type = str, help ='Phenocode',required= True)
    parser.add_argument('--tag',type = str, help ='tag')
    parser.add_argument('--onset_suffix',type = str, help ='suffix to the pheno string to get age of onset',default="_FU_AGE")
    parser.add_argument('--scores',type = file_exists,help = 'score file')
    parser.add_argument('--test',action = 'store_true')

    args = parser.parse_args()
    main(args)
