import pandas as pd
import pickle
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import seaborn as sns


def sci(x, digits=1):
    # mathtext "3.2 \times 10^{-5}" instead of "3.2e-05"
    mant, exp = f'{x:.{digits}e}'.split('e')
    return rf'{mant} \times 10^{{{int(exp)}}}'


with open('../data/coding_repeat_reset_2.pkl','rb') as file:
    coding_repeat=pickle.load(file)

with open('../data/non_coding_repeat_reset_2.pkl','rb') as file:
    non_coding_repeat=pickle.load(file)

    
with open('../data/cds_dict_human_reset.pkl', 'rb') as file:
    cds_dict = pickle.load(file)

with open('../data/noncoding_dict_human_reset.pkl', 'rb') as file:
    non_cds_dict = pickle.load(file)

with open('../data/gene_tpm.pkl','rb') as f:
    gene_tpm=pickle.load(f)

with open('../data/fs_ratio.pkl','rb') as f:
    fs_ratio=pickle.load(f)




codontable = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'
}



def plot_freq_exp(df_gene,column2,ylim=None):
    #pcc=stats.pearsonr(df_gene['tpm'],df_gene[column2])
    scc=stats.spearmanr(df_gene['tpm'],df_gene[column2],alternative='less')
    print(scc)
    plt.figure(figsize=(8, 6))
    plt.scatter(df_gene['tpm'],df_gene[column2], c='black', alpha=0.7, s=5)
    plt.ylabel('Codon repeat frequency in coding regions\nminus that in noncoding regions',fontsize=12)
    plt.xlabel('Gene expression level',fontsize=12)
    if ylim is None:
        y_min, y_max = df_gene[column2].min(), df_gene[column2].max()
        y_pad = 0.05 * (y_max - y_min)
        ylim = (y_min - y_pad, y_max + y_pad)
    plt.ylim(ylim)
    #plt.figtext(0.15, 0.83, f'pcc={pcc[0]:.2f}, p_val={pcc[1]:.3e}', fontsize=10)
    plt.figtext(0.55, 0.3, f'ρ = {scc[0]:.2f}, $P$ = ${sci(scc[1], 3)}$', fontsize=12)
    plt.savefig(f'../figures/freq_kb_exp/frequency_tpm_{column2}_reset_rescaley_2.png')
    plt.savefig(f'../figures/freq_kb_exp/frequency_tpm_{column2}_reset_rescaley_2.svg')




def add_freq_5(df_gene, threshold=5):
    df_gene[f'coding_{threshold}%']=0
    df_gene[f'noncoding_{threshold}%']=0
    for key in coding_repeat.keys():
        count=0
        for motif in coding_repeat[key].keys():
            if codontable[motif[:3]]!='*' and fs_ratio[motif[:3]]>threshold:
                if fs_ratio[motif[:3]]>threshold:
                    count+=coding_repeat[key][motif]
        df_gene.loc[key,f'coding_{threshold}%']=count*1000/len(cds_dict[key])
    for key in non_coding_repeat.keys():
        count=0
        for motif in non_coding_repeat[key].keys():
            if codontable[motif[:3]]!='*' and fs_ratio[motif[:3]]>threshold:
                if fs_ratio[motif[:3]]>threshold:
                    count+=non_coding_repeat[key][motif]
        df_gene.loc[key,f'noncoding_{threshold}%']=count*1000/len(non_cds_dict[key])
    df_gene[f'minus_{threshold}%']=df_gene[f'coding_{threshold}%']-df_gene[f'noncoding_{threshold}%']
    df_gene.replace([np.inf, -np.inf], np.nan,inplace=True)
    df_gene.dropna(inplace=True)
    return df_gene



def add_count(df_gene):
    df_gene['coding_count']=0
    df_gene['noncoding_count']=0
    for key in coding_repeat.keys():
        count=0
        for motif in coding_repeat[key].keys():
            if codontable[motif[:3]]!='*':
                count+=coding_repeat[key][motif]
        df_gene.loc[key,'coding_count']=count
    for key in non_coding_repeat.keys():
        count=0
        for motif in non_coding_repeat[key].keys():
            if codontable[motif[:3]]!='*':
                count+=non_coding_repeat[key][motif]
        df_gene.loc[key,'noncoding_count']=count
    df_gene.replace([np.inf, -np.inf], np.nan,inplace=True)
    df_gene.dropna(inplace=True)
    return df_gene


def add_count_5(df_gene, threshold=5):
    df_gene[f'coding_{threshold}%_count']=0
    df_gene[f'noncoding_{threshold}%_count']=0
    for key in coding_repeat.keys():
        count=0
        for motif in coding_repeat[key].keys():
            if codontable[motif[:3]]!='*' and fs_ratio[motif[:3]]>threshold:
                if fs_ratio[motif[:3]]>threshold:
                    count+=coding_repeat[key][motif]
        df_gene.loc[key,f'coding_{threshold}%_count']=count
    for key in non_coding_repeat.keys():
        count=0
        for motif in non_coding_repeat[key].keys():
            if codontable[motif[:3]]!='*' and fs_ratio[motif[:3]]>threshold:
                if fs_ratio[motif[:3]]>threshold:
                    count+=non_coding_repeat[key][motif]
        df_gene.loc[key,f'noncoding_{threshold}%_count']=count
    df_gene.replace([np.inf, -np.inf], np.nan,inplace=True)
    df_gene.dropna(inplace=True)
    return df_gene

# raw number
df_tpm=pd.DataFrame.from_dict(gene_tpm,orient='index',columns=['tpm'])
df_tpm['coding_len']=0
df_tpm['noncoding_len']=0
drop=[]
for key in df_tpm.index:
    if key in coding_repeat.keys():
        df_tpm.loc[key,'coding_len']=len(cds_dict[key])
    else:
        drop.append(key)

df_tpm.drop(drop,inplace=True)


drop=[]
for key in df_tpm.index:
    if key in non_coding_repeat.keys():
        df_tpm.loc[key,'noncoding_len']=len(non_cds_dict[key])
    else:
        drop.append(key)

df_tpm.drop(drop,inplace=True)


for threshold in [0.5, 1, 5]:
    # copy: add_count_5 edits in place, so without it tpm would be log10'd again every threshold
    df_5=add_count_5(df_tpm.copy(), threshold=threshold)
    df_5['tpm']=np.log10(df_5['tpm']+1)
    df_5_bin=add_count_5(df_5, threshold=threshold)
    bin_number=15
    df_5_bin['Bin']=pd.qcut(df_5_bin['tpm'],bin_number)   # equal size bins
    bins=df_5_bin['Bin'].unique()
    dfs = {i: df_5[df_5_bin['Bin'] == i].drop(columns=['Bin']) for i in bins}
    number_of_bins=len(dfs.keys())
    x=np.zeros(number_of_bins)
    y=np.zeros(number_of_bins)
    for i in range(number_of_bins):
        x[i]=np.median(dfs[list(dfs.keys())[i]]['tpm'])
        bin_count=sum(dfs[list(dfs.keys())[i]][f'coding_{threshold}%_count'])
        bin_len=sum(dfs[list(dfs.keys())[i]]['coding_len'])
        bin_count_noncoding=sum(dfs[list(dfs.keys())[i]][f'noncoding_{threshold}%_count'])
        bin_len_noncoding=sum(dfs[list(dfs.keys())[i]]['noncoding_len'])
        y[i]=bin_count*1000/bin_len-bin_count_noncoding*1000/bin_len_noncoding


    # binning
    scc=stats.spearmanr(x,y,alternative='less')
    plt.figure(figsize=(8, 6))
    plt.scatter(x,y, c='black', alpha=0.7, s=7)
    plt.ylabel('Codon repeat frequency in coding regions\nminus that in noncoding regions', fontsize=15)
    plt.xlabel(f'Gene expression level', fontsize=15)
    plt.tick_params(labelsize=13)
    plt.figtext(0.5,0.75, f'ρ = {scc[0]:.2f}, $P$ = ${sci(scc[1], 3)}$', fontsize=15)
    plt.tight_layout()
    plt.savefig(f'../figures/freq_kb_exp/minus_{threshold}_binned_{bin_number}_reset.svg')
    plt.savefig(f'../figures/freq_kb_exp/minus_{threshold}_binned_{bin_number}_reset.png')
    plt.close()


    # correlation (df_5 already has log10 tpm, same as before)
    df_plot=add_freq_5(df_5,threshold)
    plot_freq_exp(df_plot,f'minus_{threshold}%')
    plt.close()
