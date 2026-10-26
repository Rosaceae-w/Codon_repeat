import pandas as pd 
import numpy as np
from Bio import SeqIO
import pickle
import matplotlib.pyplot as plt
import re
from scipy import stats
import math 
from scipy.stats import chisquare
from scipy.stats import fisher_exact



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


with open('../data/cds_dict_human_reset.pkl', 'rb') as file:
    cds_dict = pickle.load(file)

with open('../data/noncoding_dict_human_reset.pkl', 'rb') as file:
    non_cds_dict = pickle.load(file)

with open('../data/coding_repeat_reset_2.pkl','rb') as file:
    coding_repeat=pickle.load(file)

with open('../data/non_coding_repeat_reset_2.pkl','rb') as file:
    non_coding_repeat=pickle.load(file)

with open('../data/pure_codon_repeat_coding_2.pkl','rb') as file:
    pure_codon_repeat_coding=pickle.load(file)

with open('../data/pure_codon_repeat_noncoding_2.pkl','rb') as file:
    pure_codon_repeat_noncoding=pickle.load(file)



df=pd.read_csv('../data/sup_tab_4.csv')
from collections import Counter

# Tally for coding regions
coding_tally = Counter()
for gene in coding_repeat.keys():
    for motif, count in coding_repeat[gene].items():
        if True:
            coding_tally[motif] += count

# Tally for noncoding regions
noncoding_tally = Counter()
for gene in non_coding_repeat.keys():
    for motif, count in non_coding_repeat[gene].items():
        noncoding_tally[motif] += count


df_2 = pd.read_csv('../data/frequency_len_kb_reset_2.tsv', sep='\t')
df['Avg_fs_ratio']=0
df.dropna(inplace=True)
for i in range(len(df)):
    df.loc[i,'+1 frameshifting ratio']=float(df['+1 frameshifting ratio'][i].split('±')[0])
    df.loc[i,'-1 frameshifting ratio']=float(df['-1 frameshifting ratio'][i].split('±')[0])
    df.loc[i,'Avg_fs_ratio']=(df.loc[i,'+1 frameshifting ratio']+df.loc[i,'-1 frameshifting ratio'])/2
    

df_2['+1fs_ratio'] = math.nan
df_2['-1fs_ratio'] = math.nan
df_2['Avg_fs_ratio'] = math.nan
df_2['codon_type']=df_2['codon'].str[:3]
df_2['length']=0
for i in range(len(df_2)):
    codon_repeat=df_2['codon'][i]
    df_2['coding'][i]=coding_tally[codon_repeat] if codon_repeat in coding_tally else 0
    df_2['noncoding'][i]=noncoding_tally[codon_repeat] if codon_repeat in noncoding_tally else 0
    df_2['length'][i] = len(df_2['codon'][i])
    # df_2['length'][i] = len(df_2['codon'][i])
    # df_2['GC%'][i] = (df_2['codon'][i].count('G') + df_2['codon'][i].count('C'))*100/len(df_2['codon'][i])
    if df_2['codon_type'][i] in codontable and codontable[df_2['codon_type'][i]]!='*':
        df_2['+1fs_ratio'][i]=df[df['Inserted codon'].str.contains(df_2['codon_type'][i])]['+1 frameshifting ratio'].values[0]
        df_2['-1fs_ratio'][i]=df[df['Inserted codon'].str.contains(df_2['codon_type'][i])]['-1 frameshifting ratio'].values[0]
        df_2['Avg_fs_ratio'][i]=df[df['Inserted codon'].str.contains(df_2['codon_type'][i])]['Avg_fs_ratio'].values[0]
    else:
        pass

df_2.dropna(inplace=True)

df_2 = df_2.groupby('codon_type', as_index=False).agg({'coding': 'sum', 'noncoding': 'sum', 'Avg_fs_ratio': 'first'})
df_2.sort_values('Avg_fs_ratio',inplace=True)
cds_len=34222315
non_cds_len=920316959
df_2['coding_freq_per_kb'] = df_2['coding'] / (cds_len / 1000)
df_2['noncoding_freq_per_kb'] = df_2['noncoding'] / (non_cds_len / 1000)

cds_avg_len=1758.7786514544146
non_cds_avg_len=47955.654160804545

df_select_5=df_2[df_2['Avg_fs_ratio']>5]
significance = {}
for codon in codontable.keys():
    if codontable[codon]=='*':
        continue
    else:
        coding_total=sum(pure_codon_repeat_coding[codon][i*3] for i in range(3,7))
        noncoding_total=sum(pure_codon_repeat_noncoding[codon][i*3] for i in range(3,7))
        observed=[coding_total,noncoding_total]
        expected=[(coding_total+noncoding_total)*cds_avg_len/(cds_avg_len+non_cds_avg_len),(coding_total+noncoding_total)*non_cds_avg_len/(cds_avg_len+non_cds_avg_len)]
        chi,p=chisquare(observed,expected)
        significance[codon]=p
            #table=np.array([[coding_codon_repeat[motif],non_coding_codon_repeat[motif]], [(coding_codon_repeat[motif]+non_coding_codon_repeat[motif])*cds_avg_len/(cds_avg_len+non_cds_avg_len),(coding_codon_repeat[motif]+non_coding_codon_repeat[motif])*non_cds_avg_len/(cds_avg_len+non_cds_avg_len)]])
            #res = fisher_exact(table, alternative='two-sided')
            #significance[motif]=res.pvalue

significant_codons=[]
for key in significance.keys():
    if significance[key]<0.05:
        significant_codons.append(key)

count=0
c_count = 0
n_count = 0
coding=[]
noncoding=[]

with open('../data/fs_ratio.pkl','rb') as f:
    fs_ratio=pickle.load(f)

for i in range(len(df_2)):
    if df_2['codon_type'][i] in significant_codons and df_2['Avg_fs_ratio'][i] > 0:
        count+=1
        if df_2['coding_freq_per_kb'][i] > df_2['noncoding_freq_per_kb'][i]:
                c_count+=1
                coding.append(df_2['Avg_fs_ratio'][i])
        elif df_2['noncoding_freq_per_kb'][i] > df_2['coding_freq_per_kb'][i]:
                n_count+=1
                noncoding.append(df_2['Avg_fs_ratio'][i])
        elif df_2['noncoding_freq_per_kb'][i] ==df_2['coding_freq_per_kb'][i]:
                print(df_2['codon_type'][i])


mwu=stats.mannwhitneyu(coding,noncoding,alternative='less')
#t_test=stats.ttest_ind(coding, noncoding)
print(mwu)



plt.figure(figsize=(8, 6))
is_sig = df_2['codon_type'].isin(significant_codons)
size = 35
cb=plt.scatter(df_2.noncoding_freq_per_kb[is_sig], df_2.coding_freq_per_kb[is_sig], s=size, c=np.log10(df_2['Avg_fs_ratio'][is_sig]/100), cmap='viridis_r', alpha=0.5, marker='o', label='Significant')

plt.scatter(df_2.noncoding_freq_per_kb[~is_sig], df_2.coding_freq_per_kb[~is_sig], s=size, c=np.log10(df_2['Avg_fs_ratio'][~is_sig]/100), cmap='viridis_r', alpha=0.5, marker='v', label='Non-significant',norm=cb.norm)

y_pad = 0.05 * 0.05
plt.ylim(-y_pad, 0.05)

ax = plt.gca()
xlim = ax.get_xlim()
ylim = ax.get_ylim()
lo = max(xlim[0], ylim[0])
hi = min(xlim[1], ylim[1])
ax.plot([lo, hi], [lo, hi], color='grey', linestyle=':', alpha=0.6, zorder=0)
ax.set_xlim(xlim)
ax.set_ylim(ylim)

cbar = plt.colorbar(label=r'$\log_{10}$(frameshifting ratio)')
vmin, vmax = cb.norm.vmin, cb.norm.vmax
ratios = [m * 10**e for e in range(-4, 3) for m in (1, 2, 5)]
ratios = [r for r in ratios if vmin <= np.log10(r / 100) <= vmax]
cbar.set_ticks([np.log10(r / 100) for r in ratios])
cbar.set_ticklabels([f'{np.log10(r / 100):.2g}\n({r:g}%)' for r in ratios])
cbar.ax.tick_params(labelsize=12)
cbar.set_label(r'$\log_{10}$(frameshifting ratio)', fontsize=15)
#plt.colorbar(cb)
ax.tick_params(labelsize=13)
plt.xlabel('Number of repeats per kb in noncoding regions',fontsize=15)
plt.ylabel('Number of codon repeats per kb in coding regions',fontsize=15)
#plt.title('All frameshifting repeats')
plt.legend(fontsize=13)
plt.tight_layout()
plt.savefig('../figures/fig1A/frequency_len_chi_all_kb_concatenated.svg')
plt.savefig('../figures/fig1A/frequency_len_chi_all_kb_concatenated.png')
