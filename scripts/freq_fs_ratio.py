import pandas as pd
import pickle
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


def sci(x, digits=1):
    # mathtext "3.2 \times 10^{-5}" instead of "3.2e-05"
    mant, exp = f'{x:.{digits}e}'.split('e')
    return rf'{mant} \times 10^{{{int(exp)}}}'


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


with open('../data/fs_ratio.pkl','rb') as f:
    fs_ratio=pickle.load(f)


df=pd.DataFrame.from_dict(fs_ratio,orient='index',columns=['fs_ratio'])
df_2=pd.read_csv('../data/frequency_len_kb_reset_2.tsv',sep='\t',index_col=0)
#df_2=pd.read_csv('../data/frequency_len_kb.tsv',sep='\t',index_col=0)
df_2['minus']=df_2['coding']-df_2['noncoding']
for codon in df.index:
    if codontable[codon]=='*':
        continue
    else:
        count=df_2.loc[codon*3,'minus']+df_2.loc[codon*4,'minus']+df_2.loc[codon*5,'minus']+df_2.loc[codon*6,'minus']
        #count=df_2.loc[codon*3,'coding']+df_2.loc[codon*4,'coding']+df_2.loc[codon*5,'coding']+df_2.loc[codon*6,'coding']
        df.loc[codon,'minus']=count
        #df.loc[codon,'coding']=count
    #df.loc[codon,'minus']=df_2.loc[codon*5,'minus']

scc=stats.spearmanr(df['fs_ratio']/100,df['minus'], alternative='less')

df_no_out=df[df['minus']>=-0.1]
print('excluded:',df[df['minus']<-0.1].index.tolist())
scc_no_out=stats.spearmanr(df_no_out['fs_ratio']/100,df_no_out['minus'], alternative='less')
pcc_no_out=stats.pearsonr(df_no_out['fs_ratio']/100,df_no_out['minus'], alternative='less')

print('all:',scc)
print('without outliers (n=%d): spearman'%len(df_no_out),scc_no_out,'pearson',pcc_no_out)
#pcc=stats.pearsonr(df['fs_ratio'],df['minus'])
#print(scc,pcc)

plt.figure(figsize=(8,6))
plt.scatter(np.log10(df['fs_ratio']/100),df['minus'],color='black',s=9,alpha=0.7)
plt.xlabel(r'$\log_{10}$(frameshifting ratio)',fontsize=15)
plt.ylabel('Codon repeat frequency in coding regions\nminus that in noncoding regions',fontsize=15)
plt.tick_params(labelsize=13)
#plt.text(10, 0.031,f'pcc={pcc[0]:.2f}, p_val={pcc[1]:.3f}', fontsize=10)
plt.text(-3.5, -0.2, f'ρ = {scc[0]:.2f}, $P$ = ${sci(scc[1], 2)}$', fontsize=15)
plt.tight_layout()
plt.savefig('../figures/fs_ratio/fs_ratio_freq_all_repeats_kb_reset_5.svg')
plt.savefig('../figures/fs_ratio/fs_ratio_freq_all_repeats_kb_reset_5.png')