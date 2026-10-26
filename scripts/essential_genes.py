import pandas as pd
import pickle
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from statistics import mean, median



# first to compare the frequency of repeats in essential and non-essential genes' coding region - noncoding region

essential_gene_df=pd.read_csv('/scratch/sigbio_project_root/sigbio_project19/ruoxiwan/codon_repeat_legacy/essential_gene/deg_annotation_e_human.csv',sep=';',header=None)
essential_gene_list=list(set(essential_gene_df[2].tolist()))


with open('../data/fs_ratio.pkl','rb') as f:
    fs_ratio=pickle.load(f)

with open('../data/coding_repeat_reset_2.pkl','rb') as file:
    coding_repeat=pickle.load(file)

with open('../data/non_coding_repeat_reset_2.pkl','rb') as file:
    non_coding_repeat=pickle.load(file)

    
with open('../data/cds_dict_human_reset.pkl', 'rb') as file:
    cds_dict = pickle.load(file)

with open('../data/noncoding_dict_human_reset.pkl', 'rb') as file:
    non_cds_dict = pickle.load(file)



def add_freq_5(df_gene,threshold=5):
    df_gene[f'coding_{threshold}%']=0
    df_gene[f'noncoding_{threshold}%']=0
    for key in df_gene.index:
        count=0
        try:
            for motif in coding_repeat[key].keys():
                if fs_ratio[motif[:3]]>threshold:
                    count+=coding_repeat[key][motif]
            df_gene.loc[key,f'coding_{threshold}%']=count*1000/len(cds_dict[key])
        except:
            pass
        
    for key in df_gene.index:
        count=0
        try:
            for motif in non_coding_repeat[key].keys():
                if fs_ratio[motif[:3]]>threshold:
                    count+=non_coding_repeat[key][motif]
            df_gene.loc[key,f'noncoding_{threshold}%']=count*1000/len(non_cds_dict[key])
        except:
            pass
    df_gene[f'minus_{threshold}%']=df_gene[f'coding_{threshold}%']-df_gene[f'noncoding_{threshold}%']
    df_gene.replace([np.inf, -np.inf], np.nan,inplace=True)
    df_gene.dropna(inplace=True)
    return df_gene



# df_essential=pd.DataFrame.from_list(essential_gene_list,orient='index')
for threshold in [0.5, 1, 5]:
    L_essential=[key for key in coding_repeat.keys() if key in essential_gene_list]
    df_essential=pd.DataFrame(L_essential)
    df_essential.set_index(0,inplace=True)
    df_5=add_freq_5(df_essential,threshold)

    L_nonessential=[key for key in coding_repeat.keys() if key not in essential_gene_list]
    df_nonessential=pd.DataFrame(L_nonessential)
    df_nonessential.set_index(0,inplace=True)
    df_5_non=add_freq_5(df_nonessential,threshold)

    mwu=stats.mannwhitneyu(list(df_5[f'minus_{threshold}%']),list(df_5_non[f'minus_{threshold}%']),alternative='less')
    print(f'threshold={threshold}: mean essential={mean(df_5[f"minus_{threshold}%"])}, mean nonessential={mean(df_5_non[f"minus_{threshold}%"])}, p={mwu[1]}')

    data = [list(df_5[f'minus_{threshold}%']), list(df_5_non[f'minus_{threshold}%'])]
    labels = ['Essential', 'Non-essential']

    # boxplot
    fig, ax = plt.subplots(figsize=(6, 6))
    bp = ax.boxplot(
        data,
        tick_labels=labels,
        patch_artist=True,
        notch=True,
        widths=0.5,
        showmeans=True,
        meanprops=dict(marker='D', markerfacecolor='white', markeredgecolor='black', markersize=6),
        medianprops=dict(color='black', linewidth=2),
        boxprops=dict(edgecolor='black', linewidth=1.2),
        whiskerprops=dict(color='black', linewidth=1.2),
        capprops=dict(color='black', linewidth=1.2),
        flierprops=dict(marker='o', markerfacecolor='none', markeredgecolor='gray',
                         markersize=4, alpha=0.4)
    )

    for patch in bp['boxes']:
        patch.set_facecolor('white')

    ax.axhline(y=0, color='black', linewidth=1)
    ax.set_ylabel('Codon frequency in coding regions minus that in noncoding regions', fontsize=13)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    ax.set_ylim(-3.5, 3.9)
    ax.grid(axis='y', alpha=0.3)

    y_max, h = 3.3, 0.1
    ax.plot([1, 1, 2, 2], [y_max, y_max + h, y_max + h, y_max], color='black', linewidth=1.2)

    p_val = mwu[1]
    ax.text(1.5, y_max + h + 0.03, rf'$P = {p_val:.1e}$', ha='center', va='bottom', fontsize=12)

    plt.tight_layout()
    plt.savefig(f'../figures/Fig3/boxplot_essential_nonessential_{threshold}%_3.svg')
    plt.savefig(f'../figures/Fig3/boxplot_essential_nonessential_{threshold}%_3.png', dpi=300)
    plt.close(fig)

    # violinplot
    combined = np.concatenate(data)
    lo, hi = np.min(combined), np.max(combined)
    y_low, y_high = -3.5, 3.9
    if threshold == 5:
        y_low, y_high = -1.0, 1.0
    data_trim = [[v for v in d if y_low <= v <= y_high] for d in data]

    fig, ax = plt.subplots(figsize=(6, 6))
    vp = ax.violinplot(data_trim, showextrema=False, widths=0.7, points=1000)
    for body in vp['bodies']:
        body.set_facecolor('lightgray')
        body.set_edgecolor('dimgray')
        body.set_linewidth(0.6)
        body.set_alpha(1)
    # narrow box inside each violin for median/IQR; grey circles still show fliers outside 1.5 IQR
    ax.boxplot(data, widths=0.08, whis=1.5, patch_artist=True, showfliers=True, showmeans=True,
               meanprops=dict(marker='D', markerfacecolor='white', markeredgecolor='black', markersize=5),
               medianprops=dict(color='black', linewidth=2),
               boxprops=dict(facecolor='white', edgecolor='black', linewidth=1.2),
               whiskerprops=dict(color='black', linewidth=1.2),
               capprops=dict(color='black', linewidth=1.2),
               flierprops=dict(marker='o', markerfacecolor='none', markeredgecolor='gray',
                               markersize=4, alpha=0.4))
    ax.set_xticks([1, 2], labels)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)

    ax.axhline(y=0, color='black', linewidth=1)
    ax.set_ylabel('Codon frequency in coding regions minus that in noncoding regions', fontsize=13)
    ax.set_ylim(y_low, y_high)
    ax.grid(axis='y', alpha=0.3)

    y_top = y_high - 0.15 * (y_high - y_low)
    vy, vh = y_top, 0.08 * (y_high - y_low)
    ax.plot([1, 1, 2, 2], [vy, vy + vh, vy + vh, vy], color='black', linewidth=1.2)
    ax.text(1.5, vy + vh + 0.02 * (y_high - y_low), rf'$P = {p_val:.1e}$', ha='center', va='bottom', fontsize=12)

    plt.tight_layout()
    plt.savefig(f'../figures/Fig3/violin_essential_nonessential_{threshold}%_3.svg')
    plt.savefig(f'../figures/Fig3/violin_essential_nonessential_{threshold}%_3.png', dpi=300)
    plt.close(fig)

