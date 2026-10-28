import pickle
import matplotlib.pyplot as plt
import re
from scipy import stats
import math
import pandas as pd
from collections import Counter
import numpy as np



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
    'TAC': 'Y', 'TAT': 'Y', 
    'TGC': 'C', 'TGT': 'C', 'TGG': 'W'
}

signal_df=pd.read_csv('../data/human_signal_peptide_proteins.tsv',sep='\t')
signal_peptide_genes=set(signal_df['gene_name'].tolist())


def find_codon_repeats_bins(sequence, name, threshold):
    repeats = {}
    seq_len = len(sequence)
    if seq_len < 9:
        return {}
    # Convert to NumPy array
    seq_array = np.array(list(sequence.upper()))
    # Extract triplets aligned to reading frames (only every 3rd position)
    num_codons = seq_len // 3
    triplets = np.lib.stride_tricks.sliding_window_view(seq_array, window_shape=3)[::3]
    codons = [''.join(triplet) for triplet in triplets]
    categorized_repeats = {codon: {3: 0, 4: 0, 5: 0, 6: 0} for codon in fs_ratio if fs_ratio[codon] > threshold}
    # -------------------------
    # Bin setup: 25 bins of 4% each
    # Each bin covers 4% of num_codons
    # -------------------------
    num_bins = 25  # 100% / 4% = 25 bins
    bin_size_codons = num_codons / num_bins  # bin size in codon units (can be float)
    # bin_repeat_counts[i] = total number of repeat EVENTS starting in bin i
    bin_repeat_counts = np.zeros(num_bins, dtype=float)
    # bin_length_kb[i] = length of bin i in kb (nucleotides / 1000)
    bin_lengths_kb = np.array([
        (min((b + 1) * bin_size_codons, num_codons) - b * bin_size_codons) * 3 / 1000
        for b in range(num_bins)
    ])
    prev_codon = codons[0]
    consecutive_count = 1
    def record_repeat(start_codon_idx, end_codon_idx_exclusive, codon):
        """
        Records a repeat event:
        - Updates categorized_repeats
        - Updates repeats dict with relative position
        - Updates bin_repeat_counts
        """
        repeat_len = end_codon_idx_exclusive - start_codon_idx
        rep = ''.join(codons[start_codon_idx:end_codon_idx_exclusive])
        # Categorize repeat length
        if repeat_len <= 6:
            categorized_repeats[codon][repeat_len] += 1
        else:
            categorized_repeats[codon][6] += 1
        # Relative position (0 to 1)
        rel_pos = start_codon_idx * 3 / seq_len
        if rep not in repeats:
            repeats[rep] = [rel_pos]
        else:
            repeats[rep].append(rel_pos)
        # Assign repeat event to bin based on start codon index
        bin_idx = min(int(start_codon_idx / bin_size_codons), num_bins - 1)
        bin_repeat_counts[bin_idx] += 1
    for i in range(1, len(codons)):
        if codons[i] == prev_codon:
            consecutive_count += 1
            # Handle repeat at end of sequence
            if i == len(codons) - 1:
                if consecutive_count >= 3 and fs_ratio.get(prev_codon, 0) > threshold:
                    record_repeat(i - consecutive_count + 1, i + 1, prev_codon)
        else:
            if fs_ratio.get(prev_codon, 0) > threshold:
                if consecutive_count >= 3:
                    record_repeat(i - consecutive_count, i, prev_codon)
            prev_codon, consecutive_count = codons[i], 1
    # -------------------------
    # Compute frequency per kb for each bin
    # bin_freq_per_kb[i] = repeat events in bin i / length of bin i in kb
    # -------------------------
    bin_freq_per_kb = np.divide(
        bin_repeat_counts,
        bin_lengths_kb,
        out=np.zeros_like(bin_repeat_counts),
        where=bin_lengths_kb > 0
    )
    # bin_labels: midpoint % of each bin (2%, 6%, 10%, ... 98%)
    bin_labels = [f"{int((b + 0.5) * 4)}%" for b in range(num_bins)]
    bins_result = {
        "bin_labels": bin_labels,           # ['2%', '6%', ..., '98%']
        "bin_repeat_counts": bin_repeat_counts.tolist(),
        "bin_lengths_kb": bin_lengths_kb.tolist(),
        "bin_freq_per_kb": bin_freq_per_kb.tolist()
    }
    return categorized_repeats, repeats, bins_result



    
with open('../data/cds_dict_human_reset.pkl', 'rb') as file:
    cds_dict = pickle.load(file)

with open('../data/noncoding_dict_human_reset.pkl', 'rb') as file:
    non_cds_dict = pickle.load(file)



with open('../data/fs_ratio.pkl','rb') as f:
    fs_ratio=pickle.load(f)


num_bins = 25

def aggregate_bins(bins_dict):
    """
    Averages bin_freq_per_kb across all genes.
    Returns mean and SEM per bin.
    """
    all_freqs = np.zeros((len(bins_dict), num_bins))
    for idx, (gene, bins) in enumerate(bins_dict.items()):
        all_freqs[idx] = bins["bin_freq_per_kb"]
    
    mean_freq = np.mean(all_freqs, axis=0)
    sem_freq = np.std(all_freqs, axis=0) / np.sqrt(len(bins_dict))
    return mean_freq, sem_freq


def aggregate_paired_diff_bins(coding_dict, noncoding_dict):
    """
    Paired per-gene difference: for genes present in both coding_dict and
    noncoding_dict (already excluding signal-peptide genes), computes
    coding - noncoding bin_freq_per_kb per gene, then averages across genes.
    Returns mean and SEM of the per-gene difference per bin.
    """
    shared_genes = sorted(coding_dict.keys() & noncoding_dict.keys())
    all_diffs = np.zeros((len(shared_genes), num_bins))
    for idx, gene in enumerate(shared_genes):
        coding_freq = np.array(coding_dict[gene]["bin_freq_per_kb"])
        noncoding_freq = np.array(noncoding_dict[gene]["bin_freq_per_kb"])
        all_diffs[idx] = coding_freq - noncoding_freq

    mean_diff = np.mean(all_diffs, axis=0)
    sem_diff = np.std(all_diffs, axis=0) / np.sqrt(len(shared_genes))
    return mean_diff, sem_diff



for threshold in [0.5, 1, 5]:
    print('coding')
    coding_repeat = {}
    coding_bins = {}  # Store bins_result for each gene
    for gene in cds_dict.keys():
        if gene not in signal_peptide_genes:
            sequence = cds_dict[gene]
            try:
                categorized_repeats, repeats, bins_result = find_codon_repeats_bins(sequence, gene, threshold)
                coding_repeat[gene] = repeats
                coding_bins[gene] = bins_result
            except:
                print(gene)
                continue

    print('noncoding')
    non_coding_repeat = {}
    non_coding_bins = {}  # Store bins_result for each gene
    for gene in non_cds_dict.keys():
        if gene not in signal_peptide_genes:
            sequence = non_cds_dict[gene]
            try:
                categorized_repeats, repeats, bins_result = find_codon_repeats_bins(sequence, gene, threshold)
                non_coding_repeat[gene] = repeats
                non_coding_bins[gene] = bins_result
            except:
                print(gene)
                continue


    mean_coding, sem_coding = aggregate_bins(coding_bins)
    mean_noncoding, sem_noncoding = aggregate_bins(non_coding_bins)

    bin_labels = [f"{(b + 0.5) * 0.04:.2f}" for b in range(num_bins)]
    x = np.arange(num_bins)

    rho_coding, p_coding = stats.spearmanr(x, mean_coding)
    rho_noncoding, p_noncoding = stats.spearmanr(x, mean_noncoding)


    diff_freq, sem_diff = aggregate_paired_diff_bins(coding_bins, non_coding_bins)



    fig2a, ax2a = plt.subplots(figsize=(8, 6))
    ax2a.errorbar(x, diff_freq, yerr=sem_diff, fmt="none", ecolor="black", elinewidth=1, capsize=3, zorder=1)
    ax2a.scatter(x, diff_freq, color="mediumpurple", alpha=0.8, edgecolor="black", zorder=2)

    ax2a.axhline(0, color="gray", linewidth=0.8, linestyle="--")
    ax2a.set_xticks(x)
    ax2a.set_xticklabels(bin_labels, rotation=45, ha="right")
    ax2a.tick_params(labelsize=13)
    ax2a.set_xlabel("Relative position in coding sequences", fontsize=15)
    ax2a.set_ylabel("Codon repeat frequency in coding regions\nminus that in noncoding regions", fontsize=15)

    plt.tight_layout()
    plt.savefig(f'../figures/repeat_site/freq_diff_scatter_per_bin_vs_position_{threshold}p_codon_excluding_signal_peptides.svg')
    plt.savefig(f'../figures/repeat_site/freq_diff_scatter_per_bin_vs_position_{threshold}p_codon_excluding_signal_peptides.png', dpi=300)
    plt.close(fig2a)
