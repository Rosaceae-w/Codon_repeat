import pickle
import matplotlib.pyplot as plt
import re
from scipy import stats
import math
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from scipy.stats import pearsonr



def find_motif_with_regex(sequence, motif):
    #pattern = f'(?=(?P<codon>{motif}))'
    matches = re.finditer(motif, sequence)
    return [match.start() for match in matches if match.start() % 3 != 0]



with open('../data/cds_dict_human_reset.pkl', 'rb') as file:
    cds_dict=pickle.load(file)



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


# same as conservation.py: {motif: [(start_nt, end_nt), ...]}, [start, end)
def find_codon_repeats(sequence):
    repeats = {}
    sequence = sequence.upper()
    codons = [sequence[i:i + 3] for i in range(0, len(sequence) - 2, 3)]
    start = 0
    while start < len(codons):
        end = start + 1
        while end < len(codons) and codons[end] == codons[start]:
            end += 1
        if end - start >= 3 and codons[start] in codontable:
            motif = codons[start] * (end - start)
            repeats.setdefault(motif, []).append((start * 3, end * 3))
        start = end
    return repeats


# coding_repeat_reset_2.pkl only stores counts per motif, so get positions from the sequences
coding_repeat_site = {key: find_codon_repeats(cds_dict[key]) for key in cds_dict}
'''
coding_stop_sites={}
for key in cds_dict:
    sequence=cds_dict[key]
    coding_stop_sites[key]=[]
    for codon in codontable.keys():
        if codontable[codon]=='*':
            sites=find_motif_with_regex(sequence, codon)
            coding_stop_sites[key].extend(sites)



with open('../data/coding_stop_sites.pkl', 'wb') as file:
    pickle.dump(coding_stop_sites, file)
'''
with open('../data/coding_stop_sites.pkl', 'rb') as file:
    coding_stop_sites=pickle.load(file)

coding_stop_distance={}
for key in cds_dict:
    if key not in coding_stop_sites:
        pass
    else:
        coding_stop_distance[key]={}
        for motif in coding_repeat_site[key]:
            if coding_repeat_site[key][motif]==[]:
                pass
            else:
                coding_stop_distance[key][motif]={}
                for i in range(len(coding_repeat_site[key][motif])):
                    for j in range(len(coding_stop_sites[key])):
                        distance=coding_stop_sites[key][j]-coding_repeat_site[key][motif][i][1]
                        if distance>=0:
                            try:
                                if coding_stop_distance[key][motif][i]>distance:
                                    coding_stop_distance[key][motif][i]=distance
                            except KeyError:
                                coding_stop_distance[key][motif][i]=distance



all_distances = []
for key in coding_stop_distance:
    for motif in coding_stop_distance[key]:
        all_distances.extend(coding_stop_distance[key][motif].values())

mean_distance = np.mean(all_distances)
median_distance = np.median(all_distances)
print(f'mean distance to nearest out-of-frame stop codon: {mean_distance}')
print(f'median distance to nearest out-of-frame stop codon: {median_distance}')

plt.figure(figsize=(8, 6))
plt.hist(all_distances, bins=np.linspace(0, 800, 51), color='skyblue', edgecolor='black', alpha=0.7)
x_pad = 800 * 0.05
plt.xlim(-x_pad, 800 + x_pad)
plt.xlabel('Distance to nearest out-of-frame stop codon (nt)', fontsize=15)
plt.ylabel('Count', fontsize=15)
plt.tick_params(labelsize=13)
plt.tight_layout()
plt.savefig('../figures/stop_codon/coding_stop_distance_hist.png')
plt.savefig('../figures/stop_codon/coding_stop_distance_hist.svg')

