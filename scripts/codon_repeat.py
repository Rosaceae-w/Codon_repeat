import pandas as pd
import pickle
import numpy as np
import re
import os
import time
from collections import Counter




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




def find_codon_repeats(sequence,name):
    repeats=[]
    seq_len = len(sequence)
    if seq_len < 9:
        return {}
    # Convert to NumPy array
    seq_array = np.array(list(sequence.upper()))
    # Extract triplets aligned to reading frames (only every 3rd position)
    num_codons = seq_len // 3
    triplets = np.lib.stride_tricks.sliding_window_view(seq_array, window_shape=3)[::3]
    codons = [''.join(triplet) for triplet in triplets]
    categorized_repeats = {codon:{3: 0, 4: 0, 5: 0, 6: 0} for codon in set(codontable.keys())}
    prev_codon = codons[0]
    consecutive_count = 1
    for i in range(1, len(codons)):
        if codons[i] == prev_codon:
            consecutive_count += 1
            if i == len(codons)-1:
                if consecutive_count >= 3 and consecutive_count <= 6:
                    categorized_repeats[prev_codon][consecutive_count]+=1
                    rep= ''.join(codon for codon in codons[i-consecutive_count+1:i+1])
                    repeats.append(rep)
                elif consecutive_count > 6:
                    categorized_repeats[prev_codon][6]+=1
                    rep= ''.join(codon for codon in codons[i-consecutive_count+1:i+1])
                    repeats.append(rep)
            #print(i, consecutive_count)
        elif prev_codon in codontable.keys():
            if consecutive_count >= 3 and consecutive_count <= 6:
                categorized_repeats[prev_codon][consecutive_count]+=1
                rep= ''.join(codon for codon in codons[i-consecutive_count:i])
                repeats.append(rep)
            elif consecutive_count > 6:
                #print(prev_aa, amino_acids[i-1],consecutive_count)
                categorized_repeats[prev_codon][6]+=1
                rep= ''.join(codon for codon in codons[i-consecutive_count:i])
                repeats.append(rep)
            prev_codon, consecutive_count = codons[i], 1
    return categorized_repeats,Counter(repeats)



with open('../data/cds_dict_human_reset.pkl', 'rb') as file:
    cds_dict = pickle.load(file)

with open('../data/noncoding_dict_human_reset.pkl', 'rb') as file:
    non_cds_dict = pickle.load(file)

time1 = time.time()
gene_codon_repeats_coding={}
for gene in cds_dict:
    gene_codon_repeats_coding[gene]=find_codon_repeats(cds_dict[gene],gene)[1]

time2 = time.time()
print('Time taken:', time2-time1)

time1 = time.time()
gene_codon_repeats_noncoding={}
for gene in non_cds_dict:
    try:
        gene_codon_repeats_noncoding[gene]=find_codon_repeats(non_cds_dict[gene],gene)[1]
    except:
        print(gene)
        continue

time2 = time.time()
print('Time taken:', time2-time1)

'''
with open('../data/coding_repeat_reset_2.pkl', 'wb') as file:
    pickle.dump(gene_codon_repeats_coding, file)

with open('../data/non_coding_repeat_reset_2.pkl', 'wb') as file:
    pickle.dump(gene_codon_repeats_noncoding, file)
'''