import pandas as pd
import os
from Bio import SeqIO
import pickle
from collections import Counter
import time
import numpy as np



def find_segments(r):
    segments = []
    start = r[0]
    for i in range(1, len(r)):
        if r[i] != r[i - 1] + 1:
            segments.append([start, r[i - 1]])
            start = r[i]
    segments.append([start, r[-1]])
    return segments


def preprocess_gff(input_file):
    df = pd.read_csv(input_file, header=None, sep='\t', skiprows=9)
    #df.drop(2050332,inplace=True)
    df_i=df[df[2] == 'mRNA']
    df_i = df_i[df_i[8].str.contains('Name=NM_', case=False)]
    df_i=df_i.reset_index(drop=True)
    return df_i



def coding_regions(input_file):
    df_original = pd.read_csv(input_file, header=None, sep='\t', skiprows=9)
    df_original.drop(2197887,inplace=True)
    df_coding=df_original[df_original[2] == 'CDS']
    df_coding = df_coding[df_coding[8].str.contains('Name=NP_', case=False)]
    df_coding[3]=df_coding[3].astype(int)
    df_coding[4]=df_coding[4].astype(int)
    df_coding=df_coding.reset_index(drop=True)
    return df_coding


def each_gene(df_ori):
    gene={}
    for i in range(len(df_ori)):
        gene_name=df_ori[8][i].split(';gene=')[1]
        gene_name=gene_name.split(';')[0]
        try:
            gene_id=df_ori[8][i].split('ID=rna-')[1]
            gene_id=gene_id.split(';')[0]
        except:
            continue
        try:
            # print(gene_id)
            if df_ori[4]-df_ori[3]>gene[gene_name][1]:
                gene[gene_name]=[gene_id,df_ori[0][i],df_ori[3][i],df_ori[4][i]]
            else:
                continue
        except:
            gene[gene_name]=[gene_id,df_ori[0][i],df_ori[3][i],df_ori[4][i]]
    return gene



df_gene=preprocess_gff('../data/GCF_000001405.40_GRCh38.p14_genomic.gff')
gene=each_gene(df_gene)
df_NP=coding_regions('../data/GCF_000001405.40_GRCh38.p14_genomic.gff')


cds_range={}    
for i in range(len(df_NP)):
    gene_name=df_NP[8][i].split(';gene=')[1]
    gene_name=str(gene_name.split(';')[0]).upper()
    if gene_name in gene:
        if df_NP[3][i] >= gene[gene_name][2] and df_NP[4][i] <= gene[gene_name][3] and df_NP[0][i] == gene[gene_name][1]:
            try:
                cds_range[gene_name]=list(range(df_NP[3][i],df_NP[4][i]+1))+cds_range[gene_name]
            except:
                cds_range[gene_name]=list(range(df_NP[3][i],df_NP[4][i]+1))


coding_range={}
noncoding_range={}


for i in range(len(df_gene)):
    name=df_gene[8][i].split(';gene=')[1]
    name=name.split(';')[0]
    if name in gene:
        gene[name].append('')
        gene[name].append('')
        id=df_gene[8][i].split('ID=rna-')[1]
        id=id.split(';')[0]
        if gene[name][0] == id:
            gene[name][4]=df_gene[6][i]
            gene[name][5]=df_gene[8][i]


fjob1=open('../data/coding_longest_human_reset.gff','w')
fjob2=open('../data/noncoding_longest_human_reset.gff','w')

for key in cds_range:
    #if key.upper() in repeats:
    if True:
        main=find_segments(sorted(list(set(cds_range[key]))))
        coding_range[key]=main
        whole_range=list(range(int(gene[key][2]),int(gene[key][3])+1))
        noncoding_range[key]=list(set(whole_range)-set(cds_range[key]))
        if len(noncoding_range[key])!=0:
            noncoding_range[key]=find_segments(sorted(noncoding_range[key]))
            try:
                for h in range(len(noncoding_range[key])):
                    fjob2.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (gene[key][1],'noncoding',key,int(noncoding_range[key][h][0]),int(noncoding_range[key][h][1]),'.',gene[key][4],'.',gene[key][5]))
            except:
                print(key)
        if len(coding_range[key])!=0:
            for h in range(len(coding_range[key])):
                fjob1.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (gene[key][1],'coding',key,int(coding_range[key][h][0]),int(coding_range[key][h][1]),'.',gene[key][4],'.',gene[key][5]))



fjob1.close()
fjob2.close()

os.system('bedtools getfasta -fi ../data/GCF_000001405.40_GRCh38.p14_genomic.fna -fo ../data/human_coding_reset.fa -bed ../data/coding_longest_human_reset.gff -name -s')
os.system('bedtools getfasta -fi ../data/GCF_000001405.40_GRCh38.p14_genomic.fna -fo ../data/human_noncoding_reset.fa -bed ../data/noncoding_longest_human_reset.gff -name -s')

noncds=SeqIO.parse('../data/human_noncoding_reset.fa', "fasta")

cds=SeqIO.parse('../data/GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna', "fasta")
noncds_dict={}
cds_dict={}
for record in cds:
    id_raw=record.description
    if 'NP_' in id_raw:
        id1=id_raw.split('gene=')[1]
        id=id1.split(']')[0]
        if id in cds_dict and len(cds_dict[id])<len(str(record.seq)):
            cds_dict[id]=str(record.seq)
        elif id not in cds_dict:
            cds_dict[id]=str(record.seq)
            


with open('/scratch/sigbio_project_root/sigbio_project19/ruoxiwan/codon_repeat_legacy/hg38/cds_dict_human_reset.pkl', 'wb') as file:
    pickle.dump(cds_dict,file)

for record in noncds:
    id=record.id.split('::')[0]
    try:
        noncds_dict[id]+=str(record.seq)
    except:
        noncds_dict[id]=str(record.seq)

with open('/scratch/sigbio_project_root/sigbio_project19/ruoxiwan/codon_repeat_legacy/hg38/noncoding_dict_human_reset.pkl', 'wb') as file:
    pickle.dump(noncds_dict,file)