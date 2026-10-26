import pandas as pd 
from collections import defaultdict
import pyBigWig
import pickle
import numpy as np
from Bio import SeqIO
import re
import sys
from scipy import stats
import matplotlib.pyplot as plt
from statistics import mean, median
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



def find_codon_repeats(sequence, name):
    # [start, end)
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


def repeat_sites(repeats_by_motif):
    return sorted((motif[:3], start, end)
                  for motif, spans in repeats_by_motif.items()
                  for start, end in spans)


def map_detected_site(sites, motif, position):
    """Keep the existing nearest-reported-position mapping on the longest NP CDS."""
    motif = motif.upper()
    codon = motif[:3]
    if len(motif) < 9 or len(motif) % 3 or motif != codon * (len(motif) // 3):
        return None
    candidates = [site for site in sites
                  if site[0] == codon and site[2] - site[1] >= len(motif)]
    return min(candidates, key=lambda site: (abs(site[1] - position), site[1]), default=None)


def nearest_control(sites, detected_site, detected_sites):
    # Distance is between repeat starts in CDS nucleotides; upstream wins ties.
    candidates = [site for site in sites
                  if site[0] == detected_site[0] and site not in detected_sites]
    return min(candidates,
               key=lambda site: (abs(site[1] - detected_site[1]), site[1]), default=None)



def find_position(sequence,motif):
    seq_len = len(sequence)
    if seq_len < 9:
        return {}
    positions = [m.start() for m in re.finditer(f"(?={motif})", sequence)]
    return positions


def contiguous_blocks(positions):
    # a repeat's CDS span can cross an exon/intron boundary -- ranges[gene] jumps by more
    # than +-1 there, so split into genomically-contiguous sub-blocks instead of bridging
    # the gap with one interval (which would pull intron bases into the PhastCons query).
    blocks = []
    block_start = 0
    for i in range(1, len(positions)):
        if positions[i] - positions[i - 1] not in (1, -1):
            blocks.append((block_start, i))
            block_start = i
    blocks.append((block_start, len(positions)))
    return blocks


def phastcons_score(gene, start_idx, end_idx):
    # start_idx, end_idx are CDS-relative nucleotide indices (end_idx exclusive, Python-slice
    # style). ranges[gene] runs in genomic-ascending order for plus-strand genes but
    # genomic-descending order for minus-strand genes (see parse_location), so each
    # contiguous block's two endpoints must be sorted before building a valid (ascending)
    # bigWig query interval -- otherwise a minus-strand gene produces a reversed interval
    # and pyBigWig raises "Invalid interval bounds!". Blocks are combined as a
    # length-weighted mean so an exon-spanning repeat doesn't count intron bases.
    pos = ranges[gene][start_idx:end_idx]
    chrom = chr_dict[chr[gene]]
    weighted_sum, total_len = 0.0, 0
    for b_start, b_end in contiguous_blocks(pos):
        g_start, g_end = pos[b_start], pos[b_end - 1]
        lo, hi = (g_start, g_end + 1) if g_start <= g_end else (g_end, g_start + 1)
        block_score = phastcons.stats(chrom, lo, hi, type="mean")[0]
        block_len = b_end - b_start
        weighted_sum += block_score * block_len
        total_len += block_len
    return weighted_sum / total_len




df=pd.read_csv('../data/sup_tab_5.csv')
df.sort_values(by=[df.columns[0], df.columns[2]], inplace=True)
df.reset_index(drop=True, inplace=True)

for i in range(len(df)):
    df.iloc[i,0] = df.iloc[i,0].split('_')[0]

df.drop(df['type'].index[df['type'] == 'trans-frame(splice site)'], inplace=True)
df.drop(df['type'].index[df['type'] == 'annotated protein'], inplace=True)
df.drop(df['type'].index[df['type'] == 'small peptide'], inplace=True)
df.reset_index(drop=True, inplace=True)

to_drop = []
count=0
for i in range(1,len(df)):
    if df.iloc[i,0] == df.iloc[i-1,0]:
        if df.iloc[i,2] == df.iloc[i-1,2] and df.iloc[i,3] in range(int(df.iloc[i-1,3])-len(df.iloc[i-1,2]),int(df.iloc[i-1,3])+len(df.iloc[i-1,2])):
            to_drop.append(i)
        else:
            print(df.iloc[i-1])
            print(df.iloc[i])
            count+=1

# are those two different repeats or same repeat different frameshifting events
df.drop(to_drop, inplace=True)
df.reset_index(drop=True, inplace=True)


# Load phastCons scores from WIG file
# phastcons_scores = defaultdict(dict)

TRACK = sys.argv[1] if len(sys.argv) > 1 else 'phastcons'  # 'phastcons' (default) or 'phylop'
TRACK_NAME, BW_URL = {
    'phastcons': ('PhastCons', "https://hgdownload.cse.ucsc.edu/goldenPath/hg38/phastCons20way/hg38.phastCons20way.bw"),
    'phylop': ('PhyloP', "https://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP20way/hg38.phyloP20way.bw"),
}[TRACK]
phastcons = pyBigWig.open(BW_URL)


HUMAN_FASTA = '../data/GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna'
GENE_RE = re.compile(r'\[gene=([^\]]+)\]')
PID_RE = re.compile(r'\[protein_id=([^\]]+)\]')
LOC_RE = re.compile(r'\[location=([^\]]+)\]')


def parse_location(location_str):
    # returns 0-based genomic positions in the same base-order as record.seq
    # (5'->3' of the CDS), or None if the location string can't be parsed
    complement = location_str.startswith('complement(')
    inner = location_str[len('complement('):-1] if complement else location_str
    inner = inner[len('join('):-1] if inner.startswith('join(') else inner
    exons = []
    for part in inner.split(','):
        m = re.match(r'(\d+)\.\.(\d+)', part.strip())
        if not m:
            return None
        exons.append((int(m.group(1)), int(m.group(2))))  # 1-based inclusive, as in the header
    positions = []
    if complement:
        for start, end in reversed(exons):
            positions.extend(range(end - 1, start - 2, -1))
    else:
        for start, end in exons:
            positions.extend(range(start - 1, end))
    return positions


def parse_human_longest_np(path):
    longest = {}
    for record in SeqIO.parse(path, 'fasta'):
        desc = record.description
        pid_m = PID_RE.search(desc)
        if not pid_m or not pid_m.group(1).startswith('NP_'):
            continue
        gene_m = GENE_RE.search(desc)
        loc_m = LOC_RE.search(desc)
        if not gene_m or not loc_m:
            continue
        gene = gene_m.group(1)
        seq = str(record.seq).upper()
        nc_accession = record.id.split('|')[1].split('_cds_')[0]
        if gene not in longest or len(seq) > len(longest[gene][0]):
            longest[gene] = (seq, loc_m.group(1), nc_accession)
    return longest


chr_dict={'NC_000001.11':'chr1','NC_000002.12':'chr2','NC_000003.12':'chr3','NC_000004.12':'chr4','NC_000005.10':'chr5','NC_000006.12':'chr6','NC_000007.14':'chr7','NC_000008.11':'chr8','NC_000009.12':'chr9','NC_000010.11':'chr10','NC_000011.10':'chr11','NC_000012.12':'chr12','NC_000013.11':'chr13','NC_000014.9':'chr14','NC_000015.10':'chr15','NC_000016.10':'chr16','NC_000017.11':'chr17','NC_000018.10':'chr18','NC_000019.10':'chr19','NC_000020.11':'chr20','NC_000021.9':'chr21','NC_000022.11':'chr22','NC_000023.11':'chrX','NC_000024.10':'chrY','NC_012920.1':'chrM'}

human_longest = parse_human_longest_np(HUMAN_FASTA)
print(f'human genes (NP only, longest per gene): {len(human_longest)}')

sequence={}
ranges={}
repeats={}
chr={}
skipped_location=[]
skipped_scaffold=[]
for gene, (seq, location_str, nc_accession) in human_longest.items():
    if nc_accession not in chr_dict:
        # unplaced/unlocalized scaffold (NW_...) -- PhastCons has no track for these
        skipped_scaffold.append(gene)
        continue
    pos = parse_location(location_str)
    if pos is None or len(pos) != len(seq):
        skipped_location.append(gene)
        continue
    sequence[gene]=seq
    ranges[gene]=pos
    chr[gene]=nc_accession

print(f'genes with usable location: {len(sequence)}, skipped (bad location): {len(skipped_location)}, '
      f'skipped (unplaced scaffold): {len(skipped_scaffold)}')

for id in sequence.keys():
    repeats[id]=find_codon_repeats(sequence[id],id)


df['phastcons'] = np.nan
not_found = []
multiple = []
multiple_detected = []
sites_by_gene = {gene: repeat_sites(reps) for gene, reps in repeats.items()}
detected_by_gene = defaultdict(set)
mapped_rows = {}

# Map ALL detections before selecting controls, including sites that fail scoring.
for i, row in df.iterrows():
    gene = row['gene name']
    if gene not in sites_by_gene:
        not_found.append(gene)
        continue
    site = map_detected_site(sites_by_gene[gene], row['codon repeat'],
                             row['position of codon in peptide(bp)'])
    if site is not None:
        mapped_rows[i] = (gene, site)
        detected_by_gene[gene].add(site)

pair_records = []
paired_sites = set()
for i, (gene, site) in mapped_rows.items():
    codon, start, end = site
    try:
        det = phastcons_score(gene, start, end)
        if not np.isfinite(det):
            continue
        df.loc[i, 'phastcons'] = det
    except Exception as exc:
        print(gene, 'detected scoring error:', exc)
        continue
    if (gene, site) in paired_sites:
        continue  # multiple detection rows at one repeat contribute only one pair
    control = nearest_control(sites_by_gene[gene], site, detected_by_gene[gene])
    if control is None:
        continue
    _, control_start, control_end = control
    try:
        ctrl = phastcons_score(gene, control_start, control_end)
        if not np.isfinite(ctrl):
            continue
    except Exception as exc:
        print(gene, 'control scoring error:', exc)
        continue
    paired_sites.add((gene, site))
    multiple_detected.append(det)
    multiple.append(ctrl)
    pair_records.append(dict(gene=gene, codon=codon,
                             detected_start=start, detected_end=end,
                             control_start=control_start, control_end=control_end,
                             distance_nt=abs(control_start - start),
                             detected_score=det, control_score=ctrl,
                             difference=det - ctrl))


df.dropna(inplace=True)
detected=list(df['phastcons'])

not_detected=[]
for gene in repeats.keys():
    if gene in list(df['gene name']):
        for motif in repeats[gene].keys():
            for i in range(len(repeats[gene][motif])):
                try:
                    score=float(phastcons_score(gene, repeats[gene][motif][i][0], repeats[gene][motif][i][1]))
                    not_detected.append(score)
                except:
                    if chr_dict[chr[gene]]!='chrM':
                        print(gene, 'error',chr_dict[chr[gene]])



print(mean(multiple), mean(multiple_detected))
print(f'{TRACK_NAME}: {len(multiple_detected)} pairs, mean frameshift={mean(multiple_detected):.4f}, '
      f'mean control={mean(multiple):.4f}, median frameshift={median(multiple_detected):.4f}, '
      f'median control={median(multiple):.4f}')


wilcoxon=stats.wilcoxon(multiple_detected,multiple,alternative='less')  # expect frameshift sites less conserved
print(f'Wilcoxon (frameshift < control) p={wilcoxon.pvalue:.3g}')

data = [multiple_detected, multiple]
labels = ['Frameshift sites', 'Control sites']

p_val = wilcoxon[1]
fig, ax = plt.subplots(figsize=(6, 6))

parts = ax.violinplot(data, showmeans=False, showmedians=False, showextrema=False, widths=0.8)
for pc in parts['bodies']:
    pc.set_facecolor('0.93')
    pc.set_edgecolor('0.5')
    pc.set_linewidth(1)
    pc.set_alpha(1)

bp = ax.boxplot(
    data, positions=[1, 2], widths=0.15, patch_artist=True,
    showmeans=True,
    meanprops=dict(marker='D', markerfacecolor='white', markeredgecolor='black', markersize=6),
    medianprops=dict(color='black', linewidth=2.2),
    boxprops=dict(facecolor='white', edgecolor='black', linewidth=1.6),
    whiskerprops=dict(color='black', linewidth=1.6),
    capprops=dict(color='black', linewidth=1.6),
    flierprops=dict(marker='o', markerfacecolor='none', markeredgecolor='gray', markersize=3, alpha=0.4)
)

ax.set_xticks([1, 2])
ax.set_xticklabels(labels)
ax.tick_params(labelsize=14)
ax.set_ylabel(f'{TRACK_NAME} score', fontsize=16)
ax.grid(axis='y', alpha=0.3)
if TRACK == 'phastcons':  # bounded 0-1
    ax.set_ylim(-0.05, 1.22)          # extra headroom for bracket + label
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])   # keep gridlines only at real data ticks
    y_bar = 1.05
    h = 0.04
else:  # PhyloP is unbounded (negative = faster than neutral), so place the bracket from the data
    lo_y, hi_y = min(map(min, data)), max(map(max, data))
    span = hi_y - lo_y
    y_bar, h = hi_y + 0.05 * span, 0.04 * span
    ax.set_ylim(lo_y - 0.05 * span, y_bar + h + 0.12 * span)

# significance bracket, now safely below the top of the frame
ax.plot([1, 1, 2, 2], [y_bar, y_bar + h, y_bar + h, y_bar], color='black', linewidth=1.2)
ax.text(1.5, y_bar + h + 0.01, rf'$P = {p_val:.1e}$', ha='center', va='bottom', fontsize=14)

plt.tight_layout()
plt.savefig(f'../figures/PhyloP/violinplot_{TRACK}_optionB_5.png', dpi=300)
plt.savefig(f'../figures/PhyloP/violinplot_{TRACK}_optionB_5.svg')
