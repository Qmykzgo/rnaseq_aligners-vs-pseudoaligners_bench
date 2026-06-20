# subset_to_chr21.py
import os
import re

ref_dir = '/home/yer_kanat/Downloads/rnqbench/data/reference'
raw_genome = os.path.join(ref_dir, 'Homo_sapiens.GRCh38.dna.primary_assembly.fa')
raw_gtf = os.path.join(ref_dir, 'Homo_sapiens.GRCh38.110.gtf')
ercc_fa = os.path.join(ref_dir, 'ERCC92/ERCC92.fa')
ercc_gtf = os.path.join(ref_dir, 'ERCC92/ERCC92.gtf')

out_genome = os.path.join(ref_dir, 'GRCh38_ERCC92.fa')
out_gtf = os.path.join(ref_dir, 'GRCh38_ERCC92.gtf')
out_gtf_no21 = os.path.join(ref_dir, 'GRCh38_ERCC92_noChr21.gtf')
out_tx = os.path.join(ref_dir, 'GRCh38_ERCC92_transcriptome.fa')
out_tx2gene = os.path.join(ref_dir, 'tx2gene.tsv')
out_junc = os.path.join(ref_dir, 'chr21_true_junctions.tsv')

print("Subsettting genome FASTA to Chromosome 21...")
chr21_seq = []
keep = False
with open(raw_genome, 'r') as f:
    for line in f:
        if line.startswith('>'):
            # Match chr 21
            if line.strip().split()[0] == '>21':
                keep = True
                chr21_seq.append(line)
            else:
                keep = False
        elif keep:
            chr21_seq.append(line)

print("Appending ERCC sequences to genome...")
with open(ercc_fa, 'r') as f:
    ercc_seq = f.readlines()

with open(out_genome, 'w') as f:
    f.writelines(chr21_seq)
    f.writelines(ercc_seq)
print(f"Saved subsetted genome to {out_genome}")

print("Subsettting GTF to Chromosome 21...")
chr21_gtf = []
with open(raw_gtf, 'r') as f:
    for line in f:
        if line.startswith('#'):
            chr21_gtf.append(line)
        else:
            parts = line.split('\t')
            if len(parts) > 0 and parts[0] == '21':
                chr21_gtf.append(line)

print("Appending ERCC annotations to GTF...")
with open(ercc_gtf, 'r') as f:
    ercc_gtf_lines = f.readlines()

with open(out_gtf, 'w') as f:
    f.writelines(chr21_gtf)
    f.writelines(ercc_gtf_lines)
print(f"Saved subsetted GTF to {out_gtf}")

print("Creating chr21-masked GTF (ERCC only)...")
with open(out_gtf_no21, 'w') as f:
    # Comments
    f.writelines([l for l in chr21_gtf if l.startswith('#')])
    # ERCC only
    f.writelines(ercc_gtf_lines)
print(f"Saved masked GTF to {out_gtf_no21}")

# Delete old compiled files so they get regenerated
for p in [out_tx, out_tx2gene, out_junc]:
    if os.path.exists(p):
        os.remove(p)
        print(f"Removed old file to force rebuild: {p}")

# Also clean up old indices
import shutil
for d in ['star_idx', 'star_idx_noChr21', 'hisat2_idx', 'hisat2_idx_noChr21', 'salmon_idx', 'kallisto_idx.idx']:
    path = os.path.join(ref_dir, d)
    if os.path.exists(path):
        if os.path.isdir(path):
            shutil.rmtree(path)
        else:
            os.remove(path)
        print(f"Removed old index to force rebuild: {path}")

print("Subsetting completed successfully!")
