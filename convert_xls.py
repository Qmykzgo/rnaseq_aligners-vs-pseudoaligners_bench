# convert_xls.py
import pandas as pd
import re
import os

xls_path = 'GSE83402_qPCR.xls'
gtf_path = '/mnt/d/rnaseq_aligners-vs-pseudoaligners_bench/data/reference/GRCh38_ERCC92.gtf'
out_path = '/mnt/d/rnaseq_aligners-vs-pseudoaligners_bench/data/truth_tables/GSE83402_qPCR_normalized.txt'

try:
    print("Parsing GTF for Gene ID to Symbol mapping...")
    gene_map = {}
    gene_id_pat = re.compile(r'gene_id "([^"]+)"')
    gene_name_pat = re.compile(r'gene_name "([^"]+)"')
    
    if os.path.exists(gtf_path):
        with open(gtf_path, "r") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) > 8 and parts[2] in ("gene", "transcript"):
                    attrs = parts[8]
                    gid_match = gene_id_pat.search(attrs)
                    gname_match = gene_name_pat.search(attrs)
                    if gid_match:
                        gid = gid_match.group(1)
                        gname = gname_match.group(1) if gname_match else gid
                        gene_map[gid] = gname
        print(f"Mapped {len(gene_map)} gene IDs from GTF.")
    else:
        print(f"Warning: GTF not found at {gtf_path}")
        
    print("Reading Excel...")
    df = pd.read_excel(xls_path, sheet_name='MAQC_Cq')
    df.columns = ['GeneID', 'MAQCA', 'MAQCB']
    
    print("Mapping Gene IDs to Symbols...")
    if gene_map:
        df['GeneSymbol'] = df['GeneID'].map(lambda x: gene_map.get(x, x).upper())
    else:
        df['GeneSymbol'] = df['GeneID']
        
    # Reorder columns
    df_out = df[['GeneSymbol', 'MAQCA', 'MAQCB']]
    
    # Save as tab-separated
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df_out.to_csv(out_path, sep="\t", index=False)
    print("Successfully saved normalized qPCR truth table to:", out_path)
    
except Exception as e:
    print(f"Error during conversion: {e}")
