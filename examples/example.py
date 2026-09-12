import gzip
import tempfile
from pathlib import Path
from rds2py import read_rds
import numpy as np
import pandas as pd
import paocseq as pq

input_data_path="./CMV"
sample="CMV"
gex_path = [input_data_path]
vdj_path = [input_data_path + "/all_contig_annotations.csv"]

marker_gene = ["IFNG", "TNF"]

print("Step 1: CombineData (10x + VDJ import + QC + CD4/CD8 + clonotypes)")
cell_data = pq.combine_data(
    gex_path=gex_path,
    vdj_path=vdj_path,
    marker_gene=marker_gene,
    sample_name=[sample],
    n_feature_lower_q= 0.01,
    n_feature_upper_q= 0.99,
    percent_mt_upper_q = 0.99,
)
adata = cell_data[0]
print(f"  -> {adata.n_obs} cells x {adata.n_vars} genes")

print("Step 2: threshold_marker_genes")
pq.threshold_marker_genes(adata, marker_gene, threshold_cutoff=0.9)
print(f"  -> {(adata.obs['Threshold_IFNG'] == 'high').sum()} IFNG-high cells")

print("Step 3: annotate_cell_types (clonotype summary table)")
summary = pq.annotate_cell_types(cell_data, goi="IFNG")

print("Step 4: read in reference data - this was stored as a single cell experiment RDS file format")
path_to_your_RDS_file=""
r_object = read_rds(path_to_your_RDS_file)

reference_anndata=r_object.to_anndata()

marker_gene_set = #Add a list of genes here
reference = pq.make_reference(reference_anndata, gene_list=marker_gene_set, n_cells=200, random_state=0)
scored = pq.classify_cells(adata.copy(), reference, gene_list=marker_gene_set, distance=0)
print(scored.obs[["Mdist"]].describe())