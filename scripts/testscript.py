#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt


sc.settings.verbosity = 0 
sc.settings.set_figure_params(dpi=80)

MantonBM1 = sc.read_10x_h5('/demodata/cellranger_output/MantonBM1_HiSeq_1/raw_feature_bc_matrix.h5')
MantonBM1.var_names_make_unique()
MantonBM2 = sc.read_10x_h5('/demodata/cellranger_output/MantonBM2_HiSeq_1/raw_feature_bc_matrix.h5')
MantonBM2.var_names_make_unique()
MantonBM3 = sc.read_10x_h5('/demodata/cellranger_output/MantonBM3_HiSeq_1/raw_feature_bc_matrix.h5')
MantonBM3.var_names_make_unique()
MantonBM4 = sc.read_10x_h5('/demodata/cellranger_output/MantonBM4_HiSeq_1/raw_feature_bc_matrix.h5')
MantonBM4.var_names_make_unique()
MantonBM5 = sc.read_10x_h5('/demodata/cellranger_output/MantonBM5_HiSeq_1/raw_feature_bc_matrix.h5')
MantonBM5.var_names_make_unique()
MantonBM6 = sc.read_10x_h5('/demodata/cellranger_output/MantonBM6_HiSeq_1/raw_feature_bc_matrix.h5')
MantonBM6.var_names_make_unique()

# merge into one object.
adata = MantonBM1.concatenate(MantonBM2, MantonBM3, MantonBM4, MantonBM5, MantonBM6)

# and delete individual datasets to save space
#del(MantonBM1, MantonBM2, MantonBM3)
#del(MantonBM4, MantonBM5, MantonBM6)


print(adata.obs['sample'].value_counts())

#######Checkpoint 1#########

# mitochondrial genes
adata.var['mt'] = adata.var_names.str.startswith('MT-') 
# ribosomal genes
adata.var['ribo'] = adata.var_names.str.startswith(("RPS","RPL"))
# hemoglobin genes.
adata.var['hb'] = adata.var_names.str.contains(("^HB[^(P)]"))

#adata.var


sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo','hb'], percent_top=None, log1p=False, inplace=True)


mito_genes = adata.var_names.str.startswith('MT-')
# for each cell compute fraction of counts in mito genes vs. all genes
# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
np.seterr(invalid='ignore')
# for each cell compute fraction of counts in mito genes vs. all genes
# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
adata.obs['percent_mt2'] = np.sum(
    adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1).A1


sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt','pct_counts_ribo', 'pct_counts_hb'],
             jitter=0.4, groupby = 'sample', rotation= 45,show=False)



someFig = plt.figure()
someFig.savefig('./mito_counts_violin.pdf')
plt.close(someFig) # --> (Note 1)
del someFig

######CHECKPOINT 2############

#Filtering Genes

sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_cells(adata, max_genes=6000)
sc.pp.filter_genes(adata, min_cells=3)

print(adata.n_obs, adata.n_vars)

