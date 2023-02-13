#!/usr/bin/env python
# coding: utf-8


import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from warnings import simplefilter
simplefilter("ignore", UserWarning)
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
simplefilter("ignore", FutureWarning)

sc.settings.verbosity = 1
sc.settings.figdir = 'figures'

sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=200,
    facecolor = 'white', figsize=(6,6), format='png')



import os

def make_dirs(file):
    try:
       os.makedirs(file)
    except FileExistsError:
       # directory already exists
       pass
    return

make_dirs("figures")
make_dirs("figures/umapfigures")
make_dirs("results")
make_dirs("figures/rank_genes_groups_leidenfigures")

import sys
first_arg = sys.argv[1]

print(f"This worksheet was selected as input: {first_arg}")

# importing module
from pandas import *
 

#data
def read_csv_input(file):
    if file.endswith('.csv') :
       print (f'Reading in {file}')
    else:
         raise ValueError(f'Unsupported filetype: {file}')
    return

read_csv_input(first_arg)

# reading CSV file
data = read_csv(first_arg)


samplenames = data['Sample']
print("Sample Names")
print(samplenames)



#filelist = [h5 + s + '/raw_feature_bc_matrix.h5' for s in samplenames]
#len(filelist)

filelist=data['Filepath']

print("File Paths")
print(filelist)



def read_h5(file):
    if file.endswith('.h5') :
       print (f'Reading in {file}')
    else:
        raise ValueError(f'Unsupported filetype: {file}')
    return 
              
for i in range (len(filelist)):
    read_h5(filelist[i])


print("Please ignore the next warning lines about var not being unique")
print("This is a scanpy artifact.")



adatas = [] 
for i in range(len(samplenames)):
    label = samplenames[i]
    filename = filelist[i]
    adata = sc.read_10x_h5(filename)
    adata.var_names_make_unique()
    adata.obs['sample'] = label
    adatas.append(adata) 

adata = adatas[0].concatenate(adatas[1:])

print(f"Anndata matrix for all batches contains {adata.shape[0]} obs and {adata.shape[1]} var")


# ## Print the Number of Observations Per Sample


print ("Number of observations per sample")
print(adata.obs['sample'].value_counts())

# ## Filter Cells and Genes
# 
# * Remove Cells that have Have Gene Expression levels below 500 and above 6000
# * Remove Genes that are in less thatn 5% of the Cells


#Initial quality control
min_genes = 200
min_cells = 3

sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)

print("Initial QC")
print(f"Filtering cells that express less than {min_genes} and genes expressed in less than {min_cells} cells")
print(f"Anndata matrix for all batches contains {adata.shape[0]} obs and {adata.shape[1]} var")



#calculate qc metrics, add qc metric for whether gene is mitochondrial
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


# ## Filter Cells with Between a Minimum and Maximum Number of Genes Expressed


#Filter cells with
#number of genes expressed >500
#number of genes expressed < 6000

min_genes= 500
max_genes = 6000

sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_cells(adata, max_genes=max_genes)

print (f"After filtering for genes expressed greater than {min_genes} and less than {max_genes}") 
#print(adata.obs['sample'].value_counts())
print(f"Anndata matrix for all batches contains {adata.shape[0]} obs and {adata.shape[1]} var")


# ## Filter Cells with High Mitochondrial Gene Expression


max_mito=10

mito_genes = adata.var_names.str.startswith('MT-')
# for each cell compute fraction of counts in mito genes vs. all genes
# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
adata.obs['percent_mito'] = np.sum(
    adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1).A1

adata = adata[(adata.obs['pct_counts_mt'] < max_mito), :]

#print(adata.obs['sample'].value_counts())

print (f"After filtering for cells with more than {max_mito}% mitochondrial gene expression")
print(f"Anndata matrix for all batches contains {adata.shape[0]} obs and {adata.shape[1]} var")


# ## Filter out genes expressed in less than 5% of the cells

#The number of genes that are found in less than 5% of the cells are removed
min_percent=0.05

adata = adata[:, adata.var['n_cells_by_counts'] > adata.var['n_cells_by_counts'].max()*min_percent]

#print(adata.obs['sample'].value_counts())

print (f"After filtering for genes expressed in less than {min_percent*100}% of cells")
print(f"Anndata matrix for all batches contains {adata.shape[0]} obs and {adata.shape[1]} var")


# ## Filter out Technical Artifacts


malat1 = adata.var_names.str.startswith('MALAT1')
# we need to redefine the mito_genes since they were first 
# calculated on the full object before removing low expressed genes.
mito_genes = adata.var_names.str.startswith('MT-')
remove = np.add(mito_genes, malat1)
keep = np.invert(remove)

adata = adata[:,keep]
print (f"After filtering for technical artifacts")
print(f"Anndata matrix for all batches contains {adata.shape[0]} obs and {adata.shape[1]} var")


#sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt','pct_counts_ribo', 'pct_counts_hb'],
  #           jitter=0.4, groupby = 'sample', rotation = 45)



#sc.pl.highest_expr_genes(adata, n_top=20)


# ## This section is optional 
# 
# This section looks at chrY markers and determines the likelihood of the gender of the sample origin.
# 


annot = sc.queries.biomart_annotations(
        "hsapiens",
        ["ensembl_gene_id", "external_gene_name", "start_position", "end_position", "chromosome_name"],
    ).set_index("external_gene_name")
#>>> adata.var[annot.columns] = annot

chrY_genes = adata.var_names.intersection(annot.index[annot.chromosome_name == "Y"])
#chrY_genes


adata.obs['percent_chrY'] = np.sum(
    adata[:, chrY_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1 * 100


# color inputs must be from either .obs or .var, so add in XIST expression to obs.
adata.obs["XIST-counts"] = adata.X[:,adata.var_names.str.match('XIST')].toarray()

#sc.pl.scatter(adata, x='XIST-counts', y='percent_chrY', color="sample")

for i in range(len(filelist)):

    MeanChrY = adata[adata.obs['batch'] == str(i) , :].obs['percent_chrY'].mean()
    print(adata.obs['sample'][0])
    print(f"Calculated Percentage of ChrY Genes {MeanChrY:.2E}")
    if MeanChrY > 0.05:
        sex = 'M'        
    else:
        sex = 'F'
        
    print(f"Assigned sex: {sex}")
    
    

    batchname = adata[adata.obs['batch'] == str(i) , :].obs_names
    adata.obs.loc[batchname, 'sex'] = sex


#sc.pl.violin(adata, ["XIST-counts", "percent_chrY"], jitter=0.4, groupby = 'sample', rotation= 45,show=False)

#plt.savefig("demodata/gender_counts_violin.pdf")




# ## Extract Highly Expressed Genes

# ### Log Normalize Data by 100,000 Reads Per Cell


sc.pp.filter_genes(adata, min_cells=1)
sc.pp.normalize_total(adata, target_sum=1e5)
sc.pp.log1p(adata)

print("Log Normalize Data by 100,000 Reads Per Cell")
#print(adata.n_obs, adata.n_vars)


# ## Extract Highly Variable Genes

sc.pp.highly_variable_genes(adata, n_top_genes=2000, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

print("Extract Highly Variable Genes")
print(adata.var.highly_variable.head())



#sc.pl.highly_variable_genes(adata)



print("Regress Highly Variable Genes")

sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])




max_value=10

print(f"Scale Matrix at {max_value}")

sc.pp.scale(adata, max_value=max_value)


# ## Run PCA



n_pcs=50

print(f"Run PCA with number of pcs {n_pcs}")

sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True, show=False, n_pcs=n_pcs, save='') # scanpy generates the filename automatically


# ## Run Nearest Neighbor Analysis


print(f"Run Nearest Neighbor with number of pcs {n_pcs}")

sc.pp.neighbors(adata, n_neighbors=100, n_pcs=n_pcs)


# ## Run UMAP Data


print("Generating UMAP")
sc.tl.umap(adata)


# ## Run Leiden Analysis


resolution = 1.3

print(f"Run Leiden Analysis with resolution {resolution}")

sc.tl.leiden(adata,resolution=resolution)



leidennum = len(set(adata.obs.leiden))
print(f"Leiden analysis generated {leidennum} clusters")


# ## Test Color by Leiden on UMAP


#sc.pl.umap(adata, color=['leiden'],legend_loc='on data')


# ## Find marker genes for each leiden cluster using Mann-Whitney-U test

print("Outputting rank genes group figures to figures/rank_genes_groups_leidenfigures")

sc.settings.figdir = 'figures'

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', key_added = "wilcoxon")

for i in range(len(samplenames)):
    
    samplename = list(set(adata[adata.obs['batch'] == str(i), :].obs['sample']))[0]


    sc.pl.rank_genes_groups(adata[adata.obs['batch'] == str(i) , :], n_genes=25, sharey=False, show=False, key="wilcoxon",
                        save=(str(sc.settings.figdir) + '/' + samplename + '_rank_genes_groups_leiden.png'))

#sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', key_added = "wilcoxon")
#sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, key="wilcoxon",save="")
print ("The error  invalid value encountered in log2 relates to the lack of fold change data for this example")
print ("It can be ignored in this demo.")


#The suggested function from the provided tutorial does not work as expected, 
#so a dataframe was constructed of the top genes using the following text#

top_genes=[]
tg_rows = []
for i in range(leidennum):

#    print(f"Top 5 genes in Cluster {str(i)}")
#    print(sc.get.rank_genes_groups_df(adata, group=str(i), key='wilcoxon').names[0:5])
    top_genes.append(sc.get.rank_genes_groups_df(adata, group=str(i), key='wilcoxon').names[0:10])
    tg_rows.append(str(i))

df = pd.DataFrame(top_genes)
df['cluster'] = tg_rows
df = df.set_index('cluster')
df.T

#Cluster across the top, top 10 genes listed in descending order by score



df.T.to_csv('results/top_genes_by_cluster.csv',index=False)



#Rename cluster based on top gene marker

new_cluster_names=df.index.astype(str)+'_'+df[0]
#new_cluster_names = df[0] 
new_cluster_names
adata.rename_categories('leiden', new_cluster_names)
#sc.pl.umap(adata, color=['leiden'],legend_loc='on data')

#sc.pl.umap(adata, color=['leiden'],legend_loc='on data',
#show=True, palette=sns.color_palette("tab10", n_colors=25),
#    legend_fontsize=8, frameon=True, title='Leiden Colored UMAP for All Samples')



#Violin plot
#ax = sc.pl.stacked_violin(adata, df[0], groupby='leiden', rotation=90)



results_file = 'results/combined_MantonBM.h5ad'
adata.write(results_file)



print("Outputting count matrixfiles for each sample to results/")

for i in range(len(samplenames)):
    
    samplename = list(set(adata[adata.obs['batch'] == str(i), :].obs['sample']))[0]
    results_file = 'results/' + samplename + 'count_matrix.h5ad'
    adata[adata.obs['batch'] == str(i), :].write(results_file)
    



new_cluster_names=df.index.astype(str)+'_'+ df[0]
#new_cluster_names = df[0] 
new_cluster_names
adata.rename_categories('leiden', new_cluster_names)



sc.pl.umap(adata[adata.obs['batch'] == str(i) , :], color=['leiden'],legend_loc='on data',
show=True, palette=sns.color_palette("tab10", n_colors=25),
    legend_fontsize=8, frameon=True, title='Leiden Map for Sample ' + samplename  )


print("Outputting sample figures to figures/umapfigurs")

sc.settings.figdir = 'figures'

for i in range(len(samplenames)):
    
    samplename = list(set(adata[adata.obs['batch'] == str(i), :].obs['sample']))[0]
    sc.pl.umap(adata[adata.obs['batch'] == str(i) , :], color=['leiden'],legend_loc='on data',
show=False, palette=sns.color_palette("tab10", n_colors=25),
    legend_fontsize=8, frameon=True, title='Leiden Map for Sample ' + samplename,
    save=(str(sc.settings.figdir) + '/' + samplename + '_Leiden_UMAP.png'))
    
    
