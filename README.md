# test-scanpy-workflow
Testing dockerized script that preprocesses 10X_genomics outputs from cell ranger and generates plots

### Inputs

The input to the python script requires a `.csv` file with the following format:

```
Sample,Filepath
MantonBM1,demodata/cellranger_output/MantonBM1_HiSeq_1/raw_feature_bc_matrix.h5
MantonBM2,demodata/cellranger_output/MantonBM2_HiSeq_1/raw_feature_bc_matrix.h5
MantonBM3,demodata/cellranger_output/MantonBM3_HiSeq_1/raw_feature_bc_matrix.h5
MantonBM4,demodata/cellranger_output/MantonBM4_HiSeq_1/raw_feature_bc_matrix.h5
MantonBM5,demodata/cellranger_output/MantonBM5_HiSeq_1/raw_feature_bc_matrix.h5
MantonBM6,demodata/cellranger_output/MantonBM6_HiSeq_1/raw_feature_bc_matrix.h5
MantonBM7,demodata/cellranger_output/MantonBM7_HiSeq_1/raw_feature_bc_matrix.h5
MantonBM8,demodata/cellranger_output/MantonBM8_HiSeq_1/raw_feature_bc_matrix.h5

```

The filepath should reflect the relative or absolute location of the file.

### Running the Script

The python script can be run:

```
python scandemo.py sample_sheet.csv

```

### Expected Outputs

For each channel:

*  A counts matrix in `.hda5` format
*  A UMAP `.png` file colored by 'leiden' with clusters labeled by a number and the gene symbol for the most variable gene observed for the cluster
*  A `.png1 file showing the leiden rank gene group plots for each cluster 


Example results can be viewed in the github directories /figures and /results


