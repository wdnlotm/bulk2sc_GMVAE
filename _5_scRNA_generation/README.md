# Generate using saved models and test data.
For comparison purposes, synthetic scRNA-seq data is generated from the scRNA-seq data in the test data set.
```python
from scTAPE_predictor_by_modelpt import *
modelpt="../_4_model_train_scTAPE/model_640_600_15000_1200_0.85.pt" #scTAPE model saved after scTAPE training
geneleng='../_2_data_process_for_scTAPE/GeneLength.txt'             #gene length for tpm normalization 
rna_file='../_2_data_process_for_scTAPE/lung_data_to_train_scTAPE_from_lung_test.txt' #scRNA-seq in the test data
                                                                                      #It is a raw count data.
intersect_list='../_4_model_train_scTAPE/genename_640_600_15000_1200_0.85.csv' # gene list that was used in scTAPE training.

pred_prop, total_count = predicted_proportion_by_TAPE_model(modelfile=modelpt, 
                                   genelenfile=geneleng, 
                                   intersect_gene=intersect_list, 
                                   scrna=scrna_file)
```
The last command in the above script will
1. create pseudo-bulk RNA-seq data from the scRNA-seq data, `rna_file`, in the test data.
2. perform TPM-normalization using the `genelenfile`.
3. trim genes using the gene list used in `scTAPE` training.
4. produce outputs: `pred_prop`, predicted proportion and `total_count`, the number of cells in the input scRNA-seq data.

Subsequently, execution of the scripts below will generate scRNA-seq data based on the cell counts specified in `cell_counts` using the `modelfile`. The results will be saved with the following filenames: `matrix_generated_using_test_data.mtx` and `cluster_generated_using_test_data.csv`.
```python
from GMVAE_generate_by_modelpt import *
cell_counts = [round(x*total_count) for x in pred_prop]  #cell counts in each cell types
modelfile="../_3_model_train_GMVAE/models/lung_GMVAE_zinb_ep4000_wholemodel.pt" #model saved after the GMVAE training
generate_from_pt_cellcounts(modelfile, cell_counts, suffix="generated_using_test_data")
```
## Generate from one bulk RNA-seq data.
This is demonstrated in `_7_quick_demo`
