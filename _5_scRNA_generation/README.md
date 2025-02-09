# Generate using saved models and the test data.
First, functions in the source codes are loaded.
```
from scTAPE_predictor_by_modelpt import *
from GMVAE_generate_by_modelpt import *
```
The scripts below will
1. make a pseudo-bulk RNA-seq data from scRNA-seq data, `scrna_file`, in the test data.
2. TPM-normalize using `genelenfile`.
3. trim genes using the gene list that is used in `scTAPE` training.
4. produce outputs, `pred_prop`: predicted proportion, `cell_count`: the number of cells in the input scRNA-seq data.
```
modelpt="../_4_model_train_scTAPE/model_640_600_15000_1200_0.85_ver3.pt"
geneleng='../_4_model_train_scTAPE/GeneLength.txt'

scrna_file='../_2_data_process_for_scTAPE/lung_data_to_train_scTAPE_from_lung_test.txt'
intersect_list='../_4_model_train_scTAPE/genename_640_600_15000_1200_0.85_ver3.csv'

pred_prop, cell_count = predicted_proportion_by_TAPE_model(modelfile=modelpt, 
                                   genelenfile=geneleng, 
                                   intersect_gene=intersect_list, 
                                   scrna=scrna_file, cellcounts=None)
```
```
cell_counts = [round(x) for x in pred_prop*cell_count]
print(cell_counts)
```
The result of printing `cell_counts` is [347, 194, 86, 156, 300, 78, 421, 153, 77, 77, 198, 49, 45, 151, 82, 54, 80].

Executing the scripts below will generate scRNA-seq data following the cell counts in `cell_counts` using the `modelfile`. 
The results will be saved as `matrix_generated_using_test_data.mtx` and `cluster_generated_using_test_data.csv` following the `suffix`.
```
modelfile="../_3_model_train_GMVAE/saved_model_for_quick_demo/lung_GMVAE_zinb_ep4000_wholemodel.pt"
generate_from_pt_cellcounts(modelfile, cell_counts, suffix="generated_using_test_data")
```
## Generate from one bulk RNA-seq data.
### Coming soon.
