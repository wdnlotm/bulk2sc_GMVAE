## Generate using saved models and the test data.
```
from scTAPE_predictor_by_modelpt import *
from GMVAE_generate_by_modelpt import *
```
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
