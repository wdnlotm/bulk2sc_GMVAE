# Processing data for `scTAPE` training
The Python script provided reads raw count data from the `_1_data_process_for_GMAVE` directory and 
1. reformats those into dense matrix text files meeting all requirements for `scTAPE`,
2. take five samples from the test data, 
3. and make five TPM normalized pseudo-bulk RNA-seq data with their cell type proportions.
```
python data_prep_for_scTAPE_training.py
```
