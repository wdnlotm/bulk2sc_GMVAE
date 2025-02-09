# Processing data for `scTAPE` training
The Python script provided reads raw count data from the `_1_data_process_for_GMAVE` directory and performs the following tasks:
1. Reformats the data into dense matrix text files that meet all requirements for `scTAPE`.
2. Takes five samples from the test data.
3. Creates five TPM-normalized pseudo-bulk RNA-seq datasets with their respective cell type proportions.
```
python data_prep_for_scTAPE_training.py
```
