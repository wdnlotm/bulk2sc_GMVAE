# Processing data for `scTAPE` training
The Python script provided reads raw count data from the `_1_data_process_for_GMAVE` directory and 
1. reformats those into dense matrix text files meeting all requirements for `scTAPE`,
2. take five samples from the test data, 
3. prepares a train-test split for GMVAE training. The data is split into 75% for training and 25% for testing. The resulting datasets are saved into separate folders named `train` and `test`, respectively.
```
python data_prep_for_scTAPE_training.py
```
