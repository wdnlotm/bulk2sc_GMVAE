# Processing data for `scTAPE` training
The Python script provided reads raw count data from the `_0_data_download` directory and prepares a train-test split for GMVAE training. The data is split into 75% for training and 25% for testing. The resulting datasets are saved into separate folders named `train` and `test`, respectively.
```
python process_train_test_data.py
```
