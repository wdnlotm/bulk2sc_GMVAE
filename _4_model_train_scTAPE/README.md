## scTAPE 
`scTAPE` is used to predict cell type proportions in bulk RNA-seq. However, many hyperparameters are inaccessible through the default interface of the `scTAPE` package. To address this, the source code of `scTAPE` was modified to expose additional hyperparameters, and a new package, `scTAPE_for_bulk2scGMVAE-0.1-py3-none-any.whl`, was created based on the modified code. 
```
pip install dist/scTAPE_for_bulk2scGMVAE-0.1-py3-none-any.whl # run this in _4_model_train_scTAPE
```
