## scTAPE 
`scTAPE` is used to predict cell type proportions in bulk RNA-seq. However, many hyperparameters are inaccessible through the default interface of the `scTAPE` package. To address this, the source code of `scTAPE` was modified to expose additional hyperparameters, and a new package, `scTAPE_for_bulk2scGMVAE-0.1-py3-none-any.whl`, was created based on the modified code. The package can be installed by running the following command in `_4_model_train_scTAPE` folder. 
```
pip install dist/scTAPE_for_bulk2scGMVAE-0.1-py3-none-any.whl
```
Inputs for `scTAPE` training are
```
import scTAPE_for_bulk2scGMVAE
from scTAPE_for_bulk2scGMVAE import Deconvolution
import pandas as pd
import os

Sigm, Pred = Deconvolution(data_dir+train_file, data_dir+bulk_sample_file, sep='\t',scaler='mms',
                                   datatype='TPM', genelenfile=data_dir+gene_len_file, 
                                   mode='overall', adaptive=False, variance_threshold=variance_threshold_input,
                                   save_model_name = model_name, batch_size=batch_size_input, epochs=epoch_input,
                                  samplenumber=samplenumber_input, sampling_num=sampling_num_input, 
                                   seed=epoch_input+samplenumber_input)



```
