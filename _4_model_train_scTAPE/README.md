## scTAPE 
`scTAPE` is used to predict cell type proportions in bulk RNA-seq. However, many hyperparameters are inaccessible through the default interface of the `scTAPE` package. To address this, the source code of `scTAPE` was modified to expose additional hyperparameters, and a new package, `scTAPE_for_bulk2scGMVAE-0.1-py3-none-any.whl`, was created based on the modified code. The package can be installed by running the following command in `_4_model_train_scTAPE` folder. 
```
pip install dist/scTAPE_for_bulk2scGMVAE-0.1-py3-none-any.whl
```
Inputs for `scTAPE` training are
```
import scTAPE_for_bulk2scGMVAE
from scTAPE_for_bulk2scGMVAE import Deconvolution

Sigm, Pred = Deconvolution(scRNA_train_data, bulk_data, sep='\t',scaler='mms',
                                   datatype='TPM', genelenfile='Genelength.txt', 
                                   mode='overall', adaptive=False, variance_threshold=0.85,
                                   save_model_name = "_lung_scRNA", batch_size=600, epochs=640,
                                  samplenumber=15000, sampling_num=1200, seed=15640)



```
