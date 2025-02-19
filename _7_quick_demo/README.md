# Quick demo
The script files and data files in this directory will demonstrate how trained models work together. 

The `data_for_quick_demo` directory has all the data and model files for the demonstration. `lung_GMVAE_wholePTmodel.pt` has the trained GMVAE model and `model_scTAPE_lung.pt` has the trained scTAPE model. They are trained on the lung train data. The `bulk_data_to_test.txt` is a pseudo bulk RNA-seq data aggregated from the lung test data. So, it is unseen data from the model's point of view, but it comes with the ground truth cell type proportions.  
