# Train GMVAE
The main source code for `bulk2sc` is located in the `_model_source_codes` directory. Executing the provided Python script trains the GMVAE model using the data from the `_1_data_process_for_GMVAE` directory for 4000 epochs. The model is saved every 1000 epochs. During training, saved models are stored in the `models` directory, while the TensorBoard event file is saved in the `runs` directory.
```python
epoch=4001; savefreq=1000; h_dim1=128; h_dim2=64; z_dim=32

traindir=../_1_data_process_for_GMVAE/train/
testdir=../_1_data_process_for_GMVAE/test/

python main_zinb.py --epochs ${epoch} --h_dim1 ${h_dim1} --h_dim2 ${h_dim2} --z_dim ${z_dim} \
--savefreq ${savefreq} -dl ${traindir} -tdl ${testdir} --dataset lung
```

