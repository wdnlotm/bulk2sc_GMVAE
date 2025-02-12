# Train the GMVAE
The main source code for `bulk2sc` is located in the `_model_source_codes` directory. Executing the provided Python script will train the GMVAE model using the data in the `_1_data_process_for_GMVAE` directory for 4000 epochs. The model is saved every 1000 epochs in a designated directory.
```bash
epoch=4001; savefreq=1000; h_dim1=128; h_dim2=64; z_dim=32

traindir=../_1_data_process_for_GMVAE/train/
testdir=../_1_data_process_for_GMVAE/test/

python main_zinb.py --epochs ${epoch} --h_dim1 ${h_dim1} --h_dim2 ${h_dim2} --z_dim ${z_dim} \
--savefreq ${savefreq} -dl ${traindir} -tdl ${testdir} --dataset lung
```
Training will save models in the `models` directory and the Tensorboard event file will be saved in the `runs` directory.
