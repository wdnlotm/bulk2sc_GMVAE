# Train the GMVAE
```bash
datadir=../_II_data_process_4_GMVAE/train/
testdir=../_II_data_process_4_GMVAE/test/

python main_zinb.py --epochs 4001 --h_dim1 128 --h_dim2 64 --z_dim 32 -dl ${datadir} -tdl ${testdir} --dataset lung > output.txt
```
