# bulk2sc_GMVAE
### Fig.1 Model schematic diagram
<img src="fig/Figure1_scRNAGMVAE.png" width="80%" alt="Image description">


bulk2sc is a machine learning model based on a Gaussian mixed variational autoencoder that generates scRNA-seq data from bulk RNA-seq data. 
This repository aims to demonstrate bulk2sc operation using scRNA-seq data [GSE130148](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130148).
GMVAE model uses scRNA-seq raw count table in Matrix Market exchange format (mtx file) and corresponding cluster (cell type) data (csv file). scTAPE uses scRNA-seq raw count table as a dense matrix with the cluster data as the row index. TPM-normalized bulk data needs to accompany the single-cell data.

The bulk2sc model is a machine learning tool based on a Gaussian mixed variational autoencoder (GMVAE) that generates single-cell RNA sequencing (scRNA-seq) data from bulk RNA sequencing (bulk RNA-seq) data. This repository demonstrates the operation of bulk2sc using the scRNA-seq dataset GSE130148.

The GMVAE model requires a scRNA-seq raw count table in Matrix Market exchange format (MTX file) along with the corresponding cluster (cell type) data in a CSV file. scTAPE utilizes the scRNA-seq raw count table as a dense matrix, with the cluster data as the row index. Additionally, TPM-normalized bulk data is needed to accompany the single-cell data.

<span style="color:blue">some *blue* text</span>.
<span style="background-color:gray">gray background</span>

### Fig.2 Procedures and operation order
<img src="fig/Figure2_operation_procedure_dependency.png" width="55%" alt="Procedures and dependency">
