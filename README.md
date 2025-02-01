# bulk2sc_GMVAE
### Fig.1 Model schematic diagram
<img src="fig/Figure1_scRNAGMVAE.png" width="80%" alt="Image description">

The `bulk2sc` model is a machine learning tool based on a Gaussian mixed variational autoencoder (GMVAE) that generates single-cell RNA sequencing (scRNA-seq) data from bulk RNA sequencing (bulk RNA-seq) data. This repository demonstrates the operation of `bulk2sc` using the scRNA-seq dataset [GSE130148](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130148).

The `bulk2sc` model is made up of three main components (Fig. 1):
* I. scRNA Gaussian Mixture Variational Autoencoder (GMVAE) – Learns patterns in scRNA-seq data across different cell types and generates synthetic scRNA-seq data.
* II. Bulk RNA-seq Deconvolution – Predicts cell type proportions from bulk RNA-seq data using the training from scTAPE-seq data. The [`scTAPE`](https://sctape.readthedocs.io/) is adopted.
* III. scRNA-seq Generator – Integrates components I and II to generate scRNA-seq data from bulk RNA-seq input.

The GMVAE model (Part I in Fig. 1) requires a scRNA-seq raw count table in Matrix Market exchange format (MTX file) along with the corresponding cluster (cell type) data in a CSV file. The [`scTAPE`](https://sctape.readthedocs.io/)(Part II in Fig. 1) utilizes the scRNA-seq raw count table as a dense matrix, with the cluster data as the row index. Additionally, TPM-normalized bulk data is needed to accompany the single-cell data.

### Fig.2 Procedures and operation order
<img src="fig/Figure2_operation_procedure_dependency.png" width="55%" alt="Procedures and dependency">
