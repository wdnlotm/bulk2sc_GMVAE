library(Matrix)
library(R.utils)

download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE130148&format=file&file=GSE130148%5Fraw%5Fcounts%2ERData%2Egz", "lung_GSE130148.RData.gz") 
gunzip("lung_GSE130148.RData.gz", remove=FALSE)
load("lung_GSE130148.RData")

writeMM(raw_counts, "matrix_raw_counts.mtx")
write.table(rownames(raw_counts), file = "genes_from_raw.tsv", row.names=FALSE, sep="\t", col.names=FALSE, quote=FALSE)
write.table(colnames(raw_counts), file = "barcodes_from_raw.tsv", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
