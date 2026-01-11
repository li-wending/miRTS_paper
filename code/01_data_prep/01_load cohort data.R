#' Raw sequencing data can be downloaded from SRA as described in the manuscript,
#' here include the following:
#' 1) TMM-normalized data used as input for the miR-TS,
#' 2) and the output miR-TS scores.

# TMM-norm data:
SY_miRNA_TMMcpm <- read.csv("./data/processed_miRNA_expression/SY_miRNA_TMMcpm.csv")
DFTJ_miRNA_TMMcpm <- read.csv("./data/processed_miRNA_expression/DFTJ_miRNA_TMMcpm.csv")
NAS_miRNA_TMMcpm <- read.csv("./data/processed_miRNA_expression/NAS_miRNA_TMMcpm.csv")

# miR-TS scores:
SY_miRTS_scores.257 <- read.csv("./data/miR_TS_scores/SY_miRTS_scores.257.csv")
DFTJ_miRTS_scores.257 <- read.csv("./data/miR_TS_scores/DFTJ_miRTS_scores.257.csv")
NAS_miRTS_scores.257 <- read.csv("./data/miR_TS_scores/NAS_miRTS_scores.257.csv")
