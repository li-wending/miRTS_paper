#' This script explains how to use the miRTS package to obtain the miR-TS scores for:
#' an example dataset same as in the miRTS package;
#' SY, DFTJ and NAS cohorts


# Example: Hepatitis C dataset----
# Example code from: https://github.com/li-wending/miRTS
# install.packages("remotes")
# remotes::install_github("li-wending/miRTS") # , build_vignettes = TRUE

# Alternatively, download the latest .tar.gz from the GitHub Releases page and install:
# install.packages("C:/path/to/miRTS_1.0.0.tar.gz", repos = NULL, type = "source")

library(miRTS)
# vignette(topic = "Intro_to_miRTS", package = "miRTS") # available only when vignettes was built
# ?miRTS_score
# CIBERSORT_download()
# if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
# BiocManager::install("edgeR")
# BiocManager::install("preprocessCore")
# install.packages('e1071')
# install.packages('parallel')
# miR_TS.output <- miRTS_score(Input_df = example_counts) # miRNAs on the rows (e.g., 'hsa-miR-1-3p') and samples on the columns

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("This vignette requires ggplot2. Please run 'install.packages('ggplot2')' to install it.")
}
if (!requireNamespace("scales", quietly = TRUE)) {
  stop("This vignette requires scales. Please run 'install.packages('scales')' to install it.")
}

data(example_counts)
data(miRTS_signature_v1)

# Deconvolute miRNA expression data to obtain miR-TS scores.
# CIBERSORT is essential for optimal miR-TS performance (see below):
#
source('C:/Users/wl2957/CIBERSORT.R') # path to your downloaded file

miR_TS.output <- miRTS_score(
  Input_df = example_counts,
  signature_matrix = miRTS_signature_v1,
  method = "cibersort"  #"cibersort", "xCell2", "MCP-counter"
)
# all.equal(miR_TS.output, miR_TS.output.bkup)

# Inspect returned objects####
# miR_TS.output is a list containing:
#  (1) a samples × tissues matrix of miR-TS scores,
#  (2) scaled mixture matrix (the `Input_df`),
#  (3) scaled signature matrix (the `signature_matrix`)
str(miR_TS.output)

plt_df <- stack(as.data.frame(miR_TS.output$proportions))
names(plt_df) <- c("value", "tissue")

library(ggplot2)
ggplot(plt_df) + aes(tissue, value) + geom_boxplot() +
  xlab("") + ylab("score") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 2),n.breaks = 5) +
  theme(axis.text.x = element_text(angle = 45,size = 12, hjust=1))

data(hepatitis_C.meta)

hepatitis_C.output <- Add_metadata(miR_TS.output, hepatitis_C.meta)

hepatitis_C.output$sample_type <- factor(
  hepatitis_C.output$sample_type,
  levels = c(
    "healthy control",
    "Hepatitis C, \nno fibrosis",
    "Hepatitis C, \nwith fibrosis"
  )
)

ggplot(hepatitis_C.output) + aes(sample_type, liver) + geom_boxplot() +
  xlab("") + ylab("Liver score") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 2), n.breaks = 5) +
  theme(axis.text = element_text(size = 12))



#----*******************************************************----
# Cohorts: SY, DFTJ, NAS----
## first, load the bulk data:
source("./code/01_data_prep/01_load cohort data.R")
library(tidyverse)
SY_miRNA_TMMcpm <- SY_miRNA_TMMcpm %>% filter(grepl("miR|let", miRNA)) %>% column_to_rownames("miRNA")
DFTJ_miRNA_TMMcpm <- DFTJ_miRNA_TMMcpm %>% filter(grepl("miR|let", miRNA)) %>% column_to_rownames("miRNA")
NAS_miRNA_TMMcpm <- NAS_miRNA_TMMcpm %>% filter(grepl("miR|let", miRNA)) %>% column_to_rownames("miRNA")

# the 17 tissues with >50% detect rate in plasma:
# Run the miRTS function:
# Note: for purpose of testing efficiency, SY_miRTS_scores.257 only uses the 257 miRNAs in signature matrix;
# The miR-TS scores using full miRNA data sets are stored in the SY_miRTS_scores, DFTJ_miRTS_scores and NAS_miRTS_scores.
temp <- miRTS_score(Input_df = SY_miRNA_TMMcpm, Norm = "none")
View(temp$proportions[,tissue_list.use_17])
plot(SY_miRTS_scores.257$heart ,  temp$proportions[,tissue_list.use_17]$heart)
all.equal(SY_miRTS_scores.257 %>% column_to_rownames("ID") %>% select(all_of(tissue_list.use_17)),
          temp$proportions[,tissue_list.use_17])
SY_miRTS_scores.257$ID <- rownames(SY_miRTS_scores.257)

temp <- miRTS_score(Input_df = DFTJ_miRNA_TMMcpm)
DFTJ_miRTS_scores.257 <- temp$proportions[,tissue_list.use_17]
DFTJ_miRTS_scores.257$ID <- rownames(DFTJ_miRTS_scores.257)

temp <- miRTS_score(Input_df = NAS_miRNA_TMMcpm)
NAS_miRTS_scores.257 <- temp$proportions[,tissue_list.use_17]
NAS_miRTS_scores.257$ID <- rownames(NAS_miRTS_scores.257)


