#' This script shows how to apply the miRTS package to generate
#'   miR-TS scores for 11 public datasets.
#' The scripts for plotting are provided separately
#'   in `04_figure_generation\Fig3_validation using public datasets.R`
#'
#' All data are publicly available from NCBI GEO, with accession numbers listed.
#' All rights reserved with the publishing journals, authors and affiliations for each dataset.
#'

# set the right path to easy-load data:
current.path <- getwd()
setwd("./data/other_data/")

# hepatitis_C #####
# comparing: elevated or normal levels of alanine aminotransferase (ALT) or focal nodular hyperplasia (FNH).
miR_method="Microarray"
temp = read.csv(file = "GSE74872.csv", header = TRUE)
test <- read.table("GPL19117-74051.txt", sep = "\t")
test <- test %>% select(V1, V4)
temp <- temp%>%
  inner_join(., test, by="V1")
temp <- temp %>%
  rename(miRNA=V4) %>%
  filter(grepl("hsa", miRNA)) %>%
  select(-V1)
temp <- temp %>%
  column_to_rownames("miRNA")
temp <- 2^temp
temp[is.na(temp)] <- 0
hepatitis_C <- temp
dim(hepatitis_C)

# get the miR-TS scores:

Input_df="hepatitis_C"
assign(paste0(Input_df, ".miRTS"), miRTS_score(get(Input_df)))
View(hepatitis_C.miRTS$proportions)

#
hepatitis_C.meta <- read.csv(file = "GSE74872_meta.csv", header = TRUE)
table(hepatitis_C.meta$sample_type)
hepatitis_C.meta$sample_type <- case_when(
  hepatitis_C.meta$ishak.fiborosis.score=="control" ~ "healthy control",
  hepatitis_C.meta$ishak.fiborosis.score=="0" ~ "Hepatitis C, \nno fibrosis",
  .default = "Hepatitis C, \nwith fibrosis"

)
table(hepatitis_C.meta$sample_type)

# combine with metadata ('sample_type'):
Public_datasets.miR_TS <- assign(paste0(Input_df, "__Abs.public"),
                                 Assign_miRTS_output(
                                   Input_mix = Input_df,
                                   miRTS.obj = paste0(Input_df, ".miRTS"),
                                   convert_name = T))

dim(Public_datasets.miR_TS)
str(Public_datasets.miR_TS)



# Liver_AReject ####
miR_method="NGS"

temp = read.csv(file = "GSE69579.csv", header = TRUE)
temp$miRNA.prev <- gsub("-", "_", gsub("\\(.*\\)", "", temp$ID))
temp <- temp %>%
  inner_join(., miR.alias %>% select(miRNA.prev, miRNA), by="miRNA.prev")
temp[is.na(temp)] <- 0
temp <- temp%>% select(-ID, -miRNA.prev) %>%
  column_to_rownames("miRNA")

Liver_AReject <- temp
dim(Liver_AReject)

# create the miRNA matrix:

Input_df="Liver_AReject"
assign(paste0(Input_df, ".miRTS"), miRTS_score(get(Input_df)))

#
Liver_AReject.meta <- read.csv(file = "GSE69579_meta.csv", header = TRUE)
table(Liver_AReject.meta$sample_type)
Liver_AReject.meta$sample_type <- factor(Liver_AReject.meta$sample_type,
                                         levels = c("No rejection", "Acute rejection"),
                                         labels = c("liver allograft,\nno rejection", "liver allograft,\nacute rejection"))

# combine with metadata ('sample_type'):
Public_datasets.miR_TS <- assign(paste0(Input_df, "__Abs.public"),
                                 Assign_miRTS_output(
                                   Input_mix = Input_df,
                                   miRTS.obj = paste0(Input_df, ".miRTS"),
                                   convert_name = T))

dim(Public_datasets.miR_TS)
str(Public_datasets.miR_TS)




# liver_acetaminophen-overdose ####
miR_method="NGS"

#
temp = read.table(file = "GSE59565.txt", header = TRUE)

temp = temp %>%
  mutate(miRNA=gsub("mir", "miR", miRNA))
temp <- temp%>%   column_to_rownames("miRNA")

liver_acetaminophen <- temp
dim(liver_acetaminophen)

# create the miRNA matrix:

Input_df="liver_acetaminophen"
assign(paste0(Input_df, ".miRTS"), miRTS_score(get(Input_df)))

#
liver_acetaminophen.meta <- read.csv(file = "GSE59565_meta.csv", header = TRUE)
table(liver_acetaminophen.meta$sample_type)

# combine with metadata ('sample_type'):
Public_datasets.miR_TS <- assign(paste0(Input_df, "__Abs.public"),
                                 Assign_miRTS_output(
                                   Input_mix = Input_df,
                                   miRTS.obj = paste0(Input_df, ".miRTS"),
                                   convert_name = T))

dim(Public_datasets.miR_TS)
str(Public_datasets.miR_TS)



# COVID_sev_165 ####
miR_method="NGS"

temp = read.csv(file = "GSE178246.csv", header = TRUE)

temp <- temp %>%
  mutate(miRNA=paste0("hsa-", miRNA)) %>%
  filter(grepl("miR|let", miRNA)) %>%
  column_to_rownames("miRNA")

COVID_sev_165 <- temp
dim(COVID_sev_165)
COVID_sev_165 <- COVID_sev_165[,!is.na(COVID_sev_165[1,])]



# create the miRNA matrix:

Input_df="COVID_sev_165"
assign(paste0(Input_df, ".miRTS"), miRTS_score(get(Input_df)))

#
COVID_sev_165.meta <- read.csv("GSE178246_meta.csv")
COVID_sev_165.meta <- COVID_sev_165.meta %>%
  mutate(Col_names=paste0("X", Col_names)) %>%
  filter(!grepl("Negative|Timepoint", sample_type)) %>%
  mutate(sample_type=tolower(sample_type))

table(COVID_sev_165.meta$sample_type)

# combine with metadata ('sample_type'):
Public_datasets.miR_TS <- assign(paste0(Input_df, "__Abs.public"),
                                 Assign_miRTS_output(
                                   Input_mix = Input_df,
                                   miRTS.obj = paste0(Input_df, ".miRTS"),
                                   convert_name = T))

dim(Public_datasets.miR_TS)
str(Public_datasets.miR_TS)



# postmenopausal osteoporosis ####
# Bone metabolism-related serum miRNAs to diagnose postmenopausal osteoporosis in middle-aged and elderly women
# Agilent Human miRNA Mircoarray
miR_method="Microarray"
temp = read.csv(file = "GSE201543_RAW.csv", header = TRUE)
temp <- temp %>% filter(!is.na(GSM6067330)) #[, colSums(temp)>=100000] # Following Cell2019 workflow remove TmiR<100,000. Otherwise many samples have too few reads for the sig miRs.
colnames(temp)

temp <- temp %>%
  filter(grepl("(hsa-miR)|(hsa-let)", miRNA)) %>%
  column_to_rownames("miRNA")
# temp <- as.data.frame(apply(temp, 2, function(x) x - min(x)))
temp <- 2^temp
temp1 <- temp
osteoporosis <- temp1 %>% filter(!is.na(GSM6067330))

# create the miRNA matrix:

Input_df="osteoporosis"
assign(paste0(Input_df, ".miRTS"), miRTS_score(get(Input_df)))

#
osteoporosis.meta <- read.csv("GSE201543_meta.csv")
table(osteoporosis.meta$sample_type)
osteoporosis.meta$sample_type <- ifelse(
  grepl("without", osteoporosis.meta$sample_type), "No", "Yes"
)

# combine with metadata ('sample_type'):
Public_datasets.miR_TS <- assign(paste0(Input_df, "__Abs.public"),
                                 Assign_miRTS_output(
                                   Input_mix = Input_df,
                                   miRTS.obj = paste0(Input_df, ".miRTS"),
                                   convert_name = T))

dim(Public_datasets.miR_TS)
str(Public_datasets.miR_TS)



# T2DKD ####
# The use of plasma microRNAs (miRs) as diagnostic and prognostic biomarkers of chronic kidney disease (CKD)
miR_method="NGS"
temp = read.csv("GSE262414_mature_counts.csv")
temp <- temp#[, colSums(temp)>=100000] # Following Cell2019 workflow remove TmiR<100,000. Otherwise many samples have too few reads for the sig miRs.
colnames(temp)

temp <- temp %>% column_to_rownames("X")
rownames(temp) <- gsub("-","_",rownames(temp))
T2DKD <- temp


# create the miRNA matrix:
Input_df="T2DKD"
assign(paste0(Input_df, ".miRTS"), miRTS_score(get(Input_df)))

#
T2DKD.meta <- read.csv("GSE262414_series_matrix.csv")
table(T2DKD.meta$sample_type)
T2DKD.meta$sample_type <- case_when(
  T2DKD.meta$sample_type=="Healthy" ~ "healthy control",
  T2DKD.meta$sample_type=="T2D" ~ "T2D",
  T2DKD.meta$sample_type=="T2D_DKD" ~ "T2D, with\ndiabetic kidney disease",

)

# combine with metadata ('sample_type'):
Public_datasets.miR_TS <- assign(paste0(Input_df, "__Abs.public"),
                                 Assign_miRTS_output(
                                   Input_mix = Input_df,
                                   miRTS.obj = paste0(Input_df, ".miRTS"),
                                   convert_name = T))

dim(Public_datasets.miR_TS)
str(Public_datasets.miR_TS)



# IBD-related miR__Genome-Wide Maps of Circulating miRNA Biomarkers for Ulcerative Colitis####
# Affymetrix Genechip miRNA array
miR_method="Microarray"
temp = read.csv("IBD-related miR__Genome-Wide Maps of Circulating miRNA Biomarkers for Ulcerative Colitis.csv")
temp <- temp#[, colSums(temp)>=100000] # Following Cell2019 workflow remove TmiR<100,000. Otherwise many samples have too few reads for the sig miRs.
temp$miRNA.prev <- gsub("_star", "\\*", gsub("-", "_", temp$miRNA))
length(intersect(unique(temp$miRNA.prev), miR.alias$miRNA.prev))
temp <- temp %>%
  select(miRNA.prev, contains("Microvesicles")) %>%
  left_join(., miR.alias %>%select(miRNA, miRNA.prev), by="miRNA.prev") %>%
  select(-miRNA.prev)
temp <- temp %>% filter(!is.na(miRNA)) %>% column_to_rownames("miRNA")
rownames(temp) <- gsub("-","_",rownames(temp))
IBD.PlosOne2012 <- 2^temp
# IBD.PlosOne2012.exp <- 2^temp1


# create the miRNA matrix:

Input_df="IBD.PlosOne2012"
assign(paste0(Input_df, ".miRTS"), miRTS_score(get(Input_df)))

#
IBD.PlosOne2012.meta <- data.frame(
  Col_names=colnames(IBD.PlosOne2012),
  sample_type=ifelse(grepl("control", colnames(IBD.PlosOne2012)), "control", "UC")
)
IBD.PlosOne2012.meta$sample_type <- case_when(
  IBD.PlosOne2012.meta$sample_type=="control" ~ "healthy control",
  IBD.PlosOne2012.meta$sample_type=="UC" ~ "ulcerative colitis"
)

table(IBD.PlosOne2012.meta$sample_type)

# combine with metadata ('sample_type'):
Public_datasets.miR_TS <- assign(paste0(Input_df, "__Abs.public"),
                                 Assign_miRTS_output(
                                   Input_mix = Input_df,
                                   miRTS.obj = paste0(Input_df, ".miRTS"),
                                   convert_name = T))

dim(Public_datasets.miR_TS)
str(Public_datasets.miR_TS)




# dermatitis ####
miR_method="Microarray"
temp = read.csv(file = "GSE247297.csv", header = TRUE)

temp <- temp %>%
  filter(!is.na(GSM7886764)) %>%
  filter(grepl("(hsa-miR)|(hsa-let)", miRNA)) %>%
  column_to_rownames("miRNA")
temp <- as.data.frame(apply(temp, 2, function(x) x - min(x)))
dermatitis <- temp

# create the miRNA matrix:

Input_df="dermatitis"
assign(paste0(Input_df, ".miRTS"), miRTS_score(get(Input_df)))

#
dermatitis.meta <- read.csv("GSE247297_meta.csv")
# AD_PlosOne2015.meta$sample_type <- paste(AD_PlosOne2015.meta$fracture, AD_PlosOne2015.meta$T2DM, sep = " | ")
table(dermatitis.meta$sample_type)

dermatitis.meta$sample_type <- case_when(
  dermatitis.meta$sample_type=="patient, acute-stage" ~ "dermatitis,\nacute-stage",
  dermatitis.meta$sample_type=="patient, healed" ~ "dermatitis,\nhealed",
  .default = dermatitis.meta$sample_type
)
dermatitis.meta$sample_type <-
  factor(dermatitis.meta$sample_type,
         levels = c("healthy control", "dermatitis,\nacute-stage", "dermatitis,\nhealed"))

table(dermatitis.meta$sample_type)


# combine with metadata ('sample_type'):
Public_datasets.miR_TS <- assign(paste0(Input_df, "__Abs.public"),
                                 Assign_miRTS_output(
                                   Input_mix = Input_df,
                                   miRTS.obj = paste0(Input_df, ".miRTS"),
                                   convert_name = T))

dim(Public_datasets.miR_TS)
str(Public_datasets.miR_TS)



# adipose tissue inflammation ####
# visceral adipose tissue (VAT) inflammation was classified into low/high based on an expression score derived from the mRNA levels of TNFA, IL6 and CCL2 (determined by rtPCR).
miR_method="NGS"
temp = read.csv(file = "GSE240273.csv", header = TRUE)

adipose_inflammation <- temp %>%
  filter(!is.na(miRNA)) %>%
  column_to_rownames("miRNA")

# create the miRNA matrix:

Input_df="adipose_inflammation"
assign(paste0(Input_df, ".miRTS"), miRTS_score(get(Input_df)))


adipose_inflammation.meta <- read.csv("GSE240273_meta.csv")

# AD_PlosOne2015.meta$sample_type <- paste(AD_PlosOne2015.meta$fracture, AD_PlosOne2015.meta$T2DM, sep = " | ")
table(adipose_inflammation.meta$sample_type)

adipose_inflammation.meta$sample_type <-
  factor(adipose_inflammation.meta$sample_type,
         levels = c("low", "intermediate", "high"))

table(adipose_inflammation.meta$sample_type)


# combine with metadata ('sample_type'):
Public_datasets.miR_TS <- assign(paste0(Input_df, "__Abs.public"),
                                 Assign_miRTS_output(
                                   Input_mix = Input_df,
                                   miRTS.obj = paste0(Input_df, ".miRTS"),
                                   convert_name = T))

dim(Public_datasets.miR_TS)
str(Public_datasets.miR_TS)


# fulminant myocarditis ####
miR_method="qPCR"

temp = read.csv(file = "GSE148153.csv", header = TRUE)

temp <- temp %>%
  filter(grepl("(hsa-miR)|(hsa-let)", miRNA)) %>%
  column_to_rownames("miRNA")
temp <- as.data.frame(apply(temp, 2, function(x) x - min(x, na.rm = T)))
temp <- 2^temp - 1
temp[is.na(temp)] <- 0

fulminant_myocarditis <- temp
dim(fulminant_myocarditis)

# create the miRNA matrix:

Input_df="fulminant_myocarditis"
assign(paste0(Input_df, ".miRTS"), miRTS_score(get(Input_df)))

#
fulminant_myocarditis.meta <- read.csv(file = "GSE148153_meta.csv", header = TRUE)
fulminant_myocarditis.meta <- fulminant_myocarditis.meta %>%
  filter(!grepl("myocardial infarction", sample_type))
fulminant_myocarditis.meta$sample_type <- factor(case_when(
  fulminant_myocarditis.meta$sample_type=="FM patient" ~ "fulminant myocarditis",
  fulminant_myocarditis.meta$sample_type=="healthy control" ~ "healthy control",
  .default = NA

), levels=c("healthy control",  "fulminant myocarditis" ) )

table(fulminant_myocarditis.meta$sample_type)



# combine with metadata ('sample_type'):
Public_datasets.miR_TS <- assign(paste0(Input_df, "__Abs.public"),
                                 Assign_miRTS_output(
                                   Input_mix = Input_df,
                                   miRTS.obj = paste0(Input_df, ".miRTS"),
                                   convert_name = T))

dim(Public_datasets.miR_TS)
str(Public_datasets.miR_TS)







# traumatic brain injuries (TBI) ####
# MicroRNA sequencing of rat hippocampus and human biofluids identifies acute, chronic, focal and diffuse traumatic brain injuries
miR_method="NGS"

temp = read.csv(file = "GSE131695.csv", header = TRUE)

temp <- temp %>%
  column_to_rownames("miRNA")

TBI <- temp
dim(TBI)


# create the miRNA matrix:

Input_df="TBI"
assign(paste0(Input_df, ".miRTS"), miRTS_score(get(Input_df)))


#
TBI.meta <- read.csv(file = "GSE131695_meta.csv", header = TRUE)
TBI.meta <- TBI.meta %>%
  filter(Biofluid != "CSF"
  ) %>% mutate(
    sample_type=Timepoint,
    sample_type=ifelse(Injury.Type=="Control", "healthy control", sample_type)

  )
table(TBI.meta$sample_type)

# combine with metadata ('sample_type'):
Public_datasets.miR_TS <- assign(paste0(Input_df, "__Abs.public"),
                                 Assign_miRTS_output(
                                   Input_mix = Input_df,
                                   miRTS.obj = paste0(Input_df, ".miRTS"),
                                   convert_name = F))

dim(Public_datasets.miR_TS)
str(Public_datasets.miR_TS)




# ^_^ EOF####
setwd(current.path)
