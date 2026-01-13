#
library(dplyr)
library(ggplot2)
library(tibble)
library(reshape2)
library(scales)  # For number_format()
library(ggpubr)

#Useful vectors:----
tissue_list.use_17 <- c(
  "adipocyte", "artery", "bladder", "bone", "bowel", "brain", "heart", "kidney",
  "liver", "lung", "lymph_node", "nerve", "pancreas", "salivary_gland",
  "skin", "spleen", "stomach"
)
tissue_list.use_17.rm_ <- gsub("_", " ", tissue_list.use_17) # underline sign removed

Plasma_miR.all.1_intersect.2 <-  read.csv("./data/other_data/miRNAs_detectable_in_plasma.csv")[[1]] # the list of plasma detectable miRNAS above 20%


# Useful functions:----
#
Transform_cibersort <- function(result1_1) { # remove all 0 and pct_of_0 > 0.5 columns, replace the 0 with /2.
  result1_1 <- result1_1[, colSums(result1_1)!=0] %>% select(-any_of(c("P-value","Correlation", "RMSE")))
  result1_1 <- result1_1[, colSums(result1_1!=0)>nrow(result1_1)*0] # 0.5
  min_tissue <- apply(result1_1,2,function(x){ifelse(min(x)==0,unique(sort(x))[2], min(x)) }) #get the second minimum value
  for (i in 1:ncol(result1_1)) {result1_1[result1_1[,i]==0,i] <- min_tissue[i]/2}
  result1_1
}

To_1st_upper <- function(mm) {
  mm <- capitalize_first(mm)
  mm <- gsub("mir-ts vs","miR-TS vs",mm)
  mm <- gsub("^ast","AST",mm,ignore.case = T)
  mm <- gsub("alt","ALT",mm,ignore.case = T)
  mm <- gsub("egfr","eGFR",mm,ignore.case = T)
  mm <- gsub("bmi","BMI",mm,ignore.case = T)
  mm <- gsub("a1c","A1c",mm,ignore.case = T)
  mm
}

replace_outliers_with_na <- function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)  # First quartile
  Q3 <- quantile(x, 0.75, na.rm = TRUE)  # Third quartile
  IQR_value <- Q3 - Q1  # Interquartile range

  lower_bound <- Q1 - 1.5 * IQR_value  # Lower bound
  upper_bound <- Q3 + 1.5 * IQR_value  # Upper bound

  x[x < lower_bound | x > upper_bound] <- NA  # Replace outliers with NA
  return(x)
}

capitalize_first <- function(x) {
  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
}

label_fn <- function(variable, value) {
  return(custom_labels[value])
}


theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size) #, base_family=base_family
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))

}
