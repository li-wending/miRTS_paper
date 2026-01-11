# Code to regenerate Figs. S1 and S2 are provided in Fig_simluation and benchmarking analyses.R

load_FigS_csvfile <- function(filename){
  temp <- read.csv(paste0("./data/Fig_results/",filename, ".csv"))
  temp
}

#Fig. S3----
## tissue specific miRNAs----
sig_mat <- miRTS::miRTS_signature_v1
dim(sig_mat)
#
# install.packages("pheatmap")
library(pheatmap)

pheatmap(
  mat = sig_mat[, tissue_list.use_17],
  scale = "row",                  # z-score per col
  color = colorRampPalette(c("navy","white","firebrick"))(100),
  clustering_method = "complete",
  show_rownames = TRUE, show_colnames = TRUE, fontsize_col = 18, fontsize_row = 5, angle_col = 45,
  # annotation_col = anno,
  #
  filename = "signature matrix - 257 miRNAs specific for 17 tissues 1.1.2026.png",   # extension controls device: .png/.pdf/.tiff
  width = 12, height = 12   ,
  border_color = NA
)

# Tissue Specificity Index (TSI) using Yanai's tau
# x: numeric vector of expression across tissues for one miRNA
tsi_tau <- function(x, na.rm = TRUE, eps = 1e-12) {
  if (na.rm) x <- x[!is.na(x)]
  if (length(x) < 2) return(NA_real_)
  # Ensure non-negative expression
  x <- pmax(x, 0)
  xmax <- max(x, na.rm = TRUE)
  if (!is.finite(xmax) || xmax <= eps) return(NA_real_)
  n <- length(x)
  sum(1 - (x / xmax)) / (n - 1)
}

# Apply to a matrix/data.frame: rows = miRNAs, cols = tissues
# Returns a named numeric vector of TSI per miRNA
calc_tsi_by_row <- function(expr_mat, na.rm = TRUE, eps = 1e-12) {
  expr_mat <- as.matrix(expr_mat)
  apply(expr_mat, 1, tsi_tau, na.rm = na.rm, eps = eps)
}
TA2_tissue.TMM <- read.csv( paste0("./data/processed_miRNA_expression/TA2_tissue.TMM.csv"))
TA2_tissue.TMM <- TA2_tissue.TMM %>% column_to_rownames("miRNA")
tsi_vals <- calc_tsi_by_row(TA2_tissue.TMM)  # rows=miRNAs, cols=tissues
# keep_miRNAs <- names(tsi_vals)[tsi_vals >= 0.8]
# sig_mat is just a subset of TA2_tissue.TMM
all.equal(TA2_tissue.TMM[rownames(sig_mat), colnames(sig_mat) ], sig_mat)
tsi_vals.sel <- calc_tsi_by_row(sig_mat)  # rows=miRNAs, cols=tissues

head(tsi_vals)

# plot:
library(ggplot2)

df <- rbind(
  data.frame(tsi = as.numeric(tsi_vals),     group = "All miRNAs"),
  data.frame(tsi = as.numeric(tsi_vals.sel), group = "Signature miRNAs")
)
df <- df[is.finite(df$tsi), , drop = FALSE]
# write.csv(df, file = "./data/Fig_results/Fig.S3-TSI_CDF.csv", row.names = F)

p <- ggplot(df, aes(x = tsi, color = group)) +
  stat_ecdf(linewidth = 1.0) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0.01, 0.01)) +
  labs(
    x = "Tissue Specificity Index (TSI)",
    y = "Empirical CDF",
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")

p

## Enriched Pathways:----

Pathway.top10 <- read.csv(paste0("./data/Fig_results/Fig.S3-Pathway enriched.csv"))

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}
scale_y_reordered <- function(..., sep = "___") {
  ggplot2::scale_y_discrete(labels = function(x) gsub(paste0(sep, ".+$"), "", x), ...)
}


Pathway.top10 <- Pathway.top10 %>%
  # plot_tbl <- tbl_reactome %>%
  group_by(tissue) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(term = reorder_within(Description, -log10(p.adjust), tissue))

ggplot(Pathway.top10 %>% filter(
  tissue %in% tissue_list.use_17

  # tissue %in% tissues[1:10] #c("heart", "brain")
  # tissue %in% tissues[11:22] #c("heart", "brain")
), aes(x = -log10(p.adjust), y = term, size = Count)) +
  geom_point(alpha = 0.85) + scale_size(range = c(0.6, 3))+
  facet_wrap(~tissue, scales = "free_y",ncol = 2) +
  scale_y_reordered() +
  labs(x = expression(-log[10]("FDR-adjusted P")), y = NULL, size = "Gene count") +
  theme_classic(base_size = 7) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

#Fig. S4----
pcm <- load_FigS_csvfile("Fig.S4 - comp to another atlas")
xx =
  ggplot(pcm, aes(x = Sample, y = variable,label=value)) +
  geom_point(aes(size = value, fill = both),
             alpha = 1, shape = 21, color="white") +
  scale_size_continuous( limits = c(0.000001, 100), range = c(1,15) )+ #,breaks = c(1,25,50), #guide = "none") +
  scale_fill_manual(values = c("red","grey"), guide = "none")+ #, guide = guide_legend(override.aes = list(size=5))
  geom_text()+ #hjust=0, vjust=0
  guides(
    size = guide_legend(
      override.aes = list(fill = "black", limits=c(1,50), color = "black") # Change legend color
    )
  )+
  labs( x= "miRNAs selected for signature matrix using miRTissueAtlas2",
        y = "miRNAs selected for \nsignature matrix using RNAAtlas",
        size = "No. of miRNAs", fill = "")  +
  theme(
    legend.key=element_blank(),
    axis.text.x = element_text(colour = "black", size = 14, angle = 45, hjust=1),
    axis.text.y = element_text(colour = "black", size = 12),
    axis.title.x = element_text(colour = "black", size = 14, face = "bold"),
    axis.title.y = element_text(colour = "black", size = 13, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 12,),
    panel.background = element_blank(),
    # plot.background = element_rect(fill = "transparent"),
    legend.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
    legend.position = "right",plot.margin = margin(0,0,0,0),
    panel.grid.major = element_line(colour = "grey90",linetype = "solid") )

xx

#Fig. S5----
# Spearman correlation between miR-TS in 3 cohorts:
# df = "SY"
library(miceafter)
for (df in c("DFTJ", "SY", "NAS")) {
  deconv.res <- get(paste0(df, "_miRTS_scores"))
  rownames(deconv.res) <- deconv.res$ID
  deconv.res$ID <- NULL
  result1_1 <- Transform_cibersort(deconv.res)
  result1_1 <- result1_1[, sort(colnames(result1_1))]
  # Calculate Spearman correlation matrix
  cor_matrix <- cor(result1_1, method = "spearman")

  # Melt the correlation matrix for visualization
  cor_melt <- reshape2::melt(cor_matrix)

  colnames(cor_melt) <- c("x_name"  ,  "y_name"   , "rxyi")
  cor_melt <- cor_melt %>%
    mutate(n=nrow(result1_1), sample_id=df)
  assign(paste0("corr_tissues.", df) , cor_melt)
}
#
corr_tissues.All <- rbind(corr_tissues.SY, corr_tissues.DFTJ, corr_tissues.NAS)
corr_tissues.All$x_name <- gsub("_", " ", corr_tissues.All$x_name)
corr_tissues.All$y_name <- gsub("_", " ", corr_tissues.All$y_name)


library(psychmeta)
coding_sheet <- corr_tissues.All %>%
  arrange(x_name, y_name)
levels(corr_tissues.All$x_name)
head(coding_sheet)
ma_res <- ma_r(rxyi = rxyi,
               n = n,
               construct_x = x_name,
               construct_y = y_name,
               sample_id = sample_id,
               data = coding_sheet
)

ma_res
summary(ma_res)
names(get_metatab(ma_res))
corr_pooled <- as.data.frame(get_metatab(ma_res))
# set values>1 to 1 on diagonal:
corr_pooled$mean_r <- ifelse(
  corr_pooled$mean_r >1, 1,
  corr_pooled$mean_r
)
corr_pooled$construct_x <- as.character(corr_pooled$construct_x)
corr_pooled$construct_y <- as.character(corr_pooled$construct_y)
corr_pooled$construct_x <- factor(corr_pooled$construct_x, rev(sort(unique(corr_pooled$construct_x))))
corr_pooled$construct_y <- factor(corr_pooled$construct_y, rev(sort(unique(corr_pooled$construct_y))))

library(tidyr)
library(dplyr)

# select one cohort:
cohort_sel = "NAS"
corr_tissues.All$x_name <- factor(corr_tissues.All$x_name,
                                  sort(gsub("_", " ", tissue_list.use_17)))
corr_tissues.All$y_name <- factor(corr_tissues.All$y_name,
                                  sort(gsub("_", " ", tissue_list.use_17)))

corr_plt <- corr_tissues.All %>% filter(sample_id == cohort_sel) %>%
  mutate(construct_x = x_name, construct_y = y_name, mean_r = rxyi) %>%
  arrange(construct_x, construct_y)

corr_matrix <- corr_plt %>%
  select(construct_x, construct_y, mean_r) %>%  # Remove unnecessary columns
  pivot_wider(names_from = construct_y, values_from = mean_r) %>%
  column_to_rownames(var = "construct_x") %>%
  as.matrix()
# View(corr_matrix)

# select all meta res:
cohort_sel = "ALL"
corr_matrix <- corr_pooled %>%
  arrange(desc(construct_x), desc(construct_y)) %>%
  select(construct_x, construct_y, mean_r) %>%  # Remove unnecessary columns
  pivot_wider(names_from = construct_y, values_from = mean_r) %>%
  column_to_rownames(var = "construct_x") %>%
  as.matrix()
# View(corr_matrix)

diag(corr_matrix) <- 1

# Use with corrplot()
library(corrplot)
col_palette <- colorRampPalette(c("blue", "white", "red"))(200)

# corrplot(corr_matrix, method = "color", type = "upper", tl.cex = 0.8)
corrplot(corr_matrix, type = "upper",# order = "hclust",
         col = col_palette,
         tl.col = "black", tl.srt = 45)


#Fig. S6----
Forest_ALL.using <- load_FigS_csvfile("Fig.S6_ForestPlot_Artery and heart")
Forest_ALL.using$study <- factor(Forest_ALL.using$study,
                       levels = rev(c("SY", "DFTJ", "NAS", "pooled")))


Organ_marker.using <-  c(
    "artery miR-TS vs. systolic blood pressure", "heart miR-TS vs. systolic blood pressure",
    "artery miR-TS vs. diastolic blood pressure", "heart miR-TS vs. diastolic blood pressure",
    "artery miR-TS vs. mean arterial pressure", "heart miR-TS vs. mean arterial pressure",
    "artery miR-TS vs. arterial stiffness (peripheral)", "heart miR-TS vs. arterial stiffness (peripheral)",
    "artery miR-TS vs. arterial stiffness (central)", "heart miR-TS vs. arterial stiffness (central)"
  )

dummy <- as.data.frame(rbind(
  c(0.6,"NAS","artery miR-TS vs. arterial stiffness (peripheral)"),
  c(-0.6,"NAS","artery miR-TS vs. arterial stiffness (peripheral)"),
  c(0.3,"NAS","artery miR-TS vs. arterial stiffness (central)"),
  c(-0.3,"NAS","artery miR-TS vs. arterial stiffness (central)"),
  c(0.5,"NAS","heart miR-TS vs. arterial stiffness (peripheral)"),
  c(-0.5,"NAS","heart miR-TS vs. arterial stiffness (peripheral)"),
  c(0.5,"NAS","heart miR-TS vs. arterial stiffness (central)"),
  c(-0.5,"NAS","heart miR-TS vs. arterial stiffness (central)"),

  c(0.5,"NAS","artery miR-TS vs. mean arterial pressure"),
  c(-0.5,"NAS","artery miR-TS vs. mean arterial pressure"),
  c(0.5,"NAS","artery miR-TS vs. systolic blood pressure"),
  c(-0.5,"NAS","artery miR-TS vs. systolic blood pressure"),
  c(0.5,"NAS","artery miR-TS vs. diastolic blood pressure"),
  c(-0.5,"NAS","artery miR-TS vs. diastolic blood pressure"),
  c(0.3,"NAS","heart miR-TS vs. mean arterial pressure"),
  c(-0.3,"NAS","heart miR-TS vs. mean arterial pressure"),
  c(0.3,"NAS","heart miR-TS vs. systolic blood pressure"),
  c(-0.3,"NAS","heart miR-TS vs. systolic blood pressure"),
  c(0.3,"NAS","heart miR-TS vs. diastolic blood pressure"),
  c(-0.3,"NAS","heart miR-TS vs. diastolic blood pressure")
))
colnames(dummy) <- c("Est", "study", "Organ_marker")
dummy$Est <- as.numeric(dummy$Est)
dummy$Organ_marker <- factor(dummy$Organ_marker, levels = Organ_marker.using)
Forest_ALL.using <- Forest_ALL.using %>%
  mutate(
    Organ_marker = factor(as.character(Organ_marker), levels = Organ_marker.using)
  )
p1 <-
  ggplot(Forest_ALL.using,aes(y = study, x = Est))+
  geom_segment(aes(x = CI_l, xend = CI_h, color=study, yend = study))+
  geom_point(aes(size=weights/10, shape=study, color=study, alpha=1))+
  theme_bw() +
  scale_alpha_identity()+
  scale_size_area()+
  facet_wrap(~Organ_marker,ncol=2
             ,scales="free_x", labeller = label_wrap_gen(width = 30, multi_line = TRUE) #
  )+ #labeller(category = label_fn)
  scale_shape_manual(values = rev(c(15, 15, 15, 18))) +
  scale_color_manual(values = rev(c("#2d89c9", "#39b592", "#e6a23e", "#1e3135" )))+
  geom_vline(lty=2, aes(xintercept=ref_line), colour = 'red') +
  geom_blank(data=dummy) +
  theme(strip.text.x = element_text(size = 14),
        axis.title = element_blank(),
        legend.position = "none")

p1
ggsave("test.png", width = 6, height = 8,units = "in",  dpi = 300)

#####
corr_tissues.All <- load_FigS_csvfile("Fig.S5-corr_tissues.All")
corr_pooled <- load_FigS_csvfile("Fig.S5-corr_pooled")


# select one cohort:
cohort_sel = "SY"
corr_tissues.All$x_name <- factor(corr_tissues.All$x_name,
                                  sort(gsub("_", " ", tissue_list.use_17)))
corr_tissues.All$y_name <- factor(corr_tissues.All$y_name,
                                  sort(gsub("_", " ", tissue_list.use_17)))

corr_plt <- corr_tissues.All %>% filter(sample_id == cohort_sel) %>%
  mutate(construct_x = x_name, construct_y = y_name, mean_r = rxyi) %>%
  arrange(construct_x, construct_y)

corr_matrix <- corr_plt %>%
  select(construct_x, construct_y, mean_r) %>%  # Remove unnecessary columns
  pivot_wider(names_from = construct_y, values_from = mean_r) %>%
  column_to_rownames(var = "construct_x") %>%
  as.matrix()



library(corrplot)
col_palette <- colorRampPalette(c("blue", "white", "red"))(200)
corrplot(corr_matrix, type = "upper",# order = "hclust",
         col = col_palette,
         tl.col = "black", tl.srt = 45)



# select one cohort:
cohort_sel = "DFTJ"
corr_tissues.All$x_name <- factor(corr_tissues.All$x_name,
                                  sort(gsub("_", " ", tissue_list.use_17)))
corr_tissues.All$y_name <- factor(corr_tissues.All$y_name,
                                  sort(gsub("_", " ", tissue_list.use_17)))

corr_plt <- corr_tissues.All %>% filter(sample_id == cohort_sel) %>%
  mutate(construct_x = x_name, construct_y = y_name, mean_r = rxyi) %>%
  arrange(construct_x, construct_y)

corr_matrix <- corr_plt %>%
  select(construct_x, construct_y, mean_r) %>%  # Remove unnecessary columns
  pivot_wider(names_from = construct_y, values_from = mean_r) %>%
  column_to_rownames(var = "construct_x") %>%
  as.matrix()



library(corrplot)
col_palette <- colorRampPalette(c("blue", "white", "red"))(200)
corrplot(corr_matrix, type = "upper",# order = "hclust",
         col = col_palette,
         tl.col = "black", tl.srt = 45)



# select one cohort:
cohort_sel = "NAS"
corr_tissues.All$x_name <- factor(corr_tissues.All$x_name,
                                  sort(gsub("_", " ", tissue_list.use_17)))
corr_tissues.All$y_name <- factor(corr_tissues.All$y_name,
                                  sort(gsub("_", " ", tissue_list.use_17)))

corr_plt <- corr_tissues.All %>% filter(sample_id == cohort_sel) %>%
  mutate(construct_x = x_name, construct_y = y_name, mean_r = rxyi) %>%
  arrange(construct_x, construct_y)

corr_matrix <- corr_plt %>%
  select(construct_x, construct_y, mean_r) %>%  # Remove unnecessary columns
  pivot_wider(names_from = construct_y, values_from = mean_r) %>%
  column_to_rownames(var = "construct_x") %>%
  as.matrix()



library(corrplot)
col_palette <- colorRampPalette(c("blue", "white", "red"))(200)
corrplot(corr_matrix, type = "upper",# order = "hclust",
         col = col_palette,
         tl.col = "black", tl.srt = 45)



#
# # select all meta res:
cohort_sel = "ALL"
corr_matrix <- corr_pooled %>%
  arrange((construct_x), (construct_y)) %>%
  select(construct_x, construct_y, mean_r) %>%  # Remove unnecessary columns
  pivot_wider(names_from = construct_y, values_from = mean_r) %>%
  column_to_rownames(var = "construct_x") %>%
  as.matrix()
# View(corr_matrix)


library(corrplot)
col_palette <- colorRampPalette(c("blue", "white", "red"))(200)
corrplot(corr_matrix, type = "upper",# order = "hclust",
         col = col_palette,
         tl.col = "black", tl.srt = 45)

#Fig. S7----
Pooled.ht.TSI.forbioimarkers <- load_FigS_csvfile("Fig.S7-Pooled.ht.TSI.forbioimarkers")
table(Pooled.ht.TSI.forbioimarkers$study)
Pooled.ht.TSI.formiRTS <- load_FigS_csvfile("Fig.S7-Pooled.ht.TSI.formiRTS")

## -with clinical biomarkers----
cohort_sel = "SY"
Pooled.ht.TSI.htmp <- Pooled.ht.TSI.forbioimarkers %>%
  filter(study==cohort_sel) %>%
  dplyr::mutate(
    x=gsub(" miR-TS vs. .*", "", Organ_marker),
    y=gsub("Lung", "lung", gsub(".* miR-TS vs. ", "", Organ_marker)),
    y=ifelse(y=="lung decline", "lung function decline", y),
    y=To_1st_upper(y)
  )

df_heatmap <- Pooled.ht.TSI.htmp %>%
  filter(P!=0) %>%
  group_by(y) %>%
  mutate(
    Padj=p.adjust(P, method="BH")) %>%
  mutate(
    P_trans=-log10(P),
    t_adj=case_when(
      Padj<0.01~ "111", Padj<0.05~ "11", Padj<0.2~ "1",
      .default = "") )

{
  df_heatmap_t <- as.data.frame(
    reshape2::dcast(df_heatmap, y~x, value.var = "t") %>%
      column_to_rownames("y")
  )
  temp_1 <- as.data.frame(
    reshape2::dcast(df_heatmap, y~x, value.var = "t_adj") %>%
      column_to_rownames("y")
  )
  temp_1[is.na(temp_1)] <- ""


  order_marker <- rev(c(
    "%body fat",
    "Waist/height",
    "BMI",
    "Arterial stiffness",
    "Mean arterial pressure",#
    "Cognitive impairment",
    "Coronary heart diseases",
    "eGFR", "ALT", "AST",
    "Lung function decline",
    "Lung diseases",
    "%lymphocytes",
    "Diabetes",
    "Hemoglobin A1c", #
    "Fasting blood glucose"
  ))
  order_marker.new <- intersect(order_marker, rownames(df_heatmap_t))
  df_heatmap_t <- df_heatmap_t[order_marker.new, ]
  temp_1 <- temp_1[order_marker.new, ]

  library(ComplexHeatmap)
  library(circlize)
  (t_max <- max(abs(min(df_heatmap_t, na.rm = T)), 0, abs(max(df_heatmap_t, na.rm = T))))
  newline_counts <- sapply(gregexpr("\n", rownames(df_heatmap_t)), function(x) ifelse(x[1] == -1, 0, length(x)))

  row_heights <- unit(newline_counts*3, "cm")


  ht <-
    Heatmap(as.matrix(df_heatmap_t[]),
            col = colorRamp2(c(-t_max, 0, t_max), c("blue", "white", "red")),
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            show_row_dend = F,
            row_names_side = 'left',
            row_dend_reorder = F,
            column_dend_reorder = F,
            column_names_rot = 45,
            column_names_centered = F,
            border = 1,
            row_title = "phenotype",
            column_title_gp = gpar(fontsize = 14, fontface = "bold"),
            heatmap_legend_param = list(
              title="t",
              legend_width = unit(2, "cm")),
            show_heatmap_legend = T,
            cell_fun = function(j, i, x, y, w, h, fill) {
              if(temp_1[i, j] =="***") {
                grid.text("***", x, y)
              } else if(temp_1[i, j] =="**") {
                grid.text("**", x, y)
              } else if(temp_1[i, j] =="*") {
                grid.text("*", x, y)
              } else if(temp_1[i, j] =="111") {
                grid.text("***", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else if(temp_1[i, j] =="11") {
                grid.text("**", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else if(temp_1[i, j] =="1") {
                grid.text("*", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else {
                grid.text("", x, y)
              }}
    )

  }
ht

#
cohort_sel = "DFTJ"
Pooled.ht.TSI.htmp <- Pooled.ht.TSI.forbioimarkers %>%
  filter(study==cohort_sel) %>%
  dplyr::mutate(
    x=gsub(" miR-TS vs. .*", "", Organ_marker),
    y=gsub("Lung", "lung", gsub(".* miR-TS vs. ", "", Organ_marker)),
    y=ifelse(y=="lung decline", "lung function decline", y),
    y=To_1st_upper(y)
  )

df_heatmap <- Pooled.ht.TSI.htmp %>%
  filter(P!=0) %>%
  group_by(y) %>%
  mutate(
    Padj=p.adjust(P, method="BH")) %>%
  mutate(
    P_trans=-log10(P),
    t_adj=case_when(
      Padj<0.01~ "111", Padj<0.05~ "11", Padj<0.2~ "1",
      .default = "") )

{
  df_heatmap_t <- as.data.frame(
    reshape2::dcast(df_heatmap, y~x, value.var = "t") %>%
      column_to_rownames("y")
  )
  temp_1 <- as.data.frame(
    reshape2::dcast(df_heatmap, y~x, value.var = "t_adj") %>%
      column_to_rownames("y")
  )
  temp_1[is.na(temp_1)] <- ""


  order_marker <- rev(c(
    "%body fat",
    "Waist/height",
    "BMI",
    "Arterial stiffness",
    "Mean arterial pressure",#
    "Cognitive impairment",
    "Coronary heart diseases",
    "eGFR", "ALT", "AST",
    "Lung function decline",
    "Lung diseases",
    "%lymphocytes",
    "Diabetes",
    "Hemoglobin A1c", #
    "Fasting blood glucose"
  ))
  order_marker.new <- intersect(order_marker, rownames(df_heatmap_t))
  df_heatmap_t <- df_heatmap_t[order_marker.new, ]
  temp_1 <- temp_1[order_marker.new, ]

  library(ComplexHeatmap)
  library(circlize)
  (t_max <- max(abs(min(df_heatmap_t, na.rm = T)), 0, abs(max(df_heatmap_t, na.rm = T))))
  newline_counts <- sapply(gregexpr("\n", rownames(df_heatmap_t)), function(x) ifelse(x[1] == -1, 0, length(x)))

  row_heights <- unit(newline_counts*3, "cm")


  ht <-
    Heatmap(as.matrix(df_heatmap_t[]),
            col = colorRamp2(c(-t_max, 0, t_max), c("blue", "white", "red")),
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            show_row_dend = F,
            row_names_side = 'left',
            row_dend_reorder = F,
            column_dend_reorder = F,
            column_names_rot = 45,
            column_names_centered = F,
            border = 1,
            row_title = "phenotype",
            column_title_gp = gpar(fontsize = 14, fontface = "bold"),
            heatmap_legend_param = list(
              title="t",
              legend_width = unit(2, "cm")),
            show_heatmap_legend = T,
            cell_fun = function(j, i, x, y, w, h, fill) {
              if(temp_1[i, j] =="***") {
                grid.text("***", x, y)
              } else if(temp_1[i, j] =="**") {
                grid.text("**", x, y)
              } else if(temp_1[i, j] =="*") {
                grid.text("*", x, y)
              } else if(temp_1[i, j] =="111") {
                grid.text("***", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else if(temp_1[i, j] =="11") {
                grid.text("**", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else if(temp_1[i, j] =="1") {
                grid.text("*", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else {
                grid.text("", x, y)
              }}
    )

  }
ht

#
cohort_sel = "NAS"
Pooled.ht.TSI.htmp <- Pooled.ht.TSI.forbioimarkers %>%
  filter(study==cohort_sel) %>%
  dplyr::mutate(
    x=gsub(" miR-TS vs. .*", "", Organ_marker),
    y=gsub("Lung", "lung", gsub(".* miR-TS vs. ", "", Organ_marker)),
    y=ifelse(y=="lung decline", "lung function decline", y),
    y=To_1st_upper(y)
  )

df_heatmap <- Pooled.ht.TSI.htmp %>%
  filter(P!=0) %>%
  group_by(y) %>%
  mutate(
    Padj=p.adjust(P, method="BH")) %>%
  mutate(
    P_trans=-log10(P),
    t_adj=case_when(
      Padj<0.01~ "111", Padj<0.05~ "11", Padj<0.2~ "1",
      .default = "") )

{
  df_heatmap_t <- as.data.frame(
    reshape2::dcast(df_heatmap, y~x, value.var = "t") %>%
      column_to_rownames("y")
  )
  temp_1 <- as.data.frame(
    reshape2::dcast(df_heatmap, y~x, value.var = "t_adj") %>%
      column_to_rownames("y")
  )
  temp_1[is.na(temp_1)] <- ""


  order_marker <- rev(c(
    "%body fat",
    "Waist/height",
    "BMI",
    "Arterial stiffness",
    "Mean arterial pressure",#
    "Cognitive impairment",
    "Coronary heart diseases",
    "eGFR", "ALT", "AST",
    "Lung function decline",
    "Lung diseases",
    "%lymphocytes",
    "Diabetes",
    "Hemoglobin A1c", #
    "Fasting blood glucose"
  ))
  order_marker.new <- intersect(order_marker, rownames(df_heatmap_t))
  df_heatmap_t <- df_heatmap_t[order_marker.new, ]
  temp_1 <- temp_1[order_marker.new, ]

  library(ComplexHeatmap)
  library(circlize)
  (t_max <- max(abs(min(df_heatmap_t, na.rm = T)), 0, abs(max(df_heatmap_t, na.rm = T))))
  newline_counts <- sapply(gregexpr("\n", rownames(df_heatmap_t)), function(x) ifelse(x[1] == -1, 0, length(x)))

  row_heights <- unit(newline_counts*3, "cm")


  ht <-
    Heatmap(as.matrix(df_heatmap_t[]),
            col = colorRamp2(c(-t_max, 0, t_max), c("blue", "white", "red")),
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            show_row_dend = F,
            row_names_side = 'left',
            row_dend_reorder = F,
            column_dend_reorder = F,
            column_names_rot = 45,
            column_names_centered = F,
            border = 1,
            row_title = "phenotype",
            column_title_gp = gpar(fontsize = 14, fontface = "bold"),
            heatmap_legend_param = list(
              title="t",
              legend_width = unit(2, "cm")),
            show_heatmap_legend = T,
            cell_fun = function(j, i, x, y, w, h, fill) {
              if(temp_1[i, j] =="***") {
                grid.text("***", x, y)
              } else if(temp_1[i, j] =="**") {
                grid.text("**", x, y)
              } else if(temp_1[i, j] =="*") {
                grid.text("*", x, y)
              } else if(temp_1[i, j] =="111") {
                grid.text("***", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else if(temp_1[i, j] =="11") {
                grid.text("**", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else if(temp_1[i, j] =="1") {
                grid.text("*", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else {
                grid.text("", x, y)
              }}
    )

  }
ht

#
cohort_sel = "pooled"
Pooled.ht.TSI.htmp <- Pooled.ht.TSI.forbioimarkers %>%
  filter(study==cohort_sel) %>%
  dplyr::mutate(
    x=gsub(" miR-TS vs. .*", "", Organ_marker),
    y=gsub("Lung", "lung", gsub(".* miR-TS vs. ", "", Organ_marker)),
    y=ifelse(y=="lung decline", "lung function decline", y),
    y=To_1st_upper(y)
  )

df_heatmap <- Pooled.ht.TSI.htmp %>%
  filter(P!=0) %>%
  group_by(y) %>%
  mutate(
    Padj=p.adjust(P, method="BH")) %>%
  mutate(
    P_trans=-log10(P),
    t_adj=case_when(
      Padj<0.01~ "111", Padj<0.05~ "11", Padj<0.2~ "1",
      .default = "") )

{
  df_heatmap_t <- as.data.frame(
    reshape2::dcast(df_heatmap, y~x, value.var = "t") %>%
      column_to_rownames("y")
  )
  temp_1 <- as.data.frame(
    reshape2::dcast(df_heatmap, y~x, value.var = "t_adj") %>%
      column_to_rownames("y")
  )
  temp_1[is.na(temp_1)] <- ""


  order_marker <- rev(c(
    "%body fat",
    "Waist/height",
    "BMI",
    "Arterial stiffness",
    "Mean arterial pressure",#
    "Cognitive impairment",
    "Coronary heart diseases",
    "eGFR", "ALT", "AST",
    "Lung function decline",
    "Lung diseases",
    "%lymphocytes",
    "Diabetes",
    "Hemoglobin A1c", #
    "Fasting blood glucose"
  ))
  order_marker.new <- intersect(order_marker, rownames(df_heatmap_t))
  df_heatmap_t <- df_heatmap_t[order_marker.new, ]
  temp_1 <- temp_1[order_marker.new, ]

  library(ComplexHeatmap)
  library(circlize)
  (t_max <- max(abs(min(df_heatmap_t, na.rm = T)), 0, abs(max(df_heatmap_t, na.rm = T))))
  newline_counts <- sapply(gregexpr("\n", rownames(df_heatmap_t)), function(x) ifelse(x[1] == -1, 0, length(x)))

  row_heights <- unit(newline_counts*3, "cm")


  ht <-
    Heatmap(as.matrix(df_heatmap_t[]),
            col = colorRamp2(c(-t_max, 0, t_max), c("blue", "white", "red")),
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            show_row_dend = F,
            row_names_side = 'left',
            row_dend_reorder = F,
            column_dend_reorder = F,
            column_names_rot = 45,
            column_names_centered = F,
            border = 1,
            row_title = "phenotype",
            column_title_gp = gpar(fontsize = 14, fontface = "bold"),
            heatmap_legend_param = list(
              title="t",
              legend_width = unit(2, "cm")),
            show_heatmap_legend = T,
            cell_fun = function(j, i, x, y, w, h, fill) {
              if(temp_1[i, j] =="***") {
                grid.text("***", x, y)
              } else if(temp_1[i, j] =="**") {
                grid.text("**", x, y)
              } else if(temp_1[i, j] =="*") {
                grid.text("*", x, y)
              } else if(temp_1[i, j] =="111") {
                grid.text("***", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else if(temp_1[i, j] =="11") {
                grid.text("**", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else if(temp_1[i, j] =="1") {
                grid.text("*", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else {
                grid.text("", x, y)
              }}
    )

  }
ht

## -with miR-TS scores----

#
cohort_sel = "SY"
Pooled.ht.TSI.htmp <- Pooled.ht.TSI.formiRTS %>%
  filter(study==cohort_sel) %>%
  mutate(
    x=gsub(" miR-TS vs. .*", "", Organ_marker),
    y=gsub("Lung", "lung", gsub(".* miR-TS vs. ", "", Organ_marker)),
  )

df_heatmap <- Pooled.ht.TSI.htmp %>%
  filter(P!=0) %>%
  group_by(y) %>%
  mutate(
    Padj=p.adjust(P, method="BH")) %>%
  mutate(
    P_trans=-log10(P),
    P_trans=pmin(P_trans, 10),
    t_adj=case_when(#P<0.001~ "***", P<0.01~ "**", P<0.05~ "*",
      Padj<0.01~ "111", Padj<0.05~ "11", Padj<0.2~ "1",
      .default = "") )

{
  df_heatmap_t <- as.data.frame(
    reshape2::dcast(df_heatmap, y~x, value.var = "t") %>%
      column_to_rownames("y")
  )
  temp_1 <- as.data.frame(
    reshape2::dcast(df_heatmap, y~x, value.var = "t_adj") %>%
      column_to_rownames("y")
  )
  temp_1[is.na(temp_1)] <- ""
  rownames.order <- rev(rownames(df_heatmap_t))
  df_heatmap_t <- df_heatmap_t[rownames.order, ]

  temp_1 <- temp_1[rownames.order, ]

  library(ComplexHeatmap)
  library(circlize)
  (t_max <- max(abs(min(df_heatmap_t, na.rm = T)), 0, abs(max(df_heatmap_t, na.rm = T))))
  newline_counts <- sapply(gregexpr("\n", rownames(df_heatmap_t)), function(x) ifelse(x[1] == -1, 0, length(x)))

  row_heights <- unit(newline_counts*3, "cm")


  ht <-
    Heatmap(as.matrix(df_heatmap_t[]),
            col = colorRamp2(c(-t_max, 0, t_max), c("blue", "white", "red")),
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            show_row_dend = F,
            row_names_side = 'left',
            row_dend_reorder = F,
            column_dend_reorder = F,
            column_names_rot = 45,
            column_names_centered = F,
            border = 1, row_title = "miR-TS",
            column_title_gp = gpar(fontsize = 14, fontface = "bold"),
            heatmap_legend_param = list(
              title="t",
              legend_width = unit(2, "cm")),
            show_heatmap_legend = T,
            cell_fun = function(j, i, x, y, w, h, fill) {
              if(temp_1[i, j] =="***") {
                grid.text("***", x, y)
              } else if(temp_1[i, j] =="**") {
                grid.text("**", x, y)
              } else if(temp_1[i, j] =="*") {
                grid.text("*", x, y)
              } else if(temp_1[i, j] =="111") {
                grid.text("***", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else if(temp_1[i, j] =="11") {
                grid.text("**", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else if(temp_1[i, j] =="1") {
                grid.text("*", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else {
                grid.text("", x, y)
              }}
    )

  }
ht


#
cohort_sel = "DFTJ"
Pooled.ht.TSI.htmp <- Pooled.ht.TSI.formiRTS %>%
  filter(study==cohort_sel) %>%
  mutate(
    x=gsub(" miR-TS vs. .*", "", Organ_marker),
    y=gsub("Lung", "lung", gsub(".* miR-TS vs. ", "", Organ_marker)),
  )

df_heatmap <- Pooled.ht.TSI.htmp %>%
  filter(P!=0) %>%
  group_by(y) %>%
  mutate(
    Padj=p.adjust(P, method="BH")) %>%
  mutate(
    P_trans=-log10(P),
    P_trans=pmin(P_trans, 10),
    t_adj=case_when(#P<0.001~ "***", P<0.01~ "**", P<0.05~ "*",
      Padj<0.01~ "111", Padj<0.05~ "11", Padj<0.2~ "1",
      .default = "") )

{
  df_heatmap_t <- as.data.frame(
    reshape2::dcast(df_heatmap, y~x, value.var = "t") %>%
      column_to_rownames("y")
  )
  temp_1 <- as.data.frame(
    reshape2::dcast(df_heatmap, y~x, value.var = "t_adj") %>%
      column_to_rownames("y")
  )
  temp_1[is.na(temp_1)] <- ""
  rownames.order <- rev(rownames(df_heatmap_t))
  df_heatmap_t <- df_heatmap_t[rownames.order, ]

  temp_1 <- temp_1[rownames.order, ]

  library(ComplexHeatmap)
  library(circlize)
  (t_max <- max(abs(min(df_heatmap_t, na.rm = T)), 0, abs(max(df_heatmap_t, na.rm = T))))
  newline_counts <- sapply(gregexpr("\n", rownames(df_heatmap_t)), function(x) ifelse(x[1] == -1, 0, length(x)))

  row_heights <- unit(newline_counts*3, "cm")


  ht <-
    Heatmap(as.matrix(df_heatmap_t[]),
            col = colorRamp2(c(-t_max, 0, t_max), c("blue", "white", "red")),
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            show_row_dend = F,
            row_names_side = 'left',
            row_dend_reorder = F,
            column_dend_reorder = F,
            column_names_rot = 45,
            column_names_centered = F,
            border = 1, row_title = "miR-TS",
            column_title_gp = gpar(fontsize = 14, fontface = "bold"),
            heatmap_legend_param = list(
              title="t",
              legend_width = unit(2, "cm")),
            show_heatmap_legend = T,
            cell_fun = function(j, i, x, y, w, h, fill) {
              if(temp_1[i, j] =="***") {
                grid.text("***", x, y)
              } else if(temp_1[i, j] =="**") {
                grid.text("**", x, y)
              } else if(temp_1[i, j] =="*") {
                grid.text("*", x, y)
              } else if(temp_1[i, j] =="111") {
                grid.text("***", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else if(temp_1[i, j] =="11") {
                grid.text("**", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else if(temp_1[i, j] =="1") {
                grid.text("*", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else {
                grid.text("", x, y)
              }}
    )

  }
ht

#
cohort_sel = "NAS"
Pooled.ht.TSI.htmp <- Pooled.ht.TSI.formiRTS %>%
  filter(study==cohort_sel) %>%
  mutate(
    x=gsub(" miR-TS vs. .*", "", Organ_marker),
    y=gsub("Lung", "lung", gsub(".* miR-TS vs. ", "", Organ_marker)),
  )

df_heatmap <- Pooled.ht.TSI.htmp %>%
  filter(P!=0) %>%
  group_by(y) %>%
  mutate(
    Padj=p.adjust(P, method="BH")) %>%
  mutate(
    P_trans=-log10(P),
    P_trans=pmin(P_trans, 10),
    t_adj=case_when(#P<0.001~ "***", P<0.01~ "**", P<0.05~ "*",
      Padj<0.01~ "111", Padj<0.05~ "11", Padj<0.2~ "1",
      .default = "") )

{
  df_heatmap_t <- as.data.frame(
    reshape2::dcast(df_heatmap, y~x, value.var = "t") %>%
      column_to_rownames("y")
  )
  temp_1 <- as.data.frame(
    reshape2::dcast(df_heatmap, y~x, value.var = "t_adj") %>%
      column_to_rownames("y")
  )
  temp_1[is.na(temp_1)] <- ""
  rownames.order <- rev(rownames(df_heatmap_t))
  df_heatmap_t <- df_heatmap_t[rownames.order, ]

  temp_1 <- temp_1[rownames.order, ]

  library(ComplexHeatmap)
  library(circlize)
  (t_max <- max(abs(min(df_heatmap_t, na.rm = T)), 0, abs(max(df_heatmap_t, na.rm = T))))
  newline_counts <- sapply(gregexpr("\n", rownames(df_heatmap_t)), function(x) ifelse(x[1] == -1, 0, length(x)))

  row_heights <- unit(newline_counts*3, "cm")


  ht <-
    Heatmap(as.matrix(df_heatmap_t[]),
            col = colorRamp2(c(-t_max, 0, t_max), c("blue", "white", "red")),
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            show_row_dend = F,
            row_names_side = 'left',
            row_dend_reorder = F,
            column_dend_reorder = F,
            column_names_rot = 45,
            column_names_centered = F,
            border = 1, row_title = "miR-TS",
            column_title_gp = gpar(fontsize = 14, fontface = "bold"),
            heatmap_legend_param = list(
              title="t",
              legend_width = unit(2, "cm")),
            show_heatmap_legend = T,
            cell_fun = function(j, i, x, y, w, h, fill) {
              if(temp_1[i, j] =="***") {
                grid.text("***", x, y)
              } else if(temp_1[i, j] =="**") {
                grid.text("**", x, y)
              } else if(temp_1[i, j] =="*") {
                grid.text("*", x, y)
              } else if(temp_1[i, j] =="111") {
                grid.text("***", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else if(temp_1[i, j] =="11") {
                grid.text("**", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else if(temp_1[i, j] =="1") {
                grid.text("*", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else {
                grid.text("", x, y)
              }}
    )

  }
ht

#
cohort_sel = "pooled"
Pooled.ht.TSI.htmp <- Pooled.ht.TSI.formiRTS %>%
  filter(study==cohort_sel) %>%
  mutate(
    x=gsub(" miR-TS vs. .*", "", Organ_marker),
    y=gsub("Lung", "lung", gsub(".* miR-TS vs. ", "", Organ_marker)),
  )

df_heatmap <- Pooled.ht.TSI.htmp %>%
  filter(P!=0) %>%
  group_by(y) %>%
  mutate(
    Padj=p.adjust(P, method="BH")) %>%
  mutate(
    P_trans=-log10(P),
    P_trans=pmin(P_trans, 10),
    t_adj=case_when(#P<0.001~ "***", P<0.01~ "**", P<0.05~ "*",
      Padj<0.01~ "111", Padj<0.05~ "11", Padj<0.2~ "1",
      .default = "") )

{
  df_heatmap_t <- as.data.frame(
    reshape2::dcast(df_heatmap, y~x, value.var = "t") %>%
      column_to_rownames("y")
  )
  temp_1 <- as.data.frame(
    reshape2::dcast(df_heatmap, y~x, value.var = "t_adj") %>%
      column_to_rownames("y")
  )
  temp_1[is.na(temp_1)] <- ""
  rownames.order <- rev(rownames(df_heatmap_t))
  df_heatmap_t <- df_heatmap_t[rownames.order, ]

  temp_1 <- temp_1[rownames.order, ]

  library(ComplexHeatmap)
  library(circlize)
  (t_max <- max(abs(min(df_heatmap_t, na.rm = T)), 0, abs(max(df_heatmap_t, na.rm = T))))
  newline_counts <- sapply(gregexpr("\n", rownames(df_heatmap_t)), function(x) ifelse(x[1] == -1, 0, length(x)))

  row_heights <- unit(newline_counts*3, "cm")


  ht <-
    Heatmap(as.matrix(df_heatmap_t[]),
            col = colorRamp2(c(-t_max, 0, t_max), c("blue", "white", "red")),
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            show_row_dend = F,
            row_names_side = 'left',
            row_dend_reorder = F,
            column_dend_reorder = F,
            column_names_rot = 45,
            column_names_centered = F,
            border = 1, row_title = "miR-TS",
            column_title_gp = gpar(fontsize = 14, fontface = "bold"),
            heatmap_legend_param = list(
              title="t",
              legend_width = unit(2, "cm")),
            show_heatmap_legend = T,
            cell_fun = function(j, i, x, y, w, h, fill) {
              if(temp_1[i, j] =="***") {
                grid.text("***", x, y)
              } else if(temp_1[i, j] =="**") {
                grid.text("**", x, y)
              } else if(temp_1[i, j] =="*") {
                grid.text("*", x, y)
              } else if(temp_1[i, j] =="111") {
                grid.text("***", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else if(temp_1[i, j] =="11") {
                grid.text("**", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else if(temp_1[i, j] =="1") {
                grid.text("*", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else {
                grid.text("", x, y)
              }}
    )

  }
ht


# png(paste("cancer_plot/test_cancer -eachCancer", "ALL", "linear", "lm", Cancer_type_using, "_", logTrans, " .png", sep = "_"), # As and OCM&FA metabolites__batchAdjusted.baseline_PBO.png
# png(paste("Heatmap-Score_AllOrgan_FDR_0.2 less_noAbbr", "20241122 .png", sep = "_"), # As and OCM&FA metabolites__batchAdjusted.baseline_PBO.png
png(paste("Heatmap-TSI0.8 - with miR-TS tissue- ",cohort_sel, "-20260104 .png", sep = "_"), # As and OCM&FA metabolites__batchAdjusted.baseline_PBO.png
    width=10,height=5,units="in",res=300)
draw(ht, padding =  unit(c(2, 45, 2, 2), "mm"))

dev.off()



#Fig. S8----
Pooled.ht.top1.forbioimarkers <- load_FigS_csvfile("Fig.S8-df_heatmap.biomarkers.using_top1")
Pooled.ht.top1.formiRTSscore <- load_FigS_csvfile("Fig.S8-df_heatmap.miRTSscore.using_top1")
Pooled.ht.top2.forbioimarkers <- load_FigS_csvfile("Fig.S8-df_heatmap.biomarkers.using_top2")
Pooled.ht.top2.formiRTSscore <- load_FigS_csvfile("Fig.S8-df_heatmap.miRTSscore.using_top2")
Pooled.ht.top3.forbioimarkers <- load_FigS_csvfile("Fig.S8-df_heatmap.biomarkers.using_top3")
Pooled.ht.top3.formiRTSscore <- load_FigS_csvfile("Fig.S8-df_heatmap.miRTSscore.using_top3")


## -with clinical biomarkers----
for (plot.df.use in c(
  "Pooled.ht.top1.forbioimarkers", "Pooled.ht.top2.forbioimarkers", "Pooled.ht.top3.forbioimarkers"
)) {
  df_heatmap <- get(plot.df.use)
  df_heatmap_t <- as.data.frame(
    reshape2::dcast(df_heatmap, y~x, value.var = "P_trans") %>%
      column_to_rownames("y")
  )
  temp_1 <- as.data.frame(
    reshape2::dcast(df_heatmap, y~x, value.var = "t_adj") %>%
      column_to_rownames("y")
  )
  temp_1[is.na(temp_1)] <- ""
  rownames(df_heatmap_t)

  df_heatmap_t <- df_heatmap_t[order_marker, ]

  temp_1 <- temp_1[order_marker, ]

  library(ComplexHeatmap)
  library(circlize)
  (t_max <- max(abs(min(df_heatmap_t, na.rm = T)), 0, abs(max(df_heatmap_t, na.rm = T))))
  newline_counts <- sapply(gregexpr("\n", rownames(df_heatmap_t)), function(x) ifelse(x[1] == -1, 0, length(x)))

  row_heights <- unit(newline_counts*3, "cm")


  ht1 <-
    Heatmap(as.matrix(df_heatmap_t[]),
            col = colorRamp2(c( 0, t_max), c( "white", "red")),
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            show_row_dend = F,
            row_names_side = 'left',
            row_dend_reorder = F,
            column_dend_reorder = F,
            column_names_rot = 45,
            column_names_centered = F,
            border = 1,
            column_title_gp = gpar(fontsize = 14, fontface = "bold"),
            heatmap_legend_param = list(
              title="-log10(P)",
              legend_width = unit(2, "cm")),
            show_heatmap_legend = T,
            cell_fun = function(j, i, x, y, w, h, fill) {
              if(temp_1[i, j] =="***") {
                grid.text("***", x, y)
              } else if(temp_1[i, j] =="**") {
                grid.text("**", x, y)
              } else if(temp_1[i, j] =="*") {
                grid.text("*", x, y)
              } else if(temp_1[i, j] =="111") {
                grid.text("***", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else if(temp_1[i, j] =="11") {
                grid.text("**", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else if(temp_1[i, j] =="1") {
                grid.text("*", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else {
                grid.text("", x, y)
              }}
    )
  print(ht1)

}

## -with miR-TS scores----

for (plot.df.use in c(
  "Pooled.ht.top1.formiRTSscore", "Pooled.ht.top2.formiRTSscore", "Pooled.ht.top3.formiRTSscore"
)) {
  df_heatmap <- get(plot.df.use)
  df_heatmap$y <- gsub("_", " ", df_heatmap$y)
  df_heatmap_t <- as.data.frame(
    reshape2::dcast(df_heatmap, y~x, value.var = "P_trans") %>%
      column_to_rownames("y")
  )
  temp_1 <- as.data.frame(
    reshape2::dcast(df_heatmap, y~x, value.var = "t_adj") %>%
      column_to_rownames("y")
  )
  temp_1[is.na(temp_1)] <- ""
  rownames.order <- rev(rownames(df_heatmap_t))
  df_heatmap_t <- df_heatmap_t[rownames.order, ]

  temp_1 <- temp_1[rownames.order, ]

  library(ComplexHeatmap)
  library(circlize)
  (t_max <- max(abs(min(df_heatmap_t, na.rm = T)), 0, abs(max(df_heatmap_t, na.rm = T))))
  newline_counts <- sapply(gregexpr("\n", rownames(df_heatmap_t)), function(x) ifelse(x[1] == -1, 0, length(x)))

  row_heights <- unit(newline_counts*3, "cm")


  ht2 <-
    Heatmap(as.matrix(df_heatmap_t[]),
            col = colorRamp2(c( 0, t_max), c( "white", "red")),
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            show_row_dend = F,
            row_names_side = 'left',
            row_dend_reorder = F,
            column_dend_reorder = F,
            column_names_rot = 45,
            column_names_centered = F,
            border = 1, row_title = "miR-TS",
            column_title_gp = gpar(fontsize = 14, fontface = "bold"),
            heatmap_legend_param = list(
              title="-log10(P)",
              labels=c("0", "5", "≥10"),
              legend_width = unit(2, "cm")),
            show_heatmap_legend = T,
            cell_fun = function(j, i, x, y, w, h, fill) {
              if(temp_1[i, j] =="***") {
                grid.text("***", x, y)
              } else if(temp_1[i, j] =="**") {
                grid.text("**", x, y)
              } else if(temp_1[i, j] =="*") {
                grid.text("*", x, y)
              } else if(temp_1[i, j] =="111") {
                grid.text("***", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else if(temp_1[i, j] =="11") {
                grid.text("**", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else if(temp_1[i, j] =="1") {
                grid.text("*", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else {
                grid.text("", x, y)
              }}
    )

  print(ht2)
}

#Fig. S9----
# load data for cutoff in 0.6, 0.8, 0.9:
Pooled.ht.cutoff0.6.forbioimarkers <- load_FigS_csvfile("Fig.S9-df_heatmap.biomarkers.using_TSI_cutoff0.6")
Pooled.ht.cutoff0.6.formiRTSscore <- load_FigS_csvfile("Fig.S9-df_heatmap.miRTSscore.using_TSI_cutoff0.6")
Pooled.ht.cutoff0.8.forbioimarkers <- load_FigS_csvfile("Fig.S9-df_heatmap.biomarkers.using_TSI_cutoff0.8")
Pooled.ht.cutoff0.8.formiRTSscore <- load_FigS_csvfile("Fig.S9-df_heatmap.miRTSscore.using_TSI_cutoff0.8")
Pooled.ht.cutoff0.9.forbioimarkers <- load_FigS_csvfile("Fig.S9-df_heatmap.biomarkers.using_TSI_cutoff0.9")
Pooled.ht.cutoff0.9.formiRTSscore <- load_FigS_csvfile("Fig.S9-df_heatmap.miRTSscore.using_TSI_cutoff0.9")

## -with clinical biomarkers----
for (plot.df.use in c(
  "Pooled.ht.cutoff0.6.forbioimarkers", "Pooled.ht.cutoff0.8.forbioimarkers", "Pooled.ht.cutoff0.9.forbioimarkers"
)) {
  df_heatmap <- get(plot.df.use)
  df_heatmap_t <- as.data.frame(
    reshape2::dcast(df_heatmap, y~x, value.var = "P_trans") %>%
      column_to_rownames("y")
  )
  temp_1 <- as.data.frame(
    reshape2::dcast(df_heatmap, y~x, value.var = "t_adj") %>%
      column_to_rownames("y")
  )
  temp_1[is.na(temp_1)] <- ""
  rownames(df_heatmap_t)

  df_heatmap_t <- df_heatmap_t[order_marker, ]

  temp_1 <- temp_1[order_marker, ]

  library(ComplexHeatmap)
  library(circlize)
  (t_max <- max(abs(min(df_heatmap_t, na.rm = T)), 0, abs(max(df_heatmap_t, na.rm = T))))
  newline_counts <- sapply(gregexpr("\n", rownames(df_heatmap_t)), function(x) ifelse(x[1] == -1, 0, length(x)))

  row_heights <- unit(newline_counts*3, "cm")


  ht1 <-
    Heatmap(as.matrix(df_heatmap_t[]),
            col = colorRamp2(c( 0, t_max), c( "white", "red")),
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            show_row_dend = F,
            row_names_side = 'left',
            row_dend_reorder = F,
            column_dend_reorder = F,
            column_names_rot = 45,
            column_names_centered = F,
            border = 1,
            column_title_gp = gpar(fontsize = 14, fontface = "bold"),
            heatmap_legend_param = list(
              title="-log10(P)",
              legend_width = unit(2, "cm")),
            show_heatmap_legend = T,
            cell_fun = function(j, i, x, y, w, h, fill) {
              if(temp_1[i, j] =="***") {
                grid.text("***", x, y)
              } else if(temp_1[i, j] =="**") {
                grid.text("**", x, y)
              } else if(temp_1[i, j] =="*") {
                grid.text("*", x, y)
              } else if(temp_1[i, j] =="111") {
                grid.text("***", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else if(temp_1[i, j] =="11") {
                grid.text("**", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else if(temp_1[i, j] =="1") {
                grid.text("*", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else {
                grid.text("", x, y)
              }}
    )
  print(ht1)

}

## -with miR-TS scores----

for (plot.df.use in c(
  "Pooled.ht.cutoff0.6.formiRTSscore", "Pooled.ht.cutoff0.8.formiRTSscore", "Pooled.ht.cutoff0.9.formiRTSscore"
)) {
  df_heatmap <- get(plot.df.use)
  df_heatmap$y <- gsub("_", " ", df_heatmap$y)
  df_heatmap_t <- as.data.frame(
    reshape2::dcast(df_heatmap, y~x, value.var = "P_trans") %>%
      column_to_rownames("y")
  )
  temp_1 <- as.data.frame(
    reshape2::dcast(df_heatmap, y~x, value.var = "t_adj") %>%
      column_to_rownames("y")
  )
  temp_1[is.na(temp_1)] <- ""
  rownames.order <- rev(rownames(df_heatmap_t))
  df_heatmap_t <- df_heatmap_t[rownames.order, ]

  temp_1 <- temp_1[rownames.order, ]

  library(ComplexHeatmap)
  library(circlize)
  (t_max <- max(abs(min(df_heatmap_t, na.rm = T)), 0, abs(max(df_heatmap_t, na.rm = T))))
  newline_counts <- sapply(gregexpr("\n", rownames(df_heatmap_t)), function(x) ifelse(x[1] == -1, 0, length(x)))

  row_heights <- unit(newline_counts*3, "cm")


  ht2 <-
    Heatmap(as.matrix(df_heatmap_t[]),
            col = colorRamp2(c( 0, t_max), c( "white", "red")),
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            show_row_dend = F,
            row_names_side = 'left',
            row_dend_reorder = F,
            column_dend_reorder = F,
            column_names_rot = 45,
            column_names_centered = F,
            border = 1, row_title = "miR-TS",
            column_title_gp = gpar(fontsize = 14, fontface = "bold"),
            heatmap_legend_param = list(
              title="-log10(P)",
              labels=c("0", "5", "≥10"),
              legend_width = unit(2, "cm")),
            show_heatmap_legend = T,
            cell_fun = function(j, i, x, y, w, h, fill) {
              if(temp_1[i, j] =="***") {
                grid.text("***", x, y)
              } else if(temp_1[i, j] =="**") {
                grid.text("**", x, y)
              } else if(temp_1[i, j] =="*") {
                grid.text("*", x, y)
              } else if(temp_1[i, j] =="111") {
                grid.text("***", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else if(temp_1[i, j] =="11") {
                grid.text("**", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else if(temp_1[i, j] =="1") {
                grid.text("*", x, y, vjust = 0.7,gp = gpar(cex=1.2))
              } else {
                grid.text("", x, y)
              }}
    )

  print(ht2)
}


#Fig. S10----

# The K-M curves are generated using the same code as in: Fig4_longitudinal analysis.R

#Fig. S11----
Smok_lung <- load_FigS_csvfile("Fig.S11-Smoking years and Lung miR-TS change")
x="smokquit2"
# x="ALL_smokyrs_imputed"
y="lung_change_rate"
ggplot(Smok_lung, aes(x = get(x), y = get(y))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 1) +  # Horizontal line at y = 0
  geom_jitter(width = 0.5, height = 0.01, size=1) +
  theme_minimal()+
  labs(
    # x="smoking years",
    x="quit smoking years",
    y="lung miR-TS\nchange rate"

  )+
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 12, family = "Arial", color="black"),  # font size and family
    axis.title = element_text(size = 14, face = "bold", color="black"),
    axis.text.y = element_text(size = 12, color="black"),
    axis.text.x = element_text(size = 8, color="black"),
    panel.grid.major = element_blank(),  # remove major gridlines
    panel.grid.minor = element_blank(),  # remove minor gridlines
    panel.background = element_blank(),  # remove background
    axis.ticks.x=element_line(),
    axis.line = element_line(colour = "black")  # add axis lines
  )  +
  # scale_y_continuous(limits = c(0, 2.9), breaks = c(0,1,2)) +  #trans = 'log10',  log-transformed y-axis and set limits
  # scale_x_continuous(limits = c(-10, 10),
  #                    breaks = c(-10, -5, 0, 5,  10),
  #                    labels = paste0(c(-10, -5, 0,5,  10), "y"))+
  geom_smooth(method = "lm",         # Add LOESS curve
              se = TRUE,                # Display confidence interval
              # span = 1,              # Adjust the smoothness of the curve (default is 0.75)
              color = "blue",           # Color of the LOESS curve
              fill = "lightblue")

#Fig. S12----
# No data used.
