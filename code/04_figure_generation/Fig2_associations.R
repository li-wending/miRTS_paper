# Reproduce Figures - Fig.2
source("./code/01_data_prep/00_loadbasicfunctions.R")

# Fig. 2a:####
Forestplot.df <- read.csv("./data/Fig_results/Fig.2_Forest plot.csv")

sort(unique(Forestplot.df$Organ_marker))

Organ_marker.using <-  c(
  "Artery miR-TS vs. Arterial Stiffness", "Artery miR-TS vs. Mean Arterial Pressure", "Heart miR-TS vs. Coronary Heart Diseases", "Brain miR-TS vs. Cognitive Impairment",
  "Pancreas miR-TS vs. Fasting Blood Glucose", "Pancreas miR-TS vs. Hemoglobin A1C", "Pancreas miR-TS vs. Diabetes", "Kidney miR-TS vs. eGFR",
  "Lung miR-TS vs. Lung Diseases", "Lung miR-TS vs. Lung Function Decline","Liver miR-TS vs. ALT","Liver miR-TS vs. AST",
  "Lymph Node miR-TS vs. %Lymphocytes", "Adipocyte miR-TS vs. BMI", "Adipocyte miR-TS vs. Waist/Height", "Adipocyte miR-TS vs. %Body Fat"
)
Organ_marker.using[!(Organ_marker.using %in% unique(Forestplot.df$Organ_marker))]



Forestplot.df$Organ_marker <- factor(Forestplot.df$Organ_marker, levels = Organ_marker.using)


dummy <- as.data.frame(rbind(
  c(0.6,"NAS","artery miR-TS vs. arterial stiffness"),
  c(-0.6,"NAS","artery miR-TS vs. arterial stiffness"),
  c(0.5,"NAS","artery miR-TS vs. mean arterial pressure"),
  c(-0.5,"NAS","artery miR-TS vs. mean arterial pressure"),
  c(0.5,"NAS","liver miR-TS vs. ALT"),
  c(-0.5,"NAS","liver miR-TS vs. ALT"),
  c(0.5,"NAS","liver miR-TS vs. AST"),
  c(-0.5,"NAS","liver miR-TS vs. AST"),

  c(0.5,"NAS","heart miR-TS vs. coronary heart diseases"),
  c(-0.5,"NAS","heart miR-TS vs. coronary heart diseases"),
  c(1,"NAS","lung miR-TS vs. lung diseases"),
  c(-1,"NAS","lung miR-TS vs. lung diseases"),
  c(1.5,"NAS","lung miR-TS vs. lung function decline"),
  c(-1.5,"NAS","lung miR-TS vs. lung function decline"),
  c(0.5,"NAS","kidney miR-TS vs. eGFR"),
  c(-0.5,"NAS","kidney miR-TS vs. eGFR"),

  c(4,"NAS","pancreas miR-TS vs. diabetes"),
  c(-4,"NAS","pancreas miR-TS vs. diabetes"),
  c(0.5,"NAS","pancreas miR-TS vs. fasting blood glucose"),
  c(-0.5,"NAS","pancreas miR-TS vs. fasting blood glucose"),
  c(2,"NAS","pancreas miR-TS vs. hemoglobin A1C"),
  c(-2,"NAS","pancreas miR-TS vs. hemoglobin A1C"),
  c(4,"NAS","brain miR-TS vs. cognitive impairment"),
  c(-4,"NAS","brain miR-TS vs. cognitive impairment"),

  c(0.4,"NAS","adipocyte miR-TS vs. Waist/Height"),
  c(-0.4,"NAS","adipocyte miR-TS vs. Waist/Height"),
  c(0.3,"NAS","adipocyte miR-TS vs. BMI"),
  c(-0.3,"NAS","adipocyte miR-TS vs. BMI"),
  c(0.08,"NAS","adipocyte miR-TS vs. %body fat"),
  c(-0.08,"NAS","adipocyte miR-TS vs. %body fat"),
  c(0.3,"NAS","lymph node miR-TS vs. %lymphocytes"),
  c(-0.3,"NAS","lymph node miR-TS vs. %lymphocytes")



  # c(0.5,"NAS","heart miR-TS vs. arterial stiffness"),
  # c(-0.5,"NAS","heart miR-TS vs. arterial stiffness"),
  # c(0.3,"NAS","heart miR-TS vs. SBP"),
  # c(-0.3,"NAS","heart miR-TS vs. SBP"),
  # c(0.5,"NAS",""),
  # c(-0.5,"NAS",""),


  # c(0.3,"NAS","lung miR-TS vs. fev1/fvc%"),
  # c(-0.3,"NAS","lung miR-TS vs. fev1/fvc%"),
  # c(0.3,"NAS","pleurae miR-TS vs. fev1/fvc"),
  # c(-0.3,"NAS","pleurae miR-TS vs. fev1/fvc"),
  # c(1,"NAS","pleurae miR-TS vs. Lung diseases"),
  # c(-1,"NAS","pleurae miR-TS vs. Lung diseases"),
))
colnames(dummy) <- c("Est", "study", "Organ_marker")
dummy$Est <- as.numeric(dummy$Est)
dummy <- dummy %>%
  mutate(
    Organ_marker=stringr::str_to_title(Organ_marker),
    Organ_marker=gsub("Mir-Ts Vs","miR-TS vs",Organ_marker),
    Organ_marker=gsub("Alt","ALT",Organ_marker),
    Organ_marker=gsub("Ast","AST",Organ_marker),
    Organ_marker=gsub("Fev1/Fvc","FEV1/FVC",Organ_marker),
    Organ_marker=gsub("Egfr","eGFR",Organ_marker),
    Organ_marker=gsub("A1c","A1C",Organ_marker),
    Organ_marker=gsub("Bmi","BMI",Organ_marker),
    Organ_marker=gsub("","",Organ_marker),
    Organ_marker=gsub("","",Organ_marker),
    Organ_marker=gsub("","",Organ_marker),
  )
dummy$Organ_marker <- To_1st_upper(dummy$Organ_marker)
unique(dummy$Organ_marker)
dummy$Organ_marker <- factor(dummy$Organ_marker, levels = To_1st_upper(Organ_marker.using))

custom_labels <- c(
  "brain-cognitive decline" = "brain\ncognitive decline",
  "kidney-eGFR" =  "kidney\nestimated glomerular filtration rate",
  "liver-ALT" =  "liver\nALT",
  "lung-fev1/fvc" =  "lung\nfev1/fvc",
  "lung-lung diseases" =   "lung\nlung diseases",
  "lymph node- %lymphocytes" =  "lymph node\n %lymphocytes",
  "artery-arterial stiffness" =  "artery\narterial stiffness",
  "artery-mean arterial pressure" =   "artery\nmean arterial pressure",
  "pancreas-diabetes" =  "pancreas\ndiabetes"
)


Forestplot.df$Organ_marker <- factor(To_1st_upper(as.character(Forestplot.df$Organ_marker)),
                                     levels = To_1st_upper(Organ_marker.using))

Forestplot.df$study <- factor(
  Forestplot.df$study,
  levels = c( "pooled", "NAS","DFTJ", "SY")
)

# Forestplot.df_final.ALL <- Forestplot.df
p1 <-
  ggplot(Forestplot.df,aes(y = study, x = Est))+
  geom_segment(aes(x = CI_l, xend = CI_h, color=study, yend = study))+
  geom_point(aes(size=weights/10, shape=study, color=study, alpha=1))+
  theme_bw() +
  scale_alpha_identity()+
  scale_size_area()+
  facet_wrap(~Organ_marker,ncol=4,
             ,scales="free_x", labeller = label_wrap_gen(width = 28, multi_line = TRUE) #
  )+ #labeller(category = label_fn)
  scale_shape_manual(values = rev(c(15, 15, 15, 18))) +
  scale_color_manual(values = rev(c("#2d89c9", "#39b592", "#e6a23e", "#1e3135" )))+
  geom_vline(lty=2, aes(xintercept=ref_line), colour = 'red') +
  geom_blank(data=dummy) +
  theme(strip.text.x = element_text(size = 10),
        strip.background = element_rect(fill = "grey95", color = "black"),  # Set background color
        axis.title = element_blank(),
        legend.position = "none")

p1


ggsave("figure/Fig.1a.png", width = 8, height =6,units = "in",scale = 1,  dpi = 300)

# Fig. 2b:####
df_heatmap <- read.csv("./data/Fig_results/Fig.2-Assoc_Score_allOrgans.heatmap.csv")


order_marker <- rev(c(
  "%Body Fat",
  "Waist/Height", #
  "BMI", #
  "Arterial Stiffness",
  "Mean Arterial Pressure",
  "Cognitive Impairment",
  # "Cognitive Decline",
  # "mmse30", #
  "Coronary Heart Diseases",
  # "Arterial Stiffness (c)", #
  "eGFR", "ALT", "AST",
  # "airflow limitation", #
  # "FEV1/FVC", #
  "Lung Function Decline",
  "Lung Diseases",
  "%Lymphocytes",
  "Diabetes",
  "Hemoglobin A1C", #
  "Fasting Blood Glucose"
))

# df_heatmap <- Heatmap.AllOrgans %>%
#   filter(P!=0,
#          !grepl("muscle|testis|esophagus|thyroid|vein|pleurae", x),
#          y %in% order_marker,
#   ) %>%
#   # mutate(ICD_9=as.numeric(gsub("_benign|_malig", "", y))) %>%
#   # left_join(.,Cancer_class, by="ICD_9") %>%
#   mutate(
#     Padj=p.adjust(P, method="BH")) %>%
#   mutate(
#     t_adj=case_when(#P<0.001~ "***", P<0.01~ "**", P<0.05~ "*",
#       Padj<0.05~ "11",
#       Padj<0.2~ "1",
#       # Padj>=0.05 & P<0.05~ "1",
#       .default = "") )
# t_adj=case_when(P<0.001 & t>0 ~ "***", P<0.01 & t>0 ~ "**", P<0.05 & t>0 ~ "*", .default = "")) %>%
# filter(Cancer_type==Cancer_type_using, Trans==logTrans)
# df_heatmap$y <- gsub("(.{1,40})(\\s|$)", "\\1\n", df_heatmap$y)
{
  capitalize_first <- function(x) {
    paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
  }
  df_heatmap$x <- capitalize_first(df_heatmap$x)
  df_heatmap$y <- capitalize_first(df_heatmap$y)
  df_heatmap$y <- case_when(
    df_heatmap$y == "Ast" ~ "AST",
    df_heatmap$y == "Alt" ~ "ALT",
    df_heatmap$y == "Egfr" ~ "eGFR",
    df_heatmap$y == "Bmi" ~ "BMI",
    df_heatmap$y == "Hemoglobin a1c" ~ "Hemoglobin A1c",
    df_heatmap$y == "" ~ "",
    T ~ df_heatmap$y
  )
  sort(unique(df_heatmap$y))

  }
{
  df_heatmap_t <- as.data.frame(
    reshape2::dcast(df_heatmap, y~x, value.var = "t") %>%
      column_to_rownames("y")
    # rename( #UMFA=log_UMFA,
    #   `5-mTHF`=log_5MTHF,
    #   SAM=sam_nm, SAH=sah_nm, Homocysteine=log_HCys.wk0, Cysteine=Cys.wk0, Cystathionine=log_cystathionine_nm, Methionine=methionine_um, B12=log_pB12.wk0, Choline=choline_um, Betaine=betaine_um, Dimethylglycine=log_pDMG.wk0, TMAO=log_pTMAO.wk0)
  )
  temp_1 <- as.data.frame(
    reshape2::dcast(df_heatmap, y~x, value.var = "t_adj") %>%
      column_to_rownames("y")
  )
  temp_1[is.na(temp_1)] <- ""
  order_marker.new <- rev(c(
    "%body fat",
    "Waist/height",
    "BMI",
    "Arterial stiffness",
    "Mean arterial pressure",#
    "Cognitive impairment",
    # "Cognitive decline",
    "Coronary heart diseases",
    # "Arterial Stiffness (c)",
    "eGFR", "ALT", "AST",
    # "FEV1/FVC",
    "Lung function decline",
    "Lung diseases",
    "%lymphocytes",
    "Diabetes",
    "Hemoglobin A1c", #
    "Fasting blood glucose"
  ))
  df_heatmap_t <- df_heatmap_t[order_marker.new, ]

  temp_1 <- temp_1[order_marker.new, ]
  rownames(df_heatmap_t) <- order_marker.new
  rownames(temp_1) <- order_marker.new

  library(ComplexHeatmap)
  library(circlize)
  # dup_names <- gsub(" :.*", "", rownames(df_heatmap_t))
  (t_max <- max(abs(min(df_heatmap_t, na.rm = T)), 0, abs(max(df_heatmap_t, na.rm = T))))
  # Count the number of newline characters in each element
  newline_counts <- sapply(gregexpr("\n", rownames(df_heatmap_t)), function(x) ifelse(x[1] == -1, 0, length(x)))

  row_heights <- unit(newline_counts*3, "cm")
  # colnames(df_heatmap_t) <- stringr::str_to_title(colnames(df_heatmap_t))

  ht <-
    Heatmap(as.matrix(df_heatmap_t[]),
            col = colorRamp2(c(-t_max, 0, t_max), c("blue", "white", "red")),
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            show_row_dend = F,
            # show_column_dend = T,
            # show_row_names = T,
            row_names_side = 'left',
            row_dend_reorder = F,
            column_dend_reorder = F,
            # top_annotation=colAnn,
            column_names_rot = 45,
            column_names_centered = F,
            # row_names_gp = gpar(col = ifelse(dup_names %in% dup_names[duplicated(dup_names)], "red", "black")),
            # left_annotation = rowAnn,
            # column_km = 2,
            border = 1,
            # column_title = "Arsenic exposure (bAs)", column_title_side = "bottom",
            column_names_gp = gpar(fontsize = 12),
            row_names_gp = gpar(fontsize = 12),
            column_title_gp = gpar(fontsize = 14, fontface = "bold"),
            heatmap_legend_param = list(
              title="t",
              legend_width = unit(2, "cm")),
            show_heatmap_legend = T,
            # cell_fun = function(j, i, x, y, w, h, fill) {
            #   if(temp_1[i, j] =="***") {
            #     grid.text("✱✱✱", x, y)
            #   } else if(temp_1[i, j] =="**") {
            #     grid.text("✱✱", x, y)
            #   } else if(temp_1[i, j] =="*") {
            #     grid.text("✱", x, y)
            #   } else {
            #     grid.text("", x, y)
            #   }}
            # height= row_heights,
            cell_fun = function(j, i, x, y, w, h, fill) {
              if(temp_1[i, j] =="***") {
                grid.text("***", x, y)
              } else if(temp_1[i, j] =="**") {
                grid.text("**", x, y)
              } else if(temp_1[i, j] =="*") {
                grid.text("*", x, y)
              } else if(temp_1[i, j] =="11") {
                grid.text("**", x, y, vjust = 0.7,gp = gpar(cex=1.5))
              } else if(temp_1[i, j] =="1") {
                grid.text("*", x, y, vjust = 0.7,gp = gpar(cex=1.5))
              } else {
                grid.text("", x, y)
              }}
    )

}
ht

png(paste("Heatmap-Score_AllOrgan_FDR_0.2 ALL_Cap", "2025Apr18 .png", sep = "_"), # As and OCM&FA metabolites__batchAdjusted.baseline_PBO.png
    width=8,height=4,units="in",res=300)

draw(ht, padding =  unit(c(2, 35, 2, 2), "mm"))

dev.off()
