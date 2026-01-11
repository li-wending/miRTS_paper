
load_Fig_sim_csvfile <- function(filename){
  temp <- read.csv(paste0("./data/Fig_results/Fig.S1-2 simulation/",filename, ".csv"))
  temp
}



# Fig.S1-Box plot , by tissues----

df_plot.all_tissues <- load_Fig_sim_csvfile("Fig.S1-Simulated data_miRTS final parameters")
ggplot(df_plot.all_tissues, aes(x = as.factor(true), y = estimated,
                                color = spike_color)) +
  geom_boxplot(outliers = F) +
  facet_wrap(~tissue_spiked,nrow = 3) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "black"),
                     name = "Spiked tissue?") +
  scale_x_discrete(
    labels = function(x) as.numeric(x)*100,
    breaks = as.character(seq(0,0.9, 0.1))
  ) +
  scale_y_continuous(
    labels = function(x) as.numeric(x)*100,
    breaks = seq(0,0.8, 0.2)
  ) +
  labs(
    x = "Ground truth spike-in proportion (%)",
    y = "Deconvoluted fraction (%)",
    # title = paste("CIBERSORT performance for", tissue_of_interest)
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    # axis.text.x = element_text(angle = 90),
    axis.text = element_text(color="black")
  )

# Fig.S2-correlation and RMSE----

cor_df.TA2_pipeline0       <- load_Fig_sim_csvfile("cor_df.TA2_pipeline0")
cor_df.RNAAtlas_pipeline.intsect       <- load_Fig_sim_csvfile("cor_df.RNAAtlas_pipeline.intsect")

cor_df.TA2_pipeline       <- load_Fig_sim_csvfile("cor_df.TA2_pipeline")
cor_df.TA2_TSI.published_pipeline       <- load_Fig_sim_csvfile("cor_df.TA2_TSI.published_pipeline")
cor_df.TA2.TSI0.9_pipeline       <- load_Fig_sim_csvfile("cor_df.TA2.TSI0.9_pipeline")
cor_df.TA2.CIBERSORT.sig_pipeline       <- load_Fig_sim_csvfile("cor_df.TA2.CIBERSORT.sig_pipeline")
cor_df.TA2.AutoGeneS.sig_pipeline       <- load_Fig_sim_csvfile("cor_df.TA2.AutoGeneS.sig_pipeline")

cor_df.final_pipeline       <- load_Fig_sim_csvfile("cor_df.final_pipeline")
cor_df.AutoGeneS_deconv       <- load_Fig_sim_csvfile("cor_df.AutoGeneS_deconv")
cor_df.xCell2_all       <- load_Fig_sim_csvfile("cor_df.xCell2_all")
cor_df.MCPcounter       <- load_Fig_sim_csvfile("cor_df.MCPcounter")

cor_df.RLE       <- load_Fig_sim_csvfile("cor_df.RLE")
cor_df.TMM       <- load_Fig_sim_csvfile("cor_df.TMM")

cor_df.TMM.2       <- load_Fig_sim_csvfile("cor_df.TMM.2")
cor_df.TMM.1       <- load_Fig_sim_csvfile("cor_df.TMM.1")
cor_df.TMM.5       <- load_Fig_sim_csvfile("cor_df.TMM.5")
cor_df.TMM.8       <- load_Fig_sim_csvfile("cor_df.TMM.8")

# Bind them together
cor_all <- bind_rows(
  # cor_df.TA2_pipeline0,
  # cor_df.RNAAtlas_pipeline.intsect,  # cor_df.RNAAtlas_TSI_0.9,

  # cor_df.TA2_pipeline,
  # cor_df.TA2_TSI.published_pipeline,
  # cor_df.TA2.TSI0.9_pipeline,
  # cor_df.TA2.CIBERSORT.sig_pipeline,
  # cor_df.TA2.AutoGeneS.sig_pipeline,
  #
  # cor_df.final_pipeline,
  # cor_df.AutoGeneS_deconv,
  # cor_df.xCell2_all, # cor_df.xCell2_sel257
  # cor_df.MCPcounter,
  #
  # cor_df.RLE,
  # cor_df.TMM,
  #
  cor_df.TMM.2 %>% dplyr::mutate(setting=as.character(setting)),
  cor_df.TMM.1 %>% dplyr::mutate(setting=as.character(setting)),
  cor_df.TMM.5 %>% dplyr::mutate(setting=as.character(setting)),
  cor_df.TMM.8 %>% dplyr::mutate(setting=as.character(setting))
)


# Make sure tissues and settings are nicely ordered
# (You can change the level order if you want a specific tissue ordering)
tissue_levels  <- sort(unique(cor_all$tissue))
# setting_levels <- c("miRTA2*", "RNAAtlas")
# setting_levels <- c("miRTA2.DE*",
#                     "miRTA2.source",
#                     "miRTA2.TSI_0.9",
#                     "CIBERSORT_sig",  "AutoGeneS_sig")
# setting_levels <- c("CIBERSORT*", "AutoGeneS", "xCell", "MCP-counter")

# setting_levels <- c("TMM*", "RLE")
setting_levels <- c("0.1", "0.2*", "0.5", "0.8")

cor_all <- cor_all %>% filter(tissue %in% tissue_list.use_17)%>%
  mutate(
    tissue = gsub("_", " ", tissue)
  )
cor_all <- cor_all %>%
  full_join(
    .,
    expand.grid(setting = unique(cor_all$setting), tissue = tissue_list.use_17.rm_)
  )
cor_all <- cor_all %>%
  mutate(
    tissue  = factor(tissue, levels = rev(sort(tissue_list.use_17.rm_))),
    setting = factor(setting, levels = setting_levels)
  )

# Two separate long data frames, one for Spearman, one for RMSE
df_spear <- cor_all %>%
  select(tissue, setting, value = cor_spearman)

df_rmse <- cor_all %>%
  select(tissue, setting, value = RMSE)



## Packages
library(dplyr)
library(ggplot2)
library(viridis)      # for viridis color palette (Nature-friendly)
library(stringr)

theme_nature_heatmap <- function(base_size = 8) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x  = element_text(
        angle = 45, hjust = 1, vjust = 1,colour = "black",
        size = base_size, face = "plain"
      ),
      axis.text.y  = element_text(colour = "black",
                                  size = base_size + 2, face = "plain"
      ),
      panel.grid   = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
      legend.title = element_text(size = base_size,colour = "black", face = "plain"),
      legend.text  = element_text(size = base_size - 1, colour = "black"),
      legend.key.height = unit(3, "mm"),
      legend.key.width  = unit(3, "mm"),
      panel.background = element_rect(fill = NA, colour = NA),
      plot.background  = element_rect(fill = NA, colour = NA),
      legend.background      = element_rect(fill = NA, colour = NA),
      legend.box.background  = element_rect(fill = NA, colour = NA),

      plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt")
    )
}


## Heatmap for Spearman correlation ----------------------------------

# If your Spearman correlations are mostly [0, 1], you can clamp the scale there.
# If you have negatives, set limits = c(min(df_spear$value), max(df_spear$value)).
spearman_limits <- c(0, 1) # range(df_spear$value, na.rm = TRUE) #c(0, max(df_spear$value, na.rm = TRUE)) #

p_spear <- ggplot(df_spear, aes(x = setting, y = tissue, fill = value)) +
  geom_tile(color = "white", linewidth = 0.3) +
  coord_fixed() +
  scale_fill_viridis(
    option = "magma",
    direction = 1,
    limits = spearman_limits,
    na.value = "grey40",
    oob = scales::squish,
    name = "Spearman's rho"
  ) +
  theme_nature_heatmap(base_size = 9)

p_spear


## Heatmap for RMSE ---------------------------------------------------

# For RMSE, lower = better. A light-to-dark sequential palette works well.
rmse_limits <- c(0.1, 0.4) # range(df_rmse$value, na.rm = TRUE)

p_rmse <- ggplot(df_rmse, aes(x = setting, y = tissue, fill = value)) +
  geom_tile(color = "white", linewidth = 0.3) +
  coord_fixed() +
  scale_fill_viridis(
    option = "viridis",
    direction = -1,                 # dark = low error, light = high error (or flip if you prefer)
    limits = rmse_limits,
    na.value = "grey40",
    oob = scales::squish,
    name = "RMSE"
  ) +
  theme_nature_heatmap(base_size = 9)

p_rmse


#  Fig.S2-Aitchison Distance_by methods----
metrics_collapsed <- load_Fig_sim_csvfile("Fig.S2-Aitchison Distance_by methods")
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(colorspace)
#

# across methods:####
{
  cols17  <- qualitative_hcl(17, palette = "Dark 3")
  tissue_list.use_17.rm_ <- gsub("_", " ", tissue_list.use_17)
  col.ColorShape <- c(qualitative_hcl(9, palette = "Dark 3"), qualitative_hcl(9, palette = "Dark 3")[1:8])
  names(col.ColorShape) <- sort(tissue_list.use_17.rm_)

  swatchplot(c(col.ColorShape))

  shape.ColorShape <- c(rep(16, 9), rep(17, 8))
  names(shape.ColorShape) <- sort(tissue_list.use_17.rm_)

}

# build the plot
# get the median of the miR-TS method:
temp <- as.numeric(metrics_collapsed %>% filter(methods=="miR-TS selected") %>% summarise(median(aitchison_median)))

metrics_collapsed.plt <- metrics_collapsed %>%
  mutate(
    aitchison_median = ifelse(aitchison_median>4, 4, aitchison_median)
  )

p <- ggplot(
  metrics_collapsed.plt ,
  aes(
    x = aitchison_median,
    y = methods,
    # color = tissue_spiked
  )
) +
  geom_vline(xintercept = temp , colour = "darkblue",linewidth = 0.2, linetype = "dashed")+
  geom_boxplot(linewidth = 0.5,median.linewidth =0.5,
               fill = c( "lightblue", rep("grey90", 12)),
               color = c("darkblue", rep("grey40", 12)),
               outlier.shape = NA,
               width = 0.4,
               alpha = 0.4
  ) +
  geom_jitter(
    aes(color = tissue_spiked, shape = tissue_spiked),
    height = 0.1,         # vertical jitter so points don't overlap too much
    width  = 0.02,           # keep x exact, no horizontal jitter
    size   = 1,
    alpha  = 0.6
  ) +
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10),limits = c(0,4.5),
                     breaks = c(0,1,2,4),labels = c(0,1,2,4)) +
  labs(
    x = "Aitchison Distance",
    y = NULL,
  ) +
  scale_color_manual(name   = "Tissue spiked", values = col.ColorShape) +
  scale_shape_manual(name   = "Tissue spiked", values = shape.ColorShape) +

  theme_bw(base_size = 14) +
  guides(color = guide_legend(nrow = 2),
         shape = guide_legend(nrow = 2)) + # Make legend horizontal
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(size = 12, colour = "black"),
    axis.text.y        = element_text(size = 13, colour = "black", face = c("bold", rep("plain", 13) )),
    axis.title         = element_text(face = "bold", colour = "black"),
  )
p
p1 <- p +
  theme(
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.text = element_text(margin = margin(l = 2, unit = "pt"), size = 10),
    legend.title = element_blank(),
    legend.key.size = grid::unit(0.1, "cm"),
    legend.spacing.x = grid::unit(0, "cm"),
    legend.spacing.y = grid::unit(-0.5, "cm"),
    legend.justification = "left",
    legend.box.just = "left",
    legend.box.spacing = unit(0, "pt"),
    legend.position = c(-0.2, 1.1),              # top-left *inside*
    plot.margin = margin(t = 40, r = 5.5, b = 5.5, l = 5.5),  # make room
    legend.margin        = margin(t = 0, r = 0, b = 0, l = 0),
    legend.box.margin    = margin(t = 0, r = 0, b = 0, l = 0),

  ) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE, title.position = "top"),
    shape = guide_legend(nrow = 2, byrow = TRUE, title.position = "top")  )


p1


#  Fig.S2-Aitchison Distance_heatmap----
metrics_diff_hm <- load_Fig_sim_csvfile("Fig.S2-Aitchison Distance_heatmap")


p_aitch <- ggplot(
  metrics_diff_hm,
  aes(x = tissue_spiked, y = methods, fill = diff_aitchison)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  # coord_fixed() +
  scale_fill_gradient2(
    low       = "darkgreen",
    mid       = "white",
    high      = "#A50026",
    midpoint  = 0,
    limits    = c(-0.31, 0.31),      # cap at -0.2 and 0.2
    oob       = scales::squish,    # squish values outside limits
    na.value  = "grey40",          # NAs shown as grey
    name      = "Δ Aitchison\n  Distance"
  ) +
  theme_nature_heatmap(base_size = 9) +
  labs(
    x = "Method",
    y = "Tissue spiked"
  ) %>%
  theme(
    plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5),  # make room

  )+
  guides(
    fill = guide_colorbar(
      # barwidth  = unit(4, "cm"),
      barheight = unit(4, "cm")
    )
  )

p_aitch
