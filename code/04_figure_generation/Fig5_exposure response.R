#

load_Fig5_csvfile <- function(filename){
  temp <- read.csv(paste0("./data/Fig_results/",filename, ".csv"))
}



# Forest plot:----
Organ_marker.using <- c(
  "smoking years",
  "quit smoking years"
)
plt_df <- load_Fig5_csvfile("Fig.5_Forestplot")
plt_df$study <- factor(plt_df$study,
                       levels = rev(c("SY", "DFTJ", "NAS", "pooled")))

dummy <- as.data.frame(rbind(
  #   c(0.05,"NAS","pack years"),
  # c(-0.05,"NAS","pack years"),
  c(0.05,"NAS","smoking years"),
  c(-0.1,"NAS","smoking years"),

  c(0.1,"NAS","quit smoking years"),
  c(-0.1,"NAS","quit smoking years")

  # c(0.1,"NAS","smoking years vs. lung OS"),
  # c(-0.1,"NAS","smoking years vs. lung OS"),
  #
  # c(0.1,"NAS","quit smoking years vs. lung OS"),
  # c(-0.1,"NAS","quit smoking years vs. lung OS")
))
colnames(dummy) <- c("Est", "study", "Organ_marker")
dummy$Est <- as.numeric(dummy$Est)
dummy$Organ_marker <- factor(dummy$Organ_marker, levels = Organ_marker.using)

ggplot(plt_df %>%
         filter(Organ_marker %in% Organ_marker.using) %>%
         mutate(
           Organ_marker=factor(as.character(Organ_marker),levels = Organ_marker.using ))
       ,aes(y = study, x = Est))+
  geom_segment(aes(x = CI_l, xend = CI_h, color=study, yend = study))+
  geom_point(aes(size=weights/10, shape=study, color=study, alpha=1))+
  theme_bw() +
  scale_alpha_identity()+
  scale_size_area()+
  facet_wrap(~Organ_marker, # +Trans+df,
             ncol=5
             ,scales="free_x", labeller = label_wrap_gen(width = 22, multi_line = TRUE))+ #labeller(category = label_fn)
  scale_shape_manual(values = rev(c(15, 15, 15, 18))) +
  scale_color_manual(values = rev(c("#2d89c9", "#39b592", "#e6a23e", "#1e3135" )))+
  geom_vline(lty=2, aes(xintercept=ref_line), colour = 'red') +
  geom_blank(data=dummy) +
  theme(strip.text.x = element_text(size = 14),

        text = element_text(size = 14, family = "Arial", color="black"),  # font size and family
        axis.title = element_blank(),
        axis.text.y = element_text(size = 14, color="black"),
        axis.text.x = element_text(size = 12, color="black"),
        # panel.grid.major = element_blank(),  # remove major gridlines
        # panel.grid.minor = element_blank(),  # remove minor gridlines
        # panel.background = element_blank(),  # remove background
        axis.ticks.x=element_line(),
        axis.line = element_line(colour = "black"),  # add axis lines
        panel.spacing = unit(2, "lines"),
        legend.position = "none")

# box plot:----
plot_df <- load_Fig5_csvfile("Fig.5_boxplot")

if (length(unique(plot_df$sample_type)) == 2){
  color_str <- c("grey70", "#d4242a")
} else if (length(unique(plot_df$sample_type)) == 3){
  color_str <- c("grey70", "#FF4400", "#d4242a")
} else if (length(unique(plot_df$sample_type)) == 4){
  color_str <- c("grey70", "#FEA600", "#FF4400", "#d4242a")
} else if (length(unique(plot_df$sample_type)) == 6){
  color_str <- c("grey70", rev(c( "#FEC75E", "#FEA600", "#FF7A30", "#FF4400", "#d4242a")))
} else {
  color_str <- "Nothing"
}
plot_df$sample_type <- factor(
  plot_df$sample_type,
  levels = c(
    "healthy control", "day0", "day1", "day2", "day3", "day4+"
  )
)
ggplot(plot_df, aes(x = sample_type, y = y ))  +
  geom_boxplot(aes(fill = sample_type),
               width = 0.65,
               outlier.shape = NA) +
  labs(
    x = "", y = tissue_type) +
  scale_fill_manual(values = color_str) +
  theme_Publication()+
  theme(legend.position = "none",
        axis.title.y = element_text(size=14, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(
          # angle = 45, hjust = 1,
          size=12, color = "black"
        ),
        plot.margin = margin(rep(5,4))
  )

# Heatmap:----
df_heatmap <- load_Fig5_csvfile("Fig.5_heatmap")

df_heatmap_t <- as.data.frame(
  reshape2::dcast(df_heatmap, y ~ x, value.var = "t") %>%
    column_to_rownames("y")
)

temp_1 <- as.data.frame(
  reshape2::dcast(df_heatmap, y ~ x, value.var = "t_adj") %>%
    column_to_rownames("y")
)
temp_1[is.na(temp_1)] <- ""
order_marker <- rev(c(  "ord", "quit smoking years",  "smoking years"))
df_heatmap_t <- df_heatmap_t[order_marker, ]
temp_1 <- temp_1[order_marker, ]


library(ComplexHeatmap)
library(circlize)
(t_max <- max(abs(min(df_heatmap_t, na.rm = T)), 0, abs(max(df_heatmap_t, na.rm = T))))
newline_counts <- sapply(gregexpr("\n", rownames(df_heatmap_t)), function(x) ifelse(x[1] == -1, 0, length(x)))

row_heights <- unit(newline_counts*3, "cm")
ht1 <-
  Heatmap(as.matrix(df_heatmap_t[1:2,]),
          col = colorRamp2(c(-t_max, 0, t_max), c("blue", "white", "red")),
          cluster_rows = F,
          cluster_columns = FALSE,
          show_row_dend = F,
          row_names_side = 'left',
          row_dend_reorder = F,
          column_dend_reorder = F,
          # top_annotation=colAnn,
          column_names_rot = 45,
          column_names_centered = F,
          column_names_gp = gpar(fontsize = 20, fontface = "bold"),
          border = 1,
          column_title_gp = gpar(fontsize = 14, fontface = "bold"),
          heatmap_legend_param = list(
            title="t",
            legend_width = unit(1, "cm")),
          show_heatmap_legend = T,
          cell_fun = function(j, i, x, y, w, h, fill) {
            if(temp_1[i, j] =="***") {
              grid.text("***", x, y - unit(2, "mm"), gp = gpar(fontsize = 30, col = "black"))
            } else if(temp_1[i, j] =="**") {
              grid.text("**", x, y - unit(2, "mm"), gp = gpar(fontsize = 30, col = "black"))
            } else if(temp_1[i, j] =="*") {
              grid.text("*", x, y - unit(2, "mm"), gp = gpar(fontsize = 30, col = "black"))
            } else {
              grid.text("", x, y)
            }}
  )

temp_2 <- temp_1[3,]

ht2 <-
  Heatmap(as.matrix(df_heatmap_t[3, ]),
          col = colorRamp2(c(-t_max, 0, t_max), c("blue", "white", "red")),
          cluster_rows = F,
          cluster_columns = FALSE,
          show_row_dend = F,
          row_names_side = 'left',
          row_dend_reorder = F,
          column_dend_reorder = F,
          column_names_rot = 45,
          column_names_centered = F,
          column_names_gp = gpar(fontsize = 20, fontface = "bold"),
          border = 1,
          column_title_gp = gpar(fontsize = 14, fontface = "bold"),
          heatmap_legend_param = list(
            title="t",
            legend_width = unit(1, "cm")),
          show_heatmap_legend = F,
          cell_fun = function(j, i, x, y, w, h, fill) {
            if(temp_2[i, j] =="***") {
              grid.text("***", x, y - unit(2, "mm"), gp = gpar(fontsize = 30, col = "black"))
            } else if(temp_2[i, j] =="**") {
              grid.text("**", x, y - unit(2, "mm"), gp = gpar(fontsize = 30, col = "black"))
            } else if(temp_2[i, j] =="*") {
              grid.text("*", x, y - unit(2, "mm"), gp = gpar(fontsize = 30, col = "black"))
            } else {
              grid.text("", x, y)
            }}
  )


(ht <- ht1 %v% ht2)
