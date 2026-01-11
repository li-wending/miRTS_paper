# Reproduce Figures - Fig.2
source("./code/01_data_prep/00_loadbasicfunctions.R")

load_Fig3_csvfile <- function(filename){
  temp <- read.csv(paste0(
    "./data/Fig_results/",
    "Fig.3_Public_Data_Validation-", filename, ".csv"
  ))
  cat.sample_type <- as.character(unique(temp$sample_type))
  if (sum(grepl("control|low$|no|Healthy", cat.sample_type, ignore.case = T))>0 ){
    temp$sample_type <-
      factor(temp$sample_type,
             levels=
               unique(c(cat.sample_type[grepl("control|low$|no|Healthy", cat.sample_type, ignore.case = T)],
                        sort(cat.sample_type))))
  }

  print(table(temp$sample_type))
  temp
}

# TBI####
tissue_type="brain"
plot_df <- load_Fig3_csvfile("TBI_miRTSscores")
cat.sample_type <- as.character(unique(plot_df$sample_type))
if (sum(grepl("control|low$|no|Healthy", cat.sample_type, ignore.case = T))>0 ){
  plot_df$sample_type <-
    factor(plot_df$sample_type,
           levels=
             unique(c(cat.sample_type[grepl("control|low$|no|Healthy", cat.sample_type, ignore.case = T)],
                      sort(cat.sample_type))))
}

if (length(unique(plot_df$sample_type)) == 2){
  color_str <- c("grey70", "#d4242a")
} else if (length(unique(plot_df$sample_type)) == 3){
  color_str <- c("grey70", "#FF4400", "#d4242a")
} else if (length(unique(plot_df$sample_type)) == 4){
  color_str <- c("grey70", "#d4242a", "#FF4400", "#FEA600")
} else if (length(unique(plot_df$sample_type)) == 5){
  color_str <- c("grey70", rev(c( "#FEC75E", "#FEA600", "#FF7A30", "#d4242a")))
} else {
  color_str <- "Nothing"
}

plot_df[,tissue_type] <- plot_df[,tissue_type] - min(plot_df[,tissue_type])
plot_df$y <- replace_outliers_with_na(log(plot_df[,tissue_type] + 1))
ggplot(plot_df, aes(x = sample_type, y = y ))  +
  geom_boxplot(aes(fill = sample_type),
               width = 0.65,
               outlier.shape = NA) +
  labs(x = "", y = tissue_type) +
  scale_fill_manual(values = color_str) +
  theme_Publication()+
  scale_y_continuous(breaks = seq(0, 0.2, 0.1)) +
  theme(legend.position = "none",
        axis.title.y = element_text(size=14, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(
          angle = 45, hjust = 1,
          size=12, color = "black"
        ),
        plot.margin = margin(rep(5,4))
  )

# fulminant_myocarditis####
tissue_type="heart"
plot_df <- load_Fig3_csvfile("fulminant_myocarditis")

cat.sample_type <- as.character(unique(plot_df$sample_type))
if (sum(grepl("control|low$|no|Healthy", cat.sample_type, ignore.case = T))>0 ){
  plot_df$sample_type <-
    factor(plot_df$sample_type,
           levels=
             unique(c(cat.sample_type[grepl("control|low$|no|Healthy", cat.sample_type, ignore.case = T)],
                      sort(cat.sample_type))))
}
if (length(unique(plot_df$sample_type)) == 2){
  color_str <- c("grey70", "#d4242a")
} else if (length(unique(plot_df$sample_type)) == 3){
  color_str <- c("grey70", "#FF4400", "#d4242a")
} else if (length(unique(plot_df$sample_type)) == 4){
  color_str <- c("grey70", "#FEA600", "#FF4400", "#d4242a")
} else {
  color_str <- "Nothing"
}

plot_df[,tissue_type] <- plot_df[,tissue_type] - min(plot_df[,tissue_type]); plot_df$y <- replace_outliers_with_na(log(plot_df[,tissue_type] + 1))
ggplot(plot_df
       , aes(x = sample_type, y = y ))  +
  geom_boxplot(aes(fill = sample_type),
               width = 0.65,
               outlier.shape = NA) +
  labs(
    x = "", y = tissue_type) +
  scale_fill_manual(values = color_str) +
  theme_Publication()+
  scale_x_discrete(labels = c(1:length(unique(plot_df$sample_type)))) +
  theme(legend.position = "none",
        axis.title.y = element_text(size=14, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(
          angle = 45, hjust = 1,
          size=12, color = "black"
        ),
        plot.margin = margin(rep(5,4))
  )


# adipose_inflammation####
tissue_type="adipocyte"
plot_df <- load_Fig3_csvfile("adipose_inflammation")

cat.sample_type <- as.character(unique(plot_df$sample_type))
if (sum(grepl("control|low$|no|Healthy", cat.sample_type, ignore.case = T))>0 ){
  plot_df$sample_type <-
    factor(plot_df$sample_type,
           levels=
             unique(c(cat.sample_type[grepl("control|low$|no|Healthy", cat.sample_type, ignore.case = T)],
                      sort(cat.sample_type))))
}

plot_df$sample_type <-
  factor(plot_df$sample_type,
         levels=c("low", "intermediate", "high"))

if (length(unique(plot_df$sample_type)) == 2){
  color_str <- c("grey70", "#d4242a")
} else if (length(unique(plot_df$sample_type)) == 3){
  color_str <- c("grey70", "#FF4400", "#d4242a")
} else if (length(unique(plot_df$sample_type)) == 4){
  color_str <- c("grey70", "#FEA600", "#FF4400", "#d4242a")
} else {
  color_str <- "Nothing"
}

plot_df[,tissue_type] <- plot_df[,tissue_type] - min(plot_df[,tissue_type]); plot_df$y <- replace_outliers_with_na(log(plot_df[,tissue_type] + 1))
ggplot(plot_df
       , aes(x = sample_type, y = y ))  +
  geom_boxplot(aes(fill = sample_type),
               width = 0.65,
               outlier.shape = NA) +
  labs(
    x = "", y = tissue_type) +
  scale_fill_manual(values = color_str) +
  theme_Publication()+
  scale_y_continuous(labels = number_format(accuracy = 0.1)) + # breaks = seq(0, 0.2, 0.1)
  # scale_x_discrete(labels = c(1:length(unique(plot_df$sample_type)))) +
  theme(legend.position = "none",
        axis.title.y = element_text(size=14, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(
          angle = 45, hjust = 1,
          size=12, color = "black"
        ),
        plot.margin = margin(rep(5,4))
  )

# dermatitis####
tissue_type="skin"
plot_df <- load_Fig3_csvfile("dermatitis")



if (length(unique(plot_df$sample_type)) == 2){
  color_str <- c("grey70", "#d4242a")
} else if (length(unique(plot_df$sample_type)) == 3){
  color_str <- c("grey70", "#d4242a", "#FEA600")
} else if (length(unique(plot_df$sample_type)) == 4){
  color_str <- c("grey70", "#FEA600", "#FF4400", "#d4242a")
} else {
  color_str <- "Nothing"
}

# pairwise_combinations <- pairwise_combinations[grepl("control", pairwise_combinations)]
# Create a Boxplot and Add p-value Indicators
plot_df[,tissue_type] <- plot_df[,tissue_type] - min(plot_df[,tissue_type])
plot_df$y <- replace_outliers_with_na(log(plot_df[,tissue_type] + 1))
ggplot(plot_df
       , aes(x = sample_type, y = y ))  +
  # geom_jitter(width = 0.2, color = "blue", size = 2)+
  geom_boxplot(aes(fill = sample_type),
               width = 0.65,
               outlier.shape = NA) +
  # stat_compare_means(
  #   comparisons = pairwise_combinations,
  #   # comparisons = list(c("NC", "AD"), c("NC", "VCID"), c("NC", "NVCID")), # Specify pairs to compare
  #   # comparisons = list(c("NC", "AD")), # Specify pairs to compare
  #   method = "t.test"                                         # Method for comparison
  # ) +
  labs(
    x = "", y = tissue_type) +
  scale_fill_manual(values = color_str) +
  # theme_minimal()+
  # theme_pubclean()+
  # ylim(0,0.25) +
  # scale_colour_Publication()+
  theme_Publication()+
  scale_y_continuous(breaks = seq(0, 0.2, 0.1), labels = number_format(accuracy = 0.1)) + #
  scale_x_discrete(labels = c(1:length(unique(plot_df$sample_type)))) +
  theme(legend.position = "none",
        axis.title.y = element_text(size=14, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(
          angle = 45, hjust = 1,
          size=12, color = "black"
        ),
        plot.margin = margin(rep(5,4))
  )

# IBD####
tissue_type="bowel"
plot_df <- load_Fig3_csvfile("IBD")


if (length(unique(plot_df$sample_type)) == 2){
  color_str <- c("grey70", "#d4242a")
} else if (length(unique(plot_df$sample_type)) == 3){
  color_str <- c("grey70", "#FF4400", "#d4242a")
} else if (length(unique(plot_df$sample_type)) == 4){
  color_str <- c("grey70", "#FEA600", "#FF4400", "#d4242a")
} else {
  color_str <- "Nothing"
}

plot_df[,tissue_type] <- plot_df[,tissue_type] - min(plot_df[,tissue_type]); plot_df$y <- replace_outliers_with_na(log(plot_df[,tissue_type] + 1))
ggplot(plot_df
       , aes(x = sample_type, y = y ))  +
  geom_boxplot(aes(fill = sample_type),
               width = 0.65,
               outlier.shape = NA) +
  labs(
    x = "", y = tissue_type) +
  scale_fill_manual(values = color_str) +
  theme_Publication()+
  scale_y_continuous(limits = c(0, 0.28), breaks = seq(0, 0.3, 0.1), labels = number_format(accuracy = 0.1)) + #
  scale_x_discrete(labels = c(1:length(unique(plot_df$sample_type)))) +
  theme(legend.position = "none",
        axis.title.y = element_text(size=14, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(
          angle = 45, hjust = 1,
          size=12, color = "black"
        ),
        plot.margin = margin(rep(5,4))
  )


# T2D_T2DKD####
tissue_type="kidney"
plot_df <- load_Fig3_csvfile("T2DKD")


if (length(unique(plot_df$sample_type)) == 2){
  color_str <- c("grey70", "#d4242a")
} else if (length(unique(plot_df$sample_type)) == 3){
  color_str <- c("grey70", "#FF4400", "#d4242a")
} else if (length(unique(plot_df$sample_type)) == 4){
  color_str <- c("grey70", "#FEA600", "#FF4400", "#d4242a")
} else {
  color_str <- "Nothing"
}
plot_df[,tissue_type] <- plot_df[,tissue_type] - min(plot_df[,tissue_type]); plot_df$y <- replace_outliers_with_na(log(plot_df[,tissue_type] + 1))
ggplot(plot_df
       , aes(x = sample_type, y = y ))  +
  # geom_jitter(width = 0.2, color = "blue", size = 2)+
  geom_boxplot(aes(fill = sample_type),
               width = 0.65,
               outlier.shape = NA) +
  # stat_compare_means(
  #   comparisons = pairwise_combinations,
  #   # comparisons = list(c("NC", "AD"), c("NC", "VCID"), c("NC", "NVCID")), # Specify pairs to compare
  #   # comparisons = list(c("NC", "AD")), # Specify pairs to compare
  #   method = "t.test"                                         # Method for comparison
  # ) +
  labs(
    x = "", y = tissue_type) +
  scale_fill_manual(values = color_str) +
  # theme_minimal()+
  # theme_pubclean()+
  # ylim(0,0.25) +
  # scale_colour_Publication()+
  theme_Publication()+
  scale_y_continuous(breaks = seq(0, 0.2, 0.05), labels = number_format(accuracy = 0.05)) + #
  scale_x_discrete(labels = c(1:length(unique(plot_df$sample_type)))) +
  theme(legend.position = "none",
        axis.title.y = element_text(size=14, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(
          angle = 45, hjust = 1,
          size=12, color = "black"
        ),
        plot.margin = margin(rep(5,4))
  )

# osteoporosis####
tissue_type="bone"
plot_df <- load_Fig3_csvfile("osteoporosis")

  if (length(unique(plot_df$sample_type)) == 2){
    color_str <- c("grey70", "#d4242a")
  } else if (length(unique(plot_df$sample_type)) == 3){
    color_str <- c("grey70", "#FF4400", "#d4242a")
  } else if (length(unique(plot_df$sample_type)) == 4){
    color_str <- c("grey70", "#FEA600", "#FF4400", "#d4242a")
  } else {
    color_str <- "Nothing"
  }
plot_df[,tissue_type] <- plot_df[,tissue_type] - min(plot_df[,tissue_type]); plot_df$y <- replace_outliers_with_na(log(plot_df[,tissue_type] + 1))
ggplot(plot_df %>%
         mutate(
           # sample_type=ifelse(sample_type=="VCID", "AD", sample_type),
           # sample_type=ifelse(sample_type=="NVCID", "NC", sample_type)
         )%>%
         filter(
           !is.na(y)
           # brain!=min(plot_df$brain)
         )
       , aes(x = sample_type, y = y ))  +
  # geom_jitter(width = 0.2, color = "blue", size = 2)+
  geom_boxplot(aes(fill = sample_type),
               width = 0.65,
               outlier.shape = NA) +
  # stat_compare_means(
  #   comparisons = pairwise_combinations,
  #   # comparisons = list(c("NC", "AD"), c("NC", "VCID"), c("NC", "NVCID")), # Specify pairs to compare
  #   # comparisons = list(c("NC", "AD")), # Specify pairs to compare
  #   method = "t.test"                                         # Method for comparison
  # ) +
  labs(
    x = "", y = tissue_type) +
  scale_fill_manual(values = color_str) +
  # theme_minimal()+
  # theme_pubclean()+
  # ylim(0,0.25) +
  # scale_colour_Publication()+
  theme_Publication()+
  scale_y_continuous( breaks = seq(0,2,1), labels = number_format(accuracy = 0.1)) + #
  scale_x_discrete(labels = c(1:length(unique(plot_df$sample_type)))) +
  theme(legend.position = "none",
        axis.title.y = element_text(size=14, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(
          angle = 45, hjust = 1,
          size=12, color = "black"
        ),
        plot.margin = margin(rep(5,4))
  )


# COVID_sev_165####
tissue_type="lung"
plot_df <- load_Fig3_csvfile("COVID_severity")


  if (length(unique(plot_df$sample_type)) == 2){
    color_str <- c("grey70", "#d4242a")
  } else if (length(unique(plot_df$sample_type)) == 3){
    color_str <- c("grey70", "#FF4400", "#d4242a")
  } else if (length(unique(plot_df$sample_type)) == 4){
    color_str <- c("grey70", "#FEA600", "#FF4400", "#d4242a")
  } else {
    color_str <- "Nothing"
  }

plot_df0 <- plot_df %>% mutate(y=get(tissue_type))
plot_df <- plot_df0 %>%
  ungroup() %>%
  mutate(
    # y starts at 0:
    # y = y - min(y)
  )
# plot_df$y <- replace_outliers_with_na(log(plot_df$y + 1))
ggplot(plot_df %>%
         mutate(
           sample_type=as.character(sample_type)
           # sample_type=ifelse(sample_type=="VCID", "AD", sample_type),
           # sample_type=ifelse(sample_type=="NVCID", "NC", sample_type)
         )
       , aes(x = sample_type, y = y ))  +
  # geom_jitter(width = 0.2, color = "blue", size = 2)+
  geom_boxplot(aes(fill = sample_type),
               width = 0.65,
               outlier.shape = NA) +
labs(
  x = "", y = tissue_type) +
  scale_fill_manual(values = color_str) +
  # theme_minimal()+
  # theme_pubclean()+
  # ylim(0,0.25) +
  # scale_colour_Publication()+
  theme_Publication()+
  # scale_y_continuous(limits = c(0.02, 0.1), breaks = seq(0, 0.1, 0.02), labels = number_format(accuracy = 0.02),
  #                    expand = c(0, 0.00)
  #                    ) + #
  scale_x_discrete(labels = c(1:length(unique(plot_df$sample_type)))) +
  theme(legend.position = "none",
        axis.title.y = element_text(size=14, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(
          angle = 45, hjust = 1,
          size=12, color = "black"
        ),
        plot.margin = margin(rep(5,4))
  )


# Liver_AReject####
tissue_type="liver"
plot_df <- load_Fig3_csvfile("liver_allograft_reject")


  if (length(unique(plot_df$sample_type)) == 2){
    color_str <- c("grey70", "#d4242a")
  } else if (length(unique(plot_df$sample_type)) == 3){
    color_str <- c("grey70", "#FF4400", "#d4242a")
  } else if (length(unique(plot_df$sample_type)) == 4){
    color_str <- c("grey70", "#FEA600", "#FF4400", "#d4242a")
  } else {
    color_str <- "Nothing"
  }

plot_df[,tissue_type] <- plot_df[,tissue_type] - min(plot_df[,tissue_type]); plot_df$y <- replace_outliers_with_na(log(plot_df[,tissue_type] + 1))
ggplot(plot_df
       , aes(x = sample_type, y = y ))  +
  # geom_jitter(width = 0.2, color = "blue", size = 2)+
  geom_boxplot(aes(fill = sample_type),
               width = 0.65,
               outlier.shape = NA) +
  # stat_compare_means(
  #   comparisons = pairwise_combinations,
  #   # comparisons = list(c("NC", "AD"), c("NC", "VCID"), c("NC", "NVCID")), # Specify pairs to compare
  #   # comparisons = list(c("NC", "AD")), # Specify pairs to compare
  #   method = "t.test"                                         # Method for comparison
  # ) +
  labs(
    x = "", y = tissue_type) +
  scale_fill_manual(values = color_str) +
  # theme_minimal()+
  # theme_pubclean()+
  # ylim(0,0.25) +
  # scale_colour_Publication()+
  theme_Publication()+
  scale_y_continuous(labels = number_format(accuracy = 0.1)) + # breaks = seq(0, 0.2, 0.1)
  scale_x_discrete(labels = c(1:length(unique(plot_df$sample_type)))) +
  theme(legend.position = "none",
        axis.title.y = element_text(size=14, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(
          angle = 45, hjust = 1,
          size=12, color = "black"
        ),
        plot.margin = margin(rep(5,4))
  )



# hepatitis_C####
tissue_type="liver"
plot_df <- load_Fig3_csvfile("hepatitis_C")


  if (length(unique(plot_df$sample_type)) == 2){
    color_str <- c("grey70", "#d4242a")
  } else if (length(unique(plot_df$sample_type)) == 3){
    color_str <- c("grey70", "#FF4400", "#d4242a")
  } else if (length(unique(plot_df$sample_type)) == 4){
    color_str <- c("grey70", "#FEA600", "#FF4400", "#d4242a")
  } else {
    color_str <- "Nothing"
  }
plot_df[,tissue_type] <- plot_df[,tissue_type] - min(plot_df[,tissue_type]); plot_df$y <- replace_outliers_with_na(log(plot_df[,tissue_type] + 1))
ggplot(plot_df
       , aes(x = sample_type, y = y ))  +
  # geom_jitter(width = 0.2, color = "blue", size = 2)+
  geom_boxplot(aes(fill = sample_type),
               width = 0.65,
               outlier.shape = NA) +
  # stat_compare_means(
  #   comparisons = pairwise_combinations,
  #   # comparisons = list(c("NC", "AD"), c("NC", "VCID"), c("NC", "NVCID")), # Specify pairs to compare
  #   # comparisons = list(c("NC", "AD")), # Specify pairs to compare
  #   method = "t.test"                                         # Method for comparison
  # ) +
  labs(
    x = "", y = tissue_type) +
  scale_fill_manual(values = color_str) +
  # theme_minimal()+
  # theme_pubclean()+
  # ylim(0,0.25) +
  # scale_colour_Publication()+
  theme_Publication()+
  scale_y_continuous(labels = number_format(accuracy = 0.1)) + # breaks = seq(0, 0.2, 0.1)
  scale_x_discrete(labels = c(1:length(unique(plot_df$sample_type)))) +
  theme(legend.position = "none",
        axis.title.y = element_text(size=14, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(
          angle = 45, hjust = 1,
          size=12, color = "black"
        ),
        plot.margin = margin(rep(5,4))
  )




