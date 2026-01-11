#

load_Fig4_csvfile <- function(filename){
  temp <- read.csv(paste0("./data/Fig_results/",filename, ".csv"))
}


# Myocardial infarction----
## lineplot----
plt_df <- load_Fig4_csvfile("Fig.4_Myocardial infarction_lineplot")
plt <- ggplot(plt_df, aes(x = yr_dif, y = heart)) +
  # ggplot(test2, aes(x = yr_dif, y = heart)) +
  # geom_point(size=1) +
  # theme_pubclean()+
  theme_minimal() +  # use a publication-ready theme from ggpubr
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
  scale_x_continuous(breaks = c(-10, -3, 0,  10),labels = paste0(c(-10, -3, 0,  10), "y"))+
  labs(x="")

plt + geom_line(aes(x = yr_dif, y = fit), data=plt_df %>% filter(yr_dif<=0),size=2,
                color = '#1f77b4'
) + geom_line(aes(x = yr_dif, y = fit), data=plt_df %>% filter(yr_dif>=0),size=2,
              color = '#ff7f0e'
) + theme_void()+
  geom_ribbon( aes(ymin = lwr, ymax = upr), alpha = 0.2)

## boxplot----

plt_df <- load_Fig4_csvfile("Fig.4_Myocardial infarction_boxplot")

ggplot(plt_df ) +
  aes(x = yr_to_CHD.cat, y = log(heart+1), fill=yr_to_CHD.cat) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +  # avoid drawing outliers separately
  scale_y_continuous(limits = c(0, 3)) +  #trans = 'log10',  log-transformed y-axis and set limits
  scale_fill_manual(values = c("grey90", "#1f77b4", "#1f77b4", "#ff7f0e", "#ff7f0e", "#ff7f0e", "#ff7f0e")) +  # custom fill colors
  theme_minimal() +  # use a publication-ready theme from ggpubr
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 12, family = "Arial", color="black"),  # font size and family
    axis.title = element_text(size = 14, face = "bold", color="black"),
    axis.text = element_text(size = 12, color="black"),
    legend.position = "none",  # remove legend
    panel.grid.major = element_blank(),  # remove major gridlines
    panel.grid.minor = element_blank(),  # remove minor gridlines
    panel.background = element_blank(),  # remove background
    axis.line = element_line(colour = "black")  # add axis lines
  ) +
  labs(
    y = "heart"
  )

## K-M plot----
km_plot_df <- load_Fig4_csvfile("Fig.4_Myocardial infarction - KM_plot_data")


km_df_for_plot <- km_plot_df %>%
  transmute(
    time    = as.numeric(time),
    strata  = as.character(strata),
    surv    = as.numeric(1 - event),   # or if you already have surv column, use it directly
    n.risk  = as.numeric(n.risk),
    n.event = as.numeric(n.event),
    n.censor= as.numeric(n.censor)
  ) %>%
  as.data.frame()   # important: avoid tibble/list-col surprises

library(survminer)

surv_plot <- ggsurvplot_df(
  fit = km_df_for_plot,
  fun = "event",
  palette = c("#cda12d", "#D73027"),
  risk.table = FALSE,
  conf.int = FALSE,
  break.x.by = 5,
  ggtheme = theme_pubr(),
  legend = "none",
  legend.title = "                                        ",
  legend.labs = c("No", "Yes"),
  xlab = "Years",
  ylab = "Cumulative risk",
  font.x = c(20, "bold", "black"),
  font.y = c(20, "bold", "black"),
  font.tickslab = c(16, "plain", "black"),
  font.legend = c(14, "plain", "black")
)

print(surv_plot)

# Lung diseases----
## lineplot----
plt_df <- load_Fig4_csvfile("Fig.4_Lung diseases_lineplot")
plt <- ggplot(plt_df, aes(x = yr_df.1st, y = lung)) +
  theme_minimal() +  # use a publication-ready theme from ggpubr
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
  scale_x_continuous(breaks = c(-10, -3, 0,  10),labels = paste0(c(-10, -3, 0,  10), "y"))+
  labs(x="")
plt + geom_line(aes(x = yr_df.1st, y = fit), data=plt_df %>% filter(yr_df.1st<=0),size=2,
                color = '#1f77b4'
) + geom_line(aes(x = yr_df.1st, y = fit), data=plt_df %>% filter(yr_df.1st>=0),size=2,
              color = '#ff7f0e'
) + theme_void()+
  geom_ribbon( aes(ymin = lwr, ymax = upr), alpha = 0.2)

## boxplot----

plt_df <- load_Fig4_csvfile("Fig.4_Lung diseases_boxplot")

{  order_ploting <- c( "always healthy","in the future (>5y)","in the future (5~10y)","in the future (<5y)", "newly diagnosed (<1y)", "in the past (<5y)","in the past (5~10y)","in the past (>5y)"  )
  plt_df$Obstr_lung.cat <- factor(plt_df$Obstr_lung.cat, levels = order_ploting)
}

ggplot(plt_df) +
  aes(x = Obstr_lung.cat, y = log(lung+1), fill=Obstr_lung.cat) +
  xlab("") + #ylab(x) +
  theme_minimal()+
  geom_boxplot(outlier.shape = NA, width = 0.8) +  # avoid drawing outliers separately
  scale_y_continuous(limits = c(0, 5.5), breaks = c(seq(0,4, 2))) +  #trans = 'log10',  log-transformed y-axis and set limits
  scale_fill_manual(values = c("grey90", "#1f77b4", "#1f77b4", "#ff7f0e", "#ff7f0e", "#ff7f0e", "#ff7f0e")) +  # custom fill colors
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), #
    text = element_text(size = 12, family = "Arial", color="black"),  # font size and family
    axis.title = element_text(size = 14, face = "bold", color="black"),
    axis.text = element_text(size = 12, color="black"),
    legend.position = "none",  # remove legend
    panel.grid.major = element_blank(),  # remove major gridlines
    panel.grid.minor = element_blank(),  # remove minor gridlines
    panel.background = element_blank(),  # remove background
    axis.line = element_line(colour = "black")  # add axis lines
  ) +
  labs(
    y = "lung"
  )


## K-M plot----
km_plot_df <- load_Fig4_csvfile("Fig.4_Lung diseases - KM_plot_data")


km_df_for_plot <- km_plot_df %>%
  transmute(
    time    = as.numeric(time),
    strata  = as.character(strata),
    surv    = as.numeric(1 - event),   # or if you already have surv column, use it directly
    n.risk  = as.numeric(n.risk),
    n.event = as.numeric(n.event),
    n.censor= as.numeric(n.censor)
  ) %>%
  as.data.frame()   # important: avoid tibble/list-col surprises

library(survminer)

surv_plot <- ggsurvplot_df(
  fit = km_df_for_plot,
  fun = "event",
  palette = c("#cda12d", "#D73027"),
  risk.table = FALSE,
  conf.int = FALSE,
  break.x.by = 5,
  ggtheme = theme_pubr(),
  legend = "none",
  legend.title = "                                        ",
  legend.labs = c("No", "Yes"),
  xlab = "Years",
  ylab = "Cumulative risk",
  font.x = c(20, "bold", "black"),
  font.y = c(20, "bold", "black"),
  font.tickslab = c(16, "plain", "black"),
  font.legend = c(14, "plain", "black")
)

print(surv_plot)


# Cardiovascular accident(CVA)----
## lineplot----
plt_df <- load_Fig4_csvfile("Fig.4_Cardiovascular accident(CVA)_lineplot")
plt <- ggplot(plt_df, aes(x = yr_dif, y = brain)) +
  theme_minimal() +  # use a publication-ready theme from ggpubr
  theme(
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
  scale_x_continuous(breaks = c(-10, -3, 0,  10),labels = paste0(c(-10, -3, 0,  10), "y"))+
  labs(x="")
plt + geom_line(aes(x = yr_dif, y = fit), data=plt_df %>% filter(yr_dif<=0),size=2,
                color = '#1f77b4'
) + geom_line(aes(x = yr_dif, y = fit), data=plt_df %>% filter(yr_dif>=0),size=2,
              color = '#ff7f0e'
) + theme_void()+
  geom_ribbon( aes(ymin = lwr, ymax = upr), alpha = 0.2)

## boxplot----

plt_df <- load_Fig4_csvfile("Fig.4_Cardiovascular accident(CVA)_boxplot")

{  order_ploting <- c(
  "No stroke", "happen in 5to10yrs","happenin5yrs","acute","subacute 1to5yrs"  ,"subacute 5to10yrs"
)
  plt_df$yr_to_stroke.cat <- factor(plt_df$yr_to_stroke.cat, levels = order_ploting)
}

ggplot(plt_df ) +
  aes(x = yr_to_stroke.cat, y = log(brain+1), fill=yr_to_stroke.cat) +
  geom_boxplot(outlier.shape = NA, width = 0.8) +  # avoid drawing outliers separately
  scale_fill_manual(values = c("grey90", "#1f77b4","#1f77b4", "#ff7f0e", "#ff7f0e", "#ff7f0e", "#ff7f0e")) +  # custom fill colors
  theme_minimal() +  # use a publication-ready theme from ggpubr
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 12, family = "Arial", color="black"),  # font size and family
    axis.title = element_text(size = 14, face = "bold", color="black"),
    axis.text = element_text(size = 12, color="black"),
    legend.position = "none",  # remove legend
    panel.grid.major = element_blank(),  # remove major gridlines
    panel.grid.minor = element_blank(),  # remove minor gridlines
    panel.background = element_blank(),  # remove background
    axis.line = element_line(colour = "black")  # add axis lines
  ) +
  labs(
    y = "brain"
  )


## K-M plot----
km_plot_df <- load_Fig4_csvfile("Fig.4_Cardiovascular accident(CVA) - KM_plot_data")


km_df_for_plot <- km_plot_df %>%
  transmute(
    time    = as.numeric(time),
    strata  = as.character(strata),
    surv    = as.numeric(1 - event),   # or if you already have surv column, use it directly
    n.risk  = as.numeric(n.risk),
    n.event = as.numeric(n.event),
    n.censor= as.numeric(n.censor)
  ) %>%
  as.data.frame()   # important: avoid tibble/list-col surprises

library(survminer)

surv_plot <- ggsurvplot_df(
  fit = km_df_for_plot,
  fun = "event",
  palette = c("#cda12d", "#D73027"),
  risk.table = FALSE,
  conf.int = FALSE,
  break.x.by = 5,
  ggtheme = theme_pubr(),
  legend = "none",
  legend.title = "                                        ",
  legend.labs = c("No", "Yes"),
  xlab = "Years",
  ylab = "Cumulative risk",
  font.x = c(20, "bold", "black"),
  font.y = c(20, "bold", "black"),
  font.tickslab = c(16, "plain", "black"),
  font.legend = c(14, "plain", "black")
)

print(surv_plot)

# Cognitive impairment----
## lineplot----
plt_df <- load_Fig4_csvfile("Fig.4_Cognitive impairment_lineplot")
plt <- ggplot(plt_df, aes(x = year_to_1st_CogImp, y = brain)) +
  theme_minimal() +  # use a publication-ready theme from ggpubr
  theme(
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
  scale_x_continuous(breaks = c(-10, -3, 0,  10),labels = paste0(c(-10, -3, 0,  10), "y"))+
  labs(x="")
plt + geom_line(aes(x = year_to_1st_CogImp, y = fit), data=plt_df %>% filter(year_to_1st_CogImp<=0),size=2,
                color = '#1f77b4'
) + geom_line(aes(x = year_to_1st_CogImp, y = fit), data=plt_df %>% filter(year_to_1st_CogImp>=0),size=2,
              color = '#ff7f0e'
) + theme_void()+
  geom_ribbon( aes(ymin = lwr, ymax = upr), alpha = 0.2)


## boxplot----

plt_df <- load_Fig4_csvfile("Fig.4_Cognitive impairment_boxplot")

{  order_ploting <- c( "always healthy",  "in the future (5~10y)","in the future (<5y)", "newly diagnosed (<1y)", "in the past (<5y)","in the past (5~10y)","in the past (>10y)"  )
  plt_df$Cog_imp.cat2 <- factor(plt_df$Cog_imp.cat, levels = order_ploting)
}

ggplot(plt_df) +
  aes(x = Cog_imp.cat2, y = log(brain+1), fill=Cog_imp.cat2) +
  xlab("Time of diagnosis") + #ylab(x) +
  ggtitle("Cognitive impairment")+
  theme_minimal()+
  geom_boxplot(outlier.shape = NA, width = 0.8) +  # avoid drawing outliers separately
  scale_fill_manual(values = c("grey90", "#1f77b4", "#1f77b4", "#ff7f0e", "#ff7f0e", "#ff7f0e", "#ff7f0e")) +  # custom fill colors
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), #
    text = element_text(size = 12, family = "Arial", color="black"),  # font size and family
    axis.title = element_text(size = 14, face = "bold", color="black"),
    axis.text = element_text(size = 12, color="black"),
    legend.position = "none",  # remove legend
    panel.grid.major = element_blank(),  # remove major gridlines
    panel.grid.minor = element_blank(),  # remove minor gridlines
    panel.background = element_blank(),  # remove background
    axis.line = element_line(colour = "black")  # add axis lines
  ) +
  labs(
    y = "brain"
  )


## K-M plot----
km_plot_df <- load_Fig4_csvfile("Fig.4_Cognitive impairment - KM_plot_data")


km_df_for_plot <- km_plot_df %>%
  transmute(
    time    = as.numeric(time),
    strata  = as.character(strata),
    surv    = as.numeric(1 - event),   # or if you already have surv column, use it directly
    n.risk  = as.numeric(n.risk),
    n.event = as.numeric(n.event),
    n.censor= as.numeric(n.censor)
  ) %>%
  as.data.frame()   # important: avoid tibble/list-col surprises

library(survminer)

surv_plot <- ggsurvplot_df(
  fit = km_df_for_plot,
  fun = "event",
  palette = c("#cda12d", "#D73027"),
  risk.table = FALSE,
  conf.int = FALSE,
  break.x.by = 5,
  ggtheme = theme_pubr(),
  legend = "none",
  legend.title = "                                        ",
  legend.labs = c("No", "Yes"),
  xlab = "Years",
  ylab = "Cumulative risk",
  font.x = c(20, "bold", "black"),
  font.y = c(20, "bold", "black"),
  font.tickslab = c(16, "plain", "black"),
  font.legend = c(14, "plain", "black")
)

print(surv_plot)
