
## Spearman's Rho and RMSE:####
Output_deconv <- deconv.res_____sig_mat_____sim_mat.null_tissueDepleted8_notUseRowSD.qc_TMMcpm.dtct0.v2$proportions[, 1:23]

# Make sure sample IDs line up
metadata_use <- metadata.null_tissueDepleted8_notUseRowSD
all.equal(rownames(Output_deconv), metadata_use$sample_id)

# Reorder metadata in the same order as Output_deconv rows
meta_match <- metadata_use[match(rownames(Output_deconv), metadata_use$sample_id), ]

# save it to the output object:
Output_deconv.all <- Output_deconv
meta_match.all <- meta_match

all.equal(rownames(Output_deconv.all), meta_match.all$sample_id)


# Tissues = columns of the deconvolution output
tissues <- colnames(Output_deconv.all)

# Initialize ground-truth matrix: all zeros
truth_mat <- matrix(
  0,
  nrow = nrow(Output_deconv.all),
  ncol = length(tissues),
  dimnames = list(rownames(Output_deconv.all), tissues)
)

# Fill in the true spike proportion for the spiked tissue in each sample
for (i in seq_len(nrow(Output_deconv.all))) {
  t_spike <- meta_match.all$tissue_spiked[i]
  if (!is.na(t_spike) && t_spike != "none" && t_spike %in% tissues) {
    truth_mat[i, t_spike] <- meta_match.all$spike_prop[i]
  }
}
# View(truth_mat)

rmse <- function(tru, est) {
  sqrt(mean((est - tru)^2, na.rm = TRUE))
}



cor_df <- data.frame(
  tissue        = tissues,
  cor_pearson   = NA_real_,
  cor_spearman  = NA_real_,
  R_sq = NA_real_,
  RMSE = NA_real_,
  stringsAsFactors = FALSE
)

for (j in seq_along(tissues)) {
  sample_id.sel <- meta_match.all %>% filter(tissue_spiked %in% c("none", tissues[j])) %>%
    pull(sample_id)

  est <- Output_deconv.all[ sample_id.sel, tissues[j]]
  tru <- truth_mat      [ sample_id.sel, tissues[j]]

  # If all true values are 0 (e.g., never spiked), correlation is undefined
  if (sum(tru != 0) > 1) {
    fit <- lm(tru ~ est)
    cor_df$cor_pearson[j]  <- cor(est, tru, method = "pearson")
    cor_df$cor_spearman[j] <- cor(est, tru, method = "spearman")
    cor_df$R_sq[j] <- summary(fit)$r.squared
    cor_df$RMSE[j] <- rmse(tru, est)
  }
}

cor_df$cor_spearman^2
cor_df$cor_pearson^2
cor_df$R_sq
cor_df$RMSE
cor_df.final_pipeline <- cor_df

# Aitchison Distance:####

aitchison <- function(p, q, pseudocount = 1e-12) {
  p <- Sum_to_1(p)
  q <- Sum_to_1(q)

  ## Handle zeros via pseudocount
  p[p == 0] <- pseudocount
  q[q == 0] <- pseudocount
  p <- p / sum(p)
  q <- q / sum(q)

  clr <- function(x) {
    g <- exp(mean(log(x)))     # geometric mean
    log(x / g)
  }

  clr_p <- clr(p)
  clr_q <- clr(q)

  sqrt(sum((clr_p - clr_q)^2))
}

trimedian <- function(x) {
  qs <- stats::quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  (qs[1] + 2 * qs[2] + qs[3]) / 4
}


#begins here:
metrics_by_tissue_spike.all <- data.frame()
Input_Signature = "sig_mat"
Input_mix = "sim_mat.null_tissueDepleted8_notUseRowSD"
Output_deconv <- get(paste0("deconv.res_____",
                            Input_Signature,
                            "_____",
                            Input_mix, ".qc_TMMcpm.dtct0", ".v2"))$proportions
dim(Output_deconv)
Output_deconv <- Output_deconv[, 1:(ncol(Output_deconv)-3)]

# Check sample tissue type:
metadata_use <- metadata.null_tissueDepleted8_notUseRowSD

# create into long format:
sim_long <- as.data.frame(Output_deconv) %>% rownames_to_column("sample_id") %>%
  reshape2::melt(
    id.var = "sample_id"
  ) %>% dplyr::rename(tissue = variable, estimated = value)
sim_long <- left_join(sim_long, metadata_use, by = "sample_id")%>%
  dplyr::rename(true = spike_prop) %>%
  select("sample_id" ,  "tissue"  ,  "tissue_spiked" ,  "true" ,  "estimated")


# reduce each sample to 2d, as I only spiked in one tissue type for one sample:
metrics_2d <- sim_long %>%
  group_by(sample_id, tissue_spiked) %>%
  summarise(
    p_true     = true[tissue == tissue_spiked][1],
    p_hat      = estimated[tissue == tissue_spiked][1],
    p_hat_rest = 1 - p_hat,
    aitchison_2d = aitchison(c(p_true, 1 - p_true),
                             c(p_hat,  p_hat_rest)),
    .groups = "drop"
  )


metrics_by_tissue_spike <- metrics_2d %>%
  group_by(tissue_spiked, p_true) %>%
  summarise(
    aitchison_median = trimedian(aitchison_2d),
    .groups = "drop"
  )
metrics_by_tissue_spike.all <- rbind(
  metrics_by_tissue_spike.all,
  metrics_by_tissue_spike %>%
    mutate(method = Input_Signature)
)


