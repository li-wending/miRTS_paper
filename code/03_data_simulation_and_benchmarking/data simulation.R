# The signature matrix used to generate simulated data:
sig_mat <- miRTS::miRTS_signature_v1

# Useful functions for data sim:----

simulate_bulk_from_signature <- function(
    sig_mat,
    spike_props      = c(0.05, 0.10, 0.20),
    n_reps           = 50,          # per tissue × spike_prop
    n_null           = 50,          # number of null-only individuals *per tissue*
    global_phi       = 0.2,         # overdispersion
    lib_size_factor  = 1,           # optional global scaling of means
    tissues          = colnames(sig_mat),  # which tissues (columns) to use
    seed             = NULL         # set for reproducibility
) {
  if (!is.null(seed)) set.seed(seed)

  # Basic checks
  if (is.null(rownames(sig_mat))) {
    stop("sig_mat must have rownames (miRNA IDs).")
  }
  if (is.null(colnames(sig_mat))) {
    stop("sig_mat must have colnames (tissue names).")
  }
  if (!all(tissues %in% colnames(sig_mat))) {
    stop("Some specified 'tissues' are not columns of sig_mat.")
  }

  # Subset to requested tissues and apply library size scaling
  sig_sub <- sig_mat[, tissues, drop = FALSE]
  sig_sub <- sig_sub * lib_size_factor

  # 1. Compute global mean across tissues for each miRNA (used for dispersion)
  miRNA_means_global <- rowMeans(sig_sub)

  # 2. Compute per-miRNA dispersion (size parameter) if requested
  {
    # Global overdispersion: Var = mu + phi * mu^2
    # This corresponds to NB size = 1/phi (constant across miRNAs)
    size_vec <- rep(1 / global_phi, length(miRNA_means_global))
  }

  names(miRNA_means_global) <- rownames(sig_mat)
  names(size_vec)           <- rownames(sig_mat)

  sim_list <- list()
  metadata <- data.frame()
  id_counter <- 1L

  # Helper function: simulate one sample from a mean vector
  simulate_one <- function(mu_vec, size_vec) {
    # rnbinom requires non-negative mu; ensure numeric and non-negative
    mu_vec[mu_vec < 0] <- 0
    rnbinom(
      n    = length(mu_vec),
      mu   = mu_vec,
      size = size_vec
    )
  }

  # 3 & 4. For each tissue, construct tissue-specific null (tissue-depleted)
  #        and spike-in mixtures with that tissue
  for (t in tissues) {
    # All other tissues (depleted set)
    others <- setdiff(tissues, t)
    if (length(others) == 0L) {
      stop("Need at least 2 tissues to construct a tissue-depleted null for: ", t)
    }

    # Tissue-depleted null mean: mean of all *other* tissues (for this tissue t)
    null_mean_t <- rowMeans(sig_sub[, others, drop = FALSE])

    # 3. Simulate NULL individuals for this tissue (no spike-in; p = 0)
    if (n_null > 0) {
      for (r in seq_len(n_null)) {
        counts_vec <- simulate_one(null_mean_t, size_vec)
        sim_list[[id_counter]] <- counts_vec

        metadata <- rbind(
          metadata,
          data.frame(
            sample_id     = paste0("sample_", id_counter),
            tissue_spiked = t,   # tissue under test (p = 0 means null)
            spike_prop    = 0,
            replicate     = r,
            stringsAsFactors = FALSE
          )
        )
        id_counter <- id_counter + 1L
      }
    }

    # 4. Simulate spike-in individuals for this tissue at each spike_prop
    tissue_profile <- sig_sub[, t]

    for (p in spike_props) {
      for (r in seq_len(n_reps)) {
        # Expected mean for this spiked sample:
        # (1 - p) * (mean of other tissues) + p * (tissue column)
        mu_spike <- (1 - p) * null_mean_t + p * tissue_profile

        counts_vec <- simulate_one(mu_spike, size_vec)
        sim_list[[id_counter]] <- counts_vec

        metadata <- rbind(
          metadata,
          data.frame(
            sample_id     = paste0("sample_", id_counter),
            tissue_spiked = t,
            spike_prop    = p,
            replicate     = r,
            stringsAsFactors = FALSE
          )
        )
        id_counter <- id_counter + 1L
      }
    }
  }

  # 5. Combine into a bulk matrix: rows = miRNAs, columns = samples
  sim_mat <- do.call(cbind, sim_list)
  rownames(sim_mat) <- rownames(sig_mat)
  colnames(sim_mat) <- metadata$sample_id

  # Return both the simulated matrix and metadata, plus the means/size for reference
  list(
    sim_mat        = sim_mat,
    metadata       = metadata,
    miRNA_means    = miRNA_means_global,  # global mean across tissues (for reference)
    size_vec       = size_vec,
    params = list(
      spike_props     = spike_props,
      n_reps          = n_reps,
      n_null          = n_null,
      global_phi      = global_phi,
      lib_size_factor = lib_size_factor,
      tissues         = tissues
    )
  )
}

run_miRTS_score_parallel <- function(X, Y,
                                     perm = 0,
                                     QN   = FALSE,
                                     n_cores   = 1,
                                     chunk_size = 10
) {
  # X: signature matrix (features x tissues)
  # Y: mixture matrix (features x samples)

  stopifnot(all(rownames(X) %in% rownames(Y)))

  # Align rows
  Y <- Y[rownames(X), , drop = FALSE]

  # Split sample indices into chunks
  sample_idx <- seq_len(ncol(Y))
  chunks <- split(sample_idx, ceiling(sample_idx / chunk_size))

  message("run_miRTS_score_parallel: using ", n_cores, " cores and ",
          length(chunks), " chunks.")

  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  ## 1) Export objects that live in the *current* function environment
  parallel::clusterExport(
    cl,
    varlist = c("X", "Y", "miRTS_score", "perm", "QN", "CIBERSORT"),
    envir = environment()
  )

  ## 2) Export CoreAlg (and doPerm if you use permutations) from the global env
  ##    (adjust envir if CoreAlg lives somewhere else)
  parallel::clusterExport(
    cl,
    varlist = c("CoreAlg"),   # add "doPerm" here too if needed
    envir = .GlobalEnv
  )

  ## load libraries on workers
  parallel::clusterEvalQ(cl, {
    library(preprocessCore)
    library(parallel)
    library(e1071)
    NULL
  })

  # Run miRTS on each chunk of samples in parallel
  chunk_results <- parallel::parLapply(cl, chunks, function(idx) {
    Y_chunk <- Y[, idx, drop = FALSE]
    # IMPORTANT: cores = 1 here to avoid nested parallelism
    miRTS_score(signature_matrix  = X, Input_df = Y_chunk, Dtct_cutoff  = 0, QN = QN, Norm = "none", useAbsolute = FALSE)
  })

  # Combine results: assume miRTS_score returns a list with $proportions etc.
  prop_list      <- lapply(chunk_results, function(res) res$proportions)
  props_combined <- do.call(rbind, prop_list)

  out <- chunk_results[[1]]
  out$proportions <- props_combined
  out
}

# deconvolution function and create the dataframe:

miR_TS.createDF.doParallel <- function(Mix_matrix,
                                       Sig_matrix,
                                       N_cores_use   = 6,
                                       chunk_size_use = 10,
                                       do_plots      = TRUE) {
  {
    temp <- get(Mix_matrix)
    rownames(temp) <- gsub("-","_",rownames(temp))
    df_using.qc <- temp

    # Normally this is not required but in the rare cases where
    # certain sampled draws are skewed with extreme spike-in levels,
    # it is safer to re-normalize the simulated mix data:
    matrix <- as.data.frame(df_using.qc)
    matrix <- edgeR::DGEList(counts = matrix)
    matrix <- edgeR::calcNormFactors(matrix, method = "TMM")
    df_using.qc_TMMcpm <- edgeR::cpm(matrix)

    if (do_plots) {
      plot(colSums(df_using.qc_TMMcpm))
    }

    assign(paste0(Mix_matrix, ".qc_TMMcpm.dtct", Dtct_cutoff),
           envir = .GlobalEnv,
           as.data.frame(df_using.qc_TMMcpm))

    ## Deconvolution:
    version_mark        <- ".v2"

    Input_mix = paste0(Mix_matrix, ".qc_TMMcpm.dtct", Dtct_cutoff)
    print(Input_mix)
    Y <- get(Input_mix)
    X <- get(Sig_matrix)

    rownames(X) <- gsub("-","_",rownames(X))


    rownames(X) <- gsub("-","_",rownames(X))
    rownames(Y) <- gsub("-","_",rownames(Y))

    # Parallelized my_CIBERSORT
    deconv.res <- run_miRTS_score_parallel(
      X          = X,
      Y          = Y,
      perm       = 0,
      QN         = FALSE,
      n_cores    = N_cores_use,
      chunk_size = chunk_size_use
    )

    assign(
      paste0("deconv.res_____",
             Sig_matrix,
             "_____",
             Input_mix, version_mark),
      envir = .GlobalEnv,
      deconv.res
    )
  }

}





# Simulation:----
# at: 0.05, 0.1, 0.2, 0.40, 0.60, 0.80
sim_res.null_tissueDepleted8_notUseRowSD <- simulate_bulk_from_signature(
  sig_mat        = sig_mat,
  spike_props    = seq(0.05, 0.95, 0.05),
  n_reps         = 200,
  n_null         = 200,
  lib_size_factor = 1,
  seed           = 123
)
# load("./data/other_data/sim_res.null_tissueDepleted8_notUseRowSD.rds")
sim_mat.null_tissueDepleted8_notUseRowSD  <- sim_res.null_tissueDepleted8_notUseRowSD$sim_mat      # miRNAs × simulated individuals
metadata.null_tissueDepleted8_notUseRowSD <- sim_res.null_tissueDepleted8_notUseRowSD$metadata     # ground-truth info per individual
dim(sim_mat.null_tissueDepleted8_notUseRowSD)



# Deconvolute----
Input_df <- "sim_mat.null_tissueDepleted7_notUseRowSD"
Dtct_cutoff <- 0

## use parallel:

# The following will create two objects:
# 1) deconv.res_____sig_mat_____sim_mat.null_tissueDepleted8_notUseRowSD.qc_TMMcpm.dtct0.v2 - The deconvoluted results;
# 2) sim_mat.null_tissueDepleted8_notUseRowSD.qc_TMMcpm.dtct0 - The TMM-normalized mixture matrix

{
  date()
  miR_TS.createDF.doParallel(Mix_matrix = "sim_mat.null_tissueDepleted8_notUseRowSD",
                             Sig_matrix = "sig_mat",
                             N_cores_use   = 6,
                             chunk_size_use = 50,
                             do_plots      = FALSE);      beepr::beep(8)
  date()

}
# load("./data/other_data/deconv.res_____sig_mat_____sim_mat.null_tissueDepleted8_notUseRowSD.qc_TMMcpm.dtct0.v2.rds")

str(deconv.res_____sig_mat_____sim_mat.null_tissueDepleted8_notUseRowSD.qc_TMMcpm.dtct0.v2)
dim(sim_mat.null_tissueDepleted8_notUseRowSD.qc_TMMcpm.dtct0)
