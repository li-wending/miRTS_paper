miRTS_score <- function (Input_df, signature_matrix = miRTS::miRTS_signature_v1,
          method = c("cibersort", "xCell2", "MCP-counter"), Dtct_cutoff = 0.1,
          Norm = c("TMM", "RLE", "none"), ts_miR.df = ts_miRNA, outlier_sample.cutoff = NA,
          useAbsolute = TRUE, Use_blood = TRUE, ...)
{
  method <- match.arg(method)
  Norm <- match.arg(Norm)
  rownames(Input_df) <- gsub("_", "-", rownames(Input_df))
  if (!all(grepl("miR-|let-", rownames(Input_df)))) {
    stop("Error: all row names must be miRNA names (e.g., hsa-miR-1-3p)!")
  }
  if (max(as.matrix(Input_df)) < 20) {
    warning("Max read count detected is less than 20! \n -The input data should be in non-log space.")
  }
  if (!is.na(outlier_sample.cutoff)) {
    libSizes <- colSums(as.matrix(Input_df))
    filtersamples <- function(filterParam, times = outlier_sample.cutoff) {
      samplesToRemove <- which(filterParam > median(filterParam) +
                                 times * mad(filterParam) | filterParam < median(filterParam) -
                                 times * mad(filterParam))
      samplesToRemove
    }
    samplesToRemove <- unique(unlist(lapply(list(libSizes = libSizes),
                                            filtersamples)))
    tt <- colnames(Input_df)[samplesToRemove]
    if (length(tt) > 0) {
      print("Samples removed:")
      print(paste(tt, collapse = ", "))
    }
    Input_df <- Input_df[, -samplesToRemove]
  }
  keep <- which(Matrix::rowSums(Input_df > 0) >= round(Dtct_cutoff *
                                                         ncol(Input_df)))
  Input_df = Input_df[keep, ]
  if (Norm == "none"){
    Y <- Input_df
    print("not normalizing Input_df!")
  } else{
    dge <- edgeR::DGEList(counts = as.matrix(Input_df))
    dge <- edgeR::calcNormFactors(dge, method = Norm)
    Y <- edgeR::cpm(dge)

  }
  X = signature_matrix
  rownames(X) <- gsub("-", "_", rownames(X))
  rownames(Y) <- gsub("-", "_", rownames(Y))
  Num_inSig <- sum(rownames(Y) %in% rownames(X))
  Prop_inSig <- round(Num_inSig/nrow(X), digits = 2)
  print(paste("Number of miRNAs using:", Num_inSig, "(", Prop_inSig,
              ")"))
  if (method == "cibersort") {
    if (!exists("CIBERSORT")) {
      stop("CIBERSORT function not detected! Run: \n1) ?CIBERSORT_download; \n2) navigate: Menu -> CS Archive -> CS Download -> Download CIBERSORT source code; \n3) source('CIBERSORT.R')")
    }
    X <- X[order(rownames(X)), , drop = FALSE]
    Y <- Y[order(rownames(Y)), , drop = FALSE]
    tmp_sig <- tempfile("miRTS_sig_", fileext = ".txt")
    tmp_mix <- tempfile("miRTS_mix_", fileext = ".txt")
    on.exit({
      if (file.exists(tmp_sig)) unlink(tmp_sig)
      if (file.exists(tmp_mix)) unlink(tmp_mix)
    }, add = TRUE)
    write.table(X, file = tmp_sig, sep = "\t", quote = FALSE,
                col.names = NA)
    write.table(Y, file = tmp_mix, sep = "\t", quote = FALSE,
                col.names = NA)
    CIBERSORT_output <- CIBERSORT(sig_matrix = tmp_sig, mixture_file = tmp_mix,
                                  perm = 0, QN = FALSE)
    deconv.res <- list(proportions = CIBERSORT_output, mix = Y,
                       signatures = X)
  }
  else if (method == "xCell2" & requireNamespace("xCell2",
                                                 quietly = T)) {
    min_overlap_miR = (length(intersect(rownames(Y), rownames(X))) -
                         1)/nrow(X)
    xcell2_results <- xCell2::xCell2Analysis(mix = Input_df,
                                             xcell2object = TA2_miR.raw_qc.257...DICE_demo.xCell2Ref,
                                             minSharedGenes = min_overlap_miR)
    deconv.res <- list(proportions = as.matrix(t(xcell2_results)),
                       mix = Y, signatures = X)
  }
  else if (method == "MCP-counter" & exists("compute_mcp_mirna_scores")) {
    ts_df = try(ts_miR.df[, c("miRNA", "enriched_organ")])
    if (inherits(ts_df, "try-error")) {
      stop("ts_miR.df is not correctly specified!")
    }
    mcp_scores <- compute_mcp_mirna_scores(expr_tpm = Input_df,
                                           ts_df = ts_df)
    deconv.res <- list(proportions = mcp_scores, mix = Y,
                       signatures = X)
  }
  if (method == "cibersort" & useAbsolute == TRUE) {
    X = X[rownames(X) %in% rownames(Y), ]
    print("Using absolute mode!")
    CIBERSORT.conv_ABS <- function(Input_mix, Input_Signature,
                                   Output_rel, Suppress_note = FALSE) {
      if (Suppress_note == FALSE) {
        print("Input mixture: rownames= miRs, colnames= sample names!")
        print("Input Signature: rownames= miRs, colnames= tissue types! -- note: these miRs should overlap with Input mixture!")
        print("Cibersort relative mode output to be converted: rownames= sample names, colnames= tissue types!")
      }
      Output_rel <- as.data.frame(Output_rel)
      drop_cols <- intersect(c("P-value", "Correlation",
                               "RMSE"), colnames(Output_rel))
      if (length(drop_cols) > 0) {
        Output_rel <- Output_rel[, setdiff(colnames(Output_rel),
                                           drop_cols), drop = FALSE]
      }
      Input_Signature = Input_Signature[, colnames(Output_rel)]
      if (!identical(colnames(Input_Signature), colnames(Output_rel)) ||
          !identical(rownames(Output_rel), colnames(Input_mix))) {
        stop("Error: Check sample names (mix vs CIBERSORT output) and tissue names (sig vs output).")
      }
      miR_overlap <- intersect(rownames(Input_Signature),
                               rownames(Input_mix))
      if (length(miR_overlap) == 0)
        stop("Error: no overlapping miRs!")
      sample_mean <- apply(Input_mix[miR_overlap, ], 2,
                           function(x) mean(x))
      sample_avg = ifelse(median(as.matrix(Input_mix)) ==
                            0, mean(as.matrix(Input_mix)), median(as.matrix(Input_mix)))
      scaling_factor <- as.numeric(sample_mean/sample_avg)
      sweep(Output_rel, MARGIN = 1, scaling_factor, `*`)
    }
    res.abs <- CIBERSORT.conv_ABS(Input_mix = Y, Input_Signature = X,
                                  Output_rel = deconv.res$proportions, Suppress_note = TRUE)
    orig <- as.data.frame(deconv.res$proportions)
    keep_diag <- intersect(c("P-value", "Correlation", "RMSE"),
                           colnames(orig))
    diag <- if (length(keep_diag) > 0)
      orig[, keep_diag, drop = FALSE]
    else NULL
    deconv.res$proportions <- cbind(res.abs, diag)
  }
  Transform_cibersort <- function(df) {
    df <- as.data.frame(df)
    nz <- colSums(as.matrix(df), na.rm = TRUE) != 0
    df <- df[, nz, drop = FALSE]
    drop_cols <- intersect(c("P-value", "Correlation", "RMSE"),
                           colnames(df))
    if (length(drop_cols) > 0) {
      df <- df[, setdiff(colnames(df), drop_cols), drop = FALSE]
    }
    min_tissue <- vapply(df, function(x) {
      x <- x[is.finite(x)]
      if (!length(x))
        return(NA_real_)
      ux <- sort(unique(x))
      if (length(ux) == 1)
        return(ux[1])
      if (ux[1] == 0)
        ux[2]
      else ux[1]
    }, numeric(1))
    for (i in seq_len(ncol(df))) {
      if (!is.na(min_tissue[i]))
        df[df[[i]] == 0, i] <- min_tissue[i]/2
    }
    df
  }
  if (method == "cibersort") {
    deconv.res$proportions <- Transform_cibersort(as.data.frame(deconv.res$proportions))
    if (Use_blood == TRUE) {
      blood_detect <- intersect(colnames(deconv.res$proportions),
                                miRTS::Tissues.blood_detect_17)
      deconv.res$proportions <- deconv.res$proportions[,
                                                       blood_detect]
    }
  }
  deconv.res
}
