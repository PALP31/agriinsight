#' @title Multi-Environment Trials (MET), GxE Interaction, and Stability Biplots
#' @description Comprehensive biometrical modeling of Genotype-by-Environment (GxE) interaction 
#'   including Additive Main Effects and Multiplicative Interaction (AMMI), Genotype plus Genotype-by-Environment (GGE) biplots, 
#'   AMMI Stability Value (ASV), Yield-Stability Index (YSI), Wricke's Ecovalence, Shukla's stability variance, 
#'   and Finlay-Wilkinson joint regression.
#' @param data A data.frame containing multi-environment trial data or a two-way matrix (Genotypes x Environments).
#' @param genotype Name of genotype column (default: \code{"Genotype"}).
#' @param environment Name of environment column (default: \code{"Environment"}).
#' @param yield Name of response/yield column (default: \code{"Yield"}).
#' @param rep Name of replication column (optional).
#' @param n_pc Number of principal interaction components to extract (default: \code{2}).
#' @name gxe_stability
NULL

#' @title Fit AMMI Model for Multi-Environment Trials
#' @description Fits the Additive Main Effects and Multiplicative Interaction (AMMI) model via SVD decomposition.
#' @param data A data.frame or two-way mean table.
#' @param genotype Genotype factor column name.
#' @param environment Environment factor column name.
#' @param yield Yield/trait numeric column name.
#' @param rep Replication column name (optional).
#' @param n_pc Number of IPCA components to retain (default: 2).
#' @return An S3 object of class \code{c("agri_ammi", "agri_gxe_stability")}.
#' @export
gxe_ammi <- function(data, genotype = "Genotype", environment = "Environment", yield = "Yield", rep = NULL, n_pc = 2) {
  # Build 2-way mean matrix
  if (is.matrix(data)) {
    y_mat <- data
    g_names <- rownames(data) %||% paste0("G", seq_len(nrow(data)))
    e_names <- colnames(data) %||% paste0("E", seq_len(ncol(data)))
    rownames(y_mat) <- g_names
    colnames(y_mat) <- e_names
    r_reps <- 1
  } else {
    df <- as.data.frame(data)
    if (!all(c(genotype, environment, yield) %in% colnames(df))) {
      stop(sprintf("Columns '%s', '%s', and '%s' must be present in data.", genotype, environment, yield))
    }
    
    # Calculate cell means
    fmla <- stats::as.formula(paste(yield, "~", genotype, "+", environment))
    mean_tab <- stats::aggregate(fmla, data = df, FUN = mean, na.rm = TRUE)
    
    # Reshape to matrix
    g_levels <- sort(unique(df[[genotype]]))
    e_levels <- sort(unique(df[[environment]]))
    
    y_mat <- matrix(NA_real_, nrow = length(g_levels), ncol = length(e_levels),
                    dimnames = list(g_levels, e_levels))
    
    for (i in seq_len(nrow(mean_tab))) {
      g_i <- as.character(mean_tab[[genotype]][i])
      e_j <- as.character(mean_tab[[environment]][i])
      y_mat[g_i, e_j] <- mean_tab[[yield]][i]
    }
    
    # Check for missing cells and impute with row+col mean
    if (any(is.na(y_mat))) {
      grand_m <- mean(y_mat, na.rm = TRUE)
      r_m <- rowMeans(y_mat, na.rm = TRUE)
      c_m <- colMeans(y_mat, na.rm = TRUE)
      for (r in seq_len(nrow(y_mat))) {
        for (c in seq_len(ncol(y_mat))) {
          if (is.na(y_mat[r, c])) {
            y_mat[r, c] <- r_m[r] + c_m[c] - grand_m
          }
        }
      }
    }
    
    r_reps <- if (!is.null(rep) && rep %in% colnames(df)) length(unique(df[[rep]])) else 1
  }
  
  G <- nrow(y_mat)
  E <- ncol(y_mat)
  max_pc <- min(G - 1, E - 1, n_pc)
  
  # Main effects
  mu <- mean(y_mat)
  alpha_i <- rowMeans(y_mat) - mu
  beta_j <- colMeans(y_mat) - mu
  
  # AMMI Interaction Residuals Matrix R
  R_mat <- y_mat - outer(rowMeans(y_mat), rep(1, E)) - outer(rep(1, G), colMeans(y_mat)) + mu
  
  # Singular Value Decomposition
  svd_res <- svd(R_mat)
  singular_values <- svd_res$d[seq_len(max_pc)]
  ss_ge_total <- sum(R_mat^2) * r_reps
  ss_ipc <- (singular_values^2) * r_reps
  exp_var_pct <- (singular_values^2) / sum(svd_res$d^2) * 100
  
  # Symmetrical biplot coordinates: G_score = U * sqrt(D), E_score = V * sqrt(D)
  D_sqrt <- diag(sqrt(singular_values), nrow = max_pc, ncol = max_pc)
  g_scores <- svd_res$u[, seq_len(max_pc), drop = FALSE] %*% D_sqrt
  e_scores <- svd_res$v[, seq_len(max_pc), drop = FALSE] %*% D_sqrt
  
  rownames(g_scores) <- rownames(y_mat)
  colnames(g_scores) <- paste0("IPCA", seq_len(max_pc))
  rownames(e_scores) <- colnames(y_mat)
  colnames(e_scores) <- paste0("IPCA", seq_len(max_pc))
  
  # Stability Metrics
  # 1. ASV (AMMI Stability Value)
  asv <- if (max_pc >= 2) {
    w_factor <- (exp_var_pct[1] / pmax(1e-6, exp_var_pct[2]))
    sqrt((w_factor * g_scores[, 1])^2 + (g_scores[, 2])^2)
  } else {
    abs(g_scores[, 1])
  }
  
  # 2. Wricke's Ecovalence (W_i^2)
  wricke_w2 <- rowSums(R_mat^2)
  
  # 3. Shukla's Stability Variance (Shukla 1972 Exact)
  if (G > 2 && E > 1) {
    ss_ge_mat <- sum(R_mat^2)
    shukla_var <- (G / ((G - 2) * (E - 1))) * wricke_w2 - (ss_ge_mat / ((G - 1) * (G - 2) * (E - 1)))
  } else {
    shukla_var <- rep(NA_real_, G)
  }
  
  # 4. Finlay-Wilkinson Joint Regression
  env_index <- colMeans(y_mat) - mu
  fw_slopes <- numeric(G)
  fw_ms_dev <- numeric(G)
  for (i in seq_len(G)) {
    lm_fw <- stats::lm(y_mat[i, ] ~ env_index)
    fw_slopes[i] <- coef(lm_fw)[2]
    fw_ms_dev[i] <- sum(residuals(lm_fw)^2) / pmax(1, E - 2)
  }
  
  # 5. Yield Stability Index (YSI)
  yield_ranks <- rank(-rowMeans(y_mat))
  asv_ranks <- rank(asv)
  ysi <- yield_ranks + asv_ranks
  
  stability_df <- data.frame(
    Genotype = rownames(y_mat),
    Mean_Yield = rowMeans(y_mat),
    Yield_Rank = yield_ranks,
    IPCA1 = g_scores[, 1],
    IPCA2 = if (max_pc >= 2) g_scores[, 2] else NA_real_,
    ASV = asv,
    ASV_Rank = asv_ranks,
    YSI = ysi,
    YSI_Rank = rank(ysi),
    Wricke_Ecovalence = wricke_w2,
    Shukla_Variance = pmax(0, shukla_var),
    Finlay_Wilkinson_b = fw_slopes,
    FW_MS_Dev = fw_ms_dev,
    stringsAsFactors = FALSE
  )
  stability_df <- stability_df[order(stability_df$YSI), ]
  
  pca_summary <- data.frame(
    Component = paste0("IPCA", seq_len(max_pc)),
    Singular_Value = singular_values,
    Sum_of_Squares = ss_ipc,
    Variance_Percent = exp_var_pct,
    Cumulative_Percent = cumsum(exp_var_pct),
    stringsAsFactors = FALSE
  )
  
  out <- list(
    model_type = "AMMI",
    grand_mean = mu,
    genotype_means = rowMeans(y_mat),
    environment_means = colMeans(y_mat),
    two_way_table = y_mat,
    interaction_residuals = R_mat,
    genotype_scores = as.data.frame(g_scores),
    environment_scores = as.data.frame(e_scores),
    pca_summary = pca_summary,
    stability = stability_df,
    call_info = list(genotype = genotype, environment = environment, yield = yield, n_pc = max_pc)
  )
  class(out) <- c("agri_ammi", "agri_gxe_stability", "list")
  out
}

#' @title Fit GGE Biplot Model for Multi-Environment Trials
#' @description Fits Genotype plus Genotype-by-Environment (GGE) biplot model (environment-centered SVD).
#' @param data A data.frame or two-way mean table.
#' @param genotype Genotype factor column name.
#' @param environment Environment factor column name.
#' @param yield Yield/trait numeric column name.
#' @param n_pc Number of PC components (default: 2).
#' @return An S3 object of class \code{c("agri_gge", "agri_gxe_stability")}.
#' @export
gxe_gge <- function(data, genotype = "Genotype", environment = "Environment", yield = "Yield", n_pc = 2) {
  # Build matrix
  df <- as.data.frame(data)
  fmla <- stats::as.formula(paste(yield, "~", genotype, "+", environment))
  mean_tab <- stats::aggregate(fmla, data = df, FUN = mean, na.rm = TRUE)
  
  g_levels <- sort(unique(df[[genotype]]))
  e_levels <- sort(unique(df[[environment]]))
  
  y_mat <- matrix(NA_real_, nrow = length(g_levels), ncol = length(e_levels),
                  dimnames = list(g_levels, e_levels))
  for (i in seq_len(nrow(mean_tab))) {
    y_mat[as.character(mean_tab[[genotype]][i]), as.character(mean_tab[[environment]][i])] <- mean_tab[[yield]][i]
  }
  
  # GGE Environment-Centered Matrix
  e_means <- colMeans(y_mat, na.rm = TRUE)
  R_gge <- sweep(y_mat, 2, e_means, "-")
  
  G <- nrow(y_mat)
  E <- ncol(y_mat)
  max_pc <- min(G - 1, E - 1, n_pc)
  
  svd_res <- svd(R_gge)
  singular_values <- svd_res$d[seq_len(max_pc)]
  exp_var_pct <- (singular_values^2) / sum(svd_res$d^2) * 100
  
  D_sqrt <- diag(sqrt(singular_values), nrow = max_pc, ncol = max_pc)
  g_scores <- svd_res$u[, seq_len(max_pc), drop = FALSE] %*% D_sqrt
  e_scores <- svd_res$v[, seq_len(max_pc), drop = FALSE] %*% D_sqrt
  
  rownames(g_scores) <- rownames(y_mat)
  colnames(g_scores) <- paste0("PC", seq_len(max_pc))
  rownames(e_scores) <- colnames(y_mat)
  colnames(e_scores) <- paste0("PC", seq_len(max_pc))
  
  pca_summary <- data.frame(
    Component = paste0("PC", seq_len(max_pc)),
    Singular_Value = singular_values,
    Variance_Percent = exp_var_pct,
    Cumulative_Percent = cumsum(exp_var_pct),
    stringsAsFactors = FALSE
  )
  
  out <- list(
    model_type = "GGE",
    genotype_means = rowMeans(y_mat),
    environment_means = e_means,
    two_way_table = y_mat,
    genotype_scores = as.data.frame(g_scores),
    environment_scores = as.data.frame(e_scores),
    pca_summary = pca_summary
  )
  class(out) <- c("agri_gge", "agri_gxe_stability", "list")
  out
}

#' @title Extract Biplot Coordinates and Explained Variance
#' @description Unified extraction of IPCA1, IPCA2 scores for AMMI and GGE biplots across packages.
#' @param model A fitted AMMI, GGE, or GxE model object (from agriinsight or agricolae).
#' @param ... Additional arguments.
#' @export
get_biplot_scores <- function(model, ...) {
  UseMethod("get_biplot_scores")
}

#' @export
get_biplot_scores.agri_gxe_stability <- function(model, ...) {
  list(
    genotypes = model$genotype_scores,
    environments = model$environment_scores,
    pca_summary = model$pca_summary,
    stability = model$stability %||% NULL
  )
}

#' @export
get_biplot_scores.AMMI <- function(model, ...) {
  # Support agricolae::AMMI objects
  g_scores <- model$biplot[model$biplot$type == "GEN", ]
  e_scores <- model$biplot[model$biplot$type == "ENV", ]
  
  list(
    genotypes = g_scores,
    environments = e_scores,
    pca_summary = model$analysis
  )
}

#' @export
get_biplot_scores.default <- function(model, ...) {
  if (is.matrix(model) || is.data.frame(model)) {
    fit <- gxe_ammi(model)
    return(get_biplot_scores(fit))
  }
  stop(sprintf("get_biplot_scores is not supported for class '%s'.", class(model)[1]))
}

#' @export
print.agri_gxe_stability <- function(x, ...) {
  cat(sprintf("=== %s Multi-Environment Stability Model ===\n", x$model_type))
  cat(sprintf("Genotypes: %d | Environments: %d\n", nrow(x$genotype_scores), nrow(x$environment_scores)))
  cat("\nInteraction Variance Explained:\n")
  print(x$pca_summary, row.names = FALSE)
  if (!is.null(x$stability)) {
    cat("\nTop 5 Most Stable & High Yielding Genotypes (by Yield Stability Index - YSI):\n")
    print(head(x$stability[, c("Genotype", "Mean_Yield", "ASV", "YSI", "Finlay_Wilkinson_b")], 5), row.names = FALSE)
  }
  invisible(x)
}
