#' @title Automated 2-Stage Multi-Environment Trial (MET) Weighting Bridge
#' @description Extracts Stage 1 spatial BLUEs, prediction error variances, and inverse-variance weighting 
#'   matrices (Smith et al. 2001, Piepho et al. 2012) to construct ready-to-run Stage 2 MET mixed models.
#' @param stage1_models A list of fitted single-trial models (\code{SpATS}, \code{lmerMod}, \code{mmer}, \code{lm}).
#' @param genotype_term Name of genotype factor in Stage 1 models (default: \code{"genotype"}).
#' @param trial_names Optional character vector of trial/location names.
#' @param method Weighting method: \code{"diagonal"} (\eqn{w_{ij} = 1/\text{SE}^2}) or \code{"normalized"} (scaled by trial harmonic mean variance).
#' @name stage1_to_stage2
NULL

#' @rdname stage1_to_stage2
#' @export
stage1_to_stage2 <- function(stage1_models, genotype_term = "genotype", trial_names = NULL, method = c("diagonal", "normalized")) {
  method <- match.arg(method)
  
  if (!is.list(stage1_models) || length(stage1_models) == 0) {
    stop("stage1_models must be a non-empty list of fitted single-trial models.")
  }
  
  n_trials <- length(stage1_models)
  t_names <- trial_names %||% names(stage1_models) %||% paste0("Location_", seq_len(n_trials))
  
  stage2_list <- list()
  trial_summary_list <- list()
  
  for (k in seq_len(n_trials)) {
    mod_k <- stage1_models[[k]]
    t_name <- t_names[k]
    
    # Extract BLUEs
    blues_k <- tryCatch({
      get_blues(mod_k, term = genotype_term)
    }, error = function(e) {
      # Fallback to coef
      cf <- coef(summary(mod_k))
      idx <- grep(genotype_term, rownames(cf), ignore.case = TRUE)
      if (length(idx) == 0) idx <- seq_len(nrow(cf))
      data.frame(
        Genotype = rownames(cf)[idx],
        BLUE = cf[idx, 1],
        SE = if (ncol(cf) >= 2) cf[idx, 2] else rep(1.0, length(idx)),
        stringsAsFactors = FALSE
      )
    })
    
    if (is.null(blues_k) || nrow(blues_k) == 0) {
      warning(sprintf("Could not extract BLUEs from trial '%s'. Skipping.", t_name))
      next
    }
    
    # Check SEs
    se_vec <- blues_k$SE
    if (all(is.na(se_vec)) || any(se_vec <= 0, na.rm = TRUE)) {
      se_vec[is.na(se_vec) | se_vec <= 0] <- mean(se_vec[!is.na(se_vec) & se_vec > 0], na.rm = TRUE) %||% 1.0
    }
    
    # Diagonal weights: w_i = 1 / SE_i^2
    var_vec <- se_vec^2
    raw_weights <- 1 / pmax(1e-6, var_vec)
    
    # Normalized weights by trial mean variance
    harm_mean_var <- 1 / mean(1 / pmax(1e-6, var_vec))
    norm_weights <- raw_weights * harm_mean_var
    
    final_weights <- if (method == "normalized") norm_weights else raw_weights
    
    # Heritability for trial
    h2_k <- tryCatch(get_heritability(mod_k)$Estimate[1], error = function(e) NA_real_)
    
    df_k <- data.frame(
      Trial = t_name,
      Genotype = as.character(blues_k$Genotype),
      BLUE = as.numeric(blues_k$BLUE),
      SE = as.numeric(se_vec),
      Variance = as.numeric(var_vec),
      Weight_Omega_Inv = as.numeric(final_weights),
      stringsAsFactors = FALSE
    )
    
    stage2_list[[length(stage2_list) + 1]] <- df_k
    
    trial_summary_list[[length(trial_summary_list) + 1]] <- data.frame(
      Trial = t_name,
      Genotypes_Count = nrow(df_k),
      Mean_BLUE = mean(df_k$BLUE, na.rm = TRUE),
      Mean_SE = mean(df_k$SE, na.rm = TRUE),
      Heritability_H2 = h2_k,
      Mean_Weight = mean(final_weights, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  
  full_stage2_data <- do.call(rbind, stage2_list)
  trial_summary <- do.call(rbind, trial_summary_list)
  
  # Suggested model syntax
  suggested_formula <- "BLUE ~ Trial + (1 | Genotype) + (1 | Genotype:Trial)"
  
  out <- list(
    data = full_stage2_data,
    trial_summary = trial_summary,
    total_trials = length(trial_summary_list),
    total_records = nrow(full_stage2_data),
    weight_method = method,
    suggested_lmer_code = sprintf("lmer(%s, weights = Weight_Omega_Inv, data = stage2_obj$data)", suggested_formula)
  )
  class(out) <- c("agri_stage2_bridge", "list")
  out
}

#' @export
print.agri_stage2_bridge <- function(x, ...) {
  cat("=== Stage 1 to Stage 2 MET Weighting Bridge (Smith & Piepho) ===\n")
  cat(sprintf("Total Environments/Trials: %d | Total Adjusted Records: %d\n", x$total_trials, x$total_records))
  cat(sprintf("Weighting Strategy: %s (Inverse Prediction Error Variance)\n\n", x$weight_method))
  cat("Trial Summary Overview:\n")
  print(x$trial_summary, row.names = FALSE)
  cat(sprintf("\nSuggested Stage 2 Mixed Model Execution:\n  %s\n", x$suggested_lmer_code))
  invisible(x)
}
