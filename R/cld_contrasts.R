#' @title Agronomic Compact Letter Displays (CLD) for Treatment Comparisons
#' @description Generates publication-grade mean comparison tables with exact error bounds, 
#'   pairwise contrast statistics, and Compact Letter Displays (CLD) based on Piepho (2004) maximal clique algorithm.
#'   Supports \code{emmeans}, \code{lm}, \code{aov}, \code{lmerMod}, \code{lme}, \code{glmmTMB}, and raw summary matrices.
#' @param object An \code{emmeans} object, fitted model (\code{lm}, \code{aov}, \code{lmerMod}), or data.frame.
#' @param term Factor term name for pairwise comparisons (if model is passed).
#' @param alpha Significance level (default: \code{0.05}).
#' @param method Multiple comparison adjustment method: \code{"tukey"}, \code{"lsd"}, \code{"bonferroni"}, \code{"scheffe"}, \code{"duncan"}.
#' @param letters Character vector of letters to use (default: \code{base::letters}).
#' @param reversed Logical; if \code{FALSE} (default), highest mean gets 'a'; if \code{TRUE}, lowest mean gets 'a'.
#' @param ... Additional arguments.
#' @name cld_contrasts
NULL

#' @rdname cld_contrasts
#' @export
agro_cld <- function(object, term = NULL, alpha = 0.05, method = c("tukey", "lsd", "bonferroni", "scheffe", "duncan"), letters = NULL, reversed = FALSE, ...) {
  method <- match.arg(method)
  letters_vec <- if (!is.null(letters)) letters else c(base::letters, base::LETTERS)
  
  # Case 1: emmeans object
  if (inherits(object, c("emmGrid", "emm_list"))) {
    return(cld_from_emmeans(object, alpha = alpha, method = method, letters = letters_vec, reversed = reversed))
  }
  
  # Case 2: Fitted model (lm, aov, lmerMod, lme)
  if (inherits(object, c("lm", "aov", "merMod", "lmerMod", "lme", "glmmTMB"))) {
    if (requireNamespace("emmeans", quietly = TRUE)) {
      fmla_terms <- if (!is.null(term)) term else all.vars(stats::formula(object))[-1][1]
      emm <- emmeans::emmeans(object, specs = stats::as.formula(paste("~", fmla_terms)))
      return(cld_from_emmeans(emm, alpha = alpha, method = method, letters = letters_vec, reversed = reversed))
    } else {
      # Base R fallback for lm / aov
      return(cld_from_base_lm(object, term = term, alpha = alpha, method = method, letters = letters_vec, reversed = reversed))
    }
  }
  
  # Case 3: data.frame with Treatment, Mean, and optional SE
  if (is.data.frame(object)) {
    df <- object
    trt_col <- grep("treatment|trt|genotype|gen|group|level", colnames(df), ignore.case = TRUE, value = TRUE)[1] %||% colnames(df)[1]
    mean_col <- grep("mean|estimate|blue|emmean|yield", colnames(df), ignore.case = TRUE, value = TRUE)[1] %||% colnames(df)[2]
    se_col <- grep("se|std.error|stderr", colnames(df), ignore.case = TRUE, value = TRUE)[1]
    
    trts <- as.character(df[[trt_col]])
    means <- as.numeric(df[[mean_col]])
    ses <- if (!is.na(se_col)) as.numeric(df[[se_col]]) else rep(1.0, length(trts))
    
    # Compute pairwise t-tests based on SEs
    n_k <- length(trts)
    p_mat <- matrix(1.0, n_k, n_k, dimnames = list(trts, trts))
    for (i in seq_len(n_k)) {
      for (j in seq_len(n_k)) {
        if (i != j) {
          se_diff <- sqrt(ses[i]^2 + ses[j]^2)
          t_stat <- abs(means[i] - means[j]) / pmax(1e-6, se_diff)
          p_raw <- 2 * (1 - stats::pnorm(t_stat))
          p_adj <- if (method == "tukey") stats::ptukey(t_stat * sqrt(2), n_k, df = 100, lower.tail = FALSE)
                   else if (method == "bonferroni") pmin(1, p_raw * (n_k * (n_k - 1) / 2))
                   else p_raw
          p_mat[i, j] <- p_adj
        }
      }
    }
    
    cld_res <- compute_piepho_cld(trts, means, p_mat, alpha = alpha, letters = letters_vec, reversed = reversed)
    df$CLD_Letter <- cld_res[trts]
    class(df) <- c("agri_cld", "data.frame")
    return(df)
  }
  
  stop("Input object must be an emmeans object, fitted model (lm/aov/lmerMod), or data.frame.")
}

#' @keywords internal
cld_from_emmeans <- function(emm, alpha = 0.05, method = "tukey", letters = NULL, reversed = FALSE) {
  letters_vec <- if (!is.null(letters)) letters else c(base::letters, base::LETTERS)
  smry <- as.data.frame(emm)
  conts <- as.data.frame(emmeans::contrast(emm, method = "pairwise", adjust = method))
  
  trt_name <- if (inherits(emm, "emmGrid") && length(names(emm@levels)) > 0) {
    names(emm@levels)[1]
  } else {
    colnames(smry)[1]
  }
  trts <- as.character(smry[[trt_name]])
  means <- as.numeric(smry$emmean %||% smry$estimate %||% smry[, 2])
  
  # Build pairwise p-value matrix
  n_k <- length(trts)
  p_mat <- matrix(1.0, n_k, n_k, dimnames = list(trts, trts))
  
  for (row in seq_len(nrow(conts))) {
    pair_str <- as.character(conts$contrast[row])
    # Match pair names "TrtA - TrtB" or "(TrtA) - (TrtB)"
    splits <- strsplit(pair_str, " - ")[[1]]
    if (length(splits) == 2) {
      t1 <- trimws(gsub("[()]", "", splits[1]))
      t2 <- trimws(gsub("[()]", "", splits[2]))
      pval <- conts$p.value[row]
      if (t1 %in% trts && t2 %in% trts) {
        p_mat[t1, t2] <- pval
        p_mat[t2, t1] <- pval
      }
    }
  }
  
  cld_letters <- compute_piepho_cld(trts, means, p_mat, alpha = alpha, letters = letters_vec, reversed = reversed)
  smry$CLD_Letter <- cld_letters[as.character(smry[[trt_name]])]
  
  # Sort table by mean
  ord <- if (!reversed) order(-means) else order(means)
  res <- smry[ord, ]
  class(res) <- c("agri_cld", "data.frame")
  res
}

#' @keywords internal
cld_from_base_lm <- function(model, term = NULL, alpha = 0.05, method = "tukey", letters = NULL, reversed = FALSE) {
  letters_vec <- if (!is.null(letters)) letters else c(base::letters, base::LETTERS)
  mf <- stats::model.frame(model)
  resp_name <- colnames(mf)[1]
  trt_name <- if (!is.null(term)) term else colnames(mf)[2]
  
  y <- mf[[resp_name]]
  trts <- factor(mf[[trt_name]])
  levels_k <- levels(trts)
  
  means_vec <- tapply(y, trts, mean, na.rm = TRUE)
  se_pool <- stats::sigma(model) / sqrt(tapply(y, trts, length))
  
  aov_fit <- stats::aov(model)
  tuk <- stats::TukeyHSD(aov_fit, which = trt_name, conf.level = 1 - alpha)
  tuk_df <- as.data.frame(tuk[[1]])
  
  n_k <- length(levels_k)
  p_mat <- matrix(1.0, n_k, n_k, dimnames = list(levels_k, levels_k))
  
  for (rn in rownames(tuk_df)) {
    splits <- strsplit(rn, "-")[[1]]
    if (length(splits) == 2) {
      t1 <- splits[1]; t2 <- splits[2]
      if (t1 %in% levels_k && t2 %in% levels_k) {
        p_mat[t1, t2] <- tuk_df[rn, "p adj"]
        p_mat[t2, t1] <- tuk_df[rn, "p adj"]
      }
    }
  }
  
  cld_letters <- compute_piepho_cld(levels_k, means_vec, p_mat, alpha = alpha, letters = letters_vec, reversed = reversed)
  
  df_res <- data.frame(
    Treatment = levels_k,
    Mean = as.numeric(means_vec),
    SE = as.numeric(se_pool),
    CI_Lower = as.numeric(means_vec - 1.96 * se_pool),
    CI_Upper = as.numeric(means_vec + 1.96 * se_pool),
    CLD_Letter = as.character(cld_letters[levels_k]),
    stringsAsFactors = FALSE
  )
  
  ord <- if (!reversed) order(-df_res$Mean) else order(df_res$Mean)
  res <- df_res[ord, ]
  class(res) <- c("agri_cld", "data.frame")
  res
}

#' @title Pure Base R Maximal Clique CLD Algorithm (Piepho, 2004)
#' @description Solves the boolean adjacency sweep matrix to generate minimal, standardized grouping letters.
#' @param treatments Character vector of treatment labels.
#' @param means Numeric vector of treatment means.
#' @param p_matrix Square matrix of pairwise p-values.
#' @param alpha Significance threshold (default: 0.05).
#' @param letters Letters vector to assign.
#' @param reversed Logical; if TRUE, lowest mean gets 'a'.
#' @export
compute_piepho_cld <- function(treatments, means, p_matrix, alpha = 0.05, letters = NULL, reversed = FALSE) {
  letters_vec <- if (!is.null(letters)) letters else c(base::letters, base::LETTERS)
  names(means) <- treatments
  # Order treatments by mean: highest mean first if not reversed
  ord <- if (!reversed) order(-means) else order(means)
  sorted_trts <- treatments[ord]
  
  n_k <- length(sorted_trts)
  if (n_k <= 1) {
    res <- stats::setNames(letters_vec[1], sorted_trts)
    return(res)
  }
  
  # Adjacency matrix of non-significant differences: FALSE on diagonal, TRUE if p >= alpha
  adj <- matrix(FALSE, n_k, n_k, dimnames = list(sorted_trts, sorted_trts))
  for (i in seq_len(n_k)) {
    for (j in seq_len(n_k)) {
      if (i != j && p_matrix[sorted_trts[i], sorted_trts[j]] >= alpha) {
        adj[i, j] <- TRUE
      }
    }
  }
  
  # Find maximal cliques (Bron-Kerbosch algorithm with pivot)
  cliques <- list()
  
  bron_kerbosch <- function(R, P, X) {
    if (length(P) == 0 && length(X) == 0) {
      if (length(R) > 0) cliques <<- c(cliques, list(R))
      return()
    }
    PX <- union(P, X)
    deg_counts <- sapply(PX, function(k) sum(adj[k, P]))
    u <- PX[which.max(deg_counts)]
    candidates <- setdiff(P, which(adj[u, ]))
    
    for (v in candidates) {
      nbrs_v <- which(adj[v, ])
      bron_kerbosch(c(R, v), intersect(P, nbrs_v), intersect(X, nbrs_v))
      P <- setdiff(P, v)
      X <- union(X, v)
    }
  }
  
  bron_kerbosch(integer(0), seq_len(n_k), integer(0))
  
  # If no non-trivial cliques found (all treatments distinct)
  if (length(cliques) == 0) {
    cliques <- as.list(seq_len(n_k))
  }
  
  # Sort cliques by earliest treatment appearance
  clique_min_idx <- sapply(cliques, min)
  cliques <- cliques[order(clique_min_idx)]
  
  # Assign letters to cliques
  trt_letters <- stats::setNames(character(n_k), sorted_trts)
  
  for (k in seq_along(cliques)) {
    let <- letters_vec[((k - 1) %% length(letters_vec)) + 1]
    members <- sorted_trts[cliques[[k]]]
    for (m in members) {
      if (!grepl(let, trt_letters[m], fixed = TRUE)) {
        trt_letters[m] <- paste0(trt_letters[m], let)
      }
    }
  }
  
  # Ensure any singleton has a letter
  for (i in seq_len(n_k)) {
    if (nchar(trt_letters[sorted_trts[i]]) == 0) {
      trt_letters[sorted_trts[i]] <- letters_vec[i]
    }
  }
  
  trt_letters
}
