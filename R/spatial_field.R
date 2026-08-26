#' @title Spatial Biometrics and Directional Semivariograms
#' @description Empirical Matheron semivariogram calculation across spatial field trials, 
#'   directional anisotropic filtering (0, 45, 90, 135 degrees), and theoretical variogram fitting (nugget, sill, range).
#' @name spatial_field
NULL

#' @title Calculate Empirical Directional Semivariograms from Field Trial Residuals
#' @description Computes Matheron or Cressie-Hawkins empirical semivariances across distance lags and directions.
#' @param model A fitted spatial or mixed model (\code{SpATS}, \code{lmerMod}, \code{lme}, \code{lm}), or data.frame.
#' @param x_col Column name for X coordinate (range / col).
#' @param y_col Column name for Y coordinate (row / bed).
#' @param max_lag Maximum distance lag (default: 15).
#' @param lag_width Width of each distance bin (default: 1.0).
#' @param directions Vector of directional angles in degrees (default: \code{c(0, 45, 90, 135)}).
#' @param angle_tol Angular tolerance in degrees (default: 22.5).
#' @param robust Logical; whether to use Cressie-Hawkins robust estimator (default: \code{FALSE}).
#' @return An S3 object of class \code{"agri_semivariogram"}.
#' @export
get_semivariogram <- function(model, x_col = NULL, y_col = NULL, max_lag = 15, lag_width = 1.0, 
                              directions = c(0, 45, 90, 135), angle_tol = 22.5, robust = FALSE) {
  # Extract residuals and coordinates
  if (is.data.frame(model)) {
    df <- model
    resids <- df$residuals %||% df$resids %||% df[[1]]
    x_coords <- if (!is.null(x_col)) df[[x_col]] else df[[2]]
    y_coords <- if (!is.null(y_col)) df[[y_col]] else df[[3]]
  } else {
    resids <- as.numeric(residuals(model))
    mf <- if (inherits(model, "merMod")) model@frame else stats::model.frame(model)
    
    x_name <- x_col %||% grep("col|range|x|coord_x|bench_col", colnames(mf), ignore.case = TRUE, value = TRUE)[1] %||% colnames(mf)[1]
    y_name <- y_col %||% grep("row|bed|y|coord_y|pad_distance", colnames(mf), ignore.case = TRUE, value = TRUE)[1] %||% colnames(mf)[2]
    
    x_coords <- as.numeric(mf[[x_name]])
    y_coords <- as.numeric(mf[[y_name]])
  }
  
  n_pts <- length(resids)
  if (n_pts < 5) stop("At least 5 spatial points are required for semivariogram calculation.")
  
  # Distance lags
  lags <- seq(lag_width, max_lag, by = lag_width)
  
  # Compute pairwise differences, distances, and azimuth angles
  comb <- utils::combn(n_pts, 2)
  i_idx <- comb[1, ]
  j_idx <- comb[2, ]
  
  dx <- x_coords[i_idx] - x_coords[j_idx]
  dy <- y_coords[i_idx] - y_coords[j_idx]
  dist_vec <- sqrt(dx^2 + dy^2)
  diff_vec <- resids[i_idx] - resids[j_idx]
  
  # Azimuth angle in degrees: 0 to 180
  angle_rad <- atan2(dy, dx)
  angle_deg <- (angle_rad * 180 / pi) %% 180
  
  res_list <- list()
  
  # Omnidirectional
  for (l in lags) {
    in_bin <- dist_vec >= (l - lag_width / 2) & dist_vec < (l + lag_width / 2)
    n_pairs <- sum(in_bin)
    
    if (n_pairs >= 2) {
      d_vals <- diff_vec[in_bin]
      gamma_val <- if (!robust) {
        sum(d_vals^2) / (2 * n_pairs)
      } else {
        (mean(sqrt(abs(d_vals))))^4 / (2 * (0.457 + 0.494 / n_pairs))
      }
      mean_dist <- mean(dist_vec[in_bin])
    } else {
      gamma_val <- NA_real_
      mean_dist <- l
    }
    
    res_list[[length(res_list) + 1]] <- data.frame(
      Direction = "Omnidirectional",
      Angle_Deg = NA_real_,
      Lag_Nominal = l,
      Mean_Distance = mean_dist,
      Semivariance = gamma_val,
      Pairs_Count = n_pairs,
      stringsAsFactors = FALSE
    )
  }
  
  # Directional
  for (dir_deg in directions) {
    for (l in lags) {
      # Distance condition
      in_dist <- dist_vec >= (l - lag_width / 2) & dist_vec < (l + lag_width / 2)
      # Angle condition with tolerance
      d_angle <- abs(angle_deg - dir_deg)
      d_angle <- pmin(d_angle, 180 - d_angle)
      in_angle <- d_angle <= angle_tol
      
      in_bin <- in_dist & in_angle
      n_pairs <- sum(in_bin)
      
      if (n_pairs >= 2) {
        d_vals <- diff_vec[in_bin]
        gamma_val <- if (!robust) {
          sum(d_vals^2) / (2 * n_pairs)
        } else {
          (mean(sqrt(abs(d_vals))))^4 / (2 * (0.457 + 0.494 / n_pairs))
        }
        mean_dist <- mean(dist_vec[in_bin])
      } else {
        gamma_val <- NA_real_
        mean_dist <- l
      }
      
      res_list[[length(res_list) + 1]] <- data.frame(
        Direction = paste0("Dir_", dir_deg, "deg"),
        Angle_Deg = dir_deg,
        Lag_Nominal = l,
        Mean_Distance = mean_dist,
        Semivariance = gamma_val,
        Pairs_Count = n_pairs,
        stringsAsFactors = FALSE
      )
    }
  }
  
  full_df <- do.call(rbind, res_list)
  
  # Theoretical Spherical Variogram Fit (Omnidirectional)
  omni_sub <- full_df[full_df$Direction == "Omnidirectional" & !is.na(full_df$Semivariance), ]
  
  fitted_params <- if (nrow(omni_sub) >= 3) {
    total_var <- stats::var(resids)
    c0_init <- min(omni_sub$Semivariance)
    c_init <- pmax(1e-4, total_var - c0_init)
    a_init <- max_lag * 0.4
    
    obj_sph <- function(p) {
      c0 <- p[1]; c_s <- p[2]; a <- p[3]
      h <- omni_sub$Mean_Distance
      pred <- ifelse(h <= a, c0 + c_s * (1.5 * (h / a) - 0.5 * (h / a)^3), c0 + c_s)
      sum(omni_sub$Pairs_Count * (omni_sub$Semivariance - pred)^2)
    }
    
    opt_sph <- tryCatch(
      stats::optim(c(c0 = c0_init, c = c_init, a = a_init), obj_sph,
                   method = "L-BFGS-B", lower = c(0, 1e-4, 0.5), upper = c(total_var * 2, total_var * 2, max_lag * 2)),
      error = function(e) NULL
    )
    
    if (!is.null(opt_sph)) {
      c0_est <- max(0, as.numeric(opt_sph$par["c0"]))
      c_est <- max(1e-4, as.numeric(opt_sph$par["c"]))
      a_est <- max(0.1, as.numeric(opt_sph$par["a"]))
      data.frame(
        Model = "Spherical",
        Nugget_c0 = c0_est,
        Partial_Sill_c = c_est,
        Total_Sill = c0_est + c_est,
        Range_a = a_est,
        Spatial_Dependence_Ratio = c_est / (c0_est + c_est),
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
  } else {
    NULL
  }
  
  out <- list(
    semivariogram = full_df,
    theoretical_fit = fitted_params,
    total_sample_variance = stats::var(resids),
    sample_size = n_pts
  )
  class(out) <- c("agri_semivariogram", "list")
  out
}

#' @export
print.agri_semivariogram <- function(x, ...) {
  cat("=== Empirical Spatial Semivariogram (Matheron Estimator) ===\n")
  cat(sprintf("Spatial sample size: %d | Total Variance: %.4f\n", x$sample_size, x$total_sample_variance))
  if (!is.null(x$theoretical_fit)) {
    cat("\nFitted Theoretical Model:\n")
    print(x$theoretical_fit, row.names = FALSE)
  }
  cat("\nOmnidirectional Lags:\n")
  omni <- x$semivariogram[x$semivariogram$Direction == "Omnidirectional", ]
  print(head(omni, 6), row.names = FALSE)
  invisible(x)
}
