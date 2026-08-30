#' @title Crop Nutrition, Soil Fertility, and Yield Plateau Models
#' @description Critical Soil Test Values (x_c), linear/quadratic plateau models, Mitscherlich-Bray non-linear saturation,
#'   and Economically Optimum Fertilizer Rate (EOFR) with price-ratio sensitivity.
#' @param nutrient Fertilizer application rate or soil test nutrient concentration (vector).
#' @param yield Crop yield response (vector).
#' @param method Model type: \code{"linear_plateau"}, \code{"quadratic_plateau"}, \code{"mitscherlich_bray"}, \code{"quadratic"}.
#' @param conf_level Confidence level for parameters (default: \code{0.95}).
#' @name nutrition_soils
NULL

#' @title Estimate Critical Soil Nutrient Threshold (Linear/Quadratic Plateau)
#' @description Fits non-linear yield plateau models to extract the critical soil test value (\eqn{x_c}) or optimum agronomic rate.
#' @param nutrient Fertilizer rate or soil test nutrient concentration.
#' @param yield Crop yield response.
#' @param method Model type: \code{"linear_plateau"}, \code{"quadratic_plateau"}, \code{"mitscherlich_bray"}, \code{"quadratic"}.
#' @param conf_level Confidence level for intervals (default: \code{0.95}).
#' @return An S3 object of class \code{"agri_plateau_model"}.
#' @export
get_critical_soil_value <- function(nutrient, yield, method = c("linear_plateau", "quadratic_plateau", "mitscherlich_bray", "quadratic"), conf_level = 0.95) {
  method <- match.arg(method)
  
  clean_df <- stats::na.omit(data.frame(x = as.numeric(nutrient), y = as.numeric(yield)))
  x <- clean_df$x
  y <- clean_df$y
  n <- length(x)
  
  if (n < 4) stop("At least 4 data points are required to fit plateau response models.")
  
  ss_tot <- sum((y - mean(y))^2)
  
  if (method == "linear_plateau") {
    # Objective function
    obj_lp <- function(par) {
      a <- par[1]; b <- par[2]; xc <- par[3]
      pred <- ifelse(x < xc, a + b * x, a + b * xc)
      sum((y - pred)^2)
    }
    
    init_xc <- stats::median(x)
    init_a <- min(y)
    init_b <- pmax(1e-4, (max(y) - min(y)) / pmax(1e-4, init_xc - min(x)))
    
    opt <- stats::optim(c(a = init_a, b = init_b, xc = init_xc), obj_lp,
                        method = "L-BFGS-B",
                        lower = c(-Inf, 0, min(x) + 1e-4),
                        upper = c(Inf, Inf, max(x) - 1e-4))
    
    a_est <- opt$par["a"]
    b_est <- opt$par["b"]
    xc_est <- opt$par["xc"]
    y_max <- a_est + b_est * xc_est
    ss_res <- opt$value
    r2 <- pmax(0, 1 - (ss_res / ss_tot))
    rmse <- sqrt(ss_res / (n - 3))
    
    # Delta-method approximate SE for xc
    se_xc <- rmse / (b_est * sqrt(n))
    
    res_df <- data.frame(
      Model = "Linear_Plateau",
      Critical_Threshold_xc = as.numeric(xc_est),
      SE_xc = as.numeric(se_xc),
      Plateau_Yield = as.numeric(y_max),
      Intercept_a = as.numeric(a_est),
      Marginal_Efficiency_b = as.numeric(b_est),
      R_Squared = as.numeric(r2),
      RMSE = as.numeric(rmse),
      AIC = n * log(ss_res / n) + 2 * 3,
      stringsAsFactors = FALSE
    )
    
    out <- list(
      model_type = "Linear_Plateau",
      parameters = res_df,
      coefficients = opt$par,
      fitted_values = ifelse(x < xc_est, a_est + b_est * x, y_max),
      residuals = y - ifelse(x < xc_est, a_est + b_est * x, y_max),
      data = clean_df
    )
    class(out) <- c("agri_plateau_model", "list")
    return(out)
  }
  
  if (method == "quadratic_plateau") {
    # Smooth Quadratic Plateau: y = a + b*x + c*x^2 (x < xc), y = a - b^2/(4c) (x >= xc) where xc = -b/(2c)
    obj_qp <- function(par) {
      a <- par[1]; b <- par[2]; c_val <- par[3]
      if (c_val >= 0) return(1e12) # must be concave downward
      xc <- -b / (2 * c_val)
      y_plat <- a - (b^2) / (4 * c_val)
      pred <- ifelse(x < xc, a + b * x + c_val * (x^2), y_plat)
      sum((y - pred)^2)
    }
    
    # Starting values from regular quadratic fit
    lm_init <- stats::lm(y ~ x + I(x^2))
    cf_init <- coef(lm_init)
    init_a <- cf_init[1]
    init_b <- pmax(1e-3, cf_init[2])
    init_c <- pmin(-1e-4, cf_init[3])
    
    opt <- stats::optim(c(a = init_a, b = init_b, c = init_c), obj_qp,
                        method = "L-BFGS-B",
                        lower = c(-Inf, 1e-6, -Inf),
                        upper = c(Inf, Inf, -1e-8),
                        hessian = TRUE)
    
    a_est <- opt$par["a"]
    b_est <- opt$par["b"]
    c_est <- opt$par["c"]
    xc_est <- -b_est / (2 * c_est)
    y_max <- a_est - (b_est^2) / (4 * c_est)
    ss_res <- opt$value
    r2 <- pmax(0, 1 - (ss_res / ss_tot))
    rmse <- sqrt(ss_res / (n - 3))
    
    # Delta-method exact SE for xc = -b / (2c)
    se_xc <- NA_real_
    if (!is.null(opt$hessian)) {
      cov_mat <- tryCatch(solve(opt$hessian) * (ss_res / (n - 3)), error = function(e) NULL)
      if (!is.null(cov_mat) && all(diag(cov_mat) > 0)) {
        grad_g <- c(-1 / (2 * c_est), b_est / (2 * c_est^2))
        cov_bc <- cov_mat[2:3, 2:3]
        var_xc <- as.numeric(t(grad_g) %*% cov_bc %*% grad_g)
        if (var_xc > 0) se_xc <- sqrt(var_xc)
      }
    }
    if (is.na(se_xc)) {
      se_xc <- sqrt(abs(xc_est) / (2 * n))
    }
    
    res_df <- data.frame(
      Model = "Quadratic_Plateau",
      Critical_Threshold_xc = as.numeric(xc_est),
      SE_xc = as.numeric(se_xc),
      Plateau_Yield = as.numeric(y_max),
      Intercept_a = as.numeric(a_est),
      Linear_Coef_b = as.numeric(b_est),
      Quadratic_Coef_c = as.numeric(c_est),
      R_Squared = as.numeric(r2),
      RMSE = as.numeric(rmse),
      AIC = n * log(ss_res / n) + 2 * 3,
      stringsAsFactors = FALSE
    )
    
    out <- list(
      model_type = "Quadratic_Plateau",
      parameters = res_df,
      coefficients = c(opt$par, xc = xc_est, y_max = y_max),
      fitted_values = ifelse(x < xc_est, a_est + b_est * x + c_est * (x^2), y_max),
      residuals = y - ifelse(x < xc_est, a_est + b_est * x + c_est * (x^2), y_max),
      data = clean_df
    )
    class(out) <- c("agri_plateau_model", "list")
    return(out)
  }
  
  if (method == "mitscherlich_bray") {
    # Mitscherlich saturation: y = A * (1 - exp(-c * (x + d)))
    obj_mb <- function(par) {
      A <- par[1]; c_val <- par[2]; d <- par[3]
      pred <- A * (1 - exp(-c_val * (x + d)))
      sum((y - pred)^2)
    }
    
    init_A <- max(y) * 1.05
    init_c <- 0.05
    init_d <- pmax(0, min(x))
    
    opt <- stats::optim(c(A = init_A, c = init_c, d = init_d), obj_mb,
                        method = "L-BFGS-B",
                        lower = c(max(y) * 0.8, 1e-4, 0),
                        upper = c(max(y) * 2.0, 1.0, max(x)))
    
    A_est <- opt$par["A"]
    c_est <- opt$par["c"]
    d_est <- opt$par["d"]
    # 95% sufficiency threshold: x_95 = ln(1 / 0.05) / c - d
    xc_95 <- log(20) / c_est - d_est
    ss_res <- opt$value
    r2 <- pmax(0, 1 - (ss_res / ss_tot))
    rmse <- sqrt(ss_res / (n - 3))
    
    res_df <- data.frame(
      Model = "Mitscherlich_Bray",
      Critical_Threshold_xc = as.numeric(pmax(0, xc_95)),
      SE_xc = NA_real_,
      Plateau_Yield = as.numeric(A_est),
      Asymptote_A = as.numeric(A_est),
      Curvature_c = as.numeric(c_est),
      Indigenous_Nutrient_d = as.numeric(d_est),
      R_Squared = as.numeric(r2),
      RMSE = as.numeric(rmse),
      AIC = n * log(ss_res / n) + 2 * 3,
      stringsAsFactors = FALSE
    )
    
    out <- list(
      model_type = "Mitscherlich_Bray",
      parameters = res_df,
      coefficients = opt$par,
      fitted_values = A_est * (1 - exp(-c_est * (x + d_est))),
      residuals = y - (A_est * (1 - exp(-c_est * (x + d_est)))),
      data = clean_df
    )
    class(out) <- c("agri_plateau_model", "list")
    return(out)
  }
  
  if (method == "quadratic") {
    lm_fit <- stats::lm(y ~ x + I(x^2))
    cf <- coef(lm_fit)
    a_est <- cf[1]; b_est <- cf[2]; c_est <- cf[3]
    xc_est <- if (c_est < 0) -b_est / (2 * c_est) else max(x)
    y_max <- a_est + b_est * xc_est + c_est * (xc_est^2)
    ss_res <- sum(residuals(lm_fit)^2)
    
    res_df <- data.frame(
      Model = "Quadratic_Polynomial",
      Critical_Threshold_xc = as.numeric(xc_est),
      SE_xc = NA_real_,
      Plateau_Yield = as.numeric(y_max),
      Intercept_a = as.numeric(a_est),
      Linear_Coef_b = as.numeric(b_est),
      Quadratic_Coef_c = as.numeric(c_est),
      R_Squared = summary(lm_fit)$r.squared,
      RMSE = summary(lm_fit)$sigma,
      AIC = stats::AIC(lm_fit),
      stringsAsFactors = FALSE
    )
    
    out <- list(
      model_type = "Quadratic_Polynomial",
      parameters = res_df,
      coefficients = cf,
      fitted_values = fitted(lm_fit),
      residuals = residuals(lm_fit),
      data = clean_df
    )
    class(out) <- c("agri_plateau_model", "list")
    return(out)
  }
}

#' @title Economically Optimum Fertilizer Rate (EOFR)
#' @description Computes the Economically Optimum Fertilizer Rate (EOFR) and Agronomic Optimum Rate (AOFR) 
#'   based on nutrient:crop price ratios.
#' @param model A fitted \code{agri_plateau_model} object, or numeric linear coefficient \code{beta_linear}.
#' @param price_nutrient Price per unit nutrient ($/kg).
#' @param price_crop Sale price per unit crop yield ($/kg).
#' @param beta_quadratic Quadratic coefficient if model is not supplied directly.
#' @export
get_fertilizer_optimum <- function(model, price_nutrient, price_crop, beta_quadratic = NULL) {
  price_ratio <- price_nutrient / price_crop
  
  if (inherits(model, "agri_plateau_model")) {
    params <- model$parameters
    mtype <- model$model_type
    
    if (mtype == "Linear_Plateau") {
      b_val <- params$Marginal_Efficiency_b
      xc_val <- params$Critical_Threshold_xc
      eofr <- if (b_val >= price_ratio) xc_val else 0
      aofr <- xc_val
      yield_opt <- if (eofr >= xc_val) params$Plateau_Yield else params$Intercept_a + b_val * eofr
    } else if (mtype == "Quadratic_Plateau") {
      b_val <- params$Linear_Coef_b
      c_val <- params$Quadratic_Coef_c
      xc_val <- params$Critical_Threshold_xc
      raw_eofr <- (price_ratio - b_val) / (2 * c_val)
      eofr <- pmin(xc_val, pmax(0, raw_eofr))
      aofr <- xc_val
      yield_opt <- if (eofr >= xc_val) params$Plateau_Yield else params$Intercept_a + b_val * eofr + c_val * (eofr^2)
    } else if (mtype == "Mitscherlich_Bray") {
      A_val <- params$Asymptote_A
      c_val <- params$Curvature_c
      d_val <- params$Indigenous_Nutrient_d
      raw_eofr <- log(pmax(1e-6, (A_val * c_val) / price_ratio)) / c_val - d_val
      eofr <- pmax(0, raw_eofr)
      aofr <- params$Critical_Threshold_xc
      yield_opt <- A_val * (1 - exp(-c_val * (eofr + d_val)))
    } else {
      b_val <- params$Linear_Coef_b
      c_val <- params$Quadratic_Coef_c
      raw_eofr <- (price_ratio - b_val) / (2 * c_val)
      eofr <- pmax(0, raw_eofr)
      aofr <- params$Critical_Threshold_xc
      yield_opt <- params$Intercept_a + b_val * eofr + c_val * (eofr^2)
    }
    
    net_return <- yield_opt * price_crop - eofr * price_nutrient
    
    return(data.frame(
      Model = mtype,
      Price_Ratio = as.numeric(price_ratio),
      EOFR_Economic_Optimum = as.numeric(eofr),
      AOFR_Agronomic_Optimum = as.numeric(aofr),
      Yield_at_EOFR = as.numeric(yield_opt),
      Net_Return_per_ha = as.numeric(net_return),
      stringsAsFactors = FALSE
    ))
  }
  
  # Fallback for direct beta coefficients
  b_linear <- as.numeric(model)
  b_quad <- as.numeric(beta_quadratic %||% -0.001)
  raw_eofr <- (price_ratio - b_linear) / (2 * b_quad)
  
  data.frame(
    Price_Ratio = price_ratio,
    EOFR_Optimum_Dose = pmax(0, raw_eofr),
    stringsAsFactors = FALSE
  )
}

#' @export
print.agri_plateau_model <- function(x, ...) {
  cat(sprintf("=== %s Soil Fertility & Yield Response Model ===\n", x$model_type))
  print(x$parameters, row.names = FALSE)
  invisible(x)
}
