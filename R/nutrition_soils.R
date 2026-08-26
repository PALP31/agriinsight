#' @title Crop Nutrition, Soil Fertility, and Yield Plateau Models
#' @description Critical Soil Test Values (x_c), linear/quadratic plateau models, and Economically Optimum Fertilizer Rate (EOFR).
#' @name nutrition_soils
NULL

#' @title Estimate Critical Soil Nutrient Threshold (Linear/Quadratic Plateau)
#' @param nutrient Fertilizer rate or soil test nutrient concentration.
#' @param yield Crop yield response.
#' @param method Model type ("linear_plateau", "quadratic_plateau").
#' @export
get_critical_soil_value <- function(nutrient, yield, method = c("linear_plateau", "quadratic_plateau")) {
  method <- match.arg(method)
  
  if (method == "linear_plateau") {
    obj_lp <- function(par) {
      a <- par[1]; b <- par[2]; xc <- par[3]
      pred <- ifelse(nutrient < xc, a + b * nutrient, a + b * xc)
      sum((yield - pred)^2, na.rm = TRUE)
    }
    
    init_xc <- stats::median(nutrient)
    init_a <- min(yield)
    init_b <- (max(yield) - min(yield)) / init_xc
    
    opt <- stats::optim(c(a = init_a, b = init_b, xc = init_xc), obj_lp,
                        method = "L-BFGS-B", lower = c(-Inf, 0, min(nutrient)), upper = c(Inf, Inf, max(nutrient)))
    
    return(data.frame(
      Model = "Linear_Plateau",
      Critical_Threshold_xc = opt$par["xc"],
      Plateau_Yield = opt$par["a"] + opt$par["b"] * opt$par["xc"],
      Marginal_Efficiency_b = opt$par["b"],
      stringsAsFactors = FALSE
    ))
  }
}

#' @title Economically Optimum Fertilizer Rate (EOFR)
#' @param model A fitted quadratic or Mitscherlich yield model.
#' @param price_nutrient Price per unit nutrient ($/kg).
#' @param price_crop Sale price per unit crop yield ($/kg).
#' @export
get_fertilizer_optimum <- function(beta_linear, beta_quadratic, price_nutrient, price_crop) {
  # dY/dX = beta_linear + 2 * beta_quadratic * X = price_nutrient / price_crop
  price_ratio <- price_nutrient / price_crop
  eofr <- (price_ratio - beta_linear) / (2 * beta_quadratic)
  
  data.frame(
    Price_Ratio = price_ratio,
    EOFR_Optimum_Dose = pmax(0, eofr),
    stringsAsFactors = FALSE
  )
}
