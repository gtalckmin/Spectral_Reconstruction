library(pracma)

# --- Helper Functions ---

calc_lambda_i <- function(r670, r700, r740, r780) {
  rho_i <- (r670 + r780) / 2
  denom <- r740 - r700
  if (abs(denom) < 1e-6) denom <- ifelse(denom < 0, -1, 1) * 1e-6
  700 + 40 * ((rho_i - r700) / denom)
}

predict_logistic <- function(lambda, C, r_min, r_max, lambda_i) {
  r_min + (r_max - r_min) / (1 + exp(C * (lambda_i - lambda)))
}

# --- Main Reconstruction Function ---

reconstruct_spectrum <- function(obs_r, wl, params = NULL) {
  
  required_bands <- c(400, 550, 670, 700, 740, 750, 780, 980, 1100)
  anchors <- numeric(length(required_bands))
  names(anchors) <- paste0("R", required_bands)
  
  # Ensure obs_r is a simple named vector
  if(is.list(obs_r)) obs_r <- unlist(obs_r)
  
  for(b in required_bands) {
    name <- paste0("R", b)
    if(name %in% names(obs_r)) {
      anchors[name] <- as.numeric(obs_r[name])
    } else {
      warning(paste("Band", b, "missing. Using nearest neighbor."))
    }
  }
  
  # --- Segment 1: Visible (400 - 680 nm) ---
  slope_base <- (anchors["R750"] - anchors["R400"]) / (750 - 400)
  intercept_base <- anchors["R400"] - slope_base * 400
  calc_baseline <- function(l) slope_base * l + intercept_base
  
  mu <- 670 
  base_670 <- calc_baseline(670)
  A <- base_670 - anchors["R670"]
  
  base_550 <- calc_baseline(550)
  delta_550 <- base_550 - anchors["R550"]
  ratio <- delta_550 / A
  
  if (!is.null(params$sigma)) {
    sigma <- as.numeric(params$sigma)
  } else {
    if (is.na(ratio) || ratio <= 0 || ratio >= 1) {
      sigma <- 40 
    } else {
      sigma <- sqrt( - (550 - 670)^2 / (2 * log(ratio)) )
    }
  }
  
  mu_blue <- 450
  A_blue <- if(!is.null(params$A_blue)) as.numeric(params$A_blue) else 0
  sigma_blue <- if(!is.null(params$sigma_blue)) as.numeric(params$sigma_blue) else 20
  
  vis_fun <- function(l) {
    calc_baseline(l) - A * exp( - (l - mu)^2 / (2 * sigma^2) ) - A_blue * exp( - (l - mu_blue)^2 / (2 * sigma_blue^2) )
  }
  
  # --- Segment 2: Red Edge (680 - 780 nm) ---
  r670 <- as.numeric(anchors["R670"])
  r700 <- as.numeric(anchors["R700"])
  r740 <- as.numeric(anchors["R740"])
  r780 <- as.numeric(anchors["R780"])
  
  lambda_i <- calc_lambda_i(r670, r700, r740, r780)
  
  if (!is.null(params$C)) {
    C_opt <- as.numeric(params$C)
  } else {
    rss_logistic <- function(C) {
      C <- as.numeric(C)
      preds <- predict_logistic(c(700, 740), C, r670, r780, lambda_i)
      sum((preds - c(r700, r740))^2)
    }
    opt_res <- optim(0.05, rss_logistic, method="L-BFGS-B", lower=0.001, upper=2.0)
    C_opt <- as.numeric(opt_res$par)
  }
  
  re_fun <- function(l) {
    predict_logistic(as.numeric(l), C_opt, r670, r780, lambda_i)
  }
  
  # --- Segment 3: NIR Plateau (780 - 1100 nm) ---
  slope_nir <- (anchors["R1100"] - anchors["R780"]) / (1100 - 780)
  intercept_nir <- anchors["R780"] - slope_nir * 780
  calc_nir_base <- function(l) slope_nir * l + intercept_nir
  
  mu_water <- 980
  base_980 <- calc_nir_base(980)
  A_water <- base_980 - anchors["R980"]
  
  sigma_water <- if(!is.null(params$sigma_water)) as.numeric(params$sigma_water) else 40
  
  nir_fun <- function(l) {
    calc_nir_base(l) - A_water * exp( - (l - mu_water)^2 / (2 * sigma_water^2) )
  }
  
  # --- Combine ---
  out_spec <- numeric(length(wl))
  
  mask_vis <- wl <= 680
  mask_re  <- wl > 680 & wl <= 780
  mask_nir <- wl > 780
  
  if(any(mask_vis)) out_spec[mask_vis] <- vis_fun(wl[mask_vis])
  if(any(mask_re))  out_spec[mask_re]  <- re_fun(wl[mask_re])
  if(any(mask_nir)) out_spec[mask_nir] <- nir_fun(wl[mask_nir])
  
  auc_val <- trapz(wl[mask_nir], out_spec[mask_nir])

  return(list(
    spectrum = out_spec,
    features = list(
      A = A, 
      sigma = sigma, 
      lambda_i = lambda_i, 
      C = C_opt, 
      A_water = A_water, 
      AUC = auc_val,
      A_blue = A_blue,
      sigma_blue = sigma_blue,
      sigma_water = sigma_water
    )
  ))
}
