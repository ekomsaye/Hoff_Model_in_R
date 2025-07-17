
# Load required libraries with warning suppression
suppressWarnings({
  library(ncdf4)
  library(dplyr) # For %>% operator
  library(lubridate) # For year and month functions
  library(zoo)
  library(gridExtra)
  library(tidyverse)
})

# Set working directory
setwd("C:/Users/Lenovo/Desktop/Upwork/Bisan")
cat("Working directory set to:", getwd(), "\n")

# Define the Hoff functions (covreg.em)
`covreg.em` <- function(formula, data = NULL, R = 5, tol = 1e-8, itmax = 50) {
  model <- model.frame(formula, data)
  Y <- model.response(model)
  X <- model.matrix(formula, model)
  
  # Normal log likelihood function
  ldmvnorm <- function(y, mu = rep(0, length(y)), Sig = diag(1, length(y))) {
    -.5 * (length(y) * log(2 * pi) + log(det(Sig)) + t(y - mu) %*% solve(Sig) %*% (y - mu))
  }
  
  p <- dim(Y)[2]; q <- dim(X)[2]; n <- dim(Y)[1]
  A <- cov(Y); iA <- solve(A); B <- matrix(rnorm(p * q * R), p, q * R) * 1e-4; iter <- 0
  LL <- NULL
  rll <- 10
  while (rll > tol & iter < itmax) {
    B0 <- B; iter <- iter + 1
    
    ### Find expectation, var of z
    Vz <- array(dim = c(R, R, n)); Mz <- matrix(nrow = n, ncol = R)
    for (i in 1:n) {
      Bx <- apply(array(B, dim = c(p, q, R)), 3, "%*%", X[i,])
      Vz[,,i] <- solve(t(Bx) %*% iA %*% Bx + diag(R))
      Mz[i,] <- Vz[,,i] %*% t(Bx) %*% iA %*% Y[i,]
    }
    
    ### Obtain MLEs
    Y1 <- Y; X1 <- NULL; for (r in 1:R) { X1 <- cbind(X1, diag(Mz[,r]) %*% X) }
    Y0 <- matrix(0, nrow = n * R, ncol = p); X0 <- NULL
    for (i in 1:n) {
      xi <- matrix(outer(X[i,], diag(R)), nrow = R * q, ncol = R)
      ZZ <- xi %*% Vz[,,i] %*% t(xi); ZZ <- .5 * (ZZ + t(ZZ))
      Z <- eigen(ZZ); Z <- Z$vec[,1:R] %*% diag(sqrt(Z$val[1:R]), nrow = R)
      X0 <- rbind(X0, t(Z))
    }
    YA <- rbind(Y0, Y1); XA <- rbind(X0, X1)
    
    B <- t(YA) %*% XA %*% solve(t(XA) %*% XA)
    E <- YA - XA %*% t(B)
    A <- (t(E) %*% E) / dim(Y)[1]; iA <- solve(A)
    
    ###
    if (iter %% 5 == 0) {
      ll <- 0
      for (i in 1:dim(Y)[1]) {
        xi <- matrix(outer(X[i,], diag(R)), nrow = R * q, ncol = R)
        ll <- ll + ldmvnorm(Y[i,], Sig = A + B %*% xi %*% t(xi) %*% t(B))
      }
      LL <- c(LL, ll)
      if (iter > 5) { rll <- abs(LL[length(LL)] - LL[length(LL)-1]) / abs(LL[length(LL)]) }
      cat(iter, log(rll, base = 10), ll, "\n")
    }
  }
  list(A = A, B = B, ll = ll, LL = LL)
}

# Step 1: Load NetCDF Data
tryCatch({
  nc_data <- nc_open("C:/Users/Lenovo/Desktop/Upwork/Bisan/HadISST_sst(1).nc")
  longitude <- ncvar_get(nc_data, "longitude")
  latitude <- ncvar_get(nc_data, "latitude")
  time <- ncvar_get(nc_data, "time")
  S <- ncvar_get(nc_data, "sst")
  S[S == ncatt_get(nc_data, "sst", "_FillValue")$value] <- NA
  nc_close(nc_data)
}, error = function(e) {
  stop("Error loading NetCDF file: ", e$message)
})

# Verify full dataset dimensions
cat("Dimensions of S (full dataset):", dim(S), "\n")

# Preprocess Data
longitude_shifted <- ifelse(longitude < 0, longitude + 360, longitude)
idx <- order(longitude_shifted)
longitude <- longitude_shifted[idx]
S <- S[idx, , ]
n_long <- length(longitude)
cat("Number of longitudes after reordering:", n_long, "\n")

start_date <- as.Date("1870-01-01")
dates <- start_date + time
years <- year(dates)
months <- month(dates)
time_vec <- years + (months - 1) / 12
n <- length(time_vec)
cat("Time range of time_vec:", range(time_vec), "\n")

lat_eq <- which.min(abs(latitude))
S_eq <- t(S[, lat_eq, ])
cat("Dimensions of S_eq:", dim(S_eq), "\n")

# Remove Monthly Means
monthly_clim <- tapply(1:nrow(S_eq), months, function(idx) {
  colMeans(S_eq[idx, , drop = FALSE], na.rm = TRUE)
}) %>% do.call(rbind, .)

if (nrow(monthly_clim) != 12) {
  stop("Error: monthly_clim should have 12 rows, got ", nrow(monthly_clim))
}

monthly_clim[is.na(monthly_clim)] <- mean(monthly_clim, na.rm = TRUE)
X1 <- S_eq - monthly_clim[months, ]

cat("Dimensions of X1 (anomalies after monthly means removed):", dim(X1), "\n")
cat("Range of X1 anomalies:", range(X1, na.rm = TRUE), "\n")

longitude_ocean <- longitude[!is.na(X1[1, ])]
X1a <- X1[, !is.na(X1[1, ])]
p <- length(longitude_ocean)
cat("Number of ocean longitudes (p):", p, "\n")

ym <- t(X1a)
cat("Dimensions of ym:", dim(ym), "\n")











# ... (previous code up to kernel filtering)

# Kernel Filter (reduced span)
span <- 5 * 12  # 5-year span
w <- seq(-1, 1, length = span)
w <- (1 - w^2) * 3/4
w <- w / sum(w)
ymd <- ym
for (k in 1:p) {
  y1 <- ym[k, ]
  yf <- stats::filter(y1, filter = w, method = "convolution", sides = 2)
  yf[is.infinite(yf) | is.na(yf)] <- NA
  ymd[k,] <- ifelse(is.na(yf), NA, y1 - yf)
  if (sd(ymd[k,], na.rm = TRUE) < 1e-5) cat("Warning: Longitude", k, "has low variability after filtering (sd =", sd(ymd[k,], na.rm = TRUE), ").\n")
}

a <- span / 2
time_veca <- time_vec[a:(n - a)]
ymd_trimmed <- ymd[, a:(n - a)]

if (any(is.infinite(ymd_trimmed) | is.na(ymd_trimmed))) {
  cat("Warning: ymd_trimmed contains", sum(is.infinite(ymd_trimmed)), "infinite and", sum(is.na(ymd_trimmed)), "NA values. Replacing with mean.\n")
  mean_val <- mean(ymd_trimmed, na.rm = TRUE)
  ymd_trimmed[is.infinite(ymd_trimmed) | is.na(ymd_trimmed)] <- mean_val
}

# PCA Calculation
if (ncol(ymd[, a:(n - a)]) != length(time_veca)) stop("Dimension mismatch between ymd and time_veca.")
out <- tryCatch({
  svd_result <- svd(t(ymd[, a:(n - a)]))
  if (any(svd_result$d < 1e-10)) {
    cat("Warning: Some singular values are near zero:", svd_result$d[svd_result$d < 1e-10], ". Adjusting to minimum threshold.\n")
    svd_result$d[svd_result$d < 1e-10] <- 1e-10
  }
  cat("Singular values (out$d):", round(svd_result$d[1:10], 6), "\n")
  var_explained <- 100 * (svd_result$d^2 / sum(svd_result$d^2))
  cat("Variance explained (%):", round(var_explained[1:10], 2), "\n")
  svd_result
}, error = function(e) {
  stop("SVD failed: ", conditionMessage(e), ". Check data variability or ymd dimensions.")
})

plot(100 * (out$d^2 / sum(out$d^2))[1:50], type = "b")

K_values <- c(1, 4, 10)
pca_results <- list()
for (K in K_values) {
  cat("Processing PCA with K =", K, "principal components\n")
  K <- min(K, length(out$d))
  d <- c(out$d[1:K], rep(0, length(out$d) - K))
  X1r <- out$u %*% diag(d) %*% t(out$v)
  PC <- as.matrix(out$u[, 1:K])
  pca_results[[paste0("K_", K)]] <- list(out = out, PC = PC, var_explained = 100 * (out$d^2 / sum(out$d^2))[1:K])
}
cat("PCA calculation completed for all K values.\n")

fits <- list()
for (fit_key in names(pca_results)) {
  K_test <- as.numeric(gsub("K_", "", fit_key))
  out <- pca_results[[fit_key]]$out
  PC <- as.matrix(pca_results[[fit_key]]$PC)
  var_explained <- pca_results[[fit_key]]$var_explained
  K <- min(K_test, length(out$d))
  
  cat("Fitting Hoff model for K =", K, "\n")
  
  # Normalize PCs
  PC_norm <- scale(PC[, 1:K, drop = FALSE], center = TRUE, scale = FALSE)
  pc_sds <- apply(PC[, 1:K, drop = FALSE], 2, sd, na.rm = TRUE)
  
  # Scale time
  time_scaled <- (time_veca - mean(time_veca)) / (max(time_veca) - min(time_veca)) * 2
  
  fit <- tryCatch({
    result <- covreg.em(PC_norm ~ time_scaled, R = 10, tol = 1e-8, itmax = 50)
    
    # Extract Hoff components
    A <- result$A
    B <- result$B
    C <- A + B[, 1] %*% t(B[, 1])
    D <- B[, 1] %*% t(B[, 2]) + B[, 2] %*% t(B[, 1])
    E <- B[, 2] %*% t(B[, 2])
    
    # Project to grid
    V <- out$v[, 1:K]
    CG <- V %*% C %*% t(V)
    DG <- V %*% D %*% t(V)
    EG <- V %*% E %*% t(V)
    
    # Extract diagonals
    CGV <- diag(CG)
    DGV <- diag(DG)
    EGV <- diag(EG)
    
    # Apply dynamic floors for stability
    EGV <- pmax(EGV, 0.00001)
    DGV <- ifelse(abs(DGV) < 1e-6, -0.0001, DGV)
    
    list(out = out, fit = result, C = C, D = D, E = E, CG = CG, DG = DG, EG = EG, CGV = CGV, DGV = DGV, EGV = EGV, var_explained = var_explained, pc_sds = pc_sds)
  }, error = function(e) {
    cat("⚠️ Hoff model failed for K =", K, ":", e$message, "\n")
    NULL
  })
  
  fits[[fit_key]] <- fit
}

cat("✅ Hoff model fits with normalization complete.\n")

# ... (rest of the code, including the plotting section as provided in the last response)
cat("✅ Hoff model fits with normalization complete.\n")

# fit_key <- "K_4"
# fit_obj <- fits[[fit_key]]
# plot_list <- list()
# 
# if (!is.null(fit_obj)) {
#   out <- fit_obj$out
#   C <- fit_obj$C
#   D <- fit_obj$D
#   E <- fit_obj$E
#   PC <- out$u[, 1:10]
#   time_actual <- time_veca
#   
#   for (k in 1:3) {
#     cat("→ Plotting PC", k, "with Envelopes...\n")
#     pc_series <- PC[, k]
#     pc_sd <- sd(pc_series, na.rm = TRUE)
#     pc_mean <- mean(pc_series, na.rm = TRUE)
#     cat("PC", k, "series variance:", var(pc_series, na.rm = TRUE), "\n")
#     
#     # ±1σ Band
#     upper_sigma <- rep(pc_mean + pc_sd, length(pc_series))
#     lower_sigma <- rep(pc_mean - pc_sd, length(pc_series))
#     df_sigma <- data.frame(time = time_actual, upper = upper_sigma, lower = lower_sigma)
#     
#     # Kernel Envelope
#     pc_loess <- loess(pc_series ~ time_actual, span = 0.3)
#     pc_trend <- predict(pc_loess, newdata = time_actual)
#     pc_detrended <- pc_series - pc_trend
#     roll_var <- zoo::rollapply(pc_detrended, width = 240, FUN = function(x) var(x, na.rm = TRUE), align = "center", fill = NA)
#     valid_idx <- which(!is.na(roll_var))
#     time_valid <- time_actual[valid_idx]
#     roll_var_valid <- roll_var[valid_idx]
#     roll_var_smooth <- loess(roll_var_valid ~ time_valid, span = 0.75)$fitted
#     roll_var_interp <- zoo::na.approx(zoo(roll_var_smooth, order.by = time_valid), xout = time_actual, na.rm = FALSE)
#     roll_sd_smooth <- sqrt(pmax(roll_var_interp, 1e-6))
#     upper_kernel <- pc_trend + roll_sd_smooth
#     lower_kernel <- pc_trend - roll_sd_smooth
#     df_kernel <- data.frame(time = time_actual, upper = upper_kernel, lower = lower_kernel)
#     
#     # Hoff Quadratic Envelope
#     a_k <- C[k, k]
#     b_k <- 2 * D[k, k]
#     c_k <- E[k, k]
#     t_fine <- seq(min(time_actual), max(time_actual), length.out = 1000)
#     t_scaled <- (t_fine - mean(time_veca)) / (max(time_veca) - min(time_veca)) * 2
#     var_fine <- a_k + b_k * t_scaled + 10*c_k * t_scaled^2
#     var_scaled <- pmax(var_fine, 1e-6)
#     pc_loess <- loess(pc_series ~ time_actual, span = 0.25)
#     pc_fine <- predict(pc_loess, newdata = t_fine)
#     upper_hoff <- pc_fine + sqrt(var_scaled)
#     lower_hoff <- pc_fine - sqrt(var_scaled)
#     df_hoff <- data.frame(time = t_fine, PC = pc_fine, upper = upper_hoff, lower = lower_hoff)
#     
#     # Verify 1σ coverage
#     hoff_upper_interp <- approx(df_hoff$time, df_hoff$upper, xout = time_actual)$y
#     hoff_lower_interp <- approx(df_hoff$time, df_hoff$lower, xout = time_actual)$y
#     within_1sigma <- (pc_series >= hoff_lower_interp) & (pc_series <= hoff_upper_interp)
#     coverage <- mean(within_1sigma, na.rm = TRUE)
#     cat("PC", k, "Hoff 1σ coverage:", coverage * 100, "%\n")
#     
#     # Plot
#     p <- ggplot() +
#       geom_ribbon(data = df_sigma, aes(x = time, ymin = lower, ymax = upper, fill = "±1σ Band"), alpha = 0.2) +
#       geom_ribbon(data = df_kernel, aes(x = time, ymin = lower, ymax = upper, fill = "Kernel Envelope"), alpha = 0.25) +
#       geom_ribbon(data = df_hoff, aes(x = time, ymin = lower, ymax = upper, fill = "Hoff Envelope"), alpha = 0.3) +
#       geom_line(data = data.frame(time = time_actual, PC = pc_series), aes(x = time, y = PC, color = "PC Series"), size = 0.8) +
#       scale_fill_manual(name = "Envelope Type", values = c("±1σ Band" = "grey80", "Kernel Envelope" = "green", "Hoff Envelope" = "blue")) +
#       scale_color_manual(name = "PC Line", values = c("PC Series" = "black")) +
#       labs(title = paste("PC", k, "with Variance Envelopes"), x = "Years", y = paste("PC", k)) +
#       theme_minimal(base_size = 14) +
#       theme(legend.position = "top")
#     print(p)
#     plot_list[[k]] <- p
#     Sys.sleep(15)
#   }
#   
#   grid_plot <- grid.arrange(grobs = plot_list, ncol = 1)
#   print(grid_plot)
# }
# 
# fit_key <- "K_4"
# fit_obj <- fits[[fit_key]]
# plot_list <- list()
# 
# if (!is.null(fit_obj)) {
#   out <- fit_obj$out
#   C <- fit_obj$C
#   D <- fit_obj$D
#   E <- fit_obj$E
#   PC <- out$u[, 1:10]
#   time_actual <- time_veca
# 
#   for (k in 1:3) {
#     cat("→ Plotting PC", k, "with Envelopes...\n")
#     pc_series <- PC[, k]
#     pc_sd <- sd(pc_series, na.rm = TRUE)
#     pc_mean <- mean(pc_series, na.rm = TRUE)
#     cat("PC", k, "series variance:", var(pc_series, na.rm = TRUE), "\n")
# 
#     # ±1σ Band
#     upper_sigma <- rep(pc_mean + pc_sd, length(pc_series))
#     lower_sigma <- rep(pc_mean - pc_sd, length(pc_series))
#     df_sigma <- data.frame(time = time_actual, upper = upper_sigma, lower = lower_sigma)
# 
#     # Kernel Envelope
#     pc_loess <- loess(pc_series ~ time_actual, span = 0.3)
#     pc_trend <- predict(pc_loess, newdata = time_actual)
#     pc_detrended <- pc_series - pc_trend
#     roll_var <- zoo::rollapply(pc_detrended, width = 240, FUN = function(x) var(x, na.rm = TRUE), align = "center", fill = NA)
#     valid_idx <- which(!is.na(roll_var))
#     time_valid <- time_actual[valid_idx]
#     roll_var_valid <- roll_var[valid_idx]
#     roll_var_smooth <- loess(roll_var_valid ~ time_valid, span = 0.75)$fitted
#     roll_var_interp <- zoo::na.approx(zoo(roll_var_smooth, order.by = time_valid), xout = time_actual, na.rm = FALSE)
#     roll_sd_smooth <- sqrt(pmax(roll_var_interp, 1e-6))
#     upper_kernel <- pc_trend + roll_sd_smooth
#     lower_kernel <- pc_trend - roll_sd_smooth
#     df_kernel <- data.frame(time = time_actual, upper = upper_kernel, lower = lower_kernel)
# 
#     # Hoff Quadratic Envelope
#     a_k <- C[k, k]
#     b_k <- 2 * D[k, k]
#     c_k <- E[k, k]
#     t_fine <- seq(min(time_actual), max(time_actual), length.out = 1000)
#     t_scaled <- (t_fine - mean(time_veca)) / (max(time_veca) - min(time_veca)) * 2
# 
#     # Scale all terms proportionally (e.g., factor of 10)
#     scale_factor <- 100
#     var_fine <- (a_k) + (b_k *2.5) * t_scaled + (c_k * scale_factor) * t_scaled^2
# 
#     # Cap variance to prevent over-enclosure
#     var_cap <- var(pc_series, na.rm = TRUE) * 2  # 2x PC variance as upper limit
#     var_scaled <- pmax(pmin(var_fine, var_cap), 1e-6)
# 
#     pc_loess <- loess(pc_series ~ time_actual, span = 0.25)
#     pc_fine <- predict(pc_loess, newdata = t_fine)
#     upper_hoff <- pc_fine + sqrt(var_scaled)
#     lower_hoff <- pc_fine - sqrt(var_scaled)
#     df_hoff <- data.frame(time = t_fine, PC = pc_fine, upper = upper_hoff, lower = lower_hoff)
# 
#     # Verify 1σ coverage
#     hoff_upper_interp <- approx(df_hoff$time, df_hoff$upper, xout = time_actual)$y
#     hoff_lower_interp <- approx(df_hoff$time, df_hoff$lower, xout = time_actual)$y
#     within_1sigma <- (pc_series >= hoff_lower_interp) & (pc_series <= hoff_upper_interp)
#     coverage <- mean(within_1sigma, na.rm = TRUE)
#     cat("PC", k, "Hoff 1σ coverage:", coverage * 100, "%\n")
# 
#     # Plot
#     p <- ggplot() +
#       geom_ribbon(data = df_sigma, aes(x = time, ymin = lower, ymax = upper, fill = "±1σ Band"), alpha = 0.2) +
#       geom_ribbon(data = df_kernel, aes(x = time, ymin = lower, ymax = upper, fill = "Kernel Envelope"), alpha = 0.25) +
#       geom_ribbon(data = df_hoff, aes(x = time, ymin = lower, ymax = upper, fill = "Hoff Envelope"), alpha = 0.3) +
#       geom_line(data = data.frame(time = time_actual, PC = pc_series), aes(x = time, y = PC, color = "PC Series"), size = 0.8) +
#       scale_fill_manual(name = "Envelope Type", values = c("±1σ Band" = "grey80", "Kernel Envelope" = "green", "Hoff Envelope" = "blue")) +
#       scale_color_manual(name = "PC Line", values = c("PC Series" = "black")) +
#       labs(title = paste("PC", k, "with Variance Envelopes"), x = "Years", y = paste("PC", k)) +
#       theme_minimal(base_size = 14) +
#       theme(legend.position = "top")
# 
#     plot_list[[k]] <- p
#     print(p)
#     Sys.sleep(10)
#   }
# 
#   grid_plot <- grid.arrange(grobs = plot_list, ncol = 1)
#   #print(grid_plot)
# }
# 
# 
# fit_key <- "K_4"
# fit_obj <- fits[[fit_key]]
# plot_list <- list()
# 
# if (!is.null(fit_obj)) {
#   out <- fit_obj$out
#   C <- fit_obj$C
#   D <- fit_obj$D
#   E <- fit_obj$E
#   PC <- out$u[, 1:10]
#   time_actual <- time_veca
# 
#   for (k in 1:3) {
#     cat("→ Plotting PC", k, "with Envelopes...\n")
#     pc_series <- PC[, k]
#     pc_sd <- sd(pc_series, na.rm = TRUE)
#     pc_mean <- mean(pc_series, na.rm = TRUE)
#     cat("PC", k, "series variance:", var(pc_series, na.rm = TRUE), "\n")
# 
#     # ±1σ Band
#     upper_sigma <- rep(pc_mean + pc_sd, length(pc_series))
#     lower_sigma <- rep(pc_mean - pc_sd, length(pc_series))
#     df_sigma <- data.frame(time = time_actual, upper = upper_sigma, lower = lower_sigma)
# 
#     # Kernel Envelope
#     pc_loess <- loess(pc_series ~ time_actual, span = 0.3)
#     pc_trend <- predict(pc_loess, newdata = time_actual)
#     pc_detrended <- pc_series - pc_trend
#     roll_var <- zoo::rollapply(pc_detrended, width = 240, FUN = function(x) var(x, na.rm = TRUE), align = "center", fill = NA)
#     valid_idx <- which(!is.na(roll_var))
#     time_valid <- time_actual[valid_idx]
#     roll_var_valid <- roll_var[valid_idx]
#     roll_var_smooth <- loess(roll_var_valid ~ time_valid, span = 0.75)$fitted
#     roll_var_interp <- zoo::na.approx(zoo(roll_var_smooth, order.by = time_valid), xout = time_actual, na.rm = FALSE)
#     roll_sd_smooth <- sqrt(pmax(roll_var_interp, 1e-6))
#     upper_kernel <- pc_trend + roll_sd_smooth
#     lower_kernel <- pc_trend - roll_sd_smooth
#     df_kernel <- data.frame(time = time_actual, upper = upper_kernel, lower = lower_kernel)
# 
#     # Hoff Quadratic Envelope
#     a_k <- C[k, k]
#     b_k <- 2 * D[k, k]
#     c_k <- E[k, k]
#     t_fine <- seq(min(time_actual), max(time_actual), length.out = 1000)
#     t_scaled <- (t_fine - mean(time_veca)) / (max(time_veca) - min(time_veca)) * 2
# 
#     # Scale c_k selectively to emphasize quadratic term
#     c_scale_factor <- 1000  # Adjusted to enhance parabolic shape
#     var_fine <- a_k + b_k * t_scaled + (c_k * c_scale_factor) * t_scaled^2
# 
#     # Tighter variance cap to prevent over-enclosure
#     var_cap <- var(pc_series, na.rm = TRUE) * 1.5  # Reduced to 1.5x PC variance
#     var_scaled <- pmax(pmin(var_fine, var_cap), 1e-6)
# 
#     # Diagnostic: Print variance range
#     cat("PC", k, "Hoff variance range:", range(var_scaled), "\n")
# 
#     pc_loess <- loess(pc_series ~ time_actual, span = 0.25)
#     pc_fine <- predict(pc_loess, newdata = t_fine)
#     upper_hoff <- pc_fine + sqrt(var_scaled)
#     lower_hoff <- pc_fine - sqrt(var_scaled)
#     df_hoff <- data.frame(time = t_fine, PC = pc_fine, upper = upper_hoff, lower = lower_hoff)
# 
#     # Verify 1σ coverage
#     hoff_upper_interp <- approx(df_hoff$time, df_hoff$upper, xout = time_actual)$y
#     hoff_lower_interp <- approx(df_hoff$time, df_hoff$lower, xout = time_actual)$y
#     within_1sigma <- (pc_series >= hoff_lower_interp) & (pc_series <= hoff_upper_interp)
#     coverage <- mean(within_1sigma, na.rm = TRUE)
#     cat("PC", k, "Hoff 1σ coverage:", coverage * 100, "%\n")
# 
#     # Plot
#     p <- ggplot() +
#       geom_ribbon(data = df_sigma, aes(x = time, ymin = lower, ymax = upper, fill = "±1σ Band"), alpha = 0.2) +
#       geom_ribbon(data = df_kernel, aes(x = time, ymin = lower, ymax = upper, fill = "Kernel Envelope"), alpha = 0.25) +
#       #geom_ribbon(data = df_hoff, aes(x = time, ymin = lower, ymax = upper, fill = "Hoff Envelope"), alpha = 0.3) +
#       geom_line(data = data.frame(time = time_actual, PC = pc_series), aes(x = time, y = PC, color = "PC Series"), size = 0.8) +
#       #geom_line(data = df_sigma, aes(x = time, y = upper, color = "+1σ Band"), size = 0.6)+
#       #geom_line(data = df_sigma, aes(x = time, y = lower, color = "-1σ  Band"), size = 0.6)+
#       scale_fill_manual(name = "Envelope Type", values = c("±1σ Band" = "grey80", "Kernel Envelope" = "green", "Hoff Envelope" = "blue")) +
#       scale_color_manual(name = "PC Line", values = c("PC Series" = "black")) +
#       labs(title = paste("PC", k, "with Variance Envelopes"), x = "Years", y = paste("PC", k)) +
#       theme_minimal(base_size = 14) +
#       theme(legend.position = "top")
#     print(p)
#     plot_list[[k]] <- p
#     Sys.sleep(10)
#   }
# 
#   #grid_plot <- grid.arrange(grobs = plot_list, ncol = 1)
#   #print(grid_plot)
# }


#############  Whole series ######################
# Extract equatorial SST
lat_eq <- which.min(abs(latitude))
S_eq <- t(S[, lat_eq, ])

# Remove monthly climatology
monthly_clim <- tapply(1:nrow(S_eq), months, function(idx) {
  colMeans(S_eq[idx, , drop = FALSE], na.rm = TRUE)
}) %>% do.call(rbind, .)
monthly_clim[is.na(monthly_clim)] <- mean(monthly_clim, na.rm = TRUE)
X1 <- S_eq - monthly_clim[months, ]

# Rename for compatibility with teacher's code
longitude_vec <- longitude
X1 <- X1


# Read in SST anomalies made by subtracting long term monthly means

n=length(time_vec)



# Remove the land columns where SST are all NA
#Â This makes the matrices smaller, speeds up calculations, and allows us to do SVD later
longitude_ocean=longitude_vec[!is.na(X1[1,])] # Select just longitudes that are over ocean
p=length(longitude_ocean)
data=expand.grid(X=longitude_ocean,T=time_vec)
X1a=X1[,!is.na(X1[1,])] # Select columns over ocean
ym=t(X1a)
data$y=as.vector(ym)


# Define the weights for an Epanechnikov (parabolic) kernel
span=21 # 21 years
span=span*12 # convert to months
w=seq(-1,1,length=span)
w=(1-w^2)*3/4 # the parabolic weights
w=w/sum(w) # standardise so weigths sum to one



# Test the kernel out on a univariate time series in the middle of Nino-3.4 region near 140W (220E)
# > longitude_vec[c(41,60)]
# [1] 220.5 239.5
#
y1=ym[70,]
#y1=ym[,1]*(1+(1:n)/n) with increasing variance
plot(time_vec, y1,type="l",#ylim=c(-4,4),
     xlab = "Years",ylab = "Variance",main = "Whole series Variance at 1 sd") # plot the time series
yf=stats::filter(y1,filter=w,method="convolution",sides=2)
lines(time_vec,yf,lwd=2)
a=y1-yf
va=stats::filter(a^2,filter=w,method="convolution",sides=2)
# Note: this approach has less rounding error than estimating v from filter(y1^2)-(filter(y1))^2
lines(time_vec, yf+1*sqrt(va)) # +/- 2sd intervals should include ~95% of the data points
lines(time_vec, yf-1*sqrt(va))
abline(h=0,lty=2)



############################ Standard Deviation for PC's ######################

fit_key <- "K_4"
fit_obj <- fits[[fit_key]]
plot_list <- list()

if (!is.null(fit_obj)) {
  out <- fit_obj$out
  C <- fit_obj$C
  D <- fit_obj$D
  E <- fit_obj$E
  PC <- out$u[, 1:10]
  time_actual <- time_veca
  
  for (k in 1:3) {
    cat("→ Plotting PC", k, "with Envelopes...\n")
    pc_series <- PC[, k]
    pc_sd <- sd(pc_series, na.rm = TRUE)
    pc_mean <- mean(pc_series, na.rm = TRUE)
    cat("PC", k, "series variance:", var(pc_series, na.rm = TRUE), "\n")
    
    # ±1σ Band
    upper_sigma <- rep(pc_mean + pc_sd, length(pc_series))
    lower_sigma <- rep(pc_mean - pc_sd, length(pc_series))
    df_sigma <- data.frame(time = time_actual, upper = upper_sigma, lower = lower_sigma)
    
    # Kernel Envelope
    pc_loess <- loess(pc_series ~ time_actual, span = 0.3)
    pc_trend <- predict(pc_loess, newdata = time_actual)
    pc_detrended <- pc_series - pc_trend
    roll_var <- zoo::rollapply(pc_detrended, width = 240, FUN = function(x) var(x, na.rm = TRUE), align = "center", fill = NA)
    valid_idx <- which(!is.na(roll_var))
    time_valid <- time_actual[valid_idx]
    roll_var_valid <- roll_var[valid_idx]
    roll_var_smooth <- loess(roll_var_valid ~ time_valid, span = 0.75)$fitted
    roll_var_interp <- zoo::na.approx(zoo(roll_var_smooth, order.by = time_valid), xout = time_actual, na.rm = FALSE)
    roll_sd_smooth <- sqrt(pmax(roll_var_interp, 1e-6))
    upper_kernel <- pc_trend + roll_sd_smooth
    lower_kernel <- pc_trend - roll_sd_smooth
    df_kernel <- data.frame(time = time_actual, upper = upper_kernel, lower = lower_kernel)
    
    # Hoff Quadratic Envelope
    a_k <- C[k, k]
    b_k <- 2 * D[k, k]
    c_k <- E[k, k]
    t_fine <- seq(min(time_actual), max(time_actual), length.out = 1000)
    t_scaled <- (t_fine - mean(time_veca)) / (max(time_veca) - min(time_veca)) * 2
    
    # Scale all terms proportionally (e.g., factor of 10)
    scale_factor <- 10
    var_fine <- (a_k) + (b_k ) * t_scaled + (c_k * scale_factor) * t_scaled^2
    
    # Cap variance to prevent over-enclosure
    var_cap <- var(pc_series, na.rm = TRUE) * 2  # 2x PC variance as upper limit
    var_scaled <- pmax(pmin(var_fine, var_cap), 1e-6)
    
    pc_loess <- loess(pc_series ~ time_actual, span = 0.25)
    pc_fine <- predict(pc_loess, newdata = t_fine)
    upper_hoff <- pc_fine + sqrt(var_scaled)
    lower_hoff <- pc_fine - sqrt(var_scaled)
    df_hoff <- data.frame(time = t_fine, PC = pc_fine, upper = upper_hoff, lower = lower_hoff)
    
    # Verify 1σ coverage
    hoff_upper_interp <- approx(df_hoff$time, df_hoff$upper, xout = time_actual)$y
    hoff_lower_interp <- approx(df_hoff$time, df_hoff$lower, xout = time_actual)$y
    within_1sigma <- (pc_series >= hoff_lower_interp) & (pc_series <= hoff_upper_interp)
    coverage <- mean(within_1sigma, na.rm = TRUE)
    cat("PC", k, "Hoff 1σ coverage:", coverage * 100, "%\n")
    
    # Plot
    p <- ggplot() +
      geom_ribbon(data = df_sigma, aes(x = time, ymin = lower, ymax = upper, fill = "±1σ Band"), alpha = 0.2) +
      #geom_ribbon(data = df_kernel, aes(x = time, ymin = lower, ymax = upper, fill = "Kernel Envelope"), alpha = 0.25) +
      #geom_ribbon(data = df_hoff, aes(x = time, ymin = lower, ymax = upper, fill = "Hoff Envelope"), alpha = 0.3) +
      geom_line(data = data.frame(time = time_actual, PC = pc_series), aes(x = time, y = PC, color = "PC Series"), size = 0.8) +
      
      geom_line(data = df_hoff, aes(x = time, y = upper, color = "+1σ Band"), size = 0.6)+
      geom_line(data = df_hoff, aes(x = time, y = lower, color = "-1σ  Band"), size = 0.6)+
      
      scale_fill_manual(name = "Envelope Type", values = c("±1σ Band" = "grey80", "Kernel Envelope" = "green", "Hoff Envelope" = "blue")) +
      scale_color_manual(name = "PC Line", values = c("PC Series" = "black")) +
      labs(title = paste("PC", k, "with Variance Envelopes"), x = "Years", y = paste("PC", k)) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "top")
    
    #plot_list[[k]] <- p
    print(p)
    Sys.sleep(10)
  }
  
  #grid_plot <- grid.arrange(grobs = plot_list, ncol = 1)
  #print(grid_plot)
}


############################ Kernel envelope PC's ######################

fit_key <- "K_4"
fit_obj <- fits[[fit_key]]
plot_list <- list()

if (!is.null(fit_obj)) {
  out <- fit_obj$out
  C <- fit_obj$C
  D <- fit_obj$D
  E <- fit_obj$E
  PC <- out$u[, 1:10]
  time_actual <- time_veca
  
  for (k in 1:3) {
    cat("→ Plotting PC", k, "with Envelopes...\n")
    pc_series <- PC[, k]
    pc_sd <- sd(pc_series, na.rm = TRUE)
    pc_mean <- mean(pc_series, na.rm = TRUE)
    cat("PC", k, "series variance:", var(pc_series, na.rm = TRUE), "\n")
    
    # ±1σ Band
    upper_sigma <- rep(pc_mean + pc_sd, length(pc_series))
    lower_sigma <- rep(pc_mean - pc_sd, length(pc_series))
    df_sigma <- data.frame(time = time_actual, upper = upper_sigma, lower = lower_sigma)
    
    # Kernel Envelope
    pc_loess <- loess(pc_series ~ time_actual, span = 0.3)
    pc_trend <- predict(pc_loess, newdata = time_actual)
    pc_detrended <- pc_series - pc_trend
    roll_var <- zoo::rollapply(pc_detrended, width = 240, FUN = function(x) var(x, na.rm = TRUE), align = "center", fill = NA)
    valid_idx <- which(!is.na(roll_var))
    time_valid <- time_actual[valid_idx]
    roll_var_valid <- roll_var[valid_idx]
    roll_var_smooth <- loess(roll_var_valid ~ time_valid, span = 0.75)$fitted
    roll_var_interp <- zoo::na.approx(zoo(roll_var_smooth, order.by = time_valid), xout = time_actual, na.rm = FALSE)
    roll_sd_smooth <- sqrt(pmax(roll_var_interp, 1e-6))
    upper_kernel <- pc_trend + roll_sd_smooth
    lower_kernel <- pc_trend - roll_sd_smooth
    df_kernel <- data.frame(time = time_actual, upper = upper_kernel, lower = lower_kernel)
    
    # Hoff Quadratic Envelope
    a_k <- C[k, k]
    b_k <- 2 * D[k, k]
    c_k <- E[k, k]
    t_fine <- seq(min(time_actual), max(time_actual), length.out = 1000)
    t_scaled <- (t_fine - mean(time_veca)) / (max(time_veca) - min(time_veca)) * 2
    
    # Scale all terms proportionally (e.g., factor of 10)
    scale_factor <- 1
    var_fine <- (a_k) + (b_k ) * t_scaled + (c_k * scale_factor) * t_scaled^2
    
    # Cap variance to prevent over-enclosure
    var_cap <- var(pc_series, na.rm = TRUE) * 2  # 2x PC variance as upper limit
    var_scaled <- pmax(pmin(var_fine, var_cap), 1e-6)
    
    pc_loess <- loess(pc_series ~ time_actual, span = 0.25)
    pc_fine <- predict(pc_loess, newdata = t_fine)
    upper_hoff <- pc_fine + sqrt(var_scaled)
    lower_hoff <- pc_fine - sqrt(var_scaled)
    df_hoff <- data.frame(time = t_fine, PC = pc_fine, upper = upper_hoff, lower = lower_hoff)
    
    # Verify 1σ coverage
    hoff_upper_interp <- approx(df_hoff$time, df_hoff$upper, xout = time_actual)$y
    hoff_lower_interp <- approx(df_hoff$time, df_hoff$lower, xout = time_actual)$y
    within_1sigma <- (pc_series >= hoff_lower_interp) & (pc_series <= hoff_upper_interp)
    coverage <- mean(within_1sigma, na.rm = TRUE)
    cat("PC", k, "Hoff 1σ coverage:", coverage * 100, "%\n")
    
    # Plot
    p <- ggplot() +
      #geom_ribbon(data = df_sigma, aes(x = time, ymin = lower, ymax = upper, fill = "±1σ Band"), alpha = 0.2) +
      geom_ribbon(data = df_kernel, aes(x = time, ymin = lower, ymax = upper, fill = "Kernel Envelope"), alpha = 0.25) +
      #geom_ribbon(data = df_hoff, aes(x = time, ymin = lower, ymax = upper, fill = "Hoff Envelope"), alpha = 0.3) +
      geom_line(data = data.frame(time = time_actual, PC = pc_series), aes(x = time, y = PC, color = "PC Series"), size = 0.8) +
      
      #geom_line(data = df_hoff, aes(x = time, y = upper, color = "+1σ Band"), size = 0.6)+
      #geom_line(data = df_hoff, aes(x = time, y = lower, color = "-1σ  Band"), size = 0.6)+
      
      scale_fill_manual(name = "Envelope Type", values = c("±1σ Band" = "grey80", "Kernel Envelope" = "green", "Hoff Envelope" = "blue")) +
      scale_color_manual(name = "PC Line", values = c("PC Series" = "black")) +
      labs(title = paste("PC", k, "with Variance Envelopes"), x = "Years", y = paste("PC", k)) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "top")
    
    #plot_list[[k]] <- p
    print(p)
    Sys.sleep(10)
  }
  
  #grid_plot <- grid.arrange(grobs = plot_list, ncol = 1)
  #print(grid_plot)
}

############################ Hoff envelope PC's ######################

fit_key <- "K_4"
fit_obj <- fits[[fit_key]]
plot_list <- list()

if (!is.null(fit_obj)) {
  out <- fit_obj$out
  C <- fit_obj$C
  D <- fit_obj$D
  E <- fit_obj$E
  PC <- out$u[, 1:10]
  time_actual <- time_veca
  
  for (k in 1:3) {
    cat("→ Plotting PC", k, "with Envelopes...\n")
    pc_series <- PC[, k]
    pc_sd <- sd(pc_series, na.rm = TRUE)
    pc_mean <- mean(pc_series, na.rm = TRUE)
    cat("PC", k, "series variance:", var(pc_series, na.rm = TRUE), "\n")
    
    # ±1σ Band
    upper_sigma <- rep(pc_mean + pc_sd, length(pc_series))
    lower_sigma <- rep(pc_mean - pc_sd, length(pc_series))
    df_sigma <- data.frame(time = time_actual, upper = upper_sigma, lower = lower_sigma)
    
    # Kernel Envelope
    pc_loess <- loess(pc_series ~ time_actual, span = 0.3)
    pc_trend <- predict(pc_loess, newdata = time_actual)
    pc_detrended <- pc_series - pc_trend
    roll_var <- zoo::rollapply(pc_detrended, width = 240, FUN = function(x) var(x, na.rm = TRUE), align = "center", fill = NA)
    valid_idx <- which(!is.na(roll_var))
    time_valid <- time_actual[valid_idx]
    roll_var_valid <- roll_var[valid_idx]
    roll_var_smooth <- loess(roll_var_valid ~ time_valid, span = 0.75)$fitted
    roll_var_interp <- zoo::na.approx(zoo(roll_var_smooth, order.by = time_valid), xout = time_actual, na.rm = FALSE)
    roll_sd_smooth <- sqrt(pmax(roll_var_interp, 1e-6))
    upper_kernel <- pc_trend + roll_sd_smooth
    lower_kernel <- pc_trend - roll_sd_smooth
    df_kernel <- data.frame(time = time_actual, upper = upper_kernel, lower = lower_kernel)
    
    # Hoff Quadratic Envelope
    a_k <- C[k, k]
    b_k <- 2 * D[k, k]
    c_k <- E[k, k]
    t_fine <- seq(min(time_actual), max(time_actual), length.out = 1000)
    t_scaled <- (t_fine - mean(time_veca)) / (max(time_veca) - min(time_veca)) * 2
    
    # Scale all terms proportionally (e.g., factor of 10)
    scale_factor <- 200
    var_fine <- (a_k) + (b_k ) * t_scaled + (c_k * scale_factor) * t_scaled^2
    
    # Cap variance to prevent over-enclosure
    var_cap <- var(pc_series, na.rm = TRUE) * 2  # 2x PC variance as upper limit
    var_scaled <- pmax(pmin(var_fine, var_cap), 1e-6)
    
    pc_loess <- loess(pc_series ~ time_actual, span = 0.25)
    pc_fine <- predict(pc_loess, newdata = t_fine)
    upper_hoff <- pc_fine + sqrt(var_scaled)
    lower_hoff <- pc_fine - sqrt(var_scaled)
    df_hoff <- data.frame(time = t_fine, PC = pc_fine, upper = upper_hoff, lower = lower_hoff)
    
    # Verify 1σ coverage
    hoff_upper_interp <- approx(df_hoff$time, df_hoff$upper, xout = time_actual)$y
    hoff_lower_interp <- approx(df_hoff$time, df_hoff$lower, xout = time_actual)$y
    within_1sigma <- (pc_series >= hoff_lower_interp) & (pc_series <= hoff_upper_interp)
    coverage <- mean(within_1sigma, na.rm = TRUE)
    cat("PC", k, "Hoff 1σ coverage:", coverage * 100, "%\n")
    
    # Plot
    p <- ggplot() +
      #geom_ribbon(data = df_sigma, aes(x = time, ymin = lower, ymax = upper, fill = "±1σ Band"), alpha = 0.2) +
      #geom_ribbon(data = df_kernel, aes(x = time, ymin = lower, ymax = upper, fill = "Kernel Envelope"), alpha = 0.25) +
      geom_ribbon(data = df_hoff, aes(x = time, ymin = lower, ymax = upper, fill = "Hoff Envelope"), alpha = 0.3) +
      geom_line(data = data.frame(time = time_actual, PC = pc_series), aes(x = time, y = PC, color = "PC Series"), size = 0.8) +
      
      #geom_line(data = df_hoff, aes(x = time, y = upper, color = "+1σ Band"), size = 0.6)+
      #geom_line(data = df_hoff, aes(x = time, y = lower, color = "-1σ  Band"), size = 0.6)+
      
      scale_fill_manual(name = "Envelope Type", values = c("±1σ Band" = "grey80", "Kernel Envelope" = "green", "Hoff Envelope" = "blue")) +
      scale_color_manual(name = "PC Line", values = c("PC Series" = "black")) +
      labs(title = paste("PC", k, "with Variance Envelopes"), x = "Years", y = paste("PC", k)) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "top")
    
    #plot_list[[k]] <- p
    print(p)
    Sys.sleep(10)
  }
  
  #grid_plot <- grid.arrange(grobs = plot_list, ncol = 1)
  #print(grid_plot)
}

