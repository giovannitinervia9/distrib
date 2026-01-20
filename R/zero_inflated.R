# library(distrib)
# distrib <- poisson_distrib()
# link_zi <- logit_link()
# zi <- function(distrib, link_zi = logit_link()) {
#   if (!inherits(distrib, "distrib")) {
#     stop("Input must be a 'distrib' object.")
#   }
#
#   o <- distrib
#
#   o$distrib_name <- paste0("zero-inflated_", distrib$distrib_name)
#   o$type <- distrib$type
#   o$dimension <- distrib$dimension
#   if (o$bounds[1] > 0) {
#     o$bounds[1] <- 0
#   }
#   o$params <- c(distrib$params, "zi")
#   o$n_params <- distrib$n_params + 1
#   o$params_interpretation <- c(distrib$params_interpretation, "prob. zero")
#   o$params_bounds$zi <- c(0, 1)
#   o$link_params$zi <- link_zi
#
#   split_theta <- function(theta) {
#     list(
#       orig = theta[1:(o$n_params - 1)],
#       zi = theta[o$n_params]
#     )
#   }
#
#   o$pdf <- function(y, theta, log = FALSE) {
#     pars <- split_theta(theta)
#     th_orig <- pars$orig
#     zi <- pars$zi
#
#     res <- (1 - zi) * distrib$pdf(y, th_orig, log = FALSE)
#
#     is_zero <- y == 0
#     if (any(is_zero)) {
#       res[is_zero] <- zi + res[is_zero]
#     }
#
#     if (log) {
#       log(res)
#     } else {
#       res
#     }
#   }
#
#   o$cdf <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {
#     pars <- split_theta(theta)
#     th_orig <- pars$orig
#     zi <- pars$zi
#     res <- (1 - zi) * distrib$cdf(q, pars$orig, lower.tail = TRUE, log.p = FALSE)
#     res[q >= 0] <- res[q >= 0] + zi
#     if (!lower.tail) {
#       res <- 1 - res
#     }
#     res <- pmin(pmax(res, 0), 1)
#     if (log.p) {
#       log(res)
#     } else {
#       res
#     }
#   }
#
#   o$quantile <- function(p, theta, lower.tail = TRUE, log.p = FALSE) {
#     pars <- split_theta(theta)
#     th_orig <- pars$orig
#     zi <- pars$zi
#     if (log.p) {
#       p <- exp(p)
#     }
#     if (!lower.tail) {
#       p <- 1 - p
#     }
#     p <- pmin(pmax(p, 0), 1)
#
#     F0 <- distrib$cdf(0, pars$orig)
#     p_low <- (1 - zi) * F0
#     p_high <- p_low + zi
#
#     q_vals <- numeric(length(p))
#     p_trans <- numeric(length(p))
#
#     idx_neg <- (p < p_low)
#     idx_pos <- (p > p_high)
#
#     if (any(idx_neg)) {
#       p_trans[idx_neg] <- p[idx_neg] / (1 - zi)
#     }
#
#     if (any(idx_pos)) {
#       p_trans[idx_pos] <- (p[idx_pos] - zi) / (1 - zi)
#     }
#
#     mask_compute <- idx_neg | idx_pos
#
#     if (any(mask_compute)) {
#       p_safe <- pmin(pmax(p_trans[mask_compute], 0), 1)
#       q_vals[mask_compute] <- distrib$quantile(p_safe, pars$orig)
#     }
#
#     q_vals
#   }
#
#   o$rng <- function(n, theta) {
#     pars <- split_theta(theta)
#     zi <- pars$zi
#     is_structural_zero <- stats::runif(n) < zi
#     y <- distrib$rng(n, pars$orig)
#     y[is_structural_zero] <- 0
#     y
#   }
#
#   o$loglik <- function(y, theta) {
#     o$pdf(y, theta, log = TRUE)
#   }
#
#   o$gradient <- function(y, theta, par = NULL) {
#
#   }
#
#   o$hessian <- function(y, theta, expected = FALSE) {
#
#   }
#
#   o$kernel <- function(y, theta, log = TRUE) {
#     o$pdf(y, theta, log = log)
#   }
#
#   o$normalization_constant <- function(theta, log = TRUE) {
#     if (log) 0 else 1
#   }
#
#   o$mean <- function(theta) {
#     pars <- split_theta(theta)
#     (1 - pars$zi) * distrib$mean(pars$orig)
#   }
#
#   o$mode <- function(theta) {
#     # La moda è 0 se P(0) > P(y) per ogni y > 0
#     # Altrimenti è la moda originale (spesso).
#     # Semplificazione euristica:
#     pars <- split_theta(theta)
#     p0 <- o$pdf(0, theta)
#     mode_orig <- distrib$mode(pars$orig)
#     p_mode_orig <- o$pdf(mode_orig, theta)
#     ifelse(p0 > p_mode_orig, 0, mode_orig)
#   }
#
#   o$median <- function(theta) {
#     o$quantile(0.5, theta, lower.tail = TRUE, log.p = FALSE)
#   }
#
#   o$variance <- function(theta) {
#     pars <- split_theta(theta)
#     zi <- pars$zi
#     mu_orig <- distrib$mean(pars$orig)
#     var_orig <- distrib$variance(pars$orig)
#     (1 - zi) * (var_orig + mu_orig^2) - ((1 - zi) * mu_orig)^2
#   }
#
#   o$skewness <- function(theta) {
#     pars <- split_theta(theta)
#     zi <- pars$zi
#
#     mu_x <- distrib$mean(pars$orig)
#     var_x <- distrib$variance(pars$orig)
#     skew_x <- distrib$skewness(pars$orig)
#     m3_raw_x <- skew_x * var_x^1.5 + 3 * mu_x * var_x + mu_x^3
#
#     mu_zi <- (1 - zi) * mu_x
#     var_zi <- (1 - zi) * (var_x + mu_x^2) - mu_zi^2
#     m3_raw_zi <- (1 - zi) * m3_raw_x
#
#     num <- m3_raw_zi - 3 * mu_zi * var_zi - mu_zi^3
#     den <- var_zi^1.5
#
#     ifelse(den == 0, 0, num / den)
#   }
#
#   o$kurtosis <- function(theta) {
#     pars <- split_theta(theta)
#     zi <- pars$zi
#
#     mu_x <- distrib$mean(pars$orig)
#     var_x <- distrib$variance(pars$orig)
#     skew_x <- distrib$skewness(pars$orig)
#     kurt_excess_x <- distrib$kurtosis(pars$orig)
#
#     m2_x <- var_x + mu_x^2
#     m3_x <- skew_x * var_x^1.5 + 3 * mu_x * var_x + mu_x^3
#     m4_cen_x <- (kurt_excess_x + 3) * var_x^2
#     m4_x <- m4_cen_x + 4 * mu_x * m3_x - 6 * mu_x^2 * m2_x + 3 * mu_x^4
#
#     mu_zi <- (1 - zi) * mu_x
#     var_zi <- (1 - zi) * m2_x - mu_zi^2
#
#     m4_raw_zi <- (1 - zi) * m4_x
#     m3_raw_zi <- (1 - zi) * m3_x
#     m2_raw_zi <- (1 - zi) * m2_x
#
#     m4_cen_zi <- m4_raw_zi - 4 * mu_zi * m3_raw_zi + 6 * mu_zi^2 * m2_raw_zi - 3 * mu_zi^4
#
#     den <- var_zi^2
#     ifelse(den == 0, 0, (m4_cen_zi / den) - 3)
#   }
#
#   o$initialize <- function(y) {
#     res <- distrib$initialize(y[y != 0])
#     res$zi <- min(max(mean(y == 0), 0.01), 0.99)
#     res
#   }
# }
