# pseudohuber_distrib <- function(link_mu = identity_link(), link_sigma = log_link(), link_nu = log_link()) {
#   o <- list()
#   class(o) <- c("distrib")
#
#   o$distrib_name <- "pseudo huber"
#   o$type <- "continuous"
#   o$dimension <- 1
#   o$bounds <- c(-Inf, Inf)
#
#   o$params <- c("mu", "sigma", "nu")
#   o$params_interpretation <- c(mu = "location", sigma = "scale", nu = "shape")
#   o$n_params <- 3
#   o$params_bounds <- list(
#     mu = c(-Inf, Inf),
#     sigma = c(0, Inf),
#     nu = c(0, Inf)
#   )
#   o$link_params <- list(
#     mu = link_mu,
#     sigma = link_sigma,
#     nu = link_nu
#   )
#
#
#   o$kernel <- function(y, theta, log = TRUE) {
#     k <- -sqrt(theta[[3]] + ((y - theta[[1]]) / theta[[2]])^2)
#     if (log) {
#       k
#     } else {
#       exp(k)
#     }
#   }
#
#   o$normalization_constant <- function(theta, log = TRUE) {
#     nu <- theta[[3]]
#     z <- log(2) + log(theta[[2]]) + .5 * log(nu) + log(besselK(sqrt(nu), 1))
#     if (log) {
#       z
#     } else {
#       exp(z)
#     }
#   }
#
#
#   o$pdf <- function(y, theta, log = FALSE) {
#     r <- o$kernel(y, theta, log = TRUE) - o$normalization_constant(theta, log = TRUE)
#     if (log) {
#       r
#     } else {
#       exp(r)
#     }
#   }
#
#   o$cdf <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {
#     distrib:::cdf.distrib(o, q, theta, lower.tail, log.p)
#   }
#
#   o$quantile <- function(p, theta, lower.tail = FALSE, log.p = FALSE) {
#     distrib:::quantile.distrib(o, p, theta, lower.tail, log.p)
#   }
#
#   o$loglik <- function(y, theta) {
#     o$pdf(y, theta, log = TRUE)
#   }
#
#   o$gradient <- function(y, theta) {
#     mu <- theta[[1]]
#     sigma <- theta[[2]]
#     nu <- theta[[3]]
#     res <- y - mu
#     res2 <- res * res
#     sigma2 <- sigma * sigma
#     den <- sqrt(nu + res2 / sigma2)
#     sq_nu <- sqrt(nu)
#     dmu <- res / (sigma2 * sqrt(nu + res2 / sigma2))
#     list(
#       mu = dmu,
#       sigma = (res * dmu - 1) / sigma,
#       nu = -.5 / nu - .5 / den - .5 * (dbesselK(sq_nu, 1) / besselK(sq_nu, 1)) / sq_nu
#     )
#   }
#
#   o$hessian <- function(y, theta, expected = FALSE) {
#     mu <- theta[[1]]
#     sigma <- theta[[2]]
#     nu <- theta[[3]]
#     res <- y - mu
#     res2 <- res * res
#     sigma2 <- sigma * sigma
#     sigma4 <- sigma2 * sigma2
#     den <- sqrt(nu + res2 / sigma2)
#     sq_nu <- sqrt(nu)
#     bk <- besselK(sq_nu, 1, expon.scaled = TRUE)
#     r1 <- dbesselK(sq_nu, 1, deriv = 1, expon.scaled = TRUE, mode = "standard") / bk
#     r2 <- dbesselK(sq_nu, 1, deriv = 2, expon.scaled = TRUE, mode = "standard") / bk
#
#     if (expected) {
#       list(
#         mu_mu = ,
#         sigma_sigma = ,
#         nu_nu = ,
#         mu_sigma = 0,
#         mu_nu = 0,
#         sigma_nu =
#         )
#     } else {
#       list(
#         mu_mu = -nu / (sigma2 * den^3),
#         sigma_sigma = (sigma4 - 3 * sigma2 * res2 / den + res2 * res2 / den^3) / (sigma4 * sigma2),
#         nu_nu = 0.25 / (den^3) + 0.5 / (nu^2) + 0.25 * (r1 * nu^(-1.5) + r1^2 / nu - r2 / nu),
#         mu_sigma = (-2 * nu * sigma2 * res - res^3) / (sigma2 * (nu * sigma2 + res^2)^1.5),
#         mu_nu = (-res / (2 * sigma2)) / (((res2 + nu * sigma2) / sigma2)^1.5),
#         sigma_nu = (-res^2) / (2 * sigma2 * sigma * den^3)
#       )
#     }
#   }
#
#
#
#   o$mean <- o$median <- o$mode <- function(theta) {
#     theta[[1]]
#   }
#
#   o$variance <- function(theta) {
#     sigma <- theta[[2]]
#     sqrt_nu <- sqrt(theta[[3]])
#     ratio_bessel <- besselK(sqrt_nu, 2) / besselK(sqrt_nu, 1)
#     sigma^2 * sqrt_nu * ratio_bessel
#   }
#
#   o$skewness <- function(theta) {
#     0
#   }
#
#   o$kurtosis <- function(theta) {
#     sqrt_nu <- sqrt(theta[[3]])
#     3 * (besselK(sqrt_nu, 3) * besselK(sqrt_nu, 1)) / (besselK(sqrt_nu, 2)^2)
#   }
#
#   o
# }
# #
