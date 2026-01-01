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
#     cdf.distrib(o, q, theta, lower.tail, log.p)
#   }
#
#   o$qf <- function(p, theta, lower.tail = FALSE, log.p = FALSE) {
#     if (log.p) {
#       p <- exp(p)
#     }
#
#     if (lower.tail) {
#       p <- 1 - p
#     }
#
#     newton_step <- function(y, p, theta) {
#       (o$cdf(y, theta) - p) / o$pdf(y, theta)
#     }
#   }
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
