gaussian_distrib <- function(link_mu = identity_link(), link_sigma = log_link()) {
  o <- list()
  class(o) <- c("distrib")

  o$distrib_name <- "gaussian"
  o$type <- "continuous"
  o$dimension <- 1
  o$bounds <- c(-Inf, Inf)

  o$params <- c("mu", "sigma")
  o$n_params <- 2
  o$params_bounds <- list(
    mu = c(-Inf, Inf),
    sigma = c(0, Inf)
  )
  o$link_params <- list(
    mu = link_mu,
    sigma = link_sigma
  )

  o$pdf <- function(y, theta, log = FALSE) {
    stats::dnorm(
      x = y,
      mean = theta[["mu"]],
      sd = theta[["sigma"]],
      log = log
    )
  }

  o$cdf <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {
    stats::pnorm(
      q = q,
      mean = theta[["mu"]],
      sd = theta[["sigma"]],
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$qf <- function(p, theta, lower.tail = TRUE, log.p = FALSE) {
    stats::qnorm(
      p = p,
      mean = theta[["mu"]],
      sd = theta[["sigma"]],
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$rng <- function(n, theta) {
    stats::rnorm(
      n = n,
      mean = theta[["mu"]],
      sd = theta[["sigma"]]
    )
  }

  o$loglik <- function(y, theta) {
    o$pdf(y, theta, log = TRUE)
  }

  o$gradient <- function(y, theta) {
    sigma <- theta[["sigma"]]
    sigma2 <- sigma^2
    residuals <- y - theta[["mu"]]
    list(
      mu = residuals / sigma2,
      sigma = (residuals^2 - sigma2) / sigma^3
    )
  }

  o$hessian <- function(y, theta, expected = FALSE) {
    sigma <- theta[["sigma"]]
    sigma2 <- sigma^2
    if (expected) {
      list(
        mu_mu = -1 / sigma2,
        sigma_sigma = -2 / sigma2,
        mu_sigma = 0
      )
    } else {
      residuals <- y - theta[["mu"]]
      list(
        mu_mu = -1 / sigma2,
        sigma_sigma = (sigma2 - 3 * residuals^2) / (sigma2 * sigma2),
        mu_sigma = -2 * residuals / (sigma2 * sigma)
      )
    }
  }

  o$kernel <- function(y, theta) {
    exp(-.5 * ((y - theta[["mu"]]) / theta[["sigma"]])^2)
  }

  o$normalization_constant <- function(y, theta) {
    sqrt(2 * pi) * theta[["sigma"]]
  }

  invisible(o)
}
