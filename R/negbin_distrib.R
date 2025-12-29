#' Negative Binomial `distrib` Object (NB2)
#'
#' @description
#' Creates a distribution object for the Negative Binomial distribution (NB2)
#' parameterized by mean (\eqn{\mu}) and dispersion (\eqn{\theta}).
#'
#' @param link_mu A link function object for the mean parameter \eqn{\mu}.
#'   Defaults to \code{\link[linkfunctions]{log_link}}.
#' @param link_theta A link function object for the dispersion parameter \eqn{\theta}.
#'   Defaults to \code{\link[linkfunctions]{log_link}}.
#'
#' @details
#' The PMF is: \deqn{P(Y=y; \mu, \theta) = \dfrac{\Gamma(y+\theta)}{y!\Gamma(\theta)} \left(\dfrac{\theta}{\theta+\mu}\right)^\theta \left(\dfrac{\mu}{\theta+\mu}\right)^y}
#'
#' \strong{Moments:}
#' \itemize{
#'   \item Expected value: \eqn{\mathbb{E}(y) = \mu}
#'   \item Variance: \eqn{\mathbb{V}(y) = \mu + \dfrac{\mu^2}{\theta}}
#'   \item Skewness: \eqn{\gamma_1 = \dfrac{\theta + 2\mu}{\sqrt{\mu\theta(\theta + \mu)}}}
#'   \item Excess Kurtosis: \eqn{\gamma_2 = \dfrac{6}{\theta} + \dfrac{(\theta + \mu)^2}{\mu\theta(\theta + \mu)}}
#' }
#'
#' @return A list of class \code{"distrib"}.
#' @export
negbin_distrib <- function(link_mu = log_link(), link_theta = log_link()) {
  o <- list()
  class(o) <- c("distrib")

  o$distrib_name <- "negative binomial"
  o$type <- "discrete"
  o$dimension <- 1
  o$bounds <- c(0, Inf)

  o$params <- c("mu", "theta")
  o$params_interpretation <- c(mu = "mean", theta = "dispersion")
  o$n_params <- 2
  o$params_bounds <- list(mu = c(0, Inf), theta = c(0, Inf))
  o$link_params <- list(mu = link_mu, theta = link_theta)

  o$pdf <- function(y, theta, log = FALSE) {
    stats::dnbinom(
      x = y,
      mu = theta[["mu"]],
      size = theta[["theta"]],
      log = log)
  }

  o$cdf <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {
    stats::pnbinom(
      q = q,
      mu = theta[["mu"]],
      size = theta[["theta"]],
      lower.tail = lower.tail,
      log.p = log.p)
  }

  o$qf <- function(p, theta, lower.tail = TRUE, log.p = FALSE) {
    stats::qnbinom(
      p = p,
      mu = theta[["mu"]],
      size = theta[["theta"]],
      lower.tail = lower.tail,
      log.p = log.p)
  }

  o$rng <- function(n, theta) {
    stats::rnbinom(
      n = n,
      mu = theta[["mu"]],
      size = theta[["theta"]])
  }

  o$loglik <- function(y, theta) {
    o$pdf(y, theta, log = TRUE)
  }

  o$gradient <- function(y, theta) {
    mu <- theta[["mu"]]
    th <- theta[["theta"]]
    list(
      mu = th * (y - mu) / (mu * (mu + th)),
      theta = digamma(y + th) - digamma(th) + log(th / (mu + th)) + (mu - y) / (mu + th)
    )
  }

  o$hessian <- function(y, theta, expected = FALSE) {
    mu <- theta[["mu"]]
    th <- theta[["theta"]]

    if (expected) {
      list(
        mu_mu = -th / (mu * (mu + th)),
        theta_theta = -trigamma(th) + 1/th - 1/(mu + th),
        mu_theta = -1 / (mu + th)
      )
    } else {
      # Hessiana osservata (seconda derivata della log-verosimiglianza)
      list(
        mu_mu = -th * (y / mu^2 + th / (mu + th)^2),
        theta_theta = -trigamma(y + th) + trigamma(th) - 1/th + 1/(mu + th) + (mu - y) / (mu + th)^2,
        mu_theta = -(y - mu) / (mu + th)^2
      )
    }
  }

  o$mean <- function(theta) {
    theta[["mu"]]
  }

  o$variance <- function(theta) {
    mu <- theta[["mu"]]
    th <- theta[["theta"]]
    mu + (mu^2 / th)
  }

  o$skewness <- function(theta) {
    mu <- theta[["mu"]]
    th <- theta[["theta"]]
    (th + 2*mu) / sqrt(mu * th * (th + mu))
  }

  o$kurtosis <- function(theta) {
    mu <- theta[["mu"]]
    th <- theta[["theta"]]
    6/th + (th + mu)^2 / (mu * th * (th + mu))
  }

  o
}


# rm(list = ls())
# gc()
# y <- rpois(1, 5)
# mu <- rpois(1, 5)
# theta <- rexp(1, .1)
#
#
# # lpdf
# log(Gamma(y+theta))-log(Gamma(theta))+theta*log(theta/(theta+mu))+y*log(mu/(theta+mu))
#
# # constants
# th_plus_mu <- theta + mu
# th_frac_th_plus_mu <- theta/th_plus_mu
#
#
# # d_mu
# th_frac_th_plus_mu*(y/mu - 1)
#
# # d_theta
# -digamma(theta)
#
