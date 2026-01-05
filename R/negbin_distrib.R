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
#' @importFrom stats dnbinom pnbinom qnbinom rnbinom
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
      mu = theta[[1]],
      size = theta[[2]],
      log = log
    )
  }

  o$cdf <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {
    stats::pnbinom(
      q = q,
      mu = theta[[1]],
      size = theta[[2]],
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$quantile <- function(p, theta, lower.tail = TRUE, log.p = FALSE) {
    stats::qnbinom(
      p = p,
      mu = theta[[1]],
      size = theta[[2]],
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$rng <- function(n, theta) {
    stats::rnbinom(
      n = n,
      mu = theta[[1]],
      size = theta[[2]]
    )
  }

  o$loglik <- function(y, theta) {
    o$pdf(y, theta, log = TRUE)
  }

  o$gradient <- function(y, theta) {
    mu <- theta[[1]]
    theta <- theta[[2]]
    y_plus_theta <- y + theta
    th_plus_mu <- theta + mu
    th_frac_th_plus_mu <- theta / th_plus_mu
    list(
      mu = th_frac_th_plus_mu * (y / mu - 1),
      theta = -digamma(theta) + digamma(y_plus_theta) + log(th_frac_th_plus_mu) - (y - mu) / (th_plus_mu)
    )
  }

  # function to compute expectation of trigamma(y + theta) for expected hessian
  o$E_trigamma <- function(mu, theta, p = .999) {
    mapply(\(mu, theta){
      Y <- 0:pmax(100, qnbinom(p = p, size = theta, mu = mu))
      sum(trigamma(Y + theta) * dnbinom(Y, size = theta, mu = mu))
    }, mu = mu, theta = theta)
  }

  o$hessian <- function(y, theta, expected = FALSE) {
    mu <- theta[[1]]
    theta <- theta[[2]]
    y_plus_theta <- y + theta
    th_plus_mu <- theta + mu
    th_frac_th_plus_mu <- theta / th_plus_mu

    if (expected) {
      list(
        mu_mu = -theta / (mu * (th_plus_mu)),
        theta_theta = mu / (theta * (th_plus_mu)) - trigamma(theta) + o$E_trigamma(mu, theta),
        mu_theta = 0
      )
    } else {
      list(
        mu_mu = y_plus_theta / (th_plus_mu)^2 - y / mu^2,
        theta_theta = (y * theta + mu^2) / (theta * th_plus_mu^2) - trigamma(theta) + trigamma(y_plus_theta),
        mu_theta = (y - mu) / (th_plus_mu)^2
      )
    }
  }

  o$mean <- function(theta) {
    theta[[1]]
  }

  o$variance <- function(theta) {
    mu <- theta[[1]]
    th <- theta[[2]]
    mu + (mu^2 / th)
  }

  o$skewness <- function(theta) {
    mu <- theta[[1]]
    th <- theta[[2]]
    (th + 2 * mu) / sqrt(mu * th * (th + mu))
  }

  o$kurtosis <- function(theta) {
    mu <- theta[[1]]
    th <- theta[[2]]
    6 / th + (th + mu)^2 / (mu * th * (th + mu))
  }

  o$kernel <- function(y, theta, log = TRUE) {
    mu <- theta[[1]]
    th <- theta[[2]]
    k <- lgamma(y + th) - lfactorial(y) + y * log(mu / (mu + th))
    if (log) {
      k
    } else {
      exp(k)
    }
  }

  o$normalization_constant <- function(theta, log = TRUE) {
    mu <- theta[[1]]
    th <- theta[[2]]
    z <- lgamma(th) + th * log((mu + th) / th)

    if (log) {
      z
    } else {
      exp(z)
    }
  }

  o$median <- function(theta) {
    stats::qnbinom(
      p = 0.5,
      mu = theta[[1]],
      size = theta[[2]]
    )
  }

  o$mode <- function(theta) {
    mu <- theta[[1]]
    th <- theta[[2]]
    ifelse(th > 1, floor(mu * (th - 1) / th), 0)
  }

  o
}
