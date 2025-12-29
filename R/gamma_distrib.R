#' Gamma `distrib` Object (Mean-Variance Parameterization)
#'
#' @description
#' Creates a distribution object for the Gamma distribution parameterized by mean (\eqn{\mu}) and variance (\eqn{\sigma^2}).
#'
#' @param link_mu A link function object for the mean parameter \eqn{\mu}.
#'   Defaults to \code{\link[linkfunctions]{log_link}} to ensure positivity.
#' @param link_sigma2 A link function object for the variance parameter \eqn{\sigma^2}.
#'   Defaults to \code{\link[linkfunctions]{log_link}} to ensure positivity.
#'
#' @details
#' The Gamma distribution is reparameterized from the standard shape \eqn{(\alpha)} and rate \eqn{(\lambda)} parameters using:
#' \deqn{\alpha = \dfrac{\mu^2}{\sigma^2}, \quad \lambda = \dfrac{\mu}{\sigma^2}}
#'
#' The probability density function is:
#' \deqn{f(y; \mu, \sigma^2) = \dfrac{1}{\Gamma\left(\dfrac{\mu^2}{\sigma^2}\right)} \left(\dfrac{\mu}{\sigma^2}\right)^{\dfrac{\mu^2}{\sigma^2}} y^{\dfrac{\mu^2}{\sigma^2}-1} \exp\left(-\dfrac{\mu}{\sigma^2} y\right)}
#' for \eqn{y > 0}.
#'
#' \strong{Moments:}
#' \itemize{
#'   \item Expected value: \deqn{\mathbb{E}(y) = \mu}
#'   \item Variance: \deqn{\mathbb{V}(y) = \sigma^2}
#'   \item Skewness: \deqn{\gamma_1 = \dfrac{2\sqrt{\sigma^2}}{\mu}}
#'   \item Excess Kurtosis: \deqn{\gamma_2 = \dfrac{6\sigma^2}{\mu^2}}
#' }
#'
#' \strong{Parameter Domains:}
#' \itemize{
#'   \item \eqn{\mu \in (0, +\infty)}
#'   \item \eqn{\sigma^2 \in (0, +\infty)}
#' }
#'
#' \strong{Explicit Formulas:}
#'
#' \emph{Gradient:}
#'
#' Let \eqn{\psi(z)} be the digamma function evaluated at \eqn{z = \mu^2/\sigma^2}. The derivatives of the log-pdf \eqn{\ell} are:
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{-2\mu\psi\left(\dfrac{\mu^2}{\sigma^2}\right) + 2\mu\log\left(\dfrac{\mu}{\sigma^2}\right) + \mu + 2\mu\log(y) - y}{\sigma^2}}
#' \deqn{\dfrac{\partial \ell}{\partial \sigma^2} = -\dfrac{\mu\left[-\mu\psi\left(\dfrac{\mu^2}{\sigma^2}\right) + \mu + \mu\left(\log\left(\dfrac{\mu}{\sigma^2}\right) + \log(y)\right) - y\right]}{(\sigma^2)^2}}
#'
#' \emph{Expected Hessian:}
#'
#' When \code{expected = TRUE}, substituting \eqn{\mathbb{E}[y]=\mu} and \eqn{\mathbb{E}[\log(y)] = \psi(\mu^2/\sigma^2) - \log(\mu/\sigma^2)} and letting \eqn{\psi_1(z)} be the trigamma function evaluated at \eqn{z = \mu^2/\sigma^2}, the expected second derivatives of the log-pdf are:
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = \dfrac{3\sigma^2 - 4\mu^2\psi_1\left(\dfrac{\mu^2}{\sigma^2}\right)}{(\sigma^2)^2}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial (\sigma^2)^2}\right] = -\dfrac{\mu^2\left(\mu^2\psi_1\left(\dfrac{\mu^2}{\sigma^2}\right) - \sigma^2\right)}{(\sigma^2)^4}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma^2}\right] = \dfrac{2\mu\left(\mu^2\psi_1\left(\dfrac{\mu^2}{\sigma^2}\right) - \sigma^2\right)}{(\sigma^2)^3}}
#'
#' \emph{Observed Hessian:}
#'
#' When \code{expected = FALSE}, the full analytical second derivatives depend on observed \eqn{y}.
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{1}{\sigma^2}\left[ -\dfrac{4\mu^2}{\sigma^2}\psi_1\left(\dfrac{\mu^2}{\sigma^2}\right) - 2\psi(\dfrac{\mu^2}{\sigma^2}) + 2\log\left(\dfrac{\mu}{\sigma^2}\right) + 2\log(y) + 3 \right]}
#' \deqn{\dfrac{\partial^2 \ell}{\partial (\sigma^2)^2} = -\dfrac{\mu}{(\sigma^2)^4} \left[ 2\mu\sigma^2\psi(\dfrac{\mu^2}{\sigma^2}) + \mu^3\psi_1\left(\dfrac{\mu^2}{\sigma^2}\right) + \sigma^2\left( -2\mu\log\left(\dfrac{\mu}{\sigma^2}\right) - 3\mu - 2\mu\log(y) + 2y \right) \right]}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma^2} = \dfrac{1}{(\sigma^2)^3} \left[ 2\mu\sigma^2\psi(\dfrac{\mu^2}{\sigma^2}) + 2\mu^3\psi_1\left(\dfrac{\mu^2}{\sigma^2}\right) + \sigma^2\left( -2\mu\log\left(\dfrac{\mu}{\sigma^2}\right) - 3\mu - 2\mu\log(y) + y \right) \right]}
#'
#' @return A list of class \code{"distrib"} containing the components for the Gamma distribution.
#'
#' @importFrom linkfunctions log_link
#' @importFrom stats dgamma pgamma qgamma rgamma
#' @export
gamma_distrib <- function(link_mu = log_link(), link_sigma2 = log_link()) {
  o <- list()
  class(o) <- c("distrib")

  o$distrib_name <- "gamma"
  o$type <- "continuous"
  o$dimension <- 1
  o$bounds <- c(0, Inf)

  o$params <- c("mu", "sigma2")
  o$params_interpretation <- c(mu = "mean", sigma2 = "variance")
  o$n_params <- 2
  o$params_bounds <- list(
    mu = c(0, Inf),
    sigma2 = c(0, Inf)
  )
  o$link_params <- list(
    mu = link_mu,
    sigma2 = link_sigma2
  )

  o$pdf <- function(y, theta, log = FALSE) {
    mu <- theta[["mu"]]
    sigma2 <- theta[["sigma2"]]
    stats::dgamma(
      x = y,
      shape = mu^2 / sigma2,
      rate = mu / sigma2,
      log = log
    )
  }

  o$cdf <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {
    mu <- theta[["mu"]]
    sigma2 <- theta[["sigma2"]]
    stats::pgamma(
      q = q,
      shape = mu^2 / sigma2,
      rate = mu / sigma2,
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$qf <- function(p, theta, lower.tail = TRUE, log.p = FALSE) {
    mu <- theta[["mu"]]
    sigma2 <- theta[["sigma2"]]
    stats::qgamma(
      p = p,
      shape = mu^2 / sigma2,
      rate = mu / sigma2,
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$rng <- function(n, theta) {
    mu <- theta[["mu"]]
    sigma2 <- theta[["sigma2"]]
    stats::rgamma(
      n = n,
      shape = mu^2 / sigma2,
      rate = mu / sigma2
    )
  }

  o$loglik <- function(y, theta) {
    o$pdf(y, theta, log = TRUE)
  }

  o$gradient <- function(y, theta) {
    mu <- theta[["mu"]]
    sigma2 <- theta[["sigma2"]]
    mu_sigma2 <- mu / sigma2
    mu2_sigma2 <- mu_sigma2 * mu
    digamma_mu2_sigma2 <- digamma(mu2_sigma2)

    list(
      mu = (-2 * mu * digamma_mu2_sigma2 + 2 * mu * log(mu_sigma2) + mu + 2 * mu * log(y) - y) / sigma2,
      sigma2 = -(mu * (-mu * digamma_mu2_sigma2 + mu + mu * (log(mu_sigma2) + log(y)) - y)) / sigma2^2
    )
  }

  o$hessian <- function(y, theta, expected = FALSE) {
    mu <- theta[["mu"]]
    sigma2 <- theta[["sigma2"]]
    mu_sigma2 <- mu / sigma2
    mu2_sigma2 <- mu_sigma2 * mu
    digamma_mu2_sigma2 <- digamma(mu2_sigma2)
    trigamma_mu2_sigma2 <- trigamma(mu2_sigma2)
    log_y <- log(y)

    if (expected) {
      list(
        mu_mu = (3 * sigma2 - 4 * mu^2 * trigamma_mu2_sigma2) / sigma2^2,
        sigma2_sigma2 = -mu^2 * (mu^2 * trigamma_mu2_sigma2 - sigma2) / sigma2^4,
        mu_sigma2 = 2 * mu * (mu^2 * trigamma_mu2_sigma2 - sigma2) / sigma2^3
      )
    } else {
      list(
        mu_mu = (-(4 * mu^2 * trigamma_mu2_sigma2) / sigma2 - 2 * digamma_mu2_sigma2 + 2 * log(mu_sigma2) + 2 * log_y + 3) / sigma2,
        sigma2_sigma2 = -(mu * (2 * mu * sigma2 * digamma_mu2_sigma2 + mu^3 * trigamma_mu2_sigma2 + sigma2 * (-2 * mu * log(mu_sigma2) - 3 * mu - 2 * mu * log_y + 2 * y))) / sigma2^4,
        mu_sigma2 = (2 * mu * sigma2 * digamma_mu2_sigma2 + 2 * mu^3 * trigamma_mu2_sigma2 + sigma2 * (-2 * mu * log(mu_sigma2) - 3 * mu - 2 * mu * log_y + y)) / sigma2^3
      )
    }
  }

  o$kernel <- function(y, theta) {
    mu <- theta[["mu"]]
    sigma2 <- theta[["sigma2"]]
    mu_sigma2 <- mu / sigma2
    mu2_sigma2 <- mu_sigma2 * mu
    exp((mu2_sigma2 - 1) * log(y) - y * mu_sigma2)
  }

  o$normalization_constant <- function(y, theta) {
    mu <- theta[["mu"]]
    sigma2 <- theta[["sigma2"]]
    mu_sigma2 <- mu / sigma2
    mu2_sigma2 <- mu_sigma2 * mu
    gamma(mu2_sigma2) / ((mu_sigma2)^(mu2_sigma2))
  }

  o$mean <- function(theta) {
    theta[["mu"]]
  }

  o$variance <- function(theta) {
    theta[["sigma2"]]
  }

  o$skewness <- function(theta) {
    2 * sqrt(theta[["sigma2"]]) / theta[["mu"]]
  }

  o$kurtosis <- function(theta) {
    6 * theta[["sigma2"]] / theta[["mu"]]^2
  }

  o
}
