#' Poisson `distrib` Object
#'
#' @description
#' Creates a distribution object for the Poisson distribution parameterized by the mean parameter \eqn{\mu}.
#'
#' @param link_mu A link function object for the mean parameter \eqn{\mu}.
#'   Defaults to \code{\link[linkfunctions]{log_link}} to ensure positivity.
#'
#' @details
#' The Poisson distribution has the following probability mass function (PMF):
#' \deqn{P(Y=y; \mu) = \dfrac{\mu^y e^{-\mu}}{y!}}
#' for \eqn{y = 0, 1, 2, \dots}
#'
#' \strong{Moments:}
#' \itemize{
#'   \item Expected value: \eqn{\mathbb{E}(y) = \mu}
#'   \item Variance: \eqn{\mathbb{V}(y) = \mu}
#'   \item Skewness: \eqn{\gamma_1 = \dfrac{1}{\sqrt{\mu}}}
#'   \item Excess Kurtosis: \eqn{\gamma_2 = \dfrac{1}{\mu}}
#' }
#'
#' \strong{Parameter Domains:}
#' \itemize{
#'   \item \eqn{\mu \in (0, +\infty)}
#' }
#'
#' \strong{Explicit Formulas:}
#'
#' \emph{Gradient:}
#' The derivatives of the log-probability mass function \eqn{\ell} with respect to \eqn{\mu} are:
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y}{\mu} - 1 = \dfrac{y - \mu}{\mu}}
#'
#' \emph{Expected Hessian:}
#' When \code{expected = TRUE}, substituting \eqn{\mathbb{E}[y]=\mu}, the expected second derivative of the log-probability mass function is:
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\mu}}
#'
#' \emph{Observed Hessian:}
#' The full analytical second derivative (dependent on observed \eqn{y}) is computed when \code{expected = FALSE}.
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{y}{\mu^2}}
#'
#' @return A list of class \code{"distrib"} containing the components for the Poisson distribution.
#'
#' @importFrom linkfunctions log_link
#' @importFrom stats dpois ppois qpois rpois
#' @export
poisson_distrib <- function(link_mu = log_link()) {
  o <- list()
  class(o) <- c("distrib")

  o$distrib_name <- "poisson"
  o$type <- "discrete"
  o$dimension <- 1
  o$bounds <- c(0, Inf)

  o$params <- c("mu")
  o$params_interpretation <- c(mu = "mean")
  o$n_params <- 1
  o$params_bounds <- list(
    mu = c(0, Inf)
  )
  o$link_params <- list(
    mu = link_mu
  )

  o$pdf <- function(y, theta, log = FALSE) {
    stats::dpois(
      x = y,
      lambda = theta[["mu"]],
      log = log
    )
  }

  o$cdf <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {
    stats::ppois(
      q = q,
      lambda = theta[["mu"]],
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$qf <- function(p, theta, lower.tail = TRUE, log.p = FALSE) {
    stats::qpois(
      p = p,
      lambda = theta[["mu"]],
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$rng <- function(n, theta) {
    stats::rpois(
      n = n,
      lambda = theta[["mu"]]
    )
  }

  o$loglik <- function(y, theta) {
    o$pdf(y, theta, log = TRUE)
  }

  o$gradient <- function(y, theta) {
    mu <- theta[["mu"]]
    list(
      mu = (y - mu) / mu
    )
  }

  o$hessian <- function(y, theta, expected = FALSE) {
    if (expected) {
      list(
        mu_mu = -1 / theta[["mu"]]
      )
    } else {
      list(
        mu_mu = -y / (theta[["mu"]]^2)
      )
    }
  }

  o$kernel <- function(y, theta) {
    exp(y * log(theta[["mu"]]) - lfactorial(y))
  }

  o$normalization_constant <- function(y, theta) {
    exp(theta[["mu"]])
  }

  o$mean <- function(theta) {
    theta[["mu"]]
  }

  o$variance <- function(theta) {
    theta[["mu"]]
  }

  o$skewness <- function(theta) {
    1 / sqrt(theta[["mu"]])
  }

  o$kurtosis <- function(theta) {
    1 / theta[["mu"]]
  }

  o
}
