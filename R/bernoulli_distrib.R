#' Bernoulli `distrib` Object
#'
#' @description
#' Creates a distribution object for the Bernoulli distribution parameterized by the probability of success \eqn{\mu}.
#'
#' @param link_mu A link function object for the mean parameter \eqn{\mu} (probability).
#'   Defaults to \code{\link[linkfunctions]{logit_link}} to ensure the parameter stays within (0, 1).
#'
#' @details
#' The Bernoulli distribution has the following probability mass function (PMF):
#' \deqn{P(Y=y; \mu) = \mu^y (1-\mu)^{1-y}}
#' for \eqn{y \in \{0, 1\}}.
#'
#' \strong{Moments:}
#' \itemize{
#'   \item Expected value: \eqn{\mathbb{E}(y) = \mu}
#'   \item Variance: \eqn{\mathbb{V}(y) = \mu(1-\mu)}
#'   \item Skewness: \eqn{\gamma_1 = \dfrac{1-2\mu}{\sqrt{\mu(1-\mu)}}}
#'   \item Excess Kurtosis: \eqn{\gamma_2 = \dfrac{1-6\mu(1-\mu)}{\mu(1-\mu)}}
#' }
#'
#' \strong{Parameter Domains:}
#' \itemize{
#'   \item \eqn{\mu \in (0, 1)}
#' }
#'
#' \strong{Explicit Formulas:}
#'
#' \emph{Gradient:}
#' The derivative of the log-probability mass function \eqn{\ell} with respect to \eqn{\mu} is:
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y}{\mu} - \dfrac{1-y}{1-\mu} = \dfrac{y - \mu}{\mu(1-\mu)}}
#'
#' \emph{Expected Hessian:}
#' When \code{expected = TRUE}, the expected second derivative of the log-pmf is:
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\mu(1-\mu)}}
#'
#' \emph{Observed Hessian:}
#' The full analytical second derivative (dependent on observed \eqn{y}) is computed when \code{expected = FALSE}.
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{y}{\mu^2} - \dfrac{1-y}{(1-\mu)^2}}
#'
#' @return A list of class \code{"distrib"} containing the components for the Bernoulli distribution.
#'
#' @importFrom linkfunctions logit_link
#' @importFrom stats dbinom pbinom qbinom rbinom
#' @export
bernoulli_distrib <- function(link_mu = logit_link()) {
  o <- list()
  class(o) <- c("distrib")

  o$distrib_name <- "bernoulli"
  o$type <- "discrete"
  o$dimension <- 1
  o$bounds <- c(0, 1)

  o$params <- c("mu")
  o$params_interpretation <- c(mu = "probability")
  o$n_params <- 1
  o$params_bounds <- list(
    mu = c(0, 1)
  )
  o$link_params <- list(
    mu = link_mu
  )

  o$pdf <- function(y, theta, log = FALSE) {
    stats::dbinom(
      x = y,
      size = 1,
      prob = theta[["mu"]],
      log = log
    )
  }

  o$cdf <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {
    stats::pbinom(
      q = q,
      size = 1,
      prob = theta[["mu"]],
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$qf <- function(p, theta, lower.tail = TRUE, log.p = FALSE) {
    stats::qbinom(
      p = p,
      size = 1,
      prob = theta[["mu"]],
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$rng <- function(n, theta) {
    stats::rbinom(
      n = n,
      size = 1,
      prob = theta[["mu"]]
    )
  }

  o$loglik <- function(y, theta) {
    o$pdf(y, theta, log = TRUE)
  }

  o$gradient <- function(y, theta) {
    mu <- theta[["mu"]]
    list(
      mu = (y - mu) / (mu * (1 - mu))
    )
  }

  o$hessian <- function(y, theta, expected = FALSE) {
    mu <- theta[["mu"]]
    if (expected) {
      list(
        mu_mu = -1 / (mu * (1 - mu))
      )
    } else {
      list(
        mu_mu = -(1 - y) / (1 - mu)^2 - y / mu^2
      )
    }
  }

  o$kernel <- function(y, theta) {
    mu <- theta[["mu"]]
    exp(y * log(mu / (1 - mu)) + log(1 - mu))
  }

  o$normalization_constant <- function(y, theta) {
    1
  }

  o$mean <- function(theta) {
    theta[["mu"]]
  }

  o$variance <- function(theta) {
    mu <- theta[["mu"]]
    mu * (1 - mu)
  }

  o$skewness <- function(theta) {
    mu <- theta[["mu"]]
    (1 - 2 * mu) / sqrt(mu * (1 - mu))
  }

  o$kurtosis <- function(theta) {
    mu <- theta[["mu"]]
    (1 - 6 * mu * (1 - mu)) / (mu * (1 - mu))
  }

  o
}
