#' Binomial `distrib` Object
#'
#' @description
#' Creates a distribution object for the Binomial distribution parameterized by the probability of success \eqn{\mu} and a number of trials \eqn{n} (size).
#'
#' @param link_mu A link function object for the mean parameter \eqn{\mu} (probability).
#'   Defaults to \code{\link[linkfunctions]{logit_link}} to ensure the parameter stays within (0, 1).
#' @param size Integer or Numeric Vector. The number of trials \eqn{n}.
#'   Can be a single scalar (default 1) or a vector of the same length as the observations \eqn{y}.
#'
#' @details
#' The Binomial distribution has the following probability mass function (PMF):
#' \deqn{P(Y=y; \mu, n) = \binom{n}{y} \mu^y (1-\mu)^{n-y}}
#' for \eqn{y \in \{0, 1, \dots, n\}}.
#'
#' \strong{Moments:}
#' \itemize{
#'   \item Expected value: \eqn{\mathbb{E}(y) = n\mu}
#'   \item Variance: \eqn{\mathbb{V}(y) = n\mu(1-\mu)}
#'   \item Skewness: \eqn{\gamma_1 = \dfrac{1-2\mu}{\sqrt{n\mu(1-\mu)}}}
#'   \item Excess Kurtosis: \eqn{\gamma_2 = \dfrac{1-6\mu(1-\mu)}{n\mu(1-\mu)}}
#' }
#'
#' \strong{Parameter Domains:}
#' \itemize{
#'   \item \eqn{\mu \in (0, 1)}
#'   \item \eqn{n \in \mathbb{Z}^+} (fixed in constructor, can vary per observation)
#' }
#'
#' \strong{Explicit Formulas:}
#'
#' \emph{Gradient):}
#' The derivative of the log-probability mass function \eqn{\ell} with respect to \eqn{\mu} is:
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y}{\mu} - \dfrac{n-y}{1-\mu} = \dfrac{y - n\mu}{\mu(1-\mu)}}
#'
#' \emph{Expected Hessian:}
#' When \code{expected = TRUE}, substituting \eqn{\mathbb{E}[y]=n\mu}, the expected second derivative of the log-pmf is:
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{n}{\mu(1-\mu)}}
#'
#' \emph{Observed Hessian:}
#' The full analytical second derivative (dependent on observed \eqn{y}) is computed when \code{expected = FALSE}.
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{y}{\mu^2} - \dfrac{n-y}{(1-\mu)^2}}
#'
#' @return A list of class \code{"distrib"} containing the components for the Binomial distribution.
#'
#' @importFrom linkfunctions logit_link
#' @importFrom stats dbinom pbinom qbinom rbinom
#' @export
binomial_distrib <- function(link_mu = logit_link(), size = 1) {
  o <- list()
  class(o) <- c("distrib")

  o$distrib_name <- "binomial"
  o$type <- "discrete"
  o$dimension <- 1

  o$size <- size

  o$bounds <- c(0, max(size))

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
      size = o$size,
      prob = theta[[1]],
      log = log
    )
  }

  o$cdf <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {
    stats::pbinom(
      q = q,
      size = o$size,
      prob = theta[[1]],
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$quantile <- function(p, theta, lower.tail = TRUE, log.p = FALSE) {
    stats::qbinom(
      p = p,
      size = o$size,
      prob = theta[[1]],
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$rng <- function(n, theta) {
    stats::rbinom(
      n = n,
      size = o$size,
      prob = theta[[1]]
    )
  }

  o$loglik <- function(y, theta) {
    o$pdf(y, theta, log = TRUE)
  }

  o$gradient <- function(y, theta) {
    mu <- theta[[1]]
    n_trials <- o$size
    list(
      mu = (y - n_trials * mu) / (mu * (1 - mu))
    )
  }

  o$hessian <- function(y, theta, expected = FALSE) {
    mu <- theta[[1]]
    n_trials <- o$size
    if (expected) {
      list(
        mu_mu = -n_trials / (mu * (1 - mu))
      )
    } else {
      list(
        mu_mu = -y / (mu^2) - (n_trials - y) / ((1 - mu)^2)
      )
    }
  }

  o$kernel <- function(y, theta, log = TRUE) {
    mu <- theta[[1]]
    n_trials <- o$size
    k <- lchoose(n_trials, y) + y * log(mu / (1 - mu)) + n_trials * log(1 - mu)

    if (log) {
      k
    } else {
      exp(k)
    }
  }

  o$normalization_constant <- function(theta, log = TRUE) {
    if (log) {
      0
    } else {
      1
    }
  }

  o$mean <- function(theta) {
    o$size * theta[[1]]
  }

  o$mode <- function(theta) {
    floor((o$size + 1) * theta[[1]])
  }

  o$median <- function(theta) {
    ceiling(o$size * theta[[1]])
  }

  o$variance <- function(theta) {
    mu <- theta[[1]]
    o$size * mu * (1 - mu)
  }

  o$skewness <- function(theta) {
    mu <- theta[[1]]
    (1 - 2 * mu) / sqrt(o$size * mu * (1 - mu))
  }

  o$kurtosis <- function(theta) {
    mu <- theta[[1]]
    (1 - 6 * mu * (1 - mu)) / (o$size * mu * (1 - mu))
  }

  o
}
