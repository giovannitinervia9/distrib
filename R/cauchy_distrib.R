#' Cauchy `distrib` Object
#'
#' @description
#' Creates a distribution object for the Cauchy distribution, parameterized by location (\eqn{\mu})
#' and scale (\eqn{\sigma}).
#'
#' @param link_mu A link function object for the location parameter \eqn{\mu}.
#'   Defaults to \code{\link[linkfunctions]{identity_link}}.
#' @param link_sigma A link function object for the scale parameter \eqn{\sigma}.
#'   Defaults to \code{\link[linkfunctions]{log_link}} to ensure positivity.
#'
#' @details
#' The Cauchy distribution is a continuous probability distribution with heavy tails and no
#' defined mean or variance.
#'
#' \strong{Density function:}
#' \deqn{f(y; \mu, \sigma) = \dfrac{1}{\pi \sigma \left[1 + \left(\dfrac{y-\mu}{\sigma}\right)^2\right]}}
#'
#' \strong{Moments:}
#' \itemize{
#'   \item Expected value: \eqn{\mathbb{E}(y) = \text{Undefined}}
#'   \item Variance: \eqn{\mathbb{V}(y) = \text{Undefined}}
#'   \item Skewness: \eqn{\gamma_1 = \text{Undefined}}
#'   \item Excess Kurtosis: \eqn{\gamma_2 = \text{Undefined}}
#'   \item Mode and Median: \eqn{\mu}
#' }
#'
#' \strong{Explicit Derivatives:}
#'
#' \emph{Gradient:}
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{2(y-\mu)}{\sigma^2 + (y-\mu)^2}}
#' \deqn{\dfrac{\partial \ell}{\partial \sigma} = \dfrac{(y-\mu)^2 - \sigma^2}{\sigma(\sigma^2 + (y-\mu)^2)}}
#'
#' \emph{Observed Hessian:}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{2(y-\mu)^2 - 2\sigma^2}{(\sigma^2 + (y-\mu)^2)^2}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{\sigma^4 - 4\sigma^2 (y-\mu)^2 - (y-\mu)^4}{\sigma^2(\sigma^2 + (y-\mu)^2)^2}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma} = -\dfrac{4\sigma (y-\mu)}{(\sigma^2 + (y-\mu)^2)^2}}
#'
#' \emph{Expected Hessian:}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{2\sigma^2}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right] = -\dfrac{1}{2\sigma^2}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma}\right] = 0}
#'
#' @return A list of class \code{"distrib"} containing the components for the Cauchy distribution.
#'
#' @importFrom stats dcauchy pcauchy qcauchy rcauchy median IQR
#' @export
cauchy_distrib <- function(link_mu = identity_link(), link_sigma = log_link()) {
  o <- list()
  class(o) <- c("distrib")

  o$distrib_name <- "cauchy"
  o$type <- "continuous"
  o$dimension <- 1
  o$bounds <- c(-Inf, Inf)

  o$params <- c("mu", "sigma")
  o$params_interpretation <- c(mu = "location", sigma = "scale")
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
    stats::dcauchy(
      x = y,
      location = theta[[1]],
      scale = theta[[2]],
      log = log
    )
  }

  o$cdf <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {
    stats::pcauchy(
      q = q,
      location = theta[[1]],
      scale = theta[[2]],
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$quantile <- function(p, theta, lower.tail = TRUE, log.p = FALSE) {
    stats::qcauchy(
      p = p,
      location = theta[[1]],
      scale = theta[[2]],
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$rng <- function(n, theta) {
    stats::rcauchy(
      n = n,
      location = theta[[1]],
      scale = theta[[2]]
    )
  }

  o$loglik <- function(y, theta) {
    o$pdf(y, theta, log = TRUE)
  }

  o$gradient <- function(y, theta, par = NULL) {
    if (is.null(par)) {
      par <- o$params
    }

    invalid_pars <- setdiff(par, o$params)
    if (length(invalid_pars) > 0) {
      warning(sprintf(
        "Invalid parameter(s) specified: %s. Available parameters are: %s.",
        paste(sQuote(invalid_pars), collapse = ", "),
        paste(sQuote(o$params), collapse = ", ")
      ))
    }

    sigma <- theta[[2]]
    res <- y - theta[[1]]
    den <- sigma^2 + res^2
    g <- list()
    if ("mu" %in% par) {
      g$mu <- 2 * res / den
    }
    if ("sigma" %in% par) {
      g$sigma <- (res^2 - sigma^2) / (sigma * den)
    }
    g
  }

  o$hessian <- function(y, theta, expected = FALSE) {
    sigma <- theta[[2]]
    if (expected) {
      val <- -.5 / sigma^2
      list(
        mu_mu = val,
        sigma_sigma = val,
        mu_sigma = 0
      )
    } else {
      res <- y - theta[[1]]
      res2 <- res^2
      sigma2 <- sigma * sigma
      den2 <- (sigma2 + res2)^2

      list(
        mu_mu = (2 * res2 - 2 * sigma2) / den2,
        sigma_sigma = (sigma2 * sigma2 - 4 * sigma2 * res2 - res2 * res2) / (sigma2 * den2),
        mu_sigma = -4 * sigma * res / (den2)
      )
    }
  }

  o$kernel <- function(y, theta, log = TRUE) {
    k <- -log(1 + ((y - theta[[1]]) / theta[[2]])^2)
    if (log) {
      k
    } else {
      exp(k)
    }
  }

  o$normalization_constant <- function(theta, log = TRUE) {
    z <- pi * theta[[2]]

    if (log) {
      log(z)
    } else {
      z
    }
  }

  o$mean <- function(theta) {
    NA
  }
  o$median <- o$mode <- function(theta) {
    theta[[1]]
  }

  o$variance <- function(theta) {
    NA
  }

  o$skewness <- function(theta) {
    NA
  }

  o$kurtosis <- function(theta) {
    NA
  }

  o$initialize <- function(y) {
    list(
      mu = median(y),
      sigma = .5 * IQR(y)
    )
  }

  o
}
