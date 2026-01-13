#' Lognormal `distrib` Object (Log-Scale Parameterization)
#'
#' @description
#' Creates a distribution object for the Lognormal distribution parameterized by the
#' mean (\eqn{\mu}) and the variance (\eqn{\sigma^2}) of the log-transformed variable.
#'
#' @param link_mu A link function object for the log-mean parameter \eqn{\mu}.
#'   Defaults to \code{\link[linkfunctions]{identity_link}}.
#' @param link_sigma2 A link function object for the log-variance parameter \eqn{\sigma^2}.
#'   Defaults to \code{\link[linkfunctions]{log_link}} to ensure positivity.
#'
#' @details
#' The Lognormal distribution has the following probability density function:
#' \deqn{f(y; \mu, \sigma^2) = \dfrac{1}{y\sqrt{2\pi\sigma^2}} \exp\left\{-\dfrac{(\log y - \mu)^2}{2\sigma^2}\right\}}
#' for \eqn{y > 0}.
#'
#' \strong{Moments:}
#' \itemize{
#'   \item Expected value: \deqn{\mathbb{E}(y) = \exp\left(\mu + \dfrac{\sigma^2}{2}\right)}
#'   \item Variance: \deqn{\mathbb{V}(y) = \left(\exp(\sigma^2) - 1\right)\exp(2\mu + \sigma^2)}
#'   \item Skewness: \deqn{\gamma_1 = \left(\exp(\sigma^2) + 2\right)\sqrt{\exp(\sigma^2) - 1}}
#'   \item Excess Kurtosis: \deqn{\gamma_2 = \exp(4\sigma^2) + 2\exp(3\sigma^2) + 3\exp(2\sigma^2) - 6}
#' }
#'
#' \strong{Parameter Domains:}
#' \itemize{
#'   \item \eqn{\mu \in (-\infty, +\infty)}
#'   \item \eqn{\sigma^2 \in (0, +\infty)}
#' }
#'
#' \strong{Explicit Formulas:}
#'
#' \emph{Gradient:}
#'
#' The derivatives of the log-pdf \eqn{\ell} are:
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{\log(y) - \mu}{\sigma^2}}
#' \deqn{\dfrac{\partial \ell}{\partial \sigma^2} = \dfrac{(\log(y) - \mu)^2 - \sigma^2}{2(\sigma^2)^2}}
#'
#' \emph{Expected Hessian:}
#'
#' When \code{expected = TRUE}, the expected second derivatives of the log-pdf are:
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\sigma^2}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial (\sigma^2)^2}\right] = -\dfrac{1}{2(\sigma^2)^2}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma^2}\right] = 0}
#'
#' \emph{Observed Hessian:}
#'
#' When \code{expected = FALSE}, the full analytical second derivatives depend on observed \eqn{y}.
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{1}{\sigma^2}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial (\sigma^2)^2} = -\dfrac{(\log(y) - \mu)^2}{(\sigma^2)^3} + \dfrac{1}{2(\sigma^2)^2}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma^2} = -\dfrac{\log(y) - \mu}{(\sigma^2)^2}}
#'
#' @return A list of class \code{"distrib"} containing the components for the Lognormal distribution.
#'
#' @importFrom linkfunctions identity_link log_link
#' @importFrom stats dlnorm plnorm qlnorm rlnorm
#' @export
lognormal_distrib <- function(link_mu = identity_link(), link_sigma2 = log_link()) {
  o <- list()
  class(o) <- c("distrib")

  o$distrib_name <- "lognormal"
  o$type <- "continuous"
  o$dimension <- 1
  o$bounds <- c(0, Inf)

  o$params <- c("mu", "sigma2")
  o$params_interpretation <- c(mu = "mean (log scale)", sigma2 = "variance (log scale)")
  o$n_params <- 2
  o$params_bounds <- list(
    mu = c(-Inf, Inf),
    sigma2 = c(0, Inf)
  )
  o$link_params <- list(
    mu = link_mu,
    sigma2 = link_sigma2
  )

  o$pdf <- function(y, theta, log = FALSE) {
    stats::dlnorm(
      x = y,
      meanlog = theta[[1]],
      sdlog = sqrt(theta[[2]]),
      log = log
    )
  }

  o$cdf <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {
    stats::plnorm(
      q = q,
      meanlog = theta[[1]],
      sdlog = sqrt(theta[[2]]),
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$quantile <- function(p, theta, lower.tail = TRUE, log.p = FALSE) {
    stats::qlnorm(
      p = p,
      meanlog = theta[[1]],
      sdlog = sqrt(theta[[2]]),
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$rng <- function(n, theta) {
    stats::rlnorm(
      n = n,
      meanlog = theta[[1]],
      sdlog = sqrt(theta[[2]])
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
      stop(sprintf(
        "Invalid parameter(s) specified: %s. Available parameters are: %s.",
        paste(sQuote(invalid_pars), collapse = ", "),
        paste(sQuote(o$params), collapse = ", ")
      ))
    }
    mu <- theta[[1]]
    sigma2 <- theta[[2]]
    log_y <- log(y)
    g <- list()

    if ("mu" %in% par) {
      g$mu <- (log_y - mu) / sigma2
    }

    if ("sigma2" %in% par) {
      g$sigma2 <- ((log_y - mu)^2 - sigma2) / (2 * sigma2^2)
    }
    g
  }

  o$hessian <- function(y, theta, expected = FALSE) {
    mu <- theta[[1]]
    sigma2 <- theta[[2]]
    sigma22 <- sigma2 * sigma2

    if (expected) {
      list(
        mu_mu = -1 / sigma2,
        sigma2_sigma2 = -.5 / sigma22,
        mu_sigma2 = 0
      )
    } else {
      log_y <- log(y)
      res <- (log_y - mu)
      list(
        mu_mu = -1 / sigma2,
        sigma2_sigma2 = -res^2 / (sigma22 * sigma2) + .5 / sigma22,
        mu_sigma2 = -res / sigma22
      )
    }
  }

  o$kernel <- function(y, theta, log = TRUE) {
    mu <- theta[[1]]
    sigma2 <- theta[[2]]

    k <- -.5 * (log(y) - mu)^2 / sigma2 - log(y)

    if (log) {
      k
    } else {
      exp(k)
    }
  }

  o$normalization_constant <- function(theta, log = TRUE) {
    z <- sqrt(2 * pi * theta[[2]])

    if (log) {
      log(z)
    } else {
      z
    }
  }

  o$mean <- function(theta) {
    exp(theta[[1]] + theta[[2]] / 2)
  }

  o$variance <- function(theta) {
    (exp(theta[[2]]) - 1) * exp(2 * theta[[1]] + theta[[2]])
  }

  o$skewness <- function(theta) {
    sigma2 <- theta[[2]]
    (exp(sigma2) + 2) * sqrt(exp(sigma2) - 1)
  }

  o$kurtosis <- function(theta) {
    sigma2 <- theta[[2]]
    exp(4 * sigma2) + 2 * exp(3 * sigma2) + 3 * exp(2 * sigma2) - 6
  }

  o$mode <- function(theta) {
    exp(theta[[1]] - theta[[2]])
  }

  o$median <- function(theta) {
    exp(theta[[1]])
  }

  o$initialize <- function(y) {
    log_y <- log(y)
    list(
      mu = mean(log_y),
      sigma2 = var(log_y)
    )
  }

  o
}
