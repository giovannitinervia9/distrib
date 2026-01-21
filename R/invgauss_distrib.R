#' Inverse-Gaussian `distrib` Object (Mean-Dispersion Parameterization)
#'
#' @description
#' Creates a distribution object for the Inverse-Gaussian distribution parameterized by mean (\eqn{\mu}) and dispersion (\eqn{\phi}).
#'
#' @param link_mu A link function object for the mean parameter \eqn{\mu}.
#'   Defaults to \code{\link[linkfunctions]{log_link}} to ensure positivity.
#' @param link_phi A link function object for the dispersion parameter \eqn{\phi}.
#'   Defaults to \code{\link[linkfunctions]{log_link}} to ensure positivity.
#'
#' @details
#' The Inverse-Gaussian distribution has the following density function:
#' \deqn{f(y; \mu, \phi) = \sqrt{\dfrac{1}{2\pi\phi y^3}} \exp\left\{-\dfrac{(y-\mu)^2}{2\phi\mu^2 y}\right\}}
#' for \eqn{y > 0}.
#'
#' The probability density, cumulative distribution, quantile, and random generation functions
#' are provided by the \code{\link[statmod]{dinvgauss}}, \code{\link[statmod]{pinvgauss}},
#' \code{\link[statmod]{qinvgauss}}, and \code{\link[statmod]{rinvgauss}} functions
#' from the \pkg{statmod} package.
#'
#' \strong{Moments:}
#' \itemize{
#'   \item Expected value: \deqn{\mathbb{E}(y) = \mu}
#'   \item Variance: \deqn{\mathbb{V}(y) = \phi \mu^3}
#'   \item Skewness: \deqn{\gamma_1 = 3\sqrt{\phi \mu}}
#'   \item Excess Kurtosis: \deqn{\gamma_2 = 15 \phi \mu}
#' }
#'
#' \strong{Parameter Domains:}
#' \itemize{
#'   \item \eqn{\mu \in (0, +\infty)}
#'   \item \eqn{\phi \in (0, +\infty)}
#' }
#'
#' \strong{Explicit Formulas:}
#'
#' \emph{Gradient:}
#'
#' The derivatives of the log-pdf \eqn{\ell} with respect to the parameters are:
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\phi\mu^3}}
#' \deqn{\dfrac{\partial \ell}{\partial \phi} = \dfrac{(y - \mu)^2 - y\mu^2\phi}{2y\phi^2\mu^2}}
#'
#' \emph{Expected Hessian:}
#'
#' When \code{expected = TRUE}, the expected second derivatives of the log-pdf are:
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\phi\mu^3}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \phi^2}\right] = -\dfrac{1}{2\phi^2}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu \partial \phi}\right] = 0}
#'
#' \emph{Observed Hessian:}
#'
#' When \code{expected = FALSE}, the second derivatives of the log-pdf depend on observed \eqn{y}:
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{3y - 2\mu}{\phi\mu^4}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \phi^2} = \dfrac{\phi - 2(y-\mu)^2/(\mu^2 y)}{2\phi^3}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \partial \phi} = -\dfrac{y - \mu}{\phi^2\mu^3}}
#' \strong{Note:} The observed Hessian with respect to \eqn{\phi} is not guaranteed to be negative for all
#' observed values of \eqn{y}. Specifically, \eqn{\partial^2 \ell/\partial \phi^2 < 0} only when
#' \eqn{\phi < 2(y-\mu)^2/(\mu^2 y)}. This condition may be violated when observations are far from the mean
#' or when the dispersion parameter is large, potentially causing numerical instability in optimization
#' algorithms that rely on the observed Hessian (e.g., Newton-Raphson). In such cases, using the expected
#' Hessian (\code{expected = TRUE}) is recommended for more stable convergence.
#' @return A list of class \code{"distrib"} containing the components for the Inverse-Gaussian distribution.
#'
#' @references
#' Giner, G., and Smyth, G. K. (2016). statmod: Probability calculations for the inverse Gaussian distribution.
#' \emph{R Journal} \strong{8}(1), 339-351. \url{https://journal.r-project.org/articles/RJ-2016-024/}
#'
#' @importFrom statmod dinvgauss pinvgauss qinvgauss rinvgauss
#' @export
invgauss_distrib <- function(link_mu = log_link(), link_phi = log_link()) {
  o <- list()
  class(o) <- c("distrib")

  o$distrib_name <- "inverse gaussian"
  o$type <- "continuous"
  o$dimension <- 1
  o$bounds <- c(0, Inf)

  o$params <- c("mu", "phi")
  o$params_interpretation <- c(mu = "mean", phi = "dispersion")
  o$n_params <- 2
  o$params_bounds <- list(
    mu = c(0, Inf),
    phi = c(0, Inf)
  )
  o$link_params <- list(
    mu = link_mu,
    phi = link_phi
  )

  o$pdf <- function(y, theta, log = FALSE) {
    statmod::dinvgauss(
      x = y,
      mean = theta[[1]],
      dispersion = theta[[2]],
      log = log
    )
  }

  o$cdf <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {
    statmod::pinvgauss(
      q = q,
      mean = theta[[1]],
      dispersion = theta[[2]],
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$quantile <- function(p, theta, lower.tail = TRUE, log.p = FALSE) {
    statmod::qinvgauss(
      p = p,
      mean = theta[[1]],
      dispersion = theta[[2]],
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$rng <- function(n, theta) {
    statmod::rinvgauss(
      n = n,
      mean = theta[[1]],
      dispersion = theta[[2]]
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
    mu <- theta[[1]]
    phi <- theta[[2]]
    res <- y - mu
    mu2 <- mu * mu
    g <- list()

    if ("mu" %in% par) {
      g$mu <- res / (phi * mu2 * mu)
    }

    if ("phi" %in% par) {
      g$phi <- (res^2 - y * mu2 * phi) / (2 * y * phi * phi * mu2)
    }
    g
  }

  o$hessian <- function(y, theta, expected = FALSE) {
    mu <- theta[[1]]
    phi <- theta[[2]]
    mu2 <- mu * mu
    phi2 <- phi * phi

    if (expected) {
      list(
        mu_mu = -1 / (phi * mu2 * mu),
        phi_phi = -.5 / phi2,
        mu_phi = 0
      )
    } else {
      res <- y - mu
      list(
        mu_mu = -(3 * y - 2 * mu) / (phi * mu2 * mu2),
        phi_phi = (phi - 2 * res^2 / (mu2 * y)) / (2 * phi2 * phi),
        mu_phi = -res / (phi2 * mu2 * mu)
      )
    }
  }

  o$kernel <- function(y, theta, log = TRUE) {
    mu <- theta[[1]]
    phi <- theta[[2]]
    k <- -.5 * (y - mu)^2 / (phi * mu^2 * y) - 1.5 * log(y)

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
    theta[[1]]
  }

  o$variance <- function(theta) {
    theta[[1]]^3 * theta[[2]]
  }

  o$skewness <- function(theta) {
    3 * sqrt(theta[[1]] * theta[[2]])
  }

  o$kurtosis <- function(theta) {
    15 * theta[[1]] * theta[[2]]
  }

  o$mode <- function(theta) {
    mu <- theta[[1]]
    phi <- theta[[2]]
    k <- 1.5 * mu * phi
    mu * (sqrt(1 + k^2) - k)
  }

  o$median <- function(theta) {
    o$quantile(.5, theta)
  }

  o$initialize <- function(y) {
    mu <- mean(y)
    list(
      mu = mu,
      phi = var(y) / mu^3
    )
  }

  o
}
