#' Student's t `distrib` Object (Location-Scale Parameterization)
#'
#' @description
#' Creates a distribution object for the Student's t distribution parameterized by location (\eqn{\mu}), scale (\eqn{\sigma}), and degrees of freedom (\eqn{\nu}).
#'
#' @param link_mu A link function object for the location parameter \eqn{\mu}.
#'   Defaults to \code{\link[linkfunctions]{identity_link}}.
#' @param link_sigma A link function object for the scale parameter \eqn{\sigma}.
#'   Defaults to \code{\link[linkfunctions]{log_link}} to ensure positivity.
#' @param link_nu A link function object for the degrees of freedom parameter \eqn{\nu}.
#'   Defaults to \code{\link[linkfunctions]{log_link}} to ensure positivity.
#'
#' @details
#' The probability density function is given by:
#' \deqn{f(y; \mu, \sigma, \nu) = \dfrac{\Gamma\left(\dfrac{\nu+1}{2}\right)}{\sigma\sqrt{\nu\pi}\,\Gamma\left(\dfrac{\nu}{2}\right)} \left(1 + \dfrac{(y-\mu)^2}{\nu\sigma^2}\right)^{-\dfrac{\nu+1}{2}}}
#' for \eqn{y \in (-\infty, +\infty)}.
#'
#' \strong{Moments:}
#' \itemize{
#'   \item Expected value: \deqn{\mathbb{E}(y) = \mu \quad (\text{if } \nu > 1)}
#'   \item Variance: \deqn{\mathbb{V}(y) = \sigma^2 \dfrac{\nu}{\nu-2} \quad (\text{if } \nu > 2)}
#'   \item Skewness: \deqn{\gamma_1 = 0 \quad (\text{if } \nu > 3)}
#'   \item Excess Kurtosis: \deqn{\gamma_2 = \dfrac{6}{\nu-4} \quad (\text{if } \nu > 4)}
#' }
#'
#' \strong{Parameter Domains:}
#' \itemize{
#'   \item \eqn{\mu \in (-\infty, +\infty)}
#'   \item \eqn{\sigma \in (0, +\infty)}
#'   \item \eqn{\nu \in (0, +\infty)}
#' }
#'
#' \strong{Explicit Derivatives:}
#'
#' Let \eqn{\psi(z)} be the digamma function and \eqn{\psi_1(z)} be the trigamma function.
#'
#' \emph{Gradient (First Derivatives):}
#'
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{(\nu+1)(y-\mu)}{\nu\sigma^2 + (y-\mu)^2}}
#'
#' \deqn{\dfrac{\partial \ell}{\partial \sigma} = \dfrac{\nu\left[(y-\mu)^2 - \sigma^2\right]}{\sigma\left[\nu\sigma^2 + (y-\mu)^2\right]}}
#'
#' \deqn{\dfrac{\partial \ell}{\partial \nu} = \dfrac{1}{2}\left[ -\dfrac{1}{\nu} - \psi\left(\dfrac{\nu}{2}\right) + \psi\left(\dfrac{\nu+1}{2}\right) + \dfrac{(\nu+1)(y-\mu)^2}{\nu\left[\nu\sigma^2 + (y-\mu)^2\right]} - \log\left(1 + \dfrac{(y-\mu)^2}{\nu\sigma^2}\right) \right]}
#'
#' \emph{Expected Hessian:}
#'
#' When \code{expected = TRUE}, the formulas depend only on the parameters (not on \eqn{y}).
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{\nu+1}{\sigma^2(\nu+3)}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right] = -\dfrac{2\nu}{\sigma^2(\nu+3)}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \nu^2}\right] = \dfrac{1}{4}\left[\psi_1\left(\dfrac{\nu+1}{2}\right) - \psi_1\left(\dfrac{\nu}{2}\right)\right] + \dfrac{\nu+5}{2\nu(\nu+1)(\nu+3)}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma \partial \nu}\right] = \dfrac{2}{\sigma(\nu+1)(\nu+3)}}
#'
#' The parameter \eqn{\mu} is orthogonal to \eqn{\sigma} and \eqn{\nu} (mixed expected derivatives are 0).
#'
#' \emph{Observed Hessian (Second Derivatives):}
#'
#' When \code{expected = FALSE}, the analytical second derivatives depend on observed \eqn{y}.
#'
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = \dfrac{(\nu+1)\left[(y-\mu)^2 - \nu\sigma^2\right]}{\left[\nu\sigma^2 + (y-\mu)^2\right]^2}}
#'
#' \deqn{\dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{\nu\left[\nu\sigma^4 - (3\nu+1)\sigma^2(y-\mu)^2 - (y-\mu)^4\right]}{\sigma^2\left[\nu\sigma^2 + (y-\mu)^2\right]^2}}
#'
#' \deqn{\dfrac{\partial^2 \ell}{\partial \nu^2} = \dfrac{1}{4}\left[ -\psi_1\left(\dfrac{\nu}{2}\right) + \psi_1\left(\dfrac{\nu+1}{2}\right) + \dfrac{2\left(\nu\sigma^4 + (y-\mu)^4\right)}{\nu\left[\nu\sigma^2 + (y-\mu)^2\right]^2} \right]}
#'
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma} = -\dfrac{2\nu(\nu+1)\sigma(y-\mu)}{\left[\nu\sigma^2 + (y-\mu)^2\right]^2}}
#'
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \partial \nu} = \dfrac{(y-\mu)\left[(y-\mu)^2 - \sigma^2\right]}{\left[\nu\sigma^2 + (y-\mu)^2\right]^2}}
#'
#' \deqn{\dfrac{\partial^2 \ell}{\partial \sigma \partial \nu} = \dfrac{(y-\mu)^2\left[(y-\mu)^2 - \sigma^2\right]}{\sigma\left[\nu\sigma^2 + (y-\mu)^2\right]^2}}
#'
#' @return A list of class \code{"distrib"} containing the components for the Student's t distribution.
#'
#' @importFrom linkfunctions identity_link log_link
#' @importFrom stats dt pt qt rt
#' @export
student_t_distrib <- function(
    link_mu = identity_link(),
    link_sigma = log_link(),
    link_nu = log_link()) {
  o <- list()
  class(o) <- c("distrib")

  o$distrib_name <- "student t"
  o$type <- "continuous"
  o$dimension <- 1
  o$bounds <- c(-Inf, Inf)

  o$params <- c("mu", "sigma", "nu")
  o$params_interpretation <- c(mu = "location", sigma = "scale", nu = "shape")
  o$n_params <- 3
  o$params_bounds <- list(
    mu = c(-Inf, Inf),
    sigma = c(0, Inf),
    nu = c(0, Inf)
  )
  o$link_params <- list(
    mu = link_mu,
    sigma = link_sigma,
    nu = link_nu
  )

  o$pdf <- function(y, theta, log = FALSE) {
    mu <- theta[["mu"]]
    sigma <- theta[["sigma"]]
    nu <- theta[["nu"]]
    val <- stats::dt(
      x = (y - mu) / sigma,
      df = nu,
      log = log
    )

    if (log) {
      val - log(sigma)
    } else {
      val / sigma
    }
  }

  o$cdf <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {
    mu <- theta[["mu"]]
    sigma <- theta[["sigma"]]
    nu <- theta[["nu"]]
    stats::pt(
      q = (q - mu) / sigma,
      df = nu,
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$qf <- function(p, theta, lower.tail = TRUE, log.p = FALSE) {
    mu <- theta[["mu"]]
    sigma <- theta[["sigma"]]
    nu <- theta[["nu"]]
    mu + sigma * stats::qt(
      p = p,
      df = nu,
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$rng <- function(n, theta) {
    mu <- theta[["mu"]]
    sigma <- theta[["sigma"]]
    nu <- theta[["nu"]]
    mu + sigma * stats::rt(
      n = n,
      df = nu
    )
  }

  o$loglik <- function(y, theta) {
    o$pdf(y, theta, log = TRUE)
  }

  o$gradient <- function(y, theta) {
    mu <- theta[["mu"]]
    sigma <- theta[["sigma"]]
    nu <- theta[["nu"]]
    res <- y - mu
    res2 <- res * res
    sigma2 <- sigma * sigma
    list(
      mu = ((nu + 1) * res) / (nu * sigma2 + res2),
      sigma = (nu * (res2 - sigma2)) / (sigma * (nu * sigma2 + res2)),
      nu = 0.5 * (-1 / nu - digamma(nu / 2) + digamma((nu + 1) / 2) + ((nu + 1) * (y - mu)^2) / (nu * ((y - mu)^2 + nu * sigma^2)) - log(((y - mu)^2) / (nu * sigma^2) + 1))
    )
  }

  o$hessian <- function(y, theta, expected = FALSE) {
    mu <- theta[["mu"]]
    sigma <- theta[["sigma"]]
    nu <- theta[["nu"]]
    sigma2 <- sigma * sigma

    if (expected) {
      list(
        mu_mu = -(nu + 1) / (sigma2 * (nu + 3)),
        sigma_sigma = -(2 * nu) / (sigma2 * (nu + 3)),
        nu_nu = .25 * (trigamma(.5 * (nu + 1)) - trigamma(.5 * nu)) + (nu + 5) / (2 * nu * (nu + 1) * (nu + 3)),
        mu_sigma = 0,
        mu_nu = 0,
        sigma_nu = 2 / ((nu + 1) * (nu + 3) * sigma)
      )
    } else {
      res <- y - mu
      res2 <- res * res
      res4 <- res2 * res2
      sigma4 <- sigma2 * sigma2
      list(
        mu_mu = (nu + 1) * (res2 - nu * sigma2) / (nu * sigma2 + res2)^2,
        sigma_sigma = (nu * (nu * sigma4 - (3 * nu + 1) * sigma2 * res2 - res4)) / (nu * sigma2 * sigma + sigma * res2)^2,
        nu_nu = .25 * (-trigamma(.5 * nu) + trigamma(.5 * (nu + 1)) + (2 * (nu * sigma4 + res4)) / (nu * (nu * sigma2 + res2)^2)),
        mu_sigma = -(2 * nu * (nu + 1) * sigma * res) / (nu * sigma^2 + res2)^2,
        mu_nu = (-sigma2 * res + res2 * res) / (nu * sigma2 + res2)^2,
        sigma_nu = (res4 - sigma2 * res2) / (sigma * (nu * sigma2 + res2)^2)
      )
    }
  }

  o$kernel <- function(y, theta, log = TRUE) {
    mu <- theta[["mu"]]
    sigma <- theta[["sigma"]]
    nu <- theta[["nu"]]
    z <- (y - mu) / sigma
    k <- -.5 * (nu + 1) * log(1 + z^2 / nu)

    if (log) {
      k
    } else {
      exp(k)
    }
  }

  o$normalization_constant <- function(theta, log = TRUE) {
    nu <- theta[["nu"]]
    z <- lgamma(.5 * nu) + .5 * log(nu * pi) - lgamma(.5 * (nu + 1)) + log(theta[["sigma"]])

    if (log) {
      z
    } else {
      exp(z)
    }
  }

  o$mean <- o$median <- o$mode <- function(theta) {
    ifelse(theta[["nu"]] > 1, theta[["mu"]], NA)
  }

  o$variance <- function(theta) {
    nu <- theta[["nu"]]
    ifelse(nu > 2,
      (theta[["sigma"]]^2) * nu / (nu - 2),
      ifelse(nu > 1 & nu <= 2, Inf, NA)
    )
  }

  o$skewness <- function(theta) {
    ifelse(theta[["nu"]] > 3, 0, NA)
  }

  o$kurtosis <- function(theta) {
    nu <- theta[["nu"]]
    ifelse(nu > 4,
      6 / (nu - 4),
      ifelse(nu > 2 & nu <= 4, Inf, NA)
    )
  }

  o
}
