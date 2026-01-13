#' Pseudo-Huber `distrib` Object (Location-Scale Parameterization)
#'
#' @description
#' Creates a `distrib` object for the Pseudo-Huber distribution, which corresponds to the probability density
#' defined by the Pseudo-Huber loss function kernel. It is a specific case of the Generalized Hyperbolic distribution.
#'
#' @param link_mu A link function object for the location parameter \eqn{\mu}.
#'   Defaults to \code{\link[linkfunctions]{identity_link}}.
#' @param link_sigma A link function object for the scale parameter \eqn{\sigma}.
#'   Defaults to \code{\link[linkfunctions]{log_link}} to ensure positivity.
#' @param link_nu A link function object for the shape parameter \eqn{\nu}.
#'   Defaults to \code{\link[linkfunctions]{log_link}} to ensure positivity.
#'
#' @details
#' The probability density function is proportional to the exponential of the negative Pseudo-Huber loss:
#' \deqn{f(y; \mu, \sigma, \nu) = \dfrac{1}{2 \sigma \sqrt{\nu} K_1(\sqrt{\nu})} \exp\left( - \sqrt{\nu + \left(\dfrac{y-\mu}{\sigma}\right)^2} \right)}
#'
#' Where the normalization constant involves the modified Bessel function of the second kind \eqn{K_1}.
#'
#' \strong{Moments:}
#' \itemize{
#'   \item Expected value: \deqn{\mathbb{E}(y) = \mu}
#'   \item Variance: \deqn{\mathbb{V}(y) = \sigma^2 \sqrt{\nu} \dfrac{K_2(\sqrt{\nu})}{K_1(\sqrt{\nu})}}
#'   \item Kurtosis: Depends on ratios of \eqn{K_3, K_1} and \eqn{K_2}.
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
#' Let \eqn{r = y - \mu} and \eqn{D = \sqrt{\nu + (r/\sigma)^2}}.
#'
#' \emph{Gradient (First Derivatives):}
#'
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{r}{\sigma^2 D}}
#'
#' \deqn{\dfrac{\partial \ell}{\partial \sigma} = \dfrac{1}{\sigma} \left( \dfrac{r^2}{\sigma^2 D} - 1 \right)}
#'
#' \deqn{\dfrac{\partial \ell}{\partial \nu} = -\dfrac{1}{2} \left[ \dfrac{1}{\nu} + \dfrac{1}{D} + \dfrac{K_1'(\sqrt{\nu})}{\sqrt{\nu} K_1(\sqrt{\nu})} \right]}
#'
#' Where \eqn{K_1'(z)} is the derivative of the modified Bessel function of the second kind with respect to its argument.
#'
#' \emph{Observed Hessian (Second Derivatives):}
#'
#' The following formulas utilize the observed data (via \eqn{r} and \eqn{D}). Let \eqn{R_1 = \dfrac{K_1'(\sqrt{\nu})}{K_1(\sqrt{\nu})}} and \eqn{R_2 = \dfrac{K_1''(\sqrt{\nu})}{K_1(\sqrt{\nu})}}.
#'
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{\nu}{\sigma^2 D^3}}
#'
#' \deqn{\dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{\sigma^4 - 3\sigma^2 r^2 D^{-1} + r^4 D^{-3}}{\sigma^6}}
#'
#' \deqn{\dfrac{\partial^2 \ell}{\partial \nu^2} = \dfrac{1}{4D^3} + \dfrac{1}{2\nu^2} + \dfrac{1}{4}\left(\dfrac{R_1}{\nu^{3/2}} + \dfrac{R_1^2}{\nu} - \dfrac{R_2}{\nu}\right)}
#'
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma} = \dfrac{-2\nu\sigma^2 r - r^3}{\sigma^2(\nu\sigma^2 + r^2)^{3/2}}}
#'
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \partial \nu} = -\dfrac{r}{2\sigma^2 D^3}}
#'
#' \deqn{\dfrac{\partial^2 \ell}{\partial \sigma \partial \nu} = -\dfrac{r^2}{2\sigma^3 D^3}}
#'
#'
#' \emph{Expected Hessian:}
#'
#' There are no closed-form analytical formulas for the expected Hessian of the Pseudo-Huber distribution.
#' When `expected = TRUE` inside `pseudohuber_distrib()$hessian()`, the expected Hessian is approximated using Monte Carlo simulation
#' via the function \code{\link{mc_expected_hessian}}.
#' Note that the location parameter \eqn{\mu} is orthogonal to the scale parameter \eqn{\sigma} and the shape parameter \eqn{\nu},
#' since \eqn{\mathbb{E}(r) = \mathbb{E}(r^3) = 0}.
#' @return A list of class \code{"distrib"} containing the components for the Pseudo-Huber distribution.
#'
#' @export
pseudohuber_distrib <- function(link_mu = identity_link(), link_sigma = log_link(), link_nu = log_link()) {
  o <- list()
  class(o) <- c("distrib")

  o$distrib_name <- "pseudo huber"
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


  o$kernel <- function(y, theta, log = TRUE) {
    k <- -sqrt(theta[[3]] + ((y - theta[[1]]) / theta[[2]])^2)
    if (log) {
      k
    } else {
      exp(k)
    }
  }

  o$normalization_constant <- function(theta, log = TRUE) {
    nu <- theta[[3]]
    z <- log(2) + log(theta[[2]]) + .5 * log(nu) + log(besselK(sqrt(nu), 1))
    if (log) {
      z
    } else {
      exp(z)
    }
  }


  o$pdf <- function(y, theta, log = FALSE) {
    r <- o$kernel(y, theta, log = TRUE) - o$normalization_constant(theta, log = TRUE)
    if (log) {
      r
    } else {
      exp(r)
    }
  }

  o$cdf <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {
    distrib:::cdf.distrib(o, q, theta, lower.tail, log.p)
  }

  o$quantile <- function(p, theta, lower.tail = TRUE, log.p = FALSE) {
    distrib:::quantile.distrib(o, p, theta, lower.tail, log.p)
  }

  o$rng <- function(n, theta) {
    distrib:::rng.distrib(o, n, theta)
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
    sigma <- theta[[2]]
    nu <- theta[[3]]
    res <- y - mu
    res2 <- res * res
    sigma2 <- sigma * sigma
    den <- sqrt(nu + res2 / sigma2)
    sq_nu <- sqrt(nu)
    g <- list()
    if ("mu" %in% par) {
      g$mu <- res / (sigma2 * den)
    }

    if ("sigma" %in% par) {
      if ("mu" %in% par) {
        g$sigma <- (res * g$mu - 1) / sigma
      } else {
        g$sigma <- ((res * res) / (sigma2 * den) - 1) / sigma
      }
    }

    if ("nu" %in% par) {
      g$nu <- -.5 * (1 / nu + 1 / den + (dbesselK(sq_nu, 1) / besselK(sq_nu, 1)) / sq_nu)
    }
    g
  }

  o$hessian <- function(y, theta, expected = FALSE) {
    if (expected) {
      h <- mc_expected_hessian(o, y, theta)
      h$mu_sigma <- h$mu_nu <- 0
      h
    } else {
      mu <- theta[[1]]
      sigma <- theta[[2]]
      nu <- theta[[3]]
      res <- y - mu
      res2 <- res * res
      sigma2 <- sigma * sigma
      sigma4 <- sigma2 * sigma2
      den <- sqrt(nu + res2 / sigma2)
      sq_nu <- sqrt(nu)
      bk <- besselK(sq_nu, 1)
      r1 <- dbesselK(sq_nu, 1, deriv = 1, mode = "standard") / bk
      r2 <- dbesselK(sq_nu, 1, deriv = 2, mode = "standard") / bk
      list(
        mu_mu = -nu / (sigma2 * den^3),
        sigma_sigma = (sigma4 - 3 * sigma2 * res2 / den + res2 * res2 / den^3) / (sigma4 * sigma2),
        nu_nu = 0.25 / (den^3) + 0.5 / (nu^2) + 0.25 * (r1 * nu^(-1.5) + r1^2 / nu - r2 / nu),
        mu_sigma = (-2 * nu * sigma2 * res - res^3) / (sigma2 * (nu * sigma2 + res^2)^1.5),
        mu_nu = (-res / (2 * sigma2)) / (((res2 + nu * sigma2) / sigma2)^1.5),
        sigma_nu = (-res^2) / (2 * sigma2 * sigma * den^3)
      )
    }
  }


  o$mean <- o$median <- o$mode <- function(theta) {
    theta[[1]]
  }

  o$variance <- function(theta) {
    sigma <- theta[[2]]
    sqrt_nu <- sqrt(theta[[3]])
    ratio_bessel <- besselK(sqrt_nu, 2) / besselK(sqrt_nu, 1)
    sigma^2 * sqrt_nu * ratio_bessel
  }

  o$skewness <- function(theta) {
    0
  }

  o$kurtosis <- function(theta) {
    sqrt_nu <- sqrt(theta[[3]])
    3 * (besselK(sqrt_nu, 3) * besselK(sqrt_nu, 1)) / (besselK(sqrt_nu, 2)^2) - 3
  }

  o
}
