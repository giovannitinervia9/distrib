#' Logistic `distrib` Object
#'
#' @description
#' Creates a distribution object for the Logistic distribution parameterized by location (\eqn{\mu})
#' and scale (\eqn{\sigma}).
#'
#' @param link_mu A link function object for the location parameter \eqn{\mu}.
#'   Defaults to \code{\link[linkfunctions]{identity_link}}.
#' @param link_sigma A link function object for the scale parameter \eqn{\sigma}.
#'   Defaults to \code{\link[linkfunctions]{log_link}} to ensure positivity.
#'
#' @details
#' The Logistic distribution has the following probability density function:
#' \deqn{f(y; \mu, \sigma) = \dfrac{\exp\left(-\frac{y-\mu}{\sigma}\right)}{\sigma \left[1 + \exp\left(-\frac{y-\mu}{\sigma}\right)\right]^2}}
#'
#' \strong{Moments:}
#' \itemize{
#'   \item Expected value: \eqn{\mathbb{E}(y) = \mu}
#'   \item Variance: \eqn{\mathbb{V}(y) = \dfrac{\sigma^2 \pi^2}{3}}
#'   \item Skewness: \eqn{\gamma_1 = 0}
#'   \item Excess Kurtosis: \eqn{\gamma_2 = 1.2}
#' }
#'
#' \strong{Explicit Formulas:}
#'
#' \emph{Gradient:}
#' The derivatives of the log-likelihood \eqn{\ell} with respect to the parameters are:
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{1}{\sigma} \tanh\left(\dfrac{y-\mu}{2\sigma}\right)}
#' \deqn{\dfrac{\partial \ell}{\partial \sigma} = -\dfrac{1}{\sigma} \left[ 1 - \dfrac{y-\mu}{\sigma} \tanh\left(\dfrac{y-\mu}{2\sigma}\right) \right]}
#'
#' \emph{Observed Hessian:}
#' The second order derivatives of the log-likelihood are:
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{1}{2\sigma^2} \text{sech}^2\left(\dfrac{y-\mu}{2\sigma}\right)}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{1}{\sigma^2} \left[ 1 - \dfrac{(y-\mu)^2}{2\sigma^2} \text{sech}^2\left(\dfrac{y-\mu}{2\sigma}\right) - \dfrac{y-\mu}{\sigma} \tanh\left(\dfrac{y-\mu}{2\sigma}\right) \right]}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma} = -\dfrac{1}{\sigma^2} \left[ \tanh\left(\dfrac{y-\mu}{2\sigma}\right) + \dfrac{y-\mu}{2\sigma} \text{sech}^2\left(\dfrac{y-\mu}{2\sigma}\right) \right]}
#'
#' \emph{Expected Hessian:}
#' The Fisher Information matrix has a closed constant form:
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{3\sigma^2}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right] = -\dfrac{3+\pi^2}{9\sigma^2}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma}\right] = 0}
#'
#' @return A list of class \code{"distrib"} containing the components for the Logistic distribution.
#'
#' @importFrom linkfunctions identity_link log_link
#' @importFrom stats dlogis plogis qlogis rlogis var sd
#' @export
logistic_distrib <- function(link_mu = identity_link(), link_sigma = log_link()) {
  o <- list()
  class(o) <- c("distrib")

  o$distrib_name <- "logistic"
  o$type <- "continuous"
  o$dimension <- 1
  o$bounds <- c(-Inf, Inf)

  o$params <- c("mu", "sigma")
  o$params_interpretation <- c(mu = "mean", sigma = "scale")
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
    stats::dlogis(
      x = y,
      location = theta[[1]],
      scale = theta[[2]],
      log = log
    )
  }

  o$cdf <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {
    stats::plogis(
      q = q,
      location = theta[[1]],
      scale = theta[[2]],
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$quantile <- function(p, theta, lower.tail = TRUE, log.p = FALSE) {
    stats::qlogis(
      p = p,
      location = theta[[1]],
      scale = theta[[2]],
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$rng <- function(n, theta) {
    stats::rlogis(
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

    mu <- theta[[1]]
    sigma <- theta[[2]]
    res <- y - mu
    tanh_z <- tanh(.5 * res / sigma)

    g <- list()
    if ("mu" %in% par) {
      g$mu <- tanh_z / sigma
    }
    if ("sigma" %in% par) {
      g$sigma <- (res * tanh_z - sigma) / sigma^2
    }
    g
  }

  o$hessian <- function(y, theta, expected = FALSE) {
    mu <- theta[[1]]
    sigma <- theta[[2]]
    sigma2 <- sigma * sigma

    if (expected) {
      list(
        mu_mu = -1 / (3 * sigma2),
        sigma_sigma = -(3 + pi^2) / (9 * sigma2),
        mu_sigma = 0
      )
    } else {
      z <- (y - mu) / sigma
      z_half <- .5 * z
      tanh_z <- tanh(z_half)
      sech2_z <- 1 - tanh_z^2

      list(
        mu_mu = -sech2_z / (2 * sigma2),
        sigma_sigma = (1 - 4 * z_half * tanh_z - 2 * z_half^2 * sech2_z) / sigma2,
        mu_sigma = -(tanh_z + z_half * sech2_z) / sigma2
      )
    }
  }

  o$kernel <- function(y, theta, log = TRUE) {
    sigma <- theta[[2]]
    z <- (y - theta[[1]]) / sigma
    k <- -z - 2 * log(1 + exp(-z))
    if (log) {
      k
    } else {
      exp(k)
    }
  }

  o$normalization_constant <- function(theta, log = TRUE) {
    if (log) {
      log(theta[[2]])
    } else {
      theta[[2]]
    }
  }

  o$mean <- o$median <- o$mode <- function(theta) {
    theta[[1]]
  }


  o$variance <- function(theta) {
    ((theta[[2]] * pi)^2) / 3
  }

  o$skewness <- function(theta) {
    0
  }

  o$kurtosis <- function(theta) {
    1.2
  }

  o$initialize <- function(y) {
    list(
      mu = mean(y),
      sigma = sd(y) * sqrt(3) / pi
    )
  }

  o
}
