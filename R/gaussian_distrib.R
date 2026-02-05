#' Gaussian `distrib` Object (Standard Deviation Parametrization)
#'
#' @description
#' Creates a distribution object for the Gaussian distribution parameterized by mean (\eqn{\mu}) and standard deviation (\eqn{\sigma}).
#'
#' @param link_mu A link function object for the location parameter \eqn{\mu}.
#'   Defaults to \code{\link[linkfunctions]{identity_link}}.
#' @param link_sigma A link function object for the scale parameter \eqn{\sigma}.
#'   Defaults to \code{\link[linkfunctions]{log_link}} to ensure positivity.
#'
#' @details
#' The Gaussian distribution has the following density function:
#' \deqn{f(y; \mu, \sigma) = \dfrac{1}{\sqrt{2\pi}\sigma} \exp\left\{-\dfrac{1}{2}\left(\dfrac{y-\mu}{\sigma}\right)^2\right\}}
#'
#' \strong{Moments:}
#' \itemize{
#'   \item Expected value: \eqn{\mathbb{E}(y) = \mu}
#'   \item Variance: \eqn{\mathbb{V}(y) = \sigma^2}
#'   \item Skewness: \eqn{\gamma_1 = 0}
#'   \item Excess Kurtosis: \eqn{\gamma_2 = 0}
#' }
#'
#' \strong{Parameter Domains:}
#' \itemize{
#'   \item \eqn{\mu \in (-\infty, +\infty)}
#'   \item \eqn{\sigma \in (0, +\infty)}
#' }
#'
#' \strong{Explicit Formulas:}
#'
#' \emph{Cumulative Distribution Function (CDF):}
#' \deqn{F(y; \mu, \sigma) = \Phi\left(\dfrac{y-\mu}{\sigma}\right)}
#' where \eqn{\Phi} is the standard normal CDF.
#'
#' \emph{Quantile Function (QF):}
#' \deqn{Q(p; \mu, \sigma) = \mu + \sigma \Phi^{-1}(p)}
#'
#' \emph{Gradient:}
#' The derivatives of the log-pdf \eqn{\ell} with respect to the parameters are:
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\sigma^2}}
#' \deqn{\dfrac{\partial \ell}{\partial \sigma} = \dfrac{(y - \mu)^2 - \sigma^2}{\sigma^3}}
#'
#' \emph{Observed Hessian:}
#' The second derivatives of the log-pdf are:
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{1}{\sigma^2}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \sigma^2} = \dfrac{\sigma^2 - 3(y - \mu)^2}{\sigma^4}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma} = -\dfrac{2(y - \mu)}{\sigma^3}}
#'
#' \emph{Expected Hessian:}
#' The expected values of the negative Hessian (used when \code{expected = TRUE}) are:
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\sigma^2}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \sigma^2}\right] = -\dfrac{2}{\sigma^2}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma}\right] = 0}
#'
#' \emph{Derivatives with respect to \eqn{y}:}
#' \deqn{\dfrac{\partial \ell}{\partial y} = -\dfrac{y-\mu}{\sigma^2}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial y^2} = -\dfrac{1}{\sigma^2}}
#'
#' @return A list of class \code{"distrib"} containing:
#' \item{distrib_name}{The name of the distribution ("gaussian").}
#' \item{type}{The type of distribution ("continuous").}
#' \item{params}{The names of the parameters ("mu", "sigma").}
#' \item{params_bounds}{The bounds of the parameters.}
#' \item{pdf, cdf, quantile, rng}{Functions for density, cumulative probability, quantiles, and random generation.}
#' \item{gradient}{A function to calculate the gradient of the log-pdf.}
#' \item{hessian}{A function to calculate the Hessian matrix (observed or expected).}
#'
#' @importFrom linkfunctions identity_link log_link
#' @importFrom stats dnorm pnorm qnorm rnorm
#' @export
gaussian_distrib <- function(link_mu = identity_link(), link_sigma = log_link()) {
  o <- list()
  class(o) <- c("distrib")

  o$distrib_name <- "gaussian"
  o$type <- "continuous"
  o$dimension <- 1
  o$bounds <- c(-Inf, Inf)

  o$params <- c("mu", "sigma")
  o$params_interpretation <- c(mu = "mean", sigma = "standard deviation")
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
    stats::dnorm(
      x = y,
      mean = theta[[1]],
      sd = theta[[2]],
      log = log
    )
  }

  o$cdf <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {
    stats::pnorm(
      q = q,
      mean = theta[[1]],
      sd = theta[[2]],
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$quantile <- function(p, theta, lower.tail = TRUE, log.p = FALSE) {
    stats::qnorm(
      p = p,
      mean = theta[[1]],
      sd = theta[[2]],
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$rng <- function(n, theta) {
    stats::rnorm(
      n = n,
      mean = theta[[1]],
      sd = theta[[2]]
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
    sigma2 <- sigma * sigma
    residuals <- y - theta[[1]]
    g <- list()
    if ("mu" %in% par) {
      g$mu <- residuals / sigma2
    }
    if ("sigma" %in% par) {
      g$sigma <- (residuals * residuals - sigma2) / (sigma2 * sigma)
    }
    g
  }

  o$hessian <- function(y, theta, expected = FALSE) {
    sigma <- theta[[2]]
    sigma2 <- sigma * sigma
    if (expected) {
      list(
        mu_mu = -1 / sigma2,
        sigma_sigma = -2 / sigma2,
        mu_sigma = 0
      )
    } else {
      residuals <- y - theta[[1]]
      list(
        mu_mu = -1 / sigma2,
        sigma_sigma = (sigma2 - 3 * residuals * residuals) / (sigma2 * sigma2),
        mu_sigma = -2 * residuals / (sigma2 * sigma)
      )
    }
  }

  o$grad_y <- function(y, theta) {
    -(y - theta[[1]]) / (theta[[2]]^2)
  }

  o$hess_y <- function(y, theta) {
    -1 / (theta[[2]]^2)
  }

  o$kernel <- function(y, theta, log = TRUE) {
    k <- -.5 * ((y - theta[[1]]) / theta[[2]])^2
    if (log) {
      k
    } else {
      exp(k)
    }
  }

  o$normalization_constant <- function(theta, log = TRUE) {
    z <- sqrt(2 * pi) * theta[[2]]

    if (log) {
      log(z)
    } else {
      z
    }
  }

  o$mean <- o$median <- o$mode <- function(theta) {
    theta[[1]]
  }

  o$variance <- function(theta) {
    theta[[2]]^2
  }

  o$skewness <- function(theta) {
    0
  }

  o$kurtosis <- function(theta) {
    0
  }

  o$initialize <- function(y) {
    list(
      mu = mean(y),
      sigma = sd(y)
    )
  }

  o
}





#' Gaussian `distrib` Object (Variance Parameterization)
#'
#' @description
#' Creates a distribution object for the Gaussian distribution parameterized by mean (\eqn{\mu}) and variance (\eqn{\sigma^2}).
#'
#' @param link_mu A link function object for the location parameter \eqn{\mu}.
#'   Defaults to \code{\link[linkfunctions]{identity_link}}.
#' @param link_sigma2 A link function object for the variance parameter \eqn{\sigma^2}.
#'   Defaults to \code{\link[linkfunctions]{log_link}} to ensure positivity.
#'
#' @details
#' The Gaussian distribution with variance parameterization has the following density function:
#' \deqn{f(y; \mu, \sigma^2) = \dfrac{1}{\sqrt{2\pi\sigma^2}} \exp\left\{-\dfrac{(y-\mu)^2}{2\sigma^2}\right\}}
#'
#' \strong{Moments:}
#' \itemize{
#'   \item Expected value: \eqn{\mathbb{E}(y) = \mu}
#'   \item Variance: \eqn{\mathbb{V}(y) = \sigma^2}
#'   \item Skewness: \eqn{\gamma_1 = 0}
#'   \item Excess Kurtosis: \eqn{\gamma_2 = 0}
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
#' The derivatives of the log-pdf \eqn{\ell} with respect to the parameters are:
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \dfrac{y - \mu}{\sigma^2}}
#' \deqn{\dfrac{\partial \ell}{\partial \sigma^2} = \dfrac{(y - \mu)^2 - \sigma^2}{2(\sigma^2)^2}}
#'
#' \emph{Observed Hessian:}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{1}{\sigma^2}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial (\sigma^2)^2} = \dfrac{\sigma^2 - 2(y - \mu)^2}{2(\sigma^2)^3}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma^2} = -\dfrac{y - \mu}{(\sigma^2)^2}}
#'
#' \emph{Expected Hessian:}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\dfrac{1}{\sigma^2}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial (\sigma^2)^2}\right] = -\dfrac{1}{2(\sigma^2)^2}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu \partial \sigma^2}\right] = 0}
#'
#' \emph{Derivatives with respect to \eqn{y}:}
#' \deqn{\dfrac{\partial \ell}{\partial y} = -\dfrac{y-\mu}{\sigma^2}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial y^2} = -\dfrac{1}{\sigma^2}}
#'
#' @return A list of class \code{"distrib"} containing the components for the Gaussian distribution (variance parameterization).
#'
#' @importFrom linkfunctions identity_link log_link
#' @importFrom stats dnorm pnorm qnorm rnorm
#' @export
gaussian2_distrib <- function(link_mu = identity_link(), link_sigma2 = log_link()) {
  o <- list()
  class(o) <- c("distrib")

  o$distrib_name <- "gaussian"
  o$type <- "continuous"
  o$dimension <- 1
  o$bounds <- c(-Inf, Inf)

  o$params <- c("mu", "sigma2")
  o$params_interpretation <- c(mu = "mean", sigma2 = "variance")
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
    stats::dnorm(
      x = y,
      mean = theta[[1]],
      sd = sqrt(theta[[2]]),
      log = log
    )
  }

  o$cdf <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {
    stats::pnorm(
      q = q,
      mean = theta[[1]],
      sd = sqrt(theta[[2]]),
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$quantile <- function(p, theta, lower.tail = TRUE, log.p = FALSE) {
    stats::qnorm(
      p = p,
      mean = theta[[1]],
      sd = sqrt(theta[[2]]),
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$rng <- function(n, theta) {
    stats::rnorm(
      n = n,
      mean = theta[[1]],
      sd = sqrt(theta[[2]])
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
    sigma2 <- theta[[2]]
    residuals <- y - theta[[1]]
    g <- list()
    if ("mu" %in% par) {
      g$mu <- residuals / sigma2
    }
    if ("sigma2" %in% par) {
      g$sigma2 <- (residuals^2 - sigma2) / (2 * sigma2 * sigma2)
    }
    g
  }

  o$hessian <- function(y, theta, expected = FALSE) {
    sigma2 <- theta[[2]]
    sigma4 <- sigma2 * sigma2
    if (expected) {
      list(
        mu_mu = -1 / sigma2,
        sigma2_sigma2 = -1 / (2 * sigma4),
        mu_sigma2 = 0
      )
    } else {
      residuals <- y - theta[[1]]
      list(
        mu_mu = -1 / sigma2,
        sigma2_sigma2 = (sigma2 - 2 * residuals^2) / (2 * sigma4 * sigma2),
        mu_sigma2 = -residuals / (sigma4)
      )
    }
  }

  o$grad_y <- function(y, theta) {
    -(y - theta[[1]]) / theta[[2]]
  }

  o$hess_y <- function(y, theta) {
    -1 / theta[[2]]
  }

  o$kernel <- function(y, theta, log = TRUE) {
    k <- -.5 * (y - theta[[1]])^2 / theta[[2]]
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

  o$mean <- o$median <- o$mode <- function(theta) {
    theta[[1]]
  }

  o$variance <- function(theta) {
    theta[[2]]
  }

  o$skewness <- function(theta) {
    0
  }

  o$kurtosis <- function(theta) {
    0
  }

  o$initialize <- function(y) {
    list(
      mu = mean(y),
      sigma2 = var(y)
    )
  }

  o
}





#' Gaussian `distrib` Object (Precision Parameterization)
#'
#' @description
#' Creates a distribution object for the Gaussian distribution parameterized by mean (\eqn{\mu}) and precision (\eqn{\tau = 1/\sigma^2}).
#'
#' @param link_mu A link function object for the location parameter \eqn{\mu}.
#'   Defaults to \code{\link[linkfunctions]{identity_link}}.
#' @param link_tau A link function object for the precision parameter \eqn{\tau}.
#'   Defaults to \code{\link[linkfunctions]{log_link}} to ensure positivity.
#'
#' @details
#' The Gaussian distribution with precision parameterization has the following density function:
#' \deqn{f(y; \mu, \tau) = \sqrt{\dfrac{\tau}{2\pi}} \exp\left\{-\dfrac{\tau}{2}(y-\mu)^2\right\}}
#'
#' \strong{Moments:}
#' \itemize{
#'   \item Expected value: \eqn{\mathbb{E}(y) = \mu}
#'   \item Variance: \eqn{\mathbb{V}(y) = \dfrac{1}{\tau}}
#'   \item Skewness: \eqn{\gamma_1 = 0}
#'   \item Excess Kurtosis: \eqn{\gamma_2 = 0}
#' }
#'
#' \strong{Parameter Domains:}
#' \itemize{
#'   \item \eqn{\mu \in (-\infty, +\infty)}
#'   \item \eqn{\tau \in (0, +\infty)}
#' }
#'
#' \strong{Explicit Formulas:}
#'
#' \emph{Gradient:}
#' The derivatives of the log-pdf \eqn{\ell} with respect to the parameters are:
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \tau(y - \mu)}
#' \deqn{\dfrac{\partial \ell}{\partial \tau} = \dfrac{1}{2\tau} - \dfrac{(y - \mu)^2}{2}}
#'
#' \emph{Observed Hessian:}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu^2} = -\tau}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \tau^2} = -\dfrac{1}{2\tau^2}}
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \partial \tau} = y - \mu}
#'
#' \emph{Expected Hessian:}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\tau}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \tau^2}\right] = -\dfrac{1}{2\tau^2}}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu \partial \tau}\right] = 0}
#'
#' \emph{Derivatives with respect to \eqn{y}:}
#' \deqn{\dfrac{\partial \ell}{\partial y} = -\tau(y - \mu)}
#' \deqn{\dfrac{\partial^2 \ell}{\partial y^2} = -\tau}
#'
#' @return A list of class \code{"distrib"} containing the components for the Gaussian distribution (precision parameterization).
#'
#' @importFrom linkfunctions identity_link log_link
#' @importFrom stats dnorm pnorm qnorm rnorm var sd
#' @export
gaussian3_distrib <- function(link_mu = identity_link(), link_tau = log_link()) {
  o <- list()
  class(o) <- c("distrib")

  o$distrib_name <- "gaussian"
  o$type <- "continuous"
  o$dimension <- 1
  o$bounds <- c(-Inf, Inf)

  o$params <- c("mu", "tau")
  o$params_interpretation <- c(mu = "mean", tau = "precision")
  o$n_params <- 2
  o$params_bounds <- list(
    mu = c(-Inf, Inf),
    tau = c(0, Inf)
  )
  o$link_params <- list(
    mu = link_mu,
    tau = link_tau
  )

  o$pdf <- function(y, theta, log = FALSE) {
    stats::dnorm(
      x = y,
      mean = theta[[1]],
      sd = 1 / sqrt(theta[[2]]),
      log = log
    )
  }

  o$cdf <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {
    stats::pnorm(
      q = q,
      mean = theta[[1]],
      sd = 1 / sqrt(theta[[2]]),
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$quantile <- function(p, theta, lower.tail = TRUE, log.p = FALSE) {
    stats::qnorm(
      p = p,
      mean = theta[[1]],
      sd = 1 / sqrt(theta[[2]]),
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$rng <- function(n, theta) {
    stats::rnorm(
      n = n,
      mean = theta[[1]],
      sd = 1 / sqrt(theta[[2]])
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
    tau <- theta[[2]]
    residuals <- y - theta[[1]]
    g <- list()
    if ("mu" %in% par) {
      g$mu <- tau * residuals
    }
    if ("tau" %in% par) {
      g$tau <- (1 / (2 * tau)) - (residuals^2 / 2)
    }
    g
  }

  o$hessian <- function(y, theta, expected = FALSE) {
    tau <- theta[[2]]
    if (expected) {
      list(
        mu_mu = -tau,
        tau_tau = -1 / (2 * tau^2),
        mu_tau = 0
      )
    } else {
      residuals <- y - theta[[1]]
      list(
        mu_mu = -tau,
        tau_tau = -1 / (2 * tau^2),
        mu_tau = residuals
      )
    }
  }

  o$grad_y <- function(y, theta) {
    -theta[[2]] * (y - theta[[1]])
  }

  o$hess_y <- function(y, theta) {
    -theta[[2]]
  }

  o$kernel <- function(y, theta, log = TRUE) {
    k <- -.5 * (y - theta[[1]])^2 * theta[[2]]
    if (log) {
      k
    } else {
      exp(k)
    }
  }

  o$normalization_constant <- function(theta, log = TRUE) {
    z <- sqrt(2 * pi / theta[[2]])
    if (log) {
      log(z)
    } else {
      z
    }
  }

  o$mean <- o$median <- o$mode <- function(theta) {
    theta[[1]]
  }

  o$variance <- function(theta) {
    1 / theta[[2]]
  }

  o$skewness <- function(theta) {
    0
  }

  o$kurtosis <- function(theta) {
    0
  }

  o$initialize <- function(y) {
    list(
      mu = mean(y),
      tau = 1 / var(y)
    )
  }

  o
}
