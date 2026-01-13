#' Beta `distrib` Object (Mean-Precision Parameterization)
#'
#' @description
#' Creates a distribution object for the Beta distribution parameterized by mean (\eqn{\mu}) and precision (\eqn{\phi}).
#'
#' @param link_mu A link function object for the mean parameter \eqn{\mu}.
#'   Defaults to \code{\link[linkfunctions]{logit_link}} to ensure the parameter stays within (0, 1).
#' @param link_phi A link function object for the precision parameter \eqn{\phi}.
#'   Defaults to \code{\link[linkfunctions]{log_link}} to ensure positivity.
#'
#' @details
#' The Beta distribution is reparameterized from the standard shape parameters \eqn{\alpha} and \eqn{\beta} using:
#' \deqn{\alpha = \mu\phi, \quad \beta = (1-\mu)\phi}
#' where \eqn{\phi = \alpha + \beta} acts as a precision parameter (or sample size proxy).
#'
#' The probability density function is:
#' \deqn{f(y; \mu, \phi) = \dfrac{\Gamma(\phi)}{\Gamma(\mu\phi)\Gamma((1-\mu)\phi)} y^{mu\phi-1} (1-y)^{(1-\mu)\phi-1}}
#' for \eqn{y \in (0, 1)}.
#'
#' \strong{Moments:}
#' \itemize{
#'   \item Expected value: \eqn{\mathbb{E}(y) = \mu}
#'   \item Variance: \eqn{\mathbb{V}(y) = \dfrac{\mu(1-\mu)}{1+\phi}}
#'   \item Skewness: \eqn{\gamma_1 = \dfrac{2(1-2\mu)\sqrt{1+\phi}}{(2+\phi)\sqrt{\mu(1-\mu)}}}
#'   \item Excess Kurtosis: \eqn{\gamma_2 = \dfrac{6[\mu(1-\mu)(1+\phi) + \phi(1+\phi)(1-2\mu)^2]}{\mu(1-\mu)(2+\phi)(3+\phi)} - 3}
#' }
#'
#' \strong{Parameter Domains:}
#' \itemize{
#'   \item \eqn{\mu \in (0, 1)}
#'   \item \eqn{\phi \in (0, +\infty)}
#' }
#'
#' \strong{Explicit Formulas:}
#'
#' Let \eqn{\psi(z)} be the digamma function and \eqn{\psi_1(z)} be the trigamma function.
#'
#' \emph{Gradient:}
#' \deqn{\dfrac{\partial \ell}{\partial \mu} = \phi \left[ \log\left(\dfrac{y}{1-y}\right) - \psi(\mu\phi) + \psi((1-\mu)\phi) \right]}
#' \deqn{\dfrac{\partial \ell}{\partial \phi} = \psi(\phi) - \mu\psi(\mu\phi) - (1-\mu)\psi((1-\mu)\phi) + \mu \log(y) + (1-\mu) \log(1-y)}
#'
#' \emph{Expected Hessian:}
#' When \code{expected = TRUE}, the expected second derivatives (Fisher Information) are:
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu^2}\right] = -\phi^2 \left[ \psi_1(\mu\phi) + \psi_1((1-\mu)\phi) \right]}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \phi^2}\right] = \psi_1(\phi) - \mu^2\psi_1(\mu\phi) - (1-\mu)^2\psi_1((1-\mu)\phi)}
#' \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \mu \partial \phi}\right] = -\phi \left[ \mu\psi_1(\mu\phi) - (1-\mu)\psi_1((1-\mu)\phi) \right]}
#'
#' \emph{Observed Hessian:}
#' The diagonal terms are identical to the expected Hessian.
#' The mixed term depends on \eqn{y}:
#' \deqn{\dfrac{\partial^2 \ell}{\partial \mu \partial \phi} = \log\left(\dfrac{y}{1-y}\right) - \psi(\mu\phi) + \psi((1-\mu)\phi) - \phi \left[ \mu\psi_1(\mu\phi) - (1-\mu)\psi_1((1-\mu)\phi) \right]}
#'
#' @return A list of class \code{"distrib"} containing the components for the Beta distribution.
#'
#' @importFrom linkfunctions logit_link log_link
#' @importFrom stats dbeta pbeta qbeta rbeta
#' @export
beta_distrib <- function(link_mu = logit_link(), link_phi = log_link()) {
  o <- list()
  class(o) <- c("distrib")

  o$distrib_name <- "beta"
  o$type <- "continuous"
  o$dimension <- 1
  o$bounds <- c(0, 1)

  o$params <- c("mu", "phi")
  o$params_interpretation <- c(mu = "mean", phi = "precision")
  o$n_params <- 2
  o$params_bounds <- list(
    mu = c(0, 1),
    phi = c(0, Inf)
  )
  o$link_params <- list(
    mu = link_mu,
    phi = link_phi
  )

  o$pdf <- function(y, theta, log = FALSE) {
    mu <- theta[[1]]
    phi <- theta[[2]]
    stats::dbeta(
      x = y,
      shape1 = mu * phi,
      shape2 = (1 - mu) * phi,
      log = log
    )
  }

  o$cdf <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {
    mu <- theta[[1]]
    phi <- theta[[2]]
    stats::pbeta(
      q = q,
      shape1 = mu * phi,
      shape2 = (1 - mu) * phi,
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$quantile <- function(p, theta, lower.tail = TRUE, log.p = FALSE) {
    mu <- theta[[1]]
    phi <- theta[[2]]
    stats::qbeta(
      p = p,
      shape1 = mu * phi,
      shape2 = (1 - mu) * phi,
      lower.tail = lower.tail,
      log.p = log.p
    )
  }

  o$rng <- function(n, theta) {
    mu <- theta[[1]]
    phi <- theta[[2]]
    stats::rbeta(
      n = n,
      shape1 = mu * phi,
      shape2 = (1 - mu) * phi
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
    phi <- theta[[2]]
    alpha <- mu * phi
    beta <- (1 - mu) * phi

    digamma_alpha <- digamma(alpha)
    digamma_beta <- digamma(beta)

    log_y <- log(y)
    log_1_y <- log(1 - y)
    log_ratio <- log_y - log_1_y # log(y / (1-y))

    g <- list()

    if ("mu" %in% par) {
      g$mu <- phi * (log_ratio - digamma_alpha + digamma_beta)
    }
    if ("phi" %in% par) {
      g$phi <- digamma(phi) - mu * digamma_alpha - (1 - mu) * digamma_beta + mu * log_y + (1 - mu) * log_1_y
    }
    g
  }

  o$hessian <- function(y, theta, expected = FALSE) {
    mu <- theta[[1]]
    phi <- theta[[2]]

    alpha <- mu * phi
    beta <- (1 - mu) * phi

    digamma_alpha <- digamma(alpha)
    digamma_beta <- digamma(beta)
    trigamma_alpha <- trigamma(alpha)
    trigamma_beta <- trigamma(beta)
    trigamma_phi <- trigamma(phi)

    # Common terms for diagonal
    hess_mu_mu <- -phi^2 * (trigamma_alpha + trigamma_beta)
    hess_phi_phi <- trigamma_phi - mu^2 * trigamma_alpha - (1 - mu)^2 * trigamma_beta

    if (expected) {
      # Expected mixed term simplifies greatly as E[log(y/(1-y))] cancels the digammas
      hess_mu_phi <- -phi * (mu * trigamma_alpha - (1 - mu) * trigamma_beta)
    } else {
      log_ratio <- log(y / (1 - y))
      # Observed mixed term
      term1 <- log_ratio - digamma_alpha + digamma_beta
      term2 <- phi * (mu * trigamma_alpha - (1 - mu) * trigamma_beta)
      hess_mu_phi <- term1 - term2
    }

    list(
      mu_mu = hess_mu_mu,
      phi_phi = hess_phi_phi,
      mu_phi = hess_mu_phi
    )
  }

  o$kernel <- function(y, theta, log = TRUE) {
    mu <- theta[[1]]
    phi <- theta[[2]]
    alpha <- mu * phi
    beta <- (1 - mu) * phi

    k <- (alpha - 1) * log(y) + (beta - 1) * log(1 - y)

    if (log) {
      k
    } else {
      exp(k)
    }
  }

  o$normalization_constant <- function(theta, log = TRUE) {
    mu <- theta[[1]]
    phi <- theta[[2]]
    z <- lbeta(mu * phi, (1 - mu) * phi)
    if (log) {
      z
    } else {
      exp(z)
    }
  }

  o$mean <- function(theta) {
    theta[[1]]
  }

  o$median <- function(theta) {
    o$quantile(.5, theta)
  }

  o$mode <- function(theta) {
    mu <- theta[[1]]
    phi <- theta[[2]]
    alpha <- mu * phi
    beta <- (1 - mu) * phi
    ifelse(alpha > 1 & beta > 1, (alpha - 1) / (alpha + beta - 2),
      ifelse(alpha <= 1 & beta > 1, 0,
        ifelse(alpha > 1 & beta <= 1, 1, NA)
      )
    )
  }

  o$variance <- function(theta) {
    mu <- theta[[1]]
    phi <- theta[[2]]
    (mu * (1 - mu)) / (1 + phi)
  }

  o$skewness <- function(theta) {
    mu <- theta[[1]]
    phi <- theta[[2]]
    num <- 2 * (1 - 2 * mu) * sqrt(1 + phi)
    den <- (2 + phi) * sqrt(mu * (1 - mu))
    num / den
  }

  o$kurtosis <- function(theta) {
    # Excess kurtosis
    mu <- theta[[1]]
    phi <- theta[[2]]
    term1 <- mu * (1 - mu) * (1 + phi)
    term2 <- phi * (1 + phi) * (1 - 2 * mu)^2
    denom <- mu * (1 - mu) * (2 + phi) * (3 + phi)
    (6 * (term1 + term2) / denom) - 3
  }

  o$initialize <- function(y) {
    mu <- mean(y)
    list(
      mu = mu,
      phi = mu * (1 - mu) / var(y) - 1
    )
  }

  o
}
