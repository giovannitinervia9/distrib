#' Zero-Adjusted Distribution Constructor (Discrete)
#'
#' @description
#' Creates a Zero-Adjusted (Hurdle) version of an existing **discrete** distribution object.
#' Unlike a Zero-Inflated distribution (which *adds* probability mass to zero via a mixture),
#' a Zero-Adjusted distribution *replaces* the probability mass at zero.
#'
#' @param distrib An object of class \code{"distrib"}. It must satisfy two conditions:
#'    1. \code{type == "discrete"} (e.g., Poisson, Negative Binomial).
#'    2. The support must include 0 (i.e., \code{bounds[1] <= 0}).
#' @param link_za A link function object for the zero probability parameter \eqn{\pi}.
#'    Defaults to \code{\link[linkfunctions]{logit_link}}.
#'
#' @details
#' **1. Definition and PDF:**
#' A Zero-Adjusted (Hurdle) distribution defines a process where the response \eqn{Y} is 0 with probability \eqn{\pi},
#' and follows a **zero-truncated** distribution \eqn{f_{tr}(y)} with probability \eqn{1-\pi}.
#'
#' Let \eqn{f(y; \theta)} be the PMF of the original distribution. The truncated PMF is:
#' \deqn{f_{tr}(y; \theta) = \dfrac{f(y; \theta)}{1 - f(0; \theta)} \quad \text{for } y > 0}
#'
#' The resulting Zero-Adjusted PMF is:
#' \deqn{
#' P(Y=y; \theta, \pi) =
#' \begin{cases}
#' \pi & \text{if } y = 0 \\
#' (1 - \pi) \dfrac{f(y; \theta)}{1 - f(0; \theta)} & \text{if } y > 0
#' \end{cases}
#' }
#'
#' **2. Cumulative Distribution Function (CDF):**
#' \deqn{
#' F_{ZA}(q) =
#' \begin{cases}
#' \pi + (1 - \pi) \dfrac{F(q; \theta) - f(0; \theta)}{1 - f(0; \theta)} & \text{if } q \ge 0 \\
#' 0 & \text{if } q < 0
#' \end{cases}
#' }
#'
#' **3. Quantile Function:**
#' Inversion involves solving for the truncated distribution when \eqn{p > \pi}.
#' \deqn{
#' Q_{ZA}(p) =
#' \begin{cases}
#' 0 & \text{if } p \le \pi \\
#' Q\left( u(1-f(0)) + f(0) ; \theta \right) & \text{if } p > \pi
#' \end{cases}
#' }
#' where \eqn{u = \dfrac{p - \pi}{1 - \pi}}.
#'
#' **4. Moments:**
#' The moments are derived from the moments of the zero-truncated distribution.
#' Since \eqn{0^k = 0} for \eqn{k \ge 1}, the zero mass contributes nothing to raw moments.
#' \deqn{\mathbb{E}_{ZA}[Y^k] = (1 - \pi) \mathbb{E}_{tr}[Y^k] = (1 - \pi) \dfrac{\mathbb{E}_{orig}[Y^k]}{1 - f(0; \theta)}}
#' Note: The expectation \eqn{\mathbb{E}_{orig}} in the numerator technically sums from \eqn{y=0}, but since \eqn{0^k=0}, it equals the sum from \eqn{y=1}.
#'
#' **5. Log-Likelihood and Gradient:**
#' The log-likelihood is separable into a binomial part for zero/non-zero, and a truncated part:
#' \deqn{\ell = \mathbb{I}(y=0)\log(\pi) + \mathbb{I}(y>0) \left[ \log(1-\pi) + \log f(y; \theta) - \log(1 - f(0; \theta)) \right]}
#'
#' The score vector \eqn{\nabla \ell} components are:
#'
#' * **For \eqn{\pi} (za):**
#'   \deqn{\dfrac{\partial \ell}{\partial \pi} = \mathbb{I}(y=0)\dfrac{1}{\pi} - \mathbb{I}(y>0)\dfrac{1}{1-\pi}}
#'
#' * **For \eqn{\theta} (original params):**
#'   For \eqn{y > 0}, a correction term arises from the normalization factor:
#'   \deqn{\dfrac{\partial \ell}{\partial \theta} = \nabla_\theta \log f(y; \theta) + \dfrac{f'(0; \theta)}{1 - f(0; \theta)}}
#'   For \eqn{y = 0}, the derivative w.r.t \eqn{\theta} is 0.
#'
#' **6. Hessian Matrix:**
#' Due to the separation of the likelihood, the parameters \eqn{\pi} and \eqn{\theta} are orthogonal (mixed derivatives are zero).
#'
#' * **Block \eqn{\theta\theta} (Correction):**
#'   The truncation introduces a correction term to the Hessian for \eqn{y > 0}:
#'   \deqn{H_{corr} = \dfrac{(1-f(0))f''(0) + (f'(0))^2}{(1-f(0))^2}}
#'   The observed Hessian for \eqn{y>0} is \eqn{H_{orig}(y) + H_{corr}}.
#'
#' @return A new object of class \code{"distrib"} representing the Zero-Adjusted distribution.
#'
#' @export
zero_adjusted_discrete <- function(distrib, link_za = logit_link()) {
  # --- 1. Validation ---
  if (!inherits(distrib, "distrib")) stop("Input must be a 'distrib' object.")
  if (distrib$type != "discrete") stop("This function requires a discrete distribution.")
  # Note: Strictly speaking, we need f(0) > 0 for the math to hold as written below
  if (distrib$bounds[1] > 0) stop("Original distribution must include 0 in support.")

  o <- list()
  class(o) <- "distrib"

  # --- 2. Metadata ---
  o$distrib_name <- paste0("zero-adjusted_", distrib$distrib_name)
  o$type <- "discrete"
  o$dimension <- distrib$dimension
  o$bounds <- c(0, distrib$bounds[2])

  o$params <- c(distrib$params, "za")
  o$n_params <- distrib$n_params + 1
  o$params_interpretation <- c(distrib$params_interpretation, za = "prob. zero")

  o$params_bounds <- distrib$params_bounds
  o$params_bounds$za <- c(0, 1)

  o$link_params <- distrib$link_params
  o$link_params$za <- link_za

  split_theta <- function(theta) {
    n <- length(theta)
    list(orig = theta[1:(n - 1)], za = theta[[n]])
  }

  # --- 3. PMF ---
  # P(Y=0) = za
  # P(Y=y) = (1-za) * f(y) / (1 - f(0))  per y > 0
  o$pdf <- function(y, theta, log = FALSE) {
    pars <- split_theta(theta)
    za <- pars$za

    pdf0 <- distrib$pdf(0, pars$orig, log = FALSE)
    # log(1 - f(0)) using log1p(-x) for precision when f(0) is small
    log_den <- log1p(-pdf0)

    # Log-PMF for positive part: log(1-za) + log(f(y)) - log(1-f(0))
    log_res_pos <- log(1 - za) + distrib$pdf(y, pars$orig, log = TRUE) - log_den

    n <- NROW(y)
    log_res <- numeric(n)

    is_zero <- (y == 0)
    if (any(is_zero)) log_res[is_zero] <- log(za)
    if (any(!is_zero)) log_res[!is_zero] <- log_res_pos[!is_zero]

    if (log) log_res else exp(log_res)
  }

  # --- 4. CDF ---
  o$cdf <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {
    pars <- split_theta(theta)
    za <- pars$za

    f0 <- distrib$pdf(0, pars$orig, log = FALSE)
    F_orig <- distrib$cdf(q, pars$orig, lower.tail = TRUE, log.p = FALSE)

    # Truncated CDF: (F(q) - f0) / (1 - f0)
    # Note: F_orig contains f0 mass, so we subtract it to re-normalize the rest
    cdf_trunc <- (F_orig - f0) / (1 - f0)
    cdf_trunc <- pmax(0, cdf_trunc)

    res <- za + (1 - za) * cdf_trunc

    res[q < 0] <- 0
    res <- pmin(res, 1)

    if (!lower.tail) res <- 1 - res
    if (log.p) log(res) else res
  }

  # --- 5. Quantile Function ---
  o$quantile <- function(p, theta, lower.tail = TRUE, log.p = FALSE) {
    pars <- split_theta(theta)
    za <- pars$za

    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1 - p
    p <- pmin(pmax(p, 0), 1)

    f0 <- distrib$pdf(0, pars$orig, log = FALSE)
    q_vals <- numeric(length(p))
    is_pos <- (p > za)

    if (any(is_pos)) {
      p_curr <- p[is_pos]
      za_curr <- if (length(za) > 1) za[is_pos] else za
      f0_curr <- if (length(f0) > 1) f0[is_pos] else f0
      pars_orig_curr <- lapply(pars$orig, function(x) {
        if (length(x) > 1) x[is_pos] else x
      })

      # Inversion logic:
      # p = za + (1-za) * F_trunc(q)
      # u = F_trunc(q) = (p - za)/(1 - za)
      # (F_orig(q) - f0)/(1-f0) = u
      # F_orig(q) = u*(1-f0) + f0
      u <- (p_curr - za_curr) / (1 - za_curr)
      target_prob <- u * (1 - f0_curr) + f0_curr

      # numerical safety clamp
      target_prob <- pmin(target_prob, 1 - 1e-10)

      q_vals[is_pos] <- distrib$quantile(target_prob, pars_orig_curr)
    }
    q_vals
  }

  # --- 6. Random Number Generator ---
  o$rng <- function(n, theta) {
    pars <- split_theta(theta)
    za <- pars$za

    # 1. Decide if zero
    is_zero <- stats::runif(n) < za
    y <- numeric(n)

    # 2. If not zero, generate from Truncated Distribution
    if (any(!is_zero)) {
      n_pos <- sum(!is_zero)
      f0_all <- distrib$pdf(0, pars$orig, log = FALSE)
      f0_curr <- if (length(f0_all) > 1) f0_all[!is_zero] else f0_all

      pars_orig_curr <- lapply(pars$orig, function(x) {
        if (length(x) > 1) x[!is_zero] else x
      })

      # Inverse Transform Sampling on Truncated Dist:
      # Generate U ~ Unif(0, 1)
      # Scale U to range [f(0), 1] of original CDF
      u <- stats::runif(n_pos)
      u_scaled <- f0_curr + u * (1 - f0_curr)

      y[!is_zero] <- distrib$quantile(u_scaled, pars_orig_curr)
    }

    y
  }

  o$loglik <- function(y, theta) {
    o$pdf(y, theta, log = TRUE)
  }

  # --- 7. Gradient ---
  o$gradient <- function(y, theta, par = NULL) {
    if (is.null(par)) par <- o$params
    pars <- split_theta(theta)
    za <- pars$za
    par_or <- par[par %in% distrib$params]

    # Auxiliary calculations at 0
    f0 <- distrib$pdf(0, pars$orig, log = FALSE)
    score_0 <- distrib$gradient(0, pars$orig, par_or) # Score S(0) = f'(0)/f(0)

    # Correction term for truncation
    # d/dTheta [ -log(1 - f(0)) ] = f'(0) / (1 - f(0))
    # = f(0) / (1 - f(0))
    correction_factor <- f0 / (1 - f0)

    grad_orig <- distrib$gradient(y, pars$orig, par_or)
    res_grad <- list()

    # A. Gradient w.r.t Original Params
    for (nm in names(grad_orig)) {
      if (nm %in% par) {
        # If y=0: gradient is 0 (params don't affect za)
        # If y>0: S(y) + Correction * S(0)
        term_pos <- grad_orig[[nm]] + correction_factor * score_0[[nm]]
        res_grad[[nm]] <- ifelse(y == 0, 0, term_pos)
      }
    }

    # B. Gradient w.r.t ZA
    if ("za" %in% par) {
      res_grad$za <- ifelse(y == 0, 1 / za, -1 / (1 - za))
    }

    res_grad
  }

  # --- 8. Hessian ---
  o$hessian <- function(y, theta, expected = FALSE) {
    pars <- split_theta(theta)
    th_orig <- pars$orig
    za <- pars$za

    # --- Pre-calculations (Independent of y) ---
    f0 <- distrib$pdf(0, th_orig, log = FALSE)
    grad_0 <- distrib$gradient(0, th_orig) # f'(0)/f(0)
    hess_0_obs <- distrib$hessian(0, th_orig, expected = FALSE) # H(0)
    denom <- 1 - f0
    denom2 <- denom^2

    res_hess <- list()

    # --- A. Block ZA-ZA ---
    if (expected) {
      res_hess[["za_za"]] <- -1 / (za * (1 - za))
    } else {
      # Observed
      val_0 <- -1 / (za^2)
      val_pos <- -1 / ((1 - za)^2)
      res_hess[["za_za"]] <- ifelse(y == 0, val_0, val_pos)
    }

    # --- B. Block Mixed (Always 0) ---
    for (nm in names(th_orig)) {
      res_hess[[paste0(nm, "_za")]] <- rep(0, NROW(y))
    }

    # --- C. Block Theta-Theta ---
    for (nm in names(hess_0_obs)) {
      parts <- strsplit(nm, "_")[[1]]
      p1 <- parts[1]
      p2 <- parts[length(parts)]

      # 1. Calculate Truncation Correction H_corr
      # H_corr = [ (1-f)*f'' + (f')^2 ] / (1-f)^2

      s1 <- grad_0[[p1]]
      s2 <- grad_0[[p2]]
      h_log_0 <- hess_0_obs[[nm]]

      f_prime_1 <- f0 * s1
      f_prime_2 <- f0 * s2
      f_prime_prime <- f0 * (h_log_0 + s1 * s2)
      hess_correction <- ((1 - f0) * f_prime_prime + f_prime_1 * f_prime_2) / denom2

      if (expected) {
        # E_za[H] = (1-za) * E_trunc[H_pos]
        # E_trunc[H_pos] = (E_orig[H] - f0*H_orig(0)) / (1-f0) + H_corr
        I_orig_full <- distrib$hessian(y, th_orig, expected = TRUE)[[nm]]
        E_trunc_H_orig <- (I_orig_full - f0 * h_log_0) / denom
        val_expected <- (1 - za) * (E_trunc_H_orig + hess_correction)

        res_hess[[nm]] <- val_expected
      } else {
        # Observed: 0 if y=0, H_orig + H_corr if y>0
        h_orig_y <- distrib$hessian(y, th_orig, expected = FALSE)[[nm]]
        term_pos <- h_orig_y + hess_correction
        res_hess[[nm]] <- ifelse(y == 0, 0, term_pos)
      }
    }

    expand_params(res_hess[hess_names(o$params)], length(y))
  }

  # --- 9. Moments ---
  o$raw_moment <- function(n, theta) {
    if (n == 0) {
      return(1)
    }
    pars <- split_theta(theta)
    za <- pars$za

    # E_orig[Y^n] sum starts at 0, but 0^n = 0, so it equals sum from 1.
    raw_orig <- moment(distrib, pars$orig, p = n, central = FALSE)
    f0 <- distrib$pdf(0, pars$orig, log = FALSE)

    # E_trunc = E_orig / (1-f0)
    trunc_moment <- raw_orig / (1 - f0)

    # E_za = (1-za) * E_trunc
    (1 - za) * trunc_moment
  }

  o$mean <- function(theta) {
    o$raw_moment(1, theta)
  }

  o$variance <- function(theta) {
    m1 <- o$raw_moment(1, theta)
    m2 <- o$raw_moment(2, theta)
    m2 - m1^2
  }

  o$skewness <- function(theta) {
    m1 <- o$raw_moment(1, theta)
    m2 <- o$raw_moment(2, theta)
    m3 <- o$raw_moment(3, theta)
    sigma <- sqrt(m2 - m1^2)
    (m3 - 3 * m1 * m2 + 2 * m1^3) / (sigma^3)
  }

  o$kurtosis <- function(theta) {
    m1 <- o$raw_moment(1, theta)
    m2 <- o$raw_moment(2, theta)
    m3 <- o$raw_moment(3, theta)
    m4 <- o$raw_moment(4, theta)
    var_val <- m2 - m1^2
    num <- m4 - 4 * m1 * m3 + 6 * (m1^2) * m2 - 3 * (m1^4)
    (num / (var_val^2)) - 3
  }

  o$initialize <- function(y) {
    prop_zeros <- mean(y == 0)
    za_init <- min(max(prop_zeros, 0.01), 0.99)
    y_pos <- y[y > 0]
    if (length(y_pos) == 0) y_pos <- y
    res <- distrib$initialize(y_pos)
    res$za <- za_init
    res
  }

  o
}




#' Zero-Adjusted Continuous Distribution Constructor
#'
#' @description
#' Creates a Zero-Adjusted (ZA) version of an existing **continuous** distribution object defined on \eqn{(0, \infty)}.
#' This model is often referred to as a "Zero-Inflated" continuous distribution (e.g., Zero-Inflated Gamma)
#' or "Hurdle" continuous model.
#'
#' It assumes the response variable \eqn{Y} follows a mixed distribution: a point mass at 0 with probability \eqn{p},
#' and a continuous distribution \eqn{W} (the original object) with probability \eqn{1-p}.
#'
#' @param distrib An object of class \code{"distrib"}. It must satisfy two conditions:
#'    1. \code{type == "continuous"} (e.g., Gamma, Lognormal).
#'    2. The domain is typically \eqn{(0, \infty)} (though it technically works for any continuous dist).
#' @param link_za A link function object for the zero probability parameter \eqn{p} (denoted as \code{za}).
#'    Defaults to \code{\link[linkfunctions]{logit_link}}.
#'
#' @details
#' **1. Definition and PDF:**
#' Let \eqn{Y \sim \text{ZAW}} denote the zero-adjusted variable.
#' Let \eqn{f_W(y; \theta)} be the PDF of the original continuous distribution \eqn{W} (where \eqn{W} has support \eqn{y > 0}).
#' Let \eqn{p \in [0, 1]} be the probability of zero (\code{za}).
#'
#' The PDF of \eqn{Y} is defined as (Rigby et al., Eq 9.1):
#' \deqn{
#' f_Y(y) =
#' \begin{cases}
#' p & \text{if } y = 0 \\
#' (1 - p) f_W(y; \theta) & \text{if } y > 0
#' \end{cases}
#' }
#'
#' **2. Cumulative Distribution Function (CDF):**
#' \deqn{
#' F_Y(q) =
#' \begin{cases}
#' 0 & \text{if } q < 0 \\
#' p + (1 - p) F_W(q; \theta) & \text{if } q \ge 0
#' \end{cases}
#' }
#' Note that \eqn{F_Y(0) = p} (assuming \eqn{F_W(0) = 0}, which is standard for Gamma/LogNormal).
#'
#' **3. Quantile Function:**
#' \deqn{
#' Q_Y(u) =
#' \begin{cases}
#' 0 & \text{if } u \le p \\
#' Q_W\left( \dfrac{u - p}{1 - p}; \theta \right) & \text{if } u > p
#' \end{cases}
#' }
#'
#' **4. Moments (GAMLSS Eq 9.2):**
#' The moments of \eqn{Y} are derived analytically from the raw and central moments of the original distribution \eqn{W}.
#' Let \eqn{\mu'_1 = \mathbb{E}[W]} be the mean of \eqn{W}, and \eqn{\mu_k = \mathbb{E}[(W - \mu'_1)^k]} be the \eqn{k}-th central moment of \eqn{W}.
#'
#' \itemize{
#'   \item **Mean:** \eqn{\mathbb{E}[Y] = (1-p)\mu'_1}
#'   \item **Variance:** \eqn{\mathbb{V}(Y) = (1-p)\mu_2 + p(1-p)(\mu'_1)^2}
#'   \item **Skewness:**
#'     \deqn{\mu_{3Y} = (1-p)[\mu_3 + 3p\mu_2\mu'_1 + p(2p-1)(\mu'_1)^3]}
#'     \eqn{\gamma_{1Y} = \mu_{3Y} / \mathbb{V}(Y)^{1.5}}
#'   \item **Excess Kurtosis:**
#'     \deqn{\mu_{4Y} = (1-p)[\mu_4 + 4p\mu_3\mu'_1 + 6p^2\mu_2(\mu'_1)^2 + p(1-3p+3p^2)(\mu'_1)^4]}
#'     \eqn{\gamma_{2Y} = \mu_{4Y} / \mathbb{V}(Y)^2 - 3}
#' }
#'
#' **5. Log-Likelihood and Gradient:**
#' The log-likelihood separates completely because \eqn{y=0} and \eqn{y>0} are disjoint events:
#' \deqn{\ell = \mathbb{I}(y=0)\log(p) + \mathbb{I}(y>0)[\log(1-p) + \log f_W(y; \theta)]}
#'
#' The parameters \eqn{p} and \eqn{\theta} are orthogonal.
#' \itemize{
#'   \item \eqn{\nabla_p \ell = \frac{1}{p}} if \eqn{y=0}, \eqn{\frac{-1}{1-p}} if \eqn{y>0}.
#'   \item \eqn{\nabla_\theta \ell = 0} if \eqn{y=0}, \eqn{\nabla_\theta \log f_W(y)} if \eqn{y>0}.
#' }
#'
#' **6. Hessian Matrix:**
#' \itemize{
#'   \item **Block \eqn{pp}**: \eqn{-1/p^2} (if \eqn{y=0}) or \eqn{-1/(1-p)^2} (if \eqn{y>0}).
#'   \item **Block \eqn{p\theta}**: Always 0.
#'   \item **Block \eqn{\theta\theta}**: 0 (if \eqn{y=0}) or \eqn{H_W(y)} (if \eqn{y>0}).
#' }
#'
#' **7. Expected Hessian:**
#' \itemize{
#'   \item **Block \eqn{pp}**: \eqn{\mathbb{E}[-\partial^2 \ell] = \frac{1}{p(1-p)}}.
#'   \item **Block \eqn{\theta\theta}**: \eqn{(1-p) \mathbb{E}[H_W]} (Fisher information of \eqn{W} scaled by probability of non-zero).
#' }
#'
#' @return A new object of class \code{"distrib"} representing the Zero-Adjusted Continuous distribution.
#'
#' @export
zero_adjusted_continuous <- function(distrib, link_za = logit_link()) {
  # --- 1. Validation ---
  if (!inherits(distrib, "distrib")) stop("Input must be a 'distrib' object.")
  if (distrib$type != "continuous") stop("This function requires a continuous distribution.")

  o <- list()
  class(o) <- "distrib"

  # --- 2. Metadata ---
  o$distrib_name <- paste0("zero-adjusted_", distrib$distrib_name)
  o$type <- "continuous" # Remains continuous (with a spike)
  o$dimension <- distrib$dimension
  o$bounds <- c(min(0, distrib$bounds[1]), distrib$bounds[2])

  o$params <- c(distrib$params, "za")
  o$n_params <- distrib$n_params + 1
  o$params_interpretation <- c(distrib$params_interpretation, za = "prob. zero")

  o$params_bounds <- distrib$params_bounds
  o$params_bounds$za <- c(0, 1)

  o$link_params <- distrib$link_params
  o$link_params$za <- link_za

  split_theta <- function(theta) {
    n <- length(theta)
    list(orig = theta[1:(n - 1)], za = theta[[n]])
  }

  # --- 3. PDF ---
  # P(Y=0) = za
  # f(y)   = (1-za) * f_orig(y)   for y != 0
  o$pdf <- function(y, theta, log = FALSE) {
    pars <- split_theta(theta)
    za <- pars$za

    # Calculate original log-pdf
    log_f_orig <- distrib$pdf(y, pars$orig, log = TRUE)

    n <- length(y)
    log_res <- numeric(n)

    is_zero <- (y == 0)

    # 1. Mass at zero
    if (any(is_zero)) {
      log_res[is_zero] <- log(za)
    }

    # 2. Continuous density scaled by (1-za)
    if (any(!is_zero)) {
      # log( (1-za) * f(y) ) = log(1-za) + log(f(y))
      log_res[!is_zero] <- log(1 - za) + log_f_orig[!is_zero]
    }

    if (log) log_res else exp(log_res)
  }

  # --- 4. CDF ---
  # F_za(q) = (1-za)F_orig(q) + za * I(q >= 0)
  o$cdf <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {
    pars <- split_theta(theta)
    za <- pars$za

    F_orig <- distrib$cdf(q, pars$orig, lower.tail = TRUE, log.p = FALSE)

    # Base component scaled
    res <- (1 - za) * F_orig

    # Add jump at 0
    # Note: For Gaussian, F_orig(0)=0.5. So jump happens from (1-za)*0.5 to (1-za)*0.5 + za.
    res[q >= 0] <- res[q >= 0] + za

    res <- pmin(pmax(res, 0), 1)

    if (!lower.tail) res <- 1 - res
    if (log.p) log(res) else res
  }

  # --- 5. Quantile Function ---
  o$quantile <- function(p, theta, lower.tail = TRUE, log.p = FALSE) {
    pars <- split_theta(theta)
    za <- pars$za

    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1 - p
    p <- pmin(pmax(p, 0), 1)

    # Calculate critical points for the jump at 0
    F0_orig <- distrib$cdf(0, pars$orig) # E.g., 0.5 for Gaussian, 0 for Gamma

    p_lower <- (1 - za) * F0_orig
    p_upper <- p_lower + za

    q_vals <- numeric(length(p))

    # Case A: Left Tail (e.g. Gaussian negatives) -> p < p_lower
    idx_left <- (p < p_lower)
    if (any(idx_left)) {
      # (1-za)F(q) = p  =>  F(q) = p / (1-za)
      p_trans <- p[idx_left] / (1 - za)

      # Parametri subsetting (per evitare warning di lunghezza)
      pars_sub <- lapply(pars$orig, function(x) if (length(x) > 1) x[idx_left] else x)
      q_vals[idx_left] <- distrib$quantile(p_trans, pars_sub)
    }

    # Case B: The Jump (Zero Mass) -> p_lower <= p <= p_upper
    # q_vals is already 0 initialized, so we do nothing.

    # Case C: Right Tail (Positives) -> p > p_upper
    idx_right <- (p > p_upper)
    if (any(idx_right)) {
      # (1-za)F(q) + za = p  =>  F(q) = (p - za) / (1 - za)
      p_curr <- p[idx_right]
      za_curr <- if (length(za) > 1) za[idx_right] else za

      p_trans <- (p_curr - za_curr) / (1 - za_curr)
      p_trans <- pmin(p_trans, 1) # Safety clamp

      pars_sub <- lapply(pars$orig, function(x) if (length(x) > 1) x[idx_right] else x)
      q_vals[idx_right] <- distrib$quantile(p_trans, pars_sub)
    }

    q_vals
  }

  # --- 6. RNG ---
  o$rng <- function(n, theta) {
    pars <- split_theta(theta)
    za <- pars$za

    is_zero <- stats::runif(n) < za
    y <- numeric(n)

    if (any(!is_zero)) {
      pars_sub <- lapply(pars$orig, function(x) if (length(x) > 1) x[!is_zero] else x)
      y[!is_zero] <- distrib$rng(sum(!is_zero), pars_sub)
    }
    y
  }

  o$loglik <- function(y, theta) {
    o$pdf(y, theta, log = TRUE)
  }

  # --- 7. Gradient ---
  o$gradient <- function(y, theta, par = NULL) {
    if (is.null(par)) par <- o$params
    pars <- split_theta(theta)
    za <- pars$za

    res_grad <- list()
    grad_orig <- distrib$gradient(y, pars$orig)

    for (nm in names(grad_orig)) {
      if (nm %in% par) {
        val <- grad_orig[[nm]]
        res_grad[[nm]] <- ifelse(y == 0, 0, val)
      }
    }
    if ("za" %in% par) {
      res_grad$za <- ifelse(y == 0, 1 / za, -1 / (1 - za))
    }

    res_grad
  }

  # --- 8. Hessian ---
  o$hessian <- function(y, theta, expected = FALSE) {
    pars <- split_theta(theta)
    za <- pars$za

    res_hess <- list()

    if (expected) {
      res_hess[["za_za"]] <- -1 / (za * (1 - za))
    } else {
      val_0 <- -1 / (za^2)
      val_pos <- -1 / ((1 - za)^2)
      res_hess[["za_za"]] <- ifelse(y == 0, val_0, val_pos)
    }

    for (nm in names(pars$orig)) {
      res_hess[[paste0(nm, "_za")]] <- rep(0, length(y))
    }

    if (expected) {
      h_exp_orig <- distrib$hessian(y, pars$orig, expected = TRUE)

      for (nm in names(h_exp_orig)) {
        res_hess[[nm]] <- (1 - za) * h_exp_orig[[nm]]
      }
    } else {
      h_obs_orig <- distrib$hessian(y, pars$orig, expected = FALSE)

      for (nm in names(h_obs_orig)) {
        res_hess[[nm]] <- ifelse(y == 0, 0, h_obs_orig[[nm]])
      }
    }

    expand_params(res_hess[hess_names(o$params)], length(y))
  }

  # --- 9. Moments ---
  o$raw_moment <- function(n, theta) {
    if (n == 0) {
      return(1)
    }
    pars <- split_theta(theta)
    (1 - pars$za) * moment(distrib, pars$orig, p = n, central = FALSE)
  }

  o$mean <- function(theta) {
    pars <- split_theta(theta)
    mu_orig <- distrib$mean(pars$orig)

    if (is.null(mu_orig)) {
      return(o$raw_moment(1, theta))
    }

    (1 - pars$za) * mu_orig
  }

  o$variance <- function(theta) {
    pars <- split_theta(theta)
    za <- pars$za # p

    var_orig <- distrib$variance(pars$orig) # mu_2 (Central moment of Y1)
    mu_orig <- distrib$mean(pars$orig) # mu'_1 (Raw moment of Y1)

    if (is.null(var_orig) || is.null(mu_orig)) {
      # Fallback numerical
      m1 <- o$raw_moment(1, theta)
      m2 <- o$raw_moment(2, theta)
      return(m2 - m1^2)
    }

    (1 - za) * var_orig + za * (1 - za) * mu_orig^2
  }

  o$skewness <- function(theta) {
    pars <- split_theta(theta)
    za <- pars$za # p

    mu1_p <- distrib$mean(pars$orig) # mu'_1
    mu2 <- distrib$variance(pars$orig) # mu_2
    skew_orig <- distrib$skewness(pars$orig)

    if (is.null(mu1_p) || is.null(mu2) || is.null(skew_orig)) {
      m1 <- o$raw_moment(1, theta)
      m2 <- o$raw_moment(2, theta)
      m3 <- o$raw_moment(3, theta)
      sigma <- sqrt(max(0, m2 - m1^2))
      if (sigma < 1e-12) {
        return(0)
      }
      return((m3 - 3 * m1 * m2 + 2 * m1^3) / sigma^3)
    }

    mu3 <- skew_orig * mu2^(1.5)

    term1 <- mu3
    term2 <- 3 * za * mu2 * mu1_p
    term3 <- za * (2 * za - 1) * mu1_p^3

    mu3_Y <- (1 - za) * (term1 + term2 + term3)

    var_Y <- o$variance(theta)
    mu3_Y / (var_Y^1.5)
  }

  o$kurtosis <- function(theta) {
    pars <- split_theta(theta)
    za <- pars$za # p

    mu1_p <- distrib$mean(pars$orig) # mu'_1
    mu2 <- distrib$variance(pars$orig) # mu_2

    skew_orig <- distrib$skewness(pars$orig)

    kurt_orig <- distrib$kurtosis(pars$orig)

    if (is.null(mu1_p) || is.null(mu2) || is.null(skew_orig) || is.null(kurt_orig)) {
      m1 <- o$raw_moment(1, theta)
      m2 <- o$raw_moment(2, theta)
      m3 <- o$raw_moment(3, theta)
      m4 <- o$raw_moment(4, theta)
      var_val <- max(0, m2 - m1^2)
      if (var_val < 1e-12) {
        return(0)
      }
      num <- m4 - 4 * m1 * m3 + 6 * (m1^2) * m2 - 3 * (m1^4)
      return((num / (var_val^2)) - 3)
    }

    mu3 <- skew_orig * mu2^(1.5)
    mu4 <- (kurt_orig + 3) * mu2^2

    term1 <- mu4
    term2 <- 4 * za * mu3 * mu1_p
    term3 <- 6 * za^2 * mu2 * mu1_p^2
    term4 <- za * (1 - 3 * za + 3 * za^2) * mu1_p^4

    mu4_Y <- (1 - za) * (term1 + term2 + term3 + term4)

    var_Y <- o$variance(theta)
    (mu4_Y / (var_Y^2)) - 3
  }

  o$initialize <- function(y) {
    prop_zeros <- mean(y == 0)
    za_init <- min(max(prop_zeros, 0.01), 0.99)
    y_pos <- y[y != 0]
    if (length(y_pos) == 0) y_pos <- y
    res <- distrib$initialize(y_pos)
    res$za <- za_init
    res
  }

  o
}




#' Zero-Adjusted Distribution Wrapper
#'
#' @description
#' A generic constructor for creating Zero-Adjusted (ZA) distributions.
#' This function automatically detects the type of the input distribution (discrete or continuous)
#' and dispatches the call to the appropriate specialized constructor:
#' \itemize{
#'   \item **Discrete:** Calls \code{\link{zero_adjusted_discrete}} to create a Hurdle model.
#'   \item **Continuous:** Calls \code{\link{zero_adjusted_continuous}} to create a Zero-Inflated continuous model (ZAW).
#' }
#'
#' @param distrib An object of class \code{"distrib"}.
#' @param link_za A link function object for the zero probability parameter.
#'   Defaults to \code{\link[linkfunctions]{logit_link}}.
#'
#' @return A new object of class \code{"distrib"} representing the zero-adjusted version of the input distribution.
#'
#' @seealso \code{\link{zero_adjusted_discrete}}, \code{\link{zero_adjusted_continuous}}
#'
#' @export
zero_adjusted <- function(distrib, link_za = logit_link()) {
  if (!inherits(distrib, "distrib")) {
    stop("Input must be a 'distrib' object.")
  }

  if (distrib$type == "discrete") {
    zero_adjusted_discrete(distrib, link_za)
  } else if (distrib$type == "continuous") {
    zero_adjusted_continuous(distrib, link_za)
  }
}
