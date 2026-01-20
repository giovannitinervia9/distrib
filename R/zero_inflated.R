#' Zero-Inflated Distribution Constructor (Discrete)
#'
#' @description
#' Creates a Zero-Inflated (ZI) version of an existing **discrete** distribution object.
#' A Zero-Inflated distribution is a mixture of a degenerate distribution at zero (point mass)
#' and a standard count distribution.
#'
#' @param distrib An object of class \code{"distrib"}. It must satisfy two conditions:
#'    1. \code{type == "discrete"} (e.g., Poisson, Negative Binomial, Binomial).
#'    2. The support must include 0 (i.e., \code{bounds[1] <= 0}).
#' @param link_zi A link function object for the zero-inflation probability parameter \eqn{\zeta}.
#'    Defaults to \code{\link[linkfunctions]{logit_link}}.
#'
#' @details
#' **1. Definition and PDF:**
#' Let \eqn{f(y; \theta)} be the probability mass function (PMF) of the original discrete distribution
#' with parameters \eqn{\theta}.
#' Let \eqn{\zeta \in [0, 1]} be the zero-inflation probability (structural zeros).
#'
#' The PMF of the Zero-Inflated variable \eqn{Y} is given by:
#' \deqn{
#' P(Y=y; \theta, \zeta) =
#' \begin{cases}
#' \zeta + (1 - \zeta)f(0; \theta) & \text{if } y = 0 \\
#' (1 - \zeta)f(y; \theta) & \text{if } y > 0
#' \end{cases}
#' }
#'
#' **2. Cumulative Distribution Function (CDF):**
#' The CDF \eqn{F_{ZI}(q)} is derived from the original CDF \eqn{F(q; \theta)}:
#' \deqn{F_{ZI}(q) = (1 - \zeta)F(q; \theta) + \mathbb{I}(q \ge 0)\zeta}
#'
#' **3. Quantile Function:**
#' The quantile function \eqn{Q_{ZI}(p)} is found by inverting the CDF.
#' Let \eqn{P(Y=0) = \zeta + (1-\zeta)f(0; \theta)} be the total probability at zero.
#' \deqn{
#' Q_{ZI}(p) =
#' \begin{cases}
#' 0 & \text{if } p \le \zeta + (1-\zeta)F(0; \theta) \\
#' Q\left(\dfrac{p - \zeta}{1 - \zeta}; \theta\right) & \text{if } p > \zeta + (1-\zeta)F(0; \theta)
#' \end{cases}
#' }
#' where \eqn{Q(\cdot; \theta)} is the quantile function of the original distribution.
#'
#' **4. Moments:**
#' The moments are derived from the raw moments of the original distribution.
#' Let \eqn{\mu'_k = \mathbb{E}_{orig}[Y^k]} be the \eqn{k}-th raw moment of the count component.
#' \itemize{
#'    \item \strong{Raw Moments:} \eqn{\mathbb{E}_{ZI}[Y^k] = (1 - \zeta)\mu'_k}
#'    \item \strong{Mean:} \eqn{\mathbb{E}[Y] = (1 - \zeta)\mathbb{E}_{orig}[Y]}
#'    \item \strong{Variance:} \eqn{\mathbb{V}(Y) = (1 - \zeta)\mathbb{E}_{orig}[Y^2] - [(1 - \zeta)\mathbb{E}_{orig}[Y]]^2}
#' }
#' Skewness and Kurtosis are computed numerically from the standardized central moments derived from these raw moments.
#'
#' **5. Log-Likelihood and Gradient (Score):**
#' Let \eqn{L_0 = \zeta + (1-\zeta)f(0; \theta)} be the likelihood contribution for an observation \eqn{y=0}.
#' The log-likelihood is:
#' \deqn{\ell = \mathbb{I}(y=0)\log(L_0) + \mathbb{I}(y>0)[\log(1-\zeta) + \log f(y; \theta)]}
#'
#' The score vector \eqn{\nabla \ell} has components:
#'
#' * **For the inflation parameter \eqn{\zeta}:**
#'      \deqn{\dfrac{\partial \ell}{\partial \zeta} = \mathbb{I}(y=0)\dfrac{1 - f(0; \theta)}{L_0} - \mathbb{I}(y>0)\dfrac{1}{1 - \zeta}}
#'
#' * **For original parameters \eqn{\theta}:**
#'      \deqn{\dfrac{\partial \ell}{\partial \theta} = \mathbb{I}(y=0)\dfrac{1 - \zeta}{L_0}\dfrac{\partial f(0; \theta)}{\partial \theta} + \mathbb{I}(y>0)\dfrac{\partial \log f(y; \theta)}{\partial \theta}}
#'
#' **6. Observed Hessian Matrix:**
#'
#' * **Block \eqn{\zeta\zeta} (Observed):**
#'      \deqn{\dfrac{\partial^2 \ell}{\partial \zeta^2} = -\mathbb{I}(y=0)\left[\dfrac{1 - f(0; \theta)}{L_0}\right]^2 - \mathbb{I}(y>0)\dfrac{1}{(1 - \zeta)^2}}
#'
#' * **Block \eqn{\zeta\theta} (Mixed Observed):**
#'      \deqn{\dfrac{\partial^2 \ell}{\partial \zeta \partial \theta} = -\mathbb{I}(y=0)\dfrac{f(0; \theta)}{L_0^2}\dfrac{\partial f(0; \theta)}{\partial \theta}}
#'      (Note: This term is 0 for \eqn{y > 0}).
#'
#' * **Block \eqn{\theta\theta} (Observed):**
#'      \deqn{
#'      \dfrac{\partial^2 \ell}{\partial \theta \partial \theta^T} =
#'      \begin{cases}
#'      \dfrac{1-\zeta}{L_0}\nabla^2 f(0) - \left(\dfrac{1-\zeta}{L_0}\right)^2 \nabla f(0) \nabla f(0)^T & \text{if } y = 0 \\
#'      \nabla^2 \log f(y; \theta) & \text{if } y > 0
#'      \end{cases}
#'      }
#'
#' **7. Expected Hessian (Fisher Information):**
#' The expected Hessian \eqn{\mathbb{E}[\mathcal{H}]} is calculated as the negative sum of score outer products across the full support.
#'
#' * **Block \eqn{\zeta\zeta}:**
#'      \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \zeta^2}\right] = -\left[ \dfrac{(1-f(0))^2}{L_0} + \dfrac{1-f(0)}{1-\zeta} \right]}
#'
#' * **Block \eqn{\zeta\theta}:**
#'      The cross-terms simplify analytically to:
#'      \deqn{\mathbb{E}\left[\dfrac{\partial^2 \ell}{\partial \zeta \partial \theta}\right] = -\dfrac{1}{L_0} \nabla_\theta f(0; \theta)}
#'
#' * **Block \eqn{\theta\theta}:**
#'      Computed by decomposing the expectation over \eqn{y=0} and \eqn{y>0}:
#'      \deqn{\mathbb{E}\left[\mathcal{H}_{\theta\theta}\right] = L_0 \mathcal{H}_{\theta\theta}^{\text{obs}}(0) + (1-\zeta) \left( \mathbb{E}_{orig}[\mathcal{H}_{orig}] - f(0)\mathcal{H}_{orig}^{\text{obs}}(0) \right)}
#'      where \eqn{\mathcal{H}_{orig}^{\text{obs}}(0)} is the observed Hessian of the count density evaluated at 0.
#'
#' @return A new object of class \code{"distrib"} representing the Zero-Inflated distribution.
#'
#' @export
zero_inflated <- function(distrib, link_zi = logit_link()) {
  # --- 1. Validation ---
  if (!inherits(distrib, "distrib")) {
    stop("Input must be a 'distrib' object.")
  }

  if (distrib$type != "discrete") {
    stop("zero_inflated() function can be used only of distrib object with type = 'discrete'")
  }

  if (distrib$bounds[1] > 0) {
    stop("zero_inflated() function can be used only of distrib object with 0 in their support")
  }

  o <- list()
  class(o) <- "distrib"

  # --- 2. Metadata ---
  o$distrib_name <- paste0("zero-inflated_", distrib$distrib_name)
  o$type <- distrib$type
  o$dimension <- distrib$dimension

  # Bounds are corrected to always include 0
  o$bounds <- c(min(0, distrib$bounds[1]), distrib$bounds[2])

  o$params <- c(distrib$params, "zi")
  o$n_params <- distrib$n_params + 1
  o$params_interpretation <- c(distrib$params_interpretation, zi = "prob. zero")

  o$params_bounds <- distrib$params_bounds
  o$params_bounds$zi <- c(0, 1)

  o$link_params <- distrib$link_params
  o$link_params$zi <- link_zi

  # Helper to split parameter vector into Original and ZI parts
  split_theta <- function(theta) {
    n <- length(theta)
    list(
      orig = theta[1:(n - 1)],
      zi = theta[[n]]
    )
  }


  # --- 3. Probability Density Function (PMF) ---
  o$pdf <- function(y, theta, log = FALSE) {
    pars <- split_theta(theta)
    th_orig <- pars$orig
    zi <- pars$zi

    # Base probability from count distribution
    res <- (1 - zi) * distrib$pdf(y, th_orig, log = FALSE)

    # Add point mass at 0
    is_zero <- (y == 0)
    if (any(is_zero)) {
      res[is_zero] <- zi + res[is_zero]
    }

    if (log) {
      log(res)
    } else {
      res
    }
  }

  # --- 4. Cumulative Distribution Function ---
  o$cdf <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {
    pars <- split_theta(theta)
    zi <- pars$zi

    # F_zi(q) = (1-zi)F_orig(q) + zi * I(q >= 0)
    res_orig <- distrib$cdf(q, pars$orig, lower.tail = TRUE, log.p = FALSE)
    res <- (1 - zi) * res_orig

    is_non_neg <- (q >= 0)
    if (any(is_non_neg)) {
      zi_vec <- if (length(zi) == 1) zi else zi[is_non_neg]
      res[is_non_neg] <- res[is_non_neg] + zi_vec
    }

    # Clamp for numerical safety
    res <- pmin(pmax(res, 0), 1)

    if (!lower.tail) {
      res <- 1 - res
    }
    if (log.p) {
      log(res)
    } else {
      res
    }
  }

  # --- 5. Quantile Function ---
  o$quantile <- function(p, theta, lower.tail = TRUE, log.p = FALSE) {
    pars <- split_theta(theta)
    th_orig <- pars$orig
    zi <- pars$zi

    # Normalize input
    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1 - p
    p <- pmin(pmax(p, 0), 1)

    # Calculate Probability Mass at Zero: P(Y=0)
    F0 <- distrib$cdf(0, th_orig)
    prob_at_zero <- zi + (1 - zi) * F0

    q_vals <- numeric(length(p))

    # If p <= prob_at_zero, the quantile is 0.
    # We only compute for p > prob_at_zero.
    idx_pos <- (p > prob_at_zero)

    if (any(idx_pos)) {
      p_curr <- p[idx_pos]
      zi_curr <- if (length(zi) == 1) zi else zi[idx_pos]

      # Invert mixture: F(q) = (p - zi) / (1 - zi)
      p_trans <- (p_curr - zi_curr) / (1 - zi_curr)

      # Clamp to avoid floating point overshoot (e.g. 1.000000001)
      p_trans <- pmin(pmax(p_trans, 0), 1)

      q_vals[idx_pos] <- distrib$quantile(p_trans, th_orig)
    }

    q_vals
  }

  # --- 6. Random Number Generator ---
  o$rng <- function(n, theta) {
    pars <- split_theta(theta)
    zi <- pars$zi
    # Structural zeros determined by Bernoulli(zi)
    is_structural_zero <- stats::runif(n) < zi
    # Generate from original distribution
    y <- distrib$rng(n, pars$orig)
    # Apply zeros
    y[is_structural_zero] <- 0
    y
  }

  o$loglik <- function(y, theta) {
    o$pdf(y, theta, log = TRUE)
  }

  # --- 7. Gradient (Score Function) ---
  o$gradient <- function(y, theta, par = NULL) {
    if (is.null(par)) par <- o$params
    pars <- split_theta(theta)
    zi <- pars$zi

    # Pre-calculate densities at 0 for weights
    f0 <- distrib$pdf(0, pars$orig, log = FALSE)
    l0 <- o$pdf(0, theta, log = FALSE) # zi + (1-zi)f0

    # Gradient of original distribution (returns unweighted scores)
    grad_orig <- distrib$gradient(y, pars$orig)

    # Weight w: Multiplier for original scores
    # If y=0:  (1-zi)f(0) / L0  (fraction of prob mass coming from count dist)
    # If y>0:  1
    w <- ifelse(y == 0, ((1 - zi) * f0) / l0, 1)

    res_grad <- list()

    # Apply weights to original gradients
    for (nm in names(grad_orig)) {
      if (nm %in% par) {
        res_grad[[nm]] <- w * grad_orig[[nm]]
      }
    }

    # Gradient for ZI parameter
    if ("zi" %in% par) {
      # If y=0: (1 - f(0)) / L0
      # If y>0: -1 / (1 - zi)
      res_grad$zi <- ifelse(y == 0, (1 - f0) / l0, -1 / (1 - zi))
    }

    res_grad
  }

  # --- 8. Hessian Matrix ---
  o$hessian <- function(y, theta, expected = FALSE) {
    pars <- split_theta(theta)
    th_orig <- pars$orig
    zi <- pars$zi

    # Probability mass at zero (original and ZI)
    f0 <- distrib$pdf(0, th_orig, log = FALSE)
    l0 <- o$pdf(0, theta, log = FALSE)

    res_hess <- list()

    # Pre-calculate derivatives at zero for corrections
    grad_0 <- distrib$gradient(0, th_orig)
    hess_0_obs <- distrib$hessian(0, th_orig, expected = FALSE)

    # Weights for y=0 case
    w0 <- ((1 - zi) * f0) / l0

    if (expected) {
      # --- CASE A: EXPECTED HESSIAN (Fisher Information) ---
      h_orig_exp <- distrib$hessian(y, th_orig, expected = TRUE)

      # 1. Block Theta-Theta
      for (nm in names(h_orig_exp)) {
        parts <- strsplit(nm, "_")[[1]]
        p1 <- parts[1]
        p2 <- parts[length(parts)]

        # Contribution at y=0 (Rank-1 update term)
        # Term: w0 * H_obs(0) + w0(1-w0) * S(0)S(0)'
        h_zi_0 <- w0 * hess_0_obs[[nm]] +
          (w0 * (1 - w0)) * grad_0[[p1]] * grad_0[[p2]]

        # Contribution from y > 0
        # E_pos = E[H | y>0] approx (E[H] - H(0)f(0))
        contrib_pos <- h_orig_exp[[nm]] - (hess_0_obs[[nm]] * f0)

        # Total: L0 * term(0) + (1-zi) * term(pos)
        res_hess[[nm]] <- l0 * h_zi_0 + (1 - zi) * contrib_pos
      }

      # 2. Block Theta-ZI (Mixed)
      # Formula: - S_orig(0) * f(0) / L0
      for (p in names(th_orig)) {
        # Note: grad_0 is Score (f'/f). So S_orig(0)*f0 = f'(0).
        # We calculate f'(0)/L0.
        # Since Fisher Info is negative of Expected Hessian, we negate the product of scores.
        # Score_theta(0) = (1-zi)f'/L0. Score_zi(0) = (1-f)/L0.
        # E[S_th S_zi] = L0 * [ (1-zi)f'/L0 * (1-f)/L0 ] + (1-zi) * E_pos[...]
        # To simplify, we implement the derived analytical form for discrete ZI:
        # H_mixed = - f'(0) / L0^2 * (probability mass L0) ... roughly
        # Precise implementation using gradient at zero:
        term_mixed <- -(f0 / l0) * grad_0[[p]]
        res_hess[[paste0(p, "_zi")]] <- rep(term_mixed, length(y)) # Constant expectation
      }

      # 3. Block ZI-ZI
      # Fisher Info for ZI parameter
      score_zi_0 <- (1 - f0) / l0
      score_zi_pos <- -1 / (1 - zi)
      fisher_zi <- l0 * (score_zi_0^2) + (1 - l0) * (score_zi_pos^2)
      res_hess[["zi_zi"]] <- -fisher_zi
    } else {
      # --- CASE B: OBSERVED HESSIAN (Newton-Raphson) ---
      h_orig_list <- distrib$hessian(y, th_orig, expected = FALSE)

      # Weight derivative for correction at 0
      w_prime <- ifelse(y == 0, w0 * (1 - w0), 0)

      # 1. Block Theta-Theta
      for (nm in names(h_orig_list)) {
        parts <- strsplit(nm, "_")[[1]]
        p1 <- parts[1]
        p2 <- parts[length(parts)]

        # Correction at 0: w*H(0) + w(1-w)*S(0)S(0)'
        val_y0 <- w0 * hess_0_obs[[nm]] + w_prime * grad_0[[p1]] * grad_0[[p2]]
        res_hess[[nm]] <- ifelse(y == 0, val_y0, h_orig_list[[nm]])
      }

      # 2. Block Theta-ZI
      # If y=0: - f'(0) / L0^2 = - (f(0)*S(0)) / L0^2
      # If y>0: 0
      factor_mix <- ifelse(y == 0, -f0 / (l0^2), 0)
      for (p in names(th_orig)) {
        res_hess[[paste0(p, "_zi")]] <- factor_mix * grad_0[[p]]
      }

      # 3. Block ZI-ZI
      val_zi_0 <- -((1 - f0)^2) / (l0^2)
      val_zi_pos <- -1 / ((1 - zi)^2)
      res_hess[["zi_zi"]] <- ifelse(y == 0, val_zi_0, val_zi_pos)
    }

    expand_params(res_hess[hess_names(o$params)], NROW(y))
  }

  # --- 9. Moments ---

  # Internal helper for raw moments: E[Y^n]
  # E_zi[Y^n] = (1-zi) * E_orig[Y^n]  (because 0^n = 0)
  o$raw_moment <- function(n, theta) {
    pars <- split_theta(theta)
    (1 - pars$zi) * moment(distrib, pars$orig, p = n, central = FALSE)
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
    # Skewness = (E[X^3] - 3uE[X^2] + 2u^3) / sigma^3
    sigma <- sqrt(m2 - m1^2)
    (m3 - 3 * m1 * m2 + 2 * m1^3) / (sigma^3)
  }

  o$kurtosis <- function(theta) {
    m1 <- o$raw_moment(1, theta)
    m2 <- o$raw_moment(2, theta)
    m3 <- o$raw_moment(3, theta)
    m4 <- o$raw_moment(4, theta)
    # Excess Kurtosis = (E[X^4] - 4uE[X^3] + 6u^2E[X^2] - 3u^4) / sigma^4 - 3
    var_val <- m2 - m1^2
    num <- m4 - 4 * m1 * m3 + 6 * (m1^2) * m2 - 3 * (m1^4)
    (num / (var_val^2)) - 3
  }

  # --- 10. Initialization ---
  o$initialize <- function(y) {
    # 1. Estimate zi from proportion of zeros
    prop_zeros <- mean(y == 0)
    zi_init <- min(max(prop_zeros, 0.01), 0.99)

    # 2. Initialize original distribution using only positive values
    # (Assumption: Positive values roughly follow the count distribution)
    y_pos <- y[y > 0]
    if (length(y_pos) == 0) y_pos <- y # Fallback if all zeros

    res <- distrib$initialize(y_pos)
    res$zi <- zi_init
    res
  }

  o
}
