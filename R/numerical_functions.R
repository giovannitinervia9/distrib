#' Derivatives of the Modified Bessel Function of the Second Kind
#'
#' Calculates the n-th order derivative of the modified Bessel function of the second kind, \eqn{K_\nu(x)}.
#'
#' @param x Numeric vector. The input argument (must be > 0).
#' @param nu Numeric vector. The order of the Bessel function.
#' @param deriv Integer. The order of the derivative to compute (must be \eqn{\ge 0}). Defaults to 1.
#' @param expon.scaled Logical. Controls exponential scaling.
#'        \itemize{
#'          \item If \code{mode = "standard"}: If \code{TRUE}, returns \eqn{e^x \cdot K^{(n)}_\nu(x)}. If \code{FALSE}, returns \eqn{K^{(n)}_\nu(x)}.
#'          \item If \code{mode = "deriv_scaled"}: This argument is ignored because the target function is inherently scaled.
#'        }
#' @param mode Character string. Determines the type of derivative to compute:
#'        \itemize{
#'          \item \code{"standard"} (Default): Computes the derivative of the Bessel function itself.
#'          \item \code{"deriv_scaled"}: Computes the derivative of the scaled Bessel function.
#'        }
#'
#' @details
#' The function implements two different mathematical approaches based on \code{mode}:
#'
#' \strong{1. Standard Mode (\code{mode = "standard"}):} \cr
#' Calculates \eqn{\frac{d^n}{dx^n} K_\nu(x)}. \cr
#' It uses the recurrence relation:
#' \deqn{K^{(n)}_\nu(x) = (-1/2)^n \sum_{k=0}^{n} \binom{n}{k} K_{\nu - n + 2k}(x)}
#' If \code{expon.scaled = TRUE}, the result is multiplied by \eqn{e^x} to prevent underflow for large \eqn{x}.
#'
#' \strong{2. Scaled Derivative Mode (\code{mode = "deriv_scaled"}):} \cr
#' Calculates the derivative of the product \eqn{e^x K_\nu(x)}.
#' \deqn{\frac{d^n}{dx^n} \left( e^x K_\nu(x) \right)}
#' Using the General Leibniz rule, this expands to:
#' \deqn{e^x \sum_{j=0}^{n} \binom{n}{j} K^{(j)}_\nu(x)}
#'
#' @return A numeric vector of the same length as \code{x} (or \code{nu}).
#'
#' @examples
#' # --- Case 1: Standard derivative K'(x) ---
#' # Analytical derivative
#' dy <- dbesselK(x = 5, nu = 1, deriv = 1, mode = "standard")
#'
#' # --- Case 2: Scaled derivative for large x ---
#' # Returns e^x * K'(x)
#' dy_scaled <- dbesselK(x = 1000, nu = 1, deriv = 1, expon.scaled = TRUE)
#'
#' # --- Case 3: Derivative of the scaled function (numDeriv compat) ---
#' # We want d/dx (e^x * K(x))
#' library(numDeriv)
#' x_val <- 2
#'
#' # Numerical result
#' f_sc <- function(x) besselK(x, 1, expon.scaled = TRUE)
#' num_res <- grad(f_sc, x_val)
#'
#' # Analytical result using mode = "deriv_scaled"
#' ana_res <- dbesselK(x_val, 1, deriv = 1, mode = "deriv_scaled")
#'
#' c(Numerical = num_res, Analytical = ana_res)
#'
#' @export
dbesselK <- function(x, nu, deriv = 1L, expon.scaled = FALSE, mode = c("standard", "deriv_scaled")) {
  if (deriv < 0) {
    stop("Derivative order must be >= 0.")
  }

  mode <- match.arg(mode)

  if (deriv == 0) {
    return(besselK(x, nu, expon.scaled = if (mode == "deriv_scaled") TRUE else expon.scaled))
  }

  n <- deriv

  # Determine loop sequence:
  # - deriv_scaled: 0 to n (Leibniz rule requires sum of all lower order derivatives).
  # - standard: just n (Direct calculation of the n-th derivative).
  j_seq <- if (mode == "deriv_scaled") {
    0:n
  } else {
    n
  }

  res <- 0

  for (j in j_seq) {
    # Compute pure scaled derivative of order j
    # Formula uses: (-0.5)^j * sum( binomial * K_nu_adjusted )
    s <- 0
    for (k in 0:j) {
      s <- s + choose(j, k) * besselK(x, nu - j + 2 * k, expon.scaled = TRUE)
    }
    val_j <- ((-0.5)^j) * s

    # If Leibniz (deriv_scaled): weighted accumulation of terms.
    # If Standard: direct assignment (loop runs only once).
    if (mode == "deriv_scaled") {
      res <- res + choose(n, j) * val_j
    } else {
      res <- val_j
    }
  }

  # Remove scaling only if mode is "standard" and user did NOT request scaling.
  # (In "deriv_scaled", the result is naturally scaled by e^x due to the Leibniz expansion).
  if (mode == "standard" && !expon.scaled) {
    res * exp(-x)
  } else {
    res
  }
}




#' Numerical Summation of Discrete Series
#'
#' Calculates the sum of a function `f(x)` over a sequence of integers from `start` to `end`.
#' The function is designed to handle finite sums, one-sided infinite series, and
#' doubly infinite series (from `-Inf` to `Inf`) by automatically adapting its summation strategy.
#'
#' @details
#' **1. Summation Strategies:**
#' The function automatically detects the domain topology based on `start` and `end`:
#' * **Forward (Standard):** If `start <= end` (e.g., `1` to `Inf`), it sums \eqn{f(start) + f(start+1) + \dots} directly.
#' * **Backward (Reflection):** If `start > end` (e.g., `-1` to `-Inf`), it reflects the domain internally.
#'   It calculates \eqn{\sum f(-x)} effectively iterating forward through the negative numbers.
#' * **Doubly Infinite (Folding):** If `start == -Inf` and `end == Inf`, it uses a "folding" strategy around zero.
#'   It first calculates \eqn{f(0)}, then iterates \eqn{x} from 1 to \eqn{\infty}, summing
#'   \eqn{f(x) + f(-x)} at each step.
#'
#' **2. Vectorization and Convergence:**
#' The summation is performed in blocks of size `step` to leverage R's vectorization capabilities.
#' For infinite bounds, the function checks for convergence based on the `reltol` parameter.
#'
#' **3. Underflow Protection:**
#' The function includes logic to prevent premature stopping if the series starts with a sequence
#' of zeros (e.g., a distribution centered far from the start). It continues searching until
#' a significant sum is accumulated (`s > tol`).
#'
#' @param f A function taking a vector of integers `x` and returning a vector of numeric values.
#'   **Must be vectorized** (able to handle multiple inputs at once).
#' @param start Numeric. The starting value of the sequence. Can be finite, `Inf`, or `-Inf`.
#' @param end Numeric. The ending value of the sequence. Can be finite, `Inf`, or `-Inf`.
#'   Defaults to `Inf`.
#' @param step Integer. The number of terms to calculate in a single vectorized batch.
#'   Defaults to 10000.
#' @param tol Numeric. Tolerance threshold for convergence. Defaults to `1e-10`.
#' @param maxit Integer. Safety limit for the maximum number of batch iterations.
#'   Defaults to 1000000.
#' @param reltol Logical. If `TRUE` (default), uses a **hybrid relative tolerance** criterion:
#'   `tol*max(abs(s), 1)`.  If `FALSE`, uses strict **absolute tolerance** (\code{tol}).
#'
#' @return A numeric scalar representing the calculated sum.
#'
#' @examples
#' # --- Case 1: Standard Convergent Series (Basel Problem) ---
#' # Sum of 1/x^2 from 1 to Inf => pi^2 / 6
#' f_basel <- function(x) 1 / x^2
#' series(f_basel, start = 1, end = Inf)
#'
#' # --- Case 2: Expectation of Poisson (Large Mean) ---
#' d <- poisson_distrib()
#' theta <- list(mu = 10000)
#' f_mean <- function(x) x * d$pdf(x, theta)
#' series(f_mean, start = 0, end = Inf, reltol = TRUE)
#'
#' @export
series <- function(f, start = 0, end = Inf, step = 1000, tol = 1e-10, maxit = 10000000, reltol = TRUE) {
  # --- Setup Range and Direction ---
  if (is.infinite(start) && start < 0 && is.infinite(end) && end > 0) {
    s <- f(0)
    start_internal <- 1
    end_internal <- Inf
    f_internal <- function(x) f(x) + f(-x)
  } else if (end < start) {
    s <- 0
    start_internal <- -start
    end_internal <- -end
    f_internal <- function(x) f(-x)
  } else {
    s <- 0
    start_internal <- start
    end_internal <- end
    f_internal <- f
  }

  # --- Loop Initialization ---
  # Initial eps set high to ensure loop entry
  eps <- 2 * tol
  climbing <- TRUE
  upper_limit <- min(start_internal + step, end_internal)
  it <- 0

  # --- Main Loop ---
  while (climbing && it < maxit) {
    it <- it + 1

    # Vectorized calculation
    x <- start_internal:upper_limit
    vals <- f_internal(x)
    eps <- sum(vals)
    s <- s + eps

    # STOPPING CRITERIA:
    # 1. Define tolerance type based on user input
    if (reltol) {
      # Hybrid/Relative: scales with the magnitude of the sum 's'
      scaled_tol <- tol * max(abs(s), 1)
    } else {
      # Absolute: strict fixed threshold
      scaled_tol <- tol
    }

    # 2. Check for convergence AND Underflow protection (s > tol)
    if (eps < scaled_tol && s > tol) {
      climbing <- FALSE
    } else {
      # Advancement
      start_internal <- upper_limit + 1

      if (start_internal > end_internal) {
        climbing <- FALSE
      } else {
        upper_limit <- min(start_internal + step, end_internal)
      }
    }
  }

  if (it >= maxit) {
    warning(sprintf(
      "Max iterations (%d) reached. Sum = %.6e, last eps = %.6e",
      maxit, s, eps
    ))
  }

  s
}



#' Calculate Moments of a Distribution
#'
#' A generic function to calculate the \eqn{p}-th moment (raw or central) of a probability distribution.
#' It automatically handles continuous distributions via numerical integration and discrete distributions
#' via series summation.
#'
#' @details
#' **Definitions:**
#' * **Raw Moment:** If `central = FALSE`, calculates \eqn{E[X^p] = \int x^p f(x) dx} (or \eqn{\sum x^p P(x)}).
#' * **Central Moment:** If `central = TRUE`, calculates \eqn{E[(X - \mu)^p] = \int (x - \mu)^p f(x) dx} (or \eqn{\sum (x - \mu)^p P(x)}),
#'     where \eqn{\mu} is the mean of the distribution.
#'
#' **Mechanism:**
#' * **Continuous Distributions:** Uses \code{\link[stats]{integrate}}. To ensure numerical stability for
#'   distributions with sharp peaks on infinite domains, the function
#'   attempts to identify the **mode** of the distribution. If found, the integration interval is split
#'   at the mode (\eqn{[lb, mode] + [mode, ub]}). This "anchoring" forces the numerical integrator to
#'   evaluate the function in the region of highest density, preventing underflow errors where the integral
#'   would otherwise be incorrectly estimated as zero.
#' * **Discrete Distributions:** Uses the internal function \code{series} to sum the terms.
#'
#' **Robustness:**
#' The integrand is protected against non-finite values (NaN/Inf) which are coerced to 0, assuming they occur
#' in the negligible tails of the distribution.
#'
#' **Vectorization:**
#' The function fully supports vectorized parameters (`theta`). If `theta` contains vectors of length \eqn{n},
#' the function returns a vector of \eqn{n} moments. The order `p` can also be a vector of length \eqn{n}
#' (matching `theta`) or a scalar.
#'
#' @param distrib An object of class \code{"distrib"} containing the distribution definition
#'   (pdf, bounds, type, etc.).
#' @param theta A named list of parameters for the distribution. Vectors are supported for batch calculation.
#' @param p Numeric. The order of the moment to calculate. Defaults to 2.
#' @param central Logical. If \code{TRUE}, calculates the central moment (about the mean).
#'   If \code{FALSE} (default), calculates the raw moment (about zero).
#' @param mu Numeric (optional). The mean value used for central moment calculations.
#'   If \code{NULL} (default) and \code{central = TRUE}, the function attempts to retrieve
#'   the mean from the `distrib` object or calculates it numerically by calling `moment(..., p=1)`.
#' @param ... Additional arguments passed to the underlying numerical engine:
#'   * For continuous: passed to \code{\link[stats]{integrate}} (e.g., `subdivisions`, `rel.tol`).
#'   * For discrete: passed to `series` (e.g., `step`, `tol`, `maxit`).
#'
#' @return A numeric vector containing the calculated moments.
#'
#' @importFrom stats integrate
#' @export
moment <- function(
    distrib,
    theta,
    p = 2,
    central = FALSE,
    mu = NULL,
    ...) {
  type <- distrib$type
  bounds <- distrib$bounds
  pdf <- distrib$pdf
  params <- distrib$params

  max_dim_theta <- max(lengths(theta))
  check_params_dim(theta)

  dim_p <- length(p)
  if (dim_p != 1 && max_dim_theta != 1 && dim_p != max_dim_theta) {
    stop("Dimension of 'p' does not match dimension of 'theta'.")
  }

  theta <- expand_params(theta)
  theta_list <- transpose_params(theta)

  if (central) {
    if (is.null(mu)) {
      try_mu <- distrib$mean(theta)
      if (is.null(try_mu)) {
        mu <- moment(distrib, theta, p = 1, central = FALSE, ...)
      } else {
        mu <- try_mu
      }
    }
  } else {
    mu <- 0
  }

  f <- function(x, mu, p) {
    (x - mu)^p
  }

  if (type == "continuous") {
    FUN <- function(theta, mu, p) {
      peak <- tryCatch(distrib$mode(theta), error = function(e) NA)
      use_split <- !is.na(peak) && peak > bounds[1] && peak < bounds[2]
      safe_integrand <- function(t) {
        val <- f(t, mu, p) * pdf(t, theta, log = FALSE)
        val[!is.finite(val)]
        val
      }

      tryCatch(
        {
          if (use_split) {
            int_left <- stats::integrate(safe_integrand, lower = bounds[1], upper = peak, ...)$value
            int_right <- stats::integrate(safe_integrand, lower = peak, upper = bounds[2], ...)$value
            return(int_left + int_right)
          } else {
            return(stats::integrate(safe_integrand, lower = bounds[1], upper = bounds[2], ...)$value)
          }
        },
        error = function(e) {
          warning(paste("Integration failed:", e$message))
          return(NaN)
        }
      )
    }
  } else {
    FUN <- function(theta, mu, p) {
      series(
        f = \(t) f(t, mu, p) * pdf(t, theta, log = FALSE),
        start = bounds[1],
        end = bounds[2],
        ...
      )
    }
  }

  mapply(
    FUN = FUN,
    theta = theta_list,
    mu = mu,
    p = p
  ) |> unname()
}




#' Calculate Mean for Distribution Objects
#'
#' Computes the expected value (mean) of a distribution object.
#'
#' @description
#' If the distribution object (`x`) contains an analytical formula for the mean,
#' it is used directly. Otherwise (or if `use_moment = TRUE`), the mean is calculated
#' numerically using integration (for continuous) or summation (for discrete) via
#' the \code{\link{moment}} function.
#'
#' @param x An object of class \code{"distrib"}.
#' @param theta A named list of parameters. Vectors are supported.
#' @param use_moment Logical. If \code{TRUE}, forces numerical calculation even if
#'   an analytical formula exists. Defaults to \code{FALSE}.
#' @param ... Additional arguments passed to \code{\link{moment}} (e.g., `subdivisions`, `step`, `tol`).
#'
#' @return A numeric vector representing the mean.
#' @export
mean.distrib <- function(x, theta, use_moment = FALSE, ...) {
  if (use_moment) {
    res <- moment(
      x,
      theta,
      p = 1,
      central = FALSE
    )
  } else {
    res <- x$mean(theta)
    if (is.null(res)) {
      res <- moment(
        x,
        theta,
        p = 1,
        central = FALSE
      )
    }
  }
  res
}




#' Calculate Variance for Distribution Objects
#'
#' Computes the variance (second central moment) of a distribution object.
#'
#' @description
#' Calculates the variance analytically if the formula is available in the object.
#' If not (or if forced), it calculates the variance numerically.
#'
#' @param x An object of class \code{"distrib"}.
#' @param theta A named list of parameters. Vectors are supported.
#' @param use_moment Logical. If \code{TRUE}, forces numerical calculation via the
#'   \code{\link{moment}} function, even if an analytical formula exists.
#' @param ... Additional arguments passed to \code{\link{moment}}.
#'
#' @return A numeric vector representing the variance.
#' @seealso \code{\link{std_dev.distrib}}
#' @export
variance.distrib <- function(x, theta, use_moment = FALSE, ...) {
  if (use_moment) {
    res <- moment(
      x,
      theta,
      p = 2,
      central = TRUE,
      mu = mean.distrib(x, theta, use_moment)
    )
  } else {
    res <- x$variance(theta)
    if (is.null(res)) {
      res <- moment(
        x,
        theta,
        p = 2,
        central = TRUE,
        mu = mean.distrib(x, theta, use_moment)
      )
    }
  }
  res
}




#' Calculate Standard Deviation for Distribution Objects
#'
#' Computes the standard deviation of a distribution object.
#'
#' @description
#' The standard deviation is calculated as the square root of the variance.
#' See \code{\link{variance.distrib}} for details on the calculation strategy.
#'
#' @param x An object of class \code{"distrib"}.
#' @param theta A named list of parameters.
#' @param use_moment Logical. If \code{TRUE}, forces numerical calculation via the
#'   \code{\link{moment}} function, even if an analytical formula exists.
#' @param ... Additional arguments passed to \code{\link{variance.distrib}}.
#'
#' @return A numeric vector representing the standard deviation.
#' @export
std_dev.distrib <- function(x, theta, use_moment = FALSE, ...) {
  sqrt(variance.distrib(x, theta, use_moment))
}




#' Calculate Skewness for Distribution Objects
#'
#' Computes the skewness (third standardized moment) of a distribution object.
#'
#' @description
#' Measures the asymmetry of the distribution. It uses the analytical formula if available.
#' Otherwise, it computes the third central moment via the \code{\link{moment}} function
#' and divides it by the cube of the standard deviation.
#'
#' @param x An object of class \code{"distrib"}.
#' @param theta A named list of parameters.
#' @param use_moment Logical. If \code{TRUE}, forces numerical calculation via the
#'   \code{\link{moment}} function, even if an analytical formula exists.
#' @param ... Additional arguments passed to \code{\link{moment}}.
#'
#' @return A numeric vector representing the kurtosis.
#' @export
skewness.distrib <- function(x, theta, use_moment = FALSE, ...) {
  if (use_moment) {
    m3 <- moment(
      x,
      theta,
      p = 3,
      central = TRUE,
      mu = mean.distrib(x, theta, use_moment)
    )
    s <- std_dev.distrib(x, theta, use_moment)
    res <- m3 / s^3
  } else {
    res <- x$skewness(theta)
    if (is.null(res)) {
      m3 <- moment(
        x,
        theta,
        p = 3,
        central = TRUE,
        mu = mean.distrib(x, theta, use_moment)
      )
      s <- std_dev.distrib(x, theta, use_moment)
      res <- m3 / s^3
    }
  }
  res
}




#' Calculate Kurtosis for Distribution Objects
#'
#' Computes the kurtosis (fourth standardized moment) of a distribution object.
#'
#' @description
#' Measures the "tailedness" of the distribution. It uses the analytical formula if available.
#' Otherwise, it computes the fourth central moment via the \code{\link{moment}} function
#' and divides it by the fourth power of the standard deviation.
#'
#' @param x An object of class \code{"distrib"}.
#' @param theta A named list of parameters.
#' @param use_moment Logical. If \code{TRUE}, forces numerical calculation via the
#'   \code{\link{moment}} function. Defaults to \code{FALSE}.
#' @param ... Additional arguments passed to \code{\link{moment}}.
#'
#' @return A numeric vector representing the kurtosis.
#' @export
kurtosis.distrib <- function(x, theta, use_moment = FALSE, ...) {
  if (use_moment) {
    m4 <- moment(
      x,
      theta,
      p = 4,
      central = TRUE,
      mu = mean.distrib(x, theta, use_moment)
    )
    s <- std_dev.distrib(x, theta, use_moment)
    res <- m4 / s^4 - 3
  } else {
    res <- x$kurtosis(theta)
    if (is.null(res)) {
      m4 <- moment(
        x,
        theta,
        p = 4,
        central = TRUE,
        mu = mean.distrib(x, theta, use_moment)
      )
      s <- std_dev.distrib(x, theta, use_moment)
      res <- m4 / s^4 - 3
    }
  }
  res
}




#' Cumulative Distribution Function for 'distrib' Objects
#'
#' @description
#' Computes the Cumulative Distribution Function (CDF) for a custom distribution
#' object of class \code{distrib}. This method relies on numerical integration (for continuous
#' distributions) or summation (for discrete distributions) of the distribution's kernel.
#'
#' @param x An object of class \code{"distrib"} specifying the distribution components
#'   (type, kernel, normalization constant, and bounds).
#' @param q A numeric vector of quantiles.
#' @param theta A named list of parameters for the distribution.
#'   The names must match the parameters defined in \code{x$params}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(X \leq x)},
#'   otherwise, \eqn{P(X > x)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
#' @param ... Additional arguments passed to the underlying numerical integration function
#'   (\code{\link[stats]{integrate}}) for continuous distributions and numerical
#'   summation function (\code{\link[distrib]{series}}) for discrete distributions..
#'
#' @details
#' The function computes the CDF based on the \code{type} of the distribution defined in \code{x}:
#'
#' \strong{Continuous Distributions:}
#' The probability is calculated via numerical integration of the kernel function divided by the
#' normalization constant:
#' \deqn{F(q) = \dfrac{1}{C(\theta)} \int_{lower}^{q} k(t, \theta) \, dt}
#'
#' \strong{Discrete Distributions:}
#' The probability is calculated via summation. The quantile \code{q} is first floored to the nearest integer
#' (\code{floor(q)}). The function sums the kernel from the lower bound up to \code{q}:
#' \deqn{F(q) = \dfrac{1}{C(\theta)} \sum_{x=lower}^{\lfloor q \rfloor} k(x, \theta)}
#'
#' \strong{Boundary Handling:}
#' \itemize{
#'   \item If \code{q} is smaller than the lower bound, the result is 0.
#'   \item If \code{q} is larger than the upper bound, the result is 1.
#' }
#'
#' \strong{Numerical Stability:}
#' The resulting probabilities are clipped to the range \eqn{[0, 1]} to prevent floating-point
#' errors from producing invalid probabilities (e.g., slightly greater than 1 or less than 0).
#'
#' @return A numeric vector of probabilities (or log-probabilities if \code{log.p = TRUE}).
#'
#' @seealso \code{\link{cdf}}, \code{\link[stats]{integrate}}
#'
#' @export
cdf.distrib <- function(x, q, theta, lower.tail = TRUE, log.p = FALSE, ...) {
  check_params_dim(theta)

  check_params_dim(theta)
  theta <- expand_params(theta)
  theta_list <- transpose_params(theta)

  type <- x$type
  kernel <- x$kernel
  normalization_constant <- x$normalization_constant
  bounds <- x$bounds

  if (type == "continuous") {
    p <- mapply(
      FUN = function(q, theta) {
        if (q < bounds[1]) {
          0
        } else if (q > bounds[2]) {
          1
        } else {
          stats::integrate(
            f = \(t) kernel(t, theta, log = FALSE),
            lower = bounds[1],
            upper = q,
            ...
          )$value
        }
      },
      q = q,
      theta = theta_list
    ) / normalization_constant(theta, log = FALSE)
  } else if (type == "discrete") {
    q <- floor(q)

    p <- mapply(
      FUN = function(q, theta) {
        if (q < bounds[1]) {
          0
        } else if (q > bounds[2]) {
          1
        } else {
          series(
            f = \(t) kernel(t, theta, log = FALSE),
            start = bounds[1],
            end = q
          )
        }
      },
      q = q,
      theta = theta_list
    ) / normalization_constant(theta, log = FALSE)
  }

  p <- pmin(pmax(p, 0), 1)

  if (!lower.tail) {
    p <- 1 - p
  }

  if (log.p) {
    p <- log(p)
  }

  p
}




#' Quantile Function for 'distrib' Objects
#'
#' @description
#' Computes the Quantile Function (Inverse CDF) for a custom distribution object.
#' Since most custom distributions do not have an analytical quantile function,
#' this function employs numerical methods to invert the Cumulative Distribution Function (CDF).
#'
#' @param x An object of class \code{"distrib"}.
#' @param p A vector of probabilities.
#' @param theta A named list of parameters. Vectors are supported (batch calculation).
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(X \le x)}, otherwise \eqn{P(X > x)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
#' @param ... Additional arguments passed to the internal CDF numerical engines.
#'
#' @details
#' **Numerical Strategy:**
#' * **Continuous Distributions:** Uses \code{\link[stats]{uniroot}} to find the root of \eqn{CDF(q) - p = 0}.
#'   It automatically expands the search interval if the root is not bracketed by the initial guess.
#' * **Discrete Distributions:**
#'   1. **Geometric Expansion:** Finds an interval \eqn{[low, high]} containing the quantile by exponentially
#'      expanding a search window around a starting guess (Mean or Mode).
#'   2. **Bisection:** Performs a binary search within that interval to find the smallest integer \eqn{k}
#'      such that \eqn{CDF(k) \ge p}.
#'
#' @return A numeric vector of quantiles.
#' @importFrom stats quantile
#' @export
quantile.distrib <- function(x, p, theta, lower.tail = TRUE, log.p = FALSE, ...) {
  get_start_value <- function(x, theta) {
    # Helper to check if a value is valid (not null, NA, or empty)
    is_valid <- function(val) {
      !is.null(val) && length(val) > 0 && !is.na(val)
    }

    # Priority 1: Analytical Mean
    start <- tryCatch(x$mean(theta), error = function(e) NULL)
    if (is_valid(start)) {
      return(start)
    }

    # Priority 2: Analytical Mode
    start <- tryCatch(x$mode(theta), error = function(e) NULL)
    if (is_valid(start)) {
      return(start)
    }

    # Priority 3: Analytical Median
    start <- tryCatch(x$median(theta), error = function(e) NULL)
    if (is_valid(start)) {
      return(start)
    }

    # Priority 4: Domain Center (Fallback)
    lb <- x$bounds[1]
    ub <- x$bounds[2]

    if (is.finite(lb) && is.finite(ub)) {
      return((lb + ub) / 2)
    } else if (is.finite(lb)) {
      return(lb + 1)
    } else if (is.finite(ub)) {
      return(ub - 1)
    }
    return(0) # Default for (-Inf, Inf)
  }

  quantile_discrete <- function(x, target_p, theta) {
    current <- round(get_start_value(x, theta))
    lb_global <- x$bounds[1]
    ub_global <- x$bounds[2]

    # Calculate CDF at initial guess
    F_curr <- x$cdf(current, theta)

    # Phase A: Bracketing (Geometric Expansion)
    # Find [low, high] such that CDF(low) < p <= CDF(high)
    step <- 1

    if (F_curr < target_p) {
      # Target is to the right
      low <- current
      high <- current + step

      while (x$cdf(high, theta) < target_p) {
        low <- high
        step <- step * 2 # Geometric expansion for speed
        high <- high + step
        if (high >= ub_global) {
          high <- ub_global
          break
        }
      }
    } else {
      # Target is to the left (or exactly here)
      high <- current
      low <- current - step

      while (x$cdf(low, theta) >= target_p) {
        high <- low
        step <- step * 2
        low <- low - step
        if (low <= lb_global) {
          low <- lb_global
          break
        }
      }
    }

    # Phase B: Bisection (Binary Search)
    # Find smallest integer k such that CDF(k) >= target_p
    ans <- high
    while (low <= high) {
      mid <- floor((low + high) / 2)

      # Boundary check safety
      if (mid < lb_global) mid <- lb_global
      if (mid > ub_global) mid <- ub_global

      val <- x$cdf(mid, theta)

      if (val >= target_p) {
        ans <- mid
        high <- mid - 1
      } else {
        low <- mid + 1
      }
    }
    ans
  }

  quantile_continuous <- function(x, target_p, theta) {
    # Objective function: CDF(q) - p = 0
    obj <- function(q) {
      x$cdf(q, theta) - target_p
    }

    start <- get_start_value(x, theta)
    lb <- max(start - 1, x$bounds[1])
    hb <- min(start + 1, x$bounds[2])

    # Expand bracket until signs differ (Root finding prerequisite)
    iter <- 0
    max_iter <- 50

    f_lb <- obj(lb)
    f_hb <- obj(hb)

    while (f_lb * f_hb > 0 && iter < max_iter) {
      step <- 2^iter

      # Expand lower bound
      if (x$bounds[1] == -Inf) {
        lb <- lb - step
      } else {
        lb <- max(lb - step, x$bounds[1])
      }
      # Expand upper bound
      if (x$bounds[2] == Inf) {
        hb <- hb + step
      } else {
        hb <- min(hb + step, x$bounds[2])
      }

      f_lb <- obj(lb)
      f_hb <- obj(hb)
      iter <- iter + 1
    }


    stats::uniroot(obj, interval = c(lb, hb), extendInt = "yes")$root
  }


  # Handle log probabilities and tail
  if (log.p) {
    p <- exp(p)
  }
  if (!lower.tail) {
    p <- 1 - p
  }

  # Validate bounds (0 to 1)
  if (any(p < 0 | p > 1)) {
    stop("Probabilities must be between 0 and 1.")
  }

  # Parameter validation and expansion
  check_params_dim(theta)

  dim_p <- length(p)
  max_dim_theta <- max(lengths(theta))

  if (dim_p != 1 && max_dim_theta != 1 && dim_p != max_dim_theta) {
    stop("Dimension of 'p' does not match dimension of 'theta'.")
  }

  theta <- expand_params(theta)
  theta_list <- transpose_params(theta)

  # Select solver based on distribution type
  solver <- if (x$type == "continuous") {
    quantile_continuous
  } else {
    quantile_discrete
  }

  # Execute vectorized calculation
  mapply(
    FUN = function(p_val, theta_val) solver(x, p_val, theta_val),
    p_val = p,
    theta_val = theta_list
  )
}




#' Generic Random Number Generator
#'
#' @description
#' Generates random samples for a distribution object. This function acts as a generic dispatcher:
#' \itemize{
#'   \item For \strong{discrete} distributions, it uses \strong{Inverse Transform Sampling} via the object's \code{quantile} function.
#'   \item For \strong{continuous} distributions, it delegates to the \strong{Generalized Ratio-of-Uniforms (GRoU)} method
#'   via \code{\link{rng_continuous_grou}}.
#' }
#'
#' @param x A distribution object of class \code{distrib}.
#' @param n Integer. The number of observations to be generated.
#'   \strong{Constraint:} If \code{theta} contains vectorized parameters (lengths > 1),
#'   \code{n} must strictly match the length of the expanded parameter vectors.
#' @param theta A list containing the parameters of the distribution.
#'   \itemize{
#'     \item \strong{Scalar:} All parameters have length 1.
#'     \item \strong{Vectorized:} Parameters can be vectors of length \code{n} or length 1 (which are recycled to length \code{n}).
#'   }
#' @param ... Additional arguments passed to the underlying generation method
#'   (e.g., tuning parameter \code{r} for \code{rng_continuous_grou}).
#'
#' @details
#' \strong{Vectorization and Consistency Strategy:}
#'
#' \itemize{
#'   \item \strong{i.i.d. Sampling (Scalar Parameters):} If all components of \code{theta} have length 1,
#'   the function generates \code{n} independent and identically distributed samples from the same distribution.
#'
#'   \item \strong{Pointwise Sampling (Vector Parameters):} If any component of \code{theta} has length > 1:
#'   \enumerate{
#'     \item The parameters are expanded and recycled to a common length (using standard R recycling rules).
#'     \item \strong{Consistency Check:} The function strictly verifies that the length of the expanded parameters
#'     equals \code{n}. If they do not match, an error is raised to prevent ambiguity.
#'     \item The function generates \strong{one} observation for each set of parameters (i.e., \eqn{y_i \sim D(\theta_i)}).
#'   }
#' }
#'
#' @return A numeric vector of length \code{n} containing the generated random deviates.
#'
#' @importFrom stats runif
#' @seealso \code{\link{rng_continuous_grou}}
#'
#' @export
rng.distrib <- function(x, n, theta, ...) {
  param_lens <- lengths(theta)

  if (!all(param_lens %in% c(1, n))) {
    dims_str <- paste(names(param_lens), param_lens, sep = ": ", collapse = ", ")
    stop(
      sprintf(
        "Parameter dimension mismatch: 'n' is %d, but parameters have dimensions [%s]. All parameters must have length 1 or %d.",
        n, dims_str, n
      )
    )
  }


  if (x$type == "discrete") {
    x$quantile(runif(n), theta)
  } else if (x$type == "continuous") {
    if (any(lengths(theta) > 1)) {
      theta <- transpose_params(expand_params(theta))
      sapply(
        theta,
        \(theta) {
          rng_continuous_grou(1, x, theta, ...)
        }
      )
    } else {
      rng_continuous_grou(n, x, theta, ...)
    }
  }
}




#' Calculate the Expected Value of a Function
#'
#' Computes the expected value of a given function \eqn{f(y)} with respect to a probability distribution defined by \code{distrib}.
#' It automatically handles continuous distributions (via numerical integration) and discrete distributions (via series summation).
#'
#' @param distrib An object of class \code{"distrib"}
#' @param f A function representing the transformation of the random variable \eqn{y}.
#'   **Signature:** It must accept arguments \code{y}, \code{theta}, and \code{...} (see Details).
#' @param theta A named list of parameters for the distribution (e.g., \code{list(mu=10, sigma=2)}).
#'   Vectors inside this list allow computing expectations for multiple distribution parametrizations at once.
#' @param ... Additional arguments passed directly to the function \code{f}.
#'   **Vectorization:** These arguments are fully vectorized. If vectors are provided, they are recycled
#'   against the parameters in \code{theta} according to standard R recycling rules.
#'
#' @details
#' The function calculates:
#' \itemize{
#'   \item \eqn{E[f(Y)] = \int_{lb}^{ub} f(y, \theta, \dots) \cdot p(y|\theta) \, dy} (Continuous)
#'   \item \eqn{E[f(Y)] = \sum_{y=lb}^{ub} f(y, \theta, \dots) \cdot P(y|\theta)} (Discrete)
#' }
#'
#' **Vectorization:**
#' The function iterates over the longest vector found among \code{theta} and \code{...}.
#' For example, if \code{theta$mu} has length 2 and you pass a vector of length 2 to \code{...},
#' the function computes the expectation for the paired values. If lengths differ, standard R recycling applies.
#'
#' **Requirements for `f`:**
#' The user-provided function \code{f} must be defined with the signature:
#' \code{f(y, theta, ...)}
#'
#' @return A numeric vector containing the expected values. The length corresponds to the
#'   maximum length among all vectors in \code{theta} and \code{...}.
#'
#' @importFrom stats integrate
#'
#' @examples
#' \dontrun{
#' distrib <- negbin_distrib()
#'
#' # Define f accepting y, theta, and extra parameter gamma
#' f_pow <- function(y, theta, gamma = 1) {
#'   y^gamma
#' }
#'
#' # --- Example 1: Basic usage ---
#' # Calculate E[y^2] for mu=10
#' expectation(distrib, f_pow, theta = list(mu = 10, theta = 1), gamma = 2)
#'
#' # --- Example 2: Vectorization over '...' ---
#' # Calculate 1st, 2nd, and 3rd raw moments (E[y^1], E[y^2], E[y^3]) simultaneously
#' # mu is fixed at 10, gamma varies (1, 2, 3)
#' expectation(distrib, f_pow, theta = list(mu = 10, theta = 1), gamma = c(1, 2, 3))
#'
#' # --- Example 3: Joint Vectorization ---
#' # Calculate E[y^1] for mu=10 and E[y^2] for mu=20
#' expectation(distrib, f_pow,
#'   theta = list(mu = c(10, 20), theta = 1),
#'   gamma = c(1, 2)
#' )
#' }
#'
#' @export
expectation <- function(distrib, f, theta, ...) {
  # Capture extra arguments and check for name collisions
  dots <- list(...)
  if (any(names(dots) %in% names(theta))) {
    stop("Arguments in '...' cannot have the same names as parameters in 'theta'.")
  }

  # Combine all parameters to handle vectorization
  all_params <- c(theta, dots)
  n_theta <- distrib$n_params

  # Define the worker function for a single set of parameters
  compute_single <- function(params) {
    # Extract specific subsets
    p_theta <- params[1:n_theta] # Strictly for the distribution

    p_dots <- params[-(1:n_theta)] # Strictly for the function f

    # Internal function to integrate/sum
    integrand <- function(y) {
      # Call f(y, theta, ...) dynamically
      val_f <- do.call(f, c(list(y = y, theta = p_theta), p_dots))

      # Call pdf(y, theta) strictly
      val_p <- distrib$pdf(y, p_theta, log = FALSE)

      val_f * val_p
    }

    # Execute based on distribution type
    if (distrib$type == "continuous") {
      integrate(integrand, lower = distrib$bounds[1], upper = distrib$bounds[2])$value
    } else {
      series(integrand, start = distrib$bounds[1], end = distrib$bounds[2])
    }
  }

  sapply(transpose_params(expand_params(all_params)), compute_single)
}




#' Monte Carlo Approximation of the Expected Hessian
#'
#' @description
#' Approximates the expected Hessian matrix (Fisher Information Matrix elements) for a given
#' distribution object using Monte Carlo simulation.
#'
#' This function relies on the identity that, under regularity conditions, the expected Hessian
#' is equal to the negative of the expected outer product of the score function (gradients):
#' \deqn{\mathbb{E}\left[\nabla^2 \ell(\theta)\right] = -\mathbb{E}\left[\nabla \ell(\theta) (\nabla \ell(\theta))^\top\right]}
#'
#' @param distrib An object of class \code{"distrib"}.
#'   Must contain \code{params}, \code{gradient}, and \code{rng} components.
#' @param y A vector or matrix of the observed response variable. Used to determine the number of
#'   observations (\code{n}).
#' @param theta A list containing the parameter values.
#' @param nsim Integer. The number of Monte Carlo simulations to perform for each observation.
#'   Defaults to 1000. Higher values provide better approximations but increase computation time.
#'
#' @details
#' For each observation \eqn{k} in \code{y}, the function:
#' \enumerate{
#'   \item Simulates \code{nsim} random values from the distribution using the parameters \eqn{\theta_k}.
#'   \item Calculates the gradient (score vector) for these simulated values.
#'   \item Computes the empirical mean of the negative cross-products of the gradients.
#' }
#'
#' The resulting approximation is useful when the analytical expected Hessian is difficult
#' or impossible to derive in closed form.
#'
#' @return A named list of length \eqn{p(p+1)/2} (where \eqn{p} is the number of parameters),
#' containing the unique elements of the expected Hessian matrix.
#'
#' The list elements are named as follows:
#' \itemize{
#'   \item Diagonal elements: \code{"par_par"} (e.g., \code{"mu_mu"}, \code{"sigma_sigma"})
#'   \item Off-diagonal elements: \code{"par1_par2"} (e.g., \code{"mu_sigma"}), stored in the order
#'     defined by the distribution object.
#' }
#' Each element of the list is a numeric vector of length \code{NROW(y)}, representing the
#' expected second derivative for each observation.
#'
#'
#' @export
mc_expected_hessian <- function(distrib, y, theta, nsim = 1000) {
  params <- distrib$params
  n_params <- distrib$n_params
  gradient <- distrib$gradient
  theta_k <- transpose_params(expand_params(theta))
  dim_hess <- .5 * n_params * (n_params + 1)
  out <- vector("list", dim_hess)
  n <- length(theta_k)

  hess_names <- character(dim_hess)
  for (i in 1:n_params) {
    hess_names[i] <- paste0(params[i], "_", params[i])
  }
  idx <- n_params + 1
  for (i in 1:(n_params - 1)) {
    for (j in (i + 1):n_params) {
      hess_names[idx] <- paste0(params[i], "_", params[j])
      idx <- idx + 1
    }
  }

  names(out) <- hess_names

  for (i in 1:dim_hess) {
    out[[i]] <- numeric(n)
  }

  for (k in 1:n) {
    g <- distrib$gradient(
      distrib$rng(nsim, theta_k[[k]]),
      theta_k[[k]]
    )
    k_mixed <- n_params + 1
    for (i in 1:n_params) {
      for (j in i:n_params) {
        if (i == j) {
          idx <- i
        } else {
          idx <- k_mixed
          k_mixed <- k_mixed + 1
        }
        out[[idx]][k] <- mean(-g[[i]] * g[[j]])
      }
    }
  }
  expand_params(out, NROW(y))
}




#' Empirical Variance (Population)
#'
#' @description
#' Calculates the empirical variance (second central moment) of a numeric vector.
#'
#' @details
#' Unlike \code{\link[stats]{var}}, which calculates the unbiased sample variance
#' (dividing by \eqn{n-1}), this function calculates the population variance
#' (dividing by \eqn{n}).
#'
#' The formula used is:
#' \deqn{s^2 = \frac{1}{n} \sum_{i=1}^{n} (x_i - \bar{x})^2}
#'
#' @param x A numeric vector.
#' @param na.rm Logical. Should missing values be removed before calculation? Defaults to \code{FALSE}.
#' @param ... Additional arguments (not used).
#'
#' @return A numeric scalar representing the population variance.
#' @seealso \code{\link[stats]{var}} for the sample variance.
#' @export
variance.default <- function(x, na.rm = FALSE, ...) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  sum((x - mean(x))^2) / NROW(x)
}

#' Empirical Standard Deviation (Population)
#'
#' @description
#' Calculates the empirical standard deviation of a numeric vector.
#'
#' @details
#' This is calculated as the square root of the population variance:
#' \deqn{s = \sqrt{\frac{1}{n} \sum_{i=1}^{n} (x_i - \bar{x})^2}}
#'
#' Note that this differs from \code{\link[stats]{sd}}, which uses \eqn{n-1} in the denominator.
#'
#' @param x A numeric vector.
#' @param na.rm Logical. Should missing values be removed? Defaults to \code{FALSE}.
#' @param ... Additional arguments (not used).
#'
#' @return A numeric scalar representing the population standard deviation.
#' @export
std_dev.default <- function(x, na.rm = FALSE, ...) {
  sqrt(variance(x, na.rm = na.rm))
}

#' Empirical Skewness
#'
#' @description
#' Calculates the coefficient of skewness (third standardized moment) for a numeric vector.
#'
#' @details
#' Computes the Fisher-Pearson coefficient of skewness:
#' \deqn{g_1 = \frac{\frac{1}{n} \sum_{i=1}^n (x_i - \bar{x})^3}{ \left( \frac{1}{n} \sum_{i=1}^n (x_i - \bar{x})^2 \right)^{3/2} }}
#'
#' @param x A numeric vector.
#' @param na.rm Logical. Should missing values be removed? Defaults to \code{FALSE}.
#' @param ... Additional arguments (not used).
#'
#' @return A numeric scalar. Negative values indicate a left tail, positive values indicate a right tail.
#' @export
skewness.default <- function(x, na.rm = FALSE, ...) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  (sum((x - mean(x))^3) / NROW(x)) / ((std_dev(x))^3)
}

#' Empirical Excess Kurtosis
#'
#' @description
#' Calculates the coefficient of **excess** kurtosis (fourth standardized moment minus 3) for a numeric vector.
#'
#' @details
#' Computes the sample excess kurtosis:
#' \deqn{g_2 = \frac{\frac{1}{n} \sum_{i=1}^n (x_i - \bar{x})^4}{ \left( \frac{1}{n} \sum_{i=1}^n (x_i - \bar{x})^2 \right)^2 } - 3}
#'
#' A value of 0 implies the distribution has the same "tailedness" as a Normal distribution.
#' Positive values indicate heavier tails (leptokurtic), negative values indicate lighter tails (platykurtic).
#'
#' @param x A numeric vector.
#' @param na.rm Logical. Should missing values be removed? Defaults to \code{FALSE}.
#' @param ... Additional arguments (not used).
#'
#' @return A numeric scalar representing the excess kurtosis.
#' @export
kurtosis.default <- function(x, na.rm = FALSE, ...) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  (sum((x - mean(x))^4) / NROW(x)) / ((std_dev(x))^4) - 3
}





#' Numerical Mode Finding for `distrib` Objects
#'
#' Calculates the mode (the value \eqn{y} that maximizes the probability density or mass function)
#' for a given distribution object. The function is generic, robust, and does not require
#' analytical derivatives or prior knowledge of the distribution's moments (e.g., mean).
#'
#' @description
#' This function employs a two-stage numerical optimization strategy to find the global maximum
#' of the distribution's kernel. It automatically detects the distribution type (continuous or discrete)
#' and applies the appropriate algorithm:
#' * **Continuous:** Exponential Search for bracketing followed by Brent's Method.
#' * **Discrete:** Exponential Search for bracketing followed by Integer Binary Search.
#'
#' @param x An object of class \code{"distrib"} containing the distribution definition
#'   (kernel, bounds, type, etc.).
#' @param theta A named list of parameters for the distribution. Vectors are supported for batch calculation;
#'   the function automatically expands and transposes parameters to handle multiple combinations.
#' @param maxit Integer. The maximum number of iterations allowed for the "hunting" phase
#'   (Exponential Search) to bracket the peak. Defaults to 10000.
#' @param tol Numeric. The convergence tolerance for the refinement phase (used only for
#'   continuous distributions). Defaults to \code{1e-10}.
#'
#' @details
#' **Algorithm Strategy:**
#' Since the distribution shape is treated as a "black box" (unimodal but potentially skewed or heavy-tailed),
#' standard gradient methods are avoided in favor of direct search methods that are robust to flat regions
#' and numerical noise.
#'
#' \strong{Continuous Distributions:}
#' \enumerate{
#'   \item \emph{Initialization:} Determines a safe starting point based on bounds.
#'   \item \emph{Gradient Check:} Probes the local density to determine the uphill direction.
#'   \item \emph{Exponential Search (The Hunt):} Steps away from the start point with exponentially increasing
#'     step sizes (\eqn{1, 2, 4, \dots}) until the density begins to drop. This guarantees bracketing the mode
#'     in \eqn{O(\log N)} steps, where \eqn{N} is the distance from the start.
#'   \item \emph{Refinement (The Zoom):} Uses Brent's Method (via \code{\link[stats]{optimize}}) within the
#'     found bracket to pinpoint the mode with precision \code{tol}.
#' }
#'
#' \strong{Discrete Distributions:}
#' \enumerate{
#'   \item \emph{Bracketing:} Similar to the continuous case, uses Exponential Search to find an integer interval
#'     containing the peak.
#'   \item \emph{Integer Binary Search:} Performs a binary search on the discrete domain. This reduces the search
#'     space logarithmically.
#'   \item \emph{Tie-Breaking:} In cases of bimodality (e.g., Poisson with integer \eqn{\lambda} where \eqn{P(k) = P(k-1)}),
#'     the algorithm adopts a "Greedy Right" approach, returning the larger integer (consistent with \code{floor} behavior
#'     for standard parameters).
#' }
#'
#' **Complexity:**
#' For both cases, if the mode is located at distance \eqn{N} from the origin (or lower bound), the time complexity is
#' \eqn{O(\log N)} in terms of kernel evaluations. This makes the function extremely efficient even for
#' distributions centered far from the origin (e.g., \eqn{\mu = 10^6}).
#'
#' @return A numeric vector containing the calculated mode(s). The length corresponds to the number of
#' parameter sets provided in `theta`.
#'
#' @importFrom stats optimize
#' @export
mode.distrib <- function(x, theta, maxit = 10000, tol = 1e-10) {
  mode_continuous <- function(x, theta, maxit = 10000, tol = 1e-10) {
    lb <- x$bounds[1]
    ub <- x$bounds[2]

    # Wrapper for the kernel function handling bounds and numerical stability
    obj_fun <- function(t) {
      if (t <= lb || t >= ub) {
        return(-Inf)
      }
      val <- x$kernel(t, theta, log = TRUE)
      if (is.na(val) || is.nan(val)) {
        return(-Inf)
      }
      val
    }

    # CASE A: Finite Bounds
    # If both bounds are finite, we can skip the hunt and use optimize directly
    if (is.finite(lb) && is.finite(ub)) {
      optimize(obj_fun, interval = c(lb, ub), maximum = TRUE, tol = tol)$maximum
    }

    # CASE B: Infinite Bounds (requires bracketing)

    # 1. Determine a safe starting point
    current_x <- 0
    if (is.finite(lb)) {
      current_x <- lb + 1 # Start slightly inside lower bound
    } else if (is.finite(ub)) {
      current_x <- ub - 1 # Start slightly inside upper bound
    }

    # 2. Determine initial direction (Local Gradient)
    f_curr <- obj_fun(current_x)
    step <- 1

    # Probe slightly to the right to see if density increases or decreases
    f_right <- obj_fun(current_x + 1e-04)

    if (f_right > f_curr) {
      direction <- 1 # Mode is likely to the right
    } else {
      direction <- -1 # Mode is likely to the left
    }

    # 3. Exponential Search (The Hunt)
    # Expand step size exponentially until the density value drops.
    # This brackets the mode between the previous point and the current point.

    prev_x <- current_x
    iter <- 0

    while (iter < maxit) {
      next_x <- current_x + (direction * step)

      # Check physical bounds (if one side is finite)
      if (next_x <= lb || next_x >= ub) {
        # If we hit a boundary, truncate search there and use it for optimization
        next_x <- if (direction == 1) min(next_x, ub) else max(next_x, lb)
        break
      }

      f_next <- obj_fun(next_x)

      # STOP CONDITION: Did we cross the peak?
      if (f_next < f_curr) {
        # Yes, the value dropped. The mode is trapped between prev_x and next_x.
        break
      }

      # If not, prepare for the next step
      step <- step * 2 # Double the step size (Exponential Search)
      prev_x <- current_x # Save the last valid "uphill" point
      current_x <- next_x
      f_curr <- f_next
      iter <- iter + 1
    }

    # 4. Refinement (Brent's Method)
    # We now have the mode bracketed in [min(prev, next), max(prev, next)]
    search_range <- sort(c(prev_x, next_x)) + c(-1e-03, 1e-03)

    # Correct for any potential bound oversteps (safety check)
    search_range[1] <- max(search_range[1], lb)
    search_range[2] <- min(search_range[2], ub)

    res <- optimize(obj_fun, interval = search_range, maximum = TRUE, tol = tol)
    res$maximum
  }

  mode_discrete <- function(x, theta, maxit = 10000) {
    # Bounds: ensure they are integers (e.g., Poisson starts at 0, not 0.5)
    lb <- ceiling(x$bounds[1])
    ub <- floor(x$bounds[2])

    # Wrapper for the kernel that handles integers, bounds, and numerical stability
    obj_fun <- function(k) {
      k <- round(k) # Enforce integer evaluation

      # Check physical bounds
      if (k < lb || k > ub) {
        return(-Inf)
      }
      val <- x$kernel(k, theta, log = TRUE)
      # Handle numerical errors (NA/NaN) by treating them as probability 0 (-Inf log)
      if (is.na(val) || is.nan(val)) {
        return(-Inf)
      }
      val
    }

    # --- 1. Determine Safe Starting Point ---
    # Default to 0, but respect bounds if they exist
    current_k <- 0

    if (is.finite(lb)) {
      current_k <- lb
    }
    # If lower bound is infinite but upper is finite, start at upper bound
    if (!is.finite(lb) && is.finite(ub)) {
      current_k <- ub
    }

    # Sanity Check: Ensure the starting point gives a valid density.
    # If it returns -Inf (e.g., due to edge cases), try to nudge it valid.
    if (obj_fun(current_k) == -Inf) {
      # Fallback: move slightly inside the domain if possible
      if (is.finite(lb)) {
        current_k <- lb + 1
      }
    }

    # --- 2. Determine Initial Direction (Gradient) ---
    val_curr <- obj_fun(current_k)
    val_right <- obj_fun(current_k + 1)

    direction <- 0

    if (val_right > val_curr) {
      direction <- 1 # Slope is positive, mode is to the right
    } else {
      # If slope is not positive to the right, check the left to be sure
      val_left <- obj_fun(current_k - 1)

      if (val_left > val_curr) {
        direction <- -1 # Slope is negative, mode is to the left
      } else {
        # If current is >= both left and right, we are already at the mode!
        return(current_k)
      }
    }

    # --- 3. Exponential Search (The Hunt) ---
    # Rapidly expand the search interval until we bracket the peak.
    step <- 1
    prev_k <- current_k
    iter <- 0

    found_bracket <- FALSE
    limit_hit <- FALSE

    while (iter < maxit) {
      next_k <- current_k + (direction * step)

      # Check hard bounds during the hunt
      if (next_k < lb) {
        next_k <- lb
        limit_hit <- TRUE
      }
      if (next_k > ub) {
        next_k <- ub
        limit_hit <- TRUE
      }

      val_next <- obj_fun(next_k)

      # STOP CONDITION: If the value drops, we have crossed the peak.
      if (val_next < val_curr) {
        found_bracket <- TRUE
        break
      }

      if (limit_hit) {
        break # We hit the edge of the domain and are still climbing
      }

      # Update for the next step (Double the step size for exponential reach)
      step <- step * 2
      prev_k <- current_k
      current_k <- next_k
      val_curr <- val_next
      iter <- iter + 1
    }

    # Define the final search bracket
    # The peak is strictly contained between prev_k and next_k
    search_low <- min(prev_k, next_k)
    search_high <- max(prev_k, next_k)

    # --- 4. Integer Binary Search (The Zoom) ---
    # Since the function is discrete and unimodal, we use binary search
    # to find the maximum in O(log N) time.

    while (search_low < search_high) {
      # Calculate integer midpoint
      mid <- floor((search_low + search_high) / 2)

      val_mid <- obj_fun(mid)
      val_mid_plus_1 <- obj_fun(mid + 1)

      if (val_mid_plus_1 >= val_mid - 1e-9) {
        # Slope is positive at 'mid', so the peak is to the right
        search_low <- mid + 1
      } else {
        # Slope is negative or flat, so the peak is at 'mid' or to the left
        search_high <- mid
      }
    }

    search_low
  }


  if (x$type == "continuous") {
    sapply(
      transpose_params(expand_params(theta)),
      \(th) {
        mode_continuous(x, theta, maxit, tol)
      }
    )
  } else if (x$type == "discrete") {
    sapply(
      transpose_params(expand_params(theta)),
      \(th) {
        mode_discrete(x, theta, maxit)
      }
    )
  }
}




#' Generalized Ratio-of-Uniforms (GRoU) Random Number Generator
#'
#' @description
#' Generates random variates from a continuous distribution object using the Generalized
#' Ratio-of-Uniforms (GRoU) method.
#'
#' This function is a generic sampler that does not require the inverse CDF. Instead, it relies on
#' the kernel of the probability density function and the mode of the distribution. It automatically
#' constructs the bounding rectangle for the acceptance-rejection method.
#'
#' @param n Integer. The number of random variates to generate.
#' @param x An object of class \code{"distrib"}.
#'   Must contain the following components: \code{kernel}, \code{mode}, \code{variance}, and \code{bounds}.
#' @param theta A named list of parameters for the distribution.
#' @param r Numeric. The tuning parameter for the transformation power. Defaults to 2.
#'   Common choices are:
#'   \itemize{
#'     \item \code{r = 1}: Standard Ratio-of-Uniforms.
#'     \item \code{r = 2}: Often more efficient for heavy-tailed distributions.
#'   }
#'
#' @details
#' \strong{The GRoU Method:}
#'
#' The Generalized Ratio-of-Uniforms method generates random variables \eqn{Y} with density proportional
#' to a kernel \eqn{K(y)} by simulating a point \eqn{(U, V)} uniformly over the region:
#' \deqn{A_r = \left\{ (u, v) : 0 < u \leq \left[ K\left( \frac{v}{u^r} \right) \right]^{\frac{1}{r+1}} \right\}}
#'
#' If \eqn{(U, V)} is uniformly distributed over \eqn{A_r}, then \eqn{Y = V / U^r} has the desired distribution.
#'
#' \strong{Algorithm Steps:}
#'
#' \enumerate{
#'   \item \strong{Bounding Rectangle Construction:} The algorithm first determines the bounding box
#'   \eqn{[0, u_{\max}] \times [v_{\min}, v_{\max}]} enclosing \eqn{A_r}.
#'   \itemize{
#'      \item \eqn{u_{\max} = \sup_y K(y)^{\frac{1}{r+1}}}. Calculated at the distribution's mode.
#'      \item \eqn{v_{\min} = \inf_y h(y)} and \eqn{v_{\max} = \sup_y h(y)}, where the auxiliary function is
#'      \eqn{h(y) = y \cdot K(y)^{\frac{r}{r+1}}}.
#'   }
#'   \item \strong{Heuristic Optimization:} To find \eqn{v_{\min}} and \eqn{v_{\max}}, the function performs
#'   a dynamic search starting from the mode and expanding outwards using the standard deviation (derived from
#'   \code{x$variance}) as a step size, until the function value decreases sufficiently. This defines the search
#'   interval for the \code{\link[stats]{optimize}} function.
#'   \item \strong{Sampling Loop:}
#'   \itemize{
#'      \item Generate \eqn{U \sim \text{Unif}(0, u_{\max})} and \eqn{V \sim \text{Unif}(v_{\min}, v_{\max})}.
#'      \item Compute candidate \eqn{Y = V / U^r}.
#'      \item Accept \eqn{Y} if \eqn{(r+1) \log(U) \leq \log K(Y)}.
#'   }
#' }
#'
#' \strong{Robustness:}
#' The function includes safeguards against infinite loops (max 10 consecutive failures in batch sampling)
#' and checks for NaN/NA values during the initialization of bounds.
#'
#' @return A numeric vector of length \code{n} containing the generated random variates.
#'
#' @importFrom stats runif optimize
#' @export
rng_continuous_grou <- function(n, x, theta, r = 2) {
  obj_fun <- function(t) {
    if (t < lb || t > ub) {
      return(0)
    }
    val <- t * x$kernel(t, theta, log = FALSE)^r_ratio
    if (is.na(val)) 0 else val
  }
  lb <- x$bounds[1]
  ub <- x$bounds[2]
  mode_val <- x$mode(theta)
  r_ratio <- r / (r + 1)
  u_max <- x$kernel(mode_val, theta, log = FALSE)^(1 / (r + 1))

  if (is.nan(u_max) | is.na(u_max)) {
    stop("NaN or NA value in u_max: the log-kernel is possibly not bounded")
  }

  found_lower <- FALSE
  delta <- tryCatch(sqrt(x$variance(theta)), error = function(e) 1)
  if (is.nan(delta) || !is.numeric(delta)) {
    delta <- 1
  }
  current_t <- mode_val
  f_curr <- obj_fun(current_t)

  while (!found_lower) {
    new_t <- max(lb, current_t - delta)
    f_new <- obj_fun(new_t)
    if (f_new < f_curr) {
      delta <- delta * 2
      f_curr <- f_new
      current_t <- new_t
      if (new_t <= lb) {
        found_lower <- TRUE
        search_limit_left <- lb
      }
    } else {
      search_limit_left <- new_t
      found_lower <- TRUE
    }
  }

  found_upper <- FALSE
  delta <- tryCatch(sqrt(x$variance(theta)), error = function(e) 1)
  if (is.nan(delta) || !is.numeric(delta)) {
    delta <- 1
  }
  current_t <- mode_val
  f_curr <- obj_fun(current_t)

  while (!found_upper) {
    new_t <- min(ub, current_t + delta)
    f_new <- obj_fun(new_t)
    if (f_new > f_curr) {
      delta <- delta * 2
      f_curr <- f_new
      current_t <- new_t
      if (new_t >= ub) {
        found_upper <- TRUE
        search_limit_right <- ub
      }
    } else {
      search_limit_right <- new_t
      found_upper <- TRUE
    }
  }

  v_min <- optimize(obj_fun, lower = search_limit_left, upper = search_limit_right, maximum = FALSE)$objective
  v_max <- optimize(obj_fun, lower = search_limit_left, upper = search_limit_right, maximum = TRUE)$objective


  sampled <- numeric(n)
  n_accepted <- 0
  failed <- 0

  while (n_accepted < n) {
    needed <- n - n_accepted
    batch_size <- max(10, ceiling(needed * 1.5))

    u <- runif(batch_size, 0, u_max)
    v <- runif(batch_size, v_min, v_max)
    y_cand <- v / (u^r)

    valid_range <- (y_cand >= lb & y_cand <= ub)
    if (sum(valid_range) <= 0) {
      failed <- failed + 1
      if (failed >= 10) {
        stop("Could not sample in 10 tentatives")
      }
      next
    } else {
      failed <- 0
      y_cand <- y_cand[valid_range]
      u <- u[valid_range]
      keep <- ((r + 1) * log(u) <= x$kernel(y_cand, theta, log = TRUE))

      new_values <- y_cand[keep]
      n_found <- length(new_values)

      if (n_found > 0) {
        failed <- 0
        take <- min(n_found, needed)
        idx_start <- n_accepted + 1
        idx_end <- n_accepted + take
        sampled[idx_start:idx_end] <- new_values[1:take]
        n_accepted <- n_accepted + take
      } else {
        failed <- failed + 1
        if (failed >= 10) {
          stop("Could not sample in 10 tentatives")
        }
      }
    }
  }
  sampled
}
