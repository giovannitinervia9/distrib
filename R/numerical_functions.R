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
#' For infinite bounds, the function checks for convergence: the loop terminates when the sum of
#' the current block (size `step`) drops below the tolerance threshold `tol`.
#'
#' @param f A function taking a vector of integers `x` and returning a vector of numeric values.
#'   **Must be vectorized** (able to handle multiple inputs at once).
#' @param start Numeric. The starting value of the sequence. Can be finite, `Inf`, or `-Inf`.
#' @param end Numeric. The ending value of the sequence. Can be finite, `Inf`, or `-Inf`.
#'   Defaults to `Inf`.
#' @param step Integer. The number of terms to calculate in a single vectorized batch.
#'   Larger values generally improve speed but increase memory usage. Defaults to 100.
#' @param tol Numeric. Tolerance threshold for convergence. The series stops when the
#'   sum of the current batch is less than `tol`. Defaults to `1e-10`.
#' @param maxit Integer. Safety limit for the maximum number of batch iterations to prevent
#'   infinite loops in non-convergent series. Defaults to 100,000.
#'
#' @return A numeric scalar representing the calculated sum.
#'
#' @examples
#' # --- Case 1: Standard Convergent Series (Basel Problem) ---
#' # Sum of 1/x^2 from 1 to Inf => pi^2 / 6
#' f_basel <- function(x) 1 / x^2
#' series(f_basel, start = 1, end = Inf)
#'
#' # --- Case 2: Finite Sum ---
#' # Sum from 1 to 5
#' series(f_basel, start = 1, end = 5)
#'
#' # --- Case 3: Backward Series (Negative Domain) ---
#' # Sum of 1/x^2 from -1 down to -Inf
#' series(f_basel, start = -1, end = -Inf)
#'
#' # --- Case 4: Doubly Infinite Series (Discretized Normal) ---
#' # Sum of Normal density over all integers Z (-Inf to +Inf)
#' # Should sum to approximately 1
#' f_norm <- function(x) dnorm(x)
#' series(f_norm, start = -Inf, end = Inf)
#'
#' @export
series <- function(f, start, end = Inf, step = 100, tol = 1e-10, maxit = 100000) {
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

  eps <- 2 * tol
  upper_limit <- min(start_internal + step, end_internal)
  x <- start_internal:upper_limit
  it <- 0
  while (eps > tol) {
    it <- it + 1

    eps <- sum(f_internal(x))
    s <- s + eps

    last_val <- x[length(x)]

    if (last_val >= end_internal) {
      break
    }

    next_start <- last_val + 1
    next_limit <- min(next_start + step, end_internal)
    x <- next_start:next_limit

    if (it >= maxit) {
      warning("Series may not be convergent")
      break
    }
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
#' * For **continuous** distributions (`type = "continuous"`), it uses \code{\link[stats]{integrate}}.
#' * For **discrete** distributions (`type = "discrete"`), it uses the internal function \code{series}.
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
      stats::integrate(
        f = \(t) f(t, mu, p) * pdf(t, theta, log = FALSE),
        lower = bounds[1],
        upper = bounds[2],
        ...
      )$value
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




#' Calculate Kurtosis for Distribution Objects
#'
#' Computes the kurtosis (fourth standardized moment) of a distribution object.
#'
#' @description
#' Measures the "tailedness" of the distribution. It uses the analytical formula if available.
#' Otherwise, it computes the fourth central moment and divides it by the fourth power of the
#' standard deviation.
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
    res <- m4 / s^4
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
      res <- m4 / s^4
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




#' Calculate the Expected Value of a Function
#'
#' Computes the expected value of a given function \eqn{f(y)} with respect to a probability distribution defined by \code{distrib}.
#' It automatically handles continuous distributions (via numerical integration) and discrete distributions (via series summation).
#'
#' @param distrib An object of class \code{"distrib"}
#' @param f A function representing the transformation of the random variable \eqn{y}.
#'   **Signature:** It must accept arguments \code{y}, \code{theta}, and \code{...} (see Details).
#' @param theta A named list of parameters for the distribution.
#' @param ... Additional arguments passed directly to the function \code{f}.
#'   These allow parametrizing the transformation function (e.g., passing an exponent or shift parameter).
#'
#' @details
#' The function calculates:
#' \itemize{
#'   \item \eqn{E[f(Y)] = \int_{lb}^{ub} f(y, \theta, \dots) \cdot p(y|\theta) \, dy} (Continuous)
#'   \item \eqn{E[f(Y)] = \sum_{y=lb}^{ub} f(y, \theta, \dots) \cdot P(y|\theta)} (Discrete)
#' }
#'
#' **Requirements for `f`:**
#' The user-provided function \code{f} must be defined with the signature:
#' \code{f(y, theta, ...)}
#' where:
#' \itemize{
#'   \item \code{y}: The integration/summation variable (numeric vector).
#'   \item \code{theta}: A named list of parameters for the distribution (e.g., \code{list(mu=2, theta=1)}).
#'   \item \code{...}: Optional additional arguments to support custom parameters passed via \code{expectation}.
#' }
#'
#' @return A numeric vector containing the expected values. The length corresponds to the
#'   longest parameter vector in \code{theta}.
#'
#' @importFrom stats integrate
#'
#' @examples
#' \dontrun{
#' # Example 1: Expected value of y^gamma (Raw Moment)
#' distrib <- negbin_distrib()
#' theta <- list(mu = 10, theta = 1)
#'
#' # Define f accepting y, theta, and extra parameter gamma
#' f_pow <- function(y, theta, gamma = 1) {
#'   y^gamma
#' }
#'
#' # Calculate E[y^2]
#' expectation(distrib, f_pow, theta, gamma = 2)
#'
#' # Example 2: Trig transformation
#' f_trig <- function(y, theta, ...) {
#'   sin(y)
#' }
#' expectation(distrib, f_trig, theta)
#' }
#'
#' @export
expectation <- function(distrib, f, theta, ...) {
  type <- distrib$type
  bounds <- distrib$bounds
  pdf <- distrib$pdf
  params <- distrib$params

  max_dim_theta <- max(lengths(theta))
  check_params_dim(theta)

  lb <- distrib$bounds[1]
  ub <- distrib$bounds[2]

  f_internal <- function(y, theta, ...) {
    f(y, theta, ...) * pdf(y, theta, log = FALSE)
  }

  if (type == "continuous") {
    FUN <- function(theta) {
      integrate(\(y) f_internal(y, theta, ...), lower = lb, upper = ub)$value
    }
  } else if (type == "discrete") {
    FUN <- function(theta) {
      series(\(y) f_internal(y, theta, ...), start = lb, end = ub)
    }
  }
  sapply(transpose_params(expand_params(theta)), FUN)
}
