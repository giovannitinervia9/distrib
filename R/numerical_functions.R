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
  if (dim_p != 1 && dim_p != max_dim_theta) {
    stop("Provided value of p has a dimension which does not match with theta")
  }

  theta_mat <- do.call(cbind, theta)
  theta_list <- split(theta_mat, row(theta_mat)) |>
    lapply(function(th) {
      th <- as.list(th)
      names(th) <- params
      th
    })

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
#' @method cdf distrib
#' @export
cdf.distrib <- function(x, q, theta, lower.tail = TRUE, log.p = FALSE, ...) {
  check_params_dim(theta)

  params <- x$params
  theta_mat <- do.call(cbind, theta)
  theta_list <- split(theta_mat, row(theta_mat)) |>
    lapply(function(th) {
      th <- as.list(th)
      names(th) <- params
      th
    })

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
