#' Transformer Objects
#'
#' @name transformer
#' @rdname transformer
#'
#' @description
#' A \code{transformer} object is a list of class \code{"transformer"} that defines
#' the mathematical rules for transforming a random variable. These objects are
#' used as input for the \code{\link{transformation}} function.
#'
#' @section Structure:
#' A \code{transformer} object must contain the following components:
#' \itemize{
#'   \item \strong{name}: A character string identifying the transformation (e.g., "log", "logit").
#'   \item \strong{trans_fun}: The forward transformation function \eqn{y = g(x)}.
#'   \item \strong{trans_inv}: The inverse transformation function \eqn{x = g^{-1}(y)}.
#'   \item \strong{trans_abs_jac}: The absolute value of the Jacobian of the inverse transformation \eqn{\left|J(y)\right| = \left|\dfrac{dX}{dY}\right|}.
#'     It must accept a \code{log} argument to return the log-absolute Jacobian.
#'   \item \strong{bounds_fun}: A function that takes the original distribution bounds and
#'     returns the transformed support bounds.
#'   \item \strong{valid_support}: A function that checks if the input distribution's support
#'     is mathematically compatible with the transformation.
#'   \item \strong{decreasing}: A logical value indicating if the transformation is
#'     monotonically decreasing (e.g., \eqn{1/X}).
#' }
#'
#' @section Available Transformers:
#' \itemize{
#'   \item \code{\link{log_transform}}
#'   \item \code{\link{exp_transform}}
#'   \item \code{\link{logit_transform}}
#'   \item \code{\link{expit_transform}}
#'   \item \code{\link{softplus_transform}}
#'   \item \code{\link{inverse_transform}}
#'   \item \code{\link{sqrt_transform}}
#'   \item \code{\link{power_transform}}
#'   \item \code{\link{bc_transform}} (Box-Cox)
#'   \item \code{\link{yj_transform}} (Yeo-Johnson)
#'   \item \code{\link{affine_transform}}
#' }
#'
#' @seealso \code{\link{transformation}}
NULL




#' Logarithmic Transformation
#'
#' @description
#' Creates a transformer object for the logarithmic transformation \eqn{Y = \log(X)}.
#' Useful for converting distributions defined on \eqn{(0, \infty)} to the real line.
#'
#' @details
#' **Transformation:** \eqn{Y = \log(X)}
#'
#' **Inverse:** \eqn{X = \exp(Y)}
#'
#' **Jacobian:** \eqn{J = \dfrac{dX}{dY} = \exp(Y)}
#'
#' **Support:** Requires \eqn{X > 0}. Maps to \eqn{Y \in (-\infty, \infty)}.
#'
#' @return A list of class \code{"transformer"}.
#' @export
log_transform <- function() {
  # Y = log(X)  =>  X = exp(Y)
  # J = dX/dY = exp(Y)
  # log(J) = Y
  o <- list()
  o$name <- "log"
  o$valid_support <- function(bounds) {
    (bounds[1] >= 0 & bounds[2] >= 0)
  }
  o$bounds_fun <- function(bounds) {
    res <- log(bounds)
    res[bounds == 0] <- -Inf
    res
  }
  o$trans_fun <- log
  o$trans_inv <- exp
  o$trans_abs_jac <- function(y, log = TRUE) {
    if (log) {
      y
    } else {
      exp(y)
    }
  }
  o$decreasing <- FALSE
  class(o) <- "transformer"
  o
}




#' Exponential Transformation
#'
#' @description
#' Creates a transformer object for the exponential transformation \eqn{Y = \exp(X)}.
#' Useful for converting real-valued variables to positive ones.
#'
#' @details
#' **Transformation:** \eqn{Y = \exp(X)}
#'
#' **Inverse:** \eqn{X = \log(Y)}
#'
#' **Jacobian:** \eqn{J = \dfrac{dX}{dY} = \dfrac{1}{Y}}
#'
#' **Support:** Defined for \eqn{X \in (-\infty, \infty)}. Maps to \eqn{Y \in (0, \infty)}.
#'
#' @return A list of class \code{"transformer"}.
#' @export
exp_transform <- function() {
  # Y = exp(X)  =>  X = log(Y)
  # J = dX/dY = 1/Y
  # log(J) = -log(Y)
  o <- list()
  o$name <- "exp"
  o$valid_support <- function(bounds) {
    TRUE
  }
  o$bounds_fun <- exp
  o$trans_fun <- exp
  o$trans_inv <- log
  o$trans_abs_jac <- function(y, log = TRUE) {
    if (log) {
      -log(y)
    } else {
      1 / y
    }
  }
  o$decreasing <- FALSE
  class(o) <- "transformer"
  o
}




#' Reciprocal (Inverse) Transformation
#'
#' @description
#' Creates a transformer object for the reciprocal transformation \eqn{Y = 1/X}.
#'
#' @details
#' **Transformation:** \eqn{Y = 1/X}
#'
#' **Inverse:** \eqn{X = 1/Y}
#'
#' **Jacobian:** \eqn{J = \dfrac{dX}{dY} = -\dfrac{1}{Y^2}} (Absolute value taken for density: \eqn{1/Y^2})
#'
#' **Support:** \eqn{X \neq 0}. Typically used for \eqn{X > 0} or \eqn{X < 0}.
#' Since the function is strictly decreasing, bounds and quantiles are automatically swapped.
#'
#' @return A list of class \code{"transformer"}.
#' @export
inverse_transform <- function() {
  # Y = 1/X  =>  X = 1/Y
  # J = dX/dY = -1/Y^2
  # J (absolute) = 1/Y^2
  # log(J) = -2 * log(|Y|)
  o <- list()
  o$name <- "inverse"
  o$valid_support <- function(bounds) {
    TRUE
  }
  o$bounds_fun <- function(bounds) {
    sort(1 / bounds)
  }
  o$trans_fun <- function(x) {
    1 / x
  }
  o$trans_inv <- function(y) {
    1 / y
  }
  o$trans_abs_jac <- function(y, log = TRUE) {
    if (log) {
      -2 * log(abs(y))
    } else {
      1 / y^2
    }
  }
  o$decreasing <- TRUE
  class(o) <- "transformer"
  o
}




#' Square Root Transformation
#'
#' @description
#' Creates a transformer object for the square root transformation \eqn{Y = \sqrt{X}}.
#'
#' @details
#' **Transformation:** \eqn{Y = \sqrt{X}}
#'
#' **Inverse:** \eqn{X = Y^2}
#'
#' **Jacobian:** \eqn{J = \dfrac{dX}{dY} = 2Y}
#'
#' **Support:** Requires \eqn{X \ge 0}.
#'
#' @return A list of class \code{"transformer"}.
#' @export
sqrt_transform <- function() {
  # Y = sqrt(X)  =>  X = Y^2
  # J = dX/dY = 2Y
  # J (absolute) = 2Y (since Y >= 0)
  # log(J) = log(2) + log(Y)
  o <- list()
  o$name <- "sqrt"
  o$valid_support <- function(bounds) {
    (bounds[1] >= 0)
  }
  o$bounds_fun <- sqrt
  o$trans_fun <- sqrt
  o$trans_inv <- function(y) {
    y^2
  }
  o$trans_abs_jac <- function(y, log = TRUE) {
    if (log) {
      log(2) + log(y)
    } else {
      2 * y
    }
  }
  o$decreasing <- FALSE
  class(o) <- "transformer"
  o
}




#' Power Transformation
#'
#' @description
#' Creates a transformer object for the general power transformation \eqn{Y = X^p}.
#'
#' @param p Numeric. The exponent. Defaults to 2.
#'
#' @details
#' **Transformation:** \eqn{Y = X^p}
#'
#' **Inverse:** \eqn{X = Y^{1/p}}
#'
#' **Jacobian:** \eqn{J = \dfrac{dX}{dY} = \dfrac{1}{p} Y^{\dfrac{1}{p} - 1}}
#'
#' **Support & Constraints:**
#' \itemize{
#'   \item If \eqn{p} is a fractional value, \eqn{X} must be non-negative.
#'   \item If \eqn{p} is an even integer, the transformation is not monotonic on mixed domains (negative to positive). Support is restricted to \eqn{X \ge 0} or \eqn{X \le 0}.
#'   \item If \eqn{p < 0}, the transformation is decreasing.
#' }
#'
#' @return A list of class \code{"transformer"}.
#' @export
power_transform <- function(p = 2) {
  # Y = X^p  =>  X = Y^(1/p)
  # J = dX/dY = (1/p) * Y^(1/p - 1)
  # log(J) = -log(|p|) + (1/p - 1) * log(|Y|)
  o <- list()
  o$name <- paste0("power_", p)
  o$valid_support <- function(bounds) {
    is_integer <- (round(p) == p)
    lb <- bounds[1]
    ub <- bounds[2]
    if (!is_integer && lb < 0) {
      # Fractional powers of negatives forbidden
      FALSE
    } else if (p < 0 && lb <= 0 && ub >= 0) {
      # Division by zero forbidden
      FALSE
    } else if (is_integer && (p %% 2 == 0) && lb < 0 && ub > 0) {
      # EVEN powers on mixed intervals (e.g. -2, 2) are not monotonic
      FALSE
    } else {
      TRUE
    }
  }
  o$bounds_fun <- function(bounds) {
    res <- bounds^p
    # Sort necessary if:
    # 1. p is negative (e.g. x^-1 reverses order)
    # 2. p is even AND we are on negative numbers (e.g. [-3, -2]^2 -> [9, 4])
    is_even_integer <- (round(p) == p && p %% 2 == 0)
    if (p < 0 || (is_even_integer && bounds[2] < 0)) {
      sort(res)
    } else {
      res
    }
  }
  o$trans_fun <- function(x) {
    x^p
  }
  o$trans_inv <- function(y) {
    # Safe handling of roots of negative numbers (e.g. cube root of -8)
    sign(y) * abs(y)^(1 / p)
  }
  o$trans_abs_jac <- function(y, log = TRUE) {
    log_J <- -log(abs(p)) + (1 / p - 1) * log(abs(y))
    if (log) {
      log_J
    } else {
      exp(log_J)
    }
  }
  o$decreasing <- (p < 0)
  class(o) <- "transformer"
  o
}




#' Inverse Hyperbolic Sine Transformation
#'
#' @description
#' Creates a transformer object for the inverse hyperbolic sine transformation \eqn{Y = \text{asinh}(X)}.
#'
#' @details
#' **Transformation:** \eqn{Y = \log(X + \sqrt{X^2 + 1})}
#'
#' **Inverse:** \eqn{X = \sinh(Y)}
#'
#' **Jacobian:** \eqn{J = \dfrac{dX}{dY} = \cosh(Y)}
#'
#' **Support:** Defined for all \eqn{X \in \mathbb{R}}.
#' This is a popular alternative to the log transformation that handles zero and negative values.
#'
#' @return A list of class \code{"transformer"}.
#' @export
asinh_transform <- function() {
  # Y = asinh(X) => X = sinh(Y)
  # J = dX/dY = cosh(Y)
  # log(J) = log(cosh(Y))
  o <- list()
  o$name <- "asinh"
  o$valid_support <- function(bounds) {
    TRUE
  }
  o$bounds_fun <- asinh
  o$trans_fun <- asinh
  o$trans_inv <- sinh
  o$trans_abs_jac <- function(y, log = TRUE) {
    if (log) {
      log(cosh(y))
    } else {
      cosh(y)
    }
  }
  o$decreasing <- FALSE
  class(o) <- "transformer"
  o
}




#' Box-Cox Transformation
#'
#' @description
#' Creates a transformer object for the one-parameter Box-Cox transformation.
#'
#' @param lambda Numeric. The transformation parameter.
#'
#' @details
#' **Transformation:**
#' \itemize{
#'   \item If \eqn{\lambda \neq 0}: \eqn{Y = \dfrac{X^\lambda - 1}{\lambda}}
#'   \item If \eqn{\lambda = 0}: \eqn{Y = \log(X)}
#' }
#'
#' **Inverse:** \eqn{X = (\lambda Y + 1)^{1/\lambda}}
#'
#' **Support:** Requires \eqn{X > 0}.
#'
#' @return A list of class \code{"transformer"}.
#' @export
bc_transform <- function(lambda) {
  # Box-Cox Transformation
  # If lambda -> 0: Log transform
  # If lambda != 0: Y = (X^lambda - 1) / lambda
  # Inverse: X = (lambda * Y + 1)^(1/lambda)

  # Handle the limiting case where lambda is close to 0
  if (abs(lambda) < 1e-10) {
    l <- log_transform()
    l$name <- "box_cox_0"
    return(l)
  }

  o <- list()
  o$name <- paste0("box_cox_", lambda)

  # Strictly defined for positive data (x > 0)
  o$valid_support <- function(bounds) {
    bounds[1] > 0
  }

  o$bounds_fun <- function(bounds) {
    res <- (bounds^lambda - 1) / lambda
    # If lambda < 0, the function is decreasing
    if (lambda < 0) {
      sort(res)
    } else {
      res
    }
  }

  o$trans_fun <- function(x) {
    (x^lambda - 1) / lambda
  }

  o$trans_inv <- function(y) {
    base <- lambda * y + 1
    # Numerical safety clip for floating point errors near zero
    base[base < 0] <- 0
    base^(1 / lambda)
  }

  o$trans_abs_jac <- function(y, log = TRUE) {
    # J = dX/dY = (lambda * Y + 1)^((1 - lambda) / lambda)
    term <- lambda * y + 1
    term[term < 1e-16] <- 1e-16 # Avoid log(0) issues

    log_J <- ((1 - lambda) / lambda) * log(term)

    if (log) {
      log_J
    } else {
      exp(log_J)
    }
  }

  o$decreasing <- (lambda < 0)
  class(o) <- "transformer"
  o
}




#' Yeo-Johnson Transformation
#'
#' @description
#' Creates a transformer object for the Yeo-Johnson transformation, which extends Box-Cox to negative values.
#'
#' @param lambda Numeric. The transformation parameter.
#'
#' @details
#' **Transformation:** Defined piecewise:
#' \itemize{
#'   \item \eqn{X \ge 0}: Box-Cox of \eqn{X+1} with parameter \eqn{\lambda}.
#'   \item \eqn{X < 0}: Box-Cox of \eqn{|X|+1} with parameter \eqn{2-\lambda}, negated.
#' }
#'
#' **Support:** Defined for all \eqn{X \in \mathbb{R}}.
#'
#' @return A list of class \code{"transformer"}.
#' @export
yj_transform <- function(lambda) {
  # Yeo-Johnson Transformation (Piecewise Box-Cox for real line)
  # Inverse logic is split based on Y >= 0 or Y < 0
  # Jacobian computed via inverse function theorem: log(J) = -log(dY/dX)

  o <- list()
  o$name <- paste0("yeo_johnson_", lambda)

  o$valid_support <- function(bounds) {
    TRUE
  }

  o$trans_fun <- function(x) {
    y <- numeric(length(x))
    pos <- x >= 0
    neg <- !pos

    # Case 1: x >= 0 (Box-Cox with lambda + 1)
    if (abs(lambda) < 1e-10) {
      y[pos] <- log1p(x[pos])
    } else {
      y[pos] <- ((x[pos] + 1)^lambda - 1) / lambda
    }

    # Case 2: x < 0 (Box-Cox on exp(-x) with 2 - lambda)
    lam2 <- 2 - lambda
    if (abs(lam2) < 1e-10) {
      y[neg] <- -log1p(-x[neg])
    } else {
      y[neg] <- -((-x[neg] + 1)^lam2 - 1) / lam2
    }
    y
  }

  o$bounds_fun <- function(bounds) {
    # Monotonic non-decreasing
    o$trans_fun(bounds)
  }

  o$trans_inv <- function(y) {
    x <- numeric(length(y))
    pos <- y >= 0
    neg <- !pos

    # Inverse for Y >= 0
    if (abs(lambda) < 1e-10) {
      x[pos] <- expm1(y[pos])
    } else {
      x[pos] <- (lambda * y[pos] + 1)^(1 / lambda) - 1
    }

    # Inverse for Y < 0
    lam2 <- 2 - lambda
    if (abs(lam2) < 1e-10) {
      x[neg] <- -expm1(-y[neg])
    } else {
      x[neg] <- 1 - (1 - lam2 * y[neg])^(1 / lam2)
    }
    x
  }

  o$trans_abs_jac <- function(y, log = TRUE) {
    # 1. Recover X to use the simpler forward derivative formula
    x <- o$trans_inv(y)

    log_dy_dx <- numeric(length(x))
    pos <- x >= 0
    neg <- !pos

    if (any(pos)) {
      log_dy_dx[pos] <- (lambda - 1) * log1p(x[pos])
    }
    if (any(neg)) {
      log_dy_dx[neg] <- (1 - lambda) * log1p(-x[neg])
    }

    # 2. Invert gradient (log(J_inv) = -log(J_fwd))
    log_J <- -log_dy_dx

    if (log) {
      log_J
    } else {
      exp(log_J)
    }
  }

  o$decreasing <- FALSE
  class(o) <- "transformer"
  o
}




#' Affine (Location-Scale) Transformation
#'
#' @description
#' Creates a transformer object for an affine linear transformation \eqn{Y = a + bX}.
#'
#' @param loc Numeric. The location shift \eqn{a}. Defaults to 0.
#' @param scale Numeric. The scale multiplier \eqn{b}. Defaults to 1.
#'
#' @details
#' **Transformation:** \eqn{Y = \text{loc} + \text{scale} \cdot X}
#'
#' **Inverse:** \eqn{X = (Y - \text{loc}) / \text{scale}}
#'
#' **Jacobian:** \eqn{J = 1 / |\text{scale}|}
#'
#' **Support:** Defined for all \eqn{X \in \mathbb{R}}.
#'
#' @return A list of class \code{"transformer"}.
#' @export
affine_transform <- function(loc = 0, scale = 1) {
  # Y = loc + scale * X   =>   X = (Y - loc) / scale
  # J = dX/dY = 1 / scale
  # J (absolute) = 1 / |scale|

  if (scale == 0) {
    stop("Scale parameter cannot be equal to zero.")
  }

  o <- list()
  o$name <- paste0("affine_loc", loc, "_scale", scale)

  o$valid_support <- function(bounds) {
    TRUE
  }

  o$bounds_fun <- function(bounds) {
    res <- loc + scale * bounds
    # If scale is negative, the bounds flip
    if (scale < 0) {
      sort(res)
    } else {
      res
    }
  }

  o$trans_fun <- function(x) {
    loc + scale * x
  }

  o$trans_inv <- function(y) {
    (y - loc) / scale
  }

  o$trans_abs_jac <- function(y, log = TRUE) {
    val <- 1 / abs(scale)
    if (log) {
      log(val)
    } else {
      val
    }
  }

  o$decreasing <- (scale < 0)
  class(o) <- "transformer"
  o
}




#' Logit Transformation
#'
#' @description
#' Creates a transformer object for the logit transformation \eqn{Y = \text{logit}(X)}.
#'
#' @details
#' **Transformation:** \eqn{Y = \log\left(\dfrac{X}{1-X}\right)}
#'
#' **Inverse:** \eqn{X = \dfrac{1}{1 + e^{-Y}}} (Sigmoid/Expit)
#'
#' **Jacobian:** \eqn{J = \dfrac{dX}{dY} = X(1-X)} (corresponds to the PDF of the Logistic distribution).
#'
#' **Support:** Requires \eqn{X \in (0, 1)}.
#'
#' @return A list of class \code{"transformer"}.
#' @export
logit_transform <- function() {
  # Transformation: Y = logit(X) = log(X / (1 - X))
  # Inverse:        X = expit(Y) = exp(Y) / (1 + exp(Y))
  # Jacobian:       J = dX/dY = exp(Y) / (1 + exp(Y))^2
  # Note: The Jacobian corresponds exactly to the PDF of the standard logistic distribution (dlogis).

  o <- list()
  o$name <- "logit"

  o$valid_support <- function(bounds) {
    # Support must be within (0, 1). Closed interval [0, 1] is allowed theoretically,
    # but 0 and 1 map to -Inf and Inf respectively.
    (bounds[1] >= 0 & bounds[2] <= 1)
  }

  o$bounds_fun <- function(bounds) {
    # qlogis correctly handles 0 -> -Inf and 1 -> Inf
    stats::qlogis(bounds)
  }

  o$trans_fun <- stats::qlogis
  o$trans_inv <- stats::plogis

  o$trans_abs_jac <- function(y, log = TRUE) {
    # The derivative of the inverse function (plogis) is the logistic density (dlogis).
    stats::dlogis(y, log = log)
  }

  o$decreasing <- FALSE
  class(o) <- "transformer"
  o
}




#' Expit (Sigmoid) Transformation
#'
#' @description
#' Creates a transformer object for the expit (inverse logit) transformation.
#' Useful for mapping real-valued variables to the \eqn{(0, 1)} interval.
#'
#' @details
#' **Transformation:** \eqn{Y = \dfrac{1}{1 + e^{-X}}}
#'
#' **Inverse:** \eqn{X = \text{logit}(Y)}
#'
#' **Jacobian:** \eqn{J = \dfrac{1}{Y(1-Y)}}
#'
#' **Support:** Defined for \eqn{X \in \mathbb{R}}. Maps to \eqn{Y \in (0, 1)}.
#'
#' @return A list of class \code{"transformer"}.
#' @export
expit_transform <- function() {
  # Transformation: Y = expit(X) = exp(X) / (1 + exp(X))
  # Inverse:        X = logit(Y) = log(Y / (1 - Y))
  # Jacobian:       J = dX/dY = 1 / (Y * (1 - Y))
  # log(J) = -log(Y * (1 - Y)) = -log(Y) - log(1 - Y)

  o <- list()
  o$name <- "expit"

  o$valid_support <- function(bounds) {
    TRUE
  }

  o$bounds_fun <- function(bounds) {
    # Maps (-Inf, Inf) to (0, 1)
    stats::plogis(bounds)
  }

  o$trans_fun <- stats::plogis
  o$trans_inv <- stats::qlogis

  o$trans_abs_jac <- function(y, log = TRUE) {
    # We are in the domain y in (0, 1).
    # J = 1 / (y * (1 - y))

    # Using log1p(-y) is more numerically stable than log(1 - y) for small y.
    # log(J) = - (log(y) + log(1 - y))
    val <- -(log(y) + log1p(-y))

    if (log) {
      val
    } else {
      exp(val)
    }
  }

  o$decreasing <- FALSE
  class(o) <- "transformer"
  o
}





#' Softplus Transformation
#'
#' @description
#' Creates a transformer object for the Softplus function.
#' This transformation is used to map a non-negative random variable \eqn{X \in (0, \infty)}
#' to the entire real line \eqn{Y \in \mathbb{R}}.
#'
#' Unlike the logarithmic/exponential transformation, the Softplus function approaches linearity
#' for large values (\eqn{X \approx Y} as \eqn{Y \to \infty}), while behaving like an exponential
#' near zero. This often results in better numerical stability for heavy-tailed distributions.
#'
#' @param a Numeric. Scale parameter. Defaults to 1.
#'
#' @details
#' **Transformation (Inverse Softplus):** \eqn{Y = \dfrac{1}{a} \log(e^{aX} - 1)}
#'
#' **Inverse (Softplus):** \eqn{X = \dfrac{1}{a} \log(1 + e^{aY})}
#'
#' **Jacobian:** \eqn{J = \text{sigmoid}(aY)}
#'
#' **Support:** Requires \eqn{X > 0}. Maps to \eqn{Y \in \mathbb{R}}.
#'
#' @return A list of class \code{"transformer"}.
#' @export
softplus_transform <- function(a = 1) {
  # Transformation: Y = 1/a * log(exp(a*X) - 1)  [Inverse Softplus]
  # Inverse:        X = 1/a * log(1 + exp(a*Y))  [Softplus]
  # Jacobian:       dX/dY = sigmoid(a*Y) = 1 / (1 + exp(-a*Y))
  #
  # This maps a positive variable X (0, Inf) to the real line Y (-Inf, Inf).
  # Often used to ensure positivity of parameters (like variance).

  if (a <= 0) {
    stop("Scale parameter 'a' must be greater than 0")
  }

  o <- list()
  o$name <- paste0("softplus_a", round(a, 5))

  # Strictly defined for non-negative data (X >= 0)
  # X=0 maps to Y=-Inf
  o$valid_support <- function(bounds) {
    bounds[1] >= 0
  }

  o$bounds_fun <- function(bounds) {
    # Forward transform: log(expm1(a * x)) / a
    # Handles 0 -> -Inf correctly
    res <- log(expm1(a * bounds)) / a
    res[bounds == 0] <- -Inf
    res
  }

  o$trans_fun <- function(x) {
    # Forward: X -> Y
    # log(expm1(z)) is stable for small z.
    # For very large z, expm1(z) ~ exp(z), so result ~ z/a.
    log(expm1(a * x)) / a
  }

  o$trans_inv <- function(y) {
    # Inverse: Y -> X (The Softplus function)
    # Analytical: log(1 + exp(a*y)) / a
    # Numerical trick: to avoid overflow when a*y is large positive,
    # we rewrite as: y + log(1 + exp(-a*y)) / a  (if y > 0)

    # We use the robust implementation provided:
    pmax(0, y) + log1p(exp(-abs(a * y))) / a
  }

  o$trans_abs_jac <- function(y, log = TRUE) {
    # Jacobian dX/dY is the Sigmoid function: 1 / (1 + exp(-a*y))
    # This corresponds exactly to the CDF of the Logistic distribution.
    # R's plogis is numerically stable and supports log.p directly.
    stats::plogis(a * y, log.p = log)
  }

  o$decreasing <- FALSE
  class(o) <- "transformer"
  o
}





#' Apply a Variable Transformation to a `distrib` Object
#'
#' @description
#' Creates a new distribution object by applying a bijective transformation \eqn{Y = g(X)}
#' to an existing continuous \code{"distrib"} object.
#'
#' @param distrib An object of class \code{"distrib"}. Must be a continuous distribution.
#' @param transformer An object of class \code{"transformer"} containing the transformation logic
#'   (forward, inverse, Jacobian, etc.), created by functions like
#'   \code{\link{log_transform}}, \code{\link{affine_transform}}, etc.
#'
#' @details
#' **Mathematical Formulation:**
#' This function implements the Change of Variable technique. If \eqn{X \sim f_X(x)} and \eqn{Y = g(X)},
#' then the density of \eqn{Y} is given by:
#' \deqn{f_Y(y) = f_X(g^{-1}(y)) \cdot \left| \dfrac{d}{dy} g^{-1}(y) \right|}
#' where \eqn{J = \dfrac{d}{dy} g^{-1}(y)} is the Jacobian of the inverse transformation.
#'
#' **Components Adjusted:**
#' The returned \code{"distrib"} object is modified as follows:
#' \itemize{
#'   \item **PDF:** Computed using the formula above.
#'   \item **CDF:** \eqn{F_Y(y) = F_X(g^{-1}(y))}. If the transformation is decreasing, \eqn{F_Y(y) = 1 - F_X(g^{-1}(y))}.
#'   \item **Quantile Function:** \eqn{Q_Y(p) = g(Q_X(p))}. If decreasing, \eqn{Q_Y(p) = g(Q_X(1-p))}.
#'   \item **Support:** The bounds are transformed via \eqn{g(\text{bounds})}.
#'   \item **Moments:** Analytical moments (mean, variance, skewness, kurtosis) are replaced by numerical integration methods
#'     (\code{\link{moment}}) since closed-form solutions are generally lost after transformation.
#'   \item **RNG:** Generates \eqn{x \sim X} and returns \eqn{g(x)}.
#'   \item **Mode:** Computed numerically via \code{\link{mode.distrib}}.
#' }
#'
#' **Derivatives of the Log-Likelihood:**
#' When a distribution is transformed, the log-likelihood for an observation \eqn{y} is:
#' \deqn{\ell(\theta; y) = \log f_X(g^{-1}(y); \theta) + \log |J(y)|}
#' where \eqn{J(y)} is the Jacobian of the inverse transformation. Since the transformation \eqn{g}
#' does not depend on the parameters \eqn{\theta}, the Jacobian term \eqn{\log |J(y)|} is a constant
#' with respect to \eqn{\theta}. Consequently:
#' \itemize{
#'   \item **Gradient (Score):** The gradient of the transformed distribution is identical to the
#'     gradient of the original distribution evaluated at the inverse-transformed point:
#'     \deqn{\nabla_\theta \ell(\theta; y) = \nabla_\theta \log f_X(x; \theta) \mid_{x=g^{-1}(y)}}
#'   \item **Hessian (Observed):** The observed Hessian follows the same logic:
#'     \deqn{\nabla^2_\theta \ell(\theta; y) = \nabla^2_\theta \log f_X(x; \theta) \mid_{x=g^{-1}(y)}}
#'   \item **Expected Hessian (Fisher Information):** Calculated via Monte Carlo integration
#'     (\code{\link{mc_expected_hessian}}) using the identity
#'     \eqn{\mathbb{E}[\nabla^2 \ell] = -\mathbb{E}[\nabla \ell (\theta)(\nabla \ell (\theta))^\top]}.
#' }
#'
#' **Numerical Stability:**
#' The function includes safeguards for calculating the PDF in the tails:
#' \itemize{
#'   \item Calculations are performed in log-space to avoid underflow.
#'   \item Infinite log-densities (singularities) are clamped to \code{.Machine$double.xmax} to prevent \code{Inf - Inf} issues
#'     during numerical integration.
#' }
#'
#' @return A new \code{"distrib"} object representing the transformed variable.
#'
#' @seealso \code{\link{log_transform}}, \code{\link{affine_transform}}
#' @export
transformation <- function(distrib, transformer) {
  d <- distrib
  tr <- transformer

  if (!inherits(d, "distrib")) {
    stop("Argument 'distrib' must be a distrib object.")
  }

  if (d$type != "continuous") {
    stop("This function currently supports only continuous distribution")
  }

  if (!inherits(tr, "transformer")) {
    stop("Argument 'transformer' must be a transformer object.")
  }
  if (!transformer$valid_support(distrib$bounds)) {
    stop(sprintf(
      "The '%s' transformation is not valid for the support of the '%s' distribution.",
      tr$name, d$distrib_name
    ))
  }

  old_pdf <- distrib$pdf
  old_cdf <- distrib$cdf
  old_quantile <- distrib$quantile
  old_rng <- distrib$rng
  old_gradient <- distrib$gradient
  old_hessian <- distrib$hessian
  old_kernel <- distrib$kernel
  old_initialize <- distrib$initialize
  old_median <- distrib$median

  d$distrib_name <- paste0(tr$name, "(", distrib$distrib_name, ")")
  d$bounds <- tr$bounds_fun(distrib$bounds)
  d$params_interpretation <- paste0(distrib$params_interpretation, " (", distrib$distrib_name, " scale)")
  names(d$params_interpretation) <- d$params

  d$pdf <- function(y, theta, log = FALSE) {
    log_pdf_x <- old_pdf(tr$trans_inv(y), theta, log = TRUE)
    log_J <- tr$trans_abs_jac(y, log = TRUE)

    if (any(is.infinite(log_pdf_x) & log_pdf_x > 0)) {
      log_pdf_x[is.infinite(log_pdf_x) & log_pdf_x > 0] <- log(.Machine$double.xmax)
    }
    log_val <- log_pdf_x + log_J

    if (log) {
      log_val
    } else {
      exp(log_val)
    }
  }

  d$cdf <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {
    if (tr$decreasing) {
      lower.tail <- !lower.tail
    }
    old_cdf(tr$trans_inv(q), theta, lower.tail = lower.tail, log.p = log.p)
  }

  d$quantile <- function(p, theta, lower.tail = TRUE, log.p = FALSE) {
    if (tr$decreasing) {
      lower.tail <- !lower.tail
    }
    tr$trans_fun(old_quantile(p, theta, lower.tail = lower.tail, log.p = log.p))
  }

  d$rng <- function(n, theta) {
    tr$trans_fun(old_rng(n, theta))
  }

  d$gradient <- function(y, theta, par = NULL) {
    old_gradient(tr$trans_inv(y), theta, par)
  }

  d$hessian <- function(y, theta, expected = FALSE) {
    if (expected) {
      mc_expected_hessian(d, y, theta, nsim = 1000)
    } else {
      old_hessian(tr$trans_inv(y), theta, expected = FALSE)
    }
  }

  d$kernel <- function(y, theta, log = TRUE) {
    k <- old_kernel(tr$trans_inv(y), theta, log = TRUE) + tr$trans_abs_jac(y, log = TRUE)
    if (log) {
      k
    } else {
      exp(k)
    }
  }

  d$mean <- function(theta) {
    mean.distrib(d, theta, use_moment = TRUE)
  }

  d$variance <- function(theta) {
    variance.distrib(d, theta, use_moment = TRUE)
  }

  d$skewness <- function(theta) {
    skewness.distrib(d, theta, use_moment = TRUE)
  }

  d$kurtosis <- function(theta) {
    kurtosis.distrib(d, theta, use_moment = TRUE)
  }

  d$mode <- function(theta) {
    mode.distrib(d, theta)
  }

  d$median <- function(theta) {
    tr$trans_fun(old_median(theta))
  }

  d$initialize <- function(y) {
    old_initialize(tr$trans_inv(y))
  }
  d
}
