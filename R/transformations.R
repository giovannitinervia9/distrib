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
#'   \item \strong{trans_abs_jac}: The absolute value of the Jacobian of the inverse transformation
#'     \eqn{\left|J(y)\right| = \left|\dfrac{dx}{dy}\right|}.
#'     It must accept a \code{log} argument to return the log-absolute Jacobian.
#'   \item \strong{trans_inv_hessian}: A function of \code{y} that returns the second derivative
#'     of the inverse transformation: \eqn{\dfrac{d^2 x}{d y^2}}.
#'   \item \strong{grad_log_jac}: A function of \code{y} that returns the first derivative of the
#'     log-absolute Jacobian with respect to \eqn{y}: \eqn{\dfrac{\partial \log |J(y)|}{\partial y}}.
#'   \item \strong{hess_log_jac}: A function of \code{y} that returns the second derivative of the
#'     log-absolute Jacobian with respect to \eqn{y}: \eqn{\dfrac{\partial^2 \log |J(y)|}{\partial y^2}}.
#'   \item \strong{bounds_fun}: A function that takes the original distribution bounds and
#'     returns the transformed support bounds.
#'   \item \strong{valid_support}: A function that checks if the input distribution's support
#'     is mathematically compatible with the transformation.
#'   \item \strong{decreasing}: A logical value indicating if the transformation is
#'     monotonically decreasing (e.g., \eqn{1/X}).
#' }
#'
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
#' **Log-Jacobian Derivatives:** #' \itemize{
#'   \item Gradient: \eqn{\dfrac{\partial \log|J|}{\partial y} = 1}
#'   \item Hessian: \eqn{\dfrac{\partial^2 \log|J|}{\partial y^2} = 0}
#' }
#'
#' **Inverse Hessian:** \eqn{\dfrac{d^2 X}{d Y^2} = \exp(Y)}
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
  o$trans_inv_hessian <- function(y) {
    exp(y)
  }
  o$grad_log_jac <- function(y) {
    rep(1, length(y))
  }

  o$hess_log_jac <- function(y) {
    rep(0, length(y))
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
#' **Log-Jacobian Derivatives:**
#' \itemize{
#'   \item Gradient: \eqn{\dfrac{\partial \log|J|}{\partial y} = -\dfrac{1}{y}}
#'   \item Hessian: \eqn{\dfrac{\partial^2 \log|J|}{\partial y^2} = \dfrac{1}{y^2}}
#' }
#'
#' **Inverse Hessian:** \eqn{\dfrac{d^2 X}{d Y^2} = -\dfrac{1}{Y^2}}
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

  o$grad_log_jac <- function(y) {
    -1 / y
  }

  o$hess_log_jac <- function(y) {
    1 / y^2
  }

  o$trans_inv_hessian <- function(y) {
    -1 / y^2
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
#' **Log-Jacobian Derivatives:**
#' \itemize{
#'   \item Gradient: \eqn{\dfrac{\partial \log|J|}{\partial y} = -\dfrac{2}{y}}
#'   \item Hessian: \eqn{\dfrac{\partial^2 \log|J|}{\partial y^2} = \dfrac{2}{y^2}}
#' }
#'
#' **Inverse Hessian:** \eqn{\dfrac{d^2 X}{d Y^2} = \dfrac{2}{Y^3}}
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

  o$grad_log_jac <- function(y) {
    -2 / y
  }

  o$hess_log_jac <- function(y) {
    2 / y^2
  }

  o$trans_inv_hessian <- function(y) {
    2 / y^3
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
#' **Log-Jacobian Derivatives:**
#' \itemize{
#'   \item Gradient: \eqn{\dfrac{\partial \log|J|}{\partial y} = \dfrac{1}{y}}
#'   \item Hessian: \eqn{\dfrac{\partial^2 \log|J|}{\partial y^2} = -\dfrac{1}{y^2}}
#' }
#'
#' **Inverse Hessian:** \eqn{\dfrac{d^2 X}{d Y^2} = 2}
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

  o$grad_log_jac <- function(y) {
    1 / y
  }

  o$hess_log_jac <- function(y) {
    -1 / y^2
  }

  o$trans_inv_hessian <- function(y) {
    # d2(y^2)/dy2 = 2
    rep(2, length(y))
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
#' **Log-Jacobian Derivatives:**
#' \itemize{
#'   \item Gradient: \eqn{\dfrac{\partial \log|J|}{\partial y} = (\frac{1}{p} - 1) \frac{1}{y}}
#'   \item Hessian: \eqn{\dfrac{\partial^2 \log|J|}{\partial y^2} = -(\frac{1}{p} - 1) \frac{1}{y^2}}
#' }
#'
#' **Inverse Hessian:** \eqn{\dfrac{d^2 X}{d Y^2} = \frac{1}{p}(\frac{1}{p} - 1) Y^{\frac{1}{p} - 2}}
#'
#' **Support & Constraints:**
#' \itemize{
#'   \item If \eqn{p} is a fractional value, \eqn{X} must be non-negative.
#'   \item If \eqn{p} is an even integer, the transformation is not monotonic on mixed domains.
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
      FALSE
    } else if (p < 0 && lb <= 0 && ub >= 0) {
      FALSE
    } else if (is_integer && (p %% 2 == 0) && lb < 0 && ub > 0) {
      FALSE
    } else {
      TRUE
    }
  }

  o$bounds_fun <- function(bounds) {
    res <- bounds^p
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
    sign(y) * abs(y)^(1 / p)
  }

  o$trans_abs_jac <- function(y, log = TRUE) {
    log_J <- -log(abs(p)) + (1 / p - 1) * log(abs(y))
    if (log) log_J else exp(log_J)
  }

  o$grad_log_jac <- function(y) {
    (1 / p - 1) / y
  }

  o$hess_log_jac <- function(y) {
    -(1 / p - 1) / y^2
  }

  o$trans_inv_hessian <- function(y) {
    # d2x/dy2 = 1/p * (1/p - 1) * y^(1/p - 2)
    (1 / p) * (1 / p - 1) * sign(y) * abs(y)^(1 / p - 2)
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
#' **Transformation:** \eqn{Y = \text{asinh}(X) = \log(X + \sqrt{X^2 + 1})}
#'
#' **Inverse:** \eqn{X = \sinh(Y)}
#'
#' **Jacobian:** \eqn{J = \dfrac{dX}{dY} = \cosh(Y)}
#'
#' **Log-Jacobian Derivatives:**
#' \itemize{
#'   \item Gradient: \eqn{\dfrac{\partial \log|J|}{\partial y} = \tanh(y)}
#'   \item Hessian: \eqn{\dfrac{\partial^2 \log|J|}{\partial y^2} = 1 - \tanh^2(y)}
#' }
#'
#' **Inverse Hessian:** \eqn{\dfrac{d^2 X}{d Y^2} = \sinh(Y)}
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

  o$grad_log_jac <- function(y) {
    tanh(y)
  }

  o$hess_log_jac <- function(y) {
    # sech^2(y) = 1 - tanh^2(y)
    1 - tanh(y)^2
  }

  o$trans_inv_hessian <- function(y) {
    # d2(sinh(y))/dy2 = sinh(y)
    sinh(y)
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
#' **Log-Jacobian Derivatives:**
#' \itemize{
#'   \item Gradient: \eqn{\dfrac{\partial \log|J|}{\partial y} = \dfrac{1-\lambda}{\lambda y + 1}}
#'   \item Hessian: \eqn{\dfrac{\partial^2 \log|J|}{\partial y^2} = -\dfrac{\lambda(1-\lambda)}{(\lambda y + 1)^2}}
#' }
#'
#' **Inverse Hessian:** \eqn{\dfrac{d^2 X}{d Y^2} = (1-\lambda)(\lambda Y + 1)^{\frac{1-2\lambda}{\lambda}}}
#'
#' **Support:** Requires \eqn{X > 0}.
#'
#' @return A list of class \code{"transformer"}.
#' @export
bc_transform <- function(lambda) {
  # Handle the limiting case where lambda is close to 0
  if (abs(lambda) < 1e-10) {
    l <- log_transform()
    l$name <- "box_cox_0"
    return(l)
  }

  o <- list()
  o$name <- paste0("box_cox_", lambda)

  o$valid_support <- function(bounds) {
    bounds[1] > 0
  }

  o$bounds_fun <- function(bounds) {
    res <- (bounds^lambda - 1) / lambda
    if (lambda < 0) sort(res) else res
  }

  o$trans_fun <- function(x) {
    (x^lambda - 1) / lambda
  }

  o$trans_inv <- function(y) {
    base <- lambda * y + 1
    base[base < 0] <- 0
    base^(1 / lambda)
  }

  o$trans_abs_jac <- function(y, log = TRUE) {
    term <- lambda * y + 1
    term[term < 1e-16] <- 1e-16
    log_J <- ((1 - lambda) / lambda) * log(term)
    if (log) log_J else exp(log_J)
  }

  o$grad_log_jac <- function(y) {
    term <- lambda * y + 1
    (1 - lambda) / term
  }

  o$hess_log_jac <- function(y) {
    term <- lambda * y + 1
    -(lambda * (1 - lambda)) / (term^2)
  }

  o$trans_inv_hessian <- function(y) {
    term <- lambda * y + 1
    term[term < 0] <- 0
    (1 - lambda) * term^((1 - 2 * lambda) / lambda)
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
#' **Log-Jacobian Derivatives:**
#' \itemize{
#'   \item Gradient: \eqn{\dfrac{\partial \log|J|}{\partial y} = \dfrac{1-\lambda}{\lambda y + 1}} for \eqn{y \ge 0}, and \eqn{\dfrac{\lambda-1}{1-(2-\lambda)y}} for \eqn{y < 0}.
#'   \item Hessian: \eqn{\dfrac{\partial^2 \log|J|}{\partial y^2} = -\dfrac{\lambda(1-\lambda)}{(\lambda y + 1)^2}} for \eqn{y \ge 0}, and \eqn{\dfrac{(\lambda-1)(2-\lambda)}{(1-(2-\lambda)y)^2}} for \eqn{y < 0}.
#' }
#'
#' **Inverse Hessian:** \eqn{\dfrac{d^2 X}{d Y^2}} is computed piecewise based on the second derivative of the inverse Box-Cox components.
#'
#' **Support:** Defined for all \eqn{X \in \mathbb{R}}.
#'
#' @return A list of class \code{"transformer"}.
#' @export
yj_transform <- function(lambda) {
  o <- list()
  o$name <- paste0("yeo_johnson_", lambda)

  o$valid_support <- function(bounds) {
    TRUE
  }

  o$trans_fun <- function(x) {
    y <- numeric(length(x))
    pos <- x >= 0
    neg <- !pos
    if (abs(lambda) < 1e-10) {
      y[pos] <- log1p(x[pos])
    } else {
      y[pos] <- ((x[pos] + 1)^lambda - 1) / lambda
    }
    lam2 <- 2 - lambda
    if (abs(lam2) < 1e-10) {
      y[neg] <- -log1p(-x[neg])
    } else {
      y[neg] <- -((-x[neg] + 1)^lam2 - 1) / lam2
    }
    y
  }

  o$bounds_fun <- function(bounds) {
    o$trans_fun(bounds)
  }

  o$trans_inv <- function(y) {
    x <- numeric(length(y))
    pos <- y >= 0
    neg <- !pos
    if (abs(lambda) < 1e-10) {
      x[pos] <- expm1(y[pos])
    } else {
      x[pos] <- (lambda * y[pos] + 1)^(1 / lambda) - 1
    }
    lam2 <- 2 - lambda
    if (abs(lam2) < 1e-10) {
      x[neg] <- -expm1(-y[neg])
    } else {
      x[neg] <- 1 - (1 - lam2 * y[neg])^(1 / lam2)
    }
    x
  }

  o$trans_abs_jac <- function(y, log = TRUE) {
    log_J <- numeric(length(y))
    pos <- y >= 0
    neg <- !pos
    if (any(pos)) {
      log_J[pos] <- ((1 - lambda) / lambda) * log(lambda * y[pos] + 1)
    }
    if (any(neg)) {
      lam2 <- 2 - lambda
      log_J[neg] <- ((lambda - 1) / lam2) * log(1 - lam2 * y[neg])
    }
    if (log) log_J else exp(log_J)
  }

  o$grad_log_jac <- function(y) {
    grad <- numeric(length(y))
    pos <- y >= 0
    neg <- !pos
    if (any(pos)) {
      grad[pos] <- (1 - lambda) / (lambda * y[pos] + 1)
    }
    if (any(neg)) {
      lam2 <- 2 - lambda
      grad[neg] <- (lambda - 1) / (1 - lam2 * y[neg])
    }
    grad
  }

  o$hess_log_jac <- function(y) {
    hess <- numeric(length(y))
    pos <- y >= 0
    neg <- !pos
    if (any(pos)) {
      hess[pos] <- -(lambda * (1 - lambda)) / (lambda * y[pos] + 1)^2
    }
    if (any(neg)) {
      lam2 <- 2 - lambda
      hess[neg] <- ((lambda - 1) * lam2) / (1 - lam2 * y[neg])^2
    }
    hess
  }

  o$trans_inv_hessian <- function(y) {
    h <- numeric(length(y))
    pos <- y >= 0
    neg <- !pos
    if (any(pos)) {
      h[pos] <- (1 - lambda) * (lambda * y[pos] + 1)^((1 - 2 * lambda) / lambda)
    }
    if (any(neg)) {
      lam2 <- 2 - lambda
      # Derivata seconda di 1 - (1 - lam2*y)^(1/lam2)
      h[neg] <- (lam2 - 1) * (1 - lam2 * y[neg])^((1 - 2 * lam2) / lam2)
    }
    h
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
#' **Jacobian:** \eqn{J = \dfrac{dX}{dY} = 1 / \text{scale}}
#'
#' **Log-Jacobian Derivatives:**
#' \itemize{
#'   \item Gradient: \eqn{\dfrac{\partial \log|J|}{\partial y} = 0}
#'   \item Hessian: \eqn{\dfrac{\partial^2 \log|J|}{\partial y^2} = 0}
#' }
#'
#' **Inverse Hessian:** \eqn{\dfrac{d^2 X}{d Y^2} = 0}
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

  # --- Nuovi slot per derivate rispetto a y ---

  o$grad_log_jac <- function(y) {
    rep(0, length(y))
  }

  o$hess_log_jac <- function(y) {
    rep(0, length(y))
  }

  o$trans_inv_hessian <- function(y) {
    rep(0, length(y))
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
#' **Log-Jacobian Derivatives:**
#' \itemize{
#'   \item Gradient: \eqn{\dfrac{\partial \log|J|}{\partial y} = 1 - 2\text{plogis}(y)}
#'   \item Hessian: \eqn{\dfrac{\partial^2 \log|J|}{\partial y^2} = -2\text{dlogis}(y)}
#' }
#'
#' **Inverse Hessian:** \eqn{\dfrac{d^2 X}{d Y^2} = \text{dlogis}(y)(1 - 2\text{plogis}(y))}
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
    (bounds[1] >= 0 & bounds[2] <= 1)
  }

  o$bounds_fun <- function(bounds) {
    stats::qlogis(bounds)
  }

  o$trans_fun <- stats::qlogis
  o$trans_inv <- stats::plogis

  o$trans_abs_jac <- function(y, log = TRUE) {
    stats::dlogis(y, log = log)
  }

  o$grad_log_jac <- function(y) {
    # Derivative of log(dlogis(y)) = 1 - 2*plogis(y)
    1 - 2 * stats::plogis(y)
  }

  o$hess_log_jac <- function(y) {
    # Derivative of (1 - 2*plogis(y)) = -2*dlogis(y)
    -2 * stats::dlogis(y)
  }

  o$trans_inv_hessian <- function(y) {
    # Second derivative of plogis(y): dlogis(y)*(1 - 2*plogis(y))
    stats::dlogis(y) * (1 - 2 * stats::plogis(y))
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
#' **Jacobian:** \eqn{J = \dfrac{dX}{dY} = \dfrac{1}{Y(1-Y)}}
#'
#' **Log-Jacobian Derivatives:**
#' \itemize{
#'   \item Gradient: \eqn{\dfrac{\partial \log|J|}{\partial y} = \dfrac{2y - 1}{y(1-y)}}
#'   \item Hessian: \eqn{\dfrac{\partial^2 \log|J|}{\partial y^2} = \dfrac{1}{y^2} + \dfrac{1}{(1-y)^2}}
#' }
#'
#' **Inverse Hessian:** \eqn{\dfrac{d^2 X}{d Y^2} = \dfrac{2y-1}{y^2(1-y)^2}}
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
    stats::plogis(bounds)
  }

  o$trans_fun <- stats::plogis
  o$trans_inv <- stats::qlogis

  o$trans_abs_jac <- function(y, log = TRUE) {
    # log(J) = - (log(y) + log(1 - y))
    val <- -(log(y) + log1p(-y))
    if (log) val else exp(val)
  }

  o$grad_log_jac <- function(y) {
    # d/dy [-log(y) - log(1-y)] = -1/y + 1/(1-y)
    -1 / y + 1 / (1 - y)
  }

  o$hess_log_jac <- function(y) {
    # d/dy [-1/y + 1/(1-y)] = 1/y^2 + 1/(1-y)^2
    1 / y^2 + 1 / (1 - y)^2
  }

  o$trans_inv_hessian <- function(y) {
    # d2/dy2 [log(y/(1-y))] = (2y - 1) / (y^2 * (1-y)^2)
    (2 * y - 1) / (y^2 * (1 - y)^2)
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
#' @param a Numeric. Scale parameter. Defaults to 1.
#'
#' @details
#' **Transformation (Inverse Softplus):** \eqn{Y = \dfrac{1}{a} \log(e^{aX} - 1)}
#'
#' **Inverse (Softplus):** \eqn{X = \dfrac{1}{a} \log(1 + e^{aY})}
#'
#' **Jacobian:** \eqn{J = \text{sigmoid}(aY)}
#'
#' **Log-Jacobian Derivatives:**
#' \itemize{
#'   \item Gradient: \eqn{\dfrac{\partial \log|J|}{\partial y} = a(1 - \text{plogis}(ay))}
#'   \item Hessian: \eqn{\dfrac{\partial^2 \log|J|}{\partial y^2} = -a^2 \text{dlogis}(ay)}
#' }
#'
#' **Inverse Hessian:** \eqn{\dfrac{d^2 X}{d Y^2} = a \cdot \text{dlogis}(ay)}
#'
#' **Support:** Requires \eqn{X > 0}. Maps to \eqn{Y \in \mathbb{R}}.
#'
#' @return A list of class \code{"transformer"}.
#' @export
softplus_transform <- function(a = 1) {
  # Transformation: Y = 1/a * log(exp(a*X) - 1)  [Inverse Softplus]
  # Inverse:        X = 1/a * log(1 + exp(a*Y))  [Softplus]
  # Jacobian:       dX/dY = sigmoid(a*Y) = 1 / (1 + exp(-a*Y))

  if (a <= 0) {
    stop("Scale parameter 'a' must be greater than 0")
  }

  o <- list()
  o$name <- paste0("softplus_a", round(a, 5))

  o$valid_support <- function(bounds) {
    bounds[1] >= 0
  }

  o$bounds_fun <- function(bounds) {
    res <- log(expm1(a * bounds)) / a
    res[bounds == 0] <- -Inf
    res
  }

  o$trans_fun <- function(x) {
    log(expm1(a * x)) / a
  }

  o$trans_inv <- function(y) {
    pmax(0, y) + log1p(exp(-abs(a * y))) / a
  }

  o$trans_abs_jac <- function(y, log = TRUE) {
    stats::plogis(a * y, log.p = log)
  }

  o$grad_log_jac <- function(y) {
    # d/dy [log(plogis(a*y))] = a * (1 - plogis(a*y))
    a * (1 - stats::plogis(a * y))
  }

  o$hess_log_jac <- function(y) {
    # d/dy [a * (1 - plogis(a*y))] = -a^2 * dlogis(a*y)
    -(a^2) * stats::dlogis(a * y)
  }

  o$trans_inv_hessian <- function(y) {
    # d/dy [plogis(a*y)] = a * dlogis(a*y)
    a * stats::dlogis(a * y)
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
#' **Derivatives with respect to \eqn{y}:**
#' While the derivatives with respect to the parameters \eqn{\theta} remain unchanged (evaluated at \eqn{g^{-1}(y)}),
#' the derivatives with respect to the observed variable \eqn{y} must account for the transformation
#' via the chain rule. Let \eqn{x = g^{-1}(y)} and \eqn{J = \dfrac{dx}{dy}}:
#' \itemize{
#'   \item **Gradient:** \deqn{\dfrac{\partial \ell_Y}{\partial y} = \dfrac{\partial \ell_X}{\partial x} \cdot J + \dfrac{\partial \log |J|}{\partial y}}
#'   \item **Hessian:** \deqn{\dfrac{\partial^2 \ell_Y}{\partial y^2} = \dfrac{\partial^2 \ell_X}{\partial x^2} \cdot J^2 + \dfrac{\partial \ell_X}{\partial x} \cdot \dfrac{d J}{d y} + \dfrac{\partial^2 \log |J|}{\partial y^2}}
#' }
#' where \eqn{\dfrac{d J}{d y}} is the second derivative of the inverse transformation (\code{trans_inv_hessian}).
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
  old_grad_y <- distrib$grad_y
  old_hess_y <- distrib$hess_y

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

  d$grad_y <- function(y, theta) {
    x <- tr$trans_inv(y)
    J <- tr$trans_abs_jac(y, log = FALSE)
    old_grad_y(x, theta) * J + tr$grad_log_jac(y)
  }

  d$hess_y <- function(y, theta) {
    x <- tr$trans_inv(y)
    J <- tr$trans_abs_jac(y, log = FALSE)
    d2x_dy2 <- tr$trans_inv_hessian(y)

    old_hess_y(x, theta) * (J^2) +
      old_grad_y(x, theta) * d2x_dy2 +
      tr$hess_log_jac(y)
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
