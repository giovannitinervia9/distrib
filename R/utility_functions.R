#' Check Analytical Derivatives against Numerical Approximations
#'
#' @description
#' This utility function verifies the accuracy of the analytical gradient and Hessian
#' implementations within a \code{"distrib"} object.
#' It uses **Relative Error** for comparison.
#'
#' @param distrib An object of class \code{"distrib"}.
#' @param n Integer. The number of random parameter configurations and observations to generate for testing. Defaults to 10.
#'
#' @details
#' The validation compares analytical results vs numerical approximations (via finite differences).
#' The error metric used is the **Relative Error**:
#' \deqn{Err = \frac{|x_{ana} - x_{num}|}{\max(|x_{ana}|, 10^{-10})}}
#' The denominator floor (\eqn{10^{-10}}) prevents division by zero when the analytical derivative is exactly 0.
#'
#' @return Invisibly returns a named numeric vector containing the maximum relative errors.
#'
#' @importFrom stats runif rnorm
#' @export
check_derivatives_distrib <- function(distrib, n = 10) {
  params <- distrib$params
  n_params <- distrib$n_params
  params_bounds <- distrib$params_bounds

  theta <- vector("list", n_params)
  names(theta) <- params

  for (i in 1:n_params) {
    pname <- params[i]
    pbound <- params_bounds[[pname]]

    if (all(is.infinite(pbound))) {
      theta[[i]] <- rnorm(n)
    } else if (!is.infinite(pbound[1]) & is.infinite(pbound[2])) {
      theta[[i]] <- runif(n, pbound[1] + 0.5, pbound[1] + 5)
    } else if (is.infinite(pbound[1]) & !is.infinite(pbound[2])) {
      theta[[i]] <- runif(n, pbound[2] - 5, pbound[2] - 0.5)
    } else {
      width <- pbound[2] - pbound[1]
      margin <- width * 0.1
      theta[[i]] <- runif(n, pbound[1] + margin, pbound[2] - margin)
    }
  }

  y <- distrib$rng(n, theta)

  ana_grad_list <- distrib$gradient(y, theta)
  ana_hess_list <- distrib$hessian(y, theta, expected = FALSE)

  grad_errs <- numeric(n)
  hess_errs <- numeric(n)


  calc_rel_error <- function(ana, num) {
    denominator <- max(abs(ana), 1e-10)
    abs(ana - num) / denominator
  }

  for (i in 1:n) {
    theta_vec <- numeric(n_params)
    names(theta_vec) <- params
    for (p in params) {
      theta_vec[p] <- theta[[p]][i]
    }

    lpdf <- function(x) {
      t_list <- as.list(x)
      names(t_list) <- params
      val <- tryCatch(distrib$loglik(y[i], t_list), warning = function(w) NaN)
      if (is.nan(val) || is.infinite(val)) {
        return(-1e100)
      }
      return(val)
    }

    num_g <- numDeriv::grad(func = lpdf, x = theta_vec)
    ana_g <- numeric(n_params)
    for (k in 1:n_params) {
      p_name <- params[k]
      ana_g[k] <- ana_grad_list[[p_name]][i]
    }

    errs_g <- calc_rel_error(ana_g, num_g)

    if (any(is.na(errs_g))) {
      grad_errs[i] <- Inf
    } else {
      grad_errs[i] <- max(errs_g)
    }

    num_h <- numDeriv::hessian(func = lpdf, x = theta_vec)
    hess_names <- names(ana_hess_list)
    current_hess_err <- 0

    for (h_name in hess_names) {
      val_ana <- ana_hess_list[[h_name]][i]

      parts <- strsplit(h_name, "_")[[1]]
      idx_row <- match(parts[1], params)
      if (length(parts) == 2) {
        idx_col <- match(parts[2], params)
      } else {
        idx_col <- match(parts[length(parts)], params)
      }

      val_num <- num_h[idx_row, idx_col]

      if (is.na(val_num) || is.nan(val_num)) {
        err_val <- Inf
      } else {
        err_val <- calc_rel_error(val_ana, val_num)
      }

      if (!is.na(err_val) && err_val > current_hess_err) {
        current_hess_err <- err_val
      }
    }
    hess_errs[i] <- current_hess_err
  }

  max_grad_err <- max(grad_errs, na.rm = TRUE)
  max_hess_err <- max(hess_errs, na.rm = TRUE)

  cat("----------------------------------------------------\n")
  cat("Distribution:", distrib$distrib_name, "\n")
  cat("Max Relative Error (Grad): ", formatC(max_grad_err, format = "e", digits = 5), "\n")
  cat("Max Relative Error (Hess): ", formatC(max_hess_err, format = "e", digits = 5), "\n")

  threshold <- 1e-4

  if (!is.infinite(max_grad_err) && !is.infinite(max_hess_err) &&
    max_grad_err < threshold && max_hess_err < threshold) {
    cat("\n\033[32m[OK] Analytical derivatives match numerical ones.\033[0m\n")
  } else {
    cat("\n\033[31m[WARNING] Discrepancies detected. Check formulas.\033[0m\n")
  }
  cat("----------------------------------------------------\n")

  invisible(c(grad_error = max_grad_err, hess_error = max_hess_err))
}





#' Check Expected Hessian via Monte Carlo Simulation
#'
#' @description
#' Verifies the Expected Hessian implementation by comparing it with the
#' mean of the Observed Hessians computed over a large sample.
#' Uses Absolute Error when the expected value is close to zero, and Relative Error otherwise.
#'
#' @param distrib An object of class \code{"distrib"}.
#' @param n_sim Integer. Number of Monte Carlo samples (default 50,000).
#'   Higher values reduce approximation error.
#' @param theta A named \code{list} containing the parameter values to use for the check.
#'   Each element must be a numeric vector of length 1 (scalar) with a value within
#'   the distribution's valid domain. The names must match the distribution's parameter
#'   names (e.g., \code{list(mu = 0.5, sigma = 1)}). If missing, appropriate values
#'   are generated automatically based on the parameter bounds.
#'
#' @export
check_expected_hessian_distrib <- function(distrib, n_sim = 50000, theta) {
  params <- distrib$params
  bounds <- distrib$params_bounds

  if (missing(theta)) {
    theta <- vector("list", length(params))
    names(theta) <- params

    for (p in params) {
      b <- bounds[[p]]
      if (all(is.infinite(b))) {
        theta[[p]] <- rnorm(1)
      } else if (b[1] == 0 && b[2] == 1) {
        theta[[p]] <- runif(1, 0.2, 0.8)
      } else if (b[1] == 0 && is.infinite(b[2])) {
        theta[[p]] <- runif(1, 1.0, 5.0)
      } else {
        mid <- if (is.infinite(b[1])) b[2] - 2 else if (is.infinite(b[2])) b[1] + 2 else mean(b)
        theta[[p]] <- runif(1, mid - 0.1, mid + 0.1)
      }
    }
  }

  cat("----------------------------------------------------\n")
  cat("Distribution:", distrib$distrib_name, "\n")
  cat("Parameters:  ", paste(names(theta), round(unlist(theta), 4), sep = "=", collapse = ", "), "\n")
  cat("Samples:     ", n_sim, "\n\n")

  y_sim <- distrib$rng(n_sim, theta)

  H_expected_ana <- distrib$hessian(y_sim[1], theta, expected = TRUE)
  H_observed_mean <- lapply(distrib$hessian(y_sim, theta, expected = FALSE), mean)

  max_err <- 0
  for (name in names(H_expected_ana)) {
    val_ana <- H_expected_ana[[name]]
    val_mc <- H_observed_mean[[name]]
    abs_diff <- abs(val_ana - val_mc)

    is_zero <- abs(val_ana) < 1e-8
    err_val <- if (is_zero) abs_diff else abs_diff / abs(val_ana)
    err_type <- if (is_zero) "AbsErr" else "RelErr"

    max_err <- max(max_err, err_val)

    cat(sprintf(
      "  %-15s | Exp: %9.5f | Mean Obs: %9.5f | %s: %.2e\n",
      name, val_ana, val_mc, err_type, err_val
    ))
  }

  if (max_err < 0.05) {
    cat("\n\033[32m[OK] Expected Hessian aligns with Mean Observed Hessian.\033[0m\n")
  } else {
    cat("\n\033[31m[WARNING] Large discrepancy detected.\033[0m\n")
  }
  cat("----------------------------------------------------\n")
}
