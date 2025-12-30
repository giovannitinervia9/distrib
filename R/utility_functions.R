#' Check Analytical Derivatives against Numerical Approximations
#'
#' @description
#' This utility function verifies the accuracy of the analytical gradient and Hessian
#' implementations within a \code{"distrib"} object.
#' It compares analytical calculations against numerical approximations (via finite differences)
#' and reports the maximum relative error **for each parameter individually**.
#'
#' @param distrib An object of class \code{"distrib"}.
#' @param n Integer. The number of random parameter configurations and observations to generate for testing. Defaults to 20.
#'
#' @details
#' The validation generates \code{n} random parameter sets and observations based on the
#' bounds defined in the distribution object. It then calculates:
#' \itemize{
#'   \item **Gradient:** The analytical gradient vs numeric gradient (`numDeriv::grad`).
#'   \item **Hessian:** The analytical observed Hessian vs numeric Hessian (`numDeriv::hessian`).
#' }
#'
#' The error metric used is the **Relative Error** (robust):
#' \deqn{Err = \frac{|x_{ana} - x_{num}|}{\max(|x_{ana}|, 10^{-10})}}
#'
#' A detailed report is printed to the console showing the worst-case scenario (Maximum Error)
#' found for each parameter across the `n` samples.
#'
#' @return Invisibly returns a list containing the error statistics.
#'
#' @importFrom stats runif rnorm
#' @importFrom numDeriv grad hessian
#' @export
check_derivatives_distrib <- function(distrib, n = 20) {
  params <- distrib$params
  n_params <- distrib$n_params
  params_bounds <- distrib$params_bounds

  # --- 1. Generazione Parametri Casuali ---
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

  # --- 2. Generazione Dati (y) e Calcolo Analitico ---
  y <- distrib$rng(n, theta)

  # Analitico (Vettorizzato)
  ana_grad_list <- distrib$gradient(y, theta)
  ana_hess_list <- distrib$hessian(y, theta, expected = FALSE)
  hess_names <- names(ana_hess_list)

  # --- 3. Strutture per tracciare il "Worst Case" ---
  # Usiamo data.frame per tenere traccia di Valore Ana, Valore Num e Errore Max

  stats_grad <- data.frame(
    name = params, max_err = -1, worst_ana = NA, worst_num = NA,
    stringsAsFactors = FALSE, row.names = params
  )

  stats_hess <- data.frame(
    name = hess_names, max_err = -1, worst_ana = NA, worst_num = NA,
    stringsAsFactors = FALSE, row.names = hess_names
  )

  calc_rel_error <- function(ana, num) {
    denominator <- max(abs(ana), 1e-10)
    abs(ana - num) / denominator
  }

  # --- 4. Loop di confronto Numerico ---
  for (i in 1:n) {
    # Vettore theta corrente
    theta_vec <- numeric(n_params)
    names(theta_vec) <- params
    for (p in params) theta_vec[p] <- theta[[p]][i]

    # Wrapper per numDeriv
    lpdf <- function(x) {
      t_list <- as.list(x)
      names(t_list) <- params
      val <- tryCatch(distrib$loglik(y[i], t_list), warning = function(w) NaN)
      if (is.nan(val) || is.infinite(val)) {
        return(-1e100)
      }
      return(val)
    }

    # Gradient Check
    num_g <- numDeriv::grad(func = lpdf, x = theta_vec)
    for (k in 1:n_params) {
      p_name <- params[k]
      val_ana <- ana_grad_list[[p_name]][i]
      val_num <- num_g[k]
      err <- calc_rel_error(val_ana, val_num)

      if (err > stats_grad[p_name, "max_err"]) {
        stats_grad[p_name, "max_err"] <- err
        stats_grad[p_name, "worst_ana"] <- val_ana
        stats_grad[p_name, "worst_num"] <- val_num
      }
    }

    # Hessian Check
    num_h <- numDeriv::hessian(func = lpdf, x = theta_vec)
    for (h_name in hess_names) {
      val_ana <- ana_hess_list[[h_name]][i]

      parts <- strsplit(h_name, "_")[[1]]
      idx_row <- match(parts[1], params)
      idx_col <- match(parts[length(parts)], params)
      val_num <- num_h[idx_row, idx_col]

      if (is.na(val_num) || is.nan(val_num)) {
        err <- Inf
      } else {
        err <- calc_rel_error(val_ana, val_num)
      }

      if (err > stats_hess[h_name, "max_err"]) {
        stats_hess[h_name, "max_err"] <- err
        stats_hess[h_name, "worst_ana"] <- val_ana
        stats_hess[h_name, "worst_num"] <- val_num
      }
    }
  }

  # --- 5. Stampa Report ---
  global_max_grad <- max(stats_grad$max_err)
  global_max_hess <- max(stats_hess$max_err)

  cat("----------------------------------------------------\n")
  cat("Distribution:", distrib$distrib_name, "\n")
  cat("Samples tested:", n, "\n")
  cat("(Displaying worst-case discrepancy found across samples)\n\n")

  cat(">>> Gradient Checks <<<\n")
  for (r in 1:nrow(stats_grad)) {
    cat(sprintf(
      "  %-10s | Ana: %10.4e | Num: %10.4e | Max RelErr: %.3e\n",
      stats_grad$name[r], stats_grad$worst_ana[r], stats_grad$worst_num[r], stats_grad$max_err[r]
    ))
  }
  cat("\n")

  cat(">>> Hessian Checks <<<\n")
  for (r in 1:nrow(stats_hess)) {
    cat(sprintf(
      "  %-12s | Ana: %10.4e | Num: %10.4e | Max RelErr: %.3e\n",
      stats_hess$name[r], stats_hess$worst_ana[r], stats_hess$worst_num[r], stats_hess$max_err[r]
    ))
  }

  cat("----------------------------------------------------\n")

  threshold <- 1e-4
  if (global_max_grad < threshold && global_max_hess < threshold) {
    message("[OK] Analytical derivatives match numerical ones.")
  } else {
    warning("[WARNING] Discrepancies detected. See detailed table above.")
  }
  cat("----------------------------------------------------\n")

  invisible(list(grad = stats_grad, hess = stats_hess))
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
    message("\n[OK] Expected Hessian aligns with Mean Observed Hessian.\n")
  } else {
    warning("\n[WARNING] Large discrepancy detected.\n")
  }
  cat("----------------------------------------------------\n")
}





#' Check Consistency of Parameter Dimensions
#'
#' Validates that all elements in the provided parameter list have compatible lengths.
#' Each parameter must have a length of either 1 (scalar) or exactly equal to
#' `max_length`. This ensures safe vector recycling and dimensional consistency.
#'
#' @param theta A named list of vectors (parameters). Each element represents
#'   a parameter of a distribution (e.g., `mu`, `sigma`).
#' @param max_length (Optional) An integer specifying the required maximum length.
#'   If not provided, it defaults to the maximum length found among the elements
#'   of `theta`. Providing this argument allows validation against an external
#'   dimension (e.g., sample size `n`).
#'
#' @return Returns `NULL` invisibly if the check passes.
#'
#' @section Errors:
#' The function throws an error (`stop`) if it detects any parameter with a length
#' that is neither 1 nor `max_length`. The error message lists the specific parameters
#' causing the mismatch.
#'
#' @examples
#' # --- Case 1: Implicit max length ---
#' # Valid: all scalars
#' check_params_dim(list(mu = 1, sigma = 2))
#'
#' # Valid: mixing scalar and vector
#' check_params_dim(list(mu = 1:5, sigma = 1))
#'
#' # Invalid: incompatible lengths (2 vs 3)
#' \dontrun{
#' check_params_dim(list(mu = 1:2, sigma = 1:3))
#' }
#'
#' # --- Case 2: Explicit max_length ---
#' # Valid: vector matches max_length (5)
#' check_params_dim(list(mu = 1:5, sigma = 1), max_length = 5)
#'
#' # Invalid: vector length (3) does not match required max_length (5)
#' # This is useful to enforce consistency with a dataset size n = 5
#' \dontrun{
#' check_params_dim(list(mu = 1:3, sigma = 1), max_length = 5)
#' }
#'
#' @export
check_params_dim <- function(theta, max_length) {
  len_theta <- lengths(theta)

  if (missing(max_length)) {
    max_length <- max(len_theta)
  }

  # Check: length must be 1 OR exactly max_length
  mismatch_idx <- which(len_theta != 1 & len_theta != max_length)

  if (length(mismatch_idx) > 0) {
    bad_params <- names(theta)[mismatch_idx]
    bad_lens <- len_theta[mismatch_idx]

    error_msg <- paste0(
      "Parameter dimension mismatch. All parameters should have length 1 or ", max_length, ".\n",
      "Found mismatches in:\n",
      paste(paste0("  - ", bad_params, ": length ", bad_lens), collapse = "\n")
    )

    stop(error_msg, call. = FALSE)
  }

  invisible(NULL)
}
