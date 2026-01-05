#' Print Method for `distrib` Objects
#'
#' @param x An object of class \code{"distrib"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return The object \code{x} invisibly.
#' @export
#' @method print distrib
print.distrib <- function(x, ...) {
  cat(sprintf("Distribution: %s\n", gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", x$distrib_name, perl = TRUE)))
  cat(sprintf("Type:         %s\n", paste0(toupper(substring(x$type, 1, 1)), substring(x$type, 2))))
  cat(sprintf("Dimensions:   %d\n", x$dimension))


  cat("\nParameters:\n")

  max_param_len <- max(nchar(x$params))

  for (i in seq_len(x$n_params)) {
    p_name <- x$params[i]

    interpretation <- if (!is.null(x$params_interpretation)) {
      paste0("(", x$params_interpretation[[p_name]], ")")
    } else {
      ""
    }

    link_obj <- x$link_params[[i]]
    link_name <- link_obj$link_name
    bounds <- x$params_bounds[[i]]
    domain_str <- paste0("(", bounds[1], ", ", bounds[2], ")")

    cat(sprintf(
      "  %-*s %-20s | Link: %-10s | Domain: %s\n",
      max_param_len, p_name,
      interpretation,
      link_name,
      domain_str
    ))
  }
  invisible(x)
}




#' Plot Method for `distrib` Objects
#'
#' @description
#' Visualizes the Probability Density Function (PDF) or Probability Mass Function (PMF)
#' of a distribution object for a specific set of parameters.
#'
#' @param x An object of class \code{"distrib"}.
#' @param theta A named list or vector of parameters. Must match \code{x$params}.
#' @param xlim Optional numeric vector of length 2 indicating the x-axis range.
#'   If \code{NULL} (default), the range is automatically calculated.
#' @param ... Additional arguments passed to the base \code{\link{plot}} function.
#' @importFrom graphics points segments
#' @export
#' @method plot distrib
plot.distrib <- function(x, theta, xlim = NULL, ...) {
  if (missing(theta)) {
    stop("Argument 'theta' is missing.")
  }

  if (is.numeric(theta) && !is.list(theta)) {
    theta <- as.list(theta)
  }

  if (!all(x$params %in% names(theta))) {
    stop(paste("Missing parameters in 'theta'. Expected:", paste(x$params, collapse = ", ")))
  }

  if (is.null(xlim)) {
    lower_q <- x$quantile(0.01, theta)
    upper_q <- x$quantile(0.99, theta)

    if (x$type == "discrete") {
      lower_q <- floor(lower_q)
      upper_q <- ceiling(upper_q)
      if (lower_q > x$bounds[1]) lower_q <- lower_q - 1
      if (upper_q < x$bounds[2]) upper_q <- upper_q + 1
    } else {
      span <- max(upper_q - lower_q, 1)
      lower_q <- lower_q - 0.1 * span
      upper_q <- upper_q + 0.1 * span
    }

    xlim <- c(
      max(lower_q, x$bounds[1]),
      min(upper_q, if (is.infinite(x$bounds[2])) Inf else x$bounds[2])
    )
  }

  if (x$type == "discrete") {
    seq_x <- seq(floor(xlim[1]), ceiling(xlim[2]), by = 1)
    seq_x <- seq_x[seq_x >= x$bounds[1] & seq_x <= x$bounds[2]]
  } else {
    seq_x <- seq(xlim[1], xlim[2], length.out = 300)
  }

  dens_y <- x$pdf(seq_x, theta, log = FALSE)

  dots <- list(...)
  if (is.null(dots$main)) {
    param_str <- paste(names(theta), round(unlist(theta), 2), sep = "=", collapse = ", ")
    dots$main <- paste0(tools::toTitleCase(x$distrib_name), " distribution\n(", param_str, ")")
  }
  if (is.null(dots$xlab)) dots$xlab <- "y"
  if (is.null(dots$ylab)) dots$ylab <- if (x$type == "discrete") "Probability (PMF)" else "Density (PDF)"
  if (is.null(dots$col)) dots$col <- "black"
  if (is.null(dots$lwd)) dots$lwd <- 2

  plot_args <- c(list(x = seq_x, y = dens_y, xlim = xlim), dots)

  if (x$type == "discrete") {
    plot_args$type <- "n"
    do.call(plot, plot_args)
    segments(x0 = seq_x, y0 = 0, x1 = seq_x, y1 = dens_y, col = dots$col, lwd = dots$lwd)
    points(seq_x, dens_y, pch = 16, col = dots$col)
  } else {
    plot_args$type <- "l"
    do.call(plot, plot_args)
  }

  invisible(x)
}
