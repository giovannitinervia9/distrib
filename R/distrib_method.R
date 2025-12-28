#' Print Method for `distrib` Objects
#'
#' @param x An object of class \code{"distrib"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return The object \code{x} invisibly.
#' @export
#' @method print distrib
print.distrib <- function(x, ...) {
  cat(sprintf("Distribution: %s\n", paste0(toupper(substring(x$distrib_name, 1, 1)), substring(x$distrib_name, 2))))
  cat(sprintf("Type:         %s\n", paste0(toupper(substring(x$type, 1, 1)), substring(x$type, 2))))
  cat(sprintf("Dimensions:   %d\n", x$dimension))

  cat("\nParameters:\n")
  max_param_len <- max(nchar(x$params))
  for (i in seq_len(x$n_params)) {
    link_obj <- x$link_params[[i]]
    link_name <- link_obj$link_name
    bounds <- x$params_bounds[[i]]
    domain_str <- paste0(bounds[1], ", ", bounds[2])
    cat(sprintf(
      "% -*s  |  Link: %-10s  |  Domain: %s\n",
      max_param_len, x$params[i],
      link_name,
      domain_str
    ))
  }
  invisible(x)
}