#' S3 Class `distrib`
#'
#' @name distrib-class
#' @aliases distrib
#'
#' @title S3 Class for Probability Distributions
#'
#' @description
#' A generic S3 object of class \code{"distrib"} that encapsulates all mathematical properties
#' of a probability distribution, including statistical functions (PDF, CDF, Quantile),
#' derivatives (Gradient, Hessian) for maximum likelihood estimation, and analytical moments.
#'
#' @details
#' Objects of class \code{"distrib"} are typically created by constructor functions such as
#' \code{\link{gaussian_distrib}}, \code{\link{gamma_distrib}}, \code{\link{poisson_distrib}}, etc.
#'
#' @section Structure:
#' An object of class \code{"distrib"} is a \code{list} with the following components:
#'
#' \describe{
#'   \item{\code{distrib_name}}{A character string indicating the name of the distribution (e.g., "gaussian", "gamma").}
#'   \item{\code{type}}{A character string indicating the type of distribution: \code{"continuous"} or \code{"discrete"}.}
#'   \item{\code{dimension}}{An integer indicating the dimension of the random variable (typically 1 for univariate distributions).}
#'   \item{\code{bounds}}{A numeric vector of length 2 indicating the support of the distribution (e.g., \code{c(0, Inf)}).}
#'
#'   \item{\code{params}}{A character vector containing the names of the parameters (e.g., \code{c("mu", "sigma")}).}
#'   \item{\code{n_params}}{An integer indicating the number of parameters.}
#'   \item{\code{params_interpretation}}{A named character vector providing a human-readable interpretation of the parameters (e.g., \code{c(mu="mean", sigma="std. dev.")}).}
#'   \item{\code{params_bounds}}{A named list containing the valid domain for each parameter.}
#'   \item{\code{link_params}}{A named list of link function objects (one for each parameter).}
#'
#'   \item{\code{pdf}}{A function \code{pdf(y, theta, log = FALSE)} that computes the probability density (or mass) function. \code{theta} must be a named list of parameters.}
#'   \item{\code{cdf}}{A function \code{cdf(q, theta, lower.tail = TRUE, log.p = FALSE)} that computes the cumulative distribution function.}
#'   \item{\code{quantile}}{A function \code{quantile(p, theta, lower.tail = TRUE, log.p = FALSE)} that computes the quantile function (inverse CDF).}
#'   \item{\code{rng}}{A function \code{rng(n, theta)} that generates \code{n} random draws from the distribution.}
#'
#'   \item{\code{loglik}}{A function \code{loglik(y, theta)} that computes the log-likelihood (wrapper for \code{pdf(..., log=TRUE)}).}
#'   \item{\code{gradient}}{A function \code{gradient(y, theta)} that computes the score vector (gradients of the log-pdf with respect to the parameters). Returns a named list.}
#'   \item{\code{hessian}}{A function \code{hessian(y, theta, expected = FALSE)} that computes the Hessian matrix of the log-pdf. If \code{expected = TRUE}, it returns the Fisher Information Matrix (with negative sign).}
#'
#'   \item{\code{mean}}{A function \code{mean(theta)} returning the expected value \eqn{\mathbb{E}[Y]}.}
#'   \item{\code{variance}}{A function \code{variance(theta)} returning the variance \eqn{\mathbb{V}[Y]}.}
#'   \item{\code{skewness}}{A function \code{skewness(theta)} returning the skewness \eqn{\gamma_1}.}
#'   \item{\code{kurtosis}}{A function \code{kurtosis(theta)} returning the excess kurtosis \eqn{\gamma_2}.}
#'
#'   \item{\code{kernel}}{A function \code{kernel(y, theta)} computing the unnormalized density.}
#'   \item{\code{normalization_constant}}{A function \code{normalization_constant(y, theta)} returning the normalizing term ensuring the density integrates to 1.}
#' }
#'
#' @seealso
#' \code{\link{gaussian_distrib}}, \code{\link{gamma_distrib}}, \code{\link{poisson_distrib}},
#' \code{\link{binomial_distrib}}, \code{\link{bernoulli_distrib}}
NULL
