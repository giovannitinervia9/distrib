#' Calculate Theoretical Variance
#'
#' A generic function to calculate the theoretical variance (second central moment)
#' of a distribution object.
#'
#' @details
#' The variance measures the spread of the distribution around its mean \eqn{\mu}.
#' It is defined as:
#' \deqn{Var(X) = E[(X - \mu)^2]}
#'
#' **Continuous Case:**
#' \deqn{Var(X) = \int_{-\infty}^{\infty} (x - \mu)^2 f(x) \, dx}
#'
#' **Discrete Case:**
#' \deqn{Var(X) = \sum_{x \in S} (x - \mu)^2 p(x)}
#'
#' where \eqn{f(x)} is the probability density function (PDF) and \eqn{p(x)} is the
#' probability mass function (PMF).
#'
#' @param x An object representing a probability distribution (e.g., of class \code{distrib}).
#' @param ... Additional arguments passed to specific methods (e.g., parameters \code{theta}).
#'
#' @return A numeric vector representing the variance.
#' @seealso \code{\link{std_dev}}, \code{\link{moment}}
#' @export
variance <- function(x, ...) {
  UseMethod("variance")
}




#' Calculate Standard Deviation
#'
#' A generic function to calculate the standard deviation of a distribution object.
#'
#' @details
#' The standard deviation is the square root of the variance:
#' \deqn{\sigma = \sqrt{Var(X)} = \sqrt{E[(X - \mu)^2]}}
#'
#' It shares the same units as the random variable \eqn{X}.
#'
#' @param x An object representing a probability distribution.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return A numeric vector representing the standard deviation.
#' @seealso \code{\link{variance}}
#' @export
std_dev <- function(x, ...) {
  UseMethod("std_dev")
}




#' Calculate Skewness
#'
#' A generic function to calculate the skewness (third standardized moment) of a
#' distribution object.
#'
#' @details
#' Skewness measures the asymmetry of the probability distribution about its mean.
#' It is defined as:
#' \deqn{\gamma_1 = E\left[\left(\dfrac{X - \mu}{\sigma}\right)^3\right] = \dfrac{\mu_3}{\sigma^3}}
#'
#' **Continuous Case:**
#' \deqn{\gamma_1 = \dfrac{1}{\sigma^3} \int_{-\infty}^{\infty} (x - \mu)^3 f(x) \, dx}
#'
#' **Discrete Case:**
#' \deqn{\gamma_1 = \dfrac{1}{\sigma^3} \sum_{x \in S} (x - \mu)^3 p(x)}
#'
#' where \eqn{\mu_3} is the third central moment and \eqn{\sigma} is the standard deviation.
#' Negative skew indicates a tail on the left side, positive skew indicates a tail on the right.
#'
#' @param x An object representing a probability distribution.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return A numeric vector representing the skewness coefficient.
#' @export
skewness <- function(x, ...) {
  UseMethod("skewness")
}




#' Calculate Kurtosis
#'
#' A generic function to calculate the kurtosis (fourth standardized moment) of a
#' distribution object.
#'
#' @details
#' Kurtosis measures the "tailedness" of the probability distribution.
#' It is defined as:
#' \deqn{\beta_2 = E\left[\left(\dfrac{X - \mu}{\sigma}\right)^4\right] = \dfrac{\mu_4}{\sigma^4}}
#'
#' **Continuous Case:**
#' \deqn{\beta_2 = \dfrac{1}{\sigma^4} \int_{-\infty}^{\infty} (x - \mu)^4 f(x) \, dx}
#'
#' **Discrete Case:**
#' \deqn{\beta_2 = \dfrac{1}{\sigma^4} \sum_{x \in S} (x - \mu)^4 p(x)}
#'
#' where \eqn{\mu_4} is the fourth central moment and \eqn{\sigma} is the standard deviation.
#' Note that this function typically returns the raw kurtosis (Normal distribution = 3).
#' Use \code{kurtosis - 3} to obtain the *excess kurtosis*.
#'
#' @param x An object representing a probability distribution.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return A numeric vector representing the kurtosis coefficient.
#' @export
kurtosis <- function(x, ...) {
  UseMethod("kurtosis")
}




#' Cumulative Distribution Function
#'
#' @description
#' A generic function to compute the Cumulative Distribution Function (CDF)
#' for a given distribution object.
#'
#' @param x An object representing a probability distribution.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return A vector of probabilities corresponding to the calculated CDF.
#' @export
cdf <- function(x, ...) {
  UseMethod("cdf")
}
