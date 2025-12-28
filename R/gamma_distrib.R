gamma_distrib <- function(link_mu = log_link(), link_sigma2 = log_link()) {
  o <- list()
  class(o) <- c("distrib")
  
  o$distrib_name <- "gamma"
  o$type <- "continuous"
  o$dimension <- 1
  o$bounds <- c(0, Inf)
  
  o$params <- c("mu", "sigma2")
  o$n_params <- 2
  o$params_bounds <- list(
    mu = c(0, Inf),
    sigma2 = c(0, Inf)
  )
  o$link_params <- list(
    mu = link_mu,
    sigma2 = link_sigma
  )
  
  o$pdf <- function(y, theta, log = FALSE) {
    mu <- theta[["mu"]]
    sigma2 <- theta[["sigma2"]]
    stats::dgamma(
      x = y,
      shape = mu^2/sigma2,
      rate = mu/sigma2,
      log = log
    )
  }
  
  o$cdf <- function(q, theta, lower.tail = TRUE, log.p = FALSE) {
    mu <- theta[["mu"]]
    sigma2 <- theta[["sigma2"]]
    stats::pgamma(
      q = q,
      shape = mu^2/sigma2,
      rate = mu/sigma2,
      lower.tail = lower.tail,
      log.p = log.p
    )
  }
  
  o$qf <- function(p, theta, lower.tail = TRUE, log.p = FALSE) {
    mu <- theta[["mu"]]
    sigma2 <- theta[["sigma2"]]
    stats::qgamma(
      p = p,
      shape = mu^2/sigma2,
      rate = mu/sigma2,
      lower.tail = lower.tail,
      log.p = log.p
    )
  }
  
  o$rng <- function(n, theta) {
    mu <- theta[["mu"]]
    sigma2 <- theta[["sigma2"]]
    stats::rgamma(
      n = n,
      shape = mu^2/sigma2,
      rate = mu/sigma2
    )
  }
  
  o$loglik <- function(y, theta) {
    o$pdf(y, theta, log = TRUE)
  }
  
  o$gradient <- function(y, theta) {
    mu <- theta[["mu"]]
    sigma2 <- theta[["sigma2"]]
    mu_sigma2 <- mu/sigma2
    mu2_sigma2 <- mu_sigma2*mu
    digamma_mu2_sigma2 <- digamma(mu2_sigma2)
    
    list(
      mu = (-2*mu*digamma_mu2_sigma2 + 2*mu*log(mu_sigma2) + mu + 2*mu*log(y) - y)/sigma2,
      sigma2 = -(mu*(-mu*digamma_mu2_sigma2 + mu + mu*(log(mu_sigma2) + log(y)) - y))/sigma2^2
    )
  }
  
  o$hessian <- function(y, theta, expected = FALSE) {
    mu <- theta[["mu"]]
    sigma2 <- theta[["sigma2"]]
    mu_sigma2 <- mu/sigma2
    mu2_sigma2 <- mu_sigma2*mu
    digamma_mu2_sigma2 <- digamma(mu2_sigma2)
    trigamma_mu2_sigma2 <- trigamma(mu2_sigma2)
    
    list(
      mu_mu = (-(4*mu^2*trigamma_mu2_sigma2)/sigma2 - 2*digamma_mu2_sigma2 + 2*log(mu_sigma2) + 2*log(y) + 3)/sigma2,
      sigma2_sigma2 = -(mu*(2*mu*sigma2*digamma_mu2_sigma2 + mu^3*trigamma_mu2_sigma2 + sigma2*(-2*mu*log(mu_sigma2) - 3*mu - 2*mu*log(y) + 2*y)))/sigma2^4,
      mu_sigma2 = (2*mu*sigma2*digamma_mu2_sigma2 + 2*mu^3*trigamma_mu2_sigma2 + sigma2*(-2*mu*log(mu_sigma2) - 3*mu - 2*mu*log(y) + y))/sigma2^3
    )
    
    
  }
  
  o$kernel <- function(y, theta) {
    mu <- theta[["mu"]]
    sigma2 <- theta[["sigma2"]]
    mu_sigma2 <- mu/sigma2
    mu2_sigma2 <- mu_sigma2*mu
    exp((mu2_sigma2 - 1)*log(y) - y*mu_sigma2)
    
  }
  
  o$normalization_constant <- function(y, theta) {
    mu <- theta[["mu"]]
    sigma2 <- theta[["sigma2"]]
    mu_sigma2 <- mu/sigma2
    mu2_sigma2 <- mu_sigma2*mu
    gamma(mu2_sigma2)/((mu_sigma2)^(mu2_sigma2))
  }
  
  invisible(o)
  
  
}