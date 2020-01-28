#' Inverse ALR for a vector
#' @param v a compositional vector (sum_i(v_i)=1)
#' @examples inverse_alr_vec(c(4, 0.5, 1))
#'
inverse_alr_vec <- function(v){
  ## v is a vector
  exp_v = exp(v)
  return(exp(v)/sum(exp_v))
}

#' Inverse ALR for a matrix
#' @param m a matrix in which rows are compositional vectors (rowSums(m_i) for all rows i)
#' @examples inverse_alr_mat(as(compositions::alr(MCMCpack::rdirichlet(1000, rep(1, 5))), 'matrix'))
inverse_alr_mat <- function(m){
  return(t(apply(m, 1, inverse_alr_vec)))
}

#' ALR (additive log-ratio) transformation
#' @param v a compositional vector (sum_i(v_i)=1)
#' @param idx index of the part of the composition which is used as baseline.
#' @examples alr_vec(c(0.2, 0.5, 0.3))
alr_vec <- function(v, idx=length(v)){
  return(log(v/v[idx]))
}

#' ALR (additive log-ratio) transformation of a matrix
#' @param m a matrix in which rows are compositional vectors (rowSums(m_i) for all rows i)
#' @param idx index of the part of the composition which is used as baseline.
#' @examples alr_mat(MCMCpack::rdirichlet(1000, rep(1, 5)))
alr_mat <- function(m, idx=ncol(m)){
  return(sweep(m, 1, m[,idx], '/'))
}

#' Fit compositional data on a logistic normal model
#' @param x a matrix of compositional instances in a shape of a matrix in which rows are compositional vectors (rowSums(m_i) for all rows i)
#' @examples fit_logistic_normal(MCMCpack::rdirichlet(1000, rep(1, 5)))
fit_logistic_normal <- function(x){
  alr_x = alr_mat(x)
  mu = colMeans(alr_x)
  covar = cov(alr_x)
  return(list(mu=mu, covariance=covar))
}

#' @param n number of instances to simulate
#' @param mu parameter mu of the logistic normal distribution, It's a vector of length n
#' @param sigma parameter sigma of the logistic normal distribution. It's a covariance matrix of dimensions n x n
rlogistic_normal <- function(n, mu, sigma){
  require(MASS)
  samples_invalr = mvrnorm(n = n, mu = mu, Sigma = sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  return(inverse_alr_mat(samples_invalr))
}

