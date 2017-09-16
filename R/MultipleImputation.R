ParafacMI <- function(Y, R, chain, num_impute = 20, method = "CCMV", j_0 = NULL, xi = NULL) {
  out <- list()
  iters <- floor(seq(from = 1, to = length(chain$loglik), length = num_impute))
  
  if(is.null(j_0) & method == "NIP") stop("Must include time to impute data for NIP")
  if(is.null(j_0) & method == "TLO") stop("Must include time to impute data for TLO")
  if(!is.null(chain$gamma) & method == "MAR") stop("This looks like a chain for MNAR, but using MAR imputation; fit the model under MAR to use this")
  if(is.null(xi)) xi <- rep(0, num_impute)
  stopifnot(method %in% c("CCMV", "TLO", "NIP", "MAR"))
  stopifnot(is.numeric(xi)); stopifnot(length(xi) == num_impute)
  
  for(n in 1:num_impute) {
    i <- iters[n]
    omega <- chain$omega[i, ]
    log_omega <- log(omega)
    gamma <- chain$gamma[i,,]
    if(method != "MAR") {
      log_gamma <- log(gamma)
      log_1_m_gamma <- log(1 - gamma) 
    }
    beta <- chain$beta[i,,]
    log_beta <- log(beta)
    log_1_m_beta <- log(1 - beta)
    
    if(method == "CCMV") {
      out[[n]] <- CCMVMI(Y, R, omega, log_omega, gamma, log_gamma, log_1_m_gamma, beta, log_beta, log_1_m_beta)
    } else if(method == "NIP") {
      out[[n]] <- NIPMI(Y, R, omega, log_omega, gamma, log_gamma, log_1_m_gamma, beta, log_beta, log_1_m_beta, j_0)
    } else if(method == "TLO") {
      out[[n]] <- TLOMI(Y, R, omega, log_omega, gamma, log_gamma, log_1_m_gamma, beta, log_beta, log_1_m_beta, xi, j_0)
    } else {
      out[[n]] <- MARMI(Y = Y, R = R, omega = omega, log_omega = log_omega, beta = beta, log_beta = log_beta, log_1_m_beta = log_1_m_beta)
    }
  }
  return(out)
}