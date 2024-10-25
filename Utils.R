###########################
###       Utils         ###
###########################

#####################################################
###       Convert prices into log-returns         ###
#####################################################
price_to_returns <- function(x){
  r <- c(NA, diff(log(x)))
  return(r)
}


#####################################################
###           Simulate Datasets                   ###
#####################################################
simulate_rbgarch_dpg_02 <- function(n = 1000, log_sigma_bar, rho, mu, steps = 100000) {
  n_tot <- n + 500
  log_sigma <- rep(NA, n_tot)
  ochl <- matrix(NA, ncol = 4, nrow = n_tot)
  e <- rnorm(n_tot)
  
  S <- 10
  ochl[1, 1] <- S
  log_sigma[1] <- log_sigma_bar 
  S <- c(ochl[1, 1], ochl[1, 1] + cumsum(rnorm(steps, 0, exp(log_sigma[1])/sqrt(steps))))
  ochl[1, 2:4] <- c(S[steps], max(S), min(S))
  
  for (i in 2:n_tot) {
    log_sigma[i] <- log_sigma_bar + rho * (log_sigma[i - 1] - log_sigma_bar) + mu * e[i - 1]
    ochl[i, 1] <- S[steps]
    S <- c(ochl[i, 1], ochl[i, 1] + cumsum(rnorm(steps, 0, exp(log_sigma[i])/sqrt(steps))))
    ochl[i, 2:4] <- c(S[steps], max(S), min(S))
  }
  out <- cbind(exp(log_sigma), ochl)
  colnames(out) <- c("sigma", "open", "close", "high", "low")
  out <- out[501:n_tot, ]
  return(out)
}



simulate_rbgarch_dpg_01 <- function(n = 1000, omega, alpha, beta, steps = 100000) {
  n_tot <- n + 500
  sigma2 <- rep(NA, n_tot)
  ochl <- matrix(NA, ncol = 4, nrow = n_tot)
  
  S <- 10
  ochl[1, 1] <- S
  sigma2[1] <- omega / (1 - alpha - beta)
  returns_intraday <- rnorm(steps, 0, sqrt(sigma2[1])/sqrt(steps))
  S <- c(ochl[1, 1], ochl[1, 1] + cumsum(returns_intraday))
  ochl[1, 2:4] <- c(S[steps], max(S), min(S))
  r <- ochl[1, 2] - ochl[1, 1]
  for (i in 2:n_tot) {
    sigma2[i] <- omega + alpha * r^2 + beta * sigma2[i - 1]
    ochl[i, 1] <- S[steps]
    returns_intraday <- rnorm(steps, 0, sqrt(sigma2[i])/sqrt(steps))
    S <- c(ochl[i, 1], ochl[i, 1] + cumsum(returns_intraday))
    ochl[i, 2:4] <- c(S[steps], max(S), min(S))
    r <- ochl[i, 2] - ochl[i, 1]
  }
  out <- cbind(sqrt(sigma2), ochl)
  colnames(out) <- c("sigma", "open", "close", "high", "low")
  out <- out[501:n_tot, ]
  return(out)
}

simulate_rbleverage_dgp_02 <- function(n = 1000, log_sigma2_ini, alpha, phi, rho, sigma_nu, steps = 100000) {
  n_tot <- n + 500
  v <- rep(NA, n_tot)
  u <- rep(NA, n_tot)
  w <- rep(NA, n_tot)
  log_sigma2 <- rep(NA, n_tot)
  r <- rep(NA, n_tot)
  M <- chol(matrix(c(1, rho, rho, 1), 2, 2))
  ochl <- matrix(NA, ncol = 4, nrow = n_tot)
  
  S <- 10
  ochl[1, 1] <- S
  log_sigma2[1] <- log_sigma2_ini
  S <- c(ochl[1, 1], ochl[1, 1] + cumsum(rnorm(steps, 0, sqrt(exp(log_sigma2[1]))/sqrt(steps))))
  ochl[1, 2:4] <- c(S[steps], max(S), min(S))
  r[1] <- ochl[1, 2] - ochl[1, 1]
  u[1] <- r[1] / sqrt(exp(log_sigma2[1]))
  aux_rho <- sqrt(1 - rho^2)
  for (i in 2:n_tot) {
    v[i] <- (c(u[i - 1], rnorm(1)) %*% chol(M))[2]
    w[i] <- (v[i] - rho * u[i - 1]) / aux_rho
    log_sigma2[i] <- alpha + phi * log_sigma2[i - 1]  + rho * sigma_nu * r[i - 1] / sqrt(exp(log_sigma2[i - 1])) + sigma_nu * aux_rho * w[i]
    ochl[i, 1] <- S[steps]
    S <- c(ochl[i, 1], ochl[i, 1] + cumsum(rnorm(steps, 0, sqrt(exp(log_sigma2[i]))/sqrt(steps))))
    ochl[i, 2:4] <- c(S[steps], max(S), min(S))
    r[i] <- ochl[i, 2] - ochl[i, 1]
    u[i] <- r[i] / sqrt(exp(log_sigma2[i]))
  }
  out <- cbind(sqrt(exp(log_sigma2)), ochl)
  colnames(out) <- c("sigma", "open", "close", "high", "low")
  out <- out[501:n_tot, ]
  return(out)
}



simulate_rbleverage_dgp_01 <- function(n = 1000, omega, alpha, beta, gama, steps = 100000)  {
  n_tot <- n + 500
  sigma2 <- rep(NA, n_tot)
  ochl <- matrix(NA, ncol = 4, nrow = n_tot)
  
  S <- 10
  ochl[1, 1] <- S
  sigma2[1] <- omega / (1 - alpha - beta - gama/2)
  returns_intraday <- rnorm(steps, 0, sqrt(sigma2[1])/sqrt(steps))
  S <- c(ochl[1, 1], ochl[1, 1] + cumsum(returns_intraday))
  ochl[1, 2:4] <- c(S[steps], max(S), min(S))
  r <- ochl[1, 2] - ochl[1, 1]
  for (i in 2:n_tot) {
    sigma2[i] <- omega + alpha * r^2 + beta * sigma2[i - 1] + gama * ifelse(r< 0, r^2, 0)
    ochl[i, 1] <- S[steps]
    returns_intraday <- rnorm(steps, 0, sqrt(sigma2[i])/sqrt(steps))
    S <- c(ochl[i, 1], ochl[i, 1] + cumsum(returns_intraday))
    ochl[i, 2:4] <- c(S[steps], max(S), min(S))
    r <- ochl[i, 2] - ochl[i, 1]
  }
  out <- cbind(sqrt(sigma2), ochl)
  colnames(out) <- c("sigma", "open", "close", "high", "low")
  out <- out[501:n_tot, ]
  return(out)
}



simulate_rbleverage_dgp_02_old <- function(n = 1000, log_sigma2_ini, alpha, phi, rho, steps = 100000) {
  n_tot <- n + 500
  v <- rep(NA, n_tot)
  u <- rep(NA, n_tot)
  log_sigma2 <- rep(NA, n_tot)
  M <- chol(matrix(c(1, rho, rho, 1), 2, 2))
  ochl <- matrix(NA, ncol = 4, nrow = n_tot)

  S <- 10
  ochl[1, 1] <- S
  log_sigma2[1] <- log_sigma2_ini
  S <- c(ochl[1, 1], ochl[1, 1] + cumsum(rnorm(steps, 0, sqrt(exp(log_sigma2[1]))/sqrt(steps))))
  ochl[1, 2:4] <- c(S[steps], max(S), min(S))
  u[1] <- (ochl[1, 2] - ochl[1, 1]) / sqrt(exp(log_sigma2[1]))
  
  for (i in 2:n_tot) {
    v[i] <- (c(u[i - 1], rnorm(1)) %*% chol(M))[2]
    log_sigma2[i] <- alpha + phi * log_sigma2[i - 1]  + 0.15*v[i]
    ochl[i, 1] <- S[steps]
    S <- c(ochl[i, 1], ochl[i, 1] + cumsum(rnorm(steps, 0, sqrt(exp(log_sigma2[i]))/sqrt(steps))))
    ochl[i, 2:4] <- c(S[steps], max(S), min(S))
    u[i] <- (ochl[i, 2] - ochl[i, 1]) / sqrt(exp(log_sigma2[i]))
  }
  out <- cbind(sqrt(exp(log_sigma2)), ochl)
  colnames(out) <- c("sigma", "open", "close", "high", "low")
  out <- out[501:n_tot, ]
  return(out)
}

#####################################################
###     Range-based volatility  proxies           ###
#####################################################

func_parkinson <- function(H, L, log_price = FALSE){
  if (log_price == FALSE) {
    P <- ((log(H/L))^2)/(4*log(2))
  } else {
    P <- ((H - L)^2)/(4*log(2))
  }
  return(P)
}

func_garman_klass <- function(H, L, C, O, log_price = FALSE){
  GK <- c()
  if (log_price == FALSE) {
    for (i in 1:length(H)) {
      GK[i] <- 0.5*((log(H[i]/L[i]))^2) - (2*log(2) - 1)*((log(C[i]/O[i]))^2)
    }
  } else {
    for (i in 1:length(H)) {
      GK[i] <- 0.5*(H[i] - L[i])^2 - (2*log(2) - 1)*(C[i] - O[i])^2
    }
  }
  return(GK)
}

func_meilijson <- function(H, L, C, O, log_price = FALSE) {
  M <- c()
  sigma1 <- c()
  sigma3 <- c()
  sigma4 <- c()
  sigmar <- c()
  if (log_price == FALSE) {
    for (i in 1:length(H)) {
      aux_c <- log(C[i]/O[i])
      aux_h <- log(H[i]/O[i])
      aux_l <- log(L[i]/O[i])
      
      if (aux_c > 0){
        c_line <- aux_c
        h_line <- aux_h
        l_line <- aux_l
      } else {
        c_line <- -aux_c
        h_line <- -aux_l
        l_line <- -aux_h
      }
      sigmar[i] <- c_line^2
      sigma1[i] <- 2*((h_line - c_line)^2 + l_line^2)
      sigma3[i] <- 2*(h_line - c_line - l_line)*c_line
      sigma4[i] <- -(h_line - c_line)*l_line/(2*log(2) - 5/4)
      
      M[i] <- 0.274*sigma1[i] + 0.16*sigmar[i] + 0.365*sigma3[i] + 0.2*sigma4[i] 
    }
  } else {
    for (i in 1:length(H)) {
      aux_c <- C[i] - O[i]
      aux_h <- H[i] - O[i]
      aux_l <- L[i] - O[i]
      if (aux_c > 0){
        c_line <- aux_c
        h_line <- aux_h
        l_line <- aux_l
      } else {
        c_line <- -aux_c
        h_line <- -aux_l
        l_line <- -aux_h
      }
      sigmar[i] <- c_line^2
      sigma1[i] <- 2*((h_line - c_line)^2 + l_line^2)
      sigma3[i] <- 2*(h_line - c_line - l_line)*c_line
      sigma4[i] <- -(h_line - c_line)*l_line/(2*log(2) - 5/4)
      
      M[i] <- 0.274*sigma1[i] + 0.16*sigmar[i] + 0.365*sigma3[i] + 0.2*sigma4[i] 
    }
  }
  
  return(M)
}

func_rogers_satchell <- function(H, L, C, O, log_price = FALSE){
  RS <- c()
  if (log_price == FALSE) {
    for (i in 1:length(H)) {
      aux_c <- log(C[i]/O[i])
      aux_h <- log(H[i]/O[i])
      aux_l <- log(L[i]/O[i])
      RS[i] <- aux_h*(aux_h - aux_c) + aux_l*(aux_l - aux_c)
    }
  } else {
    for (i in 1:length(H)) {
      aux_c <- C[i] - O[i]
      aux_h <- H[i] - O[i]
      aux_l <- L[i] - O[i]
      RS[i] <- aux_h*(aux_h - aux_c) + aux_l*(aux_l - aux_c)
    }
  }
  return(RS)
}

loss_mse <- function(h_hat, h) {
  mean((h_hat - h)^2)
}

loss_qlike <- function(h_hat, h) {
  mean(log(h) + h_hat/h)
}

loss_qlike_ult <- function(h_hat, h) {
  mean(h_hat/h - log(h_hat/h) - 1)
}

loss_mse_log <- function(h_hat, h) {
  mean((log(h_hat) - log(h))^2)
}

loss_mse_sd <- function(h_hat, h) {
  mean((sqrt(h_hat) - sqrt(h))^2)
}

loss_mse_prop <- function(h_hat, h) {
  mean((h_hat / h - 1)^2)
}

loss_mae <- function(h_hat, h) {
  mean(abs(h_hat - h))
}

loss_mae_log <- function(h_hat, h) {
  mean(abs(log(h_hat) - log(h)))
}

loss_mae_sd <- function(h_hat, h) {
  mean(abs(sqrt(h_hat) - sqrt(h)))
}

loss_mae_prop <- function(h_hat, h) {
  mean(abs(h_hat / h - 1))
}


#####################################################
###       Expected shortfall estimation           ###
#####################################################
xfx = function(x, mu_, sigma_) x*ddist(distribution = "norm", y = x, mu = mu_, sigma = sigma_, skew = 0)


#####################################################
###           Calibration test VQR                ###
### Gaglianone et al. (2011).                     ###
#####################################################

#function para teste de calibracao
VaR_VQR = function(r,VaR, alpha){
  
  fit1 = suppressWarnings(summary(rq(r ~ VaR, tau = alpha, method = "fn"), method="fn" , se="nid" , cov=TRUE))
  
  a1 = fit1$coefficients[1]
  a2 = fit1$coefficients[2]
  
  M = matrix(nrow = 2 , ncol=1)
  M[1,1] = a1
  M[2,1] = (a2-1)
  
  icov = matrix(nrow = 2 , ncol = 2)
  aa = fit1$cov[1,1]
  bb = fit1$cov[1,2]
  cc = fit1$cov[2,1]
  dd = fit1$cov[2,2]
  icov[2,1] = 1/(bb-aa*dd/cc)
  icov[2,2] = 1/(dd-cc*bb/aa)
  icov[1,1] = -icov[2,1]*dd/cc
  icov[1,2] = -icov[2,2]*bb/aa
  
  statistic = (t(M)) %*% icov %*% M 

  p.value = 1-pchisq(statistic[1,1], df=2)
  return(p.value)
}

