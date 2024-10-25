#####################################
###      RBGARCH: Estimation      ###
#####################################

## Parameter estimation
rbgarch <- function(y, proxi) {
  par_ini <-  grid_rbgarch(y, proxi)
  ra <- matrix(c(1, 0, 0,
                 0, 1, 0, 
                 0, 0, 1,
                 0, -1, -1), ncol = 3, byrow = TRUE)
  rb <- c(1e-06, 0, 0,-0.99999)
  param <- constrOptim(theta = par_ini, f = rbgarch_loglik, grad = NULL, ui = ra, ci = rb, r = y, proxi = proxi, outer.iterations = 400, outer.eps = 1e-07)$par
}


rbgarch_t_student <- function(y, proxi) {
  par_ini <-  grid_rbgarch_t_student(y, proxi)
  ra <- matrix(c(1, 0, 0, 0, 
                 0, 1, 0, 0,
                 0, 0, 1, 0, 
                 0, 0, 0, 1, 
                 0, -1, -1, 0), ncol = 4, byrow = TRUE)
  rb <- c(1e-06, 0, 0, 2.00001, -0.99999)
  param <- constrOptim(theta = par_ini, f = rbgarch_loglik_t_student, grad = NULL, ui = ra, ci = rb, r = y, proxi = proxi, outer.iterations = 400, outer.eps = 1e-07)$par
}



rbgjr2<- function(y, proxi) {
  par_ini <-  grid_rbgjr(y, proxi)
  ra <- matrix(c(1, 0, 0, 0,
                 0, 1, 0, 0, 
                 0, 0, 1, 0,
                 0, -1, -1, -0.5), ncol = 4, byrow = TRUE)
  rb <- c(1e-06, 0, 0, -0.99999)
  param <- constrOptim(theta = par_ini, f = rbgjr_loglik, grad = NULL, ui = ra, ci = rb, r = y, proxi = proxi, outer.iterations = 400, outer.eps = 1e-07)$par
  return(param)
}


rbgjr<- function(y, proxi) {
  par_ini <-  grid_rbgjr(y, proxi)
  ra <- matrix(c(1, 0, 0, 0,
                 0, 1, 0, 0, 
                 0, 0, 1, 0,
                 0, 0, 0, 1,
                 0, -1, -1, -0.5), ncol = 4, byrow = TRUE)
  rb <- c(1e-06, 0, 0, 0, -0.99999)
  param <- constrOptim(theta = par_ini, f = rbgjr_loglik, grad = NULL, ui = ra, ci = rb, r = y, proxi = proxi, outer.iterations = 400, outer.eps = 1e-07)$par
  return(param)
}



rbgjr_t_student<- function(y, proxi) {
  par_ini <-  grid_rbgjr_t_student(y, proxi)
  ra <- matrix(c(1, 0, 0, 0, 0, 
                 0, 1, 0, 0, 0, 
                 0, 0, 1, 0, 0,
                 0, 0, 0, 1, 0, 
                 0, 0, 0, 0, 1,  
                 0, -1, -1, -0.5, 0), ncol = 5, byrow = TRUE)
  rb <- c(1e-06, 0, 0, 0, 2.00001, -0.99999)
  param <- constrOptim(theta = par_ini, f = rbgjr_loglik_t_student, grad = NULL, ui = ra, ci = rb, r = y, proxi = proxi, outer.iterations = 400, outer.eps = 1e-07)$par
  return(param)
}



rbigarch <- function(y, proxi) {
  par_ini <-  grid_rbigarch(y, proxi)
  ra <- matrix(c(1, 0,
                 0, 1, 
                 0, -1), ncol = 2, byrow = TRUE)
  rb <- c(1e-06, 0,-1)
  param <- constrOptim(theta = par_ini, f = rbigarch_loglik, grad = NULL, ui = ra, ci = rb, r = y, proxi = proxi, outer.iterations = 400, outer.eps = 1e-07)$par
}


rbigarch_t_student <- function(y, proxi) {
  par_ini <-  grid_rbigarch_t_student(y, proxi)
  ra <- matrix(c(1, 0, 0, 
                 0, 1, 0,
                 0, -1, 0, 
                 0, 0, 1), ncol = 3, byrow = TRUE)
  rb <- c(1e-06, 0, -1, 2.00001)
  param <- constrOptim(theta = par_ini, f = rbigarch_loglik_t_student, grad = NULL, ui = ra, ci = rb, r = y, proxi = proxi, outer.iterations = 400, outer.eps = 1e-07)$par
}


## Volatility estimation
rbgarch_h <- function(y, proxi, params) {
  n <- length(y)
  h <- rep(NA, n)
  h[1] <- params[1] / (1 - params[2] - params[3])
  for (i in 2:(n + 1)) {
    h[i] = params[1] + params[2] * proxi[i-1] + params[3] * h[i-1]
  }
  return(h)
}

rbgjr_h <- function(y, proxi, params) {
  n <- length(y)
  h <- rep(NA, n)
  h[1] <- params[1] / (1 - params[2] - params[3] - 0.5* params[4])
  for (i in 2:(n + 1)) {
    if (y[i - 1] < 0) {
      h[i] = params[1] + params[2] * proxi[i-1] + params[3] * h[i-1] + params[4] * proxi[i-1]
    } else {
      h[i] = params[1] + params[2] * proxi[i-1] + params[3] * h[i-1]
    }
  }
  return(h)
}


rbigarch_h <- function(y, proxi, params) {
  n <- length(y)
  h <- rep(NA, n)
  h[1] <- var(y)
  for (i in 2:(n + 1)) {
    h[i] = params[1] + params[2] * proxi[i-1] + (1-params[2]) * h[i-1]
  }
  return(h)
}

