#######################################
###      Monte Carlo Simulation    ####
#######################################
args <- commandArgs(trailingOnly=TRUE)
library(rugarch)
library(Rcpp)
library(dplyr)
sourceCpp("Estimation.cpp")
source("Estimation.R")
source("Utils.R")



spec0 <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                    variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)), 
                    distribution = 'norm')


mc <- 10000
n <- as.numeric(args[1])
sigma2_hat <- matrix(NA, ncol = 6, nrow = mc)
colnames(sigma2_hat) <- c("TRUE", "GARCH", "RB-Parkinson", "RB-GK", "RB-Meilijson", "RB-RS")


for (i in 1:mc) {
  print(i)
  set.seed(i)
  aux <- simulate_rbgarch_dpg_01(n = 2501, omega = 0.01, alpha = 0.05, beta = 0.92, steps = 100000) |> tail(n) |> data.frame() |>
    mutate(ret = close - open,
           parkinson = func_parkinson(exp(high), exp(low), log_price = FALSE), 
           garman_Klass = func_garman_klass(exp(high), exp(low), exp(close), exp(open), log_price = FALSE),
           meilijson = func_meilijson(exp(high), exp(low), exp(close), exp(open), log_price = FALSE),
           rogers_satchell = func_rogers_satchell(exp(high), exp(low), exp(close), exp(open), log_price = FALSE))
  
  sigma2_hat[i, "TRUE"] <- aux[n, 1]^2
  sigma2_hat[i, "GARCH"] <- as.numeric(ugarchforecast(ugarchfit(spec0, aux$ret[1:(n - 1)], solver = 'hybrid'), n.ahead = 1)@forecast$sigmaFor)^2
  sigma2_hat[i, "RB-Parkinson"] <- tail(rbgarch_h(aux$ret[1:(n - 1)], aux$parkinson[1:(n - 1)], rbgarch(aux$ret[1:(n - 1)], aux$parkinson[1:(n - 1)])), 1)
  sigma2_hat[i, "RB-GK"] <- tail(rbgarch_h(aux$ret[1:(n - 1)], aux$garman_Klass[1:(n - 1)], rbgarch(aux$ret[1:(n - 1)], aux$garman_Klass[1:(n - 1)])), 1)
  sigma2_hat[i, "RB-Meilijson"] <- tail(rbgarch_h(aux$ret[1:(n - 1)], aux$meilijson[1:(n - 1)], rbgarch(aux$ret[1:(n - 1)], aux$meilijson[1:(n - 1)])), 1)
  sigma2_hat[i, "RB-RS"] <- tail(rbgarch_h(aux$ret[1:(n - 1)], aux$rogers_satchell[1:(n - 1)], rbgarch(aux$ret[1:(n - 1)], aux$rogers_satchell[1:(n - 1)])), 1)
}


comparison <- matrix(NA, ncol = 5, nrow = 9)
for (i in 2:ncol(sigma2_hat)) {
  comparison[1, i - 1] <- 100*round(loss_mse(sigma2_hat[, i], sigma2_hat[, 1]), 5)
  comparison[2, i - 1] <- 100*round(loss_qlike(sigma2_hat[, i], sigma2_hat[, 1]), 5)
  comparison[3, i - 1] <- 100*round(loss_mse_log(sigma2_hat[, i], sigma2_hat[, 1]), 5)
  comparison[4, i - 1] <- 100*round(loss_mse_sd(sigma2_hat[, i], sigma2_hat[, 1]), 5)
  comparison[5, i - 1] <- 100*round(loss_mse_prop(sigma2_hat[, i], sigma2_hat[, 1]), 5)
  comparison[6, i - 1] <- 100*round(loss_mae(sigma2_hat[, i], sigma2_hat[, 1]), 5)
  comparison[7, i - 1] <- 100*round(loss_mae_log(sigma2_hat[, i], sigma2_hat[, 1]), 5)
  comparison[8, i - 1] <- 100*round(loss_mae_sd(sigma2_hat[, i], sigma2_hat[, 1]), 5)
  comparison[9, i - 1] <- 100*round(loss_mae_prop(sigma2_hat[, i], sigma2_hat[, 1]), 5)
}

colnames(comparison) <- c("GARCH", "RB-Parkinson", "RB-GK", "RB-Meilijson", "RB-RS")
rownames(comparison) <- c("MSE", "QLIKE", "MSE-LOG", "MSE-SD", "MSE-prop", "MAE", "MAE-LOG", "MAE-SD", "MAE-prop")
write.csv(comparison, paste0("comparison_RBGARCH_", n - 1, ".csv"))









