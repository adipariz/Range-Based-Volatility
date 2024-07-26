###################################
###    Empirical Application    ###
###################################
library(rugarch)
library(Rcpp)
library(GAS)
library(esback)
library(esreg)
library(quantreg)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(xtable)
sourceCpp("Estimation.cpp")
source("Estimation.R")
source("Utils.R")

ibrx_names <- names(read.csv("ibrx_names.csv"))
ibrx_stocks <- read.csv("ochl_ibrx.csv", na.strings = c("NA", "-")) |>  mutate(Data = dmy(Data)) |> mutate_if(is.character, as.numeric) 

# Setup
spec_garch <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)), distribution = 'norm')
spec_gjr <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = list(model = 'gjrGARCH', garchOrder = c(1, 1)), distribution = 'norm')
janela_estimacao <- 1250

nomes <- c("GARCH", "RBGARCH-P", "RBGARCH-GK", "RBGARCH-M", "RBGARCH-RS", "GJR",  "RBGJR-P", "RBGJR-GK", "RBGJR-M", "RBGJR-RS")
sigma_fit_fore <- matrix(NA, ncol = 10, nrow = janela_estimacao + 1, dimnames = list(NULL, nomes))
e <- matrix(NA, ncol = 10, nrow = janela_estimacao)
alpha1    <- 0.01
alpha2.5  <- 0.025
alpha5    <- 0.05




for (stock_name in ibrx_names) {
  # Reading Data
  stock <- ibrx_stocks |> select(Data, contains(stock_name)) |> 
    filter(Data > "2010-01-01" & Data < "2024-07-01") |> 
    drop_na() |> 
    rename_with(~ str_sub(.x, 1, 4)) |> 
    rename("date" = "Data", "open" = "Aber", "close" = "Fech", "higher" = "Máxi", "lower" = "Míni") |> 
    mutate(ret = price_to_returns(close),
           parkinson = func_parkinson(higher, lower, log_price = FALSE), 
           garman_Klass = func_garman_klass(higher, lower, close, open, log_price = FALSE), 
           meilijson = func_meilijson(higher, lower, close, open, log_price = FALSE), 
           rogers_satchell = func_rogers_satchell(higher, lower, close, open, log_price = FALSE)) |> 
    select(date, ret, parkinson, garman_Klass, meilijson , rogers_satchell) |> drop_na()
  
  # DataFrames
  oos <- nrow(stock) - janela_estimacao
  
  var_1 <- matrix(NA, ncol = 10, nrow = oos, dimnames = list(NULL, nomes))
  var_2 <- matrix(NA, ncol = 10, nrow = oos, dimnames = list(NULL, nomes))
  var_5 <- matrix(NA, ncol = 10, nrow = oos, dimnames = list(NULL, nomes))
  
  var_1_normal <- matrix(NA, ncol = 10, nrow = oos, dimnames = list(NULL, nomes))
  var_2_normal <- matrix(NA, ncol = 10, nrow = oos, dimnames = list(NULL, nomes))
  var_5_normal <- matrix(NA, ncol = 10, nrow = oos, dimnames = list(NULL, nomes))
  
  es_1 <- matrix(NA, ncol = 10, nrow = oos, dimnames = list(NULL, nomes))
  es_2 <- matrix(NA, ncol = 10, nrow = oos, dimnames = list(NULL, nomes))
  es_5 <- matrix(NA, ncol = 10, nrow = oos, dimnames = list(NULL, nomes))
  
  es_1_normal <- matrix(NA, ncol = 10, nrow = oos, dimnames = list(NULL, nomes))
  es_2_normal <- matrix(NA, ncol = 10, nrow = oos, dimnames = list(NULL, nomes))
  es_5_normal <- matrix(NA, ncol = 10, nrow = oos, dimnames = list(NULL, nomes))
  sigma2_fit <- matrix(NA, ncol = 11, nrow = oos)
  
  
  
  # Rolling Window
  for (i in 1:oos) {
    print(i)
    stock_window <- stock[i:(i + janela_estimacao - 1), ]
    mu <- mean(stock_window$ret)
    stock_window_ret <- scale(stock_window$ret, center = TRUE, scale = FALSE)
    
    fit_garch <- ugarchfit(spec_garch, stock_window_ret, solver = 'hybrid')
    fit_gjr   <- ugarchfit(spec_gjr, stock_window_ret, solver = 'hybrid')
    parametros_rbgarch_p  <-  rbgarch(stock_window_ret, stock_window$parkinson)
    parametros_rbgarch_gk <-  rbgarch(stock_window_ret, stock_window$garman_Klass)
    parametros_rbgarch_m  <-  rbgarch(stock_window_ret, stock_window$meilijson)
    parametros_rbgarch_rs <-  rbgarch(stock_window_ret, stock_window$rogers_satchell)
    
    parametros_rbgjr_p  <-  rbgjr(stock_window_ret, stock_window$parkinson)
    parametros_rbgjr_gk <-  rbgjr(stock_window_ret, stock_window$garman_Klass)
    parametros_rbgjr_m  <-  rbgjr(stock_window_ret, stock_window$meilijson)
    parametros_rbgjr_rs <-  rbgjr(stock_window_ret, stock_window$rogers_satchell)
    
    sigma_fit_fore[, "GARCH"] <- c(fit_garch@fit$sigma, as.numeric(ugarchforecast(fit_garch, n.ahead = 1)@forecast$sigmaFor))
    sigma_fit_fore[, "GJR"] <- c(fit_gjr@fit$sigma, as.numeric(ugarchforecast(fit_gjr, n.ahead = 1)@forecast$sigmaFor))
    
    sigma_fit_fore[, "RBGARCH-P"] <- sqrt(rbgarch_h(stock_window_ret, stock_window$parkinson, parametros_rbgarch_p))
    sigma_fit_fore[, "RBGJR-P"] <- sqrt(rbgjr_h(stock_window_ret, stock_window$parkinson, parametros_rbgjr_p))
    
    sigma_fit_fore[, "RBGARCH-GK"] <- sqrt(rbgarch_h(stock_window_ret, stock_window$garman_Klass, parametros_rbgarch_gk))
    sigma_fit_fore[, "RBGJR-GK"] <- sqrt(rbgjr_h(stock_window_ret, stock_window$garman_Klass, parametros_rbgjr_gk))
    
    sigma_fit_fore[, "RBGARCH-M"] <- sqrt(rbgarch_h(stock_window_ret, stock_window$meilijson, parametros_rbgarch_m))
    sigma_fit_fore[, "RBGJR-M"] <- sqrt(rbgjr_h(stock_window_ret, stock_window$meilijson, parametros_rbgjr_m))
    
    sigma_fit_fore[, "RBGARCH-RS"] <- sqrt(rbgarch_h(stock_window_ret, stock_window$rogers_satchell, parametros_rbgarch_rs))
    sigma_fit_fore[, "RBGJR-RS"] <- sqrt(rbgarch_h(stock_window_ret, stock_window$rogers_satchell, parametros_rbgjr_rs))
    
    e <- rep(stock_window_ret, 10) / sigma_fit_fore[1:janela_estimacao, ]
    ret_one_step_ahead <- matrix(tail(sigma_fit_fore, 1), ncol = 10, nrow = janela_estimacao, byrow = TRUE) * e + mu
    vars <-  apply(ret_one_step_ahead, 2, quantile, prob = c(0.01, 0.025, 0.05))
    
    
    ess <- matrix(NA, ncol = 10, nrow = 3)
    for (j in 1:10) {
      for (k in 1:3) {
        ess[k, j] <- mean(ret_one_step_ahead[ret_one_step_ahead[, j] < vars[k, j], j])
      }
    }
    
    var_1[i, ] <- vars[1, ]
    var_2[i, ] <- vars[2, ]
    var_5[i, ] <- vars[3, ]
    es_1[i, ] <- ess[1, ]
    es_2[i, ] <- ess[2, ]
    es_5[i, ] <- ess[3, ]
    
    sigma2_fit[i,] <- c(tail(sigma_fit_fore, 1)^2, stock[i + janela_estimacao, ]$ret^2)
    
    var_1_normal[i, ] <- tail(sigma_fit_fore, 1) * qnorm(0.010) 
    var_2_normal[i, ] <- tail(sigma_fit_fore, 1) * qnorm(0.025) 
    var_5_normal[i, ] <- tail(sigma_fit_fore, 1) * qnorm(0.050)  
    
    for (j in 1:10){
      es_1_normal[i, j] <- integrate(xfx, -Inf, var_1_normal[i, j],  mu_ = 0, sigma_ = tail(sigma_fit_fore[, j], 1))$value/alpha1
      es_2_normal[i, j] <- integrate(xfx, -Inf, var_2_normal[i, j],  mu_ = 0, sigma_ = tail(sigma_fit_fore[, j], 1))$value/alpha2.5
      es_5_normal[i, j] <- integrate(xfx, -Inf, var_5_normal[i, j],  mu_ = 0, sigma_ = tail(sigma_fit_fore[, j], 1))$value/alpha5
      
    }
    
    var_1_normal[i, ] <- var_1_normal[i, ] + mu
    var_2_normal[i, ] <- var_2_normal[i, ] + mu
    var_5_normal[i, ] <- var_5_normal[i, ] + mu
    
    es_1_normal[i, ] <- es_1_normal[i, ] + mu
    es_2_normal[i, ] <- es_2_normal[i, ] + mu
    es_5_normal[i, ] <- es_5_normal[i, ] + mu
  }
  
  write.csv(sigma2_fit, paste0("sigma2_", stock_name, ".csv"))
  write.csv(var_1, paste0("var_1_", stock_name, ".csv"))
  write.csv(var_2, paste0("var_2_", stock_name, ".csv"))
  write.csv(var_5, paste0("var_5_", stock_name, ".csv"))
  write.csv(var_1_normal, paste0("var_1_normal_", stock_name, ".csv"))
  write.csv(var_2_normal, paste0("var_2_normal_", stock_name, ".csv"))
  write.csv(var_5_normal, paste0("var_5_normal_", stock_name, ".csv"))
  write.csv(es_1, paste0("es_1_", stock_name, ".csv"))
  write.csv(es_2, paste0("es_2_", stock_name, ".csv"))
  write.csv(es_5, paste0("es_5_", stock_name, ".csv"))
  write.csv(es_1_normal, paste0("es_1_normal_", stock_name, ".csv"))
  write.csv(es_2_normal, paste0("es_2_normal_", stock_name, ".csv"))
  write.csv(es_5_normal, paste0("es_5_normal_", stock_name, ".csv"))
  
}
  


