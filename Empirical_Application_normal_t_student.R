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

ibrx_names <- names(read.csv("ibrx_names.csv")[, -c(1:23)])
ibrx_stocks <- read.csv("ochl_ibrx.csv", na.strings = c("NA", "-")) |>  mutate(Data = dmy(Data)) |> mutate_if(is.character, as.numeric) 

# Setup

spec_garch_n <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)), distribution = 'norm')
spec_gjr_n <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = list(model = 'gjrGARCH', garchOrder = c(1, 1)), distribution = 'norm')
spec_garch_t <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)), distribution = 'std')
spec_gjr_t <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = list(model = 'gjrGARCH', garchOrder = c(1, 1)), distribution = 'std')
janela_estimacao <- 1250

nomes <- c("GARCH", "RBGARCH-P", "RBGARCH-GK", "RBGARCH-M", "RBGARCH-RS", "GJR",  "RBGJR-P", "RBGJR-GK", "RBGJR-M", "RBGJR-RS")
sigma_fit_fore_n <- matrix(NA, ncol = 10, nrow = janela_estimacao + 1, dimnames = list(NULL, nomes))
sigma_fit_fore_t <- matrix(NA, ncol = 10, nrow = janela_estimacao + 1, dimnames = list(NULL, nomes))

e_n <- matrix(NA, ncol = 10, nrow = janela_estimacao)
e_t <- matrix(NA, ncol = 10, nrow = janela_estimacao)



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
  
  var_1_n <- matrix(NA, ncol = 10, nrow = oos, dimnames = list(NULL, nomes))
  var_2_n <- matrix(NA, ncol = 10, nrow = oos, dimnames = list(NULL, nomes))
  
  var_1_t <- matrix(NA, ncol = 10, nrow = oos, dimnames = list(NULL, nomes))
  var_2_t <- matrix(NA, ncol = 10, nrow = oos, dimnames = list(NULL, nomes))
  
  
  
  es_1_n <- matrix(NA, ncol = 10, nrow = oos, dimnames = list(NULL, nomes))
  es_2_n <- matrix(NA, ncol = 10, nrow = oos, dimnames = list(NULL, nomes))
  
  
  es_1_t <- matrix(NA, ncol = 10, nrow = oos, dimnames = list(NULL, nomes))
  es_2_t <- matrix(NA, ncol = 10, nrow = oos, dimnames = list(NULL, nomes))
  
  var_1 <- matrix(NA, ncol = 20, nrow = oos, dimnames = list(NULL, rep(nomes, times = 2)))
  var_2 <- matrix(NA, ncol = 20, nrow = oos, dimnames = list(NULL, rep(nomes, times = 2)))
  
  es_1 <- matrix(NA, ncol = 20, nrow = oos, dimnames = list(NULL, rep(nomes, times = 2)))
  es_2 <- matrix(NA, ncol = 20, nrow = oos, dimnames = list(NULL, rep(nomes, times = 2)))
  
  
  sigma2_fit_n <- matrix(NA, ncol = 11, nrow = oos)
  sigma2_fit_t <- matrix(NA, ncol = 11, nrow = oos)
  
  result = tryCatch({
    
    # Rolling Window
    for (i in 1:oos) {
      print(i)
      stock_window <- stock[i:(i + janela_estimacao - 1), ]
      mu <- mean(stock_window$ret)
      stock_window_ret <- scale(stock_window$ret, center = TRUE, scale = FALSE)
      
      fit_garch_n <- ugarchfit(spec_garch_n, stock_window_ret, solver = 'hybrid')
      fit_gjr_n   <- ugarchfit(spec_gjr_n, stock_window_ret, solver = 'hybrid')
      
      fit_garch_t <- ugarchfit(spec_garch_t, stock_window_ret, solver = 'hybrid')
      fit_gjr_t   <- ugarchfit(spec_gjr_t, stock_window_ret, solver = 'hybrid')
      
      
      parametros_rbgarch_p  <-  rbgarch(stock_window_ret, stock_window$parkinson)
      parametros_rbgarch_gk <-  rbgarch(stock_window_ret, stock_window$garman_Klass)
      parametros_rbgarch_m  <-  rbgarch(stock_window_ret, stock_window$meilijson)
      parametros_rbgarch_rs <-  rbgarch(stock_window_ret, stock_window$rogers_satchell)
      
      parametros_rbgjr_p  <-  rbgjr(stock_window_ret, stock_window$parkinson)
      parametros_rbgjr_gk <-  rbgjr(stock_window_ret, stock_window$garman_Klass)
      parametros_rbgjr_m  <-  rbgjr(stock_window_ret, stock_window$meilijson)
      parametros_rbgjr_rs <-  rbgjr(stock_window_ret, stock_window$rogers_satchell)
      
      
      
      parametros_rbgarch_p_t  <-  rbgarch_t_student(stock_window_ret, stock_window$parkinson)
      parametros_rbgarch_gk_t <-  rbgarch_t_student(stock_window_ret, stock_window$garman_Klass)
      parametros_rbgarch_m_t  <-  rbgarch_t_student(stock_window_ret, stock_window$meilijson)
      parametros_rbgarch_rs_t <-  rbgarch_t_student(stock_window_ret, stock_window$rogers_satchell)
      
      parametros_rbgjr_p_t  <-  rbgjr_t_student(stock_window_ret, stock_window$parkinson)
      parametros_rbgjr_gk_t <-  rbgjr_t_student(stock_window_ret, stock_window$garman_Klass)
      parametros_rbgjr_m_t  <-  rbgjr_t_student(stock_window_ret, stock_window$meilijson)
      parametros_rbgjr_rs_t <-  rbgjr_t_student(stock_window_ret, stock_window$rogers_satchell)
      
      
      sigma_fit_fore_n[, "GARCH"] <- c(fit_garch_n@fit$sigma, as.numeric(ugarchforecast(fit_garch_n, n.ahead = 1)@forecast$sigmaFor))
      sigma_fit_fore_n[, "GJR"] <- c(fit_gjr_n@fit$sigma, as.numeric(ugarchforecast(fit_gjr_n, n.ahead = 1)@forecast$sigmaFor))
      
      sigma_fit_fore_n[, "RBGARCH-P"] <- sqrt(rbgarch_h(stock_window_ret, stock_window$parkinson, parametros_rbgarch_p))
      sigma_fit_fore_n[, "RBGJR-P"] <- sqrt(rbgjr_h(stock_window_ret, stock_window$parkinson, parametros_rbgjr_p))
      
      sigma_fit_fore_n[, "RBGARCH-GK"] <- sqrt(rbgarch_h(stock_window_ret, stock_window$garman_Klass, parametros_rbgarch_gk))
      sigma_fit_fore_n[, "RBGJR-GK"] <- sqrt(rbgjr_h(stock_window_ret, stock_window$garman_Klass, parametros_rbgjr_gk))
      
      sigma_fit_fore_n[, "RBGARCH-M"] <- sqrt(rbgarch_h(stock_window_ret, stock_window$meilijson, parametros_rbgarch_m))
      sigma_fit_fore_n[, "RBGJR-M"] <- sqrt(rbgjr_h(stock_window_ret, stock_window$meilijson, parametros_rbgjr_m))
      
      sigma_fit_fore_n[, "RBGARCH-RS"] <- sqrt(rbgarch_h(stock_window_ret, stock_window$rogers_satchell, parametros_rbgarch_rs))
      sigma_fit_fore_n[, "RBGJR-RS"] <- sqrt(rbgarch_h(stock_window_ret, stock_window$rogers_satchell, parametros_rbgjr_rs))
      
      
      
      sigma_fit_fore_t[, "GARCH"] <- c(fit_garch_t@fit$sigma, as.numeric(ugarchforecast(fit_garch_t, n.ahead = 1)@forecast$sigmaFor))
      sigma_fit_fore_t[, "GJR"] <- c(fit_gjr_t@fit$sigma, as.numeric(ugarchforecast(fit_gjr_t, n.ahead = 1)@forecast$sigmaFor))
      
      sigma_fit_fore_t[, "RBGARCH-P"] <- sqrt(rbgarch_h(stock_window_ret, stock_window$parkinson, parametros_rbgarch_p_t))
      sigma_fit_fore_t[, "RBGJR-P"] <- sqrt(rbgjr_h(stock_window_ret, stock_window$parkinson, parametros_rbgjr_p_t))
      
      sigma_fit_fore_t[, "RBGARCH-GK"] <- sqrt(rbgarch_h(stock_window_ret, stock_window$garman_Klass, parametros_rbgarch_gk_t))
      sigma_fit_fore_t[, "RBGJR-GK"] <- sqrt(rbgjr_h(stock_window_ret, stock_window$garman_Klass, parametros_rbgjr_gk_t))
      
      sigma_fit_fore_t[, "RBGARCH-M"] <- sqrt(rbgarch_h(stock_window_ret, stock_window$meilijson, parametros_rbgarch_m_t))
      sigma_fit_fore_t[, "RBGJR-M"] <- sqrt(rbgjr_h(stock_window_ret, stock_window$meilijson, parametros_rbgjr_m_t))
      
      sigma_fit_fore_t[, "RBGARCH-RS"] <- sqrt(rbgarch_h(stock_window_ret, stock_window$rogers_satchell, parametros_rbgarch_rs_t))
      sigma_fit_fore_t[, "RBGJR-RS"] <- sqrt(rbgarch_h(stock_window_ret, stock_window$rogers_satchell, parametros_rbgjr_rs_t))
      
      
      
      
      
      e_n <- rep(stock_window_ret, 10) / sigma_fit_fore_n[1:janela_estimacao, ]
      e_t <- rep(stock_window_ret, 10) / sigma_fit_fore_t[1:janela_estimacao, ]
      
      ret_one_step_ahead <- matrix(tail(sigma_fit_fore_n, 1), ncol = 10, nrow = janela_estimacao, byrow = TRUE) * e_n + mu
      ret_one_step_ahead_t <- matrix(tail(sigma_fit_fore_t, 1), ncol = 10, nrow = janela_estimacao, byrow = TRUE) * e_t + mu
      vars <-  apply(ret_one_step_ahead, 2, quantile, prob = c(0.01, 0.025))
      vars_t <-  apply(ret_one_step_ahead_t, 2, quantile, prob = c(0.01, 0.025))
      
      ess <- matrix(NA, ncol = 10, nrow = 2)
      ess_t <- matrix(NA, ncol = 10, nrow = 2)
      for (j in 1:10) {
        for (k in 1:2) {
          ess[k, j] <- mean(ret_one_step_ahead[ret_one_step_ahead[, j] < vars[k, j], j])
          ess_t[k, j] <- mean(ret_one_step_ahead_t[ret_one_step_ahead_t[, j] < vars_t[k, j], j])
        }
      }
      
      var_1_n[i, ] <- vars[1, ]
      var_2_n[i, ] <- vars[2, ]
      
      es_1_n[i, ] <- ess[1, ]
      es_2_n[i, ] <- ess[2, ]
      
      
      var_1_t[i, ] <- vars_t[1, ]
      var_2_t[i, ] <- vars_t[2, ]
      
      es_1_t[i, ] <- ess_t[1, ]
      es_2_t[i, ] <- ess_t[2, ]
      
      ###########################################
      var_1[i, ] <- cbind(var_1_n[i, ], var_1_t[i, ])
      var_2[i, ] <- cbind(var_2_n[i, ], var_2_t[i, ])
      
      es_1[i, ] <- cbind(es_1_n[i, ], es_1_t[i, ])
      es_2[i, ] <- cbind(es_2_n[i, ], es_2_t[i, ])
      
      
      #########################################
      
      sigma2_fit_n[i,] <- c(tail(sigma_fit_fore_n, 1)^2, stock[i + janela_estimacao, ]$ret^2)
      sigma2_fit_t[i,] <- c(tail(sigma_fit_fore_t, 1)^2, stock[i + janela_estimacao, ]$ret^2)
    }
    
    
    write.csv(sigma2_fit_n, paste0("sigma2_normal_", stock_name, ".csv"))
    #write.csv(var_1_n, paste0("var_1_normal_", stock_name, ".csv"))
    #write.csv(var_2_n, paste0("var_2_normal_", stock_name, ".csv"))
    
    
    #write.csv(es_1_n, paste0("es_1_normal_", stock_name, ".csv"))
    #write.csv(es_2_n, paste0("es_2_normal_", stock_name, ".csv"))
    
    
    
    write.csv(sigma2_fit_t, paste0("sigma2_t_student_", stock_name, ".csv"))
    #write.csv(var_1_t, paste0("var_1_t_student_", stock_name, ".csv"))
    #write.csv(var_2_t, paste0("var_2_t_student_", stock_name, ".csv"))
    
    #write.csv(es_1_t, paste0("es_1_t_student_", stock_name, ".csv"))
    #write.csv(es_2_t, paste0("es_2_t_student_", stock_name, ".csv"))
    
    
    write.csv(var_1, paste0("var_1_", stock_name, ".csv"))
    write.csv(var_2, paste0("var_2_", stock_name, ".csv"))
    write.csv(es_1, paste0("es_1_", stock_name, ".csv"))
    write.csv(es_2, paste0("es_2_", stock_name, ".csv"))
    
    print(xtable(round(sigma2_fit_n, 5)), file = paste0("Table_sigma2_fit_normal_", stock_name, ".tex"), compress = FALSE)
    print(xtable(round(sigma2_fit_t, 5)), file = paste0("Table_sigma2_fit_t_student_", stock_name, ".tex"), compress = FALSE)
    
    
    print(xtable(round(var_1, 5)), file = paste0("Table_var_1_", stock_name, ".tex"), compress = FALSE)
    print(xtable(round(var_2, 5)), file = paste0("Table_var_2_", stock_name, ".tex"), compress = FALSE)
    print(xtable(round(es_1, 5)), file =  paste0("Table_es_1_", stock_name, ".tex"), compress = FALSE)
    print(xtable(round(es_2, 5)), file = paste0("Table_es_2_", stock_name, ".tex"), compress = FALSE)
    
  }, warning = function(warning_condition) {
    print(stock_name)
  }, error = function(error_condition) {
    print(stock_name)
  })
}










