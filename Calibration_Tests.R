
#########################################
### Testes de calibracao ################
#########################################
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
source("Utils.R")

ibrx_names <- names(read.csv("../Range-Based-Volatility/Data/ibrx_names.csv"))
alpha1 <- 0.010
alpha2 <- 0.025


for (stock_name in ibrx_names) {
  files_names <- list.files(path = '/home/ctrucios/Dropbox/Students/MSc/Mestrado-Adi/Range-Based-Volatility/Results', pattern = stock_name, full.names = T)
  if (length(files_names) == 7) {
    es_1 <- read.csv(files_names[1])[, -1]
    es_2 <- read.csv(files_names[2])[, -1]
    ret_oos <- read.csv(files_names[3])[, -1]
    sigma2_normal <- read.csv(files_names[4])[, -1]
    sigma2_t <- read.csv(files_names[5])[, -1]
    var_1 <- read.csv(files_names[6])[, -1]
    var_2 <- read.csv(files_names[7])[, -1]
    sigma2 <- cbind(sigma2_normal[, 1:10], sigma2_t[, 1:10])
  
    # QLIKE
    comparison <- matrix(NA, ncol = 20, nrow = 2)
    for (i in 1:20) {
      comparison[1, i] <- 100*round(loss_qlike(sigma2[, i] + 0.00001, ret_oos^2+ 0.00001), 5)
      comparison[2, i] <- round(1000000*loss_mse(sigma2[, i], ret_oos^2), 5)
    }
    write.csv(comparison, paste0("QLIKE_", stock_name, ".csv"))
    print(xtable(round(comparison, 4)), file = paste0("Table_QLIKE_", stock_name, ".tex"), compress = FALSE)
    
  
    # Tests
    nomes_testes <- c("Hits", "UC", "CC",  "VQ", "ER", "CoC", "ESR_1", "ESR_2", "ESR_3")
    BackVaRES_1   = matrix(NA, ncol = 9, nrow = 20, dimnames = list(colnames(es_1), nomes_testes))
    BackVaRES_2.5 = matrix(NA, ncol = 9, nrow = 20, dimnames = list(colnames(es_1), nomes_testes))
  
    for (i in 1:20) { 
      print(i)
    
      BackT_1   = BacktestVaR(ret_oos, var_1[ ,i], alpha = alpha1)
      BackT_2.5 = BacktestVaR(ret_oos, var_2[ ,i], alpha = alpha2)
    
      EBackT_1   = ESTest(alpha = alpha1, ret_oos, es_1[ ,i], var_1[ ,i], conf.level = 0.95,  boot = TRUE)
      EBackT_2.5 = ESTest(alpha = alpha2, ret_oos, es_2[ ,i], var_2[ ,i], conf.level = 0.95,  boot = TRUE)
    
    
      BackVaRES_1[i,] = c(mean(ret_oos < var_1[ ,i])*100, 
                        BackT_1$LRuc[2], BackT_1$LRcc[2], VaR_VQR(ret_oos, var_1[ ,i], alpha1),
                        EBackT_1$boot.p.value,
                        cc_backtest(ret_oos, var_1[ ,i], es_1[ ,i], alpha  = alpha1)$pvalue_twosided_simple, 
                        esr_backtest(ret_oos, var_1[ ,i], es_1[ ,i],alpha  = alpha1, B = 0, version = 1)$pvalue_twosided_asymptotic,
                        esr_backtest(ret_oos, var_1[ ,i], es_1[ ,i],alpha  = alpha1, B = 0, version = 2)$pvalue_twosided_asymptotic,
                        esr_backtest(ret_oos, var_1[ ,i], es_1[ ,i],alpha  = alpha1, B = 0, version = 3)$pvalue_twosided_asymptotic)
    
    
      BackVaRES_2.5[i,] = c(mean(ret_oos < var_2[ ,i])*100,
                          BackT_2.5$LRuc[2], BackT_2.5$LRcc[2], VaR_VQR(ret_oos, var_2[ ,i], alpha2),
                          EBackT_2.5$boot.p.value,
                          cc_backtest(ret_oos, var_2[ ,i], es_2[ ,i], alpha  = alpha2)$pvalue_twosided_simple,
                          esr_backtest(ret_oos, var_2[ ,i], es_2[ ,i],alpha  = alpha2, B = 0, version = 1)$pvalue_twosided_asymptotic,
                          esr_backtest(ret_oos, var_2[ ,i], es_2[ ,i],alpha  = alpha2, B = 0, version = 2)$pvalue_twosided_asymptotic,
                          esr_backtest(ret_oos, var_2[ ,i], es_2[ ,i],alpha  = alpha2, B = 0, version = 3)$pvalue_twosided_asymptotic)
    
    
    }
    print(xtable(round(BackVaRES_1, 4)), file = paste0("Table_VaRES_1_", stock_name, ".tex"), compress = FALSE)
    print(xtable(round(BackVaRES_2.5, 4)), file = paste0("Table_VaRES_2_", stock_name, ".tex"), compress = FALSE)
  }
}





