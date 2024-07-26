#include<RcppArmadillo.h>
#include<Rmath.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP rbgarch_loglik(NumericVector theta, NumericVector r, NumericVector proxi){
  int n = r.size();
  NumericVector h(n), log_lik(n - 1);
  Rcpp::Function sum("sum");
  
  h[0] = theta[0] / (1 - theta[1] - theta[2]);
  for(int i = 1; i < n; i++){
    h[i] = theta[0] + theta[1]*proxi[i-1]+ theta[2]*h[i-1];
    log_lik[i - 1] = -0.5*pow(r[i], 2) / h[i] - 0.5* log(h[i]);
  }
  return Rcpp::wrap(sum(-log_lik));
}




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP rbgarch_loglik_t_student(NumericVector theta, NumericVector r, NumericVector proxi){
  int n = r.size();
  double pi = 3.141592653589793238462643383280;
  NumericVector h(n), log_lik(n - 1);
  Rcpp::Function sum("sum");
  
  h[0] = theta[0] / (1 - theta[1] - theta[2]);
  for(int i = 1; i < n; i++){
    h[i] = theta[0] + theta[1]*proxi[i-1]+ theta[2]*h[i-1];
    log_lik[i - 1] =  log(std::tgamma((theta[3] + 1.0)/2.0)) - log(std::tgamma(theta[3]/2))- 0.5*log(theta[3]*pi)- 0.5*log(h[i]) - ((theta[3] + 1.0 )/2.0)*log(1+(1+pow(r[i], 2.0)/h[i]));
    
  }
  return Rcpp::wrap(sum(-log_lik));
}






// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP rbgjr_loglik(NumericVector theta, NumericVector r, NumericVector proxi){
  int n = r.size();
  NumericVector h(n), log_lik(n - 1);
  Rcpp::Function sum("sum");
  
  h[0] = theta[0] / (1 - theta[1] - theta[2] - 0.5* theta[3]);
  for(int i = 1; i < n; i++){
    if (r[i - 1] < 0) {
      h[i] = theta[0] + theta[1] * proxi[i - 1] + theta[2] * h[i - 1] + theta[3] * proxi[i - 1];
    } else{
      h[i] = theta[0] + theta[1] * proxi[i - 1] + theta[2] * h[i - 1];
    }
    log_lik[i - 1] = -0.5*pow(r[i], 2) / h[i] - 0.5* log(h[i]);
  }
  return Rcpp::wrap(sum(-log_lik));
}




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP rbgjr_loglik_t_student(NumericVector theta, NumericVector r, NumericVector proxi){
  int n = r.size();
  double pi = 3.141592653589793238462643383280;
  NumericVector h(n), log_lik(n - 1);
  Rcpp::Function sum("sum");
  
  h[0] = theta[0] / (1 - theta[1] - theta[2] - 0.5* theta[3]);
  for(int i = 1; i < n; i++){
    if (r[i - 1] < 0) {
      h[i] = theta[0] + theta[1] * proxi[i - 1] + theta[2] * h[i - 1] + theta[3] * proxi[i - 1];
    } else{
      h[i] = theta[0] + theta[1] * proxi[i - 1] + theta[2] * h[i - 1];
    }
    log_lik[i - 1] = log(std::tgamma((theta[4] + 1.0)/2.0)) - log(std::tgamma(theta[4]/2))- 0.5*log(theta[4]*pi)- 0.5*log(h[i]) - ((theta[4] + 1.0 )/2.0)*log(1+(1+pow(r[i], 2.0)/h[i]));
  }
  return Rcpp::wrap(sum(-log_lik));
}






// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP grid_rbgjr(NumericVector y, NumericVector proxi){
  NumericVector coeff(4), vi(4);
  double omega, alpha, beta, gama;
  double omega_min = 0.00001, omega_max= 0.05, alpha_min = 0.005,  alpha_max = 0.2, beta_min = 0.65,  beta_max = 0.98,  gama_min = 0.001, gama_max = 0.4, n_omega = 5, n_alpha = 5,  n_beta = 5, n_gama = 5;
  double ml = 100000000, nml;
  double lm_omega = (omega_max - omega_min)/n_omega;
  double lm_alpha = (alpha_max - alpha_min)/n_alpha;
  double lm_beta = (beta_max - beta_min)/n_beta;
  double lm_gama = (gama_max - gama_min)/n_gama;
  Rcpp::Function rbgjr_loglik("rbgjr_loglik");
  
  for(int nj = 0; nj < n_alpha; nj++){
    for(int nk = 0; nk < n_beta; nk++){
      for(int ni =0; ni < n_omega; ni++){
        for(int nl =0; nl < n_gama; nl++){
          alpha = alpha_min + nj*lm_alpha;
          beta = beta_min + nk*lm_beta; 
          omega = omega_min +ni*lm_omega;
          gama = gama_min + nl*lm_gama;
          if(gama < 2*(1 - alpha- beta)){
            coeff[0] = omega;
            coeff[1] = alpha;
            coeff[2] = beta;
            coeff[3] = gama;
            nml = Rcpp::as<double>(rbgjr_loglik(coeff, y, proxi));
            if (nml < ml){
              vi[0] = coeff[0];
              vi[1] = coeff[1];
              vi[2] = coeff[2];
              vi[3] = coeff[3];
              ml=nml;
            }
          }
        }
      }
    }
  }
  return(vi);
}





// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP grid_rbgjr_t_student(NumericVector y, NumericVector proxi){
  NumericVector coeff(5), vi(5);
  double omega, alpha, beta, gama, tau;
  double omega_min = 0.00001, omega_max= 0.05, alpha_min = 0.005,  alpha_max = 0.2, beta_min = 0.65,  beta_max = 0.98,  gama_min = 0.001, gama_max = 0.4, tau_min = 3.0, tau_max = 20.0, n_omega = 5, n_alpha = 5,  n_beta = 5, n_gama = 5, n_tau = 5;
  double ml = 100000000, nml;
  double lm_omega = (omega_max - omega_min)/n_omega;
  double lm_alpha = (alpha_max - alpha_min)/n_alpha;
  double lm_beta = (beta_max - beta_min)/n_beta;
  double lm_gama = (gama_max - gama_min)/n_gama;
  double lm_tau = (tau_max - tau_min)/n_tau;
  Rcpp::Function rbgjr_loglik_t_student("rbgarch_loglik_t_student");
  
  for(int nj = 0; nj < n_alpha; nj++){
    for(int nk = 0; nk < n_beta; nk++){
      for(int ni =0; ni < n_omega; ni++){
        for(int nl =0; nl < n_gama; nl++){
          for(int nh =0; nh < n_tau; nh++){
          alpha = alpha_min + nj*lm_alpha;
          beta = beta_min + nk*lm_beta; 
          omega = omega_min +ni*lm_omega;
          gama = gama_min + nl*lm_gama;
          tau = tau_min + nh*lm_tau;
          if(gama < 2*(1 - alpha- beta)){
            coeff[0] = omega;
            coeff[1] = alpha;
            coeff[2] = beta;
            coeff[3] = gama;
            coeff[4] = tau;
            nml = Rcpp::as<double>(rbgjr_loglik_t_student(coeff, y, proxi));
            if (nml < ml){
              vi[0] = coeff[0];
              vi[1] = coeff[1];
              vi[2] = coeff[2];
              vi[3] = coeff[3];
              vi[4] = coeff[4];
              ml=nml;
            }
          }
         }
        }
      }
    }
  }
  return(vi);
}






// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP grid_rbgarch(NumericVector y, NumericVector proxi){
  NumericVector coeff(3), vi(3);
  double omega, alpha, beta;
  double omega_min = 0.00001, omega_max= 0.05, alpha_min = 0.005,  alpha_max = 0.2, beta_min = 0.65,  beta_max = 0.98,  n_omega = 5, n_alpha = 5,  n_beta = 5;
  double ml = 100000000, nml;
  double lm_omega = (omega_max - omega_min)/n_omega;
  double lm_alpha = (alpha_max - alpha_min)/n_alpha;
  double lm_beta = (beta_max - beta_min)/n_beta;
  Rcpp::Function rbgarch_loglik("rbgarch_loglik");
   
  for(int nj = 0; nj < n_alpha; nj++){
    for(int nk = 0; nk < n_beta; nk++){
      for(int ni =0; ni < n_omega; ni++){
        alpha = alpha_min + nj*lm_alpha;
        beta = beta_min + nk*lm_beta; 
        omega = omega_min +ni*lm_omega;
        if(alpha+ beta < 0.999){
         coeff[0] = omega;
         coeff[1] = alpha;
         coeff[2] = beta;
         nml = Rcpp::as<double>(rbgarch_loglik(coeff, y, proxi));
         if (nml < ml){
           vi[0] = coeff[0];
           vi[1] = coeff[1];
           vi[2] = coeff[2];
           ml=nml;
         }
       }
     }
   }
  }
  return(vi);
}




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP grid_rbgarch_t_student(NumericVector y, NumericVector proxi){
  NumericVector coeff(4), vi(4);
  double omega, alpha, beta, tau;
  double omega_min = 0.00001, omega_max= 0.05, alpha_min = 0.005,  alpha_max = 0.2, beta_min = 0.65,  beta_max = 0.98, tau_min = 3.0,  tau_max = 20.0,  n_omega = 5, n_alpha = 5,  n_beta = 5, n_tau = 5;
  double ml = 100000000, nml;
  double lm_omega = (omega_max - omega_min)/n_omega;
  double lm_alpha = (alpha_max - alpha_min)/n_alpha;
  double lm_beta = (beta_max - beta_min)/n_beta;
  double lm_tau = (tau_max - tau_min)/n_tau;
  Rcpp::Function rbgarch_loglik_t_student("rbgarch_loglik_t_student");
  
  for(int nj = 0; nj < n_alpha; nj++){
    for(int nk = 0; nk < n_beta; nk++){
      for(int ni =0; ni < n_omega; ni++){
        for(int nh =0; nh < n_tau; nh++){
        alpha = alpha_min + nj*lm_alpha;
        beta = beta_min + nk*lm_beta; 
        omega = omega_min +ni*lm_omega;
        tau = tau_min +nh*lm_tau;
        if(alpha+ beta < 0.999){
          coeff[0] = omega;
          coeff[1] = alpha;
          coeff[2] = beta;
          coeff[3] = tau;
          nml = Rcpp::as<double>(rbgarch_loglik_t_student(coeff, y, proxi));
          if (nml < ml){
            vi[0] = coeff[0];
            vi[1] = coeff[1];
            vi[2] = coeff[2];
            vi[3] = coeff[3];
            ml=nml;
          }
        }
      }
     }
    }
  }
  return(vi);
}


