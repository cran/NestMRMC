

#' Calculate the between-cases AUC estimator's theoretical variance and covariance
#' 
#' @description This function calculates between-cases AUC estimator's theoretical variance and covariance based on all the truths, namely, the ROI's truth labels, AUC values, covariance between ROI scores 
#' within same reader, scale factor that influences the covariance between ROI scores between readers and the variances for 
#' positive and negative ROI scores. Detailed formulas are available in following paper: 
#' Single Reader Between-Cases AUC Estimator with Nested Data. Statistical Methods in Medical Research. https://doi.org/10.1177/09622802221111539.
#' There is also a Rcpp version of this function in this package. The function name is 'true_AUC_var_abitrary_Rcpp', which is much faster than current version. They produce the exact same results.
#' @param numROI The number of positive and negative ROIs in all the patients.
#' @param AUC The AUC values used in simulated data.
#' @param cov The covariance used in simulating reading scores.
#' @param rho The scale factor used in simulating reading scores.
#' @param sigma_pos The variacne for positive ROI's reading score, defalut is 1.
#' @param sigma_neg The variacne for negative ROI's reading score, defalut is 1.
#'
#' @return The theoretical AUC estimator's (co)variance based on the simulation settings.
#' @export
#'
#' @importFrom mvtnorm pmvnorm
true_AUC_var_abitrary = function(numROI, AUC = 0.7, cov = 0.5,
                                 rho = 0.5, sigma_pos = 1, sigma_neg = 1){
 
  # reading the number of cases 
  ncase = ncol(numROI)
  
  # initialize the coefficients 
  c1 = 0
  c2 = 0
  c3 = 0
  c4 = 0
  c7 = 0
  c5 = 0
  c6 = 0
  c8 = 0
  c9 = 0
  c10 = 0
  c11 = 0
  
  
  # Using for loop to calculate all the coefficients
  for (i in 1:ncase) {
    for (j in 1:ncase) {
      if (j == i) {next}
      else {
        
        c1 = c1 + numROI[1,i]*numROI[2,j]
        c2 = c2 + numROI[1,i]*numROI[2,j]*(numROI[1,i] - 1)*(numROI[2,j] - 1)
        c3 = c3 + numROI[1,i]*numROI[2,j]*(numROI[1,i] - 1)
        c4 = c4 + numROI[1,i]*numROI[2,j]*(numROI[2,j] - 1)
        c7 = c7 + numROI[1,i]*numROI[2,j]*numROI[1,j]*numROI[2,i]
        
        for (h in 1:ncase) {
          
          if (h == j | h == i) {next}
          else{
            c5 = c5 + numROI[1,i]*numROI[2,j]*numROI[2,h]
            c6 = c6 + numROI[1,i]*(numROI[1,i] - 1)*numROI[2,j]*numROI[2,h]
            c8 = c8 + numROI[1,i]*numROI[1,j]*numROI[2,h]*numROI[2,i]
            c9 = c9 + numROI[1,i]*numROI[1,j]*numROI[2,j]*numROI[2,h]
            c10 = c10 + numROI[1,i]*numROI[1,j]*numROI[2,h]
            c11 = c11 + numROI[1,i]*numROI[1,j]*numROI[2,h]*(numROI[2,h] - 1)
            
          }
        }
      }
    }
  }
  
  # Calculate the last coefficient 
  c12 = (sum(numROI[1,])*sum(numROI[2,]) - sum(numROI[1,]*numROI[2,]))^2 - 
    (c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + c10 + c11)
  
  # Calculate the mean difference delta 
  delta = qnorm(AUC)*sqrt(sigma_pos+sigma_neg)
  
  # assign values
  sigma = cov
  pos_sig = sigma_pos
  neg_sig = sigma_neg
  
  ## Calculate all the theoretical moments' value for AUC variance 
  phi_sigma = pmvnorm(lower = c(0,0),sigma = matrix(c(pos_sig+neg_sig,sigma,sigma,pos_sig+neg_sig),2,2,byrow = T),mean = c(delta,delta))
  phi_negsigma = pmvnorm(lower = c(0,0),sigma = matrix(c(pos_sig+neg_sig,-sigma,-sigma,pos_sig+neg_sig),2,2,byrow = T),mean = c(delta,delta))
  phi_2sigma = pmvnorm(lower = c(0,0),sigma = matrix(c(pos_sig+neg_sig,2*sigma,2*sigma,pos_sig+neg_sig),2,2,byrow = T),mean = c(delta,delta))
  phi_neg2sigma = pmvnorm(lower = c(0,0),sigma = matrix(c(pos_sig+neg_sig,-2*sigma,-2*sigma,pos_sig+neg_sig),2,2,byrow = T),mean = c(delta,delta))
  phi_sig_neg = pmvnorm(lower = c(0,0),sigma = matrix(c(pos_sig+neg_sig,sigma+neg_sig,sigma+neg_sig,pos_sig+neg_sig),2,2,byrow = T),mean = c(delta,delta))
  phi_sig_pos = pmvnorm(lower = c(0,0),sigma = matrix(c(pos_sig+neg_sig,sigma+pos_sig,sigma+pos_sig,pos_sig+neg_sig),2,2,byrow = T),mean = c(delta,delta))
  phi_possig = pmvnorm(lower = c(0,0),sigma = matrix(c(pos_sig+neg_sig,pos_sig,pos_sig,pos_sig+neg_sig),2,2,byrow = T),mean = c(delta,delta))
  phi_negsig = pmvnorm(lower = c(0,0),sigma = matrix(c(pos_sig+neg_sig,neg_sig,neg_sig,pos_sig+neg_sig),2,2,byrow = T),mean = c(delta,delta)) 
  
  ## combine all the moments to calculate the AUC true variance
  true_var = 1/(c1+c2+c3+c4+c5+c6+c7+c8+c9+c10+c11+c12)*(c1*AUC +
                                                           c2*phi_2sigma  +
                                                           c3*phi_sig_neg +
                                                           c4*phi_sig_pos +
                                                           c5*phi_possig +
                                                           c6*phi_sigma + 
                                                           c7*phi_neg2sigma + 
                                                           c8*phi_negsigma +  
                                                           c9*phi_negsigma +     
                                                           c10*phi_negsig +
                                                           c11*phi_sigma +
                                                           c12*AUC^2) - AUC^2
  
  
  ## Calculate all the theoretical moments' value for AUC covariance
  phi_2rho = pmvnorm(lower = c(0,0),sigma = matrix(c(pos_sig+neg_sig,2*rho,2*rho,pos_sig+neg_sig),2,2,byrow = T),mean = c(delta,delta))
  phi_2rhosigma = pmvnorm(lower = c(0,0),sigma = matrix(c(pos_sig+neg_sig,2*rho*sigma,2*rho*sigma,pos_sig+neg_sig),2,2,byrow = T),mean = c(delta,delta))
  phi_neg2rhosigma = pmvnorm(lower = c(0,0),sigma = matrix(c(pos_sig+neg_sig,-2*rho*sigma,-2*rho*sigma,pos_sig+neg_sig),2,2,byrow = T),mean = c(delta,delta))
  phi_rhosig_rho= pmvnorm(lower = c(0,0),sigma = matrix(c(pos_sig+neg_sig,rho*sigma+rho,rho*sigma+rho,pos_sig+neg_sig),2,2,byrow = T),mean = c(delta,delta))
  phi_rho = pmvnorm(lower = c(0,0),sigma = matrix(c(pos_sig+neg_sig,rho,rho,pos_sig+neg_sig),2,2,byrow = T),mean = c(delta,delta))
  phi_rhosig = pmvnorm(lower = c(0,0),sigma = matrix(c(pos_sig+neg_sig,rho*sigma,rho*sigma,pos_sig+neg_sig),2,2,byrow = T),mean = c(delta,delta))
  phi_negrhosig = pmvnorm(lower = c(0,0),sigma = matrix(c(pos_sig+neg_sig,-rho*sigma,-rho*sigma,pos_sig+neg_sig),2,2,byrow = T),mean = c(delta,delta))
  
  ## combine all the moments to calculate the AUC true covariance
  true_cov = 1/(c1+c2+c3+c4+c5+c6+c7+c8+c9+c10+c11+c12)*(c1*phi_2rho +
                                                           c2*phi_2rhosigma  +
                                                           c3*phi_rhosig_rho +
                                                           c4*phi_rhosig_rho +
                                                           c5*phi_rho +
                                                           c6*phi_rhosig + 
                                                           c7*phi_neg2rhosigma + 
                                                           c8*phi_negrhosig +  
                                                           c9*phi_negrhosig +     
                                                           c10*phi_rho +
                                                           c11*phi_rhosig +
                                                         c12*AUC^2) - AUC^2

  # return the variance and covariance
  return(c(true_var,true_cov))
  
}

