

#' Simulation function 
#'
#' @param  sim.config list contains following parameters:
#'         I  num The number of patients.
#'         k  num The number of ROIs in each patient.
#'         R num The number of readers.
#'         correlation_t num The correlation for simulating truth label.
#'         potential_correlation_s num The correlation for simulating reading scores.
#'         AUC_all num The theoretical AUC values.
#'         sameclustersize boolean The binary variable to decide whether we have same number of ROIs in each patient.
#'         rho num The scale parameter that infulence the covariance matrix in multivariate normal distribution.
#'         fix_design boolean Binary variable to decide whether fix the truth label in simulation.
#'         stream num The integer control the random number generator.
#'
#' @return A list and the only element in the list is the simulated data with following columns:
#'         "clusterID","unitID","reader1",...,"truth"
#' 
#' @export
#'
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats qnorm
#' 

data_MRMC <- function(sim.config) {

  
  
  I = sim.config$case_num 
  k = sim.config$ROI_num
  R = sim.config$Reader_num
  correlation_t = sim.config$Truth_corr
  potential_correlation_s = sim.config$ROI_corr
  AUC_all = sim.config$AUC
  sameclustersize = sim.config$same_size_flag
  rho = sim.config$corr_factor
  fix_design = sim.config$design_flag
  stream = sim.config$stream_num
  initial_seed  =  sim.config$seed
  
  
  mu <- rep(0,k)
  covariance_t <- matrix(rep(correlation_t, k*k),k,k) + diag(x = (1 - correlation_t),nrow = k,ncol = k)
  if (fix_design == TRUE) {
  set.seed(initial_seed)
  truth_score <- rmvnorm(I, mean = mu, sigma = covariance_t)
  iMRMC::init.lecuyerRNG(seed = 1, stream)
  }
  else {
  iMRMC::init.lecuyerRNG(seed = 1, stream)
  truth_score <- rmvnorm(I, mean = mu, sigma = covariance_t)
  }
  #set.seed(s_seed)
  colnames(truth_score) <- as.factor(paste("unit",rep(1:k),sep = ""))
  rownames(truth_score) <- as.factor(paste("cluster",rep(1:I),sep = ""))
  
  #dim(truth_score) I k
  truth <- truth_score > 0
  truth <- matrix(as.numeric(truth),ncol = k)
  colnames(truth) <- as.factor(paste("unit",rep(1:k),sep = ""))
  rownames(truth) <- as.factor(paste("cluster",rep(1:I),sep = ""))
  
  testscore <- matrix(0,ncol = k*R,nrow = I)
  for (i in 1:I) {
    index <- sample(1:4, 1,replace = TRUE)
    correlation_s <- potential_correlation_s[index]
    
    covariance_s <- NULL
    for (r in 1:R) {
      covariance_A <- matrix(rep(correlation_s, k*k),k,k) + diag(x = (1 - correlation_s),nrow = k,ncol = k)
      covariance_AB <- matrix(rep(correlation_s*rho, k*k),k,k) + diag(x = (rho - correlation_s*rho),nrow = k,ncol = k)
      
      covariance_temp <- NULL
      if (r == 1) {
        covariance_temp <- cbind(covariance_temp,covariance_A)
        for (j in (r + 1):R) {
          covariance_temp <- cbind(covariance_temp,covariance_AB)
        }
      }else if (r == R) {
        for (j in 1:(r - 1)) {
          covariance_temp <- cbind(covariance_temp,covariance_AB)
        }
        covariance_temp <- cbind(covariance_temp,covariance_A)
      }else {
        for (j in 1:(r - 1)) {
          covariance_temp <- cbind(covariance_temp,covariance_AB)
        }
        covariance_temp <- cbind(covariance_temp,covariance_A)
        for (j in (r + 1):R) {
          covariance_temp <- cbind(covariance_temp,covariance_AB)
        }
      }
      covariance_s <- rbind(covariance_s,covariance_temp)
    }
    
    mu <- rep(0,k*R)
    testscore[i,] <- rmvnorm(1, mean = mu, sigma = covariance_s)
  }
  
  colnames(testscore) <- as.factor(paste(rep(paste("reader",rep(1:R),sep = ""),each = k), rep(paste("unit",rep(1:k),sep = ""),R)))
  rownames(testscore) <- as.factor(paste("cluster",rep(1:I),sep = ""))
  
  data = as.data.frame(matrix(0,nrow = k*I,ncol = 3 + R))
  colnames(data) = c("clusterID","unitID",as.vector(as.factor(paste("reader",rep(1:R),sep = ""))),"truth")
  
  data$clusterID = as.factor(paste("cluster",rep(1:I,each = k),sep = ""))
  data$unitID = as.factor(paste("unit",rep(1:k,I),sep = ""))
  for (r in 1:R) {
    data[,2 + r] = c(as.vector(t(testscore[,(1 + (r - 1)*k):(k + (r - 1)*k)])))
  }
  data$truth = c(as.vector(t(truth)))
  
  if (sameclustersize == FALSE) {
    delete_index <- sample(1:I*k, I*k/10,replace = FALSE)
    data <- data[-delete_index,]
  }
  
  index_positive <- data$truth > 0
  delta_all <- NULL
  data_final <- data
  for (r in 1:R) {
    delta <- qnorm(AUC_all[r],0,1)*sqrt(2) 
    data_final[index_positive,(2 + r)] <- data[index_positive,(2 + r)] + delta
    delta_all <- cbind(delta_all,delta)
  }
  
  data_final$mod = rep(0,I*k)
  data_final = data_final[,c(1,3:ncol(data_final),2)]
  colnames(data_final) = c("patient","reader1","reader2","truth","mod","region")
  
  return(list("data_final" = data_final))
}


#' Configuration function 
#'
#' @param I The number of patients.
#' @param k The number of ROIs in each patient.
#' @param R The number of readers.
#' @param correlation_t The correlation for simulating truth label.
#' @param potential_correlation_s The correlation for simulating reading scores.
#' @param AUC_all The theoretical AUC values.
#' @param sameclustersize The binary variable to decide whether we have same number of ROIs in each patient.
#' @param rho The scale parameter that influence the covariance matrix in multivariate normal distribution.
#' @param fix_design Binary variable to decide whether fix the truth label in simulation.
#' @param stream The integer control the random number generator.
#' @param initial_seed The integer control the random seed for truth label generation.
#'
#' @return A list of above parameters
#' 
#' @export
#'
#' @importFrom mvtnorm rmvnorm
simu_config = function(I = 100, k = 10, R = 2, 
                       correlation_t = 0, 
                       potential_correlation_s = rep(0.5,4),
                       AUC_all = rep(0.7, 2), 
                       sameclustersize = TRUE, rho = 0.5, 
                       fix_design = FALSE, stream = 20220210, initial_seed = 20220222) {
  return(
    list(
      case_num = I,
      ROI_num = k,
      Reader_num = R,
      Truth_corr = correlation_t,
      ROI_corr = potential_correlation_s,
      AUC = AUC_all,
      same_size_flag = sameclustersize,
      corr_factor = rho,
      design_flag = fix_design,
      stream_num = stream,
      seed = initial_seed
    )
  )
}


#' The test demo data to be included in my package
#'
#' @name expected_data
#' @docType data
#' @author Hongfei Du \email{hongfei@gwu.edu}
#' @keywords data
NULL