#' MRMC analysis in nested data problem
#'
#'
#' @description This function takes nested data as a data frame and runs a multi-reader multi-case analysis for single reader in nested
#' data problem based on modified U-statistics as described in the following paper:
#'
#'
#'
#'
#'
#'
#' @param data The nested data for analysis. This dataset should have specified columns:
#'             "patient","reader1","reader2","reader3","reader4","reader5","truth","mod","region".
#'
#' @return This function returns a [list] containing three dataframes.
#'
#'
#'         Here is a quick summary:
#'
#'             AUC_per_reader [data.frame] this data frame contains the AUC estimates for each reader under different modalities (Mod1 denotes modality 1 and Mod2 denotes modality 2). 
#'
#'             AUC_Var_per_reader [data.frame] this data frame contains the AUC variance estimates for each reader under different modalities.
#'             
#'             numROI [data.frame] this data frame contains the number of positive and negative ROIs in each case. 
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr do
#' @examples 
#' 
#'  
#' data = NestMRMC::expected_data
#' 
#' Outputs = AUC_per_reader_nest(data)
#'
#'
AUC_per_reader_nest = function(data){
  ## Calculate all the perfect squares and counts for each perfect square
  temp_output = success_score(data)
  success_array = temp_output[[1]]
  numROI = temp_output[[2]]
    
  nr = dim(success_array)[2]
  nm = dim(success_array)[1]
  ## 1st perfect square
  m1 = apply(success_array[,,,,1],c(1,2),sum)^2
  c1 = apply(success_array[,,,,2],c(1,2),sum)^2
  

  ## 2nd perfect square
  m2 = apply(success_array[,,,,1],c(1,2,3),sum)^2 %>% apply(c(1,2),sum)
  c2 = apply(success_array[,,,,2],c(1,2,3),sum)^2 %>% apply(c(1,2),sum)
  ## 3rd perfect square
  m3 = apply(success_array[,,,,1],c(1,2),sum_diag)*apply(success_array[,,,,1],c(1,2),sum)
  c3 = apply(success_array[,,,,2],c(1,2),sum_diag)*apply(success_array[,,,,2],c(1,2),sum)
  ## 4th perfect square
  m4 = ((apply(success_array[,,,,1],c(1,2),diag) %>% aperm(c(2,3,1)))*apply(success_array[,,,,1],c(1,2,3),sum)) %>% apply(c(1,2),sum)
  c4 = ((apply(success_array[,,,,2],c(1,2),diag) %>% aperm(c(2,3,1)))*apply(success_array[,,,,2],c(1,2,3),sum)) %>% apply(c(1,2),sum)
  ## 5th perfect square
  m5 = (apply(success_array[,,,,1],c(1,2,4),sum)*apply(success_array[,,,,1],c(1,2,3),sum)) %>% apply(c(1,2),sum)
  c5 = (apply(success_array[,,,,2],c(1,2,4),sum)*apply(success_array[,,,,2],c(1,2,3),sum)) %>% apply(c(1,2),sum)
  ## 6th perfect square
  m6 = ((
    apply(
    success_array[,,,,1],c(1,2),diag
    ) %>% 
    aperm(c(2,3,1))
    )*apply(success_array[,,,,1],c(1,2,4),sum)) %>% 
    apply(c(1,2),sum)
  
  c6 = ((apply(success_array[,,,,2],c(1,2),diag) %>% aperm(c(2,3,1)))*apply(success_array[,,,,2],c(1,2,4),sum)) %>% apply(c(1,2),sum)
  ## 7th perfect square
  m7 = (apply(success_array[,,,,1],c(1,2),diag) %>% aperm(c(2,3,1)))^2 %>% apply(c(1,2),sum)
  c7 = (apply(success_array[,,,,2],c(1,2),diag) %>% aperm(c(2,3,1)))^2 %>% apply(c(1,2),sum)
  ## 8th perfect square
  m8 = apply(success_array[,,,,1],c(1,2),m8_f)
  c8 = apply(success_array[,,,,2],c(1,2),m8_f)
  ## 9th perfect square

  m9 = apply(success_array[,,,,1],c(1,2),sum_diag)^2
  c9 = apply(success_array[,,,,2],c(1,2),sum_diag)^2
  ## 10th perfect square
  m10 = apply(success_array[,,,,1],c(1,2,4),sum)^2 %>% apply(c(1,2),sum)
  c10 = apply(success_array[,,,,2],c(1,2,4),sum)^2 %>% apply(c(1,2),sum)
  ## 11th perfect square
  m11 = apply(success_array[,,,,1],c(1,2),m11_f)
  c11 = apply(success_array[,,,,2],c(1,2),m11_f)

  
  ## change m1-m11 to m_112 something in the paper.
  
  ## Combine perfect squares
  first_square = 
    (m1 - m2 - 2*m3 + 4*m4 - 2*m5 + 4*m6 - 6*m7 + m8 + m9 - m10 + m11) /
    (c1 - c2 - 2*c3 + 4*c4 - 2*c5 + 4*c6 - 6*c7 + c8 + c9 - c10 + c11)
     
  
  
  ## Calculate the each reader's AUC estimate
  AUC_per_reader = apply(success_array[,,,,1],c(1,2),delete_diag)/apply(success_array[,,,,2],c(1,2),delete_diag)
  AUC_per_reader = as.data.frame(AUC_per_reader)
  
  ## Name the row and column for the output dataframe
  colnames(AUC_per_reader) = paste0("reader",1:nr)
  rownames(AUC_per_reader) = paste0("mod",1:nm)
  
  ## Calculate the each reader's AUC variance estimate
  AUC_Var_per_reader = AUC_per_reader^2 - first_square

  AUC_Var_per_reader = as.data.frame(AUC_Var_per_reader)
  
  ## Name the row and column for the output dataframe
  colnames(AUC_Var_per_reader) = paste0("reader",1:nr)
  rownames(AUC_Var_per_reader) = paste0("mod",1:nm)
  
  ## Calculate the AUC covariance between two readers
  AUC_Var_per_reader$cov = AUC_cov_2reader_nest(temp_output)
  
  ## Return the final outputs 
  return(list(AUC = AUC_per_reader,
              Var = AUC_Var_per_reader,  
              numROI = numROI))
}

## function to generate the success score and corresponding counts under different 
## modalities and readers 

#' Calculate the success score
#'
#' @param data the nested MRMC data
#'
#' @return The success score and number of ROIs in each case

success_score = function(data){

  np = length(unique(data$patient))
  nr = colnames(data) %>% grep(pattern = "reader") %>% length()
  #nROI = data$region %>% unique() %>% length()
  nm = data$mod %>% unique() %>% length()
  success_array = array(0,dim = c(2,nr,np,np,2))
  ## create index for all combination pairs
  comb_pairs = cbind(rep(1:np,each = np),rep(1:np,np))
  readers_name = colnames(data)[colnames(data) %>% grep(pattern = "reader")]
  mods_name = unique(data$mod)
  for (m in 1:nm) {

    for (r in 1:nr) {


     
      temp_data = data[,c("patient",readers_name[r],"truth","mod","region")] %>% 
        filter(data$mod == mods_name[m])
      
      
      colnames(temp_data)[2] = "reader"

      ## transfer data to two list contain np vectors each

     
      
      temp1 = split(temp_data, temp_data$patient)
      truth_list = lapply(temp1, function(x) t(x$truth))
      score_list = lapply(temp1, function(x) t(x$reader))
      
     
      
      
      
      num_posROI = sapply(truth_list,FUN = sum)
      num_negROI = sapply(truth_list,FUN = length) - num_posROI
      num_ROI = rbind(num_posROI,num_negROI)
      sum_comb = function(index){

        a = index[1]
        b = index[2]
        ma = truth_list[[a]] %>% sum()
        nb = truth_list[[b]] %>% length() - truth_list[[b]] %>% sum()
       
        if (ma == 0 | nb == 0) {

          return(c(0,0))
        }
        else{
          this_pos = score_list[[a]][truth_list[[a]] == 1]
          this_neg = score_list[[b]][truth_list[[b]] == 0]
          temp_m1 = matrix(this_neg,nrow = ma,ncol = nb,byrow = TRUE)
          temp_m2 = matrix(this_pos,nrow = ma,ncol = nb,byrow = FALSE)
          temp_dif = temp_m2 - temp_m1
          temp_dif[temp_dif > 0] = 1
          temp_dif[temp_dif == 0] = 0.5
          temp_dif[temp_dif < 0] = 0
          return(c(sum(temp_dif),ma*nb))
        }

      }



      temp_out = apply(comb_pairs,MARGIN = 1,FUN = sum_comb)

      success_array[m,r,,,1] = temp_out[1,] %>% matrix(np,np,byrow = TRUE)
      success_array[m,r,,,2] = temp_out[2,] %>% matrix(np,np,byrow = TRUE)


    }

  }



  return(list(success_array,num_ROI))
}


#' Delete diagonal term function
#'
#' @param m the input matrix for deleting diagonal term
#'
#' @return diagonal term removed matrix

delete_diag = function(m){
  return(sum(m) - sum(diag(m)))
}

# sum the diagonal terms ----
#' sum the diagonal terms
#'
#' @param m input matrix
#'
#' @return sum of diagonal terms

sum_diag = function(m){sum(diag(m))}


# function for calculating the 8th moment ----
#' function for calculating the 8th moment
#'
#' @param m input matrix
#'
#' @return the 8th moment
m8_f = function(m){
  np = dim(m)[1]
  return(sum((m %>% matrix(1,np^2))*(m %>% t() %>% matrix(1,np^2))))
}
# function for calculating the 11th moment ----
#' function for calculating the 11th moment
#'
#' @param m input matrix
#'
#' @return the 11th moment

m11_f = function(m){
  return(sum(m^2))
}
#' covariance 8th moment middle calculation part one
#'
#' @param m input matrix
#'
#' @return the middle values for calculating covariance 8th moment

cov_m8_f1 = function(m) {
  np = dim(m)[1]
  return(m %>% t() %>% matrix(1,np^2))
}
#' covariance 8th moment middle calculation part two
#'
#' @param m input matrix
#'
#' @return the middle values for calculating covariance 8th moment
cov_m8_f2 = function(m) {
  np = dim(m)[1]
  return(m %>% matrix(1,np^2))
}


# Function for calculating 2 reader AUC covariance ----

#' Function for calculating 2 reader AUC covariance
#'
#' @param success_score The success score for nested data
#'
#' @return the covariance between two readers' AUC

AUC_cov_2reader_nest = function(success_score){
  ## Calculate all the perfect squares and counts for each perfect square

  success_array = success_score[[1]]

  
  nr = dim(success_array)[2]
  nm = dim(success_array)[1]
  ## 1st cov perfect square
  m1 = apply(success_array[,1,,,1],c(1),sum)*apply(success_array[,2,,,1],c(1),sum)
  c1 = apply(success_array[,1,,,2],c(1),sum)*apply(success_array[,2,,,2],c(1),sum)
  
  
  ## 2nd perfect square
  m2 = (apply(success_array[,1,,,1],c(1,2),sum)*apply(success_array[,2,,,1],c(1,2),sum)) %>% apply(c(1),sum)
  c2 = (apply(success_array[,1,,,2],c(1,2),sum)*apply(success_array[,2,,,2],c(1,2),sum)) %>% apply(c(1),sum)
  ## 3rd 2 perfect squares
  m31 = apply(success_array[,1,,,1],c(1),sum_diag)*apply(success_array[,2,,,1],c(1),sum)
  c31 = apply(success_array[,1,,,2],c(1),sum_diag)*apply(success_array[,2,,,2],c(1),sum)
  
  m32 = apply(success_array[,2,,,1],c(1),sum_diag)*apply(success_array[,1,,,1],c(1),sum)
  c32 = apply(success_array[,2,,,2],c(1),sum_diag)*apply(success_array[,1,,,2],c(1),sum)
  ## 4th 2 perfect square
  
  m41 = ((apply(success_array[,1,,,1],c(1),diag) %>% aperm(c(2,1)))*apply(success_array[,2,,,1],c(1,2),sum)) %>% apply(c(1),sum)
  c41 = ((apply(success_array[,1,,,2],c(1),diag) %>% aperm(c(2,1)))*apply(success_array[,2,,,2],c(1,2),sum)) %>% apply(c(1),sum)
  
  m42 = ((apply(success_array[,2,,,1],c(1),diag) %>% aperm(c(2,1)))*apply(success_array[,1,,,1],c(1,2),sum)) %>% apply(c(1),sum)
  c42 = ((apply(success_array[,2,,,2],c(1),diag) %>% aperm(c(2,1)))*apply(success_array[,1,,,2],c(1,2),sum)) %>% apply(c(1),sum)
  
  ## 5th 2 perfect square
  m51 = (apply(success_array[,1,,,1],c(1,3),sum)*apply(success_array[,2,,,1],c(1,2),sum)) %>% apply(c(1),sum)
  c51 = (apply(success_array[,1,,,2],c(1,3),sum)*apply(success_array[,2,,,2],c(1,2),sum)) %>% apply(c(1),sum)
  
  m52 = (apply(success_array[,2,,,1],c(1,3),sum)*apply(success_array[,1,,,1],c(1,2),sum)) %>% apply(c(1),sum)
  c52 = (apply(success_array[,2,,,2],c(1,3),sum)*apply(success_array[,1,,,2],c(1,2),sum)) %>% apply(c(1),sum)
  
  ## 6th 2 perfect square
  
  m61 = ((apply(success_array[,1,,,1],c(1),diag) %>% aperm(c(2,1)))*apply(success_array[,2,,,1],c(1,3),sum)) %>% apply(c(1),sum)
  c61 = ((apply(success_array[,1,,,2],c(1),diag) %>% aperm(c(2,1)))*apply(success_array[,2,,,2],c(1,3),sum)) %>% apply(c(1),sum)
  
  m62 = ((apply(success_array[,2,,,1],c(1),diag) %>% aperm(c(2,1)))*apply(success_array[,1,,,1],c(1,3),sum)) %>% apply(c(1),sum)
  c62 = ((apply(success_array[,2,,,2],c(1),diag) %>% aperm(c(2,1)))*apply(success_array[,1,,,2],c(1,3),sum)) %>% apply(c(1),sum)
  
  ## 7th perfect square
  m7 = (apply(success_array[,1,,,1],c(1),diag) %>% aperm(c(2,1))*apply(success_array[,2,,,1],c(1),diag) %>% aperm(c(2,1))) %>% apply(c(1),sum)
  c7 = (apply(success_array[,1,,,2],c(1),diag) %>% aperm(c(2,1))*apply(success_array[,2,,,2],c(1),diag) %>% aperm(c(2,1))) %>% apply(c(1),sum)
 
  ## 8th perfect square
  m8 = (apply(success_array[,1,,,1],c(1),cov_m8_f1)*apply(success_array[,2,,,1],c(1),cov_m8_f2)) %>% t() %>% apply(c(1),sum)
  c8 = (apply(success_array[,1,,,2],c(1),cov_m8_f1)*apply(success_array[,2,,,2],c(1),cov_m8_f2)) %>% t() %>% apply(c(1),sum)
  ## 9th perfect square
  
  m9 = apply(success_array[,1,,,1],c(1),sum_diag)*apply(success_array[,2,,,1],c(1),sum_diag)
  c9 = apply(success_array[,1,,,2],c(1),sum_diag)*apply(success_array[,2,,,2],c(1),sum_diag)
  
  ## 10th perfect square
  m10 = (apply(success_array[,1,,,1],c(1,3),sum) * apply(success_array[,2,,,1],c(1,3),sum)) %>% apply(c(1),sum)
  c10 = (apply(success_array[,1,,,2],c(1,3),sum) * apply(success_array[,2,,,2],c(1,3),sum)) %>% apply(c(1),sum)
  ## 11th perfect square
  m11 = (success_array[,1,,,1]*success_array[,2,,,1]) %>% apply(c(1),sum)
  c11 = (success_array[,1,,,2]*success_array[,2,,,2]) %>% apply(c(1),sum)
  
  
  ## Combine perfect squares
  first_square = (m1 - m2 - m31 - m32 + 2*m41 + 2*m42 - m51 - m52 + 2*m61 + 2*m62 - 6*m7 + m8 + m9 - m10 + m11)/(c1 - 
                       c2 - c31 - c32 + 2*c41 + 2*c42 - c51 - c52 + 2*c61 + 2*c62 - 6*c7 + c8 + c9 - c10 + c11)
  AUC_per_reader = apply(success_array[,,,,1],c(1,2),delete_diag)/apply(success_array[,,,,2],c(1,2),delete_diag)

  
  AUC_COV_per_reader = AUC_per_reader[,1]*AUC_per_reader[,2] - first_square
  
 
  
  return(AUC_COV_per_reader)
}
