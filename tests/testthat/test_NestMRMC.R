library(NestMRMC)
library(testthat)
context("NestMRMC_R")








# Simulate the Demo data set or Load the Demo data  
# It serves as the baseline to be compared with

generate_new_data = FALSE

if (generate_new_data) {
  
# simulate the parameters with default values
sim.config = simu_config()
    
expected_data = data_MRMC(sim.config)$data_final

} else {
  expected_data = NestMRMC::expected_data
}




# Simulate the same data set for testing 

sim.config = simu_config()
current_data = data_MRMC(sim.config)$data_final

 

# Do the testing for simulation data
test_that(
  "NestMRMC Demo simulation data does not change", {
    expect_equal(expected_data, current_data)
  }
)


## Test the analysis results

system.time(
  {
    expected_output = AUC_per_reader_nest(expected_data) 
  }
)


system.time(
  {
    current_output = AUC_per_reader_nest(current_data)
  }
)





# Do the testing for analysis results
test_that(
  "NestMRMC analysis result does not change", {
    expect_equal(expected_output, current_output)
  }
)



# Do the testing for theoretical (co)variance results
# Comment the other arguments 
system.time(
  {expected_var = true_AUC_var_abitrary(expected_output$numROI)}
)




system.time(
  {
    current_var = true_AUC_var_abitrary(current_output$numROI)
  }
)





test_that(
  "NestMRMC theoretical (co)variance result does not change", {
    expect_equal(expected_output, current_output)
  }
)


## Test the Rcpp results with org version

system.time(
  {
   
    expected_Rcpp = true_AUC_var_abitrary_Rcpp(numROI = expected_output$numROI) 
    
  }
)

test_that(
  "NestMRMC theoretical (co)variance result does not change", {
    expect_equal(expected_var, expected_Rcpp)
  }
)


 
