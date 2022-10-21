#include <Rcpp.h>
using namespace Rcpp;
//' Calculate the each moments coefficient in variance 
//'
//' @param numROI number of ROIs in each case
//' @return all the coefficients
// [[Rcpp::export]]
NumericVector var_coef(NumericMatrix numROI){
 
  int n = numROI.ncol();
  NumericVector coef(11);
  
  for (int i=0; i<n; i++){
  	for (int j=0; j<n; j++){
        	if(j == i){
                  continue;
                 }else{
            
        coef[0] = coef[0] + numROI(0,i)*numROI(1,j);
        coef[1] = coef[1] + numROI(0,i)*numROI(1,j)*(numROI(0,i)-1)*(numROI(1,j)-1);
        coef[2] = coef[2] + numROI(0,i)*numROI(1,j)*(numROI(0,i)-1);
        coef[3] = coef[3] + numROI(0,i)*numROI(1,j)*(numROI(1,j)-1);
        coef[6] = coef[6] + numROI(0,i)*numROI(1,j)*numROI(0,j)*numROI(1,i);

        	for (int h=0; h<n; h++){
 			     if ((h==j) | (h==i)) {continue;
                         }else{

            coef[4] = coef[4]+numROI(0,i)*numROI(1,j)*numROI(1,h);
            coef[5] = coef[5]+numROI(0,i)*(numROI(0,i)-1)*numROI(1,j)*numROI(1,h);
            coef[7] = coef[7] + numROI(0,i)*numROI(0,j)*numROI(1,h)*numROI(1,i);
            coef[8] = coef[8] + numROI(0,i)*numROI(0,j)*numROI(1,j)*numROI(1,h);
            coef[9] = coef[9] + numROI(0,i)*numROI(0,j)*numROI(1,h);
            coef[10] = coef[10] + numROI(0,i)*numROI(0,j)*numROI(1,h)*(numROI(1,h)-1);
            }
            }
      } 
    }
  }
  return coef;

}
