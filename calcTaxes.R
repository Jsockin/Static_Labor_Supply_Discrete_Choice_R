#
#   This function will find total taxes paid.
#

calcTaxes <- function(h_opt,y,w,gamma,theta,alpha,tau,tax_thresh){
  
  is_taxed  <- w*h_opt > tax_thresh 
  taxes     <- is_taxed * (tau*w*h_opt)
  
  return(taxes)
}  