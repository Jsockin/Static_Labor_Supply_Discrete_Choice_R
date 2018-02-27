#
#   This function will find h that maximizes utility u(c,h)
#

calcHours <- function(epsilon,h_grid,y,w,gamma,theta,alpha,tau,tax_thresh,tax_policy,tax_setting){
  
  # Calculate taxes paid
  if (tax_policy == 1){
    is_taxed  <- w*h_grid > tax_thresh 
    taxes     <- is_taxed * (tau*w*h_grid)
  } else{
    taxes     <- tau*w*h_grid
  }
  
  # Calculate utility  
  c       <- y + w*h_grid - taxes
  utility <- gamma * ( (c^theta)/theta - alpha*h_grid) + epsilon
  
  # Find optimal h
  hstar   <- h_grid[which(utility==max(utility))]
  return(hstar)
  
}  