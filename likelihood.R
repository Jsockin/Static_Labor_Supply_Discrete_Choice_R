#
# This function will calculate the likelihood for an individual,
# which can then be summed across individuals to obtain the total
# likelihood function for the static discrete choice model.
#

likelihood <- function(h_opt,h_grid,y,w,gamma,theta,alpha,tau,tax_thresh,mu_w,sigma_w){
  
    # Calculate individual's likelihood
    pos_hstar <- which(h_grid==h_opt)
    
    if (h_opt == 0){
      
      # Monte Carlo appraoch to prob when h=0
      set.seed(10)
      draws     <- 1000
      logw_vals <- rlnorm(draws,mu_w,sigma_w)
      
      # Calculate utility for each wage and each h in h_grid
      wages     <- exp(logw_vals)
      is_taxed  <- lapply(seq_len(length(wages)),function(ii) wages[ii]*h_grid > tax_thresh )
      taxes     <- lapply(seq_len(length(wages)),function(ii) is_taxed[[ii]] * (tau*wages[ii]*h_grid) )
      c         <- lapply(seq_len(length(wages)),function(ii) y + wages[ii]*h_grid - taxes[[ii]] )
      utility   <- lapply(seq_len(length(wages)),function(ii) gamma * ( (c[[ii]]^theta)/theta - alpha*h_grid) )
      prob_h0   <- lapply(seq_len(length(wages)),function(ii) exp(utility[[ii]][1]) / sum(exp(utility[[ii]])) )
      prob      <- mean(unlist(prob_h0))
      
    } else{
      
      # Calculate utility for each h in h_grid
      is_taxed  <- w*h_grid > tax_thresh 
      taxes     <- is_taxed * (tau*w*h_grid)
      c         <- y + w*h_grid - taxes
      utility   <- gamma * ( (c^theta)/theta - alpha*h_grid)
      
      # Likelihood for P_i = 1 (hours = h, wage = w)
      prob_h    <- exp(utility[pos_hstar]) / sum(exp(utility))
      prob_w    <- dlnorm(w,mu_w,sigma_w)
      prob      <- prob_h * prob_w
      
    }
      
    return(prob)
}