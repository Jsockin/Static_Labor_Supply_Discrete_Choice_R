#
#   This program will solve problem set 1 for Economics of Labor.
#
#   Author: Jason Sockin
#   Date:   January 2018

#------------------------------
#   0. Housekeeping and Functions 
#------------------------------

rm(list=ls())
set.seed(1)
path <- file.path("C:","Users","jsock","Dropbox","Penn","Homework","Shephard","HW_1")

# Load libraries
#install.packages("gumbel")
library(gumbel)
library(stats)

# Import functions
source(paste(path,'/calcHours.R',sep=''))
source(paste(path,'/calcTaxes.R',sep=''))
source(paste(path,'/likelihood.R',sep=''))

#------------------------------
#   1. Environment 
#------------------------------

agents      <- 1000
h_grid      <- c(0,10,20,30,40)

#     U(c,h) = gamma * [c^theta / theta - alpha*h] + epsilon_h
#       c = y + wh - T(wh)
#       log(w) = mu_w + epsilon_w, epsilon_w ~ N(0,sigma_w)

theta       <- 0.3
alpha       <- 0.1
gamma       <- 2
tau         <- 0.3
tax_thresh  <- 80
mu_w        <- 1
sigma_w     <- 0.55

# Stochastic draws
y_vals    <- runif(agents,10,100)
log_w     <- mu_w + rnorm(agents,0,sigma_w)
wage      <- exp(log_w)
epsilon   <- matrix(rgumbel(agents*length(h_grid),1,1),agents,length(h_grid))

#------------------------------
#   2. Agents' Labor Supply Decisions 
#------------------------------

tax_policy <- 1
h_star      <- unlist(lapply(seq_len(agents),function(x) calcHours(epsilon[x,],h_grid,y_vals[x],wage[x],gamma,theta,alpha,tau,tax_thresh,tax_policy)))

# Plot distribution of hours worked
savefile <- paste(path,'/Distribution_h.eps',sep='')
postscript(savefile, horizontal = FALSE, onefile = FALSE)
hist(h_star,breaks=15,main='Distribution of Hours Worked (h)')
dev.off()

# Non-labor income distribution given work hours
savefile <- paste(path,'/Nonlabor_income_hours.eps',sep='')
postscript(savefile, horizontal = FALSE, onefile = FALSE)
df_y_h    <- cbind(h_star,y_vals)
dens_y_h  <- lapply(seq_len(length(h_grid)),function(x) density(df_y_h[df_y_h[,1]==h_grid[x],2]))
plot(dens_y_h[[1]],ylim=c(0,0.08),xlim=c(0,120),col=1,xlab='non-labor income (y)',ylab='density',main='Distribution of Non-labor income (y) given hours worked (h)')
for (ii in seq(2,length(h_grid))){
  par(new=TRUE)
  plot(dens_y_h[[ii]],ylim=c(0,0.08),xlim=c(0,120),col=ii,xlab='non-labor income',ylab='density',main='')
}
legend(80,0.06,paste('h =',h_grid),lty=rep(1,length(h_grid)),col=1:length(h_grid),bty='n',cex=1.25)
dev.off()

# Wage distribution given work hours
savefile <- paste(path,'/Wage_hours.eps',sep='')
postscript(savefile, horizontal = FALSE, onefile = FALSE)
df_w_h    <- cbind(h_star,wage)
dens_w_h  <- lapply(seq_len(length(h_grid)),function(x) density(df_w_h[df_w_h[,1]==h_grid[x],2]))
plot(dens_w_h[[1]],ylim=c(0,1.1),xlim=c(0,12),col=1,xlab='wage offered (w)',ylab='density',main='Distribution of wage offered (w) given hours worked (h)')
for (ii in seq(2,length(h_grid))){
  par(new=TRUE)
  plot(dens_w_h[[ii]],ylim=c(0,1.1),xlim=c(0,12),col=ii,xlab='wage offered (w)',ylab='density',main='')
}
legend(9,1,paste('h =',h_grid),lty=rep(1,length(h_grid)),col=1:length(h_grid),bty='n',cex=1.25)
dev.off()

#------------------------------
#   3. Agents' Labor Supply Decisions Per Gamma
#------------------------------

gamma_vals  <- seq(0,3,0.02)
h_star_gamma  <- lapply(seq_len(length(gamma_vals)), function(ii) unlist(lapply(seq_len(agents),function(x) calcHours(epsilon[x,],h_grid,y_vals[x],wage[x],gamma_vals[ii],theta,alpha,tau,tax_thresh,tax_policy))))
reg_df      <- lapply(seq_len(length(gamma_vals)), function(ii) as.data.frame(cbind(wage,h_star_gamma[[ii]])))

# regression for all h
reg_result  <- lapply(seq_len(length(gamma_vals)), function(ii) lm(V2 ~ wage,data=reg_df[[ii]]))
wage_coef   <-unlist(lapply(seq_len(length(gamma_vals)), function(ii) reg_result[[ii]]$coefficients[2]))    

# regression conditional h > 0
reg_result_working  <- lapply(seq_len(length(gamma_vals)), function(ii) lm(V2 ~ wage,data=reg_df[[ii]][reg_df[[ii]]$V2>0,]))
wage_coef_working   <-unlist(lapply(seq_len(length(gamma_vals)), function(ii) reg_result_working[[ii]]$coefficients[2]))    

savefile <- paste(path,'/Bias_gamma.eps',sep='')
postscript(savefile, horizontal = FALSE, onefile = FALSE)
plot(gamma_vals,wage_coef,type='l',col="blue",ylim=c(0,4),ylab="",xlab="",main="",cex.main=1.5)
par(new=TRUE)
plot(gamma_vals,wage_coef_working,type='l',col="red",ylim=c(0,4),ylab="beta",xlab="gamma",main="hours = beta * wage + error",cex.main=1.5)
legend(2,1,c('All h','Given h > 0'),lty=c(1,1),col=c('blue','red'),bty='n',cex=1.25)
dev.off()

#------------------------------
#   4. Calculate likelihood 
#------------------------------

log_like  <- function(vec) -sum(log(unlist(lapply(seq_len(agents),function(x) likelihood(h_star[x],h_grid,y_vals[x],wage[x],
                                  vec[1],vec[2],vec[3],tau,tax_thresh,vec[4],vec[5])))))
theta_0   <- c(gamma,theta,alpha,mu_w,sigma_w)
theta_hat <- optim(theta_0,log_like)

#------------------------------
#   5. Counterfactual Scenario (constant marginal tax rate = 0.2)
#------------------------------

tax_policy        <- 2
tau_prime         <- 0.2
tax_thresh_prime  <- 0
h_star_prime      <- unlist(lapply(seq_len(agents),function(x) calcHours(epsilon[x,],h_grid,y_vals[x],wage[x],gamma,theta,alpha,tau_prime,tax_thresh_prime,tax_policy)))

# Plot new distribution of hours worked
savefile <- paste(path,'/Distribution_h_counterfactual.eps',sep='')
postscript(savefile, horizontal = FALSE, onefile = FALSE)
hist(h_star_prime,breaks=15,main='Distribution of Hours Worked Under Constant Marginal Tax Rate')
dev.off()

# Calculate tax revenue under each policy
tax_revenue_original <-sum(unlist(lapply(seq_len(agents),function(x) calcTaxes(h_star[x],y_vals[x],wage[x],gamma,theta,alpha,tau,tax_thresh))))
tax_revenue_counterfactual <-sum(unlist(lapply(seq_len(agents),function(x) calcTaxes(h_star_prime[x],y_vals[x],wage[x],gamma,theta,alpha,tau_prime,tax_thresh_prime))))

# Calculate cons. marg. tax rate to keep revenue unchanged
residualRevenue <- function(tax_rate) tax_revenue_original - sum(unlist(lapply(seq_len(agents),function(x) calcTaxes(h_star_prime[x],y_vals[x],wage[x],gamma,theta,alpha,tax_rate,tax_thresh_prime))))
opt_tau         <- uniroot(residualRevenue,c(0,1))$root

