#########################
# Necessary Functions
#########################
### Estimation
uDCS_t_model_estimator <- function(dati, param){
  
  Start <- Sys.time()
  ###Take T
  T <- length(dati)
  
  ###Parameter Selections Dynamic Location
  omega <- param[1]
  phi   <- param[2]
  k     <- param[3]
  
  varsigma <- param[4]
  nu       <- param[5]
  
  ###Create a vector with the parameters
  theta_st <- c(omega, phi, k, varsigma, nu)
  
  ###Take Bounds
  lower <- c(-Inf, -0.999, -2, 1e-05, 2.099)
  upper <- c( Inf,  0.999,  2, Inf, 300)
  
  #------> Optimize uDCS_t_model Filters w/L-BFGS-B
  #optimizer <- suppressWarnings(optim(par = theta_st, fn = interprete_uDCS_t_model, 
  #                                    dati = dati, method = "L-BFGS-B", 
  #                                    control = list(trace = 1), hessian = FALSE,
  #                                    lower = lower, upper = upper))
  
  #------> Optimize every uDCS_t_model Filters w/solnp := WAY faster than optim
  #optimizer <- suppressWarnings(Rsolnp::solnp(pars = theta_st, fun = interprete_uDCS_t_model, 
  #                                            dati  = dati, control = list(trace = 1), 
  #                                            LB = lower, UB = upper)) 
  
  #------> Optimize every uDCS_t_model Filters w/nlminb := WAY faster than optim
  optimizer <- suppressWarnings(nlminb(start = theta_st, objective = interprete_uDCS_t_model, 
                                       dati  = dati, gradient = NULL, 
                                       control = list(trace = 0), hessian = NULL,
                                       lower = lower, upper = upper))
  
  #------> Save the optimized parameters Dynamic Location
  omega_opt <- optimizer$par[1]  
  phi_opt <- optimizer$par[2]
  k_opt  <- (optimizer$par[3])
  
  varsigma_opt <- optimizer$par[4]
  nu_opt       <- optimizer$par[5]
  
  ###Create a vector with ALL the optimized parameters
  theta_opt <- c(omega_opt, phi_opt, k_opt, varsigma_opt, nu_opt)
  
  ###Create a list with ALL the optimized parameters
  theta_list <- list(omega = omega_opt,
                     phi = phi_opt,
                     k  = k_opt,
                     varsigma = varsigma_opt,
                     nu    = nu_opt)
  
  ######################
  ####### OUTPUT #######
  ######################
  
  #------> Some detail
  Elapsed_Time <- Sys.time() - Start
  print(paste("Elapsed Time: ", toString(Elapsed_Time)))
  
  ###Make List
  out <- list(theta_list = theta_list,
              theta      = theta_opt,
              optimizer  = optimizer)
  
  return(out) 
}

################
##### Interprete
################
interprete_uDCS_t_model <- function(dati, param){
  
  ###Take T
  T <- length(dati)
  
  ###Parameter Selections Dynamic Location
  omega <- param[1]
  phi   <- param[2]
  k     <- param[3]
  
  varsigma <- param[4]
  nu       <- param[5]
  
  ###Create a new vector with the parameters
  theta_new <- c(omega, phi, k, varsigma, nu)
  
  #------> Fitness Functions
  fitness <- uDCS_t_model_filter(dati, theta_new)$Log_Likelihood
  
  if(is.na(fitness) | !is.finite(fitness)) fitness <- -1e10
  if(fitness != fitness) fitness <- -1e10
  
  return(-fitness)
} 


### Filtration
uDCS_t_model_filter <- function(y, theta){
  
  ###Take T
  T <- length(y)
  
  ###Define LogLikelihoods
  dloglik <- array(data = NA, dim = c(T))
  loglik  <- numeric()
  
  ###Parameter Selections Dynamic Location
  omega <- theta[1]
  phi   <- theta[2]
  k     <- theta[3]
  
  varsigma <- theta[4]
  nu       <- theta[5]
  
  ###Define Dynamic Location and Innovations
  mu_t <- array(data = NA, dim = c(T+1))
  u_t  <- array(data = NA, dim = c(T))
  
  ###Initialize Dynamic Location
  mu_t[1]   <- (omega)
  
  ###Initialize Likelihood
  dloglik[1] <- uSTDT_uDCS_t(y[1], mu_t[1], varsigma = varsigma, nu = nu, log = TRUE)
  loglik     <- dloglik[1]
  
  for(t in 2:(T+1)) {
    ###Dynamic Location Innovations
    u_t[t-1] <- martingale_diff_u_t(y[t-1], mu_t[t-1], varsigma, nu)
    ###Updating Filter                    
    mu_t[t]   <- omega + phi * (mu_t[t-1] - omega) + k * u_t[t-1]
    
    if(t < (T+1)){
      ###Updating Likelihoods
      dloglik[t] <- uSTDT_uDCS_t(y[t], mu_t = mu_t[t], varsigma = varsigma, nu = nu, log = TRUE)
      loglik     <- loglik + dloglik[t]
    }
  }
  
  ######################
  ####### OUTPUT #######
  ######################
  mu_t <- ts(mu_t, start = start(y), frequency = frequency(y))
  u_t  <- ts(u_t, start = start(y), frequency = frequency(y))
  
  ###Make List
  out <- list(Dynamic_Location = mu_t,
              Innovation_u_t   = u_t,
              Log_Densities_i  = dloglik,
              Log_Likelihood   = loglik)
  
  return(out)
}

############################################
####### ADDITIONAL FUNCTIONS #######
############################################

########################################################
martingale_diff_u_t <- function(y, mu_t, varsigma, nu){
  
  u_t <- c((1 / (1 + (y - mu_t)^2/(nu*varsigma)) * (y - mu_t)))
  
  return(u_t)
}

########################################################
uSTDT_uDCS_t <- function(y, mu_t, varsigma, nu, log = TRUE){
  
  ulpdf <- (lgamma((nu + 1) / 2) - lgamma(nu / 2) - (1/2) * log(varsigma) -
              (1/2)  * log(pi * nu) - ((nu + 1) / 2) * log(1 + (y - mu_t)^2 / (nu*varsigma) ))
  
  if(log != TRUE){
    ulpdf <- exp(ulpdf)
  } 
  
  return(ulpdf)
}


### Simulation
uDCS_t_model_simulator <- function(T, omega, phi, k, varsigma, nu){
  
  ###Define the Processes
  y   <- array(data = NA, dim = c(T) ) 
  
  ###Define Dynamic Location and Innovations
  mu_t <- array(data = NA, dim = c(T))
  u_t  <- array(data = NA, dim = c(T-1))
  
  ###Initial value for the recursion 
  mu_t[1]   <- omega
  
  ###Generate the first observations of the process
  y[1]   <- uSTDT_rnd(1, mu_t[1], varsigma, nu)
  
  ###Dynamics 
  for (t in 2:T) {
    
    ###Factor Innovations
    u_t[t-1] <- martingale_diff_u_t(y[t-1], mu_t[t-1], varsigma, nu)
    
    ###Updating Filters                    
    mu_t[t]   <- omega + phi * (mu_t[t-1] - omega) + k * u_t[t-1]
    
    ###Generate the observations of the processes
    y[t] <- uSTDT_rnd(1, mu_t[t], varsigma, nu)
  }
  ######################
  ####### OUTPUT #######
  ######################
  
  ###Make List
  out <- list(y_t_gen          = as.ts(y),
              Dynamic_Location = as.ts(mu_t),
              Innovation_u_t   = as.ts(u_t))
  
  return(out)
}

############################################################################
################# Univariate Student's t Random Generator ##################
############################################################################

uSTDT_rnd <- function(n, mu, varsigma, nu) {
  
  z <- rt(n, df = nu) 
  y <- numeric()
  for(i in 1:n){
    y[i] <- c(mu + z[i]* sqrt(varsigma) ) 
  }
  
  return(y)
}



#########################
# Question No.2
#########################
# The case with nu = 3
T <- 2000
omega <- 0
phi <- 0.8
k <- 0.8
varsigma <- 1
nu <- 3

theta <- c(omega, phi, k, varsigma, nu)

simu <- uDCS_t_model_simulator(T, omega, phi, k, varsigma, nu)
y = simu$y_t_gen
ts.plot(y)


filter <- uDCS_t_model_filter(y, theta)

v = 0 
for(t in 2:(T+1)){
  v[t-1]=y[t-1]-filter$Dynamic_Location[t-1]
}


ts.plot(y)
lines(v,col = "blue")
lines(filter$Innovation_u_t,col = "red")

# We can see that ut has no sensitivity to extreme observations while vt follows the extreme observations very closely
# and vt is very sensitive to extreme observations 


# The case with nu = 10
T <- 2000
omega <- 0
phi <- 0.8
k <- 0.8
varsigma <- 1
nu <- 10

theta <- c(omega, phi, k, varsigma, nu)

simu <- uDCS_t_model_simulator(T, omega, phi, k, varsigma, nu)
y = simu$y_t_gen
ts.plot(y)


filter <- uDCS_t_model_filter(y, theta)

v = 0 
for(t in 2:(T+1)){
  v[t-1]=y[t-1]-filter$Dynamic_Location[t-1]
}


ts.plot(y)
lines(v,col = "blue")
lines(filter$Innovation_u_t,col = "red")

# same as before, we can see that ut has no sensitivity to extreme observations while vt follows the extreme
# observations very closely and vt is very sensitive to extreme observations
# compared to the previous case, ut is more sensitive to extreme observations
# In this case we have also more extreme observations


# The case with nu = 200
T <- 2000
omega <- 0
phi <- 0.8
k <- 0.8
varsigma <- 1
nu <- 200

theta <- c(omega, phi, k, varsigma, nu)

simu <- uDCS_t_model_simulator(T, omega, phi, k, varsigma, nu)
y = simu$y_t_gen
ts.plot(y)


filter <- uDCS_t_model_filter(y, theta)

v = 0 
for(t in 2:(T+1)){
  v[t-1]=y[t-1]-filter$Dynamic_Location[t-1]
}


ts.plot(y)
lines(v,col = "blue")
lines(filter$Innovation_u_t,col = "red")

# Here we can see that the values of u are much closer to the values of v (They are almost equal) than the previous
# cases, and it can be said that vt is less sensitive to extreme observations but ut is more sensitive to extreme
# observations than the previous cases
# In general we can say as nu increases, ut and vt become essentially the same (it can be implied from the formula of
# ut where we can see that as nu goes to infinity, we will have ut = vt), and hence ut will be more sensitive to 
# extreme observations


##############################
# Question 3
##############################

T <- 2000
omega <- 0
phi <- 0.9
k <- 1
varsigma <- 1.3
nu <- 5

theta <- c(omega, phi, k, varsigma, nu)

simu <- uDCS_t_model_simulator(T, omega, phi, k, varsigma, nu)
y = simu$y_t_gen

### part a
ts.plot(y)

### part b
filter <- uDCS_t_model_filter(y, theta)
ts.plot(filter$Dynamic_Location)

ts.plot(y)
lines(filter$Dynamic_Location, col = 'green')

# we can see that the time varying parameter (location) is not following most of extreme observations and it is
# far away from very extreme observations

### part c

v = 0 
for(t in 2:(T+1)){
  v[t-1]=y[t-1]-filter$Dynamic_Location[t-1]
}


ts.plot(y)
lines(v,col = "blue")
lines(filter$Innovation_u_t,col = "red")

# we can see that vt is very sensitive to extreme observations but ut shows no sign of sensitivity to
# extreme observations

### part d

acf(y, lag.max = 100, drop.lag.0 = FALSE)
pacf(y, lag.max = 100, drop.lag.0 = FALSE)
# based on acf and pacf, yt is an AR(1) process

acf(filter$Innovation_u_t, lag.max = 100, drop.lag.0 = FALSE)
pacf(filter$Innovation_u_t, lag.max = 100, drop.lag.0 = FALSE)
acf(v, lag.max = 100, drop.lag.0 = FALSE)
pacf(v, lag.max = 100, drop.lag.0 = FALSE)
# based on acf and pacf, ut and vt are white noise processes


### part e
est  <- uDCS_t_model_estimator(y, theta)
est$theta_list

# omega = -0.149823, phi = 0.8798777, k = 1.021846, varsigma = 1.298681, and nu = 4.275446
# estimation of all of the parameters is close to their real values, except for nu, where the real value is 5 but the
# estimated value is 4.275446