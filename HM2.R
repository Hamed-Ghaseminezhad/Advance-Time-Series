################################################################################
# Question No. 1
################################################################################
# AR(1) signal plus noise model with Gaussian and Student-t measurement noise
# the case with Sigma_e = sigma_eta = 1
################################################################################

# simulate the AR(1) time varying mean or signal (or trend)

# fix the length of the series 
n = 300

# simulate the noise of the signal (in the transition equation (SSF))
sigma_eta = 1
eta = sigma_eta * rnorm(n)

# simulate the measurement noise (in the observation equation)
sigma_e = 1

# q, signal-to-noise ratio, q = sigma^2_eta/sigma^2_e

q = sigma_eta^2/sigma_e^2

e = sigma_e * rnorm(n)

# plot of eta (error affecting the signal) versus epsilon (e, measurement error)

ts.plot(e)
lines(eta, col = "red")

phi <- 0.8
mu = 0
mu[1] <- 0  # unconditional mean of an AR(1) process with no intercept
y = 0 
for(t in 1:(n-1)){
  mu[t+1] = phi * mu[t] + eta[t]
  y[t+1] =  mu[t+1] + e[t+1]
}
ts.plot(mu)

# y = mu + e

ts.plot(y)
lines(mu,col="red")

# briefly check the correlation structure of the simulated process 
# to be sure that we are properly simulating an AR(1) signal plus noise process

acf(y, lag.max = 60, drop.lag.0 = FALSE)
pacf(y, lag.max = 60)
acf(mu, lag.max = 60, drop.lag.0 = FALSE)
pacf(mu, lag.max = 60)

# indeed, the acf and pacf of mu are those of an AR(1) process, while the acf
# and pacf of y  are those of an ARMA(1,1) process, as expected 


################################################################################
# KALMAN filter recursions (see Durbin and Koopman, 2001) 
################################################################################

# we have simulated y_t and an underlying signal mu_t  
# (observed and latent or unobserved in real data, respectively).
# Now we wish to estimate mu_{t|t-1} - note that we are 
# assuming that phi, sigma_e, sigma_eta, i.e. the static parameters, are known.
# In practice, they are estimated by maximising the Gaussian likelihood 
# (we shall see it) 

# initial conditions needed 
# mu_10, P_10 that we can set at the unconditional mean and variance of mu 


# KF recursion

# allocate space 
mu_pred <- 0  # this is mu_{t|t-1}
P <- 0        # this is P_{t|t-1} 
v <-0         # this  will be the innovation error 
K <- 0        # the Kalman gain 
F <- 0        # the conditional variance of v_t 
llk<-0        # the log-likelihood value 

# initialise the recursion 

# mu_pred[1] is set equal to its unconditional mean 
# P is the conditional variance of mu_t, let us initialise it as the signal-to-noise ratio

mu_pred[1] = 0
P[1] = (sigma_eta^2)/(1-phi^2)
llk[1]=0

# the recursion 

for(t in 1:(n-1)){
  v[t] = y[t] - mu_pred[t]; 
  F[t] = P[t] + sigma_e^2;
  K[t] = (phi * P[t])/F[t]
  P[t+1] = phi^2 * P[t] + sigma_eta^2 - K[t]*F[t]*K[t]
  mu_pred[t+1] = phi * mu_pred[t] + K[t]*v[t]
  llk[t] = - 0.5 * log(F[t] + (v[t]^2/F[t]))
}

llk = sum(llk)

ts.plot(y)
lines(mu, col = 'red')
lines(mu_pred, col = 'green')

ts.plot(y)
lines(mu, col = "red")

ts.plot(y)
lines(mu_pred, col = "green")

ts.plot(mu)
lines(mu_pred, col = "blue")

# in this case we have llk = -179.81

################################################################################
# Question No. 1
################################################################################
# AR(1) signal plus noise model with Gaussian and Student-t measurement noise
# the case with Sigma_e = 0.1 and sigma_eta = 1
################################################################################

# simulate the AR(1) time varying mean or signal (or trend)

# fix the length of the series 
n = 300

# simulate the noise of the signal (in the transition equation (SSF))
sigma_eta = 1
eta = sigma_eta * rnorm(n)

# simulate the measurement noise (in the observation equation)
sigma_e = .1

# q, signal-to-noise ratio, q = sigma^2_eta/sigma^2_e

q = sigma_eta^2/sigma_e^2

e = sigma_e * rnorm(n)

# plot of eta (error affecting the signal) versus epsilon (e, measurement error)

ts.plot(e)
lines(eta, col = "red")

phi <- 0.8
mu = 0
mu[1] <- 0  # unconditional mean of an AR(1) process with no intercept
y = 0 
for(t in 1:(n-1)){
  mu[t+1] = phi * mu[t] + eta[t]
  y[t+1] =  mu[t+1] + e[t+1]
}
ts.plot(mu)

# y = mu + e

# Since sigma_e is small, y and mu are the same (approximately)
ts.plot(y)
lines(mu,col="red")

# briefly check the correlation structure of the simulated process 
# to be sure that we are properly simulating an AR(1) signal plus noise process

acf(y, lag.max = 60, drop.lag.0 = FALSE)
pacf(y, lag.max = 60)
acf(mu, lag.max = 60, drop.lag.0 = FALSE)
pacf(mu, lag.max = 60)

# indeed, the acf and pacf of mu and y are those of an AR(1) process, as expected 


################################################################################
# KALMAN filter recursions (see Durbin and Koopman, 2001) 
################################################################################

# we have simulated y_t and an underlying signal mu_t  
# (observed and latent or unobserved in real data, respectively).
# Now we wish to estimate mu_{t|t-1} - note that we are 
# assuming that phi, sigma_e, sigma_eta, i.e. the static parameters, are known.
# In practice, they are estimated by maximising the Gaussian likelihood 
# (we shall see it) 

# initial conditions needed 
# mu_10, P_10 that we can set at the unconditional mean and variance of mu 


# KF recursion

# allocate space 
mu_pred <- 0  # this is mu_{t|t-1}
P <- 0        # this is P_{t|t-1} 
v <-0         # this  will be the innovation error 
K <- 0        # the Kalman gain 
F <- 0        # the conditional variance of v_t 
llk<-0        # the log-likelihood value 

# initialise the recursion 

# mu_pred[1] is set equal to its unconditional mean 
# P is the conditional variance of mu_t, let us initialise it as the signal-to-noise ratio

mu_pred[1] = 0
P[1] = (sigma_eta^2)/(1-phi^2)
llk[1]=0

# the recursion 

for(t in 1:(n-1)){
  v[t] = y[t] - mu_pred[t]; 
  F[t] = P[t] + sigma_e^2;
  K[t] = (phi * P[t])/F[t]
  P[t+1] = phi^2 * P[t] + sigma_eta^2 - K[t]*F[t]*K[t]
  mu_pred[t+1] = phi * mu_pred[t] + K[t]*v[t]
  llk[t] = - 0.5 * log(F[t] + (v[t]^2/F[t]))
}

llk = sum(llk)

ts.plot(y)
lines(mu, col = 'red')
lines(mu_pred, col = 'green')

ts.plot(y)
lines(mu, col = "red")

ts.plot(y)
lines(mu_pred, col = "green")

ts.plot(mu)
lines(mu_pred, col = "blue")

# in this case, we have llk = -84.11 which shows great improvement over the estimation of
# mu_{t | t-1} in the previous model (the one with sigma_e = sigma_eta = 1). We can say
# that the lower value of sigma_e have led to an better estimation



################################################################################
# Question No. 1
################################################################################
# AR(1) signal plus noise model with Gaussian and Student-t measurement noise
# the case with Sigma_e = 1 and sigma_eta = 0.1
################################################################################

# simulate the AR(1) time varying mean or signal (or trend)

# fix the length of the series 
n = 300

# simulate the noise of the signal (in the transition equation (SSF))
sigma_eta = .1
eta = sigma_eta * rnorm(n)

# simulate the measurement noise (in the observation equation)
sigma_e = 1

# q, signal-to-noise ratio, q = sigma^2_eta/sigma^2_e

q = sigma_eta^2/sigma_e^2

e = sigma_e * rnorm(n)


# plot of eta (error affecting the signal) versus epsilon (e, measurement error)

ts.plot(e)
lines(eta, col = "red")

phi <- 0.8
mu = 0
mu[1] <- 0  # unconditional mean of an AR(1) process with no intercept
y = 0 
for(t in 1:(n-1)){
  mu[t+1] = phi * mu[t] + eta[t]
  y[t+1] =  mu[t+1] + e[t+1]
}
ts.plot(mu)

# y = mu + e

# Since sigma_eta is small and sigma_e and sigma_eta are distant, y has higher values
ts.plot(y)
lines(mu,col="red")

# briefly check the correlation structure of the simulated process 
# to be sure that we are properly simulating an AR(1) signal plus noise process

acf(y, lag.max = 60, drop.lag.0 = FALSE)
pacf(y, lag.max = 60)
acf(mu, lag.max = 60, drop.lag.0 = FALSE)
pacf(mu, lag.max = 60)

# here, the acf and pacf of y suggest that there is no significant (considerable) correlation between the values of y
# as all of the lags of acf (except for zero, naturally) are not significant, or we can say y exhibits properties of 
# a white noise process. On the other hand, the acf and pacf of mu display AR(1) model behaviour


################################################################################
# KALMAN filter recursions (see Durbin and Koopman, 2001) 
################################################################################

# we have simulated y_t and an underlying signal mu_t  
# (observed and latent or unobserved in real data, respectively).
# Now we wish to estimate mu_{t|t-1} - note that we are 
# assuming that phi, sigma_e, sigma_eta, i.e. the static parameters, are known.
# In practice, they are estimated by maximising the Gaussian likelihood 
# (we shall see it) 

# initial conditions needed 
# mu_10, P_10 that we can set at the unconditional mean and variance of mu 

# KF recursion

# allocate space 
mu_pred <- 0  # this is mu_{t|t-1}
P <- 0        # this is P_{t|t-1} 
v <-0         # this  will be the innovation error 
K <- 0        # the Kalman gain 
F <- 0        # the conditional variance of v_t 
llk<-0        # the log-likelihood value 

# initialise the recursion 

# mu_pred[1] is set equal to its unconditional mean 
# P is the conditional variance of mu_t, let us initialise it as the signal-to-noise ratio

mu_pred[1] = 0
P[1] = (sigma_eta^2)/(1-phi^2)
llk[1]=0

# the recursion 

for(t in 1:(n-1)){
  v[t] = y[t] - mu_pred[t]; 
  F[t] = P[t] + sigma_e^2;
  K[t] = (phi * P[t])/F[t]
  P[t+1] = phi^2 * P[t] + sigma_eta^2 - K[t]*F[t]*K[t]
  mu_pred[t+1] = phi * mu_pred[t] + K[t]*v[t]
  llk[t] = - 0.5 * log(F[t] + (v[t]^2/F[t]))
}

llk = sum(llk)

ts.plot(y)
lines(mu, col = 'red')
lines(mu_pred, col = 'green')

ts.plot(y)
lines(mu, col = "red")

ts.plot(y)
lines(mu_pred, col = "green")

ts.plot(mu)
lines(mu_pred, col = "blue")

# In this case we have llk = -78.67. We have seen that the best estimate of mu_{t|t-1} among the three cases,
# is when we have sigma_e = 1 and sigma_eta = 0.1 (which makes sense as mu[t+1] = phi * mu[t] + eta[t]), 
# and the worst estimate is for the case with sigma_e = sigma_eta = 1
# llk for all the cases respectively: (-179.81, -84.11, -78.67)
# The results are as we expected since epsilon is the measurement error and eta is the error affecting the signal (mu)
# hence, sigma_eta has more influence on the estimation that sigma_epsilon (or we can say eta has more influence
# than epsilon)


################################################################################
# Question No. 1 part b
################################################################################
# AR(1) signal plus noise model with Gaussian and Student-t measurement noise
# the case with nu = 3 and sigma_e = sigma_eta = 1
################################################################################

# simulate the AR(1) time varying mean or signal (or trend)

# fix the length of the series 
n = 300

# simulate the noise of the signal (in the transition equation (SSF))
sigma_eta = 1
eta = sigma_eta * rnorm(n)

# simulate the measurement noise (in the observation equation)
sigma_e = 1

# q, signal-to-noise ratio, q = sigma^2_eta/sigma^2_e

q = sigma_eta^2/sigma_e^2

# case when the data are generated by a heavy tailed distribution  
# note that assuming s Student t distribution for the measurement error does break
# the assumption of KF which is optimal (in MMSE sense) for linear Gaussian models 
# nu is the degrees of freedom parameter (nu -> infinity then t -> Gaussian)

nu = 3
e = sqrt((nu-2)/nu) * sigma_e * rt(n,nu)
ts.plot(e)

# plot of eta (error affecting the signal) versus epsilon (e, measurement error)

ts.plot(e)
lines(eta, col = "red")

phi <- 0.8
mu = 0
mu[1] <- 0  # unconditional mean of an AR(1) process with no intercept
y = 0 
for(t in 1:(n-1)){
  mu[t+1] = phi * mu[t] + eta[t]
  y[t+1] =  mu[t+1] + e[t+1]
}
ts.plot(mu)

# y = mu + e

# Plot of y versus mu
ts.plot(y)
lines(mu,col="red")

# briefly check the correlation structure of the simulated process 
# to be sure that we are properly simulating an AR(1) signal plus noise process

acf(y, lag.max = 60, drop.lag.0 = FALSE)
pacf(y, lag.max = 60)
acf(mu, lag.max = 60, drop.lag.0 = FALSE)
pacf(mu, lag.max = 60)

# here, the acf and pacf of y suggest the AR(2) model behaviour but the acf and pacf of mu suggest AR(1) model
# behaviour


################################################################################
# KALMAN filter recursions (see Durbin and Koopman, 2001) 
################################################################################

# we have simulated y_t and an underlying signal mu_t  
# (observed and latent or unobserved in real data, respectively).
# Now we wish to estimate mu_{t|t-1} - note that we are 
# assuming that phi, sigma_e, sigma_eta, i.e. the static parameters, are known.
# In practice, they are estimated by maximising the Gaussian likelihood 
# (we shall see it) 

# initial conditions needed 
# mu_10, P_10 that we can set at the unconditional mean and variance of mu 


# KF recursion

# allocate space 
mu_pred <- 0  # this is mu_{t|t-1}
P <- 0        # this is P_{t|t-1} 
v <-0         # this  will be the innovation error 
K <- 0        # the Kalman gain 
F <- 0        # the conditional variance of v_t 
llk<-0        # the log-likelihood value 

# initialise the recursion 

# mu_pred[1] is set equal to its unconditional mean 
# P is the conditional variance of mu_t, let us initialise it as the signal-to-noise ratio

mu_pred[1] = 0
P[1] = (sigma_eta^2)/(1-phi^2)
llk[1]=0

# the recursion 

for(t in 1:(n-1)){
  v[t] = y[t] - mu_pred[t]; 
  F[t] = P[t] + sigma_e^2;
  K[t] = (phi * P[t])/F[t]
  P[t+1] = phi^2 * P[t] + sigma_eta^2 - K[t]*F[t]*K[t]
  mu_pred[t+1] = phi * mu_pred[t] + K[t]*v[t]
  llk[t] = - 0.5 * log(F[t] + (v[t]^2/F[t]))
}

llk = sum(llk)

ts.plot(y)
lines(mu, col = 'red')
lines(mu_pred, col = 'green')

ts.plot(y)
lines(mu, col = "red")

ts.plot(y)
lines(mu_pred, col = "green")

ts.plot(mu)
lines(mu_pred, col = "blue")
# in this case we have llk = -170.63, it has a slight improvement over the case where we did not use nu
# (the first case in part a)

################################################################################
# Question No. 1 part b
################################################################################
# AR(1) signal plus noise model with Gaussian and Student-t measurement noise
# the case with nu = 6 and sigma_e = sigma_eta = 1
################################################################################

# simulate the AR(1) time varying mean or signal (or trend)

# fix the length of the series 
n = 300

# simulate the noise of the signal (in the transition equation (SSF))
sigma_eta = 1
eta = sigma_eta * rnorm(n)

# simulate the measurement noise (in the observation equation)
sigma_e = 1

# q, signal-to-noise ratio, q = sigma^2_eta/sigma^2_e

q = sigma_eta^2/sigma_e^2

# case when the data are generated by a heavy tailed distribution  
# note that assuming s Student t distribution for the measurement error does break
# the assumption of KF which is optimal (in MMSE sense) for linear Gaussian models 
# nu is the degrees of freedom parameter (nu -> infinity then t -> Gaussian)

nu = 6
e = sqrt((nu-2)/nu) * sigma_e * rt(n,nu)
ts.plot(e)

# plot of eta (error affecting the signal) versus epsilon (e, measurement error)

ts.plot(e)
lines(eta, col = "red")

phi <- 0.8
mu = 0
mu[1] <- 0  # unconditional mean of an AR(1) process with no intercept
y = 0 
for(t in 1:(n-1)){
  mu[t+1] = phi * mu[t] + eta[t]
  y[t+1] =  mu[t+1] + e[t+1]
}
ts.plot(mu)

# y = mu + e

# Plot of y versus mu
ts.plot(y)
lines(mu,col="red")

# briefly check the correlation structure of the simulated process 
# to be sure that we are properly simulating an AR(1) signal plus noise process

acf(y, lag.max = 60, drop.lag.0 = FALSE)
pacf(y, lag.max = 60)
acf(mu, lag.max = 60, drop.lag.0 = FALSE)
pacf(mu, lag.max = 60)

# here, the acf and pacf of y and mu suggest the AR(1) model behaviour


################################################################################
# KALMAN filter recursions (see Durbin and Koopman, 2001) 
################################################################################

# we have simulated y_t and an underlying signal mu_t  
# (observed and latent or unobserved in real data, respectively).
# Now we wish to estimate mu_{t|t-1} - note that we are 
# assuming that phi, sigma_e, sigma_eta, i.e. the static parameters, are known.
# In practice, they are estimated by maximising the Gaussian likelihood 
# (we shall see it) 

# initial conditions needed 
# mu_10, P_10 that we can set at the unconditional mean and variance of mu 


# KF recursion

# allocate space 
mu_pred <- 0  # this is mu_{t|t-1}
P <- 0        # this is P_{t|t-1} 
v <-0         # this  will be the innovation error 
K <- 0        # the Kalman gain 
F <- 0        # the conditional variance of v_t 
llk<-0        # the log-likelihood value 

# initialise the recursion 

# mu_pred[1] is set equal to its unconditional mean 
# P is the conditional variance of mu_t, let us initialise it as the signal-to-noise ratio

mu_pred[1] = 0
P[1] = (sigma_eta^2)/(1-phi^2)
llk[1]=0

# the recursion 

for(t in 1:(n-1)){
  v[t] = y[t] - mu_pred[t]; 
  F[t] = P[t] + sigma_e^2;
  K[t] = (phi * P[t])/F[t]
  P[t+1] = phi^2 * P[t] + sigma_eta^2 - K[t]*F[t]*K[t]
  mu_pred[t+1] = phi * mu_pred[t] + K[t]*v[t]
  llk[t] = - 0.5 * log(F[t] + (v[t]^2/F[t]))
}

llk = sum(llk)

ts.plot(y)
lines(mu, col = 'red')
lines(mu_pred, col = 'green')

ts.plot(y)
lines(mu, col = "red")

ts.plot(y)
lines(mu_pred, col = "green")

ts.plot(mu)
lines(mu_pred, col = "blue")
# in this case we have slightly better results than the case with nu = 3, as llk = -168.62

################################################################################
# Question No. 1 part b
################################################################################
# AR(1) signal plus noise model with Gaussian and Student-t measurement noise
# the case with nu = 12 and sigma_e = sigma_eta = 1
################################################################################

# simulate the AR(1) time varying mean or signal (or trend)

# fix the length of the series 
n = 300

# simulate the noise of the signal (in the transition equation (SSF))
sigma_eta = 1
eta = sigma_eta * rnorm(n)

# simulate the measurement noise (in the observation equation)
sigma_e = 1

# q, signal-to-noise ratio, q = sigma^2_eta/sigma^2_e

q = sigma_eta^2/sigma_e^2

# case when the data are generated by a heavy tailed distribution  
# note that assuming s Student t distribution for the measurement error does break
# the assumption of KF which is optimal (in MMSE sense) for linear Gaussian models 
# nu is the degrees of freedom parameter (nu -> infinity then t -> Gaussian)

nu = 12
e = sqrt((nu-2)/nu) * sigma_e * rt(n,nu)
ts.plot(e)

# plot of eta (error affecting the signal) versus epsilon (e, measurement error)

ts.plot(e)
lines(eta, col = "red")

phi <- 0.8
mu = 0
mu[1] <- 0  # unconditional mean of an AR(1) process with no intercept
y = 0 
for(t in 1:(n-1)){
  mu[t+1] = phi * mu[t] + eta[t]
  y[t+1] =  mu[t+1] + e[t+1]
}
ts.plot(mu)

# y = mu + e

# Plot of y versus mu
ts.plot(y)
lines(mu,col="red")

# briefly check the correlation structure of the simulated process 
# to be sure that we are properly simulating an AR(1) signal plus noise process

acf(y, lag.max = 60, drop.lag.0 = FALSE)
pacf(y, lag.max = 60)
acf(mu, lag.max = 60, drop.lag.0 = FALSE)
pacf(mu, lag.max = 60)

# here, the acf and pacf of y and mu suggest the AR(1) model


################################################################################
# KALMAN filter recursions (see Durbin and Koopman, 2001) 
################################################################################

# we have simulated y_t and an underlying signal mu_t  
# (observed and latent or unobserved in real data, respectively).
# Now we wish to estimate mu_{t|t-1} - note that we are 
# assuming that phi, sigma_e, sigma_eta, i.e. the static parameters, are known.
# In practice, they are estimated by maximising the Gaussian likelihood 
# (we shall see it) 

# initial conditions needed 
# mu_10, P_10 that we can set at the unconditional mean and variance of mu 


# KF recursion

# allocate space 
mu_pred <- 0  # this is mu_{t|t-1}
P <- 0        # this is P_{t|t-1} 
v <-0         # this  will be the innovation error 
K <- 0        # the Kalman gain 
F <- 0        # the conditional variance of v_t 
llk<-0        # the log-likelihood value 

# initialise the recursion 

# mu_pred[1] is set equal to its unconditional mean 
# P is the conditional variance of mu_t, let us initialise it as the signal-to-noise ratio

mu_pred[1] = 0
P[1] = (sigma_eta^2)/(1-phi^2)
llk[1]=0

# the recursion 

for(t in 1:(n-1)){
  v[t] = y[t] - mu_pred[t]; 
  F[t] = P[t] + sigma_e^2;
  K[t] = (phi * P[t])/F[t]
  P[t+1] = phi^2 * P[t] + sigma_eta^2 - K[t]*F[t]*K[t]
  mu_pred[t+1] = phi * mu_pred[t] + K[t]*v[t]
  llk[t] = - 0.5 * log(F[t] + (v[t]^2/F[t]))
}

llk = sum(llk)

ts.plot(y)
lines(mu, col = 'red')
lines(mu_pred, col = 'green')

ts.plot(y)
lines(mu, col = "red")

ts.plot(y)
lines(mu_pred, col = "green")

ts.plot(mu)
lines(mu_pred, col = "blue")
# in this case we have llk = -173.06, hence, the estimation got worse compared to the cases with nu = 3 and nu = 6


################################################################################
# Question No. 1 part b
################################################################################
# AR(1) signal plus noise model with Gaussian and Student-t measurement noise
# the case with nu = 28 and sigma_e = sigma_eta = 1
################################################################################

# simulate the AR(1) time varying mean or signal (or trend)

# fix the length of the series 
n = 300

# simulate the noise of the signal (in the transition equation (SSF))
sigma_eta = 1
eta = sigma_eta * rnorm(n)

# simulate the measurement noise (in the observation equation)
sigma_e = 1

# q, signal-to-noise ratio, q = sigma^2_eta/sigma^2_e

q = sigma_eta^2/sigma_e^2

# case when the data are generated by a heavy tailed distribution  
# note that assuming s Student t distribution for the measurement error does break
# the assumption of KF which is optimal (in MMSE sense) for linear Gaussian models 
# nu is the degrees of freedom parameter (nu -> infinity then t -> Gaussian)

nu = 28
e = sqrt((nu-2)/nu) * sigma_e * rt(n,nu)
ts.plot(e)

# plot of eta (error affecting the signal) versus epsilon (e, measurement error)

ts.plot(e)
lines(eta, col = "red")

phi <- 0.8
mu = 0
mu[1] <- 0  # unconditional mean of an AR(1) process with no intercept
y = 0 
for(t in 1:(n-1)){
  mu[t+1] = phi * mu[t] + eta[t]
  y[t+1] =  mu[t+1] + e[t+1]
}
ts.plot(mu)

# y = mu + e

# Plot of y versus mu
ts.plot(y)
lines(mu,col="red")

# briefly check the correlation structure of the simulated process 
# to be sure that we are properly simulating an AR(1) signal plus noise process

acf(y, lag.max = 60, drop.lag.0 = FALSE)
pacf(y, lag.max = 60)
acf(mu, lag.max = 60, drop.lag.0 = FALSE)
pacf(mu, lag.max = 60)

# here, the acf and pacf of y suggest the AR(2) model behaviour but the acf and pacf of mu suggest AR(1) model
# behaviour


################################################################################
# KALMAN filter recursions (see Durbin and Koopman, 2001) 
################################################################################

# we have simulated y_t and an underlying signal mu_t  
# (observed and latent or unobserved in real data, respectively).
# Now we wish to estimate mu_{t|t-1} - note that we are 
# assuming that phi, sigma_e, sigma_eta, i.e. the static parameters, are known.
# In practice, they are estimated by maximising the Gaussian likelihood 
# (we shall see it) 

# initial conditions needed 
# mu_10, P_10 that we can set at the unconditional mean and variance of mu 


# KF recursion

# allocate space 
mu_pred <- 0  # this is mu_{t|t-1}
P <- 0        # this is P_{t|t-1} 
v <-0         # this  will be the innovation error 
K <- 0        # the Kalman gain 
F <- 0        # the conditional variance of v_t 
llk<-0        # the log-likelihood value 

# initialise the recursion 

# mu_pred[1] is set equal to its unconditional mean 
# P is the conditional variance of mu_t, let us initialise it as the signal-to-noise ratio

mu_pred[1] = 0
P[1] = (sigma_eta^2)/(1-phi^2)
llk[1]=0

# the recursion 

for(t in 1:(n-1)){
  v[t] = y[t] - mu_pred[t]; 
  F[t] = P[t] + sigma_e^2;
  K[t] = (phi * P[t])/F[t]
  P[t+1] = phi^2 * P[t] + sigma_eta^2 - K[t]*F[t]*K[t]
  mu_pred[t+1] = phi * mu_pred[t] + K[t]*v[t]
  llk[t] = - 0.5 * log(F[t] + (v[t]^2/F[t]))
}

llk = sum(llk)

ts.plot(y)
lines(mu, col = 'red')
lines(mu_pred, col = 'green')

ts.plot(y)
lines(mu, col = "red")

ts.plot(y)
lines(mu_pred, col = "green")

ts.plot(mu)
lines(mu_pred, col = "blue")
# in this case we have llk = -170.43, it has the worst estimation compared to the cases with nu = 3 and nu = 6,
# and and better estimation than the case with nu = 12


################################################################################
# Question No. 1 part b
################################################################################
# AR(1) signal plus noise model with Gaussian and Student-t measurement noise
# the case with nu = 200 and sigma_e = sigma_eta = 1
################################################################################

# simulate the AR(1) time varying mean or signal (or trend)

# fix the length of the series 
n = 300

# simulate the noise of the signal (in the transition equation (SSF))
sigma_eta = 1
eta = sigma_eta * rnorm(n)

# simulate the measurement noise (in the observation equation)
sigma_e = 1

# q, signal-to-noise ratio, q = sigma^2_eta/sigma^2_e

q = sigma_eta^2/sigma_e^2

# case when the data are generated by a heavy tailed distribution  
# note that assuming s Student t distribution for the measurement error does break
# the assumption of KF which is optimal (in MMSE sense) for linear Gaussian models 
# nu is the degrees of freedom parameter (nu -> infinity then t -> Gaussian)

nu = 200
e = sqrt((nu-2)/nu) * sigma_e * rt(n,nu)
ts.plot(e)

# plot of eta (error affecting the signal) versus epsilon (e, measurement error)

ts.plot(e)
lines(eta, col = "red")

phi <- 0.8
mu = 0
mu[1] <- 0  # unconditional mean of an AR(1) process with no intercept
y = 0 
for(t in 1:(n-1)){
  mu[t+1] = phi * mu[t] + eta[t]
  y[t+1] =  mu[t+1] + e[t+1]
}
ts.plot(mu)

# y = mu + e

# Plot of y versus mu
ts.plot(y)
lines(mu,col="red")

# briefly check the correlation structure of the simulated process 
# to be sure that we are properly simulating an AR(1) signal plus noise process

acf(y, lag.max = 60, drop.lag.0 = FALSE)
pacf(y, lag.max = 60)
acf(mu, lag.max = 60, drop.lag.0 = FALSE)
pacf(mu, lag.max = 60)

# here, the acf and pacf of y suggest the AR(2) model behaviour but the acf and pacf of mu suggest AR(1) model
# behaviour


################################################################################
# KALMAN filter recursions (see Durbin and Koopman, 2001) 
################################################################################

# we have simulated y_t and an underlying signal mu_t  
# (observed and latent or unobserved in real data, respectively).
# Now we wish to estimate mu_{t|t-1} - note that we are 
# assuming that phi, sigma_e, sigma_eta, i.e. the static parameters, are known.
# In practice, they are estimated by maximising the Gaussian likelihood 
# (we shall see it) 

# initial conditions needed 
# mu_10, P_10 that we can set at the unconditional mean and variance of mu 


# KF recursion

# allocate space 
mu_pred <- 0  # this is mu_{t|t-1}
P <- 0        # this is P_{t|t-1} 
v <-0         # this  will be the innovation error 
K <- 0        # the Kalman gain 
F <- 0        # the conditional variance of v_t 
llk<-0        # the log-likelihood value 

# initialise the recursion 

# mu_pred[1] is set equal to its unconditional mean 
# P is the conditional variance of mu_t, let us initialise it as the signal-to-noise ratio

mu_pred[1] = 0
P[1] = (sigma_eta^2)/(1-phi^2)
llk[1]=0

# the recursion 

for(t in 1:(n-1)){
  v[t] = y[t] - mu_pred[t]; 
  F[t] = P[t] + sigma_e^2;
  K[t] = (phi * P[t])/F[t]
  P[t+1] = phi^2 * P[t] + sigma_eta^2 - K[t]*F[t]*K[t]
  mu_pred[t+1] = phi * mu_pred[t] + K[t]*v[t]
  llk[t] = - 0.5 * log(F[t] + (v[t]^2/F[t]))
}

llk = sum(llk)

ts.plot(y)
lines(mu, col = 'red')
lines(mu_pred, col = 'green')

ts.plot(y)
lines(mu, col = "red")

ts.plot(y)
lines(mu_pred, col = "green")

ts.plot(mu)
lines(mu_pred, col = "blue")
# in this case we have llk = -168.31, it has the best estimation than the previous cases
# we have llk = (-170.63, -168.62, -173.06, -170.43, -168.31) corresponding respectively to nu = (3, 6, 12, 28, 200)
# the case with nu = 12 has the worst estimation



################################################################################
# Question No. 1 part c
################################################################################
# LLM model with Gaussian and Student-t measurement noise
# the case with sigma_e = 1 and sigma_eta = 0.5
################################################################################

# simulate the LLM time varying mean or signal (or trend)

# fix the length of the series 
n = 300

# simulate the noise of the signal (in the transition equation (SSF))
sigma_eta = .5
eta = sigma_eta * rnorm(n)

# simulate the measurement noise (in the observation equation)
sigma_e = 1
e = sigma_e * rnorm(n)

# q, signal-to-noise ratio, q = sigma^2_eta/sigma^2_e

q = sigma_eta^2/sigma_e^2


# plot of eta (error affecting the signal) versus epsilon (e, measurement error)

ts.plot(e)
lines(eta, col = "red")

mu = 0
mu[1] = 0  # initial condition  
y = 0 

for(t in 1:(n-1)){
  mu[t+1] = mu[t] + eta[t]
  y[t+1] =  mu[t+1] + e[t+1]
}
ts.plot(mu)

# y = mu + e

# Plot of y versus mu
ts.plot(y)
lines(mu,col="red")

# briefly check the correlation structure of the simulated process 
# to be sure that we are properly simulating a LLM process

acf(y, lag.max = 80, drop.lag.0 = FALSE)
pacf(y, lag.max = 80)
acf(mu, lag.max = 80, drop.lag.0 = FALSE)
pacf(mu, lag.max = 80)

# as expected, the acf and pacf of y and mu suggest non-stationary process

dy = diff(y)
ts.plot(dy)
acf(dy,lag.max = 100, drop.lag.0 = FALSE)
pacf(dy,lag.max = 100, drop.lag.0 = FALSE)

# while the reduced form of (1-L)y_t is MA(1),  (1-L)mu_t 
# is a white noise process, actually, BY CONSTRUCTION, is eta_t

dmu = diff(mu)
ts.plot(dmu)
lines(eta, col = "red")
acf(dmu,lag.max = 100, drop.lag.0 = FALSE)
pacf(dmu,lag.max = 100, drop.lag.0 = FALSE)

# hence, our LLM process is an ARIMA(0,1,1) process

################################################################################
# KALMAN filter recursions (see Durbin and Koopman, 2001) 
################################################################################

# we have simulated y_t and an underlying signal mu_t  
# (observed and latent or unobserved in real data, respectively).
# Now we wish to estimate mu_{t|t-1} - note that we are 
# assuming that phi, sigma_e, sigma_eta, i.e. the static parameters, are known.
# In practice, they are estimated by maximising the Gaussian likelihood 
# (we shall see it) 

# initial conditions needed 
# mu_10, P_10 that we can set at the unconditional mean and variance of mu 


# KF recursion

# allocate space 
mu_pred <- 0  # this is mu_{t|t-1}
P <- 0        # this is P_{t|t-1} 
v <-0         # this  will be the innovation error 
K <- 0        # the Kalman gain 
F <- 0        # the conditional variance of v_t 
llk<-0        # the log-likelihood value 

# initialise the recursion 

# mu_pred[1] is set equal to its unconditional mean 
# P is the conditional variance of mu_t, let us initialise it as the signal-to-noise ratio

mu_pred[1] = 0
P[1] = (sigma_eta^2)/(1-phi^2)
llk[1]=0

# the recursion 

for(t in 1:(n-1)){
  v[t] = y[t] - mu_pred[t]; 
  F[t] = P[t] + sigma_e^2;
  K[t] =  P[t] / F[t]
  P[t+1] = P[t] + sigma_eta^2 - K[t]*F[t]*K[t]
  mu_pred[t+1] = mu_pred[t] + K[t]*v[t]
  llk[t] = - 0.5 * log(F[t] + (v[t]^2/F[t]))
}

llk = sum(llk)

ts.plot(y)
lines(mu, col = 'red')
lines(mu_pred, col = 'green')

ts.plot(y)
lines(mu, col = "red")

ts.plot(y)
lines(mu_pred, col = "green")

ts.plot(mu)
lines(mu_pred, col = "blue")

# in this case, we have llk = -137.1


################################################################################
# Question No. 1 part c
################################################################################
# LLM model with Gaussian and Student-t measurement noise
# the case with sigma_e = 1 and sigma_eta = 0.5 and nu = 3
################################################################################

# simulate the LLM time varying mean or signal (or trend)

# fix the length of the series 
n = 300

# simulate the noise of the signal (in the transition equation (SSF))
sigma_eta = .5
eta = sigma_eta * rnorm(n)

# simulate the measurement noise (in the observation equation)
sigma_e = 1

# q, signal-to-noise ratio, q = sigma^2_eta/sigma^2_e

q = sigma_eta^2/sigma_e^2

# case when the data are generated by a heavy tailed distribution  
# note that assuming s Student t distribution for the measurement error does break
# the assumption of KF which is optimal (in MMSE sense) for linear Gaussian models 
# nu is the degrees of freedom parameter (nu -> infinity then t -> Gaussian)

nu = 3
e = sqrt((nu-2)/nu) * sigma_e * rt(n,nu)
ts.plot(e)

# plot of eta (error affecting the signal) versus epsilon (e, measurement error)

ts.plot(e)
lines(eta, col = "red")

mu = 0
mu[1] = 0  # initial condition  
y = 0 

for(t in 1:(n-1)){
  mu[t+1] = mu[t] + eta[t]
  y[t+1] =  mu[t+1] + e[t+1]
}
ts.plot(mu)

# y = mu + e

# Plot of y versus mu
ts.plot(y)
lines(mu,col="red")

# briefly check the correlation structure of the simulated process 
# to be sure that we are properly simulating a LLM process

acf(y, lag.max = 80, drop.lag.0 = FALSE)
pacf(y, lag.max = 80)
acf(mu, lag.max = 80, drop.lag.0 = FALSE)
pacf(mu, lag.max = 80)

# as expected, the acf and pacf of y and mu suggest non-stationary process

dy = diff(y)
ts.plot(dy)
acf(dy,lag.max = 100, drop.lag.0 = FALSE)
pacf(dy,lag.max = 100, drop.lag.0 = FALSE)

# while the reduced form of (1-L)y_t is MA(1),  (1-L)mu_t 
# is a white noise process, actually, BY CONSTRUCTION, is eta_t

dmu = diff(mu)
ts.plot(dmu)
lines(eta, col = "red")
acf(dmu,lag.max = 100, drop.lag.0 = FALSE)
pacf(dmu,lag.max = 100, drop.lag.0 = FALSE)

# hence, our LLM process is an ARIMA(0,1,1) process

################################################################################
# KALMAN filter recursions (see Durbin and Koopman, 2001) 
################################################################################

# we have simulated y_t and an underlying signal mu_t  
# (observed and latent or unobserved in real data, respectively).
# Now we wish to estimate mu_{t|t-1} - note that we are 
# assuming that phi, sigma_e, sigma_eta, i.e. the static parameters, are known.
# In practice, they are estimated by maximising the Gaussian likelihood 
# (we shall see it) 

# initial conditions needed 
# mu_10, P_10 that we can set at the unconditional mean and variance of mu 


# KF recursion

# allocate space 
mu_pred <- 0  # this is mu_{t|t-1}
P <- 0        # this is P_{t|t-1} 
v <-0         # this  will be the innovation error 
K <- 0        # the Kalman gain 
F <- 0        # the conditional variance of v_t 
llk<-0        # the log-likelihood value 

# initialise the recursion 

# mu_pred[1] is set equal to its unconditional mean 
# P is the conditional variance of mu_t, let us initialise it as the signal-to-noise ratio

mu_pred[1] = 0
P[1] = (sigma_eta^2)/(1-phi^2)
llk[1]=0

# the recursion 

for(t in 1:(n-1)){
  v[t] = y[t] - mu_pred[t]; 
  F[t] = P[t] + sigma_e^2;
  K[t] =  P[t] / F[t]
  P[t+1] = P[t] + sigma_eta^2 - K[t]*F[t]*K[t]
  mu_pred[t+1] = mu_pred[t] + K[t]*v[t]
  llk[t] = - 0.5 * log(F[t] + (v[t]^2/F[t]))
}

llk = sum(llk)

ts.plot(y)
lines(mu, col = 'red')
lines(mu_pred, col = 'green')

ts.plot(y)
lines(mu, col = "red")

ts.plot(y)
lines(mu_pred, col = "green")

ts.plot(mu)
lines(mu_pred, col = "blue")

# in this case, we have llk = -125.06 which is better than the previous case (without nu)


################################################################################
# Question No. 1 part c
################################################################################
# LLM model with Gaussian and Student-t measurement noise
# the case with sigma_e = 1 and sigma_eta = 0.5 and nu = 6
################################################################################

# simulate the LLM time varying mean or signal (or trend)

# fix the length of the series 
n = 300

# simulate the noise of the signal (in the transition equation (SSF))
sigma_eta = .5
eta = sigma_eta * rnorm(n)

# simulate the measurement noise (in the observation equation)
sigma_e = 1

# q, signal-to-noise ratio, q = sigma^2_eta/sigma^2_e

q = sigma_eta^2/sigma_e^2

# case when the data are generated by a heavy tailed distribution  
# note that assuming s Student t distribution for the measurement error does break
# the assumption of KF which is optimal (in MMSE sense) for linear Gaussian models 
# nu is the degrees of freedom parameter (nu -> infinity then t -> Gaussian)

nu = 6
e = sqrt((nu-2)/nu) * sigma_e * rt(n,nu)
ts.plot(e)

# plot of eta (error affecting the signal) versus epsilon (e, measurement error)

ts.plot(e)
lines(eta, col = "red")

mu = 0
mu[1] = 0  # initial condition  
y = 0 

for(t in 1:(n-1)){
  mu[t+1] = mu[t] + eta[t]
  y[t+1] =  mu[t+1] + e[t+1]
}
ts.plot(mu)

# y = mu + e

# Plot of y versus mu
ts.plot(y)
lines(mu,col="red")

# briefly check the correlation structure of the simulated process 
# to be sure that we are properly simulating a LLM process

acf(y, lag.max = 80, drop.lag.0 = FALSE)
pacf(y, lag.max = 80)
acf(mu, lag.max = 80, drop.lag.0 = FALSE)
pacf(mu, lag.max = 80)

# as expected, the acf and pacf of y and mu suggest non-stationary process

dy = diff(y)
ts.plot(dy)
acf(dy,lag.max = 100, drop.lag.0 = FALSE)
pacf(dy,lag.max = 100, drop.lag.0 = FALSE)

# while the reduced form of (1-L)y_t is ARMA(1, 1),  (1-L)mu_t 
# is a white noise process, actually, BY CONSTRUCTION, is eta_t

dmu = diff(mu)
ts.plot(dmu)
lines(eta, col = "red")
acf(dmu,lag.max = 100, drop.lag.0 = FALSE)
pacf(dmu,lag.max = 100, drop.lag.0 = FALSE)

# hence, our LLM process is an ARIMA(1,1,1) process

################################################################################
# KALMAN filter recursions (see Durbin and Koopman, 2001) 
################################################################################

# we have simulated y_t and an underlying signal mu_t  
# (observed and latent or unobserved in real data, respectively).
# Now we wish to estimate mu_{t|t-1} - note that we are 
# assuming that phi, sigma_e, sigma_eta, i.e. the static parameters, are known.
# In practice, they are estimated by maximising the Gaussian likelihood 
# (we shall see it) 

# initial conditions needed 
# mu_10, P_10 that we can set at the unconditional mean and variance of mu 


# KF recursion

# allocate space 
mu_pred <- 0  # this is mu_{t|t-1}
P <- 0        # this is P_{t|t-1} 
v <-0         # this  will be the innovation error 
K <- 0        # the Kalman gain 
F <- 0        # the conditional variance of v_t 
llk<-0        # the log-likelihood value 

# initialise the recursion 

# mu_pred[1] is set equal to its unconditional mean 
# P is the conditional variance of mu_t, let us initialise it as the signal-to-noise ratio

mu_pred[1] = 0
P[1] = (sigma_eta^2)/(1-phi^2)
llk[1]=0

# the recursion 

for(t in 1:(n-1)){
  v[t] = y[t] - mu_pred[t]; 
  F[t] = P[t] + sigma_e^2;
  K[t] =  P[t] / F[t]
  P[t+1] = P[t] + sigma_eta^2 - K[t]*F[t]*K[t]
  mu_pred[t+1] = mu_pred[t] + K[t]*v[t]
  llk[t] = - 0.5 * log(F[t] + (v[t]^2/F[t]))
}

llk = sum(llk)

ts.plot(y)
lines(mu, col = 'red')
lines(mu_pred, col = 'green')

ts.plot(y)
lines(mu, col = "red")

ts.plot(y)
lines(mu_pred, col = "green")

ts.plot(mu)
lines(mu_pred, col = "blue")

# in this case, we have llk = -125.25 which is worse than the previous case where nu = 3 but still better than the
# original case

################################################################################
# Question No. 1 part c
################################################################################
# LLM model with Gaussian and Student-t measurement noise
# the case with sigma_e = 1 and sigma_eta = 0.5 and nu = 12
################################################################################

# fix the length of the series 
n = 300

# simulate the noise of the signal (in the transition equation (SSF))
sigma_eta = .5
eta = sigma_eta * rnorm(n)

# simulate the measurement noise (in the observation equation)
sigma_e = 1

# q, signal-to-noise ratio, q = sigma^2_eta/sigma^2_e

q = sigma_eta^2/sigma_e^2

# case when the data are generated by a heavy tailed distribution  
# note that assuming s Student t distribution for the measurement error does break
# the assumption of KF which is optimal (in MMSE sense) for linear Gaussian models 
# nu is the degrees of freedom parameter (nu -> infinity then t -> Gaussian)

nu = 12
e = sqrt((nu-2)/nu) * sigma_e * rt(n,nu)
ts.plot(e)

# plot of eta (error affecting the signal) versus epsilon (e, measurement error)

ts.plot(e)
lines(eta, col = "red")

mu = 0
mu[1] = 0  # initial condition  
y = 0 

for(t in 1:(n-1)){
  mu[t+1] = mu[t] + eta[t]
  y[t+1] =  mu[t+1] + e[t+1]
}
ts.plot(mu)

# y = mu + e

# Plot of y versus mu
ts.plot(y)
lines(mu,col="red")

# briefly check the correlation structure of the simulated process 
# to be sure that we are properly simulating a LLM process

acf(y, lag.max = 80, drop.lag.0 = FALSE)
pacf(y, lag.max = 80)
acf(mu, lag.max = 80, drop.lag.0 = FALSE)
pacf(mu, lag.max = 80)

# as expected, the acf and pacf of y and mu suggest non-stationary process

dy = diff(y)
ts.plot(dy)
acf(dy,lag.max = 100, drop.lag.0 = FALSE)
pacf(dy,lag.max = 100, drop.lag.0 = FALSE)

# while the reduced form of (1-L)y_t is ARMA(1, 1),  (1-L)mu_t 
# is a white noise process, actually, BY CONSTRUCTION, is eta_t

dmu = diff(mu)
ts.plot(dmu)
lines(eta, col = "red")
acf(dmu,lag.max = 100, drop.lag.0 = FALSE)
pacf(dmu,lag.max = 100, drop.lag.0 = FALSE)

# hence, our LLM process is an ARIMA(1,1,1) process

################################################################################
# KALMAN filter recursions (see Durbin and Koopman, 2001) 
################################################################################

# we have simulated y_t and an underlying signal mu_t  
# (observed and latent or unobserved in real data, respectively).
# Now we wish to estimate mu_{t|t-1} - note that we are 
# assuming that phi, sigma_e, sigma_eta, i.e. the static parameters, are known.
# In practice, they are estimated by maximising the Gaussian likelihood 
# (we shall see it) 

# initial conditions needed 
# mu_10, P_10 that we can set at the unconditional mean and variance of mu 


# KF recursion

# allocate space 
mu_pred <- 0  # this is mu_{t|t-1}
P <- 0        # this is P_{t|t-1} 
v <-0         # this  will be the innovation error 
K <- 0        # the Kalman gain 
F <- 0        # the conditional variance of v_t 
llk<-0        # the log-likelihood value 

# initialise the recursion 

# mu_pred[1] is set equal to its unconditional mean 
# P is the conditional variance of mu_t, let us initialise it as the signal-to-noise ratio

mu_pred[1] = 0
P[1] = (sigma_eta^2)/(1-phi^2)
llk[1]=0

# the recursion 

for(t in 1:(n-1)){
  v[t] = y[t] - mu_pred[t]; 
  F[t] = P[t] + sigma_e^2;
  K[t] =  P[t] / F[t]
  P[t+1] = P[t] + sigma_eta^2 - K[t]*F[t]*K[t]
  mu_pred[t+1] = mu_pred[t] + K[t]*v[t]
  llk[t] = - 0.5 * log(F[t] + (v[t]^2/F[t]))
}

llk = sum(llk)

ts.plot(y)
lines(mu, col = 'red')
lines(mu_pred, col = 'green')

ts.plot(y)
lines(mu, col = "red")

ts.plot(y)
lines(mu_pred, col = "green")

ts.plot(mu)
lines(mu_pred, col = "blue")

# in this case, we have llk = -129.7 which is worse than the last two cases where nu = 3 and nu = 6



################################################################################
# Question No. 1 part c
################################################################################
# LLM model with Gaussian and Student-t measurement noise
# the case with sigma_e = 1 and sigma_eta = 0.5 and nu = 200
################################################################################

# fix the length of the series 
n = 300

# simulate the noise of the signal (in the transition equation (SSF))
sigma_eta = .5
eta = sigma_eta * rnorm(n)

# simulate the measurement noise (in the observation equation)
sigma_e = 1

# q, signal-to-noise ratio, q = sigma^2_eta/sigma^2_e

q = sigma_eta^2/sigma_e^2

# case when the data are generated by a heavy tailed distribution  
# note that assuming s Student t distribution for the measurement error does break
# the assumption of KF which is optimal (in MMSE sense) for linear Gaussian models 
# nu is the degrees of freedom parameter (nu -> infinity then t -> Gaussian)

nu = 200
e = sqrt((nu-2)/nu) * sigma_e * rt(n,nu)
ts.plot(e)

# plot of eta (error affecting the signal) versus epsilon (e, measurement error)

ts.plot(e)
lines(eta, col = "red")

mu = 0
mu[1] = 0  # initial condition  
y = 0 

for(t in 1:(n-1)){
  mu[t+1] = mu[t] + eta[t]
  y[t+1] =  mu[t+1] + e[t+1]
}
ts.plot(mu)

# y = mu + e

# Plot of y versus mu
ts.plot(y)
lines(mu,col="red")

# briefly check the correlation structure of the simulated process 
# to be sure that we are properly simulating a LLM process

acf(y, lag.max = 80, drop.lag.0 = FALSE)
pacf(y, lag.max = 80)
acf(mu, lag.max = 80, drop.lag.0 = FALSE)
pacf(mu, lag.max = 80)

# as expected, the acf and pacf of y and mu suggest non-stationary process

dy = diff(y)
ts.plot(dy)
acf(dy,lag.max = 100, drop.lag.0 = FALSE)
pacf(dy,lag.max = 100, drop.lag.0 = FALSE)

# while the reduced form of (1-L)y_t is ARMA(1, 1),  (1-L)mu_t 
# is a white noise process, actually, BY CONSTRUCTION, is eta_t

dmu = diff(mu)
ts.plot(dmu)
lines(eta, col = "red")
acf(dmu,lag.max = 100, drop.lag.0 = FALSE)
pacf(dmu,lag.max = 100, drop.lag.0 = FALSE)

# hence, our LLM process is an ARIMA(1,1,1) process

################################################################################
# KALMAN filter recursions (see Durbin and Koopman, 2001) 
################################################################################

# we have simulated y_t and an underlying signal mu_t  
# (observed and latent or unobserved in real data, respectively).
# Now we wish to estimate mu_{t|t-1} - note that we are 
# assuming that phi, sigma_e, sigma_eta, i.e. the static parameters, are known.
# In practice, they are estimated by maximising the Gaussian likelihood 
# (we shall see it) 

# initial conditions needed 
# mu_10, P_10 that we can set at the unconditional mean and variance of mu 


# KF recursion

# allocate space 
mu_pred <- 0  # this is mu_{t|t-1}
P <- 0        # this is P_{t|t-1} 
v <-0         # this  will be the innovation error 
K <- 0        # the Kalman gain 
F <- 0        # the conditional variance of v_t 
llk<-0        # the log-likelihood value 

# initialise the recursion 

# mu_pred[1] is set equal to its unconditional mean 
# P is the conditional variance of mu_t, let us initialise it as the signal-to-noise ratio

mu_pred[1] = 0
P[1] = (sigma_eta^2)/(1-phi^2)
llk[1]=0

# the recursion 

for(t in 1:(n-1)){
  v[t] = y[t] - mu_pred[t]; 
  F[t] = P[t] + sigma_e^2;
  K[t] =  P[t] / F[t]
  P[t+1] = P[t] + sigma_eta^2 - K[t]*F[t]*K[t]
  mu_pred[t+1] = mu_pred[t] + K[t]*v[t]
  llk[t] = - 0.5 * log(F[t] + (v[t]^2/F[t]))
}

llk = sum(llk)

ts.plot(y)
lines(mu, col = 'red')
lines(mu_pred, col = 'green')

ts.plot(y)
lines(mu, col = "red")

ts.plot(y)
lines(mu_pred, col = "green")

ts.plot(mu)
lines(mu_pred, col = "blue")

# in this case, we have llk = -134.89 which is worse than the last three cases where nu = 3, nu = 6, and nu = 12
# we have llk = (-125.06, -125.25, -129.7, -134.89) corresponding with nu = (3, 6, 12, 200) respectively



################################################################################
# Question No. 1 part c
################################################################################
# LLM model with Gaussian and Student-t measurement noise
# the case with sigma_e = 1 and sigma_eta = 0.5 and nu = 200000
################################################################################

# fix the length of the series 
n = 300

# simulate the noise of the signal (in the transition equation (SSF))
sigma_eta = .5
eta = sigma_eta * rnorm(n)

# simulate the measurement noise (in the observation equation)
sigma_e = 1

# q, signal-to-noise ratio, q = sigma^2_eta/sigma^2_e

q = sigma_eta^2/sigma_e^2

# case when the data are generated by a heavy tailed distribution  
# note that assuming s Student t distribution for the measurement error does break
# the assumption of KF which is optimal (in MMSE sense) for linear Gaussian models 
# nu is the degrees of freedom parameter (nu -> infinity then t -> Gaussian)

nu = 200000
e = sqrt((nu-2)/nu) * sigma_e * rt(n,nu)
ts.plot(e)

# plot of eta (error affecting the signal) versus epsilon (e, measurement error)

ts.plot(e)
lines(eta, col = "red")

mu = 0
mu[1] = 0  # initial condition  
y = 0 

for(t in 1:(n-1)){
  mu[t+1] = mu[t] + eta[t]
  y[t+1] =  mu[t+1] + e[t+1]
}
ts.plot(mu)

# y = mu + e

# Plot of y versus mu
ts.plot(y)
lines(mu,col="red")

# briefly check the correlation structure of the simulated process 
# to be sure that we are properly simulating an AR(1) signal plus noise process

acf(y, lag.max = 80, drop.lag.0 = FALSE)
pacf(y, lag.max = 80)
acf(mu, lag.max = 80, drop.lag.0 = FALSE)
pacf(mu, lag.max = 80)

# as expected, the acf and pacf of y and mu suggest non-stationary process

dy = diff(y)
ts.plot(dy)
acf(dy,lag.max = 100, drop.lag.0 = FALSE)
pacf(dy,lag.max = 100, drop.lag.0 = FALSE)

# while the reduced form of (1-L)y_t is ARMA(1, 1),  (1-L)mu_t 
# is a white noise process, actually, BY CONSTRUCTION, is eta_t

dmu = diff(mu)
ts.plot(dmu)
lines(eta, col = "red")
acf(dmu,lag.max = 100, drop.lag.0 = FALSE)
pacf(dmu,lag.max = 100, drop.lag.0 = FALSE)

# hence, our LLM process is an ARIMA(1,1,1) process

################################################################################
# KALMAN filter recursions (see Durbin and Koopman, 2001) 
################################################################################

# we have simulated y_t and an underlying signal mu_t  
# (observed and latent or unobserved in real data, respectively).
# Now we wish to estimate mu_{t|t-1} - note that we are 
# assuming that phi, sigma_e, sigma_eta, i.e. the static parameters, are known.
# In practice, they are estimated by maximising the Gaussian likelihood 
# (we shall see it) 

# initial conditions needed 
# mu_10, P_10 that we can set at the unconditional mean and variance of mu 


# KF recursion

# allocate space 
mu_pred <- 0  # this is mu_{t|t-1}
P <- 0        # this is P_{t|t-1} 
v <-0         # this  will be the innovation error 
K <- 0        # the Kalman gain 
F <- 0        # the conditional variance of v_t 
llk<-0        # the log-likelihood value 

# initialise the recursion 

# mu_pred[1] is set equal to its unconditional mean 
# P is the conditional variance of mu_t, let us initialise it as the signal-to-noise ratio

mu_pred[1] = 0
P[1] = (sigma_eta^2)/(1-phi^2)
llk[1]=0

# the recursion 

for(t in 1:(n-1)){
  v[t] = y[t] - mu_pred[t]; 
  F[t] = P[t] + sigma_e^2;
  K[t] =  P[t] / F[t]
  P[t+1] = P[t] + sigma_eta^2 - K[t]*F[t]*K[t]
  mu_pred[t+1] = mu_pred[t] + K[t]*v[t]
  llk[t] = - 0.5 * log(F[t] + (v[t]^2/F[t]))
}

llk = sum(llk)

ts.plot(y)
lines(mu, col = 'red')
lines(mu_pred, col = 'green')

ts.plot(y)
lines(mu, col = "red")

ts.plot(y)
lines(mu_pred, col = "green")

ts.plot(mu)
lines(mu_pred, col = "blue")

# in this case, we have llk = -132.44 with nu = 200000


################################################################################
# Question No. 2
################################################################################
# LLT model with Gaussian and Student-t measurement noise
# the case with sigma_e = 0.4, sigma_eta = 0.8, and sigma_xi = 1.5
################################################################################

# fix the length of the series 
n = 300

# simulate the noise of the signal (in the transition equation (SSF))
sigma_eta = .4
eta = sigma_eta * rnorm(n)

# simulate the measurement noise (in the observation equation)
sigma_e = .8
e = sigma_e * rnorm(n)
# q, signal-to-noise ratio, q = sigma^2_eta/sigma^2_e

q = sigma_eta^2/sigma_e^2


sigma_xi = 1.5
xi = sigma_xi * rnorm(n)

# plot of eta (error affecting the signal) versus epsilon (e, measurement error) and xi

ts.plot(e)
lines(eta, col = "red")
lines(xi, col = "blue")


mu = 0
mu[1] <- 0   
beta = 0
beta[1] <- 0   
y = 0 
for(t in 1:(n-1)){
  beta[t+1] = beta[t] + xi[t]
  mu[t+1] = beta[t] + mu[t] + eta[t]
  y[t+1] =  mu[t+1] + e[t+1]
}

ts.plot(mu)

# y = mu + e

# Plot of y versus mu
ts.plot(y)
lines(mu,col="red")

# briefly check the correlation structure of the simulated process 
# to be sure that we are properly simulating an AR(1) signal plus noise process

acf(y, lag.max = 80, drop.lag.0 = FALSE)
pacf(y, lag.max = 80)
acf(mu, lag.max = 80, drop.lag.0 = FALSE)
pacf(mu, lag.max = 80)

# as expected, the acf and pacf of y and mu suggest non-stationary process


# Differences
dy = diff(y)
ts.plot(dy)
acf(dy,lag.max = 100, drop.lag.0 = FALSE)
pacf(dy,lag.max = 100, drop.lag.0 = FALSE)



dmu = diff(mu)
ts.plot(dmu)
lines(eta, col = "red")
acf(dmu,lag.max = 100, drop.lag.0 = FALSE)
pacf(dmu,lag.max = 100, drop.lag.0 = FALSE)



dy = diff(diff(y))
ts.plot(dy)
acf(dy,lag.max = 100, drop.lag.0 = FALSE)
pacf(dy,lag.max = 100, drop.lag.0 = FALSE)

# while the reduced form of (1-L^2)y_t is ARMA(1,1),  (1-L^2)mu_t 
# is a MA(1) process

dmu = diff(diff(mu))
ts.plot(dmu)
lines(eta, col = "red")
acf(dmu,lag.max = 100, drop.lag.0 = FALSE)
pacf(dmu,lag.max = 100, drop.lag.0 = FALSE)

# hence, our LLT process is an ARIMA(1,2,1) process

################################################################################
# KALMAN filter recursions (see Durbin and Koopman, 2001) 
################################################################################

# we have simulated y_t and an underlying signal mu_t  
# (observed and latent or unobserved in real data, respectively).
# Now we wish to estimate mu_{t|t-1} - note that we are 
# assuming that phi, sigma_e, sigma_eta, i.e. the static parameters, are known.
# In practice, they are estimated by maximising the Gaussian likelihood 
# (we shall see it) 

# initial conditions needed 
# mu_10, P_10 that we can set at the unconditional mean and variance of mu 


# KF recursion

mu_pred <- 0 # this is mu_{t|t-1}
beta_pred <- 0 # this is beta_{t|t-1}
P11 <- 0 # this is P_{t|t-1}
P21 <- 0 
P12 <- 0 
P22 <- 0 
v <- 0 # this will be the innovation error
Kmu <- 0 # the Kalman gain
Kbeta <- 0 
F <- 0 # the conditional variance of v_t
llk<- 0 # the log-likelihood value


# initialise the recursion

# mu_pred[1] is set equal to its unconditional mean
# P is the conditional variance of mu_t, let us initialise it as the signal-to-noise ratio

mu_pred[1] <- 0
beta_pred[1] <- 0
P11[1] <- (sigma_eta^2)/(sigma_e^2)
P22[1] <- (sigma_xi^2)/(sigma_e^2)
P12[1] <- 0
P21[1] <- 0
llk[1] <- 0


# the recursion
for(t in 1:(n-1)){
  v[t] = y[t] - mu_pred[t];
  F[t] = P11[t] + sigma_e^2;
  Kmu[t] = (P11[t]+P21[t])/F[t]
  Kbeta[t] = P21[t]/F[t]
  P11[t+1] = P11[t] + P12[t] + P21[t] + P22[t] + sigma_eta^2 - (P11[t] + P21[t])^2/F[t]
  P12[t+1] = P12[t] + P22[t] - (P11[t]*P21[t]+ P21[t]^2)/F[t]
  P21[t+1] = P21[t] + P22[t] - (P11[t]*P21[t]+ P21[t]^2)/F[t]
  P22[t+1] = P22[t] + sigma_xi^2 - (P21[t]^2)/F[t]
  mu_pred[t+1] = mu_pred[t] + beta_pred[t] + Kmu[t]*v[t]
  beta_pred[t+1] = beta_pred[t] + Kbeta[t]*v[t]
  llk[t] = - 0.5 * log(F[t] + (v[t]^2/F[t]))
}

llk = sum(llk)

ts.plot(y)
lines(mu, col = 'red')
lines(mu_pred, col = 'green')

ts.plot(y)
lines(mu, col = "red")

ts.plot(y)
lines(mu_pred, col = "green")

ts.plot(mu)
lines(mu_pred, col = "blue")

# in this case, we have llk = -269.28, as expected, the estimation mu for the LLT model is not as good as the LLM
# (or AR(1) + noise) as there are more parameters to estimate in LLT, also xi affects the estimation while xi is 
# not involved in LLM (or AR(1) + noise)


################################################################################
# Question No. 2
################################################################################
# LLT model with Gaussian and Student-t measurement noise
# the case with sigma_e = 1, sigma_eta = 1, and sigma_xi = 1
################################################################################

# fix the length of the series 
n = 300

# simulate the noise of the signal (in the transition equation (SSF))
sigma_eta = 1
eta = sigma_eta * rnorm(n)

# simulate the measurement noise (in the observation equation)
sigma_e = 1
e = sigma_e * rnorm(n)
# q, signal-to-noise ratio, q = sigma^2_eta/sigma^2_e

q = sigma_eta^2/sigma_e^2

sigma_xi = 1
xi = sigma_xi * rnorm(n)

# plot of eta (error affecting the signal) versus epsilon (e, measurement error) and xi

ts.plot(e)
lines(eta, col = "red")
lines(xi, col = "blue")


mu = 0
mu[1] <- 0   
beta = 0
beta[1] <- 0   
y = 0 
for(t in 1:(n-1)){
  beta[t+1] = beta[t] + xi[t]
  mu[t+1] = beta[t] + mu[t] + eta[t]
  y[t+1] =  mu[t+1] + e[t+1]
}

ts.plot(mu)

# y = mu + e

# Plot of y versus mu
ts.plot(y)
lines(mu,col="red")

# briefly check the correlation structure of the simulated process 
# to be sure that we are properly simulating an AR(1) signal plus noise process

acf(y, lag.max = 80, drop.lag.0 = FALSE)
pacf(y, lag.max = 80)
acf(mu, lag.max = 80, drop.lag.0 = FALSE)
pacf(mu, lag.max = 80)

# as expected, the acf and pacf of y and mu suggest non-stationary process


# Differences
dy = diff(y)
ts.plot(dy)
acf(dy,lag.max = 100, drop.lag.0 = FALSE)
pacf(dy,lag.max = 100, drop.lag.0 = FALSE)



dmu = diff(mu)
ts.plot(dmu)
lines(eta, col = "red")
acf(dmu,lag.max = 100, drop.lag.0 = FALSE)
pacf(dmu,lag.max = 100, drop.lag.0 = FALSE)



dy = diff(diff(y))
ts.plot(dy)
acf(dy,lag.max = 100, drop.lag.0 = FALSE)
pacf(dy,lag.max = 100, drop.lag.0 = FALSE)

# while the reduced form of (1-L^2)y_t is ARMA(2,1),  (1-L^2)mu_t 
# is ARMA (1,1)

dmu = diff(diff(mu))
ts.plot(dmu)
lines(eta, col = "red")
acf(dmu,lag.max = 100, drop.lag.0 = FALSE)
pacf(dmu,lag.max = 100, drop.lag.0 = FALSE)

# hence, our LLT process is an ARIMA(2,2,1) process

################################################################################
# KALMAN filter recursions (see Durbin and Koopman, 2001) 
################################################################################

# we have simulated y_t and an underlying signal mu_t  
# (observed and latent or unobserved in real data, respectively).
# Now we wish to estimate mu_{t|t-1} - note that we are 
# assuming that phi, sigma_e, sigma_eta, i.e. the static parameters, are known.
# In practice, they are estimated by maximising the Gaussian likelihood 
# (we shall see it) 

# initial conditions needed 
# mu_10, P_10 that we can set at the unconditional mean and variance of mu 

# KF recursion

mu_pred <- 0 # this is mu_{t|t-1}
beta_pred <- 0 # this is beta_{t|t-1}
P11 <- 0 # this is P_{t|t-1}
P21 <- 0 
P12 <- 0 
P22 <- 0 
v <- 0 # this will be the innovation error
Kmu <- 0 # the Kalman gain
Kbeta <- 0 
F <- 0 # the conditional variance of v_t
llk<- 0 # the log-likelihood value


# initialise the recursion

# mu_pred[1] is set equal to its unconditional mean
# P is the conditional variance of mu_t, let us initialise it as the signal-to-noise ratio

mu_pred[1] <- 0
beta_pred[1] <- 0
P11[1] <- (sigma_eta^2)/(sigma_e^2)
P22[1] <- (sigma_xi^2)/(sigma_e^2)
P12[1] <- 0
P21[1] <- 0
llk[1] <- 0


# the recursion
for(t in 1:(n-1)){
  v[t] = y[t] - mu_pred[t];
  F[t] = P11[t] + sigma_e^2;
  Kmu[t] = (P11[t]+P21[t])/F[t]
  Kbeta[t] = P21[t]/F[t]
  P11[t+1] = P11[t] + P12[t] + P21[t] + P22[t] + sigma_eta^2 - (P11[t] + P21[t])^2/F[t]
  P12[t+1] = P12[t] + P22[t] - (P11[t]*P21[t]+ P21[t]^2)/F[t]
  P21[t+1] = P21[t] + P22[t] - (P11[t]*P21[t]+ P21[t]^2)/F[t]
  P22[t+1] = P22[t] + sigma_xi^2 - (P21[t]^2)/F[t]
  mu_pred[t+1] = mu_pred[t] + beta_pred[t] + Kmu[t]*v[t]
  beta_pred[t+1] = beta_pred[t] + Kbeta[t]*v[t]
  llk[t] = - 0.5 * log(F[t] + (v[t]^2/F[t]))
}

llk = sum(llk)

ts.plot(y)
lines(mu, col = 'red')
lines(mu_pred, col = 'green')

ts.plot(y)
lines(mu, col = "red")

ts.plot(y)
lines(mu_pred, col = "green")

ts.plot(mu)
lines(mu_pred, col = "blue")

# in this case, we have llk = -276.66, which is worse than the last case



################################################################################
# Question No. 2
################################################################################
# LLT model with Gaussian and Student-t measurement noise
# the case with sigma_e = 0.1, sigma_eta = 1, and sigma_xi = 1
################################################################################

# fix the length of the series 
n = 300

# simulate the noise of the signal (in the transition equation (SSF))
sigma_eta = 1
eta = sigma_eta * rnorm(n)

# simulate the measurement noise (in the observation equation)
sigma_e = .1
e = sigma_e * rnorm(n)
# q, signal-to-noise ratio, q = sigma^2_eta/sigma^2_e

q = sigma_eta^2/sigma_e^2

# case when the data are generated by a heavy tailed distribution  
# note that assuming s Student t distribution for the measurement error does break
# the assumption of KF which is optimal (in MMSE sense) for linear Gaussian models 
# nu is the degrees of freedom parameter (nu -> infinity then t -> Gaussian)

#nu = 200
#e = sqrt((nu-2)/nu) * sigma_e * rt(n,nu)
#ts.plot(e)

sigma_xi = 1
xi = sigma_xi * rnorm(n)

# plot of eta (error affecting the signal) versus epsilon (e, measurement error) and xi

ts.plot(e)
lines(eta, col = "red")
lines(xi, col = "blue")


mu = 0
mu[1] <- 0   
beta = 0
beta[1] <- 0   
y = 0 
for(t in 1:(n-1)){
  beta[t+1] = beta[t] + xi[t]
  mu[t+1] = beta[t] + mu[t] + eta[t]
  y[t+1] =  mu[t+1] + e[t+1]
}

ts.plot(mu)

# y = mu + e

# Plot of y versus mu
ts.plot(y)
lines(mu,col="red")

# briefly check the correlation structure of the simulated process 
# to be sure that we are properly simulating an AR(1) signal plus noise process

acf(y, lag.max = 80, drop.lag.0 = FALSE)
pacf(y, lag.max = 80)
acf(mu, lag.max = 80, drop.lag.0 = FALSE)
pacf(mu, lag.max = 80)

# as expected, the acf and pacf of y and mu suggest non-stationary process


# Differences
dy = diff(y)
ts.plot(dy)
acf(dy,lag.max = 100, drop.lag.0 = FALSE)
pacf(dy,lag.max = 100, drop.lag.0 = FALSE)



dmu = diff(mu)
ts.plot(dmu)
lines(eta, col = "red")
acf(dmu,lag.max = 100, drop.lag.0 = FALSE)
pacf(dmu,lag.max = 100, drop.lag.0 = FALSE)



dy = diff(diff(y))
ts.plot(dy)
acf(dy,lag.max = 100, drop.lag.0 = FALSE)
pacf(dy,lag.max = 100, drop.lag.0 = FALSE)

#  (1-L^2)y_t and (1-L^2)mu_t are both ARMA (1,1)

dmu = diff(diff(mu))
ts.plot(dmu)
lines(eta, col = "red")
acf(dmu,lag.max = 100, drop.lag.0 = FALSE)
pacf(dmu,lag.max = 100, drop.lag.0 = FALSE)

# hence, our LLT process is an ARIMA(1,2,1) process

################################################################################
# KALMAN filter recursions (see Durbin and Koopman, 2001) 
################################################################################

# we have simulated y_t and an underlying signal mu_t  
# (observed and latent or unobserved in real data, respectively).
# Now we wish to estimate mu_{t|t-1} - note that we are 
# assuming that phi, sigma_e, sigma_eta, i.e. the static parameters, are known.
# In practice, they are estimated by maximising the Gaussian likelihood 
# (we shall see it) 

# initial conditions needed 
# mu_10, P_10 that we can set at the unconditional mean and variance of mu 

# KF recursion

mu_pred <- 0 # this is mu_{t|t-1}
beta_pred <- 0 # this is beta_{t|t-1}
P11 <- 0 # this is P_{t|t-1}
P21 <- 0 
P12 <- 0 
P22 <- 0 
v <- 0 # this will be the innovation error
Kmu <- 0 # the Kalman gain
Kbeta <- 0 
F <- 0 # the conditional variance of v_t
llk<- 0 # the log-likelihood value


# initialise the recursion

# mu_pred[1] is set equal to its unconditional mean
# P is the conditional variance of mu_t, let us initialise it as the signal-to-noise ratio

mu_pred[1] <- 0
beta_pred[1] <- 0
P11[1] <- (sigma_eta^2)/(sigma_e^2)
P22[1] <- (sigma_xi^2)/(sigma_e^2)
P12[1] <- 0
P21[1] <- 0
llk[1] <- 0


# the recursion
for(t in 1:(n-1)){
  v[t] = y[t] - mu_pred[t];
  F[t] = P11[t] + sigma_e^2;
  Kmu[t] = (P11[t]+P21[t])/F[t]
  Kbeta[t] = P21[t]/F[t]
  P11[t+1] = P11[t] + P12[t] + P21[t] + P22[t] + sigma_eta^2 - (P11[t] + P21[t])^2/F[t]
  P12[t+1] = P12[t] + P22[t] - (P11[t]*P21[t]+ P21[t]^2)/F[t]
  P21[t+1] = P21[t] + P22[t] - (P11[t]*P21[t]+ P21[t]^2)/F[t]
  P22[t+1] = P22[t] + sigma_xi^2 - (P21[t]^2)/F[t]
  mu_pred[t+1] = mu_pred[t] + beta_pred[t] + Kmu[t]*v[t]
  beta_pred[t+1] = beta_pred[t] + Kbeta[t]*v[t]
  llk[t] = - 0.5 * log(F[t] + (v[t]^2/F[t]))
}

llk = sum(llk)

ts.plot(y)
lines(mu, col = 'red')
lines(mu_pred, col = 'green')

ts.plot(y)
lines(mu, col = "red")

ts.plot(y)
lines(mu_pred, col = "green")

ts.plot(mu)
lines(mu_pred, col = "blue")

# in this case, we have llk = -188.44, we have massive improvement over the last 2 cases, but still the estimation is
# worse than the LLM (or AR(1) + noise) models. The improvement suggests that lower values of sigma_e leads
# to better estimation




################################################################################
# Question No. 2
################################################################################
# LLT model with Gaussian and Student-t measurement noise
# the case with sigma_e = 1, sigma_eta = 0.1, and sigma_xi = 0.1
################################################################################

# fix the length of the series 
n = 300

# simulate the noise of the signal (in the transition equation (SSF))
sigma_eta = .1
eta = sigma_eta * rnorm(n)

# simulate the measurement noise (in the observation equation)
sigma_e = 1
e = sigma_e * rnorm(n)
# q, signal-to-noise ratio, q = sigma^2_eta/sigma^2_e

q = sigma_eta^2/sigma_e^2

# case when the data are generated by a heavy tailed distribution  
# note that assuming s Student t distribution for the measurement error does break
# the assumption of KF which is optimal (in MMSE sense) for linear Gaussian models 
# nu is the degrees of freedom parameter (nu -> infinity then t -> Gaussian)

#nu = 200
#e = sqrt((nu-2)/nu) * sigma_e * rt(n,nu)
#ts.plot(e)

sigma_xi = .1
xi = sigma_xi * rnorm(n)

# plot of eta (error affecting the signal) versus epsilon (e, measurement error) and xi

ts.plot(e)
lines(eta, col = "red")
lines(xi, col = "blue")


mu = 0
mu[1] <- 0   
beta = 0
beta[1] <- 0   
y = 0 
for(t in 1:(n-1)){
  beta[t+1] = beta[t] + xi[t]
  mu[t+1] = beta[t] + mu[t] + eta[t]
  y[t+1] =  mu[t+1] + e[t+1]
}

ts.plot(mu)

# y = mu + e

# Plot of y versus mu
ts.plot(y)
lines(mu,col="red")

# briefly check the correlation structure of the simulated process 
# to be sure that we are properly simulating an AR(1) signal plus noise process

acf(y, lag.max = 80, drop.lag.0 = FALSE)
pacf(y, lag.max = 80)
acf(mu, lag.max = 80, drop.lag.0 = FALSE)
pacf(mu, lag.max = 80)

# as expected, the acf and pacf of y and mu suggest non-stationary process


# Differences
dy = diff(y)
ts.plot(dy)
acf(dy,lag.max = 100, drop.lag.0 = FALSE)
pacf(dy,lag.max = 100, drop.lag.0 = FALSE)



dmu = diff(mu)
ts.plot(dmu)
lines(eta, col = "red")
acf(dmu,lag.max = 100, drop.lag.0 = FALSE)
pacf(dmu,lag.max = 100, drop.lag.0 = FALSE)



dy = diff(diff(y))
ts.plot(dy)
acf(dy,lag.max = 100, drop.lag.0 = FALSE)
pacf(dy,lag.max = 100, drop.lag.0 = FALSE)

#  (1-L^2)y_t is MA(1),  (1-L^2)mu_t is ARMA (1,1)

dmu = diff(diff(mu))
ts.plot(dmu)
lines(eta, col = "red")
acf(dmu,lag.max = 100, drop.lag.0 = FALSE)
pacf(dmu,lag.max = 100, drop.lag.0 = FALSE)

# hence, our LLT process is an ARIMA(0,2,1) process

################################################################################
# KALMAN filter recursions (see Durbin and Koopman, 2001) 
################################################################################

# we have simulated y_t and an underlying signal mu_t  
# (observed and latent or unobserved in real data, respectively).
# Now we wish to estimate mu_{t|t-1} - note that we are 
# assuming that phi, sigma_e, sigma_eta, i.e. the static parameters, are known.
# In practice, they are estimated by maximising the Gaussian likelihood 
# (we shall see it) 

# initial conditions needed 
# mu_10, P_10 that we can set at the unconditional mean and variance of mu 

# KF recursion

mu_pred <- 0 # this is mu_{t|t-1}
beta_pred <- 0 # this is beta_{t|t-1}
P11 <- 0 # this is P_{t|t-1}
P21 <- 0 
P12 <- 0 
P22 <- 0 
v <- 0 # this will be the innovation error
Kmu <- 0 # the Kalman gain
Kbeta <- 0 
F <- 0 # the conditional variance of v_t
llk<- 0 # the log-likelihood value


# initialise the recursion

# mu_pred[1] is set equal to its unconditional mean
# P is the conditional variance of mu_t, let us initialise it as the signal-to-noise ratio

mu_pred[1] <- 0
beta_pred[1] <- 0
P11[1] <- (sigma_eta^2)/(sigma_e^2)
P22[1] <- (sigma_xi^2)/(sigma_e^2)
P12[1] <- 0
P21[1] <- 0
llk[1] <- 0


# the recursion
for(t in 1:(n-1)){
  v[t] = y[t] - mu_pred[t];
  F[t] = P11[t] + sigma_e^2;
  Kmu[t] = (P11[t]+P21[t])/F[t]
  Kbeta[t] = P21[t]/F[t]
  P11[t+1] = P11[t] + P12[t] + P21[t] + P22[t] + sigma_eta^2 - (P11[t] + P21[t])^2/F[t]
  P12[t+1] = P12[t] + P22[t] - (P11[t]*P21[t]+ P21[t]^2)/F[t]
  P21[t+1] = P21[t] + P22[t] - (P11[t]*P21[t]+ P21[t]^2)/F[t]
  P22[t+1] = P22[t] + sigma_xi^2 - (P21[t]^2)/F[t]
  mu_pred[t+1] = mu_pred[t] + beta_pred[t] + Kmu[t]*v[t]
  beta_pred[t+1] = beta_pred[t] + Kbeta[t]*v[t]
  llk[t] = - 0.5 * log(F[t] + (v[t]^2/F[t]))
}

llk = sum(llk)

ts.plot(y)
lines(mu, col = 'red')
lines(mu_pred, col = 'green')

ts.plot(y)
lines(mu, col = "red")

ts.plot(y)
lines(mu_pred, col = "green")

ts.plot(mu)
lines(mu_pred, col = "blue")

# in this case, we have llk = -128.32, we have massive improvement over the last 3 cases, like the LLM model with
# sigma_e = 1 and sigma_eta = 0.1 where llk substationally increased compared to the LLM with
# sigma_e = sigma_eta = 1. The improvement suggests that lower values of sigma_eta and sigma_xi makes the estimation
# more accurate and, eta and xi have more impact on estimation than epsilon