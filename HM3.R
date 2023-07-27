# install.packages("TSA")
library(TSA)
data = Nile

################################################################################
# Question 1
################################################################################

KF <-function(y, phi, sigma_e, sigma_eta){
  
# we here separate the arguments of the function that will be estimated - phi, sigma_e,
# sigma_eta by the arguments of the function that will not - mu_1|0 and P_1|0 
  
m10 = 0
P10 = 1
  
n = NROW(y)
  
# allocate space 
mu_pred = array(data = NA, dim = c(n))  # this is mu_{t|t-1}
P = 0                                   # this is P_{t|t-1} 
v = array(data = NA, dim = c(n))        # this  will be the innovation error 
K = 0                                   # the Kalman gain 
F = 0                                   # the conditional variance of v_t 
dllk = array(data = NA, dim = c(n))    #the log-likelihood value 
llk = 0 
  
# initialise the recursion 
  
mu_pred[1] = m10
P[1] = P10
  
# the recursion 
  
for(t in 1:(n-1)){
    v[t] = y[t] - mu_pred[t]; 
    F[t] = P[t] + sigma_e^2;
    K[t] = (phi * P[t])/F[t]
    P[t+1] = phi^2 * P[t] + sigma_eta^2 - K[t]*F[t]*K[t]
    mu_pred[t+1] = phi * mu_pred[t] + K[t]*v[t]
    dllk[t] = - 0.5 * log(F[t] + (v[t]^2/F[t]))
    llk  = llk + dllk[t]
}
  
#llk = sum(dllk)
  
out <- list(mu_pred = as.ts(mu_pred), llk = llk)
#out <- list(v, mu_pred)
  
return(out)
}



# q = sigma_eta^2/sigma_e^2 implies that sigma_eta ^ 2 = q * sigma_e ^ 2
# in the chapter, sigma_e ^ 2 = 15099 were used
likelihood = 0
sigma_e = sqrt(15099)

########################### q = 1
q = 1
sigma_eta = sqrt(q) * sigma_e

first = KF(data, 1, sigma_e, sigma_eta)
mu_pred1 <- ts(first$mu_pred, start=1871, frequency=1) # changes the type of the column to "Time-Series" (to be able to plot together with y)
likelihood[1] = first$llk


ts.plot(data)
lines(mu_pred1, col = 'green')

# llk = -523.2093

########################### q = 0.0360
q = 0.0360
sigma_eta = sqrt(q) * sigma_e

second = KF(data, 1, sigma_e, sigma_eta)
mu_pred2 <- ts(second$mu_pred, start=1871, frequency=1)
likelihood[2] = second$llk


ts.plot(data)
lines(mu_pred2, col = 'green')

# llk = -485.3073
# based on llk, we see improvement over the last case

########################### q = 0.0745
q = 0.0745
sigma_eta = sqrt(q) * sigma_e

third = KF(data, 1, sigma_e, sigma_eta)
mu_pred3 <- ts(third$mu_pred, start=1871, frequency=1)
likelihood[3] = third$llk


ts.plot(data)
lines(mu_pred3, col = 'green')

# llk = -489.3709
# based on llk, the estimation got slightly worse than the previous case

########################### q = 0.0974
q = 0.0974
sigma_eta = sqrt(q) * sigma_e

fourth = KF(data, 1, sigma_e, sigma_eta)
mu_pred4 <- ts(fourth$mu_pred, start=1871, frequency=1)
likelihood[4] = fourth$llk


ts.plot(data)
lines(mu_pred4, col = 'green')

# llk = -491.2756
# based on llk, the estimation got slightly worse than the previous case


########################### q = 0.0973
q = 0.0973
sigma_eta = sqrt(q) * sigma_e

fifth = KF(data, 1, sigma_e, sigma_eta)
mu_pred5 <- ts(fifth$mu_pred, start=1871, frequency=1)
likelihood[5] = fifth$llk


ts.plot(data)
lines(mu_pred5, col = 'green')

# llk = -491.2678
# the llk is very similar to the last case

########################### Results

likelihood

ts.plot(data)
lines(mu_pred1, col = 'green')
lines(mu_pred2, col = 'red')
lines(mu_pred3, col = 'blue')
lines(mu_pred4, col = 'yellow')
lines(mu_pred5, col = 'cyan')

# we can see that mu_pred4 is not noticeable as it is covered by mu_pred5 (they completely fall on each other)

# The best estimation is for mu_pred2 (red), and the worst is for mu_pred1 (green)


################################################################################
# Question 2
################################################################################

n = 100
# y = data 
# or y = Nile
# here we used the "data" we used in question 1 (which was Nile dataset)

q = 0.0973
sigma_e = 15099
ep <- rnorm(n)*sigma_e

sigma_eta = sqrt(q) * sqrt(sigma_e)
etap <- rnorm(n)*sigma_eta

#allocate space
yp <- rep(0,n)
alphap <- rep(0,n)
alphap[1] <- 0 # initial value
for (t in 1:(n-1)){
  alphap[t+1] = alphap[t] + etap[t]
  yp[t+1] = alphap[t+1] + ep[t+1]
}

# recursion for epsilon plus hat

# allocate space
alpha_pred <- rep(0,n) # this is mu_{t|t-1}
P <- rep(0,n) # this is P_{t|t-1}
v <- rep(0,n) # this will be the innovation error
K <- rep(0,n) # the Kalman gain
F <- rep(0,n) # the conditional variance of v_t
L <- rep(0,n)
u <- rep(0,n)
r <- rep(0,n)
ep_hat <- rep(0,n) # epsilon plus hat

# initialise the recursion
# mu_pred[1] is set equal to its unconditional mean
# P is the conditional variance of mu_t, let us initialise it as the signal-to-noise ratio

alpha_pred[1] = 0
P[1] = (sigma_eta^2)/(sigma_e^2)
r[n] = 0

for(t in 1:(n-1)){
  v[t] = yp[t] - alpha_pred[t]
  F[t] = P[t] + sigma_e^2;
  K[t] = P[t]/F[t]
  L[t] = 1 - K[t] # equation 2.30
  P[t+1] = P[t] + sigma_eta^2 - K[t]*F[t]*K[t]
  alpha_pred[t+1] = alpha_pred[t] + K[t]*v[t]
}

# to obtain the 100th (nth) element of v, F, K, and L
t = n
v[t] = yp[t] - alpha_pred[t]
F[t] = P[t] + sigma_e^2;
K[t] = P[t]/F[t]
L[t] = 1 - K[t]


for(t in n:2){
  r[t-1] = v[t]/F[t] + L[t]*r[t] # equation 2.36
}
for(t in 1:n){
  u[t] = v[t]/F[t] - K[t]*r[t] # equation 2.45
  ep_hat[t] = u[t] * sigma_e^2 # equation 2.44
}

# recursion for epsilon hat

# allocate space
alpha_pred <- rep(0,n) # this is mu_{t|t-1}
P <- rep(0,n) # this is P_{t|t-1}
v <- rep(0,n) # this will be the innovation error
K <- rep(0,n) # the Kalman gain
F <- rep(0,n) # the conditional variance of v_t
L <- rep(0,n)
u <- rep(0,n)
r <- rep(0,n)
e_hat <- rep(0,n) # epsilon hat

# initialise the recursion
# mu_pred[1] is set equal to its unconditional mean
# P is the conditional variance of mu_t, let us initialise it as the signal-to-noise ratio

alpha_pred[1] = 0
P[1] = (sigma_eta^2)/(sigma_e^2)
r[n] = 0

for(t in 1:(n-1)){
  v[t] = data[t] - alpha_pred[t]
  F[t] = P[t] + sigma_e^2;
  K[t] = P[t]/F[t]
  L[t] = 1 - K[t] # equation 2.30
  P[t+1] = P[t] + sigma_eta^2 - K[t]*F[t]*K[t]
  alpha_pred[t+1] = alpha_pred[t] + K[t]*v[t]
}

# to obtain the 100th (nth) element of v, F, K, and L
t = n
v[t] = data[t] - alpha_pred[t]
F[t] = P[t] + sigma_e^2;
K[t] = P[t]/F[t]
L[t] = 1 - K[t]


for(t in n:2){
  r[t-1] = v[t]/F[t] + L[t]*r[t] # equation 2.36
}
for(t in 1:n){
  u[t] = v[t]/F[t] - K[t]*r[t] # equation 2.45
  e_hat[t] = u[t] * sigma_e^2 # equation 2.44
}

# allocate space
etilde <- ep - ep_hat + e_hat
alpha_tilde <- rep(0,n)
eta_tilde <- rep(0,n)
eta_hat <- rep(0,n)

alpha_tilde[1]<- data[1]- etilde[1]


for(t in 2:n){
  alpha_tilde[t]<- data[t] - etilde[t]
  eta_tilde[t-1]<- alpha_tilde[t] - alpha_tilde[t-1]
}

for(t in 1:(n-1)){
  eta_hat[t]<- alpha_pred[t+1] - alpha_pred[t]
}

ts.plot(alpha_tilde)
lines(alpha_pred, col='red') # alpha hat
lines(alphap, col='blue')

ts.plot(etilde)
lines(e_hat, col='red')

ts.plot(eta_tilde)
lines(eta_hat, col='red')



