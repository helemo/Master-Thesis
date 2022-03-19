library(nimble)
library(tidyverse)
library(ggpubr)
library(reshape2)
library(DiceKriging)
library(bayesplot)
library(hexbin)
library(coda)
library(GGally)



time <-c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 
         0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 
         0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 
         0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 
         0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 
         0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 
         0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 
         0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 
         0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 
         0.99, 1)

Q <- c(1.35, 4.023, 9.972, 21.546, 41.022, 69.903, 108.666, 156.393, 
       210.825, 268.866, 327.438, 383.427, 433.791, 475.776, 507.357, 
       527.517, 536.355, 535.086, 525.663, 510.255, 490.806, 468.783, 
       445.221, 420.741, 395.667, 370.134, 344.16, 317.628, 290.286, 
       261.819, 232.029, 201.006, 169.254, 137.682, 107.568, 83.349, 
       53.757, 17.964, -22.032, -43.884, -41.922, -24.687, -7.974, 0.324,
       1.827, 0.963, 0.369, 0.954, 1.017, 0.855, 0, 0, 0, 0, 0, 0, 0, 
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
       0, 0)


# Noise in Q
#Q <- Q+rinvgamma(length(Q),11,50)
#plot(time,Q)
# R = 0.9
# C = 1.3

WK_simulate.R = function(flow01, flow_min, flow_max, time, R, C){
  #Q = flow
  Q = flow01*(flow_max-flow_min)+flow_min
  
  impedance_WK2 = function(w,R,C){
    Z=R/(1.0+1i*w*R*C)
    return(Z)
  }
  d = time[2]-time[1]
  n = length(time)
  if((n %% 2) == 0) {
    f = c(seq(0,n/2-1, by = 1), seq(-n/2, -1,  by=1))/(d*n)
  } else {
    f = c(seq(0,(n-1)/2, by = 1), seq(-(n-1)/2, -1,  by=1))/(d*n)
  }
  w = 2*pi*f
  Qfft = fft(Q) #+rnorm(length(Q),0,3)
  
  Zmod = impedance_WK2(w,R,C)
  P = Re(fft(Zmod*Qfft, inverse = TRUE))/n
  return(P)
}

# plot(time,P)

WK_nimble <- nimbleRcall(function(flow01 = double(1),flow_min=double(0), flow_max=double(0), time = double(1), R=double(0), C = double(0)){}
                         , Rfun = 'WK_simulate.R', returnType = double(1))

WK3_simulate = function(flow01,flow_max, flow_min, time, R, C, Z){
  Q = flow01*(flow_max-flow_min)+flow_min
  #Q = flow
  Z_ao = Z
  impedance_WK3 = function(w,R,C,Z_ao){
    Z=(Z_ao + R + 1i*w*R*C*Z_ao)/(1.0+1i*w*R*C)
    return(Z)
  }
  d = time[2]-time[1]
  n = length(time)
  if((n %% 2) == 0) {
    f = c(seq(0,n/2-1, by = 1), seq(-n/2, -1,  by=1))/(d*n)
  } else {
    f = c(seq(0,(n-1)/2, by = 1), seq(-(n-1)/2, -1,  by=1))/(d*n)
  }
  w = 2*pi*f
  Qfft = fft(Q)
  
  Zmod = impedance_WK3(w,R,C, Z_ao)
  P = Re(fft(Zmod*Qfft, inverse = TRUE))/n
  return(P)
}

cov_sq_exp <- nimbleFunction(     
  run = function(X = double(2), rho = double(1), sigma = double(0), delta = double(0)) {
    returnType(double(2))
    n <- nimDim(X)[1]
    p <- nimDim(X)[2]
    K <- matrix(nrow = n, ncol = n, init = FALSE)
    sigma2 <- sigma*sigma
    rho2 <- rho^2
    for(i in 1:(n-1)){
      K[i,i] <- sigma2
      for(j in (i+1):n){ 
        d <- sum((X[i,1:p]-X[j,1:p])^2/rho2[1:p])
        K[i, j] <- sigma2*exp(-d/2)
        K[j, i] <- K[i,j]
      }
      K[n,n] <- sigma2 
    }
    K = K + diag(rep(delta^2, n))
    return(K)
  })

ccov_sq_exp <- compileNimble(cov_sq_exp)




code <- nimbleCode({
  # Fix the priors
  R ~ dunif(0.5,3)
  C ~ dunif(0.5,3)
  
  # Full calib:
  # sigma ~ dunif(0,10)
  # sigma_B ~ dunif(0,500)
  # rho[1] ~ dunif(0,1) # #dunif(0,1000)
  # rho[2] ~  dunif(0,1) #dunif(0,1) #dunif(0,5000)
  # 
  # Mod 1
  sigma ~  dunif(0,10) #dinvgamma(8102,24303) # l = 2.9 u = 3.1  #dinvgamma(11,50) #l = 0 u = 10 #dunif(0,10)#dinvgamma(112.25,389.375) #l = 2.5 u = 4.5 dinvgamma(8102,24303) # #dinvgamma(8102,24303) # l = 2.9 u = 3.1 #dunif(0,6) #
  sigma_B ~ dunif(0,500)#dinvgamma(11,250) # l = 0 u = 50 dunif(0,500)#dinvgamma(11,250) # l = 0 u = 50  #dunif(0,80) #dinvgamma(11,250) # l = 0 u = 50 #dunif(0,50) #dinvgamma(146,580)
  # fix lengthscale MAP from bias est
  rho[1] <- 0.5 #0.63952674  #0.45   #0.4820121 ##0.6317292    #0.6327226
  rho[2] <- 0.07 #0.05086166  #0.62#0.6183833 # #0.6804413    #0.6814067

  # Mod 2 
  # sigma ~ dinvgamma(8102,24303) # l = 2.9 u = 3.1  dunif(0,10)# iid   #dinvgamma(11,50) #l = 0 u = 10
  # sigma_B ~ dinvgamma(11,500) #l 0 u = 100 dunif(0,500)# iid    # dinvgamma(11,250) #dinvgamma(11,2500) #l = 0 u =  500
  # # lengthscale post from bias est as prior
  # rho[1] ~  dinvgamma(11,5) #l = 0 u = 1  dinvgamma(83,36.9) # l = 0.3 u = 0.6
  # rho[2] ~  dinvgamma(11,1/2) #l = 0 u = 0.1 dinvgamma(62.84, 40.196) # l = 0.4 u = 0.9

  # adding transformed flow to covariance matrix
  X[1:N,1:p] <- matrix(c(flow01[1:N],time[1:N]),ncol=2)
  
  # NB if you want to run it with non transformed values, make sure to change flow_t and X_t
  # back to regular flow and X in P_sim and cov below.
  
  P_sim[1:N] <- WK_nimble(flow01=flow01[1:N], flow_min=Qmin, flow_max=Qmax, time = time[1:N], R=R, C=C)
  cov[1:N,1:N] <- cov_sq_exp(X = X[1:N,1:p], rho[1:p], sigma_B, sigma)
  
  # likelihood
  y[1:N] ~ dmnorm(P_sim[1:N], cov = cov[1:N,1:N])
  
})

trans_01 <- function(Q, flow){
  Qmax = max(Q); Qmin = min(Q)
  flow01 = (flow-Qmin)/(Qmax-Qmin)
  return(flow01)
}


# Transform flow Q to 01 scale
Qmax = max(Q); Qmin = min(Q)
ind = round(seq(1,101,length.out = 50))
flow = Q[ind]; time = time[ind]
flow01 = trans_01(Q,flow)

lP = 50
nrep = 3
#fixed noise for Z experimen5
fixed_noise = rnorm(nrep*lP,0,3)

# Set true values of parameters
Rtrue_WK3=0.9;Ctrue_WK3=1.3

df <- data.frame(samp = double())#, Z = double())
Z_list <- c(0.01,0.05,0.1,0.15,0.2)

name <- "Z-m1"


for(i in 1:length(Z_list)){
  
  print("sim nr")
  print(i)
  
  sd_noise <- 3
  
  
  Ztrue_WK3=Z_list[i] 
  Rtrue=(0.9+Ztrue_WK3);Ctrue=1.3
  
  P = WK3_simulate(flow=flow01,flow_max = Qmax, flow_min = Qmin, time=time, R=Rtrue_WK3, C=Ctrue_WK3, Z = Ztrue_WK3)
  #plot(time,P)
  
  set.seed(321)
 
  
  ####################################################
  # iid noise
  P_obs = rep(P,nrep) +  fixed_noise
  ####################################################
  
  N = length(P_obs)
  X = matrix(c(rep(flow01,nrep),rep(time,nrep)),ncol=2)
  constants <- list(N = N, p=2, time = rep(time,nrep), flow01 = rep(flow01,nrep), X = X, Qmax = Qmax, Qmin = Qmin)
  
  data <- list(y = P_obs)
  inits <- list(R = 0.8, C = 0.8, sigma = 2, rho  = c(0.5, 0.5), sigma_B=10)
  inits$cov <- ccov_sq_exp(X, inits$rho, inits$sigma_B, inits$sigma)
  
  model <- nimbleModel(code, constants = constants, data = data, inits = inits)
  cModel <- compileNimble(model)
  
  conf <- configureMCMC(model)
  conf$printSamplers()
  MCMC <- buildMCMC(conf)
  cMCMC <- compileNimble(MCMC, project = cModel)
  # run MCMC chain
  samples <- runMCMC(cMCMC, niter = 10000, nburnin = 5000)
  samp = mcmc(samples)
  print(effectiveSize(samp))
  df <- rbind(df,data.frame(samp))


  
  
  
}

samp

df_l <- df %>% select(colnames(df)) %>% gather(key="parameter", value="value")
df_l$Z <- rep(c(rep("0.01",length(samp[,"C"])),rep("0.05",length(samp[,"C"])), 
                rep("0.1",length(samp[,"C"])),rep("0.15",length(samp[,"C"])), rep("0.2",length(samp[,"C"]))),ncol(df))

p <- ggplot(df_l,aes(value, fill = Z)) + geom_histogram(aes( y= ..density..),bins = 25)
pl_hist = p + facet_wrap(~parameter, scales = "free")+
  theme_bw()
ggsave(stringr::str_c(name,"-post.png"))

pl_hist



#effectiveSize(samp)
#par(mfrow = c(ncol(samples)/2,2))
# for(i in 1:ncol(samples)){
#   plot(samples[,i], type="l", xlab = "iterations", ylab = colnames(samples)[i])
# }
