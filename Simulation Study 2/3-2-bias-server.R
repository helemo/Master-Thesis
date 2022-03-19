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


# Define WK2 model simulator in NIMBLE
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

# Define WK3 model simulator in NIMBLE
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

# Covariance function
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

# Define NIMBLE model
code <- nimbleCode({
  # Fix the priors
  R ~ dunif(0.5,3)
  C ~ dunif(0.5,3)
  
  # Full calib:
  sigma ~ dunif(0,10)
  sigma_B ~ dunif(0,500)
  rho[1] ~ dunif(0,1) # #dunif(0,1000)
  rho[2] ~  dunif(0,1) #dunif(0,1) #dunif(0,5000)

  # Mod 1
  # sigma ~  dunif(0,10) #dinvgamma(8102,24303) # l = 2.9 u = 3.1  #dinvgamma(11,50) #l = 0 u = 10 #dunif(0,10)#dinvgamma(112.25,389.375) #l = 2.5 u = 4.5 dinvgamma(8102,24303) # #dinvgamma(8102,24303) # l = 2.9 u = 3.1 #dunif(0,6) #
  # sigma_B ~ dunif(0,500)#dinvgamma(11,250) # l = 0 u = 50 dunif(0,500)#dinvgamma(11,250) # l = 0 u = 50  #dunif(0,80) #dinvgamma(11,250) # l = 0 u = 50 #dunif(0,50) #dinvgamma(146,580)
  # # fix lengthscale MAP from bias est
  # rho[1] <- 0.5 #0.63952674  #0.45   #0.4820121 ##0.6317292    #0.6327226
  # rho[2] <- 0.07 #0.05086166  #0.62#0.6183833 # #0.6804413    #0.6814067

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



# Set true values of parameters
Rtrue_WK3=0.9;Ctrue_WK3=1.3
Ztrue_WK3=0.1
Rtrue=(0.9+Ztrue_WK3);Ctrue=1.3

P = WK3_simulate(flow=flow01,flow_max = Qmax, flow_min = Qmin, time=time, R=Rtrue_WK3, C=Ctrue_WK3, Z = Ztrue_WK3)
#plot(time,P)

set.seed(321)
lP = length(P)
nrep = 3


# data frame to store parameter estimates 
df_params <- data.frame(C_lower = double(),
                        C = double(),
                        C_upper = double(),
                        
                        R_lower = double(),
                        R = double(),
                        R_upper = double(),
                        
                        rho1_lower = double(),
                        rho1 = double(),
                        rho1_upper = double(),
                        
                        rho2_lower = double(),
                        rho2 = double(),
                        rho2_upper = double(),
                        
                        sigma_lower = double(),
                        sigma = double(),
                        sigma_upper = double(),
                        
                        sigma_bias_lower = double(),
                        sigma_bias = double(),
                        sigma_bias_upper = double()
                        
)

# vector to store coverage for multiple runs
coverage <- c(0,0,0,0,0,0)

# vector to store noise and data
df_noise <- data.frame(noise = double(),
                       iteration = integer())
df_all_samp <- data.frame(samp = double(),
                          iteration = integer())

# Set nr of sample runs
N_sims = 1

#See name for simulation
name <- "NoName"


for(sim in 1:N_sims){
  
  print("sim nr")
  print(sim)
  
  sd_noise <- 3
  
  ####################################################
  # dependent noise 
  #NB modified for bias estimates
#   theta <- c(0.1,0.1)
#   print(theta)
#   sigma <- sd_noise
# 
#   eps_model <- km(design = data.frame(flow = flow01, time = time), response = rep(1, length(time)), covtype = "matern5_2",
#                   coef.cov = theta, coef.var = sigma^2, nugget = 2)
#   eps <- simulate(eps_model, nsim = 1, newdata = data.frame(flow = flow01, time = time))
# 
#   #plot(flow01,eps)
#   plot(time,eps)
# #
#   df_noise <- rbind(df_noise,data.frame(noise = (eps[1,]),
#                     iteration = rep(sim,length(eps[1,]))))
# 
#   write.csv(df_noise,stringr::str_c(name,"-eps.csv"))
# 
#   P_obs = rep(P,nrep) + rep(eps[1,],nrep)
#   plot(rep(time,nrep), P_obs)


  ####################################################
  # iid noise
  P_obs = rep(P,nrep) +  rnorm(nrep*lP,0,sd_noise)
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
  # obtain MAP estimate for parameters
  mean <-  apply(samp,2,mean)
  mean
  
  # store 90% credible interval
  CI <- apply(samp,2,sort)

  df_params <- rbind(df_params, data.frame( C_lower = CI[dim(samp)[1]*0.05,'C'], C = mean['C'], C_upper = CI[dim(samp)[1]*0.95,'C'],
                                            R_lower = CI[dim(samp)[1]*0.05,'R'], R = mean['R'], R_upper = CI[dim(samp)[1]*0.95,'R'],
# 
                                            # rho1_lower = CI[dim(samp)[1]*0.05,'rho[1]'], rho1 =  mean['rho[1]'],rho1_upper = CI[dim(samp)[1]*0.95,'rho[1]'],
                                            # rho2_lower = CI[dim(samp)[1]*0.05,'rho[2]'],rho2 =  mean['rho[2]'],rho2_upper = CI[dim(samp)[1]*0.95,'rho[2]'],

                                            sigma_lower = CI[dim(samp)[1]*0.05,'sigma'], sigma = mean['sigma'], sigma_upper = CI[dim(samp)[1]*0.95,'sigma'],
                                            sigma_bias_lower = CI[dim(samp)[1]*0.05,'sigma_B'], sigma_bias =  mean['sigma_B'], sigma_bias_upper = CI[dim(samp)[1]*0.95,'sigma_B']))


  # # Check coverage (mean in CI?)
  if (between(Ctrue,df_params$C_lower[sim],df_params$C_upper[sim])) coverage[1] <- coverage[1]+1
  if (between(Rtrue,df_params$R_lower[sim],df_params$R_upper[sim])) coverage[2] <- coverage[2]+1
  # if (between(mean['rho[1]'],df_params$rho1_lower[sim],df_params$rho1_upper[sim])) coverage[3] <- coverage[3]+1
  # if (between(mean['rho[2]'],df_params$rho2_lower[sim],df_params$rho2_upper[sim])) coverage[4] <- coverage[4]+1

  if (between(sd_noise,df_params$sigma_lower[sim],df_params$sigma_upper[sim])) coverage[5] <- coverage[5]+1
  if (between(mean['sigma_B'],df_params$sigma_bias_lower[sim],df_params$sigma_bias_upper[sim])) coverage[6] <- coverage[6]+1

  df_all_samp <- rbind(df_all_samp,data.frame(samp = samp, iteration = rep(sim,length(samp))))
  
  


}




##################
# Convergence check

# 
# inits2 <- list(R = 1, C = 1,sigma =5, rho  = c(0.8, 0.8), sigma_B=50)
# inits2$cov <- ccov_sq_exp(X, inits2$rho, inits2$sigma_B, inits2$sigma)
# model2 <- nimbleModel(code, constants = constants, data = data, inits = inits2)
# cModel2 <- compileNimble(model2)
# conf2 <- configureMCMC(model2)
# conf2$printSamplers()
# MCMC2 <- buildMCMC(conf2)
# cMCMC2 <- compileNimble(MCMC2, project = cModel2)
# samples2 <- runMCMC(cMCMC2, niter = 100000, nburnin = 50000)
# samp2 = mcmc(samples2)
# 
# inits3 <- list(R = 0.5, C = 0.5, sigma =1, rho  = c(0.3, 0.9), sigma_B=25)
# inits3$cov <- ccov_sq_exp(X, inits3$rho, inits3$sigma_B, inits3$sigma)
# model3 <- nimbleModel(code, constants = constants, data = data, inits = inits3)
# cModel3 <- compileNimble(model3)
# conf3 <- configureMCMC(model3)
# conf3$printSamplers()
# MCMC3 <- buildMCMC(conf3)
# cMCMC3 <- compileNimble(MCMC3, project = cModel3)
# samples3 <- runMCMC(cMCMC3, niter = 100000, nburnin = 50000)
# samp3 = mcmc(samples3)
# 
# inits4 <- list(R = 2, C = 2, sigma =8, rho  = c(0.7, 0.1), sigma_B=25)
# inits4$cov <- ccov_sq_exp(X, inits4$rho, inits4$sigma_B, inits4$sigma)
# model4 <- nimbleModel(code, constants = constants, data = data, inits = inits4)
# cModel4 <- compileNimble(model4)
# conf4 <- configureMCMC(model4)
# conf4$printSamplers()
# MCMC4 <- buildMCMC(conf4)
# cMCMC4 <- compileNimble(MCMC4, project = cModel4)
# samples4 <- runMCMC(cMCMC4, niter = 100000, nburnin = 50000)
# samp4 = mcmc(samples4)
# 
# par(mfrow = c(ncol(samples)/2,2))
# par(mfrow = c(ncol(samples)/2,2))
# for(i in 1:ncol(samples)){
#   plot(samples[,i], type="l", xlab = "iterations", ylab = colnames(samples)[i], col = 1)
#   lines(samples2[,i], type="l", xlab = "iterations", ylab = colnames(samples2)[i], col = 2)
#   lines(samples3[,i], type="l", xlab = "iterations", ylab = colnames(samples3)[i], col = 3)
#   lines(samples4[,i], type="l", xlab = "iterations", ylab = colnames(samples4)[i], col = 4)
#   legend("topleft", legend = 1:4, col=1:4, pch=1) # optional legend
#   
# }
# 
# df <- data.frame(samp)
# df$Chain = rep(1, nrow(df))
# df2 <- data.frame(samp2)
# df2$Chain = rep(2, nrow(df2))
# df3 <- data.frame(samp3)
# df3$Chain = rep(3, nrow(df3)) 
# df4 <- data.frame(samp4)
# df4$Chain = rep(4, nrow(df4))                 
# df <- rbind(df,df2,df3,df4)  
# mcmc_pairs( x = df, off_diag_args = list(size = .5), off_diag_fun="hex")
# warnings()


#################



# # ONE RUN CHECK
df <- data.frame(samp)
ggpairs(df,
        lower =  list(continuous = "density"))+theme_bw()
ggsave(stringr::str_c(name,"-corr.png"))

#write.csv(df_all_samp,stringr::str_c(name,"-data.csv"))

#effectiveSize(samp)
#par(mfrow = c(ncol(samples)/2,2))
# for(i in 1:ncol(samples)){
#   plot(samples[,i], type="l", xlab = "iterations", ylab = colnames(samples)[i])
# }

#apply(samp,2,mean)
#mcmc_pairs(samp, off_diag_args = list(size = .5), off_diag_fun="hex")


# ---------------------------
total_mean <- apply(df_params,2,mean)
#FIXED RHOS
params <- data.frame(parameter = c("C", "R", "sigma", "sigma_B"),
                     Z = c(Ctrue, Rtrue, sd_noise, total_mean['sigma_bias']),
                     E = c(total_mean['C'],
                           total_mean['R'],

                           total_mean['sigma'],
                           total_mean['sigma_bias']))
# #PRIOR RHO
# params <- data.frame(parameter = c("C", "R", "rho.1.", "rho.2.", "sigma", "sigma_B"),
#                      Z = c(Ctrue, Rtrue,total_mean['rho1'],total_mean['rho2'], sd_noise, total_mean['sigma_bias']),
#                      E = c(total_mean['C'],
#                            total_mean['R'],
#                            total_mean['rho1'],
#                            total_mean['rho2'],
#                            total_mean['sigma'],
#                            total_mean['sigma_bias']))
#----------------------------


# with mean in plot
df_l <- df %>% select(colnames(df)) %>% gather(key="parameter", value="value")
df_l$iteration = rep(1:nrow(df), ncol(df))
par_val <- params

p <- ggplot(df_l,aes(value)) + geom_histogram(aes( y= ..density..),bins = 25)
pl_hist = p + facet_wrap(~parameter, scales = "free")+
  geom_vline(data = par_val, aes(xintercept = Z), colour = "blue", size = 1.25)+
  geom_vline(data = par_val, aes(xintercept = E), colour = "red", size = 1.25)+
  theme_bw()

ps <- df_l %>% ggplot(aes(x=iteration, y = value)) + geom_line()
pl_trace = ps + facet_wrap(~parameter, scales = "free") +
  geom_hline(data = par_val, aes(yintercept = Z), colour = "blue", size = 1.25)+
  theme_bw()

ggarrange(pl_hist, pl_trace, ncol = 1)
ggsave(stringr::str_c(name,"-post.png"))

# ----------------------------
# # traceplot of CI for 100 simulations
# 
# # fixed rhos
# total_CI <- data.frame(parameter = c(rep("C",N_sims),rep("R",N_sims), rep("sigma", N_sims), rep("sigma_bias", N_sims)),
#                        lower = c(df_params$C_lower, df_params$R_lower, df_params$sigma_lower, df_params$sigma_bias_lower),
#                        upper = c(df_params$C_upper, df_params$R_upper, df_params$sigma_upper, df_params$sigma_bias_upper))
# 
# # # WITH RHO
# total_CI <- data.frame(parameter = c(rep("C",N_sims),rep("R",N_sims),rep("rho1", N_sims),rep("rho2", N_sims), rep("sigma", N_sims), rep("sigma_bias", N_sims)),
#                  lower = c(df_params$C_lower, df_params$R_lower, df_params$rho1_lower, df_params$rho2_lower, df_params$sigma_lower, df_params$sigma_bias_lower),
#                  upper = c(df_params$C_upper, df_params$R_upper, df_params$rho1_upper, df_params$rho2_upper, df_params$sigma_upper, df_params$sigma_bias_upper))
# 
# # -----------------------------
# total_CI$iteration <- rep(1:nrow(df_params), ncol(CI))
# 
# 
# C_CI <- ggplot(total_CI[1:N_sims,], aes(C))+
#   geom_segment(aes(x = lower, y = iteration, xend = upper, yend = iteration))+
#   geom_vline(data = params[1,], aes(xintercept = Z), colour = "blue")+
#   geom_vline(data = params[1,], aes(xintercept = E), colour = "red")+
#   theme_bw()
# ggsave(stringr::str_c(name,"-C_CI.png"))
# 
# R_CI <- ggplot(total_CI[(N_sims+1):(2*N_sims),], aes(R))+
#   geom_segment(aes(x = lower, y = iteration, xend = upper, yend = iteration))+
#   geom_vline(data = params[2,], aes(xintercept = Z), colour = "blue")+
#   geom_vline(data = params[2,], aes(xintercept = E), colour = "red")+
#   theme_bw()
# ggsave(stringr::str_c(name,"-R_CI.png"))
# 
# # ## FIXED RHOS -------
# total_CI
# # sigma_CI <- ggplot(total_CI[(2*N_sims+1):(3*N_sims),], aes(sigma))+
# #   geom_segment(aes(x = lower, y = iteration, xend = upper, yend = iteration))+
# #   geom_vline(data = params[3,], aes(xintercept = Z), colour = "blue")+
# #   geom_vline(data = params[3,], aes(xintercept = E), colour = "red")+
# #   theme_bw()
# # ggsave(stringr::str_c(name,"-sigma_CI.png"))
# # 
# # sigma_B_CI <- ggplot(total_CI[(3*N_sims+1):(4*N_sims),], aes(sigma_B))+
# #   geom_segment(aes(x = lower, y = iteration, xend = upper, yend = iteration))+
# #   geom_vline(data = params[4,], aes(xintercept = Z), colour = "blue")+
# #   geom_vline(data = params[4,], aes(xintercept = E), colour = "red")+
# #   theme_bw()
# # ggsave(stringr::str_c(name,"-sigma_B_CI.png"))
# 
# # ------------------
# 
# # # PRIORS RHOS ---------
# #
# rho1_CI <- ggplot(total_CI[(2*N_sims+1):(3*N_sims),], aes(rho1))+
#   geom_segment(aes(x = lower, y = iteration, xend = upper, yend = iteration))+
#   geom_vline(data = params[3,], aes(xintercept = Z), colour = "blue")+
#   geom_vline(data = params[3,], aes(xintercept = E), colour = "red")+
#   theme_bw()
# ggsave(stringr::str_c(name,"-rho1_CI.png"))
# 
# rho2_CI <- ggplot(total_CI[(3*N_sims+1):(4*N_sims),], aes(rho2))+
#   geom_segment(aes(x = lower, y = iteration, xend = upper, yend = iteration))+
#   geom_vline(data = params[4,], aes(xintercept = Z), colour = "blue")+
#   geom_vline(data = params[4,], aes(xintercept = E), colour = "red")+
#   theme_bw()
# ggsave(stringr::str_c(name,"-rho2_CI.png"))
# 
# sigma_CI <- ggplot(total_CI[(4*N_sims+1):(5*N_sims),], aes(sigma))+
#   geom_segment(aes(x = lower, y = iteration, xend = upper, yend = iteration))+
#   geom_vline(data = params[5,], aes(xintercept = Z), colour = "blue")+
#   geom_vline(data = params[5,], aes(xintercept = E), colour = "red")+
#   theme_bw()
# ggsave(stringr::str_c(name,"-sigma_CI.png"))
# 
# sigma_B_CI <- ggplot(total_CI[(5*N_sims+1):(6*N_sims),], aes(sigma_B))+
#   geom_segment(aes(x = lower, y = iteration, xend = upper, yend = iteration))+
#   geom_vline(data = params[6,], aes(xintercept = Z), colour = "blue")+
#   geom_vline(data = params[6,], aes(xintercept = E), colour = "red")+
#   theme_bw()
# ggsave(stringr::str_c(name,"-sigma_B_CI.png"))
# C_CI
# # # # --------------------------------------------------
# # 
# df_coverage <- data.frame(parameter = c("C", "R","rho1", "rho2", "sigma","sigma_B"),coverage = c(coverage[1],coverage[2],coverage[3],coverage[4],coverage[5],coverage[6]))
# 
# cover <- ggplot(df_coverage, aes(parameter, coverage))+
#   geom_point()
# ggsave(stringr::str_c(name,"-coverage.png"))
# cover
# ###################################
# 
# ################ PREDICT ######################
# 
# # 
# # # # # 
# sampl_mtrx <-  data.matrix(df)
# # smpl_size <- data_frame(parameter = c("C","R","sigma"), Z = as.double(coda::effectiveSize(samples)))
# # ggplot(smpl_size, aes(parameter,Z)) + geom_point()
# # ggsave(stringr::str_c(name,"-ess.png"))
# #
# dataNodes <- model$getNodeNames(dataOnly = TRUE)
# 
# parentNodes <- model$getParents(dataNodes, stochOnly = TRUE)
# ## Ensure we have both data nodes and deterministic intermediates (e.g., lifted nodes)
# simNodes <- model$getDependencies(parentNodes, self = FALSE)
# 
# nSamp = 5000
# samples_sub <- sampl_mtrx[1:nSamp,]
# n <- length(data$y)
# ppSamples <- matrix(0, nSamp, n)
# #
# # set.seed(0)
# for(i in 1:nSamp){
#   cModel[["R"]] <- samples_sub[i, "R"]
#   cModel[["C"]] <- samples_sub[i, "C"]
#   cModel[["rho[1]"]] <- samples_sub[i, "rho.1."]
#   cModel[["rho[2]"]] <- samples_sub[i, "rho.2."]
#   cModel[["sigma"]] <- samples_sub[i, "sigma"]
#   cModel[["sigma_B"]] <- samples_sub[i, "sigma_B"]
# 
# 
#   cModel$simulate(simNodes, includeData = TRUE)
#   ppSamples[i, ] <- cModel[["y"]]
# }
# 
# mean_pred <-  apply(ppSamples[,1:50],2,mean)
# #
# plot(time[1:50], mean_pred)
# plot(time[1:50], ppSamples[3,1:50])
# 
# 
# cov_disc <- ccov_sq_exp(X, rho = c(mean["rho[1]"],mean["rho[2]"]), sigma = mean['sigma_B'], delta = mean['sigma'])
# pred <- rmnorm_chol(1,mean_pred, chol(cov_disc))
# plot(time,pred[1:50])
# # ###
# # # Dicekrig vignette try
# # 
# # # FIXED RHOS
# # rho <- c(0.63952674,0.05086166)
# # #PRIOR RHOS
# # #rho <- c(mean[])
# # 
# # sigma_B <- mean["sigma_B"]
# # sigma <-  mean['sigma']
# # bias_cov <- ccov_sq_exp(X, rho, sigma_B, sigma)
# # chol_bias_cov <- chol(bias_cov)
# # mean_pred<- WK_simulate.R(flow01,Qmin,Qmax, time, Rtrue, Ctrue)
# # 
# # #rmnorm_chol(1,rep(0,length(P_obs)),chol_bias_cov)
# # bias_corr_pred <- rmnorm_chol(1,rep(mean_pred,3),chol_bias_cov)
# # 
# # plot(time,mean_pred[1:50])
# # 
# # 
# # # ###
# # #
# # #
# # #
# # # ppSamples = as.data.frame(ppSamples)
# # # ppSamples$sample = 1:nrow(ppSamples)
# # # #ppSamples <- ppSamples[1:50,1:50]
# # # colnames(ppSamples)[grep("V", colnames(ppSamples))] = time
# # # df_plot_pp = melt(ppSamples, id.var = "sample")
# # # colnames(df_plot_pp)[2:3] = c("time", "P")
# # # as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
# # # df_plot_pp$time = as.numeric.factor(df_plot_pp$time)
# # #
# # #
# # # df_obs <- data.frame(P_obs=P_obs, time=time)
# # # df_true <- data.frame(P=P, time=time)
# # # ggplot()+
# # #   geom_point(data =df_obs,aes(x=time, y=P_obs, colour="observed"), shape=16, size = 1)+
# # #   geom_line(data=df_plot_pp[df_plot_pp$sample==1:500,], aes(x=time,y=P,group=sample, color = "pp samples"), alpha=0.5, size=0.7)+
# # #   geom_line(data =df_true,aes(x=time, y=P, colour="true"),size = 1)+
# # #   theme_bw()+
# # #   scale_colour_manual(values = c("black","red", "red4")
# # #                       ,guide = guide_legend(override.aes = list(
# # #                         linetype = c( "blank", "solid", "solid"),shape = c(1,NA,NA))))+
# # #   theme(legend.title = element_blank(),legend.position = c(0.9, 0.8))
# # #
# # # # #ggsave(stringr::str_c(name,"-fit.png"))
# # # # #
# # #
