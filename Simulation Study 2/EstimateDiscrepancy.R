library(nimble)
library(DiceKriging)
library(ggpubr)
library(reshape2)
library(tidyverse)

library(coda)
library(bayesplot)
library(hexbin)
library(MASS)
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

WK_simulate.R = function(flow01, flow_min, flow_max, time, R, C){
  
  Q = flow01*(flow_max-flow_min)+flow_min
  #Q = flow01
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
  Qfft = fft(Q)
  
  Zmod = impedance_WK2(w,R,C)
  P = Re(fft(Zmod*Qfft, inverse = TRUE))/n
  return(P)
}

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
  sigma ~ dunif(0,50)
 
  P_sim[1:N] <- WK_nimble(flow01=flow01[1:N], flow_min=Qmin, flow_max=Qmax, time = time[1:N], R=R, C=C)
  
  for(i in 1:N) {
    y[i] ~ dnorm(P_sim[i], sd = sigma) #mvn
  }

})



Qmax = max(Q); Qmin = min(Q)
ind = round(seq(1,101,length.out = 50))
flow = Q[ind]; time = time[ind]
flow01 = (flow-Qmin)/(Qmax-Qmin)


#P = WK_simulate.R(flow01 = flow, flow_min = Qmin, flow_max = Qmax, time = time, R = 0.9, C=1.3)
P = WK3_simulate(flow=flow01,flow_max = Qmax, flow_min = Qmin, time=time, R=0.9, C=1.3, Z = 0.1)
plot(time,P)

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
                        
                        sigma_lower = double(),
                        sigma = double(),
                        sigma_upper = double()
                      
                        
)



# Set nr of sample runs
N_sims = 1


for(sim in 1:N_sims){
  print("sim nr")
  print(sim)
  
  sd_noise <- 3

  
  ####################################################
  # iid noise
  P_obs = rep(P,nrep) + rnorm(nrep*lP,0,sd_noise)
  ####################################################
  
  N = length(P_obs)
  
  X = matrix(c(rep(flow01,nrep),rep(time,nrep)),ncol=2)
  
  #NB for WK2 gjør om funksjonen for WKsimulate NÅ!!!!
  
  constants <- list(N = N, p=2, time = rep(time,nrep), flow01 = rep(flow01,nrep), X = X, Qmax = Qmax, Qmin = Qmin)
  
  data <- list(y = P_obs)
  #data <- list(y = rep(b[ind],nrep))
  
  inits <- list(R = 1.1, C = 1.1, sigma = 10)
  #inits <- list(sigma =2, rho  = c(0.5, 0.5), sigma_B=10)
  
  #inits$cov <- ccov_sq_exp(X, inits$rho, inits$sigma_B, inits$sigma)
  
  model <- nimbleModel(code, constants = constants, data = data, inits = inits)
  cModel <- compileNimble(model)
  conf <- configureMCMC(model)
  conf$printSamplers()
  MCMC <- buildMCMC(conf)
  cMCMC <- compileNimble(MCMC, project = cModel)
  samples <- runMCMC(cMCMC, niter = 10000, nburnin = 5000)
  samp = mcmc(samples)
  
  # obtain MAP estimate for parameters
  mean <-  apply(samp,2,mean)
  
  # store 90% credible interval
  CI <- apply(samp,2,sort)
  df_params <- rbind(df_params, data.frame(C_lower = CI[dim(samp)[1]*0.05,'C'], C = mean['C'], C_upper = CI[dim(samp)[1]*0.95,'C'], 
                                           R_lower = CI[dim(samp)[1]*0.05,'R'], R = mean['R'], R_upper = CI[dim(samp)[1]*0.95,'R'],
                                            sigma_lower = CI[dim(samp)[1]*0.05,'sigma'], sigma = mean['sigma'], sigma_upper = CI[dim(samp)[1]*0.95,'sigma']))

  

}
df <- data.frame(samp)


sampl_mtrx <-  data.matrix(df)
smpl_size <- data_frame(parameter = c("C","R","sigma"), Z = as.double(coda::effectiveSize(samples)))

dataNodes <- model$getNodeNames(dataOnly = TRUE)
parentNodes <- model$getParents(dataNodes, stochOnly = TRUE)  
## Ensure we have both data nodes and deterministic intermediates (e.g., lifted nodes)
simNodes <- model$getDependencies(parentNodes, self = FALSE)

nSamp = 5000
samples_sub <- sampl_mtrx[1:nSamp,]
n <- length(data$y)
ppSamples <- matrix(0, nSamp, n)

set.seed(0)
for(i in 1:nSamp){
  cModel[["R"]] <- samples_sub[i, "R"]
  cModel[["C"]] <- samples_sub[i, "C"] 
  cModel[["sigma"]] <- samples_sub[i, "sigma"]
  cModel$simulate(simNodes, includeData = TRUE)
  ppSamples[i, ] <- cModel[["y"]]
}


#plot(time[1:50], apply(ppSamples[,0:50],2,mean))

# 
# 
# ppSamples = as.data.frame(ppSamples)
# ppSamples$sample = 1:nrow(ppSamples)
# #ppSamples <- ppSamples[1:50,]
# colnames(ppSamples)[grep("V", colnames(ppSamples))] = time
# df_plot_pp = melt(ppSamples, id.var = "sample")
# 
# colnames(df_plot_pp)[2:3] = c("time", "P")
# as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
# df_plot_pp$time = as.numeric.factor(df_plot_pp$time)
# 
# 
# 
# 
# df_obs <- data.frame(P_obs=P_obs, time=time)
# df_true <- data.frame(P=P, time=time)
# ggplot()+
#   geom_point(data =df_obs,aes(x=time, y=P_obs, colour="observed"), shape=16, size = 1)+
#   geom_line(data=df_plot_pp[df_plot_pp$sample==1:50,], aes(x=time,y=P,group=sample, color = "pp samples"), alpha=0.5, size=0.5)+
#   geom_line(data =df_true,aes(x=time, y=P, colour="true"),size = 1)+
#   theme_bw()+
#   scale_colour_manual(values = c("black","red", "red4")
#                       ,guide = guide_legend(override.aes = list(
#                         linetype = c( "blank", "solid", "solid"),shape = c(1,NA,NA))))+
#   theme(legend.title = element_blank(),legend.position = c(0.9, 0.8))
# 
# #ggsave(stringr::str_c(name,"-fit.png"))




#################################################3
# # #calculate bias

pred <- apply(ppSamples,2,mean)

b <- P_obs - pred[1:150]
plot(rep(time,3),P_obs)
plot(rep(time,3),pred[1:150])
plot(rep(time,3),b)
#rep(1, length(b))


code_bias <- nimbleCode({
  # Fix the priors
  sigma ~ dunif(0,10) 
  sigma_B ~ dunif(0,500)
  
  rho[1] ~ dinvgamma(11,scale = 5) # l = 0 u = 1  
  rho[2] ~  dinvgamma(11,scale = 5) # l = 0 u = 1  
  
  zeros[1:N] <- rep(0,N)
  
  # adding transformed flow to covariance matrix
  X[1:N,1:p] <- matrix(c(flow01[1:N],time[1:N]),ncol=2)
  cov[1:N,1:N] <- cov_sq_exp(X = X[1:N,1:p], rho[1:p], sigma_B, sigma)
  # likelihood
  y[1:N] ~ dmnorm(zeros[1:N], cov = cov[1:N,1:N])
  
})

df_params <- data.frame(rho1_lower = double(),
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




data <- list(y = rep(b[ind],nrep))

inits <- list(sigma =3, rho  = c(0.5, 0.5), sigma_B=10)
inits$cov <- ccov_sq_exp(X, inits$rho, inits$sigma_B, inits$sigma)
model <- nimbleModel(code_bias, constants = constants, data = data, inits = inits)
cModel <- compileNimble(model)
conf <- configureMCMC(model)
conf$printSamplers()
MCMC <- buildMCMC(conf)
cMCMC <- compileNimble(MCMC, project = cModel)
samples <- runMCMC(cMCMC, niter = 10000, nburnin = 5000)
samp = mcmc(samples)

# ONE RUN CHECK
df <- data.frame(samp)

ggpairs(df,
        lower =  list(continuous = "density"))+theme_bw()




# ---------------------------
total_mean <- apply(df,2,mean)
total_mean['sigma']
params <- data.frame(parameter = c( "rho.1.", "rho.2.", "sigma", "sigma_B"),
                     Z = c(total_mean['rho.1.'],
                           total_mean['rho.2.'],
                           sd_noise,
                           total_mean['sigma_B']),
                     E = c(total_mean['rho.1.'],
                           total_mean['rho.2.'],
                           total_mean['sigma'],
                           total_mean['sigma_B']))
# ----------------------------


# with mean in plot
df_l <- df %>% select(colnames(df)) %>% gather(key="parameter", value="value")
df_l$iteration = rep(1:nrow(df), ncol(df))
par_val <- params

p <- ggplot(df_l,aes(value)) + geom_histogram(aes( y= ..density..),bins = 25)
pl_hist = p + facet_wrap(~parameter, scales = "free")+
  geom_vline(data = par_val, aes(xintercept = Z), colour = "blue", size = 1.25)+
  geom_vline(data = par_val, aes(xintercept = E), colour = "red", size = 1.25)+
  theme_bw()
pl_hist

ps <- df_l %>% ggplot(aes(x=iteration, y = value)) + geom_line()
pl_trace = ps + facet_wrap(~parameter, scales = "free") +
  geom_hline(data = par_val, aes(yintercept = Z), colour = "blue", size = 1.25)+
  theme_bw()

ggarrange(pl_hist,pl_trace, ncol = 1)


# store 90% credible interval
# CI <- apply(samp,2,sort)
# df_params <- rbind(df_params, data.frame(
#   rho1_lower = CI[dim(samp)[1]*0.05,'rho[1]'], rho1 =  mean['rho[1]'],rho1_upper = CI[dim(samp)[1]*0.95,'rho[1]'],
#   rho2_lower = CI[dim(samp)[1]*0.05,'rho[2]'],rho2 =  mean['rho[2]'],rho2_upper = CI[dim(samp)[1]*0.95,'rho[2]'],
# 
#   sigma_lower = CI[dim(samp)[1]*0.05,'sigma'], sigma = mean['sigma'], sigma_upper = CI[dim(samp)[1]*0.95,'sigma'],
#   sigma_bias_lower = CI[dim(samp)[1]*0.05,'sigma_B'], sigma_bias =  mean['sigma_B'], sigma_bias_upper = CI[dim(samp)[1]*0.95,'sigma_B']))
# 


# Check coverage (mean in CI?)

# if (between(mean['rho[1]'],df_params$rho1_lower[sim],df_params$rho1_upper[sim])) coverage[1] <- coverage[1]+1
# if (between(mean['rho[2]'],df_params$rho2_lower[sim],df_params$rho2_upper[sim])) coverage[2] <- coverage[2]+1
# 
# if (between(sd_noise,df_params$sigma_lower[sim],df_params$sigma_upper[sim])) coverage[3] <- coverage[3]+1
# if (between(mean['sigma_B'],df_params$sigma_bias_lower[sim],df_params$sigma_bias_upper[sim])) coverage[4] <- coverage[4]+1

#################################################
effectiveSize(samp)

par(mfrow = c(ncol(samples)/2,2))
for(i in 1:ncol(samples)){
  plot(samples[,i], type="l", xlab = "iterations", ylab = colnames(samples)[i], col = 1)
  lines(samples2[,i], type="l", xlab = "iterations", ylab = colnames(samples2)[i], col = 2)
  lines(samples3[,i], type="l", xlab = "iterations", ylab = colnames(samples3)[i], col = 3)
  lines(samples4[,i], type="l", xlab = "iterations", ylab = colnames(samples4)[i], col = 4)
  legend("topleft", legend = 1:4, col=1:4, pch=1) # optional legend
  
}
# with mean in plot


# just run
m<-apply(samp,2,mean)
m2<-apply(samp2,2,mean)
m3<-apply(samp3,2,mean)
m4<-apply(samp4,2,mean)


mcmc_pairs(samp, off_diag_args = list(size = .5), off_diag_fun="hex")
warnings()

# convergence check
df <- data.frame(samp)
df$Chain = rep(1, nrow(df))
df2 <- data.frame(samp2)
df2$Chain = rep(2, nrow(df2))
df3 <- data.frame(samp3)
df3$Chain = rep(3, nrow(df3)) 
df4 <- data.frame(samp4)
df4$Chain = rep(4, nrow(df4))                 
df <- rbind(df,df2,df3,df4)  
mcmc_pairs( x = df, off_diag_args = list(size = .5), off_diag_fun="hex")
warnings()

#--------------------------------

# total_mean <- apply(df_params,2,mean)
# 
# params <- data.frame(parameter = c("rho.1.", "rho.2.", "sigma", "sigma_B"), 
#                      Z = c(total_mean['rho1'],total_mean['rho2'], sd_noise, total_mean['sigma_bias']),
#                      E = c(
#                            total_mean['rho1'],
#                            total_mean['rho2'],
#                            total_mean['sigma'],
#                            total_mean['sigma_bias']))
#   
# # WITH RHO
# total_CI <- data.frame(parameter = c(rep("rho1", N_sims),rep("rho2", N_sims), rep("sigma", N_sims), rep("sigma_bias", N_sims)),
#                  lower = c(df_params$rho1_lower, df_params$rho2_lower, df_params$sigma_lower, df_params$sigma_bias_lower),
#                  upper = c(df_params$rho1_upper, df_params$rho2_upper, df_params$sigma_upper, df_params$sigma_bias_upper))
# 
# total_CI$iteration <- rep(1:nrow(df_params), ncol(CI))
# 
# name <- "bias-iid"
# 
# rho1_CI <- ggplot(total_CI[1:N_sims,], aes(rho1))+
#   geom_segment(aes(x = lower, y = iteration, xend = upper, yend = iteration))+
#   geom_vline(data = params[1,], aes(xintercept = Z), colour = "blue")+
#   geom_vline(data = params[1,], aes(xintercept = E), colour = "red")+
#   theme_bw()
# ggsave(stringr::str_c(name,"-rho1_CI.png"))
# 
# rho2_CI <- ggplot(total_CI[(N_sims+1):(2*N_sims),], aes(rho2))+
#   geom_segment(aes(x = lower, y = iteration, xend = upper, yend = iteration))+
#   geom_vline(data = params[2,], aes(xintercept = Z), colour = "blue")+
#   geom_vline(data = params[2,], aes(xintercept = E), colour = "red")+
#   theme_bw()
# ggsave(stringr::str_c(name,"-rho2_CI.png"))
# 
# sigma_CI <- ggplot(total_CI[(2*N_sims+1):(3*N_sims),], aes(sigma))+
#   geom_segment(aes(x = lower, y = iteration, xend = upper, yend = iteration))+
#   geom_vline(data = params[3,], aes(xintercept = Z), colour = "blue")+
#   geom_vline(data = params[3,], aes(xintercept = E), colour = "red")+
#   theme_bw()
# ggsave(stringr::str_c(name,"-sigma_CI.png"))
# 
# sigma_B_CI <- ggplot(total_CI[(3*N_sims+1):(4*N_sims),], aes(sigma_B))+
#   geom_segment(aes(x = lower, y = iteration, xend = upper, yend = iteration))+
#   geom_vline(data = params[4,], aes(xintercept = Z), colour = "blue")+
#   geom_vline(data = params[4,], aes(xintercept = E), colour = "red")+
#   theme_bw()
# ggsave(stringr::str_c(name,"-sigma_B_CI.png"))
# 
# 
# 
