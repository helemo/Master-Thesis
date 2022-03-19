library(nimble)
library(tidyverse)
library(ggpubr)
library(reshape2)
library(DiceKriging)
library(bayesplot)
library(hexbin)
library(coda)
library(GGally)


## Windkessel 3 parameters model calibration

#The 3 parameters Windkessel model is represented by the following differential equation

#$$Q(t)(1 + \frac{Z}{R})+ ZC\frac{dQ(t)}{dt} = \frac{1}{R}P(t) + C\frac{dP(t)}{dt}.$$
#We want to estimate model parameters $R$ and $C$ and $Z$ given the observed inflow $Q_{obs}(t_i)$ and the observed blood pressure $P_{obs}(t_i).$ For simplicity, we assume that we observe noisy pressure data and we denote with $P^{\text{WK2}}(Q,t,R,C,Z)$ the simulated pressure data from the WK3 model given the inflow ($Q$) in time ($t$) and model parameters $R$ and $C$ and $Z$. 

#To simulate from the model we use the following time and inflow data, $Q(t)$
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
#plot(time,Q, type="l")

#Remove before github 
# df <- data.frame(time=time,Q=Q)
# ggplot(data=df, aes(x = time, y = Q, colour = "Q"))+
#   geom_line()+
#   theme_classic()
# till here

#The WK3 model simulator is the following:
  

WK3_simulate = function(flow, time, R, C, Z = NULL){
  Q = flow
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
  Zmod = impedance_WK3(w,R,C,Z_ao)
  P = Re(fft(Zmod*Qfft, inverse = TRUE))/n
  return(P)
}

#We can simulate $P(t)$ for different values of Z as follows
# 
# Rtrue=0.9;Ctrue=0.01
# Z1=0.01;Z2= 0.05;Z3 = 0.1;Z4=0.2
# Psim1 = WK3_simulate(Q,time, Rtrue, Ctrue, Z1)
# Psim2 = WK3_simulate(Q,time, Rtrue, Ctrue, Z2)
# Psim3 = WK3_simulate(Q,time, Rtrue, Ctrue, Z3)
# Psim4 = WK3_simulate(Q,time, Rtrue, Ctrue, Z4)
# Psim_df = data.frame(time=time,P1 = Psim1, P2 = Psim2, P3=Psim3, P4 = Psim4)
# #Psim_df2 = data.frame(time= rep(time,4),P = c(Psim1,Psim2,Psim3,Psim4))
# #Psim_df2$Z = c(rep(Z1,length(Psim1)), rep(Z2,length(Psim1)),rep(Z3,length(Psim1)),rep(Z4,length(Psim1)))
# 
# ggplot(data = Psim_df)+
#   geom_line(aes(x = time, y = P1, colour = "Z = 0.01"))+
#   geom_line(aes(x = time, y = P2, colour = "Z = 0.05"))+
#   geom_line(aes(x = time, y = P3, colour = "Z = 0.1"))+
#   geom_line(aes(x = time, y = P4, colour = "Z = 0.2"))+
#   theme_classic()+
#   labs(y = "P")


#Now we want to use this function in a nimble model and in order to be compatible with nimble code we have to do the following


WK3_nimble <- nimbleRcall(function(flow = double(1), time = double(1), R=double(0), C = double(0), Z = double(0) ){}, Rfun = 'WK3_simulate',
                          returnType = double(1))


#So the last chunk of code translates the R code in nimble code. When you do that you have to specify the type of inputs (`double(0)` for scalar, `double(1)` for vector and `double(2)` for matrix). The same holds for outputs (`returnType`)




#We assume that we observe the blood pressure, $Pobs(t)$ with $i.i.d.$ noise  $\varepsilon \sim N(0,\sigma^2).$ So we can write that 
#$P_{obs}(t_i)\sim N(P^{WK2}(Q(t_i), R, C), \sigma^2)$ and lets define this model in nimble. 

#First we will define the probabilistic model. The model is similar to the linear regression model we defined before with the main difference that the linear predictor is replaced by the simulator output $P^{WK3}(Q(t_i), R, C, Z).$ As before, we have to define priors for the physical parameters $R,C,Z$ and also for the $\sigma$ parameter. 


code_WK3 <- nimbleCode({
  # Priors
  C ~ dunif(0.5,3)
  R ~ dunif(0.5,3)
  Z ~ dunif(0.01,0.2)
  sigma ~ dunif(0, 30) 
  
  # Simulator
  P_sim[1:N] <- WK3_nimble(flow = Q[1:N], time = time[1:N], R=R, C=C, Z=Z)
  
  # likelihood
  for(i in 1:N) {
    y[i] ~ dnorm(P_sim[i], sd = sigma)
  }
})



#To create observed data, we simualte data from the WK3 model and we add $i.i.d.$ noise $\varepsilon \sim N(0,3^2)$


ind = round(seq(1,101,length.out = 31))
Q = Q[ind]; time = time[ind]

Rtrue=0.9; Ctrue=1.3; Ztrue = 0.2
P = WK3_nimble(flow = Q, time = time, R = Rtrue, C = Ctrue, Z = Ztrue)
N=length(P)
set.seed(0)

# fixed noise for Z experiment
fixed_noise <-  rnorm(N, 0, 3)

df_params <- data.frame(C_lower = double(),
                        C = double(),
                        C_upper = double(),
                        
                        R_lower = double(),
                        R = double(),
                        R_upper = double(),
                        
                        Z_lower = double(),
                        Z = double(),
                        Z_upper = double(),
                        
                        sigma_lower = double(),
                        sigma = double(),
                        sigma_upper = double())

##############################

N_sims = 1
# df <- data.frame(C = double(), 
#                  R = double(), 
#                  Z = double(),
#                  sigma = double()) 

coverage <- c(0,0,0,0)

for(sim in 1:N_sims){
  
  print("sim nr")
  print(sim)
  
  sd_noise = 3
  
  ################################################
  # iid noise
  # eps = rnorm(N, 0, sd_noise)
  # plot(time, eps)

  P_obs = P + fixed_noise #eps
  plot(time, P_obs)
  
  ################################################
  # dependent noise
  # theta <- 0.4
  # sigma <- sd_noise #sd_noise*0.5
  # eps_model <- km(design = data.frame(x = time),  response = rep(1, length(time)), covtype = "matern5_2",
  #                 coef.cov = theta, coef.var = sigma^2, nugget = 2)
  # eps <- simulate(eps_model, nsim = 1, newdata = data.frame(x = time))
  # plot(time,eps)
  # 
  # P_obs = P + eps[1,]
  # plot(time, P_obs)

  ################################################
  
  
  #Then we define the nimble model
  
  constants = list(N = N, time = time, Q = Q)
  data = list(y = P_obs)
  inits = list(R = 1, C = 1,Z = 0.05, sigma = 1)
  model = nimbleModel(code_WK3, constants = constants, data = data, inits = inits)
  
  
  #It is important to check the message above for warnings and errors.
  
  
  #As before we define the sampling algorithm
  
  
  cModel = compileNimble(model)
  printErrors()
  conf = configureMCMC(model)
  conf$printSamplers()
  MCMC = buildMCMC(conf)
  cMCMC = compileNimble(MCMC, project = cModel)
  
  #and we run the MCMC algorithm
  samples = runMCMC(cMCMC, niter = 10000, nburnin = 5000)
  samp = mcmc(samples)
  
  mean <-  apply(samp,2,mean)
  
  CI <- apply(samp,2,sort)
  #samp
  df_params <- rbind(df_params, data.frame( C_lower = CI[dim(samp)[1]*0.05,'C'], C = mean['C'], C_upper = CI[dim(samp)[1]*0.95,'C'], 
                                            R_lower = CI[dim(samp)[1]*0.05,'R'], R = mean['R'], R_upper = CI[dim(samp)[1]*0.95,'R'],
                                            Z_lower = CI[dim(samp)[1]*0.05,'Z'], Z = mean['Z'], Z_upper = CI[dim(samp)[1]*0.95,'Z'],
                                            sigma_lower = CI[dim(samp)[1]*0.05,'sigma'], sigma = mean['sigma'], sigma_upper = CI[dim(samp)[1]*0.95,'sigma']))
  
  
  # Check coverage Not really needed?
  if (between(Ctrue,df_params$C_lower[sim],df_params$C_upper[sim])) coverage[1] <- coverage[1]+1
  if (between(Rtrue,df_params$R_lower[sim],df_params$R_upper[sim])) coverage[2] <- coverage[2]+1
 
  if (between(sd_noise,df_params$sigma_lower[sim],df_params$sigma_upper[sim])) coverage[3] <- coverage[3]+1
  if (between(Ztrue,df_params$Z_lower[sim],df_params$Z_upper[sim])) coverage[4] <- coverage[4]+1
  
  # 
  # df_samples <- data.frame(samples)
  # 
  # 
  # 
  # N_samples = nrow(df_samples)
  # 
  # # Calculate sample posterior mean
  # C = sum(df_samples['C'])/N_samples
  # R = sum(df_samples['R'])/N_samples
  # Z = sum(df_samples['Z'])/N_samples
  # sigma = sum(df_samples['sigma'])/N_samples
  #mean <- rbind(mean, data.frame(C_mean, R_mean, sigma_mean))
  
  # Calculate sample posterior 90%CI
  # C_lower <- sort(df_samples$C)[N_samples*0.05]; C_upper <- sort(df_samples$C)[N_samples*0.95]
  # R_lower <- sort(df_samples$R)[N_samples*0.05]; R_upper <- sort(df_samples$R)[N_samples*0.95]
  # Z_lower <- sort(df_samples$Z)[N_samples*0.05]; Z_upper <- sort(df_samples$Z)[N_samples*0.95]
  # sigma_lower <- sort(df_samples$sigma)[N_samples*0.05]; sigma_upper <- sort(df_samples$sigma)[N_samples*0.95]
  # df_params <- rbind(df_params, data.frame(C_lower, C, C_upper, 
  #                                          R_lower, R, R_upper,
  #                                          Z_lower, Z, Z_upper,
  #                                          sigma_lower, sigma, sigma_upper))
  # 
  # # Check coverage
  # if (between(Ctrue,C_lower,C_upper)) coverage[1] <- coverage[1]+1
  # if (between(Rtrue,R_lower,R_upper)) coverage[2] <- coverage[2]+1
  # if (between(sd_noise,sigma_lower,sigma_upper)) coverage[3] <- coverage[3]+1
  # if (between(Ztrue,Z_lower,Z_upper)) coverage[4] <- coverage[4]+1
  # 
  
  # df <- rbind(df,df_samples)
  
}

# ONE RUN CHECK
# posteriors after 1 run
df <- data.frame(samp)

ggpairs(df,
        lower =  list(continuous = "density"))+theme_bw()

effectiveSize(samp)
# par(mfrow = c(ncol(samples)/2,2))
# for(i in 1:ncol(samples)){
#   plot(samples[,i], type="l", xlab = "iterations", ylab = colnames(samples)[i])
# }

# apply(samp,2,mean)

#mcmc_pairs(samp, off_diag_args = list(size = .5), off_diag_fun="hex")


# ---------------------------
total_mean <- apply(df_params,2,mean)
total_mean

params <- data.frame(parameter = c("C", "R", "sigma","Z"), Z = c(Ctrue, Rtrue, sd_noise,Ztrue),
                     E = c(total_mean['C'],
                           total_mean['R'],
                           total_mean['sigma'],
                           total_mean['Z']
                     ))



# ----------------------------




# with mean in plot
df_l <- df %>% select(colnames(df)) %>% gather(key="parameter", value="value")
df_l$iteration = rep(1:nrow(df), ncol(df))
par_val <- params #data.frame(parameter = c("C", "R", "sigma", "lambda", "sigma_bias"), Z = c(Ctrue, Rtrue,sd_noise, lambda_true, sigma_bias_true), E =c(avg_C, avg_R, avg_sigma,avg_lambda,avg_sigma_bias) ) #true and estimated values

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


###############################################

CI <- data.frame(parameter = c(rep("C",N_sims),rep("R",N_sims),rep("Z",N_sims), rep("sigma", N_sims)),
                 lower = c(df_params$C_lower, df_params$R_lower,df_params$Z_lower, df_params$sigma_lower),
                 upper = c(df_params$C_upper, df_params$R_upper,df_params$Z_upper,df_params$sigma_upper))
CI$iteration <- rep(1:nrow(df_params), 4) #4 = nr of params

name <- "Spam"
CI

C_CI <- ggplot(CI[1:N_sims,], aes(C))+
  geom_segment(aes(x = lower, y = iteration, xend = upper, yend = iteration))+
  geom_vline(data = params[1,], aes(xintercept = Z), colour = "blue")+
  geom_vline(data = params[1,], aes(xintercept = E), colour = "red")+
  theme_bw()
ggsave(stringr::str_c(name,"-C_CI.png"))

R_CI <- ggplot(CI[(N_sims+1):(2*N_sims),], aes(R))+
  geom_segment(aes(x = lower, y = iteration, xend = upper, yend = iteration))+
  geom_vline(data = params[2,], aes(xintercept = Z), colour = "blue")+
  geom_vline(data = params[2,], aes(xintercept = E), colour = "red")+
  theme_bw()
ggsave(stringr::str_c(name,"-R_CI.png"))

sigma_CI <- ggplot(CI[(3*N_sims+1):(4*N_sims),], aes(sigma))+
  geom_segment(aes(x = lower, y = iteration, xend = upper, yend = iteration))+
  geom_vline(data = params[3,], aes(xintercept = Z), colour = "blue")+
  geom_vline(data = params[3,], aes(xintercept = E), colour = "red")+
  theme_bw()
ggsave(stringr::str_c(name,"-sigma_CI.png"))

Z_CI <- ggplot(CI[(2*N_sims+1):(3*N_sims),], aes(Z))+
  geom_segment(aes(x = lower, y = iteration, xend = upper, yend = iteration))+
  geom_vline(data = params[4,], aes(xintercept = Z), colour = "blue")+
  geom_vline(data = params[4,], aes(xintercept = E), colour = "red")+
  theme_bw()
ggsave(stringr::str_c(name,"-Z_CI.png"))


df_coverage <- data.frame(parameter = c("C", "R", "sigma", "Z"),coverage = c(coverage[1],coverage[2],coverage[3], coverage[4]))

cover <- ggplot(df_coverage, aes(parameter, coverage))+
  geom_point()
ggsave(stringr::str_c(name,"-coverage.png"))


###############################


#Mixing is good! The effective sample size is also ok. For better parameter estimates it might be useful to have replicates.

#sample_size <- coda::effectiveSize(samples)
#plot(c("C","R","Z","sigma"), sample_size)

# smpl_size <- data_frame(parameter = c("C","R","Z","sigma"), Z = as.double(coda::effectiveSize(samples)))


# ggplot(smpl_size, aes(parameter,Z)) + geom_point()
# ggsave(stringr::str_c(name,"-ess.png"))

#############
# 
# dataNodes <- model$getNodeNames(dataOnly = TRUE)
# parentNodes <- model$getParents(dataNodes, stochOnly = TRUE)  
# ## Ensure we have both data nodes and deterministic intermediates (e.g., lifted nodes)
# simNodes <- model$getDependencies(parentNodes, self = FALSE)
# 
# nSamp = 5000
# samples_sub <- samples[1:nSamp,]
# n <- length(data$y)
# ppSamples <- matrix(0, nSamp, n)
# 
# set.seed(0)
# for(i in 1:nSamp){
#   cModel[["R"]] <- samples_sub[i, "R"]
#   cModel[["C"]] <- samples_sub[i, "C"] 
#   cModel[["Z"]] <- samples_sub[i, "Z"] 
#   cModel[["sigma"]] <- samples_sub[i, "sigma"]
#   cModel$simulate(simNodes, includeData = TRUE)
#   ppSamples[i, ] <- cModel[["y"]]
# }
# 
# 
# 
# ppSamples = as.data.frame(ppSamples)
# ppSamples$sample = 1:nrow(ppSamples)
# ppSamples <- ppSamples[1:5,]
# 
# colnames(ppSamples)[grep("V", colnames(ppSamples))] = time
# df_plot_pp = melt(ppSamples, id.var = "sample")
# colnames(df_plot_pp)[2:3] = c("time", "P")
# as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
# df_plot_pp$time = as.numeric.factor(df_plot_pp$time)
# 
# 
# df_obs <- data.frame(P_obs=P_obs, time=time)
# df_true <- data.frame(P=P, time=time)
# ggplot()+
#   geom_line(data=df_plot_pp[df_plot_pp$sample==1:5,], aes(x=time,y=P,group=sample, color = "pp samples"), alpha=0.3, size=0.5)+
#   geom_line(data =df_true,aes(x=time, y=P, colour="true"))+
#   theme_bw()+
#   geom_point(data =df_obs,aes(x=time, y=P_obs, colour="observed"), shape=16, size = 1)+
#   scale_colour_manual(values = c("black","red", "red4")
#                       ,guide = guide_legend(override.aes = list(
#                         linetype = c( "blank", "solid", "solid"),shape = c(1,NA,NA))))+
#   theme(legend.title = element_blank(),legend.position = c(0.9, 0.8))
