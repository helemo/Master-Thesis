library(ggpubr)
library(GGally)

#Plot CI and MAP fir Simulation Studies

########################################################################
# SIMLATION STUDY 2
##########################
# C
avg_param <- data.frame(param = rep("C",6), model = rep(c("Full Calib", "Mod 1", "Mod 2"),2), 
                        avg_90_CI = c(1.25,1,1.2,1.25,1,1.2), 
                        CI_end = c(1.55,2.2,2.25,1.8,2.4,2.5), 
                        bias_MAP = c(1.33,1.37,1.6,1.35,1.53,1.62),
                        quant_bias = c(0.03, 0.07,0.3, 0.05,0.23,0.32),
                        simulation_number = c("2.3","2.5","2.7","2.4","2.6","2.8"),
                        noise = c(rep("independent",3), rep("dependent",3)))

ggplot(data = avg_param)+
  geom_segment(aes(x = avg_90_CI, xend = CI_end, y =  simulation_number, yend =  simulation_number, color = model, linetype= noise), size = 4)+
  geom_vline(xintercept = 1.3, size = 1.5)+
  geom_point(aes(x = bias_MAP, y =  simulation_number), shape = 4, size = 8)+
  theme_bw()



##########################

# R
avg_param <- data.frame(param = rep("R",6), model = rep(c("Full Calib", "Mod 1", "Mod 2"),2), 
                        avg_90_CI = c(0.8,0.93,0.85,0.8,0.8,0.75), 
                        CI_end = c(1.3,1.06,1.15,1.4,1.15,1.3), 
                        bias_MAP = c(1.05,1.005,1.01,1.06,1.01,1.03), 
                        quant_bias = c(0.05, 0.005,0.01,0.06,0.01,0.03),
                        simulation_number = c("2.3","2.5","2.7","2.4","2.6","2.8"),
                        noise = c(rep("independent",3), rep("dependent",3)))

ggplot(data = avg_param)+
  geom_segment(aes(x = avg_90_CI, xend = CI_end, y =  simulation_number, yend =  simulation_number, color = model, linetype= noise), size = 4)+
  geom_vline(xintercept = 1, size = 1.5)+
  geom_point(aes(x = bias_MAP, y =  simulation_number), shape = 4, size = 8)+
  theme_bw()



##########################

# sigma
avg_param <- data.frame(param = rep("sigma",6), model = rep(c("Full Calib", "Mod 1", "Mod 2"),2), 
                        avg_CI = c(2.75,2.8,2.75, 0.6,0.8,2.92), 
                        CI_end = c(3.25,3.3,3.3,1.1,1.2,3.02), 
                        bias_MAP = c(3.001,3.001,3.05,0.8,1.1,2.96), 
                        quant_bias = c(0.001, 0.001,0.05,-2.2,-1.9,-0.04), 
                        simulation_number = c("2.3","2.5","2.7","2.4","2.6","2.8"),
                        noise = c(rep("independent",3), rep("dependent",3)))

ggplot(data = avg_param)+
  geom_segment(aes(x = avg_CI, xend = CI_end, y = simulation_number, yend = simulation_number, color = model, linetype = noise), size = 4)+
  geom_vline(xintercept = 3, size = 1.5)+
  geom_point(aes(x = bias_MAP, y = simulation_number), shape = 4, size = 8)+
  theme_bw()


##########################

# sigmaB
avg_param <- data.frame(param = rep("sigmaB",6), model = rep(c("Full Calib", "Mod 1", "Mod 2"),2), 
                        avg_CI = c(25,12,17,20,20,26), 
                        CI_end = c(60,18,37,50,30,46), 
                        bias_MAP = c(37,14,24,30,22,34), 
                        simulation_number = c("2.3","2.5","2.7","2.4","2.6","2.8"),
                        noise = c(rep("independent",3), rep("dependent",3)))

ggplot(data = avg_param)+
  geom_segment(aes(x = avg_CI, xend = CI_end, y = simulation_number, yend = simulation_number, color = model, linetype = noise), size = 4)+
  geom_point(aes(x = bias_MAP, y = simulation_number), shape = 4, size = 8)+
  theme_bw()

##########################

# rho1
avg_param <- data.frame(param = rep("rho1",4), model = rep(c("Full Calib", "Mod 2"),2), 
                        avg_CI = c(0.55,0.49,0.3,0.55), 
                        CI_end = c(1,0.9,0.85,1.05), 
                        bias_MAP = c(0.78,0.63,0.58,0.73), 
                        simulation_number = c("2.3","2.7","2.4","2.8"),
                        noise = c(rep("independent",2), rep("dependent",2)))

ggplot(data = avg_param)+
  geom_segment(aes(x = avg_CI, xend = CI_end, y = simulation_number, yend = simulation_number, color = model, linetype = noise), size = 4)+
  geom_point(aes(x = bias_MAP, y = simulation_number), shape = 4, size = 8)+
  theme_bw()


##########################

# rho2
avg_param <- data.frame(param = rep("rho2",4), model = rep(c("Full Calib", "Mod 2"),2), 
                        avg_CI = c(0.62,0.15,0.35,0.2), 
                        CI_end = c(1,0.4,0.85,0.5), 
                        bias_MAP = c(0.83,0.24,0.45,0.3), 
                        simulation_number = c("2.3","2.7","2.4","2.8"),
                        noise = c(rep("independent",2), rep("dependent",2)))

ggplot(data = avg_param)+
  geom_segment(aes(x = avg_CI, xend = CI_end, y = simulation_number, yend = simulation_number, color = model, linetype = noise), size = 4)+
  geom_point(aes(x = bias_MAP, y = simulation_number), shape = 4, size = 8)+
  theme_bw()


########################################################################
# SIMLATION STUDY 1

##########################
# C
avg_param <- data.frame(param = rep("C",6), model = rep(c("2-2", "3-3", "3-2"),2), 
                        avg_90_CI = c(1.2,1.2,0.9,1.25,1.26,0.9), 
                        CI_end = c(1.45,1.4,1.8,1.42,1.40,1.8), 
                        bias_MAP = c(1.31,1.305,1.2,1.31,1.29,1.2),
                        quant_bias = c(0.01, 0.005,0.1, 0.01,0.01,0.1),
                        simulation_number = c("1.1","1.2","1.3","1.4","1.5","1.6"),
                        noise = c(rep("independent",3), rep("dependent",3)))

ggplot(data = avg_param)+
  geom_segment(aes(x = avg_90_CI, xend = CI_end, y =  simulation_number, yend =  simulation_number, color = model, linetype= noise), size = 4)+
  geom_vline(xintercept = 1.3, size = 1.5)+
  geom_point(aes(x = bias_MAP, y =  simulation_number), shape = 4, size = 8)+
  theme_bw()



##########################

# R
avg_param <- data.frame(param = rep("R",6), model = rep(c("2-2", "3-3", "3-2"),2), 
                        avg_90_CI = c(0.885,0.885,0.95,0.89,0.898,0.93), 
                        CI_end = c(0.913,0.918,1.04,0.91,0.905,1.07), 
                        bias_MAP = c(0.899,0.902,0.995,0.89,0.92,1.01), 
                        quant_bias = c(0.001, 0.002,0.05,0.01,0.02,0.01),
                        simulation_number = c("1.1","1.2","1.3","1.4","1.5","1.6"),
                        noise = c(rep("independent",3), rep("dependent",3)))

ggplot(data = avg_param)+
  geom_segment(aes(x = avg_90_CI, xend = CI_end, y =  simulation_number, yend =  simulation_number, color = model, linetype= noise), size = 4)+
  geom_vline(xintercept = 1, size = 1.5)+
  geom_vline(xintercept = 0.9, size = 1.5, linetype = 2)+
  geom_point(aes(x = bias_MAP, y =  simulation_number), shape = 4, size = 8)+
  theme_bw()



##########################

# sigma
avg_param <- data.frame(param = rep("sigma",6), model = rep(c("2-2", "3-3", "3-2"),2), 
                        avg_CI = c(2.5,2.3,17,1.7,1.5,16), 
                        CI_end = c(4.5,4.3,23,2.8,2.8,24), 
                        bias_MAP = c(3.2,3.1,19,2.2,1.9,18), 
                        quant_bias = c(0.001, 0.001,0.05,-2.2,-1.9,-0.04), 
                        simulation_number = c("1.1","1.2","1.3","1.4","1.5","1.6"),
                        noise = c(rep("independent",3), rep("dependent",3)))

ggplot(data = avg_param)+
  geom_segment(aes(x = avg_CI, xend = CI_end, y = simulation_number, yend = simulation_number, color = model, linetype = noise), size = 4)+
  geom_vline(xintercept = 3, size = 1.5)+
  geom_point(aes(x = bias_MAP, y = simulation_number), shape = 4, size = 8)+
  theme_bw()


##########################
