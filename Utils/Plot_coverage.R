library(ggpubr)
library(GGally)

#EXPERIMENT 1  ###################################################################
coverage_22 <- data.frame(coverage = c(0.94, 0.92, 0.86, 0.40,0.20,0.27), parameter = c(rep(c("C", "R", "sigma"),2)), 
                          noise = c(rep("independent", 3), rep("dependent",3)), model = rep("2-2",2))
coverage_33 <- data.frame(coverage = c(0.88, 0.87, 0.89, 0.50,0.20,0.17), parameter = c(rep(c("C", "R", "sigma"),2)), 
                          noise = c(rep("independent", 3), rep("dependent",3)), model = rep("3-3",2))
coverage_32 <- data.frame(coverage = c(1.0, 1.0, 0, 1.0,0.98,0), parameter = c(rep(c("C", "R", "sigma"),2)), 
                          noise = c(rep("independent", 3), rep("dependent",3)), model = rep("3-2",2))


coverage_df <- rbind(coverage_22,coverage_33, coverage_32)


ggplot(data = coverage_df, aes(x = parameter, y = coverage, color = noise))+
  geom_point(aes(shape = model),size = 4.5)+
  geom_hline(yintercept = 0.9, size = 1, linetype = "dashed")+
  scale_color_manual(values=c('#00CC66','#CC00CC'))+
  theme_bw()


# EXPERIMENT 2 ###################################################################
coverage_fc <- data.frame(coverage = c(1, 1, 0.92, 0.8,0.97,0), parameter = c(rep(c("C", "R", "sigma"),2)), 
                          noise = c(rep("independent", 3), rep("dependent",3)), model = rep("FullCal",2) )
coverage_m1 <- data.frame(coverage = c(1, 1, 0.92, 1,1,0), parameter = c(rep(c("C", "R", "sigma"),2)), 
                          noise = c(rep("independent", 3), rep("dependent",3)), model = rep("Mod1",2))
coverage_m2 <- data.frame(coverage = c(1.0, 1.0, 0.89, 1.0,1,1), parameter = c(rep(c("C", "R", "sigma"),2)), 
                          noise = c(rep("independent", 3), rep("dependent",3)), model = rep("Mod2",2))


coverage_df_2 <- rbind(coverage_fc,coverage_m1, coverage_m2)

ggplot(data = coverage_df_2, aes(x = parameter, y = coverage, color = noise))+
  geom_point(aes(shape = model),size = 4.5)+
  geom_hline(yintercept = 0.9, size = 1, linetype = "dashed")+
  scale_color_manual(values=c('#00CC66','#CC00CC'))+
  theme_bw()

# EXPERIMENT 3  ###################################################################
coverage_fc <- data.frame(coverage = c(0.94, 0.92, 0.86, 0.40,0.20,0.27), parameter = c(rep(c("C", "R", "sigma"),2)), 
                          noise = c(rep("independent", 3), rep("dependent",3)))
coverage_m1 <- data.frame(coverage = c(0.88, 0.87, 0.89, 0.50,0.20,0.17), parameter = c(rep(c("C", "R", "sigma"),2)), 
                          noise = c(rep("independent", 3), rep("dependent",3)))
coverage_m2 <- data.frame(coverage = c(1.0, 1.0, 0, 1.0,0.98,0), parameter = c(rep(c("C", "R", "sigma"),2)), 
                          noise = c(rep("independent", 3), rep("dependent",3)))


coverage_df <- rbind(coverage_22,coverage_33, coverage_32)
coverage_df$model <- c(rep("2-2",3), rep("3-3",3), rep("3-2",3))


ggplot(data = coverage_df, aes(x = parameter, y = coverage, color = noise))+
  geom_point(aes(shape = model),size = 3)+
  geom_hline(yintercept = 0.9, size = 1, linetype = "dashed")+
  scale_color_manual(values=c('#00CC66','#CC00CC'))+
  theme_bw()
  