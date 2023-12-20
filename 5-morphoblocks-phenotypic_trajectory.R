#CH.5 - Phenotypic trajectory analyses on morphoblocks PC scores

#LOAD LIBRARIES ----
#always do this first!!
library(geomorph) 
library(geiger)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(ggrepel)
library(gginnards)
library(ggphylomorpho)
library(ggfortify)
library(RColorBrewer) 
library(borealis)
library(ggthemes)
library(ggpubr)
library(ggplotify)
library(Morpho)
library(rphylopic)
library(png)
library(gridExtra)
library(phytools)
library(evomap)
library(reshape2)
library(scales)

#TRAJECTORY ANALYSIS ----
#Need multiple observations per group/vertebra and all have to be represented (all categories in each group)

##Shape changes along column - Regionalization ----
#First perform lm.rrpp to create linear model that describes what we are trying to test

###Age ----
#Shape changes along column at different ages
fit_shape_column_age <- lm.rrpp(pcscores ~ age * vertebra, iter = 999, data = gdf_all, RRPP = F)

#Check that there is a significant correlation
summary(fit_shape_column_age)

#Save results to file
sink("Output/fit_shape_column_age.txt")
summary(fit_shape_column_age)
sink() 

#Use fit to calculate trajectories
trajectory_column_age <- trajectory.analysis(fit_shape_column_age, groups = gdf_all$age, traj.pts = gdf_all$vertebra, 
                                        pca = TRUE) 

#View results
#Magnitude differences between trajectories, standard summary - are trajectories different in length?
trajectory_column_age_MD <- summary(trajectory_column_age, show.trajectories = TRUE, attribute = "MD") 
trajectory_column_age_MD 
#Trajectory correlations -  are trajectories different in angle/direction?
trajectory_column_age_TC <- summary(trajectory_column_age, show.trajectories = TRUE, attribute = "TC", angle.type = "deg")
trajectory_column_age_TC
#Trajectory shape differences - are trajectories different in shape?
trajectory_column_age_SD <- summary(trajectory_column_age, show.trajectories = TRUE, attribute = "SD") 
trajectory_column_age_SD

#Save results to file
sink("Output/trajectory_column_age.txt")
print("Magnitude difference (absolute difference between path distances) - length")
summary(trajectory_column_age, show.trajectories = F, attribute = "MD") 
print("Correlations (angles) between trajectories - direction")
summary(trajectory_column_age, show.trajectories = F, attribute = "TC", angle.type = "deg")
print("Shape differences between trajectory vectors - shape")
summary(trajectory_column_age, show.trajectories = F, attribute = "SD") 
sink() 

#Plot results - PCA of fitted values
trajectory_column_age_plot <- plot(trajectory_column_age, main = "Trajectories shape change on column by age",  pch = shapes, #title and type of point to be used
                              col = c("darkgray","lightgray"), bg = c("darkgray","lightgray"), cex = 1, font.main = 2) #improve graphics
#Add line between groups
add.trajectories(trajectory_column_age_plot, 
                 traj.pch = shapes, traj.col = 1, traj.lty = 1, traj.lwd = 1, traj.cex = 1.5, traj.bg = 1, 
                 start.bg = "green", end.bg = "red") #trajectory line graphics
#Add legend to see which trajectory belongs to each group
legend(list(x = -1.3,y = -0.6), legend = c("adult"  ,  "juvenile", "newborn"), pch =  shapes, pt.bg = 1, cex = 0.8)

##Make better PCA plot using ggplot
#Read PC scores as tibble
trajectory_column_age_pcscores <- as_tibble(trajectory_column_age_plot [["pc.points"]])

#Add group names and other attributes to tibble as columns
trajectory_column_age_pcscores <- trajectory_column_age_pcscores %>% mutate(vertebra = gdf_all$vertebra, age = gdf_all$age)
glimpse(trajectory_column_age_pcscores)

#Calculate means of PC1 and PC2 at each vertebra per group to draw trajectories
trajectory_column_age_pcscores_means <- trajectory_column_age_pcscores %>% group_by(age, vertebra) %>%
  summarise_at(vars(PC1, PC2), list(mean = mean))              #function for means, both columns
glimpse(trajectory_column_age_pcscores_means)

#Rename columns so they are easier to use for plot
trajectory_column_age_pcscores_means <- trajectory_column_age_pcscores_means %>% rename(x = PC1_mean, y = PC2_mean)
glimpse(trajectory_column_age_pcscores_means)

#Nice plot
trajectory_column_age_ggplot <- ggplot(trajectory_column_age_pcscores, aes(x = PC1, y = PC2, shape = vertebra, group = age))+
  geom_point(size = 2.2, colour = "black", fill = "darkgray", alpha = 0.4, show.legend = F)+
  geom_point(data = trajectory_column_age_pcscores_means, aes(x = x, y = y, fill = age, shape = vertebra, group = age), colour = "black",
             size = 6, inherit.aes = F, alpha = 0.8)+
  geom_path(data = trajectory_column_age_pcscores_means, aes(x = x, y = y, colour = age, group = age), inherit.aes = F, linewidth = 1,
           linejoin = 'mitre', show.legend = F)+
  scale_shape_manual(name = "Vertebra", labels = c("Thoracic", "Lumbar", "Caudal"), values = shapes)+
  scale_colour_manual(name = "Age", labels = c("Adult"  ,  "Juvenile", "Neonate"), #copy from as.factor(genera)
                      values = mypalette_age, aesthetics = c("colour", "fill"))+
  theme_bw(base_size = 14)+
  xlab(paste0("PC 1 (",round(as.numeric(trajectory_column_age[["pca"]][["sdev"]][1]^2/sum(trajectory_column_age[["pca"]][["sdev"]]^2)*100), digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(as.numeric(trajectory_column_age[["pca"]][["sdev"]][2]^2/sum(trajectory_column_age[["pca"]][["sdev"]]^2)*100), digits = 2),"%)"))+
  theme(legend.key = element_blank(), legend.title = element_text(size = 15, face = "bold"), legend.text = element_text(size = 14),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  guides(shapes = guide_legend(override.aes = list(colour = "black", fill = "darkgray", alpha = 1)),
    color = guide_legend(override.aes = list(shape = 23, colour = "black", fill = mypalette_age)))

#Visualize plot and save as PDF using menu in bar on the right
trajectory_column_age_ggplot

###Sex ----
#Shape changes along column at different sexs
fit_shape_column_sex <- lm.rrpp(pcscores ~ sex * vertebra, iter = 999, data = gdf_all, RRPP = F)

#Check that there is a significant correlation
summary(fit_shape_column_sex)

#Save results to file
sink("Output/fit_shape_column_sex.txt")
summary(fit_shape_column_sex)
sink() 

#Use fit to calculate trajectories
trajectory_column_sex <- trajectory.analysis(fit_shape_column_sex, groups = gdf_all$sex, traj.pts = gdf_all$vertebra, 
                                             pca = TRUE) 

#View results
#Magnitude differences between trajectories, standard summary - are trajectories different in length?
trajectory_column_sex_MD <- summary(trajectory_column_sex, show.trajectories = TRUE, attribute = "MD") 
trajectory_column_sex_MD
#Trajectory correlations -  are trajectories different in angle/direction?
trajectory_column_sex_TC <- summary(trajectory_column_sex, show.trajectories = TRUE, attribute = "TC", angle.type = "deg")
trajectory_column_sex_TC
#Trajectory shape differences - are trajectories different in shape?
trajectory_column_sex_SD <- summary(trajectory_column_sex, show.trajectories = TRUE, attribute = "SD") 
trajectory_column_sex_SD

#Save results to file
sink("Output/trajectory_column_sex.txt")
print("Magnitude difference (absolute difference between path distances) - length")
summary(trajectory_column_sex, show.trajectories = F, attribute = "MD") 
print("Correlations (angles) between trajectories - direction")
summary(trajectory_column_sex, show.trajectories = F, attribute = "TC", angle.type = "deg")
print("Shape differences between trajectory vectors - shape")
summary(trajectory_column_sex, show.trajectories = F, attribute = "SD") 
sink() 

#Plot results - PCA of fitted values
trajectory_column_sex_plot <- plot(trajectory_column_sex, main = "Trajectories shape change on column by sex",  pch = shapes, #title and type of point to be used
                                   col = c("darkgray","lightgray"), bg = c("darkgray","lightgray"), cex = 1, font.main = 2) #improve graphics
#Add line between groups
add.trajectories(trajectory_column_sex_plot, 
                 traj.pch = shapes, traj.col = 1, traj.lty = 1, traj.lwd = 1, traj.cex = 1.5, traj.bg = 1, 
                 start.bg = "green", end.bg = "red") #trajectory line graphics
#Add legend to see which trajectory belongs to each group
legend(list(x = -1.3,y = -0.35), legend = c( "F"   ,    "M"      , "unknown"), pch =  shapes, pt.bg = 1, cex = 0.8)

###Plot
trajectory_column_sex_pcscores <- as_tibble(trajectory_column_sex_plot [["pc.points"]])

#Add group names and other attributes to tibble as columns
trajectory_column_sex_pcscores <- trajectory_column_sex_pcscores %>% mutate(vertebra = gdf_all$vertebra, sex = gdf_all$sex)
glimpse(trajectory_column_sex_pcscores)

#Calculate means of PC1 and PC2 at each vertebra per group to draw trajectories
trajectory_column_sex_pcscores_means <- trajectory_column_sex_pcscores %>% group_by(sex, vertebra) %>%
  summarise_at(vars(PC1, PC2), list(mean = mean))              #function for means, both columns
glimpse(trajectory_column_sex_pcscores_means)

#Rename columns so they are easier to use for plot
trajectory_column_sex_pcscores_means <- trajectory_column_sex_pcscores_means %>% rename(x = PC1_mean, y = PC2_mean)
glimpse(trajectory_column_sex_pcscores_means)

#Plot including unknwon
trajectory_column_sex_ggplot <- ggplot(trajectory_column_sex_pcscores, aes(x = PC1, y = PC2, shape = vertebra, group = sex))+
  geom_point(size = 2.2, colour = "black", fill = "darkgray", alpha = 0.4, show.legend = F)+
  geom_point(data = trajectory_column_sex_pcscores_means, aes(x = x, y = y, fill = sex, shape = vertebra, group = sex), colour = "black",
             size = 6, inherit.aes = F, alpha = 0.8)+
  geom_path(data = trajectory_column_sex_pcscores_means, aes(x = x, y = y, colour = sex, group = sex), inherit.aes = F, linewidth = 1,
            linejoin = 'mitre', show.legend = F)+
  scale_shape_manual(name = "Vertebra", labels = c("Thoracic", "Lumbar", "Caudal"), values = shapes)+
  scale_colour_manual(name = "Sex", labels = c("F"  ,  "M", "Unknown"), #copy from as.factor(genera)
                      values = c(mypalette_sex,  "chocolate"), aesthetics = c("colour", "fill"))+
  theme_bw()+
  xlab(paste0("PC 1 (",round(as.numeric(trajectory_column_sex[["pca"]][["sdev"]][1]^2/sum(trajectory_column_sex[["pca"]][["sdev"]]^2)*100), digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(as.numeric(trajectory_column_sex[["pca"]][["sdev"]][2]^2/sum(trajectory_column_sex[["pca"]][["sdev"]]^2)*100), digits = 2),"%)"))+
  theme(legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  guides(shapes = guide_legend(override.aes = list(colour = "black", fill = "darkgray", alpha = 1)),
         color = guide_legend(override.aes = list(shape = 23, colour = "black", fill = c(mypalette_sex,  "chocolate"))))

#Visualize plot and save as PDF using menu in bar on the right
trajectory_column_sex_ggplot

#Plot M + F only
trajectory_column_sex_ggplot2 <- ggplot(trajectory_column_sex_pcscores[!trajectory_column_sex_pcscores$sex%in% c("unknown"),], aes(x = PC1, y = PC2, shape = vertebra, group = sex))+
  geom_point(size = 2.2, colour = "black", fill = "darkgray", alpha = 0.4, show.legend = F)+
  geom_point(data = trajectory_column_sex_pcscores_means[!trajectory_column_sex_pcscores_means$sex%in% c("unknown"),], aes(x = x, y = y, fill = sex, shape = vertebra, group = sex), colour = "black",
             size = 6, inherit.aes = F, alpha = 0.8)+
  geom_path(data = trajectory_column_sex_pcscores_means[!trajectory_column_sex_pcscores_means$sex%in% c("unknown"),], aes(x = x, y = y, colour = sex, group = sex), inherit.aes = F, linewidth = 1,
            linejoin = 'mitre', show.legend = F)+
  scale_shape_manual(name = "Vertebra", labels = c("Thoracic", "Lumbar", "Caudal"), values = shapes)+
  scale_colour_manual(name = "Sex", labels = c("Females"  ,  "Males"), #copy from as.factor(genera)
                      values = mypalette_sex, aesthetics = c("colour", "fill"))+
  theme_bw(base_size = 14)+
  xlab(paste0("PC 1 (",round(as.numeric(trajectory_column_sex[["pca"]][["sdev"]][1]^2/sum(trajectory_column_sex[["pca"]][["sdev"]]^2)*100), digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(as.numeric(trajectory_column_sex[["pca"]][["sdev"]][2]^2/sum(trajectory_column_sex[["pca"]][["sdev"]]^2)*100), digits = 2),"%)"))+
  theme(legend.key = element_blank(), legend.title = element_text(size = 15, face = "bold"), legend.text = element_text(size = 14),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  guides(shapes = guide_legend(override.aes = list(colour = "black", fill = "darkgray", alpha = 1)),
         color = guide_legend(override.aes = list(shape = 23, colour = "black", fill = mypalette_sex)))

#Visualize plot and save as PDF using menu in bar on the right
trajectory_column_sex_ggplot2


##Shape changes per vertebra type - Growth trajectories and sexual dimorphism ----
#First perform lm.rrpp to create linear model that describes what we are trying to test

###Age ----
#Shape changes along vertebra at different ages
fit_shape_vertebra_age <- lm.rrpp(pcscores ~ vertebra * age, iter = 999, data = gdf_all, RRPP = F)

#Check that there is a significant correlation
summary(fit_shape_vertebra_age)

#Save results to file
sink("Output/fit_shape_vertebra_age.txt")
summary(fit_shape_vertebra_age)
sink() 

#Use fit to calculate trajectories
trajectory_vertebra_age <- trajectory.analysis(fit_shape_vertebra_age, groups = gdf_all$vertebra, traj.pts = gdf_all$age, 
                                             pca = TRUE) 

#View results
#Magnitude differences between trajectories, standard summary - are trajectories different in length?
trajectory_vertebra_age_MD <- summary(trajectory_vertebra_age, show.trajectories = TRUE, attribute = "MD") 
trajectory_vertebra_age_MD
#Trajectory correlations -  are trajectories different in angle/direction?
trajectory_vertebra_age_TC <- summary(trajectory_vertebra_age, show.trajectories = TRUE, attribute = "TC", angle.type = "deg")
trajectory_vertebra_age_TC
#Trajectory shape differences - are trajectories different in shape?
trajectory_vertebra_age_SD <- summary(trajectory_vertebra_age, show.trajectories = TRUE, attribute = "SD") 
trajectory_vertebra_age_SD

#Save results to file
sink("Output/trajectory_vertebra_age.txt")
print("Magnitude difference (absolute difference between path distances) - length")
summary(trajectory_vertebra_age, show.trajectories = F, attribute = "MD") 
print("Correlations (angles) between trajectories - direction")
summary(trajectory_vertebra_age, show.trajectories = F, attribute = "TC", angle.type = "deg")
print("Shape differences between trajectory vectors - shape")
summary(trajectory_vertebra_age, show.trajectories = F, attribute = "SD") 
sink() 

#Plot results - PCA of fitted values
trajectory_vertebra_age_plot <- plot(trajectory_vertebra_age, main = "Trajectories shape change on vertebra by age",  pch = shapes, #title and type of point to be used
                                   col = c("darkgray","lightgray"), bg = c("darkgray","lightgray"), cex = 1, font.main = 2) #improve graphics
#Add line between groups
add.trajectories(trajectory_vertebra_age_plot, 
                 traj.pch = shapes, traj.col = 1, traj.lty = 1, traj.lwd = 1, traj.cex = 1.5, traj.bg = 1, 
                 start.bg = "green", end.bg = "red") #trajectory line graphics
#Add legend to see which trajectory belongs to each group
legend(list(x = -1.3, y = -0.6), legend = c("Thoracic" ,"Lumbar"  , "Caudal" ), pch =  shapes, pt.bg = 1, cex = 0.8)

##Make better PCA plot using ggplot
#Read PC scores as tibble
trajectory_vertebra_age_pcscores <- as_tibble(trajectory_vertebra_age_plot [["pc.points"]])

#Add group names and other attributes to tibble as vertebras
trajectory_vertebra_age_pcscores <- trajectory_vertebra_age_pcscores %>% mutate(vertebra = gdf_all$vertebra, age = gdf_all$age)
glimpse(trajectory_vertebra_age_pcscores)

#Calculate means of PC1 and PC2 at each vertebra per group to draw trajectories
trajectory_vertebra_age_pcscores_means <- trajectory_vertebra_age_pcscores %>% group_by(vertebra, age) %>%
  summarise_at(vars(PC1, PC2), list(mean = mean))              #function for means, both vertebras
glimpse(trajectory_vertebra_age_pcscores_means)

#Rename vertebras so they are easier to use for plot
trajectory_vertebra_age_pcscores_means <- trajectory_vertebra_age_pcscores_means %>% rename(x = PC1_mean, y = PC2_mean)
glimpse(trajectory_vertebra_age_pcscores_means)

#Nice plot
trajectory_vertebra_age_ggplot <- ggplot(trajectory_vertebra_age_pcscores, aes(x = PC1, y = PC2, shape = age, group = vertebra))+
  geom_point(size = 2.2, colour = "black", fill = "darkgray", alpha = 0.4, show.legend = F)+
  geom_point(data = trajectory_vertebra_age_pcscores_means, aes(x = x, y = y, fill = vertebra, shape = age, group = vertebra), colour = "black",
             size = 6, inherit.aes = F, alpha = 0.8)+
  geom_path(data = trajectory_vertebra_age_pcscores_means, aes(x = x, y = y, colour = vertebra, group = vertebra), inherit.aes = F, linewidth = 1,
            linejoin = 'mitre', show.legend = F)+
  scale_shape_manual(name = "Age", labels = c("Adult"  ,  "Juvenile", "Neonate"), values = shapes)+
  scale_colour_manual(name = "Vertebra", labels = c("Thoracic", "Lumbar", "Caudal"), #copy from as.factor(genera)
                      values = mypalette_vertebrae, aesthetics = c("colour", "fill"))+
  theme_bw()+
  xlab(paste0("PC 1 (",round(as.numeric(trajectory_vertebra_age[["pca"]][["sdev"]][1]^2/sum(trajectory_vertebra_age[["pca"]][["sdev"]]^2)*100), digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(as.numeric(trajectory_vertebra_age[["pca"]][["sdev"]][2]^2/sum(trajectory_vertebra_age[["pca"]][["sdev"]]^2)*100), digits = 2),"%)"))+
  theme(legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  guides(shapes = guide_legend(override.aes = list(colour = "black", fill = "darkgray", alpha = 1)),
         color = guide_legend(override.aes = list(shape = 23, colour = "black", fill = mypalette_vertebrae)))

#Visualize plot and save as PDF using menu in bar on the right
trajectory_vertebra_age_ggplot

###Sex ----

#Make sure to only use males and females!
rm_U_vertebra <- which(gdf_all$sex == "unknown")
  

gdf_all_2 <- geomorph.data.frame(pcscores = gdf_all$pcscores[-rm_U_vertebra,],
                                  sex = as.factor(as.vector(gdf_all$sex[-rm_U_vertebra])), 
                                 vertebra = gdf_all$vertebra[-rm_U_vertebra])
glimpse(gdf_all_2)

#Shape changes along vertebra at different ages
fit_shape_vertebra_sex <- lm.rrpp(pcscores ~ vertebra * sex, iter = 999, data = gdf_all_2, RRPP = F)

#Check that there is a significant correlation
summary(fit_shape_vertebra_sex)

#Save results to file
sink("Output/fit_shape_vertebra_sex.txt")
summary(fit_shape_vertebra_sex)
sink() 

#Use fit to calculate trajectories
trajectory_vertebra_sex <- trajectory.analysis(fit_shape_vertebra_sex, groups = gdf_all_2$vertebra, traj.pts = gdf_all_2$sex, 
                                               pca = TRUE) 

#View results
#Magnitude differences between trajectories, standard summary - are trajectories different in length?
trajectory_vertebra_sex_MD <- summary(trajectory_vertebra_sex, show.trajectories = TRUE, attribute = "MD") 
trajectory_vertebra_sex_MD
#Trajectory correlations -  are trajectories different in angle/direction?
trajectory_vertebra_sex_TC <- summary(trajectory_vertebra_sex, show.trajectories = TRUE, attribute = "TC", angle.type = "deg")
trajectory_vertebra_sex_TC
#Trajectory shape differences - are trajectories different in shape?
trajectory_vertebra_sex_SD <- summary(trajectory_vertebra_sex, show.trajectories = TRUE, attribute = "SD") 
trajectory_vertebra_sex_SD

#Save results to file
sink("Output/trajectory_vertebra_sex.txt")
print("Magnitude difference (absolute difference between path distances) - length")
summary(trajectory_vertebra_sex, show.trajectories = F, attribute = "MD") 
print("Correlations (angles) between trajectories - direction")
summary(trajectory_vertebra_sex, show.trajectories = F, attribute = "TC", angle.type = "deg")
print("Shape differences between trajectory vectors - shape")
summary(trajectory_vertebra_sex, show.trajectories = F, attribute = "SD") 
sink() 

#Plot results - PCA of fitted values
trajectory_vertebra_sex_plot <- plot(trajectory_vertebra_sex, main = "Trajectories shape change on vertebra by age", pch  = shapes,   #title and type of point to be used
                                     col = c("darkgray","lightgray"), bg = c("darkgray","lightgray"), cex = 1, font.main = 2) #improve graphics
#Add line between groups
add.trajectories(trajectory_vertebra_sex_plot, traj.pch = shapes,
                 traj.col = 1, traj.lty = 1, traj.lwd = 1, traj.cex = 1.5, traj.bg = 1, 
                 start.bg = "green", end.bg = "red") #trajectory line graphics
#Add legend to see which trajectory belongs to each group
legend(list(x = -1.4, y = -0.3), legend = c("Thoracic" ,"Lumbar"  , "Caudal" ), pch = shapes,  pt.bg = 1, cex = 0.8)

##Make better PCA plot using ggplot
#Read PC scores as tibble
trajectory_vertebra_sex_pcscores <- as_tibble(trajectory_vertebra_sex_plot [["pc.points"]])

#Add group names and other attributes to tibble as vertebras
trajectory_vertebra_sex_pcscores <- trajectory_vertebra_sex_pcscores %>% mutate(vertebra = gdf_all_2$vertebra, sex = gdf_all_2$sex)
glimpse(trajectory_vertebra_sex_pcscores)

#Calculate means of PC1 and PC2 at each vertebra per group to draw trajectories
trajectory_vertebra_sex_pcscores_means <- trajectory_vertebra_sex_pcscores %>% group_by(vertebra, sex) %>%
  summarise_at(vars(PC1, PC2), list(mean = mean))              #function for means, both vertebras
glimpse(trajectory_vertebra_sex_pcscores_means)

#Rename vertebras so they are easier to use for plot
trajectory_vertebra_sex_pcscores_means <- trajectory_vertebra_sex_pcscores_means %>% rename(x = PC1_mean, y = PC2_mean)
glimpse(trajectory_vertebra_sex_pcscores_means)

#Nice plot
trajectory_vertebra_sex_ggplot <- ggplot(trajectory_vertebra_sex_pcscores, aes(x = PC1, y = PC2, shape = sex, group = vertebra))+
  geom_point(size = 2.2, colour = "black", fill = "darkgray", alpha = 0.4, show.legend = F)+
  geom_point(data = trajectory_vertebra_sex_pcscores_means, aes(x = x, y = y, fill = vertebra, shape = sex, group = vertebra), colour = "black",
             size = 6, inherit.aes = F, alpha = 0.8)+
  geom_path(data = trajectory_vertebra_sex_pcscores_means, aes(x = x, y = y, colour = vertebra, group = vertebra), inherit.aes = F, linewidth = 1,
            linejoin = 'mitre', show.legend = F)+
  scale_shape_manual(name = "Sex", labels = c("Females", "Males"), values = shapes[c(1,2)])+
  scale_colour_manual(name = "Vertebra", labels = c("Thoracic", "Lumbar", "Caudal"), #copy from as.factor(genera)
                      values = mypalette_vertebrae, aesthetics = c("colour", "fill"))+
  theme_bw()+
  xlab(paste0("PC 1 (",round(as.numeric(trajectory_vertebra_sex[["pca"]][["sdev"]][1]^2/sum(trajectory_vertebra_sex[["pca"]][["sdev"]]^2)*100), digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(as.numeric(trajectory_vertebra_sex[["pca"]][["sdev"]][2]^2/sum(trajectory_vertebra_sex[["pca"]][["sdev"]]^2)*100), digits = 2),"%)"))+
  theme(legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  guides(shapes = guide_legend(override.aes = list(colour = "black", fill = "darkgray", alpha = 1)),
         color = guide_legend(override.aes = list(shape = 23, colour = "black", fill = mypalette_vertebrae)))

#Visualize plot and save as PDF using menu in bar on the right
trajectory_vertebra_sex_ggplot
