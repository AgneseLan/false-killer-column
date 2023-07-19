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

##Shape changes along column ----
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
summary(trajectory_column_age, show.trajectories = TRUE, attribute = "MD") 
print("Correlations (angles) between trajectories - direction")
summary(trajectory_column_age, show.trajectories = TRUE, attribute = "TC", angle.type = "deg")
print("Shape differences between trajectory vectors - shape")
summary(trajectory_column_age, show.trajectories = TRUE, attribute = "SD") 
sink() 

#Plot results - PCA of fitted values
trajectory_column_age_plot <- plot(trajectory_column_age, main = "Trajectories shape change on column by age",  pch = shapes, #title and type of point to be used
                              col = c("darkgray","lightgray"), bg = c("darkgray","lightgray"), cex = 1, font.main = 2) #improve graphics
#Add line between groups
add.trajectories(trajectory_column_age_plot, 
                 traj.pch = shapes, traj.col = 1, traj.lty = 1, traj.lwd = 1, traj.cex = 1.5, traj.bg = 1, 
                 start.bg = "green", end.bg = "red") #trajectory line graphics
#Add legend to see which trajectory belongs to each group
legend("bottomleft", legend = c("adult"  ,  "juvenile", "newborn"), pch =  shapes, pt.bg = 0.5, cex = 0.5)

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

######################################################################################

#Trajectory AGE Nice plot (color blind friendly)

mycols_age <- colors()[c("#DDCC77", "#6699CC", "#AA4499")] 


mycols_age= carto_pal(12, "Safe")
mycols_age

trajectory_column_age_ggplot <- ggplot(trajectory_column_age_pcscores, aes(x = PC1, y = PC2, shape = vertebra, group = age))+
  geom_point(size = 2.2, colour = "darkgray", fill = "darkgray", alpha = 0.5)+
  geom_point(data = trajectory_column_age_pcscores_means, aes(x = x, y = y, colour = age, shape = vertebra, group = age), fill = "white",
            size = 2.75, stroke = 2, inherit.aes = F)+
  geom_path(data = trajectory_column_age_pcscores_means, aes(x = x, y = y, colour = age, group = age), inherit.aes = F, size = 1,
           linejoin = 'mitre', show.legend = F)+
  scale_shape_manual(name = "Vertebra", labels = c("Thoracic", "Lumbar", "Caudal"), values = shapes)+
  scale_colour_manual(name = "Age", labels = c("adult"  ,  "juvenile", "neonate"), #copy from as.factor(genera)
                      values = c("#F0E442", "#6699CC", "#AA4499"), aesthetics = c("colour", "fill"))+
  theme_bw()+
  xlab("PC 1 (63.78%)")+ #copy this from standard trajectory plot
  ylab("PC 2 (20.27%)")+
  theme(legend.position = "right", legend.direction = "vertical", legend.justification = "center")

#Visualize plot and save as PDF using menu in bar on the right
trajectory_column_age_ggplot

##############################################################################

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
sink("trajectory_column_sex.txt")
print("Magnitude difference (absolute difference between path distances) - length")
summary(trajectory_column_sex, show.trajectories = TRUE, attribute = "MD") 
print("Correlations (angles) between trajectories - direction")
summary(trajectory_column_sex, show.trajectories = TRUE, attribute = "TC", angle.type = "deg")
print("Shape differences between trajectory vectors - shape")
summary(trajectory_column_sex, show.trajectories = TRUE, attribute = "SD") 
sink() 

#Plot results - PCA of fitted values
trajectory_column_sex_plot <- plot(trajectory_column_sex, main = "Trajectories shape change on column by sex",  pch = shapes, #title and type of point to be used
                                   col = c("darkgray","lightgray"), bg = c("darkgray","lightgray"), cex = 1, font.main = 2) #improve graphics
#Add line between groups
add.trajectories(trajectory_column_sex_plot, 
                 traj.pch = shapes, traj.col = 1, traj.lty = 1, traj.lwd = 1, traj.cex = 1.5, traj.bg = 1, 
                 start.bg = "green", end.bg = "red") #trajectory line graphics
#Add legend to see which trajectory belongs to each group
legend("bottomleft", legend = c( "F"   ,    "M"      , "unknown"), pch =  shapes, pt.bg = 0.5, cex = 0.5)

###
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

## nice plot
trajectory_column_sex_ggplot <- ggplot(trajectory_column_sex_pcscores, aes(x = PC1, y = PC2, shape = sex, group = sex))+
  geom_point(size = 3, colour = "darkgray", fill = "darkgray", alpha = 0.5)+
  geom_point(data = trajectory_column_sex_pcscores_means, aes(x = x, y = y, colour = sex, shape = sex, group = sex), fill = "white",
             size = 4, stroke = 2, inherit.aes = F)+
  geom_path(data = trajectory_column_sex_pcscores_means, aes(x = x, y = y, colour = sex, group = sex), inherit.aes = F, size = 1,
            linejoin = 'mitre', show.legend = F)+
  scale_shape_manual(name = "sex", labels = c("female"  ,  "male", "neonate"), values = shapes)+
  scale_colour_manual(name = "sex", labels = c("female"  ,  "male", "neonate"), #copy from as.factor(genera)
                      values = mypalette_age_sex, aesthetics = c("colour", "fill"))+
  theme_bw()+
  xlab("PC 1 (63.78%)")+ #copy this from standard trajectory plot
  ylab("PC 2 (20.27%)")+
  ggtitle("Trajectories shape change on column by sex")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14), 
        legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center")

#################################################################################################################
#Fig. 5

#Palette to use for sex PCAs
mypalette_age_sex <- c(mypalette_paired[5], mypalette_paired[2], mypalette_paired[10])
image(1:3, 1, as.matrix(1:3), col = mypalette_age_sex,
      ylab = "", yaxt = "n")

#create shape palette for sex
mycols_age= carto_pal(7, "Purp")
mycols_age

mycols_age= carto_pal(7, "BluGrn")
mycols_age

#colorblind code   "#117733" "#332288"  "#44AA99" "#999933" "#882255" "#661100"  
#[!master_df$GenderDescription %in% c("Unknown"),]
shapes <- c(19, 15, 17)

trajectory_column_sex_ggplot <- ggplot(trajectory_column_sex_pcscores, aes(x = PC1, y = PC2, shape = vertebra, group = sex))+
  geom_point(size = 2.75, colour = "gray", fill = "gray", alpha = 0.5)+
  geom_point(data = trajectory_column_sex_pcscores_means, aes(x = x, y = y, colour = sex, shape = vertebra, group = sex), fill = "white",
             size = 2, stroke = 2, inherit.aes = F)+
  geom_path(data = trajectory_column_sex_pcscores_means, aes(x = x, y = y, colour = sex, group = sex), inherit.aes = F, size = 1,
            linejoin = 'mitre', show.legend = F)+
  scale_shape_manual(name = "vertebra", labels = c("thoracic"  ,  "lumbar", "caudal"), values = shapes)+
  scale_colour_manual(name = "sex", labels = c("Female"  ,  "Male", "Unknown (neonate)"), #copy from as.factor(genera)
                      values = c("#E597B9", "#44AA99", "#826DBA"), aesthetics = c("colour", "fill"))+
  theme_bw()+
  xlab("PC 1 (63.78%)")+ #copy this from standard trajectory plot
  ylab("PC 2 (20.27%)")+
  theme(legend.position = "right", legend.direction = "vertical", legend.justification = "center")

trajectory_column_sex_ggplot

######################################
####remove neonates
###ggplot(data=master_df[!master_df$GenderDescription %in% c("Unknown"),] 

trajectory_column_sex_ggplot <- ggplot(trajectory_column_sex_pcscores[!trajectory_column_sex_pcscores$sex%in% c("unknown"),], aes(x = PC1, y = PC2, shape = vertebra, group = sex))+
  geom_point(size = 2.75, colour = "gray", fill = "gray", alpha = 0.5)+
  geom_point(data = trajectory_column_sex_pcscores_means[!trajectory_column_sex_pcscores_means$sex%in% c("unknown"),], aes(x = x, y = y, colour = sex, shape = vertebra, group = sex), fill = "white",
             size = 2, stroke = 2, inherit.aes = F)+
  geom_path(data = trajectory_column_sex_pcscores_means[!trajectory_column_sex_pcscores_means$sex%in% c("unknown"),], aes(x = x, y = y, colour = sex, group = sex), inherit.aes = F, size = 1,
            linejoin = 'mitre', show.legend = F)+
  scale_shape_manual(name = "vertebra", labels = c("thoracic"  ,  "lumbar", "caudal"), values = shapes)+
  scale_colour_manual(name = "sex", labels = c("Female"  ,  "Male"), #copy from as.factor(genera)
                      values = c("#E597B9", "#332288"), aesthetics = c("colour", "fill"))+
  theme_bw()+
  xlab("PC 1 (63.78%)")+ #copy this from standard trajectory plot
  ylab("PC 2 (20.27%)")+
  theme(legend.position = "right", legend.direction = "vertical", legend.justification = "center")

trajectory_column_sex_ggplot

###################################
#including the neonates for supplementary materials
trajectory_column_sex_ggplot <- ggplot(trajectory_column_sex_pcscores, aes(x = PC1, y = PC2, shape = vertebra, group = sex))+
  geom_point(size = 2.75, colour = "gray", fill = "gray", alpha = 0.5)+
  geom_point(data = trajectory_column_sex_pcscores_means, aes(x = x, y = y, colour = sex, shape = vertebra, group = sex), fill = "white",
             size = 2, stroke = 2, inherit.aes = F)+
  geom_path(data = trajectory_column_sex_pcscores_means, aes(x = x, y = y, colour = sex, group = sex), inherit.aes = F, size = 1,
            linejoin = 'mitre', show.legend = F)+
  scale_shape_manual(name = "vertebra", labels = c("thoracic"  ,  "lumbar", "caudal"), values = shapes)+
  scale_colour_manual(name = "sex", labels = c("Female"  ,  "Male", "Neonates (unknown)"), #copy from as.factor(genera)
                      values = c("#E597B9", "#332288","#44AA99"), aesthetics = c("colour", "fill"))+
  theme_bw()+
  xlab("PC 1 (63.78%)")+ #copy this from standard trajectory plot
  ylab("PC 2 (20.27%)")+
  theme(legend.position = "right", legend.direction = "vertical", legend.justification = "center")

trajectory_column_sex_ggplot

############################################################################################################################"""""""""""""""
##Create project palettes----
mypalette_paired <- brewer.pal(12,"Paired")
image(1:12, 1, as.matrix(1:12), col = mypalette_paired, xlab = "Paired",
      ylab = "", yaxt = "n")

mypalette_blue <- as.matrix(ggthemes_data[["tableau"]][["color-palettes"]][["ordered-sequential"]][["Blue"]][["value"]])
image(1:20, 1, as.matrix(1:20), col = mypalette_blue, xlab = "Blue",
      ylab = "", yaxt = "n")

#Palette for age and sex - early, late/new, immature, adult
mypalette_age_sex <- c(mypalette_paired[5], mypalette_paired[2], mypalette_paired[10])
image(1:3, 1, as.matrix(1:3), col = mypalette_age_sex,
      ylab = "", yaxt = "n")

#Create shape palette 3 sex/age
shapes <- c(21, 22, 24)
shapes <- c(21, 22, 23)
##Shape changes per vertebra type ----
#First perform lm.rrpp to create linear model that describes what we are trying to test

###Age ----
#Shape changes along vertebra at different ages
fit_shape_vertebra_age <- lm.rrpp(pcscores ~ vertebra * age, iter = 999, data = gdf_all, RRPP = F)

#Check that there is a significant correlation
summary(fit_shape_vertebra_age)

#Save results to file
sink("fit_shape_vertebra_age.txt")
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
sink("trajectory_vertebra_age.txt")
print("Magnitude difference (absolute difference between path distances) - length")
summary(trajectory_vertebra_age, show.trajectories = TRUE, attribute = "MD") 
print("Correlations (angles) between trajectories - direction")
summary(trajectory_vertebra_age, show.trajectories = TRUE, attribute = "TC", angle.type = "deg")
print("Shape differences between trajectory vectors - shape")
summary(trajectory_vertebra_age, show.trajectories = TRUE, attribute = "SD") 
sink() 

#Plot results - PCA of fitted values
trajectory_vertebra_age_plot <- plot(trajectory_vertebra_age, main = "Trajectories shape change on vertebra by age",  pch = shapes, #title and type of point to be used
                                   col = c("darkgray","lightgray"), bg = c("darkgray","lightgray"), cex = 1, font.main = 2) #improve graphics
#Add line between groups
add.trajectories(trajectory_vertebra_age_plot, 
                 traj.pch = shapes, traj.col = 1, traj.lty = 1, traj.lwd = 1, traj.cex = 1.5, traj.bg = 1, 
                 start.bg = "green", end.bg = "red") #trajectory line graphics
#Add legend to see which trajectory belongs to each group
legend("bottomleft", legend = c("Thoracic" ,"Lumbar"  , "Caudal" ), pch =  shapes, pt.bg = 0.5, cex = 0.5)

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
  geom_point(size = 3, colour = "darkgray", fill = "darkgray", alpha = 0.5)+
  geom_point(data = trajectory_vertebra_age_pcscores_means, aes(x = x, y = y, colour = vertebra, shape = age, group = vertebra), fill = "white",
             size = 4, stroke = 2.2, inherit.aes = F)+
  geom_path(data = trajectory_vertebra_age_pcscores_means, aes(x = x, y = y, colour = vertebra, group = vertebra), inherit.aes = F, size = 1.2,
            linejoin = 'mitre', show.legend = F)+
  scale_shape_manual(name = "Age", labels = c("adult"  ,  "juvenile", "neonate"), values = shapes)+
  scale_colour_manual(name = "Vertebra", labels = c("Thoracic", "Lumbar", "Caudal"), 
                      values = mypalette_age_sex, aesthetics = c("colour", "fill"))+
  theme_bw()+
  xlab("PC 1 (63.78%)")+ #copy this from standard trajectory plot
  ylab("PC 2 (20.27%)")


#Visualize plot and save as PDF using menu in bar on the right
trajectory_vertebra_age_ggplot

###Sex ----
#Shape changes along vertebra at different sex
fit_shape_vertebra_sex <- lm.rrpp(pcscores ~ vertebra * sex, iter = 999, data = gdf_all, RRPP = F)

#Check that there is a significant correlation
summary(fit_shape_vertebra_sex)

#Save results to file
sink("Output/fit_shape_vertebra_sex.txt")
summary(fit_shape_vertebra_sex)
sink() 

#Use fit to calculate trajectories
trajectory_vertebra_sex <- trajectory.analysis(fit_shape_vertebra_sex, groups = gdf_all$vertebra, traj.pts = gdf_all$sex, 
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
summary(trajectory_vertebra_sex, show.trajectories = TRUE, attribute = "MD") 
print("Correlations (angles) between trajectories - direction")
summary(trajectory_vertebra_sex, show.trajectories = TRUE, attribute = "TC", angle.type = "deg")
print("Shape differences between trajectory vectors - shape")
summary(trajectory_vertebra_sex, show.trajectories = TRUE, attribute = "SD") 
sink() 

#Plot results - PCA of fitted values
trajectory_vertebra_sex_plot <- plot(trajectory_vertebra_sex, main = "Trajectories shape change on vertebra by sex",  pch = shapes, #title and type of point to be used
                                   col = c("darkgray","lightgray"), bg = c("darkgray","lightgray"), cex = 1, font.main = 2) #improve graphics
#Add line between groups
add.trajectories(trajectory_vertebra_sex_plot, 
                 traj.pch = shapes, traj.col = 1, traj.lty = 1, traj.lwd = 1, traj.cex = 1.5, traj.bg = 1, 
                 start.bg = "green", end.bg = "red") #trajectory line graphics
#Add legend to see which trajectory belongs to each group
legend("bottomleft", legend = c( "Thoracic"   ,    "Lumbar"      , "Caudal"), pch =  shapes, pt.bg = 0.5, cex = 0.5)
