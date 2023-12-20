#CH.4 - Allometry analyses morphoblocks PC scores

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
library(ggfortify)
library(RColorBrewer) 
library(ggthemes)
library(ggpubr)
library(ggplotify)
library(Morpho)
library(rphylopic)
library(png)
library(gridExtra)
library(reshape2)
library(scales)


#ALLOMETRY CORRECTION ----
##Evaluate allometry and get the allometry-free shapes using LogCS, use this for analyses

#Regression pc scores on logCS size
allometry <- lm.rrpp(gdf_all$pcscores ~ gdf_all$size, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of allometry with logCS
summary(allometry) 

#Regression score of shape vs logCS - regression method with "RegScore" plotting
allometry_plot_regscore <- plot(allometry, type = "regression",predictor = gdf_all$size , reg.type = "RegScore",
                                main = "PCs vs logCS",xlab = "logCS", pch = 21, col = "chartreuse4", bg = "chartreuse4", cex = 1.2, font.main = 2)   #improve graphics
text(x = gdf_all$size , y = allometry_plot_regscore$RegScore, labels = gdf_all$specimen,
     pos = 3, offset = 0.5, cex = 0.75)    #improve appearance of labels

##Test different allometry between types of vertebrae ----
allometry_vert_comb <-  lm.rrpp(gdf_all$pcscores ~ gdf_all$size  + gdf_all$vertebra, iter=999, print.progress = TRUE) 
allometry_vert_int <-  lm.rrpp(gdf_all$pcscores ~ gdf_all$size  * gdf_all$vertebra, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of allometry with logCS
summary(allometry_vert_comb)
summary(allometry_vert_int) 

#Save results of significant regression to file
sink("Output/allometry_models.txt")
print("Null")
summary(allometry)

print("Combination +")
summary(allometry_vert_comb) 

print("Interaction *")
summary(allometry_vert_int)
sink() 

#ANOVAs - is a model significantly better than the others?
anova_allometry_models <- anova(allometry, allometry_vert_comb, allometry_vert_int)
anova_allometry_models

#Pairwise comparison for the combination and interaction model
#Helps determine if there is a significant difference in slope (int model) in the allometry trajectory on top of difference in intercept (comb model)
pairwise_allometry_vert <- pairwise(allometry_vert_int, fit.null = allometry_vert_comb,
                                                   groups = gdf_all$vertebra, 
                                                   covariate =  gdf_all$size , print.progress = FALSE) 
pairwise_allometry_vert

#Distances between slope vectors (end-points) - absolute difference between slopes of groups
#if significant means int model better than comb
pairwise_allometry_vert_dist <- summary(pairwise_allometry_vert, confidence = 0.95, test.type = "dist") 
pairwise_allometry_vert_dist

#Correlation between slope vectors (and angles) - similarity of vector orientation or angle,
#if significant means the vectors of the groups are oriented in different ways 
pairwise_allometry_vert_VC <- summary(pairwise_allometry_vert, confidence = 0.95, test.type = "VC",
        angle.type = "deg") 
pairwise_allometry_vert_VC

#Absolute difference between slope vector lengths - difference in rate of change per covariate unit (size),
#if significant means there is a significant rate of change difference in shape between groups during growth
pairwise_allometry_vert_DL <-summary(pairwise_allometry_vert, confidence = 0.95, test.type = "DL") 
pairwise_allometry_vert_DL

#Save results to file
sink("Output/pairwise_allometry_vert.txt")
print("ANOVA models")
print(anova_allometry_models)

print("1-Pairwise absolute distances slopes")
summary(pairwise_allometry_vert, confidence = 0.95, test.type = "dist") 

print("2-Distance between angles (slope directions)")
summary(pairwise_allometry_vert, confidence = 0.95, test.type = "VC",
        angle.type = "deg") 

print("3-Difference in slope vector length (difference in rate of change of shape per unit of size)")
summary(pairwise_allometry_vert, confidence = 0.95, test.type = "DL") 

sink()

###Plot regression by vertebra type ----
#Regression score of shape vs logCS and comb or int (best model)- regression method with "RegScore" plotting
allometry_vert_plot_regscore <- plot(allometry_vert_int, type = "regression",predictor = gdf_all$size , reg.type = "RegScore",
                                     main = "PCs vs logCS",xlab = "logCS", pch = 21, col = "chartreuse4", bg = "chartreuse4", cex = 1.2, font.main = 2)   #improve graphics
text(x = gdf_all$size , y = allometry_vert_plot_regscore$RegScore, labels = gdf_all$specimen,
     pos = 3, offset = 0.5, cex = 0.75)    #improve appearance of labels

##Make better allometry plot with ggplot
#Create data frame object that ggplot can read - use data from plot object you want to improve
allometry_vert_plot <- data.frame(logCS = allometry_vert_plot_regscore[["plot_args"]][["x"]], 
                                   RegScores = allometry_vert_plot_regscore[["plot_args"]][["y"]])
#Convert data frame to tibble
allometry_vert_plot <- as_tibble(allometry_vert_plot)
#Add labels and other attributes to tibble as columns
allometry_vert_plot <- allometry_vert_plot %>% 
  mutate(specimens = gdf_all$specimen,  sex = gdf_all$sex, age = gdf_all$age, 
         vertebra = gdf_all$vertebra)
glimpse(allometry_vert_plot)

#Plot allometry regression by vertebra
allometry_vert_ggplot <- ggplot(allometry_vert_plot, aes(x = logCS, y = RegScores))+
  geom_point(size = 3, aes(shape = age, fill = vertebra), alpha = 0.2, color = "black")+  
  geom_smooth(method = "lm", aes(fill = vertebra, colour = vertebra, linetype = vertebra),           #confidence intervals and reg line, before points
              linewidth = 1.0, alpha = 0.5, show.legend = F)+      #put col and other graphics OUTSIDE of aes()!!!
      #points after, so they are on top
  scale_color_manual(name = "Vertebra", labels = c( "Thoracic","Lumbar", "Caudal"),
                     values = mypalette_vertebrae, aesthetics = c("color","fill"))+         
  scale_shape_manual(name = "Age", labels =  c("Adult", "Juvenile", "Neonate"), 
                     values = shapes)+
  scale_linetype_manual(name = "Vertebra", labels = c( "Thoracic","Lumbar", "Caudal"),
                        values = c(1,2,3))+
  theme_classic(base_size = 14)+
  xlab("Size (logCS)")+
  ylab("Regression Score")+
  theme(legend.key = element_blank(), legend.title = element_text(size = 15, face = "bold"), legend.text = element_text(size = 14),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  guides(color = guide_legend(override.aes = list(shape = 23, colour = "black", fill = mypalette_vertebrae, alpha =1)),
         shape = guide_legend(override.aes = list(colour = "black", fill = "gray50", alpha =1)))
allometry_vert_ggplot
#Add lines manually to legend


