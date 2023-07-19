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
sink("allometry_models.txt")
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
sink("pairwise_allometry_vert.txt")
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
allometry_vert_plot <- data.frame(logCS = allometry_vert_plot_regscore[["plot.args"]][["x"]], 
                                   RegScores = allometry_vert_plot_regscore[["plot.args"]][["y"]])
#Convert data frame to tibble
allometry_vert_plot <- as_tibble(allometry_vert_plot)
#Add labels and other attributes to tibble as columns
allometry_vert_plot <- allometry_vert_plot %>% 
  mutate(specimens = gdf_all$specimen,  sex = gdf_all$sex, age = gdf_all$age, 
         vertebra = gdf_all$vertebra)
glimpse(allometry_vert_plot)

#Plot allometry regression by vertebra
allometry_vert_ggplot <- ggplot(allometry_vert_plot, aes(x = logCS, y = RegScores, colour = vertebra, fill = vertebra))+
  geom_point(size = 3, aes(shape = age), alpha = 0.3)+  
  geom_smooth(method = "lm", aes(fill = vertebra, colour = vertebra),           #confidence intervals and reg line, before points
              size = 1.0, alpha = 0.6, show.legend = F)+      #put col and other graphics OUTSIDE of aes()!!!
      #points after, so they are on top
  scale_color_manual(name = "Vertebra", labels = c( "Thoracic","Lumbar", "Caudal"),
                     values = mypalette_age_sex, aesthetics = c("color","fill"))+         
  scale_shape_manual(name = "Age", labels =  c("adult", "juvenile", "newborn"), 
                     values = shapes)+
  theme_classic(base_size = 12)+
  ylab("Regression Score")+
  theme(plot.title = element_text(face = "bold", hjust = 0.8, size = 11), legend.text = element_text(size = 12), 
        legend.title = element_text(size = 13), legend.position = "bottom", legend.direction = "horizontal")+
  guides(fill = guide_legend(label = F, title = NULL, override.aes = list(shape = NA)))
allometry_vert_ggplot

##Create residuals array from best model to then save as coordinates for analyses
allometry_residuals <- allometry_vert_int[["LM"]][["residuals"]]

#PCA ALLOMETRY RESIDUALS ----

##Make better PCA plot using ggplot
#Read PC scores as tibble
#Create tibble with new pc scores from allometry residuals and other concatenated variables
pcscores_res_vertebrae <- as_tibble(allometry_residuals)

#Add labels and other attributes to tibble as columns
pcscores_res_vertebrae <- pcscores_res_vertebrae %>% mutate(age = gdf_all$age, sex = gdf_all$sex, 
                                                            specimens = gdf_all$specimens, size = gdf_all$size ,
                                                            vertebra = gdf_all$vertebra)
glimpse(pcscores_res_vertebrae)

#PC1 and PC2 % values cannot be estimated so leave blank

#Make hulls for PCA plot with hulls around age
hulls_res_age_vertebrae <- pcscores_res_vertebrae %>%
  group_by(age) %>%
  slice(chull(comp1, comp2)) %>%
  rename(x = comp1, y = comp2)

#Nice PCA plot with hulls around age
PCA_res_vertebrae_age_ggplot <- ggplot(pcscores_res_vertebrae, aes(x = comp1, y = comp2, colour = age))+
  geom_point(size = 3, aes(shape = vertebra), fill = "white")+
  scale_colour_manual(name = "Age", labels =  c("adult" ,   "juvenile" ,"newborn" ), #to be ordered as they appear in tibble
                      values = mypalette_age_sex)+            #legend and color adjustments
  geom_polygon(data = hulls_res_age_vertebrae, aes(x = x, y = y, colour = age, fill = age), 
               alpha = .3, show.legend = FALSE, inherit.aes = F)+ #colored hulls with transparency
  scale_fill_manual(name = "Age", labels = c("adult"  ,  "juvenile", "newborn" ),
                    values =  mypalette_age_sex )+ #must match scale_colour_manual
  scale_shape_manual(name = "Vertebra", labels = c( "Thoracic","Lumbar", "Caudal"), #copy from as.factor(sex)
                     values = shapes)+
  theme_bw()+
  xlab("PC 1")+ #copy this from PC1_res object
  ylab("PC 2")+ #copy this from PC2_res object
  ggtitle("PCA res data vertebrae")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  #Remove legend for a scale_ using guide
  guides(fill = guide_legend(label = F, title = NULL, override.aes = list(shape = NA)))

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_res_vertebrae_age_ggplot

#Make hulls for PCA plot with hulls around sex
hulls_res_sex_vertebrae  <- pcscores_res_vertebrae  %>%
  group_by(sex) %>%
  slice(chull(comp1, comp2)) %>%
  rename(x = comp1, y = comp2)

#Nice PCA plot with hulls around sex
PCA_res_vertebrae_sex_ggplot <- ggplot(pcscores_res_vertebrae, aes(x = comp1, y = comp2, colour = sex))+
  geom_point(size = 3, aes(shape = vertebra), fill = "white")+
  scale_colour_manual(name = "Sex", labels =  c("F"   ,    "M"      , "unknown"), #to be ordered as they appear in tibble
                      values = mypalette_age_sex)+            #legend and color adjustments
  geom_polygon(data = hulls_res_sex_vertebrae, aes(x = x, y = y, colour = sex, fill = sex), 
               alpha = .3, show.legend = FALSE, inherit.aes = F)+ #colored hulls with transparency
  scale_fill_manual(name = "Sex", labels = c( "F"   ,    "M"      , "unknown"),
                    values =  mypalette_age_sex )+ #must match scale_colour_manual
  scale_shape_manual(name = "Vertebra", labels = c( "Thoracic","Lumbar", "Caudal"), #copy from as.factor(sex)
                     values = shapes)+
  theme_bw()+
  xlab("PC 1")+ #copy this from PC1_res object
  ylab("PC 2")+ #copy this from PC2_res object
  ggtitle("PCA res data vertebrae")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  #Remove legend for a scale_ using guide
  guides(fill = guide_legend(label = F, title = NULL, override.aes = list(shape = NA)))

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_res_vertebrae_sex_ggplot

#Make hulls for PCA plot with hulls around vertebra type
hulls_res_vertebra_vertebrae  <- pcscores_res_vertebrae  %>%
  group_by(vertebra) %>%
  slice(chull(comp1, comp2)) %>%
  rename(x = comp1, y = comp2)

#Nice PCA plot with hulls around sex
PCA_res_vertebrae_vertebra_ggplot <- ggplot(pcscores_res_vertebrae, aes(x = comp1, y = comp2, colour = vertebra))+
  geom_point(size = 3, aes(shape = age), fill = "white")+
  scale_colour_manual(name = "Vertebra", labels =  c("Thoracic","Lumbar", "Caudal"), #to be ordered as they appear in tibble
                      values = mypalette_age_sex )+            #legend and color adjustments
  geom_polygon(data = hulls_res_vertebra_vertebrae, aes(x = x, y = y, colour = vertebra, fill = vertebra), 
               alpha = .3, show.legend = FALSE, inherit.aes = F)+ #colored hulls with transparency
  scale_fill_manual(name = "Vertebra", labels = c("Thoracic","Lumbar", "Caudal"),
                    values =  mypalette_age_sex)+ #must match scale_colour_manual
  scale_shape_manual(name = "Age", labels = c( "adult"  ,  "juvenile", "newborn"), #copy from as.factor(age)
                     values = shapes)+
  theme_bw()+
  xlab("PC 1")+ #copy this from PC1_res object
  ylab("PC 2")+ #copy this from PC2_res object
  ggtitle("PCA all data vertebrae")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  #Remove legend for a scale_ using guide
  guides(fill = guide_legend(label = F, title = NULL, override.aes = list(shape = NA)))

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_res_vertebrae_vertebra_ggplot


#ANOVAs OF RESIDUALS PC SCORES AND SEX/AGE BY VERTEBRA TYPE CORRECTED BY SIZE ----
#Is there differences vertebra shape between sexes and ages after size is removed?
#Taking into account the shape difference based on column position
#Shape is expressed as residuals pc scores

#Calculate null model - changes in res based on vertebra type
res_vert_null <- lm.rrpp(allometry_residuals ~ pcscores_res_vertebrae$vertebra,  iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis null model
summary(res_vert_null)

#Save results to file
sink("Output/res_vertebra_model.txt")
print("res (pc scores) size corrected by vertebra type")
summary(res_vert_null) 
sink() 

##Sex ----
##Test different res between types of vertebrae ----
res_vert_sex_comb <- lm.rrpp(allometry_residuals ~ pcscores_res_vertebrae$vertebra + pcscores_res_vertebrae$sex,  iter=999, print.progress = TRUE) 
res_vert_sex_int <-  lm.rrpp(allometry_residuals ~ pcscores_res_vertebrae$vertebra * pcscores_res_vertebrae$sex,  iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of res with logCS
summary(res_vert_sex_comb)
summary(res_vert_sex_int) 

#Save results of significant regression to file
sink("Output/res_sex_models.txt")
print("Combination +")
summary(res_vert_sex_comb) 

print("Interaction *")
summary(res_vert_sex_int)
sink() 

#ANOVAs - is a model significantly better than the others?
anova_res_sex_models <- anova(res_vert_null, res_vert_sex_comb, res_vert_sex_int)
anova_res_sex_models

##Age ----
##Test different res between types of vertebrae ----
res_vert_age_comb <- lm.rrpp(allometry_residuals ~ pcscores_res_vertebrae$vertebra + pcscores_res_vertebrae$age,  iter=999, print.progress = TRUE) 
res_vert_age_int <-  lm.rrpp(allometry_residuals ~ pcscores_res_vertebrae$vertebra * pcscores_res_vertebrae$age,  iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of res with logCS
summary(res_vert_age_comb)
summary(res_vert_age_int) 

#Save results of significant regression to file
sink("Output/res_age_models.txt")
print("Combination +")
summary(res_vert_age_comb) 

print("Interaction *")
summary(res_vert_age_int)
sink() 

#ANOVAs - is a model significantly better than the others?
anova_res_age_models <- anova(res_vert_null, res_vert_age_comb, res_vert_age_int)
anova_res_age_models 

#Pairwise comparison for the combination and interaction model
#Helps determine if there is a significant difference in slope (int model) in the res_age trajectory on top of difference in intercept (comb model)
pairwise_res_age_vert <- pairwise(res_vert_age_int, fit.null = res_vert_null,
                                    groups = gdf_all$age, print.progress = FALSE) 
pairwise_res_age_vert

#Distances between slope vectors (end-points) - absolute difference between slopes of groups
#if significant means int model better than comb
pairwise_res_age_vert_dist <- summary(pairwise_res_age_vert, confidence = 0.95, test.type = "dist") 
pairwise_res_age_vert_dist

#Correlation between slope vectors (and angles) - similarity of vector orientation or angle,
#if significant means the vectors of the groups are oriented in different ways 
pairwise_res_age_vert_VC <- summary(pairwise_res_age_vert, confidence = 0.95, test.type = "VC",
                                      angle.type = "deg") 
pairwise_res_age_vert_VC

#Absolute difference between slope vector lengths - difference in rate of change per covariate unit (size),
#if significant means there is a significant rate of change difference in res between groups during growth
pairwise_res_age_vert_DL <-summary(pairwise_res_age_vert, confidence = 0.95, test.type = "DL") 
pairwise_res_age_vert_DL

#Save results to file
sink("Output/pairwise_res_age_vert.txt")
print("ANOVA models")
print(anova_res_age_models)

print("1-Pairwise absolute distances slopes")
summary(pairwise_res_age_vert, confidence = 0.95, test.type = "dist") 

print("2-Distance between angles (slope directions)")
summary(pairwise_res_age_vert, confidence = 0.95, test.type = "VC",
        angle.type = "deg") 

print("3-Difference in slope vector length (difference in rate of change of res per unit of size)")
summary(pairwise_res_age_vert, confidence = 0.95, test.type = "DL") 

sink()
