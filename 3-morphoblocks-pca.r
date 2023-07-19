#CH.3 - Run morphoblocks PCA and ANOVA sex/age

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
library(gridExtra)
library(reshape2)
library(scales)
library(morphoBlocks)

#devtools::install_github("aharmer/morphoBlocks")

#BLOCKS PC SCORES ----

#Add column with vertebra type, pc scores and raw coordinates
gdf_thoracic$vertebra <- as.factor(rep("Thoracic", times = length(gdf_thoracic$Id)))
gdf_thoracic$
gdf_thoracic$raw <- final_dataset_thoracic

#do this for each type of vertebra
gdf_lumbar$vertebra <- as.factor(rep("Lumbar", times = length(gdf_lumbar$Id)))

gdf_lumbar$raw <- final_dataset_lumbar

gdf_caudal$vertebra <- as.factor(rep("Caudal", times = length(gdf_caudal$Id)))

gdf_caudal$raw <- final_dataset_caudal

#Remove extra vertebrae to have 3 per type per specimen
#Same number of elements per block per specimen needed for morphoblock to work
lumbar_delete <-  seq(4, length(gdf_lumbar$specimens), by = 4)
caudal_delete <- c(seq(4, length(gdf_caudal$specimens), by = 4), seq(5, length(gdf_caudal$specimens), by = 5))

#Raw data with extra vertebrae removed
raw_thoracic <- final_dataset_thoracic
raw_lumbar <- final_dataset_lumbar[,,-lumbar_delete]
raw_caudal <- final_dataset_caudal[,,-caudal_delete]

#MORPHOBLOCKS#
#Create coordinate blocks
block_thoracic <- formatBlock(raw_thoracic, k = 3, gpa = TRUE)
block_lumbar <- formatBlock(raw_lumbar, k = 3, gpa = TRUE)
block_caudal <- formatBlock(raw_caudal, k = 3, gpa = TRUE)

#Create block logCS
block_thoracic_logCS <- log(block_thoracic@centroid)
block_lumbar_logCS <- log(block_lumbar@centroid)
block_caudal_logCS <- log(block_caudal@centroid)

#Create logCS for consensus shape
mean_logCS <- data.frame(t = block_thoracic_logCS, l = block_lumbar_logCS, c = block_caudal_logCS)
mean_logCS <- mean_logCS %>% mutate(consensus = rowMeans(mean_logCS))
block_consensus_logCS <- mean_logCS[,4]
  
#List blocks for analysis
blocklist <- combineBlocks(blocks = c(block_thoracic, block_lumbar, block_caudal))

#Perform adjusted PCA on block list
#See analyseBlocks vignette for standard PCA option
pca_all <- analyseBlocks(blocklist, option = "rcpca", ncomp = 95)
#Check it worked
scoresPlot(pca_all)
par(mfrow=c(1,1))

#Make gdf with all data needed for analysis
#Make sure each variable is added in the same order (if thoracic first in pca first in all other etc.)

gdf_all <- geomorph.data.frame(pcscores = rbind(pca_all[["scores"]][[1]], pca_all[["scores"]][[2]],pca_all[["scores"]][[3]]), 
                              sex = as.factor(c(gdf_thoracic$sex, gdf_lumbar$sex[-lumbar_delete], gdf_caudal$sex[-caudal_delete])),
                               specimen = as.factor(c(gdf_thoracic$specimens, gdf_lumbar$specimens[-lumbar_delete], gdf_caudal$specimens[-caudal_delete])), 
                               age = as.factor(c(gdf_thoracic$age, gdf_lumbar$age[-lumbar_delete], gdf_caudal$age[-caudal_delete])), 
                               size = c(block_thoracic_logCS, block_lumbar_logCS, block_caudal_logCS), 
                              vertebra = as.factor(c(gdf_thoracic$vertebra, gdf_lumbar$vertebra[-lumbar_delete], gdf_caudal$vertebra[-caudal_delete])))

glimpse(gdf_all)

consensus <- rep(levels(gdf_all$vertebra), times = length(levels(gdf_thoracic$specimens)))

gdf_consensus <- geomorph.data.frame(pcscores = pca_all[["scores"]][[4]], 
                               sex = as.factor(gdf_thoracic$sex),
                               specimen = as.factor(gdf_thoracic$specimens), 
                               age = as.factor(gdf_thoracic$age), 
                               size = block_consensus_logCS, 
                               vertebra = consensus)

glimpse(gdf_consensus)

#PCA ALL vertebrae ----

#Create tibble with new pc scores from gdf_all$pcscores and other concatenated variables
pcscores_all_vertebrae <- as_tibble(gdf_all$pcscores)

#Add labels and other attributes to tibble as columns
pcscores_all_vertebrae <- pcscores_all_vertebrae %>% mutate(age = gdf_all$age, sex = gdf_all$sex, 
                                        specimens = gdf_all$specimen, size = gdf_all$size,
                                       vertebra = gdf_all$vertebra)
glimpse(pcscores_all_vertebrae)

#PC1 and PC2 % values can be estimated as the average of each type of specimen

#PC1
#Calculate mean of all 3 values
PC1_all <- mean(pca_all[["result"]][["AVE"]][["AVE_X"]][[1]][["comp1"]], pca_all[["result"]][["AVE"]][["AVE_X"]][[2]][["comp1"]],
                     pca_all[["result"]][["AVE"]][["AVE_X"]][[3]][["comp1"]])*100
                     
#PC2
#Calculate mean of all 3 values
PC2_all <- mean(pca_all[["result"]][["AVE"]][["AVE_X"]][[1]][["comp2"]], pca_all[["result"]][["AVE"]][["AVE_X"]][[2]][["comp2"]],
                pca_all[["result"]][["AVE"]][["AVE_X"]][[3]][["comp2"]])*100

#Plot
PCA_all_vertebrae_ggplot <- ggplot(pcscores_all_vertebrae, aes(x = comp1, y = comp2, label = specimens, colour = age, fill = age))+
  geom_point(size = 3, aes(shape = sex))+
  geom_text_repel(colour = "black", size = 4, max.overlaps = 40)+
  scale_colour_manual(name = "Age", labels =  c("adult"  ,  "juvenile", "newborn" ), #to be ordered as they appear in tibble
                      values = mypalette_age_sex , aesthetics = c("colour","fill"))+            #legend and color adjustments
  scale_shape_manual(name = "Sex", labels = c( "F"   ,    "M"      , "unknown"), #copy from as.factor(sex)
                     values = shapes)+
  theme_bw()+
  xlab("PC 1 (28.73%)")+ #copy this from PC1_all object
  ylab("PC 2 (26.97%)")+ #copy this from PC2_all object
  ggtitle("PCA all data vertebrae")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_vertebrae_ggplot

#Make hulls for PCA plot with hulls around age
hulls_all_age_vertebrae <- pcscores_all_vertebrae %>%
  group_by(age) %>%
  slice(chull(comp1, comp2)) %>%
  rename(x = comp1, y = comp2)


##############################################

####colours
# We can create choose a palette based on the R chart as follow:
mycols_age <- colors()[c(142, 616, 548)] 


mycols_age <- colors()[c("#DDCC77", "#6699CC", "#AA4499")] 

############# color blind friendly

mycols_age= carto_pal(12, "Safe")
mycols_age


#Nice PCA plot with hulls around age COLOR BLIND FRIENDLY

PCA_all_vertebrae_age_ggplot <- ggplot(pcscores_all_vertebrae, aes(x = comp1, y = comp2, colour = age))+
  geom_point(size = 3, aes(shape = vertebra), fill = "white")+
  scale_colour_manual(name = "Age", labels =  c("adult" ,   "juvenile" ,"neonate" ), #to be ordered as they appear in tibble
                      values = c("#F0E442", "#6699CC", "#AA4499"))+            #legend and color adjustments
  geom_polygon(data = hulls_all_age_vertebrae, aes(x = x, y = y, colour = age, fill = age), 
               alpha = .3, show.legend = FALSE, inherit.aes = F)+ #colored hulls with transparency
  scale_fill_manual(name = "Age", labels = c("adult"  ,  "juvenile", "neonate" ),
                    values =  c("#F0E442", "#6699CC", "#AA4499") )+ #must match scale_colour_manual
  scale_shape_manual(name = "Vertebra", labels = c("Thoracic","Lumbar", "Caudal"), #copy from as.factor(sex)
                     values = shapes)+
  theme_bw()+
  xlab("PC 1 (28.73%)")+ #copy this from PC1_all object
  ylab("PC 2 (26.97%)")+ #copy this from PC2_all object
 
  #Remove legend for a scale_ using guide
  guides(fill = guide_legend(label = F, title = NULL, override.aes = list(shape = NA)))

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_vertebrae_age_ggplot


#Make hulls for PCA plot with hulls around sex
hulls_all_sex_vertebrae  <- pcscores_all_vertebrae  %>%
  group_by(sex) %>%
  slice(chull(comp1, comp2)) %>%
  rename(x = comp1, y = comp2)

#Nice PCA plot with hulls around sex
PCA_all_vertebrae_sex_ggplot <- ggplot(pcscores_all_vertebrae[!pcscores_all_vertebrae$sex%in% c("unknown"),], aes(x = comp1, y = comp2, colour = sex))+
  geom_point(size = 3, aes(shape = vertebra), fill = "white")+
  scale_colour_manual(name = "Sex", labels =  c("Female"   ,    "Male"      , "Unknown (neonate)"), #to be ordered as they appear in tibble
                      values = c("#E597B9", "#332288"))+            #legend and color adjustments
  geom_polygon(data = hulls_all_sex_vertebrae[!hulls_all_sex_vertebrae$sex%in% c("unknown"),], aes(x = x, y = y, colour = sex, fill = sex), 
               alpha = .3, show.legend = FALSE, inherit.aes = F)+ #colored hulls with transparency
  scale_fill_manual(name = "Sex", labels = c( "Female"   ,    "Male"      , "Unknown (neonate)"),
                    values = c("#E597B9", "#332288"))+ #must match scale_colour_manual
  scale_shape_manual(name = "Vertebra", labels = c("Thoracic","Lumbar", "Caudal"), #copy from as.factor(sex)
                     values = shapes)+
  theme_bw()+
  xlab("PC 1 (28.73%)")+ #copy this from PC1_all object
  ylab("PC 2 (26.97%)")+ #copy this from PC2_all object
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  #Remove legend for a scale_ using guide
  guides(fill = guide_legend(label = F, title = NULL, override.aes = list(shape = NA)))

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_vertebrae_sex_ggplot

##################################################

PCA_all_vertebrae_sex_ggplot <- ggplot(pcscores_all_vertebrae, aes(x = comp1, y = comp2, colour = sex))+
  geom_point(size = 3, aes(shape = vertebra), fill = "white")+
  scale_colour_manual(name = "Sex", labels =  c("Female"   ,    "Male"      , "Unknown (neonate)"), #to be ordered as they appear in tibble
                      values = c("#E597B9", "#332288", "#44AA99"))+            #legend and color adjustments
  geom_polygon(data = hulls_all_sex_vertebrae, aes(x = x, y = y, colour = sex, fill = sex), 
               alpha = .3, show.legend = FALSE, inherit.aes = F)+ #colored hulls with transparency
  scale_fill_manual(name = "Sex", labels = c( "Female"   ,    "Male"      , "Unknown (neonate)"),
                    values = c("#E597B9", "#332288", "#44AA99"))+ #must match scale_colour_manual
  scale_shape_manual(name = "Vertebra", labels = c("Thoracic","Lumbar", "Caudal"), #copy from as.factor(sex)
                     values = shapes)+
  theme_bw()+
  xlab("PC 1 (28.73%)")+ #copy this from PC1_all object
  ylab("PC 2 (26.97%)")+ #copy this from PC2_all object
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  #Remove legend for a scale_ using guide
  guides(fill = guide_legend(label = F, title = NULL, override.aes = list(shape = NA)))

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_vertebrae_sex_ggplot

############################################################################################

#Make hulls for PCA plot with hulls around vertebra type
hulls_all_vertebra_vertebrae  <- pcscores_all_vertebrae  %>%
  group_by(vertebra) %>%
  slice(chull(comp1, comp2)) %>%
  rename(x = comp1, y = comp2)

#Nice PCA plot with hulls around vertebra
PCA_all_vertebrae_vertebra_ggplot <- ggplot(pcscores_all_vertebrae, aes(x = comp1, y = comp2, colour = vertebra))+
  geom_point(size = 3, aes(shape = age), fill = "white")+
  scale_colour_manual(name = "Vertebra", labels =  c("Thoracic","Lumbar", "Caudal"), #to be ordered as they appear in tibble
                      values = mypalette_age_sex )+            #legend and color adjustments
  geom_polygon(data = hulls_all_vertebra_vertebrae, aes(x = x, y = y, colour = vertebra, fill = vertebra), 
               alpha = .3, show.legend = FALSE, inherit.aes = F)+ #colored hulls with transparency
  scale_fill_manual(name = "Vertebra", labels = c("Thoracic","Lumbar", "Caudal"),
                    values =  mypalette_age_sex)+ #must match scale_colour_manual
  scale_shape_manual(name = "Age", labels = c( "adult"  ,  "juvenile", "newborn"), #copy from as.factor(age)
                     values = shapes)+
  theme_bw()+
  xlab("PC 1 (28.73%)")+ #copy this from PC1_all object
  ylab("PC 2 (26.97%)")+ #copy this from PC2_all object
  ggtitle("PCA all data vertebrae")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  #Remove legend for a scale_ using guide
  guides(fill = guide_legend(label = F, title = NULL, override.aes = list(shape = NA)))

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_vertebrae_vertebra_ggplot

#Make hulls for PCA plot with hulls around specimen
hulls_all_specimen_vertebrae  <- pcscores_all_vertebrae  %>%
  group_by(specimens) %>%
  slice(chull(comp1, comp2)) %>%
  rename(x = comp1, y = comp2)

#Nice PCA plot with hulls around sex
PCA_all_vertebrae_specimen_ggplot <- ggplot(pcscores_all_vertebrae, aes(x = comp1, y = comp2, colour = specimens))+
  geom_point(size = 3, aes(shape = vertebra), fill = "white")+
  geom_polygon(data = hulls_all_specimen_vertebrae, aes(x = x, y = y, colour = specimens, fill = specimens), 
               alpha = .1, show.legend = FALSE, inherit.aes = F)+ #colored hulls with transparency
  scale_shape_manual(name = "Vertebra", labels = c("Thoracic","Lumbar", "Caudal"), #copy from as.factor(sex)
                     values = shapes)+
  theme_bw()+
  xlab("PC 1 (28.73%)")+ #copy this from PC1_all object
  ylab("PC 2 (26.97%)")+ #copy this from PC2_all object
  ggtitle("PCA all data vertebrae")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  #Remove legend for a scale_ using guide
  guides(fill = guide_legend(label = F, title = NULL, override.aes = list(shape = NA)))

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_vertebrae_specimen_ggplot

###Regression PC1 and PC2 ----

#Calculate regression for each component for size
reg_PC1all_vertebrae_size <- lm(comp1 ~ size, data = pcscores_all_vertebrae)
reg_PC2all_vertebrae_size <- lm(comp2 ~ size, data = pcscores_all_vertebrae)

#View results and p-value
summary(reg_PC1all_vertebrae_size)
summary(reg_PC2all_vertebrae_size)
anova(reg_PC1all_vertebrae_size)
anova(reg_PC2all_vertebrae_size)

#Save results of significant regression to file
sink("output/PC1-2all_vertebrae_size_lm.txt")
print("PC1")
summary(reg_PC1all_vertebrae_size)
anova(reg_PC1all_vertebrae_size)
print("PC2")
summary(reg_PC2all_vertebrae_size)
anova(reg_PC2all_vertebrae_size)
sink() 

#Calculate regression for each component taking sex into account
reg_PC1all_vertebrae_sex <- lm(comp1 ~ sex, data = pcscores_all_vertebrae)
reg_PC2all_vertebrae_sex <- lm(comp2 ~ sex, data = pcscores_all_vertebrae)

#View results and p-value
summary(reg_PC1all_vertebrae_sex)
summary(reg_PC2all_vertebrae_sex)
anova(reg_PC1all_vertebrae_sex)
anova(reg_PC2all_vertebrae_sex)

#Save results of significant regression to file
sink("output/PC1-2all_vertebrae_sex_lm.txt")
print("PC1")
summary(reg_PC1all_vertebrae_sex)
anova(reg_PC1all_vertebrae_sex)
print("PC2")
summary(reg_PC2all_vertebrae_sex)
anova(reg_PC2all_vertebrae_sex)
sink() 

#Calculate regression for each component taking age into account
reg_PC1all_vertebrae_age <- lm(comp1 ~ age, data = pcscores_all_vertebrae)
reg_PC2all_vertebrae_age <- lm(comp2 ~ age, data = pcscores_all_vertebrae)

#View results and p-value
summary(reg_PC1all_vertebrae_age)
summary(reg_PC2all_vertebrae_age)
anova(reg_PC1all_vertebrae_age)
anova(reg_PC2all_vertebrae_age)

#Save results of significant regression to file
sink("output/PC1-2all_vertebrae_age_lm.txt")
print("PC1")
summary(reg_PC1all_vertebrae_age)
anova(reg_PC1all_vertebrae_age)
print("PC2")
summary(reg_PC2all_vertebrae_age)
anova(reg_PC2all_vertebrae_age)
sink() 

#Calculate regression for each component taking specimens into account
reg_PC1all_vertebrae_specimens <- lm(comp1 ~ specimens, data = pcscores_all_vertebrae)
reg_PC2all_vertebrae_specimens <- lm(comp2 ~ specimens, data = pcscores_all_vertebrae)

#View results and p-value
summary(reg_PC1all_vertebrae_specimens)
summary(reg_PC2all_vertebrae_specimens)
anova(reg_PC1all_vertebrae_specimens)
anova(reg_PC2all_vertebrae_specimens)

#Save results of significant regression to file
sink("output/PC1-2all_vertebrae_specimens_lm.txt")
print("PC1")
summary(reg_PC1all_vertebrae_specimens)
anova(reg_PC1all_vertebrae_specimens)
print("PC2")
summary(reg_PC2all_vertebrae_specimens)
anova(reg_PC2all_vertebrae_specimens)
sink() 

#Calculate regression for each component taking vertebra into account
reg_PC1all_vertebrae_vertebra <- lm(comp1 ~ vertebra, data = pcscores_all_vertebrae)
reg_PC2all_vertebrae_vertebra <- lm(comp2 ~ vertebra, data = pcscores_all_vertebrae)

#View results and p-value
summary(reg_PC1all_vertebrae_vertebra)
summary(reg_PC2all_vertebrae_vertebra)
anova(reg_PC1all_vertebrae_vertebra)
anova(reg_PC2all_vertebrae_vertebra)

#Save results of significant regression to file
sink("output/PC1-2all_vertebrae_vertebra_lm.txt")
print("PC1")
summary(reg_PC1all_vertebrae_vertebra)
anova(reg_PC1all_vertebrae_vertebra)
print("PC2")
summary(reg_PC2all_vertebrae_vertebra)
anova(reg_PC2all_vertebrae_vertebra)
sink() 

#Save results of all regressions to 1 file
sink("output/PC1-2_all_vertebrae_lm.txt")
print("PC1")
anova(reg_PC1all_vertebrae_size)
anova(reg_PC1all_vertebrae_sex)
anova(reg_PC1all_vertebrae_age)
anova(reg_PC1all_vertebrae_specimens)
anova(reg_PC1all_vertebrae_vertebra)
print("PC2")
anova(reg_PC2all_vertebrae_size)
anova(reg_PC2all_vertebrae_sex)
anova(reg_PC2all_vertebrae_age)
anova(reg_PC2all_vertebrae_specimens)
anova(reg_PC2all_vertebrae_vertebra)
sink()

#ANOVAs OF PC SCORES AND SEX/AGE BY VERTEBRA TYPE AND CONSENSUS ----
##Vertebra type ----
#Is there differences vertebra shape between sexes and ages?
#Taking into account the shape difference based on column position
#Shape is expressed as pc scores

#Make unique row names
rownames(gdf_all$pcscores) <- make.names(gdf_all$specimen, unique = TRUE)

#Calculate null model - changes in shape based on vertebra type
shape_vert_null <- lm.rrpp(pcscores ~ vertebra, data = gdf_all, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis null model
summary(shape_vert_null)

#Save results onto file
sink("Output/shape_vertebra_model.txt")
print("Shape (pc scores) not size corrected by vertebra type")
summary(shape_vert_null) 
sink() 

###Sex ----
##Test different shape between types of vertebrae ----
shape_vert_sex_comb <- lm.rrpp(pcscores ~ vertebra + sex, data = gdf_all, iter=999, print.progress = TRUE) 
shape_vert_sex_int <-  lm.rrpp(pcscores ~ vertebra * sex, data = gdf_all, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of shape with logCS
summary(shape_vert_sex_comb)
summary(shape_vert_sex_int) 

#ANOVAs - is a model significantly better than the others?
anova_shape_sex_models <- anova(shape_vert_null, shape_vert_sex_comb, shape_vert_sex_int)
anova_shape_sex_models

#Save results of significant regression to file
sink("shape_sex_models.txt")
print("Combination +")
summary(shape_vert_sex_comb) 

print("Interaction *")
summary(shape_vert_sex_int)

print("ANOVA")
anova_shape_sex_models
sink() 

###Age ----
##Test different shape between types of vertebrae ----
shape_vert_age_comb <- lm.rrpp(pcscores ~ vertebra + age, data = gdf_all, iter=999, print.progress = TRUE) 
shape_vert_age_int <-  lm.rrpp(pcscores ~ vertebra * age, data = gdf_all, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of shape with logCS
summary(shape_vert_age_comb)
summary(shape_vert_age_int) 

#Save results of significant regression to file
sink("shape_age_models.txt")
print("Combination +")
summary(shape_vert_age_comb) 

print("Interaction *")
summary(shape_vert_age_int)
sink() 

#ANOVAs - is a model significantly better than the others?
anova_shape_age_models <- anova(shape_vert_null, shape_vert_age_comb, shape_vert_age_int)
anova_shape_age_models

#Pairwise comparison for the combination and interaction model
#Helps determine if there is a significant difference in slope (int model) in the shape_age trajectory on top of difference in intercept (comb model)
pairwise_shape_age_vert <- pairwise(shape_vert_age_int, fit.null = shape_vert_null,
                                    groups = gdf_all$age, print.progress = FALSE) 
pairwise_shape_age_vert

#Distances between slope vectors (end-points) - absolute difference between slopes of groups
#if significant means int model better than comb
pairwise_shape_age_vert_dist <- summary(pairwise_shape_age_vert, confidence = 0.95, test.type = "dist") 
pairwise_shape_age_vert_dist

#Correlation between slope vectors (and angles) - similarity of vector orientation or angle,
#if significant means the vectors of the groups are oriented in different ways 
pairwise_shape_age_vert_VC <- summary(pairwise_shape_age_vert, confidence = 0.95, test.type = "VC",
                                      angle.type = "deg") 
pairwise_shape_age_vert_VC

#Absolute difference between slope vector lengths - difference in rate of change per covariate unit (size),
#if significant means there is a significant rate of change difference in shape between groups during growth
pairwise_shape_age_vert_DL <-summary(pairwise_shape_age_vert, confidence = 0.95, test.type = "DL") 
pairwise_shape_age_vert_DL

#Save results to file
sink("pairwise_shape_age_vert.txt")
print("ANOVA models")
print(anova_shape_age_models)

print("1-Pairwise absolute distances slopes")
summary(pairwise_shape_age_vert, confidence = 0.95, test.type = "dist") 

print("2-Distance between angles (slope directions)")
summary(pairwise_shape_age_vert, confidence = 0.95, test.type = "VC",
        angle.type = "deg") 

print("3-Difference in slope vector length (difference in rate of change of shape per unit of size)")
summary(pairwise_shape_age_vert, confidence = 0.95, test.type = "DL") 

sink()

###Sex vs Size ----
##Test different size between types of vertebrae and sex ----
size_vert_sex_null <- lm.rrpp(size ~ vertebra, data = gdf_all, iter=999, print.progress = TRUE) 
size_vert_sex_comb <- lm.rrpp(size ~ vertebra + sex, data = gdf_all, iter=999, print.progress = TRUE) 
size_vert_sex_int <-  lm.rrpp(size ~ vertebra * sex, data = gdf_all, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of shape with logCS
summary(size_vert_sex_null)
summary(size_vert_sex_comb)
summary(size_vert_sex_int) 

#ANOVAs - is a model significantly better than the others?
anova_size_sex_models <- anova(size_vert_sex_null, size_vert_sex_comb, size_vert_sex_int)
anova_size_sex_models

#Save results of significant regression to file
sink("shape_sex_models.txt")
print("Null")
summary(size_vert_sex_null)

print("Combination +")
summary(size_vert_sex_comb) 

print("Interaction *")
summary(size_vert_sex_int)
sink() 

#Pairwise comparison for the combination and interaction model
#Helps determine if there is a significant difference in slope (int model) in the size trajectory on top of difference in intercept (comb model)
pairwise_size_sex_vert <- pairwise(size_vert_sex_int, fit.null = size_vert_sex_null,
                                    groups = gdf_all$sex, print.progress = FALSE) 
pairwise_size_sex_vert

#Distances between slope vectors (end-points) - absolute difference between slopes of groups
#if significant means int model better than comb
pairwise_size_sex_vert_dist <- summary(pairwise_size_sex_vert, confidence = 0.95, test.type = "dist") 
pairwise_size_sex_vert_dist

#Correlation between slope vectors (and angles) - similarity of vector orientation or angle,
#if significant means the vectors of the groups are oriented in different ways 
pairwise_size_sex_vert_VC <- summary(pairwise_size_sex_vert, confidence = 0.95, test.type = "VC",
                                      angle.type = "deg") 
pairwise_size_sex_vert_VC

#Absolute difference between slope vector lengths - difference in rate of change per covariate unit (size),
#if significant means there is a significant rate of change difference in shape between groups during growth
pairwise_size_sex_vert_DL <-summary(pairwise_size_sex_vert, confidence = 0.95, test.type = "DL") 
pairwise_size_sex_vert_DL

#Save results to file
sink("pairwise_size_sex_vert.txt")
print("ANOVA models")
print(anova_size_sex_models)

print("1-Pairwise absolute distances slopes")
summary(pairwise_size_sex_vert, confidence = 0.95, test.type = "dist") 

print("2-Distance between angles (slope directions)")
summary(pairwise_size_sex_vert, confidence = 0.95, test.type = "VC",
        angle.type = "deg") 

print("3-Difference in slope vector length (difference in rate of change of shape per unit of size)")
summary(pairwise_size_sex_vert, confidence = 0.95, test.type = "DL") 

sink()

####NOT RUN#####
##Consensus ----
#Is there differences in consensus column shape between sexes and ages?
#Taking into account the shape difference based on column position
#Shape is expressed as pc scores

#Make unique row names
rownames(gdf_consensus$pcscores) <- make.names(gdf_consensus$specimen, unique = TRUE)

#Calculate null model - changes in shape based on vertebra type
shape_all_vert_null <- lm.rrpp(pcscores ~ vertebra, data = gdf_consensus, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis null model
summary(shape_all_vert_null)

#Save results onto file
sink("Output/shape_all_specimen_model.txt")
print("Shape (pc scores) not size corrected by vertebra type")
summary(shape_all_vert_null) 
sink() 

###Sex ----
##Test different shape between types of vertebrae ----
shape_all_sex_comb <- lm.rrpp(pcscores ~ vertebra + sex, data = gdf_consensus, iter=999, print.progress = TRUE) 
shape_all_sex_int <-  lm.rrpp(pcscores ~ vertebra * sex, data = gdf_consensus, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of shape with logCS
summary(shape_all_sex_comb)
summary(shape_all_sex_int) 

#ANOVAs - is a model significantly better than the others?
anova_shape_all_sex_models <- anova(shape_all_vert_null, shape_all_sex_comb, shape_all_sex_int)
anova_shape_all_sex_models

#Save results of significant regression to file
sink("Output/shape_sex_models.txt")
print("Combination +")
summary(shape_all_sex_comb) 

print("Interaction *")
summary(shape_all_sex_int)

print("ANOVA")
anova_shape_all_sex_models
sink() 

###Age ----
##Test different shape between types of vertebrae ----
shape_all_age_comb <- lm.rrpp(pcscores ~ vertebra + age, data = gdf_consensus, iter=999, print.progress = TRUE) 
shape_all_age_int <-  lm.rrpp(pcscores ~ vertebra * age, data = gdf_consensus, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of shape with logCS
summary(shape_all_age_comb)
summary(shape_all_age_int) 

#Save results of significant regression to file
sink("Output/shape_age_models.txt")
print("Combination +")
summary(shape_all_age_comb) 

print("Interaction *")
summary(shape_all_age_int)
sink() 

#ANOVAs - is a model significantly better than the others?
anova_shape_all_age_models <- anova(shape_all_vert_null, shape_all_age_comb, shape_all_age_int)
anova_shape_all_age_models

#Pairwise comparison for the combination and interaction model
#Helps determine if there is a significant difference in slope (int model) in the shape_all_age trajectory on top of difference in intercept (comb model)
pairwise_shape_all_age_vert <- pairwise(shape_all_age_int, fit.null = shape_all_vert_null,
                                    groups = gdf_consensus$age, print.progress = FALSE) 
pairwise_shape_all_age_vert

#Distances between slope vectors (end-points) - absolute difference between slopes of groups
#if significant means int model better than comb
pairwise_shape_all_age_vert_dist <- summary(pairwise_shape_all_age_vert, confidence = 0.95, test.type = "dist") 
pairwise_shape_all_age_vert_dist

#Correlation between slope vectors (and angles) - similarity of vector orientation or angle,
#if significant means the vectors of the groups are oriented in different ways 
pairwise_shape_all_age_vert_VC <- summary(pairwise_shape_all_age_vert, confidence = 0.95, test.type = "VC",
                                      angle.type = "deg") 
pairwise_shape_all_age_vert_VC

#Absolute difference between slope vector lengths - difference in rate of change per covariate unit (size),
#if significant means there is a significant rate of change difference in shape between groups during growth
pairwise_shape_all_age_vert_DL <-summary(pairwise_shape_all_age_vert, confidence = 0.95, test.type = "DL") 
pairwise_shape_all_age_vert_DL

#Save results to file
sink("Output/pairwise_shape_all_age_vert.txt")
print("ANOVA models")
print(anova_shape_all_age_models)

print("1-Pairwise absolute distances slopes")
summary(pairwise_shape_all_age_vert, confidence = 0.95, test.type = "dist") 

print("2-Distance between angles (slope directions)")
summary(pairwise_shape_all_age_vert, confidence = 0.95, test.type = "VC",
        angle.type = "deg") 

print("3-Difference in slope vector length (difference in rate of change of shape per unit of size)")
summary(pairwise_shape_all_age_vert, confidence = 0.95, test.type = "DL") 

sink()

#PCA CONSENSUS ----

#Create tibble with new pc scores from gdf_consensus$pcscores and other concatenated variables
pcscores_all_consensus <- as_tibble(gdf_consensus$pcscores)

#Add labels and other attributes to tibble as columns
pcscores_all_consensus <- pcscores_all_consensus %>% mutate(age = gdf_consensus$age, sex = gdf_consensus$sex, 
                                                            specimens = gdf_consensus$specimen, size = gdf_consensus$size,
                                                            vertebra = gdf_consensus$vertebra)
glimpse(pcscores_all_consensus)

#PC1 and PC2 % values can be estimated as the average of each type of specimen

#PC1
#Calculate mean of all 3 values
PC1_consensus <- pca_all[["result"]][["AVE"]][["AVE_X"]][[4]][["comp1"]]*100

#PC2
#Calculate mean of all 3 values
PC2_consensus <- pca_all[["result"]][["AVE"]][["AVE_X"]][[4]][["comp2"]]*100

#Plot
PCA_all_consensus_ggplot <- ggplot(pcscores_all_consensus, aes(x = comp1, y = comp2, colour = age, fill = age))+
  geom_point(size = 3, aes(shape = vertebra))+
  scale_colour_manual(name = "Age", labels =  c("adult"  ,  "juvenile", "newborn" ), #to be ordered as they appear in tibble
                      values = mypalette_age_sex , aesthetics = c("colour","fill"))+            #legend and color adjustments
  scale_shape_manual(name = "Vertebra", labels = c("Thoracic","Lumbar", "Caudal"), #copy from as.factor(sex)
                     values = shapes)+
  theme_bw()+
  xlab("PC 1 (23.83%)")+ #copy this from PC1_consensus object
  ylab("PC 2 (13.43%)")+ #copy this from PC2_consensus object
  ggtitle("PCA all data consensus")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_consensus_ggplot

#Make hulls for PCA plot with hulls around age
hulls_all_age_consensus <- pcscores_all_consensus %>%
  group_by(age) %>%
  slice(chull(comp1, comp2)) %>%
  rename(x = comp1, y = comp2)

#Nice PCA plot with hulls around age
PCA_all_consensus_age_ggplot <- ggplot(pcscores_all_consensus, aes(x = comp1, y = comp2, colour = age))+
  geom_point(size = 3, aes(shape = vertebra), fill = "white")+
  scale_colour_manual(name = "Age", labels =  c("adult" ,   "juvenile" ,"neonate" ), #to be ordered as they appear in tibble
                      values = mypalette_age_sex)+            #legend and color adjustments
  geom_polygon(data = hulls_all_age_consensus, aes(x = x, y = y, colour = age, fill = age), 
               alpha = .3, show.legend = FALSE, inherit.aes = F)+ #colored hulls with transparency
  scale_fill_manual(name = "Age", labels = c("adult"  ,  "juvenile", "neonate" ),
                    values =  mypalette_age_sex )+ #must match scale_colour_manual
  scale_shape_manual(name = "Vertebra", labels = c("Thoracic","Lumbar", "Caudal"), #copy from as.factor(sex)
                     values = shapes)+
  theme_bw()+
  xlab("PC 1 (23.83%)")+ #copy this from PC1_consensus object
  ylab("PC 2 (13.43%)")+ #copy this from PC2_consensus object
 
  #Remove legend for a scale_ using guide
  guides(fill = guide_legend(label = F, title = NULL, override.aes = list(shape = NA)))

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_consensus_age_ggplot

#Make hulls for PCA plot with hulls around sex
hulls_all_sex_consensus <- pcscores_all_consensus %>%
  group_by(sex) %>%
  slice(chull(comp1, comp2)) %>%
  rename(x = comp1, y = comp2)

#Nice PCA plot with hulls around sex
PCA_all_consensus_sex_ggplot <- ggplot(pcscores_all_consensus, aes(x = comp1, y = comp2, colour = sex))+
  geom_point(size = 3, aes(shape = vertebra), fill = "white")+
  scale_colour_manual(name = "Sex", labels =  c("F" ,   "M" ,"unknown" ), #to be ordered as they appear in tibble
                      values = mypalette_age_sex)+            #legend and color adjustments
  geom_polygon(data = hulls_all_sex_consensus, aes(x = x, y = y, colour = sex, fill = sex), 
               alpha = .3, show.legend = FALSE, inherit.aes = F)+ #colored hulls with transparency
  scale_fill_manual(name = "Sex", labels =  c("F" ,   "M" ,"unknown" ),
                    values =  mypalette_age_sex )+ #must match scale_colour_manual
  scale_shape_manual(name = "Vertebra", labels = c("Thoracic","Lumbar", "Caudal"), #copy from as.factor(sex)
                     values = shapes)+
  theme_bw()+
  xlab("PC 1 (23.83%)")+ #copy this from PC1_consensus object
  ylab("PC 2 (13.43%)")+ #copy this from PC2_consensus object
  ggtitle("PCA all data consensus")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  #Remove legend for a scale_ using guide
  guides(fill = guide_legend(label = F, title = NULL, override.aes = list(shape = NA)))

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_consensus_sex_ggplot

#Make hulls for PCA plot with hulls around vertebra type
hulls_all_vertebra_consensus  <- pcscores_all_consensus  %>%
  group_by(vertebra) %>%
  slice(chull(comp1, comp2)) %>%
  rename(x = comp1, y = comp2)

#Nice PCA plot with hulls around sex
PCA_all_consensus_vertebra_ggplot <- ggplot(pcscores_all_consensus, aes(x = comp1, y = comp2, colour = vertebra))+
  geom_point(size = 3, aes(shape = age), fill = "white")+
  scale_colour_manual(name = "Vertebra", labels =  c("Thoracic","Lumbar", "Caudal"), #to be ordered as they appear in tibble
                      values = mypalette_age_sex )+            #legend and color adjustments
  geom_polygon(data = hulls_all_vertebra_consensus, aes(x = x, y = y, colour = vertebra, fill = vertebra), 
               alpha = .3, show.legend = FALSE, inherit.aes = F)+ #colored hulls with transparency
  scale_fill_manual(name = "Vertebra", labels = c("Thoracic","Lumbar", "Caudal"),
                    values =  mypalette_age_sex)+ #must match scale_colour_manual
  scale_shape_manual(name = "Age", labels = c( "adult"  ,  "juvenile", "newborn"), #copy from as.factor(age)
                     values = shapes)+
  theme_bw()+
  xlab("PC 1 (23.83%)")+ #copy this from PC1_consensus object
  ylab("PC 2 (13.43%)")+ #copy this from PC2_consensus object
  ggtitle("PCA all data consensus")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  #Remove legend for a scale_ using guide
  guides(fill = guide_legend(label = F, title = NULL, override.aes = list(shape = NA)))

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_consensus_vertebra_ggplot

#GRAPHICS CODE ----
##Create project palettes----
mypalette_paired <- brewer.pal(12,"Paired")
image(1:12, 1, as.matrix(1:12), col = mypalette_paired, xlab = "Paired",
      ylab = "", yaxt = "n")

mypalette_blue <- as.matrix(ggthemes_data[["tableau"]][["color-palettes"]][["ordered-sequential"]][["Blue"]][["value"]])
image(1:20, 1, as.matrix(1:20), col = mypalette_blue, xlab = "Blue",
      ylab = "", yaxt = "n")

install.packages("rcartocolor")
library(rcartocolor)
display_carto_all(colorblind_friendly = TRUE)
display_carto_pal(12, "Safe")


#Palette for vertebrae
mypalette_avertebrae <- c(mypalette_paired[3], mypalette_paired[8], mypalette_paired[9])
image(1:3, 1, as.matrix(1:3), col = mypalette_age_sex,
      ylab = "", yaxt = "n")

#Palette for age - early, late/new, immature, adult
mypalette_age<- c(mypalette_paired[12], mypalette_paired[11], mypalette_paired[9])
image(1:3, 1, as.matrix(1:3), col = mypalette_age_sex,
      ylab = "", yaxt = "n")

#Palette to use for sex PCAs
mypalette_sex <- c(mypalette_paired[5], mypalette_paired[2], mypalette_paired[10])
image(1:3, 1, as.matrix(1:3), col = mypalette_age_sex,
      ylab = "", yaxt = "n")

#create shape palette for sex
shapes <- c(19, 15, 17)

#Create shape palette 3 age
shapes <- c(21, 22, 24) 


#Palette for vertebrae - thoracic, lumbar, caudal, consensus
mypalette_specimen <- c(mypalette_blue[3,], mypalette_blue[7,], mypalette_blue[12,], mypalette_blue[20,])
image(1:4, 1, as.matrix(1:4), col = mypalette_specimen,
      ylab = "", yaxt = "n")

shapes2 <- c(21, 22, 23, 24)