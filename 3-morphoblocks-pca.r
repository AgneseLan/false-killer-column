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

#GPA data with extra vertebrae removed
coords_thoracic_blocks <- gdf_thoracic$coords
coords_lumbar_blocks <- gdf_lumbar$coords[,,-lumbar_delete]
coords_caudal_blocks <- gdf_caudal$coords[,,-caudal_delete]

#CS with extra vertebrae removed
block_thoracic_logCS <- gdf_thoracic$size
block_lumbar_logCS <- gdf_lumbar$size[-lumbar_delete]
block_caudal_logCS <- gdf_caudal$size[-caudal_delete]

#MORPHOBLOCKS#
#Create coordinate blocks
block_thoracic <- formatBlock(coords_thoracic_blocks, cs = block_thoracic_logCS, k = 3, gpa = F)
block_lumbar <- formatBlock(coords_lumbar_blocks, cs = block_lumbar_logCS, k = 3, gpa = F)
block_caudal <- formatBlock(coords_caudal_blocks, cs = block_caudal_logCS, k = 3, gpa = F)

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

#Visualize loadings for the first two components
loadingsPlot(pca_all, comp = 1, cex.3d = 10)
rgl.snapshot(filename = "Output/loadings_comp1.png")
clear3d()

loadingsPlot(pca_all, comp = 2, cex.3d = 7)
rgl.snapshot(filename = "Output/loadings_comp2.png")


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
  geom_text_repel(colour = "black", size = 4, max.overlaps = 200)+
  scale_colour_manual(name = "Age", labels =  c("adult"  ,  "juvenile", "Neonate" ), #to be ordered as they appear in tibble
                      values = mypalette_age , aesthetics = c("colour","fill"))+            #legend and color adjustments
  scale_shape_manual(name = "Sex", labels = c( "F"   ,    "M"      , "unknown"), #copy from as.factor(sex)
                     values = shapes)+
  theme_bw()+
  xlab(paste0("PC 1 (",round(PC1_all, digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(PC2_all, digits = 2),"%)"))+ 
  ggtitle("PCA all data vertebrae")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_vertebrae_ggplot

#Make hulls for PCA plot with hulls around age
hulls_all_age_vertebrae <- pcscores_all_vertebrae %>%
  group_by(age) %>%
  slice(chull(comp1, comp2)) %>%
  rename(x = comp1, y = comp2)


PCA_all_vertebrae_age_ggplot <- ggplot(pcscores_all_vertebrae, aes(x = comp1, y = comp2, colour = age))+
  geom_point(size = 3, aes(shape = vertebra, fill = age), color = "black") +
  geom_polygon(data = hulls_all_age_vertebrae, aes(x = x, y = y, colour = age, fill = age), 
               alpha = .3, show.legend = FALSE, inherit.aes = FALSE) +
  scale_colour_manual(name = "Age", labels = c("Adult", "Juvenile", "Neonate"),
                      values = mypalette_age) +
  scale_fill_manual(name = "Age", labels = c("Adult", "Juvenile", "Neonate"),
                    values = mypalette_age) +
  scale_shape_manual(name = "Vertebra", labels = c("Thoracic","Lumbar", "Caudal"), #copy from as.factor(sex)
                     values = shapes)+
  theme_bw(base_size = 14)+
  xlab(paste0("PC 1 (",round(PC1_all, digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(PC2_all, digits = 2),"%)"))+ 
  guides(color = guide_legend(override.aes = list(shape = 23, colour = "black", fill = mypalette_age))) +
  theme(legend.key = element_blank(), legend.title = element_text(size = 15, face = "bold"), legend.text = element_text(size = 14),
        legend.background = element_blank(), legend.box.background =  element_blank())

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
  geom_point(size = 3, aes(shape = vertebra, fill = sex), color = "black") +
  geom_polygon(data = hulls_all_sex_vertebrae[!hulls_all_sex_vertebrae$sex%in% c("unknown"),], aes(x = x, y = y, colour = sex, fill = sex), 
               alpha = .3, show.legend = FALSE, inherit.aes = F)+ #colored hulls with transparency
  scale_colour_manual(name = "Sex", labels =  c("Female"   ,    "Male"), #to be ordered as they appear in tibble
                      values= mypalette_sex)+            #legend and color adjustments
  scale_fill_manual(name = "Sex", labels = c( "Female"   ,    "Male"),
                    values =   mypalette_sex)+ #must match scale_colour_manual
  scale_shape_manual(name = "Vertebra", labels = c("Thoracic","Lumbar", "Caudal"), #copy from as.factor(sex)
                     values = shapes)+
  theme_bw(base_size = 14)+
  xlab(paste0("PC 1 (",round(PC1_all, digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(PC2_all, digits = 2),"%)"))+ 
  guides(color = guide_legend(override.aes = list(shape = 23, colour = "black", fill = mypalette_sex))) +
  theme(legend.key = element_blank(), legend.title = element_text(size = 15, face = "bold"), legend.text = element_text(size = 14),
        legend.background = element_blank(), legend.box.background =  element_blank())

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_vertebrae_sex_ggplot

#Palette for vertebrae
mypalette_paired <- brewer.pal(12,"Paired")
image(1:12, 1, as.matrix(1:12), col = mypalette_paired, xlab = "Paired",
      ylab = "", yaxt = "n")

mypalette_vertebrae <- c(mypalette_paired[3], mypalette_paired[8], mypalette_paired[9])
image(1:3, 1, as.matrix(1:3), col = mypalette_vertebrae,
      ylab = "", yaxt = "n")

#Make hulls for PCA plot with hulls around vertebra type
hulls_all_vertebra_vertebrae  <- pcscores_all_vertebrae  %>%
  group_by(vertebra) %>%
  slice(chull(comp1, comp2)) %>%
  rename(x = comp1, y = comp2)

#Nice PCA plot with hulls around vertebra
PCA_all_vertebrae_vertebra_ggplot <- ggplot(pcscores_all_vertebrae, aes(x = comp1, y = comp2, colour = vertebra))+
  geom_point(size = 3, aes(shape = age, fill = vertebra), color = "black")+
  geom_polygon(data = hulls_all_vertebra_vertebrae, aes(x = x, y = y, colour = vertebra, fill = vertebra), 
               alpha = .3, show.legend = FALSE, inherit.aes = F)+ #colored hulls with transparency
  scale_fill_manual(name = "Vertebra", labels = c("Thoracic","Lumbar", "Caudal"),
                    values =  mypalette_vertebrae)+ #must match scale_colour_manual
  scale_colour_manual(name = "Vertebra", labels =  c("Thoracic","Lumbar", "Caudal"), #to be ordered as they appear in tibble
                      values = mypalette_vertebrae)+            #legend and color adjustments
  scale_shape_manual(name = "Age", labels = c( "Adult"  ,  "Juvenile", "Neonate"), #copy from as.factor(age)
                     values = shapes)+
  theme_bw()+
  xlab(paste0("PC 1 (",round(PC1_all, digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(PC2_all, digits = 2),"%)"))+ 
  guides(color = guide_legend(override.aes = list(shape = 23, colour = "black", fill = mypalette_vertebrae))) +
  theme(legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
        legend.background = element_blank(), legend.box.background =  element_blank())

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_vertebrae_vertebra_ggplot

####Predict shapes based on each vertebra scores ----
#PC1
pcs_thoracic_1 <- pca_all[["scores"]][[1]][,1]
pcs_lumbar_1 <- pca_all[["scores"]][[2]][,1]
pcs_caudal_1 <- pca_all[["scores"]][[3]][,1]

preds_thoracic_1 <- shape.predictor(block_thoracic@gpa.3D, x= pcs_thoracic_1, Intercept = FALSE, 
                                   pred_min = min(pcs_thoracic_1), pred_max = max(pcs_thoracic_1))
preds_lumbar_1 <- shape.predictor(block_lumbar@gpa.3D, x= pcs_lumbar_1, Intercept = FALSE, 
                                    pred_min = min(pcs_lumbar_1), pred_max = max(pcs_lumbar_1))
preds_caudal_1 <- shape.predictor(block_caudal@gpa.3D, x= pcs_caudal_1, Intercept = FALSE, 
                                    pred_min = min(pcs_caudal_1), pred_max = max(pcs_caudal_1))

#3D deformation using PLY files
#Extremes PC axes

##Create warp mesh, to use as reference for visualization of analyses
PC1min_morphoblocks_thoracic_mesh <- warpRefMesh(mesh = refmesh_thoracic, mesh.coord = warp_specimen_thoracic, 
                                      ref = preds_thoracic_1$pred_min, color = "gray", centered = FALSE)
##Save warped mesh as ply (Morpho)
mesh2ply(PC1min_morphoblocks_thoracic_mesh, filename = "Output/PC1min_morphoblocks_thoracic")

##Create warp mesh, to use as reference for visualization of analyses
PC1max_morphoblocks_thoracic_mesh <- warpRefMesh(mesh = refmesh_thoracic, mesh.coord = warp_specimen_thoracic, 
                                      ref = preds_thoracic_1$pred_max, color = "gray", centered = FALSE)
##Save warped mesh as ply (Morpho)
mesh2ply(PC1max_morphoblocks_thoracic_mesh, filename = "Output/PC1max_morphoblocks_thoracic")

##Create warp mesh, to use as reference for visualization of analyses
PC1min_morphoblocks_lumbar_mesh <- warpRefMesh(mesh = refmesh_lumbar, mesh.coord = warp_specimen_lumbar, 
                                                 ref = preds_lumbar_1$pred_min, color = "gray", centered = FALSE)
##Save warped mesh as ply (Morpho)
mesh2ply(PC1min_morphoblocks_lumbar_mesh, filename = "Output/PC1min_morphoblocks_lumbar")

##Create warp mesh, to use as reference for visualization of analyses
PC1max_morphoblocks_lumbar_mesh <- warpRefMesh(mesh = refmesh_lumbar, mesh.coord = warp_specimen_lumbar, 
                                                 ref = preds_lumbar_1$pred_max, color = "gray", centered = FALSE)
##Save warped mesh as ply (Morpho)
mesh2ply(PC1max_morphoblocks_lumbar_mesh, filename = "Output/PC1max_morphoblocks_lumbar")

##Create warp mesh, to use as reference for visualization of analyses
PC1min_morphoblocks_caudal_mesh <- warpRefMesh(mesh = refmesh_caudal, mesh.coord = warp_specimen_caudal, 
                                                 ref = preds_caudal_1$pred_min, color = "gray", centered = FALSE)
##Save warped mesh as ply (Morpho)
mesh2ply(PC1min_morphoblocks_caudal_mesh, filename = "Output/PC1min_morphoblocks_caudal")

##Create warp mesh, to use as reference for visualization of analyses
PC1max_morphoblocks_caudal_mesh <- warpRefMesh(mesh = refmesh_caudal, mesh.coord = warp_specimen_caudal, 
                                                 ref = preds_caudal_1$pred_max, color = "gray", centered = FALSE)
##Save warped mesh as ply (Morpho)
mesh2ply(PC1max_morphoblocks_caudal_mesh, filename = "Output/PC1max_morphoblocks_caudal")

#PC2
pcs_thoracic_2 <- pca_all[["scores"]][[1]][,2]
pcs_lumbar_2 <- pca_all[["scores"]][[2]][,2]
pcs_caudal_2 <- pca_all[["scores"]][[3]][,2]

preds_thoracic_2 <- shape.predictor(block_thoracic@gpa.3D, x= pcs_thoracic_2, Intercept = FALSE, 
                                    pred_min = min(pcs_thoracic_2), pred_max = max(pcs_thoracic_2))
preds_lumbar_2 <- shape.predictor(block_lumbar@gpa.3D, x= pcs_lumbar_2, Intercept = FALSE, 
                                  pred_min = min(pcs_lumbar_2), pred_max = max(pcs_lumbar_2))
preds_caudal_2 <- shape.predictor(block_caudal@gpa.3D, x= pcs_caudal_2, Intercept = FALSE, 
                                  pred_min = min(pcs_caudal_2), pred_max = max(pcs_caudal_2))

#3D deformation using PLY files
#Extremes PC axes

##Create warp mesh, to use as reference for visualization of analyses
PC2min_morphoblocks_thoracic_mesh <- warpRefMesh(mesh = refmesh_thoracic, mesh.coord = warp_specimen_thoracic, 
                                                 ref = preds_thoracic_2$pred_min, color = "gray", centered = FALSE)
##Save warped mesh as ply (Morpho)
mesh2ply(PC2min_morphoblocks_thoracic_mesh, filename = "Output/PC2min_morphoblocks_thoracic")

##Create warp mesh, to use as reference for visualization of analyses
PC2max_morphoblocks_thoracic_mesh <- warpRefMesh(mesh = refmesh_thoracic, mesh.coord = warp_specimen_thoracic, 
                                                 ref = preds_thoracic_2$pred_max, color = "gray", centered = FALSE)
##Save warped mesh as ply (Morpho)
mesh2ply(PC2max_morphoblocks_thoracic_mesh, filename = "Output/PC2max_morphoblocks_thoracic")

##Create warp mesh, to use as reference for visualization of analyses
PC2min_morphoblocks_lumbar_mesh <- warpRefMesh(mesh = refmesh_lumbar, mesh.coord = warp_specimen_lumbar, 
                                               ref = preds_lumbar_2$pred_min, color = "gray", centered = FALSE)
##Save warped mesh as ply (Morpho)
mesh2ply(PC2min_morphoblocks_lumbar_mesh, filename = "Output/PC2min_morphoblocks_lumbar")

##Create warp mesh, to use as reference for visualization of analyses
PC2max_morphoblocks_lumbar_mesh <- warpRefMesh(mesh = refmesh_lumbar, mesh.coord = warp_specimen_lumbar, 
                                               ref = preds_lumbar_2$pred_max, color = "gray", centered = FALSE)
##Save warped mesh as ply (Morpho)
mesh2ply(PC2max_morphoblocks_lumbar_mesh, filename = "Output/PC2max_morphoblocks_lumbar")

##Create warp mesh, to use as reference for visualization of analyses
PC2min_morphoblocks_caudal_mesh <- warpRefMesh(mesh = refmesh_caudal, mesh.coord = warp_specimen_caudal, 
                                               ref = preds_caudal_2$pred_min, color = "gray", centered = FALSE)
##Save warped mesh as ply (Morpho)
mesh2ply(PC2min_morphoblocks_caudal_mesh, filename = "Output/PC2min_morphoblocks_caudal")

##Create warp mesh, to use as reference for visualization of analyses
PC2max_morphoblocks_caudal_mesh <- warpRefMesh(mesh = refmesh_caudal, mesh.coord = warp_specimen_caudal, 
                                               ref = preds_caudal_2$pred_max, color = "gray", centered = FALSE)
##Save warped mesh as ply (Morpho)
mesh2ply(PC2max_morphoblocks_caudal_mesh, filename = "Output/PC2max_morphoblocks_caudal")


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
sink("Output/PC1-2all_vertebrae_size_lm.txt")
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
sink("Output/PC1-2all_vertebrae_sex_lm.txt")
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
sink("Output/PC1-2all_vertebrae_age_lm.txt")
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
sink("Output/PC1-2all_vertebrae_specimens_lm.txt")
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
sink("Output/PC1-2all_vertebrae_vertebra_lm.txt")
print("PC1")
summary(reg_PC1all_vertebrae_vertebra)
anova(reg_PC1all_vertebrae_vertebra)
print("PC2")
summary(reg_PC2all_vertebrae_vertebra)
anova(reg_PC2all_vertebrae_vertebra)
sink() 

#Save results of all regressions to 1 file
sink("Output/PC1-2_all_vertebrae_lm.txt")
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
#Is there differences vertebra shape between vertebra type, sexes and ages across the data?
#Taking into account the shape difference based on column position
#Shape is expressed as pc scores

#Make unique row names
rownames(gdf_all$pcscores) <- make.names(gdf_all$specimen, unique = TRUE)

##Shape ----
###Vertebra type ----
#Calculate null model - changes in shape based on vertebra type
shape_vert_null <- lm.rrpp(pcscores ~ vertebra, data = gdf_all, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis null model
summary(shape_vert_null)

#Save results onto file
sink("Output/shape_vertebra_model.txt")
print("Shape (pc scores) not size corrected by vertebra type")
summary(shape_vert_null) 
sink() 

###Vertebra type vs Sex ----

#Use only males and females to avoid losing signal due to unknowns

rm_U_all <- which(gdf_all$sex == "unknown")

shape_vert_sex <-  lm.rrpp(gdf_all$pcscores[-c(rm_U_all),] ~  as.vector(gdf_all$vertebra[-rm_U_all]) * as.vector(gdf_all$sex[-rm_U_all]), iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of shape with logCS
summary(shape_vert_sex) 

#Save results of significant regression to file
sink("Output/shape_vertebra_sex_model.txt")
print("Shape ~ vertebra * sex")
summary(shape_vert_sex)
sink() 

#Calculate null model without unknowns - changes in shape based on vertebra type
shape_vert_null1 <- lm.rrpp(gdf_all$pcscores[-c(rm_U_all),] ~  as.vector(gdf_all$vertebra[-rm_U_all]), iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis null model
summary(shape_vert_null1)

#Pairwise comparison for the combination and interaction model
#Helps determine if there is a significant difference in slope (int model) in the shape_sex trajectory on top of difference in intercept (comb model)
pairwise_shape_sex_vert <- pairwise(shape_vert_sex, fit.null = shape_vert_null1,
                                    groups = gdf_all$sex[-rm_U_all], print.progress = FALSE) 
pairwise_shape_sex_vert

#Distances between slope vectors (end-points) - absolute difference between slopes of groups
#if significant means int model better than comb
pairwise_shape_sex_vert_dist <- summary(pairwise_shape_sex_vert, confidence = 0.95, test.type = "dist") 
pairwise_shape_sex_vert_dist

#Correlation between slope vectors (and angles) - similarity of vector orientation or angle,
#if significant means the vectors of the groups are oriented in different ways 
pairwise_shape_sex_vert_VC <- summary(pairwise_shape_sex_vert, confidence = 0.95, test.type = "VC",
                                      angle.type = "deg") 
pairwise_shape_sex_vert_VC

#Absolute difference between slope vector lengths - difference in rate of change per covariate unit (size),
#if significant means there is a significant rate of change difference in shape between groups during growth
pairwise_shape_sex_vert_DL <-summary(pairwise_shape_sex_vert, confidence = 0.95, test.type = "DL") 
pairwise_shape_sex_vert_DL

#Save results to file
sink("Output/pairwise_shape_sex_vert.txt")
print("1-Pairwise absolute distances slopes")
summary(pairwise_shape_sex_vert, confidence = 0.95, test.type = "dist") 

print("2-Distance between angles (slope directions)")
summary(pairwise_shape_sex_vert, confidence = 0.95, test.type = "VC",
        angle.type = "deg") 

print("3-Difference in slope vector length (difference in rate of change of shape per unit of size)")
summary(pairwise_shape_sex_vert, confidence = 0.95, test.type = "DL") 

sink()

###Vertebra type vs Age ----
shape_vert_age <-  lm.rrpp(pcscores ~ vertebra * age, data = gdf_all, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of shape with logCS
summary(shape_vert_age) 

#Save results of significant regression to file
sink("Output/shape_vertebra_age_model.txt")
print("Shape ~ vertebra * age")
summary(shape_vert_age)
sink() 

#Pairwise comparison for the combination and interaction model
#Helps determine if there is a significant difference in slope (int model) in the shape_age trajectory on top of difference in intercept (comb model)
pairwise_shape_age_vert <- pairwise(shape_vert_age, fit.null = shape_vert_null,
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
sink("Output/pairwise_shape_age_vert.txt")
print("1-Pairwise absolute distances slopes")
summary(pairwise_shape_age_vert, confidence = 0.95, test.type = "dist") 

print("2-Distance between angles (slope directions)")
summary(pairwise_shape_age_vert, confidence = 0.95, test.type = "VC",
        angle.type = "deg") 

print("3-Difference in slope vector length (difference in rate of change of shape per unit of size)")
summary(pairwise_shape_age_vert, confidence = 0.95, test.type = "DL") 

sink()

##Size ----
###Vertebra type ----
#Calculate null model - changes in size based on vertebra type
size_vert_null <- lm.rrpp(size ~ vertebra, data = gdf_all, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis null model
summary(size_vert_null)

#Save results onto file
sink("Output/size_vertebra_model.txt")
print("Size differences by vertebra type")
summary(size_vert_null) 
sink() 

###Vertebra type vs Sex ----

#Use only males and females to avoid losing signal due to unknowns

size_vert_sex <-  lm.rrpp(gdf_all$size[-rm_U_all] ~  as.vector(gdf_all$vertebra[-rm_U_all]) * as.vector(gdf_all$sex[-rm_U_all]), iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of size with logCS
summary(size_vert_sex) 

#Save results of significant regression to file
sink("Output/size_vertebra_sex_model.txt")
print("size ~ vertebra * sex")
summary(size_vert_sex)
sink() 

#Calculate null model without unknowns - changes in size based on vertebra type
size_vert_null1 <- lm.rrpp(gdf_all$size[-rm_U_all] ~  as.vector(gdf_all$vertebra[-rm_U_all]), iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis null model
summary(size_vert_null1)

#Pairwise comparison for the combination and interaction model
#Helps determine if there is a significant difference in slope (int model) in the size_sex trajectory on top of difference in intercept (comb model)
pairwise_size_sex_vert <- pairwise(size_vert_sex, fit.null = size_vert_null1,
                                    groups = gdf_all$sex[-rm_U_all], print.progress = FALSE) 
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
#if significant means there is a significant rate of change difference in size between groups during growth
pairwise_size_sex_vert_DL <-summary(pairwise_size_sex_vert, confidence = 0.95, test.type = "DL") 
pairwise_size_sex_vert_DL

#Save results to file
sink("Output/pairwise_size_sex_vert.txt")
print("1-Pairwise absolute distances slopes")
summary(pairwise_size_sex_vert, confidence = 0.95, test.type = "dist") 

print("2-Distance between angles (slope directions)")
summary(pairwise_size_sex_vert, confidence = 0.95, test.type = "VC",
        angle.type = "deg") 

print("3-Difference in slope vector length (difference in rate of change of size per unit of size)")
summary(pairwise_size_sex_vert, confidence = 0.95, test.type = "DL") 

sink()

###Vertebra type vs Age ----
size_vert_age <-  lm.rrpp(size ~ vertebra * age, data = gdf_all, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of size with logCS
summary(size_vert_age) 

#Save results of significant regression to file
sink("Output/size_vertebra_age_model.txt")
print("size ~ vertebra * age")
summary(size_vert_age)
sink() 

#Pairwise comparison for the combination and interaction model
#Helps determine if there is a significant difference in slope (int model) in the size_age trajectory on top of difference in intercept (comb model)
pairwise_size_age_vert <- pairwise(size_vert_age, fit.null = size_vert_null,
                                    groups = gdf_all$age, print.progress = FALSE) 
pairwise_size_age_vert

#Distances between slope vectors (end-points) - absolute difference between slopes of groups
#if significant means int model better than comb
pairwise_size_age_vert_dist <- summary(pairwise_size_age_vert, confidence = 0.95, test.type = "dist") 
pairwise_size_age_vert_dist

#Correlation between slope vectors (and angles) - similarity of vector orientation or angle,
#if significant means the vectors of the groups are oriented in different ways 
pairwise_size_age_vert_VC <- summary(pairwise_size_age_vert, confidence = 0.95, test.type = "VC",
                                      angle.type = "deg") 
pairwise_size_age_vert_VC

#Absolute difference between slope vector lengths - difference in rate of change per covariate unit (size),
#if significant means there is a significant rate of change difference in size between groups during growth
pairwise_size_age_vert_DL <-summary(pairwise_size_age_vert, confidence = 0.95, test.type = "DL") 
pairwise_size_age_vert_DL

#Save results to file
sink("Output/pairwise_size_age_vert.txt")
print("1-Pairwise absolute distances slopes")
summary(pairwise_size_age_vert, confidence = 0.95, test.type = "dist") 

print("2-Distance between angles (slope directions)")
summary(pairwise_size_age_vert, confidence = 0.95, test.type = "VC",
        angle.type = "deg") 

print("3-Difference in slope vector length (difference in rate of change of size per unit of size)")
summary(pairwise_size_age_vert, confidence = 0.95, test.type = "DL") 

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
PCA_all_consensus_ggplot <- ggplot(pcscores_all_consensus, aes(x = comp1, y = comp2, label = specimens, colour = age, fill = age))+
  geom_point(size = 3, aes(shape = sex))+
  geom_text_repel(colour = "black", size = 4, max.overlaps = 200)+
  scale_colour_manual(name = "Age", labels =  c("adult"  ,  "juvenile", "Neonate" ), #to be ordered as they appear in tibble
                      values = mypalette_age , aesthetics = c("colour","fill"))+            #legend and color adjustments
  scale_shape_manual(name = "Sex", labels = c( "F"   ,    "M"      , "unknown"), #copy from as.factor(sex)
                     values = shapes)+
  theme_bw()+
  xlab(paste0("PC 1 (",round(PC1_consensus, digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(PC2_consensus, digits = 2),"%)"))+ 
  ggtitle("PCA all data consensus")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_consensus_ggplot

#Make hulls for PCA plot with hulls around age
hulls_all_age_consensus <- pcscores_all_consensus %>%
  group_by(age) %>%
  slice(chull(comp1, comp2)) %>%
  rename(x = comp1, y = comp2)


PCA_all_consensus_age_ggplot <- ggplot(pcscores_all_consensus, aes(x = comp1, y = comp2, colour = age))+
  geom_point(size = 3, aes(shape = vertebra, fill = age), color = "black") +
  geom_polygon(data = hulls_all_age_consensus, aes(x = x, y = y, colour = age, fill = age), 
               alpha = .3, show.legend = FALSE, inherit.aes = FALSE) +
  scale_colour_manual(name = "Age", labels = c("Adult", "Juvenile", "Neonate"),
                      values = mypalette_age) +
  scale_fill_manual(name = "Age", labels = c("Adult", "Juvenile", "Neonate"),
                    values = mypalette_age) +
  scale_shape_manual(name = "Vertebra", labels = c("Thoracic","Lumbar", "Caudal"), #copy from as.factor(sex)
                     values = shapes)+
  theme_bw()+
  xlab(paste0("PC 1 (",round(PC1_consensus, digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(PC2_consensus, digits = 2),"%)"))+ 
  guides(color = guide_legend(override.aes = list(shape = 23, colour = "black", fill = mypalette_age))) +
  theme(legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
        legend.background = element_blank(), legend.box.background =  element_blank())

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_consensus_age_ggplot


#Make hulls for PCA plot with hulls around sex
hulls_all_sex_consensus  <- pcscores_all_consensus  %>%
  group_by(sex) %>%
  slice(chull(comp1, comp2)) %>%
  rename(x = comp1, y = comp2)

#Nice PCA plot with hulls around sex
PCA_all_consensus_sex_ggplot <- ggplot(pcscores_all_consensus[!pcscores_all_consensus$sex%in% c("unknown"),], aes(x = comp1, y = comp2, colour = sex))+
  geom_point(size = 3, aes(shape = vertebra, fill = sex), color = "black") +
  geom_polygon(data = hulls_all_sex_consensus[!hulls_all_sex_consensus$sex%in% c("unknown"),], aes(x = x, y = y, colour = sex, fill = sex), 
               alpha = .3, show.legend = FALSE, inherit.aes = F)+ #colored hulls with transparency
  scale_colour_manual(name = "Sex", labels =  c("Female"   ,    "Male"), #to be ordered as they appear in tibble
                      values= mypalette_sex)+            #legend and color adjustments
  scale_fill_manual(name = "Sex", labels = c( "Female"   ,    "Male"),
                    values =   mypalette_sex)+ #must match scale_colour_manual
  scale_shape_manual(name = "Vertebra", labels = c("Thoracic","Lumbar", "Caudal"), #copy from as.factor(sex)
                     values = shapes)+
  theme_bw()+
  xlab(paste0("PC 1 (",round(PC1_consensus, digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(PC2_consensus, digits = 2),"%)"))+ 
  guides(color = guide_legend(override.aes = list(shape = 23, colour = "black", fill = mypalette_sex))) +
  theme(legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
        legend.background = element_blank(), legend.box.background =  element_blank())

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_consensus_sex_ggplot

#Make hulls for PCA plot with hulls around vertebra type
hulls_all_vertebra_consensus  <- pcscores_all_consensus  %>%
  group_by(vertebra) %>%
  slice(chull(comp1, comp2)) %>%
  rename(x = comp1, y = comp2)

#Nice PCA plot with hulls around vertebra
PCA_all_consensus_vertebra_ggplot <- ggplot(pcscores_all_consensus, aes(x = comp1, y = comp2, colour = vertebra))+
  geom_point(size = 3, aes(shape = age, fill = vertebra), color = "black")+
  geom_polygon(data = hulls_all_vertebra_consensus, aes(x = x, y = y, colour = vertebra, fill = vertebra), 
               alpha = .3, show.legend = FALSE, inherit.aes = F)+ #colored hulls with transparency
  scale_fill_manual(name = "Vertebra", labels = c("Thoracic","Lumbar", "Caudal"),
                    values =  mypalette_vertebrae)+ #must match scale_colour_manual
  scale_colour_manual(name = "Vertebra", labels =  c("Thoracic","Lumbar", "Caudal"), #to be ordered as they appear in tibble
                      values = mypalette_vertebrae)+            #legend and color adjustments
  scale_shape_manual(name = "Age", labels = c( "Adult"  ,  "Juvenile", "Neonate"), #copy from as.factor(age)
                     values = shapes)+
  theme_bw()+
  xlab(paste0("PC 1 (",round(PC1_consensus, digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(PC2_consensus, digits = 2),"%)"))+ 
  guides(color = guide_legend(override.aes = list(shape = 23, colour = "black", fill = mypalette_vertebrae))) +
  theme(legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
        legend.background = element_blank(), legend.box.background =  element_blank())

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_consensus_vertebra_ggplot
