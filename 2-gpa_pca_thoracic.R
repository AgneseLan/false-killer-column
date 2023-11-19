#CH.2 - Prepare final dataset for analysis, run GPA and PCA

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
library(ggthemes)
library(ggpubr)
library(ggplotify)
library(Morpho)
library(rphylopic)
library(png)
library(gridExtra)
library(phytools)
library(evomap)
library(car)
library(Rvcg)

#devtools::install_github("JeroenSmaers/evomap")
#devtools::install_github("wabarr/ggphylomorpho")
#devtools::install_github("aphanotus/borealis")

#DATA PREP ----

###SET WD to root folder from console!! -->

#Import classifiers
classifiers_thoracic <- read_csv("Data/info_thoracic.csv")
glimpse(classifiers_thoracic)

#Make sure the specimens are in the same order before proceeding!!!
identical(dimnames(final_dataset_thoracic)[[3]], classifiers_thoracic$thoracic_vertebrae, attrib.as.set = T)

#Count number of vertebrae per specimen
count_specimens <- classifiers_thoracic %>% count(reg_no_audit)

#Check for outliers in raw data shape array, they would be displayed in red
#Might be due to absent bones, check misstable list
plotOutliers(final_dataset_thoracic)

#Plot outliers landmarks to search for possible problems - check file name to find number in classifiers_thoracic
checkLM(final_dataset_thoracic, path="Data/ply/thoracic_ply/", pt.size = 2, suffix=".ply", render = "s", begin = 1, point = "s")

#Save specimens names as object
specimens_thoracic <- classifiers_thoracic$reg_no_audit

#Save age as factor, useful for later analysis
age_thoracic <- classifiers_thoracic$age
#Check how many colors are needed
as.factor(age_thoracic)

#Save sex as factor, useful for later analysis
sex_thoracic <- classifiers_thoracic$sex
#Check how many colors are needed
as.factor(sex_thoracic)

#Vertebra number
number_thoracic <- classifiers_thoracic$number

#Vertebra Id
Id_thoracic <- classifiers_thoracic$thoracic_vertebrae

#GRAPHICS CODE ----
##Create project palettes----
mypalette_age <- c("#F0E442", "#6699CC", "#AA4499")

mypalette_sex <-  c("#E597B9", "#332288")

shapes <- c(21,22,24)

###Plot landmarks on vertebra ----

##Save mesh with plotted landmarks
#Find mean specimen raw data
findMeanSpec(final_dataset_thoracic)

#Import simplified ply
#Less faces, no holes or isolated triangles
mesh_3D_thoracic <- vcgImport("Data/ply/thoracic_ply/1961_6_14_75_2_merged_c.ply")

#Plot on surface
fixed_LMs_thoracic <- c(1:16)

shade3d(mesh_3D_thoracic, col = "white", alpha = 0.5, fastTransparency = T)
spheres3d(slidedlms_thoracic[fixed_LMs_thoracic,,70], col =  "turquoise4", type = "s",
          radius = 3, aspect = T, main = "mean",axes = F, main = F, fov = 0)
spheres3d(slidedlms_thoracic[-fixed_LMs_thoracic,,70], col =  "turquoise1", type = "s",
          radius = 2, aspect = T, main = "mean",axes = F, main = F, fov = 0)

rgl.snapshot(filename = "Output/landmarks_thoracic_dorsal.png")
rgl.snapshot(filename = "Output/landmarks_thoracic_lateral.png")
rgl.snapshot(filename = "Output/landmarks_thoracic_posterior.png")
rgl.snapshot(filename = "Output/landmarks_thoracic_anterior.png") 

#GPA ALIGNMENT ----

#Procrustes alignment, should also show mean config coordinates
gpa_thoracic <- gpagen(final_dataset_thoracic) 

#Log-transform Centroid size as object
logCsize_thoracic <- log10(gpa_thoracic$Csize) 

#Save mean shape to create links
mean_shape_thoracic <- gpa_thoracic$consensus 

#Coordinates of all specimens after gpa_lumbar alignment
#Removed mirrored landmarks - useful for gpa alignment but can insert error in analyses
coords_thoracic<- gpa_thoracic$coords[c(left.lm_lumbar, midline_lumbar),,]

#Plot all specimens with mean to check that all landmarks are ok
plotAllSpecimens(coords_thoracic, mean = TRUE, label = F, plot.param = list(pt.cex = 0.05, mean.cex = 3, mean.bg = "black"))

#Check for outliers, they would be displayed in red - most immature ones are normal as outliers
plotOutliers(coords_thoracic)

#Plot landmarks from outliers in 3D to check how they look
spheres3d(coords_thoracic[,,10], r = 0.002)

#checkLM(final_dataset_thoracic, path="Data/ply/thoracic_ply/, pt.size = 2, suffix=".ply", render="s", begin = 65) 
#to run if needed to check plotting of points on mesh

##Make data frame for analyses in geomorph
gdf_thoracic <- geomorph.data.frame(coords = coords_thoracic,
                                  sex = as.factor(sex_thoracic), Id = as.factor(Id_thoracic), 
                                  number = as.factor(number_thoracic),
                                  age = as.factor(age_thoracic), specimens = as.factor(specimens_thoracic), 
                                  size = logCsize_thoracic)
glimpse(gdf_thoracic)

#PREPARE WARP MESH  ----

#Find specimen closer to mean, useful to create warp mesh
findMeanSpec(coords_thoracic) #number below specimen name is the number of the specimen in the array

#Create object containing only that specimen coordinates
warp_specimen_thoracic <- coords_thoracic[,,21] #number displayed by findMeanSpec
warp_specimen_thoracic

#Import simplified mesh to create warp mesh on
refmesh_thoracic <- read.ply("Data/ply/thoracic_ply/1961_6_14_245_3_c.ply") #make sure NO binary encoding (ASCII)

#Check range of mesh and coordinates to make sure it has same scale
range(refmesh_thoracic$vb[1:3,]) #if this is too big/small, scale in editor and re-import
range(warp_specimen_thoracic)

#Re-import scaled mesh
#Also clean, reduce in size and remove tecture for faster processing
#Import simplified mesh to create warp mesh on
refmesh_thoracic <- read.ply("Data/ply/refmesh_thoracic.ply") #make sure NO binary encoding (ASCII)

#Check range of mesh and coordinates again to make sure it has same scale
range(refmesh_thoracic$vb[1:3,]) 
range(warp_specimen_thoracic)

#PCA COMPLETE DATASET ----

#Run PCA on complete dataset
PCA_all_thoracic <- gm.prcomp(gdf_thoracic$coords)

#List of PC components and proportion of variation
PCA_all_thoracic

#Save PCA results to file
sink("Output/PCA_all_thoracic_components.txt")
print("PCA complete dataset thoracic")
print(PCA_all_thoracic)
sink() 

#Change row names to codes to make plot readable
row.names(PCA_all_thoracic$x) <- gdf_thoracic$number

##View plot
plot(PCA_all_thoracic, main = "PCA all data thoracic - PC1-PC2",  pch = 21, #title and type of point to be used
     col = "deeppink",   #border of points
     bg = "deeppink",    #fill of points
     cex = 1,            #size of points (1=regular)
     font.main = 2)       #bold font for title
#Add quick labels to plot
text(x = PCA_all_thoracic$x[,1], y = PCA_all_thoracic$x[,2], labels = rownames(PCA_all_thoracic$x), 
     pos = 1,       #position relative to data point
     offset = 0.5,  #distance from data point
     cex = 0.75)    #font size (1=regular)

#Save PC scores as object to use later
pcscores_all_thoracic <- PCA_all_thoracic$x 

#Save shapes of extremes for axes used in plot
PC1min_all_thoracic <- PCA_all_thoracic[["shapes"]][["shapes.comp1"]][["min"]]
PC1max_all_thoracic <- PCA_all_thoracic[["shapes"]][["shapes.comp1"]][["max"]] 
PC2min_all_thoracic <- PCA_all_thoracic[["shapes"]][["shapes.comp2"]][["min"]] 
PC2max_all_thoracic <- PCA_all_thoracic[["shapes"]][["shapes.comp2"]][["max"]] 

#3D deformation using PLY files
#Extremes PC axes

##Create warp mesh, to use as reference for visualization of analyses
PC1min_all_thoracic_mesh <- warpRefMesh(mesh = refmesh_thoracic, mesh.coord = warp_specimen_thoracic, 
                                      ref = PC1min_all_thoracic, color = "gray", centered = FALSE)
##Save warped mesh as ply (Morpho)
mesh2ply(PC1min_all_thoracic_mesh, filename = "Output/PC1min_all_thoracic")

##Create warp mesh, to use as reference for visualization of analyses
PC1max_all_thoracic_mesh <- warpRefMesh(mesh = refmesh_thoracic, mesh.coord = warp_specimen_thoracic, 
                                      ref = PC1max_all_thoracic, color = "gray", centered = FALSE)
##Save warped mesh as ply (Morpho)
mesh2ply(PC1max_all_thoracic_mesh, filename = "Output/PC1max_all_thoracic")

##Create warp mesh, to use as reference for visualization of analyses
PC2min_all_thoracic_mesh <- warpRefMesh(mesh = refmesh_thoracic, mesh.coord = warp_specimen_thoracic, 
                                      ref = PC2min_all_thoracic, color = "gray", centered = FALSE)
##Save warped mesh as ply (Morpho)
mesh2ply(PC2min_all_thoracic_mesh, filename = "Output/PC2min_all_thoracic")

##Create warp mesh, to use as reference for visualization of analyses
PC2max_all_thoracic_mesh <- warpRefMesh(mesh = refmesh_thoracic, mesh.coord = warp_specimen_thoracic, 
                                      ref = PC2max_all_thoracic, color = "gray", centered = FALSE)
##Save warped mesh as ply (Morpho)
mesh2ply(PC2max_all_thoracic_mesh, filename = "Output/PC2max_all_thoracic")

##Make better PCA plot using ggplot
#Read PC scores as tibble
pcscores_all_thoracic <- as_tibble(pcscores_all_thoracic)

#Add labels and other attributes to tibble as columns
pcscores_all_thoracic <- pcscores_all_thoracic %>% mutate(Id = gdf_thoracic$number, age = gdf_thoracic$age,
                                                      sex = gdf_thoracic$sex, specimens = gdf_thoracic$specimens, size = gdf_thoracic$size,
                                                      TL = as.numeric(classifiers_thoracic$size_cm))
glimpse(pcscores_all_thoracic)


#Nice PCA plot with stages and groups
PCA_all_thoracic_ggplot <- ggplot(pcscores_all_thoracic, aes(x = Comp1, y = Comp2, label = Id))+
  geom_point(size = 3, aes(shape =  sex, fill = age), color = "black")+
  geom_text_repel(size = 4, max.overlaps = 40)+
  scale_fill_manual(name = "Age", labels =  c("adult"  ,  "juvenile", "newborn" ), #to be ordered as they appear in tibble
                    values = mypalette_age)+            #legend and color adjustments
  scale_shape_manual(name = "Sex", labels = c( "F"   ,    "M"      , "unknown"), #copy from as.factor(sex)
                     values = shapes)+
  theme_bw()+
  xlab(paste0("PC 1 (",round(as.numeric(PCA_all_thoracic$sdev[1]^2/sum(PCA_all_thoracic$sdev^2)*100), digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(as.numeric(PCA_all_thoracic$sdev[2]^2/sum(PCA_all_thoracic$sdev^2)*100), digits = 2),"%)"))+
  ggtitle("PCA all data thoracic")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  guides(fill = guide_legend(override.aes = list(shape = 23, colour = "black", fill = mypalette_age)))

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_thoracic_ggplot

#Make hulls for PCA plot with hulls around age
hulls_all_age_thoracic <- pcscores_all_thoracic %>%
  group_by(age) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

#Nice PCA plot with hulls around age
PCA_all_thoracic_age_ggplot <- ggplot(pcscores_all_thoracic, aes(x = Comp1, y = Comp2)) +
  geom_point(size = 3, aes(shape = sex, fill = age), color = "black") +
  geom_polygon(data = hulls_all_age_thoracic, aes(x = x, y = y, colour = age, fill = age), 
               alpha = .3, show.legend = FALSE, inherit.aes = FALSE) +
  scale_colour_manual(name = "Age", labels = c("Adult", "Juvenile", "Neonate"),
                      values = mypalette_age) +
  scale_fill_manual(name = "Age", labels = c("Adult", "Juvenile", "Neonate"),
                    values = mypalette_age) +
  scale_shape_manual(name = "Sex", labels = c("Female", "Male", "Neonate"),
                     values = shapes) +
  theme_bw() +
  xlab(paste0("PC 1 (", round(as.numeric(PCA_all_thoracic$sdev[1]^2/sum(PCA_all_thoracic$sdev^2)*100), digits = 2), "%)")) +
  ylab(paste0("PC 2 (", round(as.numeric(PCA_all_thoracic$sdev[2]^2/sum(PCA_all_thoracic$sdev^2)*100), digits = 2), "%)")) +
  guides(color = guide_legend(override.aes = list(shape = 23, colour = "black", fill = mypalette_age))) +
  theme(legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  annotate("text", x = 0.1, y = -0.12, label = "Thoracic",  fontface = 3, size = 12)

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_thoracic_age_ggplot

#Make hulls for PCA plot with hulls around sex
hulls_all_sex_thoracic  <- pcscores_all_thoracic  %>%
  group_by(sex) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

#Nice PCA plot with hulls around sex
PCA_all_thoracic_sex_ggplot <- ggplot(pcscores_all_thoracic[!pcscores_all_thoracic$sex%in% c("unknown"),], aes(x = Comp1, y = Comp2))+
  geom_point(size = 3, aes(shape = age, fill = sex), color = "black") +
  geom_polygon(data = hulls_all_sex_thoracic[!hulls_all_sex_thoracic$sex%in% c("unknown"),], aes(x = x, y = y, colour = sex, fill = sex), 
               alpha = .3, show.legend = FALSE, inherit.aes = F)+ #colored hulls with transparency
  scale_colour_manual(name = "Sex", labels =  c("Female"   ,    "Male"), #to be ordered as they appear in tibble
                      values= mypalette_sex)+            #legend and color adjustments
  scale_fill_manual(name = "Sex", labels = c( "Female"   ,    "Male"),
                    values =   mypalette_sex)+ #must match scale_colour_manual
  scale_shape_manual(name = "Age", labels = c( "Adult"  ,  "Juvenile", "Neonate"), #copy from as.factor(sex)
                     values = shapes)+
  theme_bw()+
  xlab(paste0("PC 1 (",round(as.numeric(PCA_all_thoracic$sdev[1]^2/sum(PCA_all_thoracic$sdev^2)*100), digits = 2),"%)"))+ #copy this from standard PCA plot (PCA_all_plot)
  ylab(paste0("PC 2 (",round(as.numeric(PCA_all_thoracic$sdev[2]^2/sum(PCA_all_thoracic$sdev^2)*100), digits = 2),"%)"))+
  guides(color = guide_legend(override.aes = list(shape = 23, colour = "black", fill = mypalette_sex))) +
  theme(legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  annotate("text", x = 0.1, y = -0.12, label = "Thoracic",  fontface = 3, size = 12)

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_thoracic_sex_ggplot

###Regression PC1 and PC2 ----

#Calculate regression for each component for size
reg_PC1all_thoracic_size <- lm(Comp1 ~ size, data = pcscores_all_thoracic)
reg_PC2all_thoracic_size <- lm(Comp2 ~ size, data = pcscores_all_thoracic)

#View results and p-value
summary(reg_PC1all_thoracic_size)
summary(reg_PC2all_thoracic_size)
anova(reg_PC1all_thoracic_size)
anova(reg_PC2all_thoracic_size)

#Save results of significant regression to file
sink("Output/PC1-2all_thoracic_size_lm.txt")
print("PC1")
summary(reg_PC1all_thoracic_size)
anova(reg_PC1all_thoracic_size)
print("PC2")
summary(reg_PC2all_thoracic_size)
anova(reg_PC2all_thoracic_size)
sink() 

#Calculate regression for each component taking sex into account
reg_PC1all_thoracic_sex <- lm(Comp1 ~ sex, data = pcscores_all_thoracic)
reg_PC2all_thoracic_sex <- lm(Comp2 ~ sex, data = pcscores_all_thoracic)

#View results and p-value
summary(reg_PC1all_thoracic_sex)
summary(reg_PC2all_thoracic_sex)
anova(reg_PC1all_thoracic_sex)
anova(reg_PC2all_thoracic_sex)

#Save results of significant regression to file
sink("Output/PC1-2all_thoracic_sex_lm.txt")
print("PC1")
summary(reg_PC1all_thoracic_sex)
anova(reg_PC1all_thoracic_sex)
print("PC2")
summary(reg_PC2all_thoracic_sex)
anova(reg_PC2all_thoracic_sex)
sink() 

#Calculate regression for each component taking age into account
reg_PC1all_thoracic_age <- lm(Comp1 ~ age, data = pcscores_all_thoracic)
reg_PC2all_thoracic_age <- lm(Comp2 ~ age, data = pcscores_all_thoracic)

#View results and p-value
summary(reg_PC1all_thoracic_age)
summary(reg_PC2all_thoracic_age)
anova(reg_PC1all_thoracic_age)
anova(reg_PC2all_thoracic_age)

#Save results of significant regression to file
sink("Output/PC1-2all_thoracic_age_lm.txt")
print("PC1")
summary(reg_PC1all_thoracic_age)
anova(reg_PC1all_thoracic_age)
print("PC2")
summary(reg_PC2all_thoracic_age)
anova(reg_PC2all_thoracic_age)
sink() 

#Calculate regression for each component taking TL into account
reg_PC1all_thoracic_TL <- lm(Comp1 ~ TL, data = pcscores_all_thoracic)
reg_PC2all_thoracic_TL <- lm(Comp2 ~ TL, data = pcscores_all_thoracic)

#View results and p-value
summary(reg_PC1all_thoracic_TL)
summary(reg_PC2all_thoracic_TL)
anova(reg_PC1all_thoracic_TL)
anova(reg_PC2all_thoracic_TL)

#Save results of significant regression to file
sink("Output/PC1-2all_thoracic_TL_lm.txt")
print("PC1")
summary(reg_PC1all_thoracic_TL)
anova(reg_PC1all_thoracic_TL)
print("PC2")
summary(reg_PC2all_thoracic_TL)
anova(reg_PC2all_thoracic_TL)
sink() 

#Save results of all regressions to 1 file
sink("Output/PC1-2_all_thoracic_lm.txt")
print("PC1")
anova(reg_PC1all_thoracic_size)
anova(reg_PC1all_thoracic_sex)
anova(reg_PC1all_thoracic_age)
anova(reg_PC1all_thoracic_TL)
print("PC2")
anova(reg_PC2all_thoracic_size)
anova(reg_PC2all_thoracic_sex)
anova(reg_PC2all_thoracic_age)
anova(reg_PC2all_thoracic_TL)
sink()


#ANOVA SHAPE AND SIZE vs SEX and AGE -----

##Shape vs Sex ----
#Use only males and females to avoid losing signal due to unknowns

rm_U_thoracic <- which(gdf_thoracic$sex == "unknown")

shape_thoracic_sex <- procD.lm(gdf_thoracic$coords[,,-c(rm_U_thoracic)] ~  as.vector(gdf_thoracic$sex[-rm_U_thoracic]), iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of shape with logCS
summary(shape_thoracic_sex)

#Pairwise comparison for the combination and interaction model
#Helps determine if there is a significant difference in slope (int model) in the size trajectory on top of difference in intercept (comb model)
pairwise_shape_thoracic_sex <- pairwise(shape_thoracic_sex, 
                                      groups =  as.vector(gdf_thoracic$sex[-rm_U_thoracic]), print.progress = FALSE) 
pairwise_shape_thoracic_sex

#Distances between slope vectors (end-points) - absolute difference between slopes of groups
#if significant means int model better than comb
pairwise_shape_thoracic_sex_dist <- summary(pairwise_shape_thoracic_sex, confidence = 0.95, test.type = "dist") 
pairwise_shape_thoracic_sex_dist

#Correlation between slope vectors (and angles) - similarity of vector orientation or angle,
#if significant means the vectors of the groups are oriented in different ways 
pairwise_shape_thoracic_sex_VC <- summary(pairwise_shape_thoracic_sex, confidence = 0.95, test.type = "VC",
                                        angle.type = "deg") 
pairwise_shape_thoracic_sex_VC

#Absolute difference between slope vector lengths - difference in rate of change per covariate unit (size),
#if significant means there is a significant rate of change difference in shape between groups during growth
pairwise_shape_thoracic_sex_DL <-summary(pairwise_shape_thoracic_sex, confidence = 0.95, test.type = "DL") 
pairwise_shape_thoracic_sex_DL

#Save results to file
sink("Output/ANOVA_shape_thoracic_sex.txt")
print("ANOVA")
summary(shape_thoracic_sex)

print("1-Pairwise absolute distances slopes")
summary(pairwise_shape_thoracic_sex, confidence = 0.95, test.type = "dist") 

print("2-Distance between angles (slope directions)")
summary(pairwise_shape_thoracic_sex, confidence = 0.95, test.type = "VC",
        angle.type = "deg") 

print("3-Difference in slope vector length (difference in rate of change of shape per unit of size)")
summary(pairwise_shape_thoracic_sex, confidence = 0.95, test.type = "DL") 
sink()

##Shape vs Age ----
shape_thoracic_age <- procD.lm(gdf_thoracic$coords ~ gdf_thoracic$age, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of shape with logCS
summary(shape_thoracic_age)

#Pairwise comparison for the combination and interaction model
#Helps determine if there is a significant difference in slope (int model) in the size trajectory on top of difference in intercept (comb model)
pairwise_shape_thoracic_age <- pairwise(shape_thoracic_age, 
                                      groups = gdf_thoracic$age, print.progress = FALSE) 
pairwise_shape_thoracic_age

#Distances between slope vectors (end-points) - absolute difference between slopes of groups
#if significant means int model better than comb
pairwise_shape_thoracic_age_dist <- summary(pairwise_shape_thoracic_age, confidence = 0.95, test.type = "dist") 
pairwise_shape_thoracic_age_dist

#Correlation between slope vectors (and angles) - similarity of vector orientation or angle,
#if significant means the vectors of the groups are oriented in different ways 
pairwise_shape_thoracic_age_VC <- summary(pairwise_shape_thoracic_age, confidence = 0.95, test.type = "VC",
                                        angle.type = "deg") 
pairwise_shape_thoracic_age_VC

#Absolute difference between slope vector lengths - difference in rate of change per covariate unit (size),
#if significant means there is a significant rate of change difference in shape between groups during growth
pairwise_shape_thoracic_age_DL <-summary(pairwise_shape_thoracic_age, confidence = 0.95, test.type = "DL") 
pairwise_shape_thoracic_age_DL

#Save results to file
sink("Output/ANOVA_shape_thoracic_age.txt")
print("ANOVA")
summary(shape_thoracic_age)

print("1-Pairwise absolute distances slopes")
summary(pairwise_shape_thoracic_age, confidence = 0.95, test.type = "dist") 

print("2-Distance between angles (slope directions)")
summary(pairwise_shape_thoracic_age, confidence = 0.95, test.type = "VC",
        angle.type = "deg") 

print("3-Difference in slope vector length (difference in rate of change of shape per unit of size)")
summary(pairwise_shape_thoracic_age, confidence = 0.95, test.type = "DL") 
sink()


##Size vs Sex ----
#Use only males and females to avoid losing signal due to unknowns
size_thoracic_sex <- lm.rrpp(gdf_thoracic$size[-rm_U_thoracic] ~ as.vector(gdf_thoracic$sex[-rm_U_thoracic]), iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of size with logCS
summary(size_thoracic_sex)

#Pairwise comparison for the combination and interaction model
#Helps determine if there is a significant difference in slope (int model) in the size trajectory on top of difference in intercept (comb model)
pairwise_size_thoracic_sex <- pairwise(size_thoracic_sex, 
                                     groups =   as.vector(gdf_thoracic$sex[-rm_U_thoracic]), print.progress = FALSE) 
pairwise_size_thoracic_sex

#Distances between slope vectors (end-points) - absolute difference between slopes of groups
#if significant means int model better than comb
pairwise_size_thoracic_sex_dist <- summary(pairwise_size_thoracic_sex, confidence = 0.95, test.type = "dist") 
pairwise_size_thoracic_sex_dist

#Correlation between slope vectors (and angles) - similarity of vector orientation or angle,
#if significant means the vectors of the groups are oriented in different ways 
pairwise_size_thoracic_sex_VC <- summary(pairwise_size_thoracic_sex, confidence = 0.95, test.type = "VC",
                                       angle.type = "deg") 
pairwise_size_thoracic_sex_VC

#Absolute difference between slope vector lengths - difference in rate of change per covariate unit (size),
#if significant means there is a significant rate of change difference in size between groups during growth
pairwise_size_thoracic_sex_DL <-summary(pairwise_size_thoracic_sex, confidence = 0.95, test.type = "DL") 
pairwise_size_thoracic_sex_DL

#Save results to file
sink("Output/ANOVA_size_thoracic_sex.txt")
print("ANOVA")
summary(size_thoracic_sex)

print("1-Pairwise absolute distances slopes")
summary(pairwise_size_thoracic_sex, confidence = 0.95, test.type = "dist") 

print("2-Distance between angles (slope directions)")
summary(pairwise_size_thoracic_sex, confidence = 0.95, test.type = "VC",
        angle.type = "deg") 

print("3-Difference in slope vector length (difference in rate of change of size per unit of size)")
summary(pairwise_size_thoracic_sex, confidence = 0.95, test.type = "DL") 
sink()

##Size vs Age ----
size_thoracic_age <- lm.rrpp(gdf_thoracic$size ~ gdf_thoracic$age, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of size with logCS
summary(size_thoracic_age)

#Pairwise comparison for the combination and interaction model
#Helps determine if there is a significant difference in slope (int model) in the size trajectory on top of difference in intercept (comb model)
pairwise_size_thoracic_age <- pairwise(size_thoracic_age, 
                                     groups = gdf_thoracic$age, print.progress = FALSE) 
pairwise_size_thoracic_age

#Distances between slope vectors (end-points) - absolute difference between slopes of groups
#if significant means int model better than comb
pairwise_size_thoracic_age_dist <- summary(pairwise_size_thoracic_age, confidence = 0.95, test.type = "dist") 
pairwise_size_thoracic_age_dist

#Correlation between slope vectors (and angles) - similarity of vector orientation or angle,
#if significant means the vectors of the groups are oriented in different ways 
pairwise_size_thoracic_age_VC <- summary(pairwise_size_thoracic_age, confidence = 0.95, test.type = "VC",
                                       angle.type = "deg") 
pairwise_size_thoracic_age_VC

#Absolute difference between slope vector lengths - difference in rate of change per covariate unit (size),
#if significant means there is a significant rate of change difference in size between groups during growth
pairwise_size_thoracic_age_DL <-summary(pairwise_size_thoracic_age, confidence = 0.95, test.type = "DL") 
pairwise_size_thoracic_age_DL

#Save results to file
sink("Output/ANOVA_size_thoracic_age.txt")
print("ANOVA")
summary(size_thoracic_age)

print("1-Pairwise absolute distances slopes")
summary(pairwise_size_thoracic_age, confidence = 0.95, test.type = "dist") 

print("2-Distance between angles (slope directions)")
summary(pairwise_size_thoracic_age, confidence = 0.95, test.type = "VC",
        angle.type = "deg") 

print("3-Difference in slope vector length (difference in rate of change of size per unit of size)")
summary(pairwise_size_thoracic_age, confidence = 0.95, test.type = "DL") 
sink()


#PHENOTYPIC TRAJECTORY SEXUAL DIMORPHISM ---- 

#Make sure to only use males and females!

#Shape changes along column at different ages
fit_shape_thoracic_age <- procD.lm(gdf_thoracic$coords[,,-c(rm_U_thoracic)] ~ as.vector(gdf_thoracic$sex[-rm_U_thoracic]) * gdf_thoracic$age[-rm_U_thoracic], iter = 999, RRPP = F)


#Check that there is a significant correlation
summary(fit_shape_thoracic_age)

#Use fit to calculate trajectories
trajectory_thoracic_age <- trajectory.analysis(fit_shape_thoracic_age, groups = as.vector(gdf_thoracic$sex[-rm_U_thoracic]), traj.pts = gdf_thoracic$age[-rm_U_thoracic], 
                                             pca = TRUE) 

#View results
#Magnitude differences between trajectories, standard summary - are trajectories different in length?
trajectory_thoracic_age_MD <- summary(trajectory_thoracic_age, show.trajectories = TRUE, attribute = "MD") 
trajectory_thoracic_age_MD 
#Trajectory correlations -  are trajectories different in angle/direction?
trajectory_thoracic_age_TC <- summary(trajectory_thoracic_age, show.trajectories = TRUE, attribute = "TC", angle.type = "deg")
trajectory_thoracic_age_TC
#Trajectory shape differences - are trajectories different in shape?
trajectory_thoracic_age_SD <- summary(trajectory_thoracic_age, show.trajectories = TRUE, attribute = "SD") 
trajectory_thoracic_age_SD

#Save results to file
sink("Output/trajectory_thoracic_age.txt")
print("Initial fit shape ~ sex * age")
summary(fit_shape_thoracic_age)
print("Magnitude difference (absolute difference between path distances) - length")
summary(trajectory_thoracic_age, show.trajectories = TRUE, attribute = "MD") 
print("Correlations (angles) between trajectories - direction")
summary(trajectory_thoracic_age, show.trajectories = TRUE, attribute = "TC", angle.type = "deg")
print("Shape differences between trajectory vectors - shape")
summary(trajectory_thoracic_age, show.trajectories = TRUE, attribute = "SD") 
sink() 


#Plot results - PCA of fitted values
trajectory_thoracic_age_plot <- plot(trajectory_thoracic_age, main = "Trajectories shape change during growth by sex", pch = shapes[c(1,2)],  #title and type of point to be used
                                   col = c("darkgray","lightgray"), bg = c("darkgray","lightgray"), cex = 1, font.main = 2) #improve graphics
#Add line between groups
add.trajectories(trajectory_thoracic_age_plot, traj.pch = shapes[c(1,2)],
                 traj.col = 1, traj.lty = 1, traj.lwd = 1, traj.cex = 1.5, traj.bg = 1, 
                 start.bg = "green", end.bg = "red") #trajectory line graphics
#Add legend to see which trajectory belongs to each group
legend("bottomleft", pch = shapes[c(1,2)],legend = c("Female"  ,  "Male"), pt.bg = 1, cex = 1)

##Make better PCA plot using ggplot
#Read PC scores as tibble
trajectory_thoracic_age_pcscores <- as_tibble(trajectory_thoracic_age_plot [["pc.points"]])

#Add group names and other attributes to tibble as columns
trajectory_thoracic_age_pcscores <- trajectory_thoracic_age_pcscores %>% mutate(sex = gdf_thoracic$sex[-rm_U_thoracic], age = gdf_thoracic$age[-rm_U_thoracic])
glimpse(trajectory_thoracic_age_pcscores)

#Calculate means of PC1 and PC2 at each vertebra per group to draw trajectories
trajectory_thoracic_age_pcscores_means <- trajectory_thoracic_age_pcscores %>% group_by(age, sex) %>%
  summarise_at(vars(PC1, PC2), list(mean = mean))              #function for means, both columns
glimpse(trajectory_thoracic_age_pcscores_means)

#Rename columns so they are easier to use for plot
trajectory_thoracic_age_pcscores_means <- trajectory_thoracic_age_pcscores_means %>% rename(x = PC1_mean, y = PC2_mean)
glimpse(trajectory_thoracic_age_pcscores_means)

#Nice trajectory plot by sex and age
trajectory_thoracic_age_ggplot <- ggplot(trajectory_thoracic_age_pcscores, aes(x = PC1, y = PC2, shape = age, group = sex))+
  geom_point(size = 2.2, colour = "black", fill = "lightgray", alpha = 0.5, show.legend = F)+
  geom_point(data = trajectory_thoracic_age_pcscores_means, aes(x = x, y = y, fill = sex, shape = age, group = sex), colour = "black",
             size = 5, inherit.aes = F, alpha = 0.8)+
  geom_path(data = trajectory_thoracic_age_pcscores_means, aes(x = x, y = y, colour = sex, group = sex), inherit.aes = F, linewidth = 1,
            linejoin = 'mitre', arrow = arrow(angle = 40, length = unit(0.03, "npc"), ends = "first", type = "open"), show.legend = F)+
  scale_colour_manual(name = "Sex", labels = c("Female", "Male"), values = mypalette_sex, aesthetics = c("colour", "fill"))+
  scale_shape_manual(name = "Age", labels = c("Adult"  ,  "Juvenile", "Neonate"), #copy from as.factor(genera)
                     values = shapes)+
  theme_bw()+
  xlab(paste0("PC 1 (",round(as.numeric(trajectory_thoracic_age[["pca"]][["sdev"]][1]^2/sum(trajectory_thoracic_age[["pca"]][["sdev"]]^2)*100), digits = 2),"%)"))+ #copy this from standard trajectory plot
  ylab(paste0("PC 2 (",round(as.numeric(trajectory_thoracic_age[["pca"]][["sdev"]][2]^2/sum(trajectory_thoracic_age[["pca"]][["sdev"]]^2)*100), digits = 2),"%)"))+
  theme(legend.key = element_blank(), legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  guides(color = guide_legend(override.aes = list(shape = 23, colour = "black", fill = mypalette_sex)))+
  annotate("text", x = 0.09, y = -0.062, label = "Thoracic",  fontface = 3, size = 12)

#Visualize plot and save as PDF using menu in bar on the right
trajectory_thoracic_age_ggplot


