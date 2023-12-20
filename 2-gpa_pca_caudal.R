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
classifiers_caudal <- read_csv("Data/info_caudal.csv")
glimpse(classifiers_caudal)

#Make sure the specimens are in the same order before proceeding!!!
identical(dimnames(final_dataset_caudal)[[3]], classifiers_caudal$caudal_vertebrae, attrib.as.set = T)

#Count number of vertebrae per specimen
count_specimens <- classifiers_caudal %>% count(reg_no_audit)

#Check for outliers in raw data shape array, they would be displayed in red
#Might be due to absent bones, check misstable list
plotOutliers(final_dataset_caudal)

#Plot outliers landmarks to search for possible problems - check file name to find number in classifiers_caudal
checkLM(final_dataset_caudal, path="Data/ply/caudal_ply/", pt.size = 2, suffix=".ply", render = "s", begin = 1, point = "s")

#Save specimens names as object
specimens_caudal <- classifiers_caudal$reg_no_audit

#Save age as factor, useful for later analysis
age_caudal <- classifiers_caudal$age
#Check how many colors are needed
as.factor(age_caudal)

#Save sex as factor, useful for later analysis
sex_caudal <- classifiers_caudal$sex
#Check how many colors are needed
as.factor(sex_caudal)

#Vertebra number
number_caudal <- classifiers_caudal$number

#Vertebra Id
Id_caudal <- classifiers_caudal$caudal_vertebrae


#GRAPHICS CODE ----
##Create project palettes----
mypalette_age <- c("#F0E442", "#6699CC", "#AA4499")

mypalette_sex <-  c("#E597B9", "#332288")

shapes <- c(21,22,24)


###Plot landmarks on vertebra ----

##Save mesh with plotted landmarks
#Find mean specimen raw data
findMeanSpec(final_dataset_caudal)

#Import simplified ply
#Less faces, no holes or isolated triangles
mesh_3D_caudal <- vcgImport("Data/ply/caudal_ply/1961_6_14_22_3_e.ply")

#Plot on surface
fixed_LMs_caudal <- c(1:13)

shade3d(mesh_3D_caudal, col = "white", alpha = 0.5, fastTransparency = T)
spheres3d(slidedlms_caudal[fixed_LMs_caudal,,30], col =  "turquoise4", type = "s",
          radius = 3, aspect = T, main = "mean",axes = F, main = F, fov = 0)
spheres3d(slidedlms_caudal[-fixed_LMs_caudal,,30], col =  "turquoise1", type = "s",
          radius = 2, aspect = T, main = "mean",axes = F, main = F, fov = 0)

rgl.snapshot(filename = "Output/landmarks_caudal_dorsal.png")
rgl.snapshot(filename = "Output/landmarks_caudal_lateral.png")
rgl.snapshot(filename = "Output/landmarks_caudal_posterior.png")
rgl.snapshot(filename = "Output/landmarks_caudal_anterior.png") 


#GPA ALIGNMENT ----

#Procrustes alignment, should also show mean config coordinates
gpa_caudal <- gpagen(final_dataset_caudal) 

#Log-transform Centroid size as object
logCsize_caudal <- log10(gpa_caudal$Csize) 

#Save mean shape to create links
mean_shape_caudal <- gpa_caudal$consensus 

#Coordinates of all specimens after gpa_lumbar alignment
#Removed mirrored landmarks - useful for gpa alignment but can insert error in analyses
coords_caudal<- gpa_caudal$coords[c(left.lm_lumbar, midline_lumbar),,]

#Plot all specimens with mean to check that all landmarks are ok
plotAllSpecimens(coords_caudal, mean = TRUE, label = F, plot.param = list(pt.cex = 0.05, mean.cex = 3, mean.bg = "black"))

#Check for outliers, they would be displayed in red - most immature ones are normal as outliers
plotOutliers(coords_caudal)

#Plot landmarks from outliers in 3D to check how they look
spheres3d(coords_caudal[,,133], r = 0.002)

#checkLM(final_dataset_caudal, path="Data/ply/caudal_ply/, pt.size = 2, suffix=".ply", render="s", begin = 65) 
#to run if needed to check plotting of points on mesh

##Make data frame for analyses in geomorph
gdf_caudal <- geomorph.data.frame(coords = coords_caudal,
                                  sex = as.factor(sex_caudal), Id = as.factor(Id_caudal), 
                                  number = as.factor(number_caudal),
                                  age = as.factor(age_caudal), specimens = as.factor(specimens_caudal), 
                                  size = logCsize_caudal)
glimpse(gdf_caudal)

#PREPARE WARP MESH  ----

#Find specimen closer to mean, useful to create warp mesh

findMeanSpec(coords_caudal) #number below specimen name is the number of the specimen in the array

#Create object containing only that specimen coordinates
warp_specimen_caudal <- coords_caudal[,,38] #number displayed by findMeanSpec
warp_specimen_caudal

#Import simplified mesh to create warp mesh on
refmesh_caudal <- read.ply("Data/ply/caudal_ply/1961_6_14_247_4_b.ply") #make sure NO binary encoding (ASCII)

#Check range of mesh and coordinates to make sure it has same scale
range(refmesh_caudal$vb[1:3,]) #if this is too big/small, scale in editor and re-import
range(warp_specimen_caudal)

#Re-import scaled mesh
#Also clean, reduce in size and remove tecture for faster processing
#Import simplified mesh to create warp mesh on
refmesh_caudal <- read.ply("Data/ply/refmesh_caudal.ply") #make sure NO binary encoding (ASCII)

#Check range of mesh and coordinates again to make sure it has same scale
range(refmesh_caudal$vb[1:3,]) 
range(warp_specimen_caudal)

#PCA COMPLETE DATASET ----

#Run PCA on complete dataset
PCA_all_caudal <- gm.prcomp(gdf_caudal$coords)

#List of PC components and proportion of variation
PCA_all_caudal

#Save PCA results to file
sink("Output/PCA_all_caudal_components.txt")
print("PCA complete dataset caudal")
print(PCA_all_caudal)
sink() 

#Change row names to codes to make plot readable
row.names(PCA_all_caudal$x) <- gdf_caudal$number

##View plot
plot(PCA_all_caudal, main = "PCA all data caudal - PC1-PC2",  pch = 21, #title and type of point to be used
     col = "deeppink",   #border of points
     bg = "deeppink",    #fill of points
     cex = 1,            #size of points (1=regular)
     font.main = 2)       #bold font for title
#Add quick labels to plot
text(x = PCA_all_caudal$x[,1], y = PCA_all_caudal$x[,2], labels = rownames(PCA_all_caudal$x), 
     pos = 1,       #position relative to data point
     offset = 0.5,  #distance from data point
     cex = 0.75)    #font size (1=regular)

#Save PC scores as object to use later
pcscores_all_caudal <- PCA_all_caudal$x 

#Save shapes of extremes for axes used in plot
PC1min_all_caudal <- PCA_all_caudal[["shapes"]][["shapes.comp1"]][["min"]]
PC1max_all_caudal <- PCA_all_caudal[["shapes"]][["shapes.comp1"]][["max"]] 
PC2min_all_caudal <- PCA_all_caudal[["shapes"]][["shapes.comp2"]][["min"]] 
PC2max_all_caudal <- PCA_all_caudal[["shapes"]][["shapes.comp2"]][["max"]] 

#3D deformation using PLY files
#Extremes PC axes

##Create warp mesh, to use as reference for visualization of analyses
PC1min_all_caudal_mesh <- warpRefMesh(mesh = refmesh_caudal, mesh.coord = warp_specimen_caudal, 
                                        ref = PC1min_all_caudal, color = "gray", centered = FALSE)
##Save warped mesh as ply (Morpho)
mesh2ply(PC1min_all_caudal_mesh, filename = "Output/PC1min_all_caudal")

##Create warp mesh, to use as reference for visualization of analyses
PC1max_all_caudal_mesh <- warpRefMesh(mesh = refmesh_caudal, mesh.coord = warp_specimen_caudal, 
                                        ref = PC1max_all_caudal, color = "gray", centered = FALSE)
##Save warped mesh as ply (Morpho)
mesh2ply(PC1max_all_caudal_mesh, filename = "Output/PC1max_all_caudal")

##Create warp mesh, to use as reference for visualization of analyses
PC2min_all_caudal_mesh <- warpRefMesh(mesh = refmesh_caudal, mesh.coord = warp_specimen_caudal, 
                                        ref = PC2min_all_caudal, color = "gray", centered = FALSE)
##Save warped mesh as ply (Morpho)
mesh2ply(PC2min_all_caudal_mesh, filename = "Output/PC2min_all_caudal")

##Create warp mesh, to use as reference for visualization of analyses
PC2max_all_caudal_mesh <- warpRefMesh(mesh = refmesh_caudal, mesh.coord = warp_specimen_caudal, 
                                        ref = PC2max_all_caudal, color = "gray", centered = FALSE)
##Save warped mesh as ply (Morpho)
mesh2ply(PC2max_all_caudal_mesh, filename = "Output/PC2max_all_caudal")

##Make better PCA plot using ggplot
#Read PC scores as tibble
pcscores_all_caudal <- as_tibble(pcscores_all_caudal)

#Add labels and other attributes to tibble as columns
pcscores_all_caudal <- pcscores_all_caudal %>% mutate(Id = gdf_caudal$number, age = gdf_caudal$age,
                                                          sex = gdf_caudal$sex, specimens = gdf_caudal$specimens, size = gdf_caudal$size,
                                                          TL = as.numeric(classifiers_caudal$size_cm))
glimpse(pcscores_all_caudal)


#Nice PCA plot with stages and groups
PCA_all_caudal_ggplot <- ggplot(pcscores_all_caudal, aes(x = Comp1, y = Comp2, label = Id))+
  geom_point(size = 3, aes(shape =  sex, fill = age), color = "black")+
  geom_text_repel(size = 4, max.overlaps = 40)+
  scale_fill_manual(name = "Age", labels =  c("adult"  ,  "juvenile", "newborn" ), #to be ordered as they appear in tibble
                      values = mypalette_age)+            #legend and color adjustments
  scale_shape_manual(name = "Sex", labels = c( "F"   ,    "M"      , "unknown"), #copy from as.factor(sex)
                     values = shapes)+
  theme_bw(base_size = 14)+
  xlab(paste0("PC 1 (",round(as.numeric(PCA_all_caudal$sdev[1]^2/sum(PCA_all_caudal$sdev^2)*100), digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(as.numeric(PCA_all_caudal$sdev[2]^2/sum(PCA_all_caudal$sdev^2)*100), digits = 2),"%)"))+
  ggtitle("PCA all data caudal")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  guides(fill = guide_legend(override.aes = list(shape = 23, colour = "black", fill = mypalette_age)))

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_caudal_ggplot

#Make hulls for PCA plot with hulls around age
hulls_all_age_caudal <- pcscores_all_caudal %>%
  group_by(age) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

#Nice PCA plot with hulls around age
PCA_all_caudal_age_ggplot <- ggplot(pcscores_all_caudal, aes(x = Comp1, y = Comp2)) +
  geom_point(size = 3, aes(shape = sex, fill = age), color = "black") +
  geom_polygon(data = hulls_all_age_caudal, aes(x = x, y = y, colour = age, fill = age), 
               alpha = .3, show.legend = FALSE, inherit.aes = FALSE) +
  scale_colour_manual(name = "Age", labels = c("Adult", "Juvenile", "Neonate"),
                      values = mypalette_age) +
  scale_fill_manual(name = "Age", labels = c("Adult", "Juvenile", "Neonate"),
                    values = mypalette_age) +
  scale_shape_manual(name = "Sex", labels = c("Female", "Male", "Neonate"),
                     values = shapes) +
  theme_bw(base_size = 14) +
  xlab(paste0("PC 1 (", round(as.numeric(PCA_all_caudal$sdev[1]^2/sum(PCA_all_caudal$sdev^2)*100), digits = 2), "%)")) +
  ylab(paste0("PC 2 (", round(as.numeric(PCA_all_caudal$sdev[2]^2/sum(PCA_all_caudal$sdev^2)*100), digits = 2), "%)")) +
  guides(color = guide_legend(override.aes = list(shape = 23, colour = "black", fill = mypalette_age))) +
  theme(legend.key = element_blank(), legend.title = element_text(size = 15, face = "bold"), legend.text = element_text(size = 14),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  annotate("text", x = 0.15, y = -0.12, label = "Caudal",  fontface = 3, size = 12)

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_caudal_age_ggplot

#Make hulls for PCA plot with hulls around sex
hulls_all_sex_caudal  <- pcscores_all_caudal  %>%
  group_by(sex) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

#Nice PCA plot with hulls around sex
PCA_all_caudal_sex_ggplot <- ggplot(pcscores_all_caudal[!pcscores_all_caudal$sex%in% c("unknown"),], aes(x = Comp1, y = Comp2))+
  geom_point(size = 3, aes(shape = age, fill = sex), color = "black") +
  geom_polygon(data = hulls_all_sex_caudal[!hulls_all_sex_caudal$sex%in% c("unknown"),], aes(x = x, y = y, colour = sex, fill = sex), 
               alpha = .3, show.legend = FALSE, inherit.aes = F)+ #colored hulls with transparency
  scale_colour_manual(name = "Sex", labels =  c("Female"   ,    "Male"), #to be ordered as they appear in tibble
                      values= mypalette_sex)+            #legend and color adjustments
  scale_fill_manual(name = "Sex", labels = c( "Female"   ,    "Male"),
                    values =   mypalette_sex)+ #must match scale_colour_manual
  scale_shape_manual(name = "Age", labels = c( "Adult"  ,  "Juvenile", "Neonate"), #copy from as.factor(sex)
                     values = shapes)+
  theme_bw(base_size = 14)+
  xlab(paste0("PC 1 (",round(as.numeric(PCA_all_caudal$sdev[1]^2/sum(PCA_all_caudal$sdev^2)*100), digits = 2),"%)"))+ #copy this from standard PCA plot (PCA_all_plot)
  ylab(paste0("PC 2 (",round(as.numeric(PCA_all_caudal$sdev[2]^2/sum(PCA_all_caudal$sdev^2)*100), digits = 2),"%)"))+
  guides(color = guide_legend(override.aes = list(shape = 23, colour = "black", fill = mypalette_sex))) +
  theme(legend.key = element_blank(), legend.title = element_text(size = 15, face = "bold"), legend.text = element_text(size = 14),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  annotate("text", x = 0.15, y = -0.12, label = "Caudal",  fontface = 3, size = 12)

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_caudal_sex_ggplot

###Regression PC1 and PC2 ----

#Calculate regression for each component for size
reg_PC1all_caudal_size <- lm(Comp1 ~ size, data = pcscores_all_caudal)
reg_PC2all_caudal_size <- lm(Comp2 ~ size, data = pcscores_all_caudal)

#View results and p-value
summary(reg_PC1all_caudal_size)
summary(reg_PC2all_caudal_size)
anova(reg_PC1all_caudal_size)
anova(reg_PC2all_caudal_size)

#Save results of significant regression to file
sink("Output/PC1-2all_caudal_size_lm.txt")
print("PC1")
summary(reg_PC1all_caudal_size)
anova(reg_PC1all_caudal_size)
print("PC2")
summary(reg_PC2all_caudal_size)
anova(reg_PC2all_caudal_size)
sink() 

#Calculate regression for each component taking sex into account
reg_PC1all_caudal_sex <- lm(Comp1 ~ sex, data = pcscores_all_caudal)
reg_PC2all_caudal_sex <- lm(Comp2 ~ sex, data = pcscores_all_caudal)

#View results and p-value
summary(reg_PC1all_caudal_sex)
summary(reg_PC2all_caudal_sex)
anova(reg_PC1all_caudal_sex)
anova(reg_PC2all_caudal_sex)

#Save results of significant regression to file
sink("Output/PC1-2all_caudal_sex_lm.txt")
print("PC1")
summary(reg_PC1all_caudal_sex)
anova(reg_PC1all_caudal_sex)
print("PC2")
summary(reg_PC2all_caudal_sex)
anova(reg_PC2all_caudal_sex)
sink() 

#Calculate regression for each component taking age into account
reg_PC1all_caudal_age <- lm(Comp1 ~ age, data = pcscores_all_caudal)
reg_PC2all_caudal_age <- lm(Comp2 ~ age, data = pcscores_all_caudal)

#View results and p-value
summary(reg_PC1all_caudal_age)
summary(reg_PC2all_caudal_age)
anova(reg_PC1all_caudal_age)
anova(reg_PC2all_caudal_age)

#Save results of significant regression to file
sink("Output/PC1-2all_caudal_age_lm.txt")
print("PC1")
summary(reg_PC1all_caudal_age)
anova(reg_PC1all_caudal_age)
print("PC2")
summary(reg_PC2all_caudal_age)
anova(reg_PC2all_caudal_age)
sink() 

#Calculate regression for each component taking TL into account
reg_PC1all_caudal_TL <- lm(Comp1 ~ TL, data = pcscores_all_caudal)
reg_PC2all_caudal_TL <- lm(Comp2 ~ TL, data = pcscores_all_caudal)

#View results and p-value
summary(reg_PC1all_caudal_TL)
summary(reg_PC2all_caudal_TL)
anova(reg_PC1all_caudal_TL)
anova(reg_PC2all_caudal_TL)

#Save results of significant regression to file
sink("Output/PC1-2all_caudal_TL_lm.txt")
print("PC1")
summary(reg_PC1all_caudal_TL)
anova(reg_PC1all_caudal_TL)
print("PC2")
summary(reg_PC2all_caudal_TL)
anova(reg_PC2all_caudal_TL)
sink() 

#Save results of all regressions to 1 file
sink("Output/PC1-2_all_caudal_lm.txt")
print("PC1")
anova(reg_PC1all_caudal_size)
anova(reg_PC1all_caudal_sex)
anova(reg_PC1all_caudal_age)
anova(reg_PC1all_caudal_TL)
print("PC2")
anova(reg_PC2all_caudal_size)
anova(reg_PC2all_caudal_sex)
anova(reg_PC2all_caudal_age)
anova(reg_PC2all_caudal_TL)
sink()


#ANOVA SHAPE AND SIZE vs SEX and AGE -----

##Shape vs Sex ----
#Use only males and females to avoid losing signal due to unknowns

rm_U_caudal <- which(gdf_caudal$sex == "unknown")

shape_caudal_sex <- procD.lm(gdf_caudal$coords[,,-c(rm_U_caudal)] ~  as.vector(gdf_caudal$sex[-rm_U_caudal]), iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of shape with logCS
summary(shape_caudal_sex)

#Pairwise comparison for the combination and interaction model
#Helps determine if there is a significant difference in slope (int model) in the size trajectory on top of difference in intercept (comb model)
pairwise_shape_caudal_sex <- pairwise(shape_caudal_sex, 
                                   groups =  as.vector(gdf_caudal$sex[-rm_U_caudal]), print.progress = FALSE) 
pairwise_shape_caudal_sex

#Distances between slope vectors (end-points) - absolute difference between slopes of groups
#if significant means int model better than comb
pairwise_shape_caudal_sex_dist <- summary(pairwise_shape_caudal_sex, confidence = 0.95, test.type = "dist") 
pairwise_shape_caudal_sex_dist

#Correlation between slope vectors (and angles) - similarity of vector orientation or angle,
#if significant means the vectors of the groups are oriented in different ways 
pairwise_shape_caudal_sex_VC <- summary(pairwise_shape_caudal_sex, confidence = 0.95, test.type = "VC",
                                     angle.type = "deg") 
pairwise_shape_caudal_sex_VC

#Absolute difference between slope vector lengths - difference in rate of change per covariate unit (size),
#if significant means there is a significant rate of change difference in shape between groups during growth
pairwise_shape_caudal_sex_DL <-summary(pairwise_shape_caudal_sex, confidence = 0.95, test.type = "DL") 
pairwise_shape_caudal_sex_DL

#Save results to file
sink("Output/ANOVA_shape_caudal_sex.txt")
print("ANOVA")
summary(shape_caudal_sex)

print("1-Pairwise absolute distances slopes")
summary(pairwise_shape_caudal_sex, confidence = 0.95, test.type = "dist") 

print("2-Distance between angles (slope directions)")
summary(pairwise_shape_caudal_sex, confidence = 0.95, test.type = "VC",
        angle.type = "deg") 

print("3-Difference in slope vector length (difference in rate of change of shape per unit of size)")
summary(pairwise_shape_caudal_sex, confidence = 0.95, test.type = "DL") 
sink()

##Shape vs Age ----
shape_caudal_age <- procD.lm(gdf_caudal$coords ~ gdf_caudal$age, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of shape with logCS
summary(shape_caudal_age)

#Pairwise comparison for the combination and interaction model
#Helps determine if there is a significant difference in slope (int model) in the size trajectory on top of difference in intercept (comb model)
pairwise_shape_caudal_age <- pairwise(shape_caudal_age, 
                                      groups = gdf_caudal$age, print.progress = FALSE) 
pairwise_shape_caudal_age

#Distances between slope vectors (end-points) - absolute difference between slopes of groups
#if significant means int model better than comb
pairwise_shape_caudal_age_dist <- summary(pairwise_shape_caudal_age, confidence = 0.95, test.type = "dist") 
pairwise_shape_caudal_age_dist

#Correlation between slope vectors (and angles) - similarity of vector orientation or angle,
#if significant means the vectors of the groups are oriented in different ways 
pairwise_shape_caudal_age_VC <- summary(pairwise_shape_caudal_age, confidence = 0.95, test.type = "VC",
                                        angle.type = "deg") 
pairwise_shape_caudal_age_VC

#Absolute difference between slope vector lengths - difference in rate of change per covariate unit (size),
#if significant means there is a significant rate of change difference in shape between groups during growth
pairwise_shape_caudal_age_DL <-summary(pairwise_shape_caudal_age, confidence = 0.95, test.type = "DL") 
pairwise_shape_caudal_age_DL

#Save results to file
sink("Output/ANOVA_shape_caudal_age.txt")
print("ANOVA")
summary(shape_caudal_age)

print("1-Pairwise absolute distances slopes")
summary(pairwise_shape_caudal_age, confidence = 0.95, test.type = "dist") 

print("2-Distance between angles (slope directions)")
summary(pairwise_shape_caudal_age, confidence = 0.95, test.type = "VC",
        angle.type = "deg") 

print("3-Difference in slope vector length (difference in rate of change of shape per unit of size)")
summary(pairwise_shape_caudal_age, confidence = 0.95, test.type = "DL") 
sink()


##Size vs Sex ----
#Use only males and females to avoid losing signal due to unknowns
size_caudal_sex <- lm.rrpp(gdf_caudal$size[-rm_U_caudal] ~ as.vector(gdf_caudal$sex[-rm_U_caudal]), iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of size with logCS
summary(size_caudal_sex)

#Pairwise comparison for the combination and interaction model
#Helps determine if there is a significant difference in slope (int model) in the size trajectory on top of difference in intercept (comb model)
pairwise_size_caudal_sex <- pairwise(size_caudal_sex, 
                                      groups =   as.vector(gdf_caudal$sex[-rm_U_caudal]), print.progress = FALSE) 
pairwise_size_caudal_sex

#Distances between slope vectors (end-points) - absolute difference between slopes of groups
#if significant means int model better than comb
pairwise_size_caudal_sex_dist <- summary(pairwise_size_caudal_sex, confidence = 0.95, test.type = "dist") 
pairwise_size_caudal_sex_dist

#Correlation between slope vectors (and angles) - similarity of vector orientation or angle,
#if significant means the vectors of the groups are oriented in different ways 
pairwise_size_caudal_sex_VC <- summary(pairwise_size_caudal_sex, confidence = 0.95, test.type = "VC",
                                        angle.type = "deg") 
pairwise_size_caudal_sex_VC

#Absolute difference between slope vector lengths - difference in rate of change per covariate unit (size),
#if significant means there is a significant rate of change difference in size between groups during growth
pairwise_size_caudal_sex_DL <-summary(pairwise_size_caudal_sex, confidence = 0.95, test.type = "DL") 
pairwise_size_caudal_sex_DL

#Save results to file
sink("Output/ANOVA_size_caudal_sex.txt")
print("ANOVA")
summary(size_caudal_sex)

print("1-Pairwise absolute distances slopes")
summary(pairwise_size_caudal_sex, confidence = 0.95, test.type = "dist") 

print("2-Distance between angles (slope directions)")
summary(pairwise_size_caudal_sex, confidence = 0.95, test.type = "VC",
        angle.type = "deg") 

print("3-Difference in slope vector length (difference in rate of change of size per unit of size)")
summary(pairwise_size_caudal_sex, confidence = 0.95, test.type = "DL") 
sink()

##Size vs Age ----
size_caudal_age <- lm.rrpp(gdf_caudal$size ~ gdf_caudal$age, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of size with logCS
summary(size_caudal_age)

#Pairwise comparison for the combination and interaction model
#Helps determine if there is a significant difference in slope (int model) in the size trajectory on top of difference in intercept (comb model)
pairwise_size_caudal_age <- pairwise(size_caudal_age, 
                                      groups = gdf_caudal$age, print.progress = FALSE) 
pairwise_size_caudal_age

#Distances between slope vectors (end-points) - absolute difference between slopes of groups
#if significant means int model better than comb
pairwise_size_caudal_age_dist <- summary(pairwise_size_caudal_age, confidence = 0.95, test.type = "dist") 
pairwise_size_caudal_age_dist

#Correlation between slope vectors (and angles) - similarity of vector orientation or angle,
#if significant means the vectors of the groups are oriented in different ways 
pairwise_size_caudal_age_VC <- summary(pairwise_size_caudal_age, confidence = 0.95, test.type = "VC",
                                        angle.type = "deg") 
pairwise_size_caudal_age_VC

#Absolute difference between slope vector lengths - difference in rate of change per covariate unit (size),
#if significant means there is a significant rate of change difference in size between groups during growth
pairwise_size_caudal_age_DL <-summary(pairwise_size_caudal_age, confidence = 0.95, test.type = "DL") 
pairwise_size_caudal_age_DL

#Save results to file
sink("Output/ANOVA_size_caudal_age.txt")
print("ANOVA")
summary(size_caudal_age)

print("1-Pairwise absolute distances slopes")
summary(pairwise_size_caudal_age, confidence = 0.95, test.type = "dist") 

print("2-Distance between angles (slope directions)")
summary(pairwise_size_caudal_age, confidence = 0.95, test.type = "VC",
        angle.type = "deg") 

print("3-Difference in slope vector length (difference in rate of change of size per unit of size)")
summary(pairwise_size_caudal_age, confidence = 0.95, test.type = "DL") 
sink()


#PHENOTYPIC TRAJECTORY SEXUAL DIMORPHISM ---- 

#Make sure to only use males and females!

#Shape changes along column at different ages
fit_shape_caudal_age <- procD.lm(gdf_caudal$coords[,,-c(rm_U_caudal)] ~ as.vector(gdf_caudal$sex[-rm_U_caudal]) * gdf_caudal$age[-rm_U_caudal], iter = 999, RRPP = F)


#Check that there is a significant correlation
summary(fit_shape_caudal_age)

#Use fit to calculate trajectories
trajectory_caudal_age <- trajectory.analysis(fit_shape_caudal_age, groups = as.vector(gdf_caudal$sex[-rm_U_caudal]), traj.pts = gdf_caudal$age[-rm_U_caudal], 
                                             pca = TRUE) 

#View results
#Magnitude differences between trajectories, standard summary - are trajectories different in length?
trajectory_caudal_age_MD <- summary(trajectory_caudal_age, show.trajectories = TRUE, attribute = "MD") 
trajectory_caudal_age_MD 
#Trajectory correlations -  are trajectories different in angle/direction?
trajectory_caudal_age_TC <- summary(trajectory_caudal_age, show.trajectories = TRUE, attribute = "TC", angle.type = "deg")
trajectory_caudal_age_TC
#Trajectory shape differences - are trajectories different in shape?
trajectory_caudal_age_SD <- summary(trajectory_caudal_age, show.trajectories = TRUE, attribute = "SD") 
trajectory_caudal_age_SD

#Save results to file
sink("Output/trajectory_caudal_age.txt")
print("Initial fit shape ~ sex * age")
summary(fit_shape_caudal_age)
print("Magnitude difference (absolute difference between path distances) - length")
summary(trajectory_caudal_age, show.trajectories = F, attribute = "MD") 
print("Correlations (angles) between trajectories - direction")
summary(trajectory_caudal_age, show.trajectories = F, attribute = "TC", angle.type = "deg")
print("Shape differences between trajectory vectors - shape")
summary(trajectory_caudal_age, show.trajectories = F, attribute = "SD") 
sink() 


#Plot results - PCA of fitted values
trajectory_caudal_age_plot <- plot(trajectory_caudal_age, main = "Trajectories shape change during growth by sex", pch = shapes[c(1,2)],  #title and type of point to be used
                                   col = c("darkgray","lightgray"), bg = c("darkgray","lightgray"), cex = 1, font.main = 2) #improve graphics
#Add line between groups
add.trajectories(trajectory_caudal_age_plot, traj.pch = shapes[c(1,2)],
                 traj.col = 1, traj.lty = 1, traj.lwd = 1, traj.cex = 1.5, traj.bg = 1, 
                 start.bg = "green", end.bg = "red") #trajectory line graphics
#Add legend to see which trajectory belongs to each group
legend("bottomleft", pch = shapes[c(1,2)],legend = c("Female"  ,  "Male"), pt.bg = 1, cex = 1)

##Make better PCA plot using ggplot
#Read PC scores as tibble
trajectory_caudal_age_pcscores <- as_tibble(trajectory_caudal_age_plot [["pc.points"]])

#Add group names and other attributes to tibble as columns
trajectory_caudal_age_pcscores <- trajectory_caudal_age_pcscores %>% mutate(sex = gdf_caudal$sex[-rm_U_caudal], age = gdf_caudal$age[-rm_U_caudal])
glimpse(trajectory_caudal_age_pcscores)

#Calculate means of PC1 and PC2 at each vertebra per group to draw trajectories
trajectory_caudal_age_pcscores_means <- trajectory_caudal_age_pcscores %>% group_by(age, sex) %>%
  summarise_at(vars(PC1, PC2), list(mean = mean))              #function for means, both columns
glimpse(trajectory_caudal_age_pcscores_means)

#Rename columns so they are easier to use for plot
trajectory_caudal_age_pcscores_means <- trajectory_caudal_age_pcscores_means %>% rename(x = PC1_mean, y = PC2_mean)
glimpse(trajectory_caudal_age_pcscores_means)

#Nice trajectory plot by sex and age
trajectory_caudal_age_ggplot <- ggplot(trajectory_caudal_age_pcscores, aes(x = PC1, y = PC2, shape = age, group = sex))+
  geom_point(size = 2.2, colour = "black", fill = "lightgray", alpha = 0.5, show.legend = F)+
  geom_point(data = trajectory_caudal_age_pcscores_means, aes(x = x, y = y, fill = sex, shape = age, group = sex), colour = "black",
            size = 5, inherit.aes = F, alpha = 0.8)+
  geom_path(data = trajectory_caudal_age_pcscores_means, aes(x = x, y = y, colour = sex, group = sex), inherit.aes = F, linewidth = 1,
            linejoin = 'mitre', arrow = arrow(angle = 40, length = unit(0.03, "npc"), ends = "first", type = "open"), show.legend = F)+
  scale_colour_manual(name = "Sex", labels = c("Female", "Male"), values = mypalette_sex, aesthetics = c("colour", "fill"))+
  scale_shape_manual(name = "Age", labels = c("Adult"  ,  "Juvenile", "Neonate"), #copy from as.factor(genera)
                     values = shapes)+
  theme_bw(base_size = 14)+
  xlab(paste0("PC 1 (",round(as.numeric(trajectory_caudal_age[["pca"]][["sdev"]][1]^2/sum(trajectory_caudal_age[["pca"]][["sdev"]]^2)*100), digits = 2),"%)"))+ 
  ylab(paste0("PC 2 (",round(as.numeric(trajectory_caudal_age[["pca"]][["sdev"]][2]^2/sum(trajectory_caudal_age[["pca"]][["sdev"]]^2)*100), digits = 2),"%)"))+
  theme(legend.key = element_blank(), legend.title = element_text(size = 15, face = "bold"), legend.text = element_text(size = 14),
        legend.background = element_blank(), legend.box.background =  element_blank())+
  guides(color = guide_legend(override.aes = list(shape = 23, colour = "black", fill = mypalette_sex)))+
  annotate("text", x = 0.09, y = -0.12, label = "Caudal",  fontface = 3, size = 12)

#Visualize plot and save as PDF using menu in bar on the right
trajectory_caudal_age_ggplot
