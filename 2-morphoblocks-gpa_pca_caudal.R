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
checkLM(final_dataset_caudal, path="Data/ply/caudal_ply/", pt.size = 2, suffix=".ply", render = "s", begin = 177, point = "s")


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

#GPA ALIGNMENT ----

#Procrustes alignment, should also show mean config coordinates
gpa_caudal <- gpagen(final_dataset_caudal) 

#Log-transform Centroid size as object
logCsize_caudal <- log10(gpa_caudal$Csize) 

#Save mean shape to create links
mean_shape_caudal <- gpa_caudal$consensus 

#Coordinates of all specimens after gpa_lumbar alignment
coords_caudal<- gpa_caudal$coords[c(left.lm_lumbar, midline_lumbar),,]

#Plot all specimens with mean to check that all landmarks are ok
plotAllSpecimens(coords_caudal, mean = TRUE, label = F, plot.param = list(pt.cex = 0.05, mean.cex = 3, mean.bg = "black"))

#Check for outliers, they would be displayed in red - most immature ones are normal as outliers
plotOutliers(coords_caudal)

#Plot landmarks from outliers in 3D to check how they look
spheres3d(coords_caudal[,,3], r = 0.002)

#checkLM(final_dataset_caudal, path="", pt.size = 2, suffix=".ply", render="s", begin = 65) 
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
warp_specimen_caudal <- coords_caudal[,,43] #number displayed by findMeanSpec
warp_specimen_caudal

#Import simplified mesh to create warp mesh on
mesh_3D_caudal <- read.ply("Data/ply/refmesh_caudal.ply") #make sure NO binary encoding (ASCII)

#Check range of mesh and coordinates to make sure it has same scale
range(mesh_3D_caudal$vb[1:3,]) #if this is too big/small, scale in editor and re-import
range(warp_specimen_caudal)

#PCA COMPLETE DATASET ----

#Run PCA on complete dataset
PCA_all_caudal <- gm.prcomp(gdf_caudal$coords)

#List of PC components and proportion of variation
PCA_all_caudal

#Save PCA results to file
sink("output/PCA_all_caudal_components.txt")
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

#Show 3D deformation from mean with points overlay, do this for all 4 extremes - using spheres3D for points
#PC1min colors
#spheres3d(mean_shape, radius=.0005, color = "gray60", alpha = 0.5, fastTransparency = T) - plot mean specimens with transparency
PC1min_all_caudal_points <- c(
  spheres3d(PC1min_all_caudal[centrum_caudal,], radius=.001, color = "red"),
  spheres3d(PC1min_all_caudal[process_caudal,], radius=.001, color = "blue"))

PC1max_all_caudal_points <- c(
  spheres3d(PC1min_all_caudal[centrum_caudal,], radius=.001, color = "red"),
  spheres3d(PC1min_all_caudal[process_caudal,], radius=.001, color = "blue"))

PC2min_all_caudal_points <- c(
  spheres3d(PC1min_all_caudal[centrum_caudal,], radius=.001, color = "red"),
  spheres3d(PC1min_all_caudal[process_caudal,], radius=.001, color = "blue"))

PC2max_all_caudal_points <- c(
  spheres3d(PC1min_all_caudal[centrum_caudal,], radius=.001, color = "red"),
  spheres3d(PC1min_all_caudal[process_caudal,], radius=.001, color = "blue"))

#3D deformation using PLY files
#Extremes PC axes

##Create warp mesh, to use as reference for visualization of analyses
PC1min_all_caudal_mesh <- warpRefMesh(mesh = mesh_3D_caudal, mesh.coord = warp_specimen_caudal, 
                                        ref = PC1min_all_caudal, color = "gray", centered = FALSE)
##Save warped mesh as ply (Morpho)
mesh2ply(PC1min_all_caudal_mesh, filename = "output/PC1min_all_caudal")

##Create warp mesh, to use as reference for visualization of analyses
PC1max_all_caudal_mesh <- warpRefMesh(mesh = mesh_3D_caudal, mesh.coord = warp_specimen_caudal, 
                                        ref = PC1max_all_caudal, color = "gray", centered = FALSE)
##Save warped mesh as ply (Morpho)
mesh2ply(PC1max_all_caudal_mesh, filename = "output/PC1max_all_caudal")

##Create warp mesh, to use as reference for visualization of analyses
PC2min_all_caudal_mesh <- warpRefMesh(mesh = mesh_3D_caudal, mesh.coord = warp_specimen_caudal, 
                                        ref = PC2min_all_caudal, color = "gray", centered = FALSE)
##Save warped mesh as ply (Morpho)
mesh2ply(PC2min_all_caudal_mesh, filename = "output/PC2min_all_caudal")

##Create warp mesh, to use as reference for visualization of analyses
PC2max_all_caudal_mesh <- warpRefMesh(mesh = mesh_3D_caudal, mesh.coord = warp_specimen_caudal, 
                                        ref = PC2max_all_caudal, color = "gray", centered = FALSE)
##Save warped mesh as ply (Morpho)
mesh2ply(PC2max_all_caudal_mesh, filename = "output/PC2max_all_caudal")

##Make better PCA plot using ggplot
#Read PC scores as tibble
pcscores_all_caudal <- as_tibble(pcscores_all_caudal)

#Add labels and other attributes to tibble as columns
pcscores_all_caudal <- pcscores_all_caudal %>% mutate(Id = gdf_caudal$number, age = gdf_caudal$age,
                                                          sex = gdf_caudal$sex, specimens = gdf_caudal$specimens, size = gdf_caudal$size,
                                                          TL = as.numeric(classifiers_caudal$size_cm))
glimpse(pcscores_all_caudal)

#Nice PCA plot with stages and groups
PCA_all_caudal_ggplot <- ggplot(pcscores_all_caudal, aes(x = Comp1, y = Comp2, label = Id, colour = age, fill = age))+
  geom_point(size = 3, aes(shape = sex))+
  geom_text_repel(colour = "black", size = 4, max.overlaps = 40)+
  scale_colour_manual(name = "Age", labels =  c("adult"  ,  "juvenile", "newborn" ), #to be ordered as they appear in tibble
                      values = mypalette_age_sex , aesthetics = c("colour","fill"))+            #legend and color adjustments
  scale_shape_manual(name = "Sex", labels = c( "F"   ,    "M"      , "unknown"), #copy from as.factor(sex)
                     values = shapes)+
  theme_bw()+
  xlab("PC 1 (55.66%)")+ #copy this from standard PCA plot (PCA_all_plot)
  ylab("PC 2 (20.74%)")+
  ggtitle("PCA all data caudal")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_caudal_ggplot

#Make hulls for PCA plot with hulls around age
hulls_all_age_caudal <- pcscores_all_caudal %>%
  group_by(age) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

#Nice PCA plot with hulls around age
PCA_all_caudal_age_ggplot <- ggplot(pcscores_all_caudal, aes(x = Comp1, y = Comp2, colour = age))+
  geom_point(size = 3, aes(shape = sex), fill = "white")+
  scale_colour_manual(name = "Age", labels =  c("adult" ,   "juvenile" ,"neonate" ), #to be ordered as they appear in tibble
                      values = c("#F0E442", "#6699CC", "#AA4499"))+            #legend and color adjustments
  geom_polygon(data = hulls_all_age_caudal, aes(x = x, y = y, colour = age, fill = age), 
               alpha = .3, show.legend = FALSE, inherit.aes = F)+ #colored hulls with transparency
  scale_fill_manual(name = "Age", labels = c("adult"  ,  "juvenile", "newborn" ),
                    values =  c("#F0E442", "#6699CC", "#AA4499"))+ #must match scale_colour_manual
  scale_shape_manual(name = "Sex", labels = c( "Female"   ,    "Male"      , "Neonate"), #copy from as.factor(sex)
                     values = shapes)+
  theme_bw()+
  xlab("PC 1 (55.66%)")+ #copy this from standard PCA plot (PCA_all_plot)
  ylab("PC 2 (20.66%)")

  #Remove legend for a scale_ using guide
  guides(fill = guide_legend(label = F, title = NULL, override.aes = list(shape = NA)))

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_caudal_age_ggplot

#Make hulls for PCA plot with hulls around sex
hulls_all_sex_caudal  <- pcscores_all_caudal  %>%
  group_by(sex) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

#Nice PCA plot with hulls around sex
PCA_all_caudal_sex_ggplot <- ggplot(pcscores_all_caudal[!pcscores_all_caudal$sex%in% c("unknown"),], aes(x = Comp1, y = Comp2, colour = sex))+
  geom_point(size = 3, aes(shape = age), fill = "white")+
  scale_colour_manual(name = "Sex", labels =  c("Female"   ,    "Male"      , "Unknown (neonate)"), #to be ordered as they appear in tibble
                      values=c("#E597B9", "#332288"))+            #legend and color adjustments
  geom_polygon(data = hulls_all_sex_caudal[!hulls_all_sex_caudal$sex%in% c("unknown"),], aes(x = x, y = y, colour = sex, fill = sex), 
               alpha = .3, show.legend = FALSE, inherit.aes = F)+ #colored hulls with transparency
  scale_fill_manual(name = "Sex", labels = c( "Female"   ,    "Male"      , "Unknown (neonate)"),
                    values =  c("#E597B9", "#332288"))+ #must match scale_colour_manual
  scale_shape_manual(name = "Age", labels = c( "adult"  ,  "juvenile", "neonate"), #copy from as.factor(sex)
                     values = shapes)+
  theme_bw()+
  xlab("PC 1 (55.66%)")+ #copy this from standard PCA plot (PCA_all_plot)
  ylab("PC 2 (20.66%)")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  #Remove legend for a scale_ using guide
  guides(fill = guide_legend(label = F, title = NULL, override.aes = list(shape = NA)))

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
sink("output/PC1-2all_caudal_size_lm.txt")
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
sink("output/PC1-2all_caudal_sex_lm.txt")
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
sink("output/PC1-2all_caudal_age_lm.txt")
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
sink("output/PC1-2all_caudal_TL_lm.txt")
print("PC1")
summary(reg_PC1all_caudal_TL)
anova(reg_PC1all_caudal_TL)
print("PC2")
summary(reg_PC2all_caudal_TL)
anova(reg_PC2all_caudal_TL)
sink() 

#Save results of all regressions to 1 file
sink("output/PC1-2_all_caudal_lm.txt")
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

#Calculate regression for each component taking TL into account
reg_PC1all_caudal_specs <- lm(Comp1 ~ specimens, data = pcscores_all_caudal)
reg_PC2all_caudal_specs <- lm(Comp2 ~ specimens, data = pcscores_all_caudal)

#View results and p-value
anova(reg_PC1all_caudal_specs)
anova(reg_PC2all_caudal_specs)

#GRAPHICS CODE ----
##Create project palettes----
mypalette_paired <- brewer.pal(12,"Paired")
image(1:12, 1, as.matrix(1:12), col = mypalette_paired, xlab = "Paired",
      ylab = "", yaxt = "n")

mypalette_blue <- as.matrix(ggthemes_data[["tableau"]][["color-palettes"]][["ordered-sequential"]][["Blue"]][["value"]])
image(1:20, 1, as.matrix(1:20), col = mypalette_blue, xlab = "Blue",
      ylab = "", yaxt = "n")

#Palette for age - early, late/new, immature, adult
mypalette_age_sex <- c(mypalette_paired[3], mypalette_paired[8], mypalette_paired[9])
image(1:3, 1, as.matrix(1:3), col = mypalette_age_sex,
      ylab = "", yaxt = "n")

#Palette to use for sex PCAs
mypalette_age_sex <- c(mypalette_paired[5], mypalette_paired[2], mypalette_paired[10])
image(1:3, 1, as.matrix(1:3), col = mypalette_age_sex,
      ylab = "", yaxt = "n")

#create shape palette for sex
shapes <- c(19, 15, 17)


#Create shape palette 3 age
shapes <- c(21, 22, 24)


###Plot landmarks on vertebra ----

##Save mesh with plotted landmarks
#Find mean specimen raw data
findMeanSpec(final_dataset_caudal)
#Get spec number
match("1961_6_14_22_3_e", dimnames(final_dataset_caudal)[[3]])

#Import simplified ply
#Less faces, no holes or isolated triangles
refmesh_caudal <- vcgImport("Data/ply/caudal_ply/1961_6_14_22_3_e.ply")

#Plot on surface
fixed_LMs <- c(1:13)

shade3d(refmesh_caudal, col = "white", alpha = 0.5, fastTransparency = T)
spheres3d(slidedlms_caudal[fixed_LMs,,30], col =  "turquoise4", type = "s",
          radius = 3, aspect = T, main = "mean",axes = F, main = F, fov = 0)
spheres3d(slidedlms_caudal[-fixed_LMs,,30], col =  "turquoise1", type = "s",
          radius = 2, aspect = T, main = "mean",axes = F, main = F, fov = 0)

rgl.snapshot(filename = "landmarks_caudal_dorsal.png")
rgl.snapshot(filename = "landmarks_caudal_lateral.png")
rgl.snapshot(filename = "landmarks_caudal_posterior.png")
rgl.snapshot(filename = "landmarks_caudal_anterior.png") 
