#CH.1 - Import landmarks, fix missing and mirror

#LOAD LIBRARIES ----
#always do this first!!
library(tidyverse)
library(Morpho)
library(geomorph)
library(Rvcg)
library(paleomorph)
library(EMMLi)
library(qgraph)
library(ape)
library(geiger)
library(abind)
library(devtools)
library(magick)
library(cli)
library(tibble)
library(rlang)
library(dplyr)
#install the two different surge package

#devtools::install_github("rnfelice/SURGE", force=TRUE) 

library(SURGE)

#DATA IMPORT caudal----

#Import LMs list - curves listed in curve_table
LM_caudal_table <-  read_csv("Data/landmark_caudal.csv")

#Import table defining curves
curve_caudal_table <- read_csv("Data/curves_caudal.csv")

#curve_table$bone[11:12] <- "premaxilla" to change bone group if needed for plots

#Identify the folder where your pts files are (I have a different one for each vertebral type)

ptsfolder_caudal <- "Data/pts/caudal_pts"

#Import the pts file names
ptslist_caudal <- dir(ptsfolder_caudal, pattern='.pts', recursive=F)

#Create curve info and indicate fixed points
my_curves_caudal <- create_curve_info(curve_caudal_table, n_fixed = 13)

##Check that your pts and ply match
#Identify ply folder
plyfolder_caudal <- "Data/ply/caudal_ply"

#Make list of ply file names
plylist_caudal <-  dir(plyfolder_caudal, pattern='\\.ply$', recursive=F)

#Make sublist to check names
ptslist11<-gsub(pattern="\\.pts$","",ptslist_caudal)
plylist11<-gsub(pattern="\\.ply$","",plylist_caudal)

#Check that names match up
identical(plylist11,ptslist11) #should be TRUE if order the same

#Set wd where the checkpoint data is
setwd(ptsfolder_caudal)

#Write file with list of specimens/vertebrae
write.csv(plylist9, file = "caudal.csv")

###RESAMPLE AND FIX MISSING DATA ----

#Import pts data
subsampled.lm_caudal <- import_chkpt_data(ptslist_caudal, my_curves_caudal, subsampl = TRUE, verbose=TRUE)

#If you have any missing points, Checkpoint will set them to x,y,z=9999
#This makes for bad plots in checkLM, so switch them to NA

subsampled.lm_caudal[subsampled.lm_caudal == 9999] <- NA


###SET WD to ply from console!! -->

#Make sure you have converted your plys from binary to ASCII
#Check to make sure your curves look okay on each specimen

specs_tofix_caudal <- checkLM(subsampled.lm_caudal, path="", pt.size = 4, suffix=".ply", render = "s", begin = 158, point = "s")
specs_tofix2caudal <- checkLM(subsampled.lm_caudal, path="", pt.size = 4, suffix=".ply", render = "s", begin = 177, point = "s")

#Write down specs numbers as you will need to clean evironment to import them again

#Check single specimen problems with spheres
checkLM(subsampled.lm_caudal, path="", pt.size = 1, suffix=".ply", render = "s", begin = 177, point = "s")

#Create object with new resampled points
newpts_caudal <- subsampled.lm_caudal

#Create missing list 
misslist_caudal <- createMissingList(dim(newpts_caudal)[3])

for (j in 1:dim(newpts_caudal)[[3]]){
  misslist_caudal[[j]]<-which(is.na(newpts_caudal[,1,j]))
} 


##Fix missing landmarks
newpts2_caudal <- fixLMtps(newpts_caudal)

checkLM(newpts2_caudal$out, path="", pt.size = 1, suffix=".ply", render="s", begin = 177)


#SLIDE ----

#================================#
#      RUN Slider3D CODE         #
#================================#

###!!SET WD to ply from console!! -->

#Slide points on surface
slided4.all_caudal <- slider3d_2(newpts2_caudal$out, SMvector= my_curves_caudal$Sliding.LMs,
                                   outlines = my_curves_caudal$Curves, 
                                   #copy ply folder path to ensure it works - set as working directory from console to print
                                   sur.path = ".",
                                   sur.name = NULL, 
                                   meshlist = paste("",dimnames(newpts2_caudal$out)[[3]],".ply",sep=""), ignore = NULL,
                                   sur.type = "ply", tol = 1e-10, deselect = FALSE, inc.check = FALSE,
                                   recursive = TRUE, iterations = 3, initproc = TRUE,
                                   pairedLM = 0, mc.cores = 1, bending=TRUE,
                                   fixRepro = F, stepsize=0.2,
                                   missingList = misslist_caudal)

#Name specimens slided data
dimnames(slided4.all_caudal[["dataslide"]])[3] <- dimnames(newpts2_caudal$out)[3]

#Re-estimate missing post sliding
slided4.all_caudal$dataslide[which(is.na(newpts_caudal))] <- NA

#Fix missing landmarks
slid.lms_caudal <- fixLMtps(slided4.all_caudal$dataslide)

#Extract landmark array from object
slidedlms_caudal <- slid.lms_caudal$out

plotOutliers(slidedlms_caudal)

#Check to see how they look (plys must be ASCII format)
specs_tofix_slid_caudal <- checkLM(slidedlms_caudal, path="", pt.size = 5, 
                                     suffix=".ply", render = "s", begin = 104, point = "s")

#Save slided LMs as R data file
save(slidedlms_caudal, file = 'slidedlms_caudal.Rdata')

#List of points and curves for different bones - useful for plots
centrum_caudal <- c(LM_caudal_table$lm[which(LM_caudal_table$bone%in%c("centra"))], 
                      my_curves_caudal$Curves[which(curve_caudal_table$bone%in%c("centra"))]) %>% unlist(.)%>%unique(.)%>%sort(.)
process_caudal <- c(LM_caudal_table$lm[which(LM_caudal_table$bone%in%c("process"))], 
                      my_curves_caudal$Curves[which(curve_caudal_table$bone%in%c("process"))]) %>% unlist(.)%>%unique(.)%>%sort(.)

#MIRROR SYMMETRICAL TRAITS ---- 

#LMs and curves that are on the midline 
midline_caudal <- as.integer(c(LM_caudal_table$lm[which(LM_caudal_table$side%in%c("M"))], 
                                 my_curves_caudal$Curves[which(curve_caudal_table$side%in%c("M"))]) %>% unlist(.)%>%unique(.)%>%sort(.))


#Generate length with the both sides - midline
mirrormatrix_caudal <- dim(slidedlms_caudal)[1]*2-length(midline_caudal)

#create new matrix with size of full dataset
m.matrix_caudal <- array(NA,dim = c(mirrormatrix_caudal,3,dim(slidedlms_caudal)[3])) #n landmarks, 3 dimesions, n of specimens
m.matrix_caudal[1:dim(slidedlms_caudal)[1],,] <- slidedlms_caudal

#Left side
left.lm_caudal <- c(1:dim(slidedlms_caudal)[1])[-midline_caudal] # exclude midline points

#Right side
right.lm_caudal <- c(c(dim(slidedlms_caudal)[1]+1):mirrormatrix_caudal) # space for the LM on the left side = last point on the right side + # of point needed on each side

#Landmark pairs
bilat.landmarks_caudal <- cbind(left.lm_caudal,right.lm_caudal) #one column is the rows of the right side LM and the other column the rows of the left side

#Mirror landmarks
mirrored_caudal <- mirrorfill(A = m.matrix_caudal, l1 = midline_caudal, l2 = bilat.landmarks_caudal) # the NA rows are now filled
mirrored_caudal

#Plot mirrored points
spheres3d(mirrored_caudal[,,5],col=2,radius=4)

spheres3d(mirrored_caudal[left.lm_caudal,,5],col=40,radius=4)
spheres3d(mirrored_caudal[midline_caudal,,5],col=205,radius=4)
spheres3d(mirrored_caudal[right.lm_caudal,,5],col=254,radius=4)

#Set specimens names
dimnames(mirrored_caudal)[3] <- dimnames(slidedlms_caudal)[3]

#Create new object for analyses with all mirrored data, include only shape data
final_dataset_caudal <- mirrored_caudal




