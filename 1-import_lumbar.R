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

#DATA IMPORT lumbar----

#Import LMs list - curves listed in curve_table
LM_lumbar_table <-  read_csv("Data/landmark_lumbar.csv")

#Import table defining curves
curve_lumbar_table <- read_csv("Data/curves_lumbar.csv")

#Identify the folder where your pts files are (I have a different one for each vertebral type)
ptsfolder_lumbar <- "Data/pts/lumbar_pts"

#Import the pts file names
ptslist_lumbar <- dir(ptsfolder_lumbar, pattern='.pts', recursive=F)

#Create curve info and indicate fixed points
my_curves_lumbar <- create_curve_info(curve_lumbar_table, n_fixed = 16)

##Check that your pts and ply match
#Identify ply folder
plyfolder_lumbar <- "Data/ply/lumbar_ply"

#Make list of ply file names
plylist_lumbar <-  dir(plyfolder_lumbar, pattern='\\.ply$', recursive=F)

#Make sublist to check names
ptslist8<-gsub(pattern="\\.pts$","",ptslist_lumbar)
plylist8<-gsub(pattern="\\.ply$","",plylist_lumbar)

#Check that names match up
identical(plylist8,ptslist8) #should be TRUE if order the same

#Write file with list of specimens/vertebrae
write.csv(plylist8, file = "Output/lumbar.csv")

#Set wd where the checkpoint data is
setwd(ptsfolder_lumbar)

###RESAMPLE AND FIX MISSING DATA ----

#Import pts data
subsampled.lm_lumbar <- import_chkpt_data(ptslist_lumbar, my_curves_lumbar, subsampl = TRUE, verbose=TRUE)

#If you have any missing points, Checkpoint will set them to x,y,z=9999
#This makes for bad plots in checkLM, so switch them to NA

subsampled.lm_lumbar[subsampled.lm_lumbar == 9999] <- NA


###SET WD to lumbar_ply from console!! -->

#Make sure you have converted your plys from binary to ASCII
#Check to make sure your curves look okay on each specimen

specs_tofix_lumbar <- checkLM(subsampled.lm_lumbar, path="", pt.size = 4, suffix=".ply", render = "s", begin = 152, point = "s")
specs_tofix2lumbar <- checkLM(subsampled.lm_lumbar, path="", pt.size = 4, suffix=".ply", render = "s", begin = 152, point = "s")

#Write down specs numbers as you will need to clean evironment to import them again

#Check single specimen problems with spheres
checkLM(subsampled.lm_lumbar, path="", pt.size = 1, suffix=".ply", render = "s", begin = 1, point = "s")

#Create object with new resampled points
newpts_lumbar <- subsampled.lm_lumbar

#Create missing list 
misslist_lumbar <- createMissingList(dim(newpts_lumbar)[3])

for (j in 1:dim(newpts_lumbar)[[3]]){
  misslist_lumbar[[j]]<-which(is.na(newpts_lumbar[,1,j]))
} 

##Fix missing landmarks
newpts2_lumbar <- fixLMtps(newpts_lumbar)

checkLM(newpts2_lumbar$out, path="", pt.size = 1, suffix=".ply", render="s", begin = 1)


#SLIDE ----

#================================#
#      RUN Slider3D CODE         #
#================================#

###!!SET WD to ply from console!! -->

#Slide points on surface
slided4.all_lumbar <- slider3d_2(newpts2_lumbar$out, SMvector= my_curves_lumbar$Sliding.LMs,
                                   outlines = my_curves_lumbar$Curves, 
                                   #copy ply folder path to ensure it works - set as working directory from console to print
                                   sur.path = ".",
                                   sur.name = NULL, 
                                   meshlist = paste("",dimnames(newpts2_lumbar$out)[[3]],".ply",sep=""), ignore = NULL,
                                   sur.type = "ply", tol = 1e-10, deselect = FALSE, inc.check = FALSE,
                                   recursive = TRUE, iterations = 3, initproc = TRUE,
                                   pairedLM = 0, mc.cores = 1, bending=TRUE,
                                   fixRepro = F, stepsize=0.2,
                                   missingList = misslist_lumbar)

#Name specimens slided data
dimnames(slided4.all_lumbar[["dataslide"]])[3] <- dimnames(newpts2_lumbar$out)[3]

#Re-estimate missing post sliding
slided4.all_lumbar$dataslide[which(is.na(newpts_lumbar))] <- NA

#Fix missing landmarks
slid.lms_lumbar <- fixLMtps(slided4.all_lumbar$dataslide)

#Extract landmark array from object
slidedlms_lumbar <- slid.lms_lumbar$out

plotOutliers(slidedlms_lumbar)

#Check to see how they look (plys must be ASCII format)
specs_tofix_slid_lumbar <- checkLM(slidedlms_lumbar, path="", pt.size = 5, 
                                     suffix=".ply", render = "s", begin = 1, point = "s")

##CHANGE WD to source folder in console --->

#Save slided LMs as R data file
save(slidedlms_lumbar, file = 'Output/slidedlms_lumbar.Rdata')

#List of points and curves for different bones - useful for plots
centrum_lumbar <- c(LM_lumbar_table$lm[which(LM_lumbar_table$bone%in%c("centra"))], 
                      my_curves_lumbar$Curves[which(curve_lumbar_table$bone%in%c("centra"))]) %>% unlist(.)%>%unique(.)%>%sort(.)
process_lumbar <- c(LM_lumbar_table$lm[which(LM_lumbar_table$bone%in%c("process"))], 
                      my_curves_lumbar$Curves[which(curve_lumbar_table$bone%in%c("process"))]) %>% unlist(.)%>%unique(.)%>%sort(.)

#MIRROR SYMMETRICAL TRAITS ---- 

#LMs and curves that are on the midline 
midline_lumbar <- as.integer(c(LM_lumbar_table$lm[which(LM_lumbar_table$side%in%c("M"))], 
                                 my_curves_lumbar$Curves[which(curve_lumbar_table$side%in%c("M"))]) %>% unlist(.)%>%unique(.)%>%sort(.))


#Generate length with the both sides - midline
mirrormatrix_lumbar <- dim(slidedlms_lumbar)[1]*2-length(midline_lumbar)

#create new matrix with size of full dataset
m.matrix_lumbar <- array(NA,dim = c(mirrormatrix_lumbar,3,dim(slidedlms_lumbar)[3])) #n landmarks, 3 dimesions, n of specimens
m.matrix_lumbar[1:dim(slidedlms_lumbar)[1],,] <- slidedlms_lumbar

#Left side
left.lm_lumbar <- c(1:dim(slidedlms_lumbar)[1])[-midline_lumbar] # exclude midline points

#Right side
right.lm_lumbar <- c(c(dim(slidedlms_lumbar)[1]+1):mirrormatrix_lumbar) # space for the LM on the left side = last point on the right side + # of point needed on each side

#Landmark pairs
bilat.landmarks_lumbar <- cbind(left.lm_lumbar,right.lm_lumbar) #one column is the rows of the right side LM and the other column the rows of the left side

#Mirror landmarks
mirrored_lumbar <- mirrorfill(A = m.matrix_lumbar, l1 = midline_lumbar, l2 = bilat.landmarks_lumbar) # the NA rows are now filled
mirrored_lumbar

#Plot mirrored points
spheres3d(mirrored_lumbar[,,5],col=2,radius=4)

spheres3d(mirrored_lumbar[left.lm_lumbar,,5],col=40,radius=4)
spheres3d(mirrored_lumbar[midline_lumbar,,5],col=205,radius=4)
spheres3d(mirrored_lumbar[right.lm_lumbar,,5],col=254,radius=4)

#Set specimens names
dimnames(mirrored_lumbar)[3] <- dimnames(slidedlms_lumbar)[3]

#Create new object for analyses with all mirrored data, include only shape data
final_dataset_lumbar <- mirrored_lumbar



