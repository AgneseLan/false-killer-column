#CH.1 - Import landmarks, fix missing and mirror

#LOAD LIBRARIES ----
#always do this first!!library(tidyverse)
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

#DATA IMPORT THORACIC----

#Import LMs list - curves listed in curve_table
LM_thoracic_table <-  read_csv("Data/landmark_thoracic.csv")

#Import table defining curves
curve_thoracic_table <- read_csv("Data/curves_thoracic.csv")

#curve_table$bone[11:12] <- "premaxilla" to change bone group if needed for plots

#Identify the folder where your pts files are (I have a different one for each vertebral type)
ptsfolder_thoracic <- "Data/pts/thoracic_pts"

#Import the pts file names
ptslist_thoracic <- dir(ptsfolder_thoracic, pattern='.pts', recursive=F)

#Create curve info and indicate fixed points
my_curves_thoracic <- create_curve_info(curve_thoracic_table, n_fixed = 16)

##Check that your pts and ply match
#Identify ply folder
plyfolder_thoracic <- "Data/ply/thoracic_ply"

#Make list of ply file names
plylist_thoracic <-  dir(plyfolder_thoracic, pattern='\\.ply$', recursive=F)

#Make sublist to check names
ptslist2<-gsub(pattern="\\.pts$","",ptslist_thoracic)
plylist2<-gsub(pattern="\\.ply$","",plylist_thoracic)

#Check that names match up
identical(plylist2,ptslist2) #should be TRUE if order the same

#Write file with list of specimens/vertebrae
write.csv(plylist2, file = "Output/thoracic.csv")

#Set wd where the checkpoint data is
setwd(ptsfolder_thoracic)

###RESAMPLE AND FIX MISSING DATA ----

#Import pts data
subsampled.lm_thoracic <- import_chkpt_data(ptslist_thoracic, my_curves_thoracic, subsampl = TRUE, verbose=TRUE)

#If you have any missing points, Checkpoint will set them to x,y,z=9999
#This makes for bad plots in checkLM, so switch them to NA

subsampled.lm_thoracic[subsampled.lm_thoracic == 9999] <- NA

###SET WD to thoracic_ply from console!! -->

#Make sure you have converted your plys from binary to ASCII
#Check to make sure your curves look okay on each specimen

specs_tofix_thoracic <- checkLM(subsampled.lm_thoracic, path="", pt.size = 4, suffix=".ply", render = "s", begin = 1, point = "s")
specs_tofix2thoracic <- checkLM(subsampled.lm_thoracic, path="", pt.size = 4, suffix=".ply", render = "s", begin = 104, point = "s")

#Write down specs numbers as you will need to clean evironment to import them again

#Check single specimen problems with spheres
checkLM(subsampled.lm_thoracic, path="", pt.size = 1, suffix=".ply", render = "s", begin = 1, point = "s")

#Create object with new resampled points
newpts_thoracic <- subsampled.lm_thoracic

#Create missing list 
misslist_thoracic <- createMissingList(dim(newpts_thoracic)[3])

for (j in 1:dim(newpts_thoracic)[[3]]){
  misslist_thoracic[[j]]<-which(is.na(newpts_thoracic[,1,j]))
} 

##Fix missing landmarks
newpts2_thoracic <- fixLMtps(newpts_thoracic)

checkLM(newpts2_thoracic$out, path="", pt.size = 1, suffix=".ply", render="s", begin = 1)

#SLIDE ----

#================================#
#      RUN Slider3D CODE         #
#================================#

###!!SET WD to ply from console!! -->

#Slide points on surface
slided4.all_thoracic <- slider3d_2(newpts2_thoracic$out, SMvector= my_curves_thoracic$Sliding.LMs,
                          outlines = my_curves_thoracic$Curves, 
                          #copy ply folder path to ensure it works - set as working directory from console to print
                          sur.path = ".",
                          sur.name = NULL, 
                          meshlist = paste("",dimnames(newpts2_thoracic$out)[[3]],".ply",sep=""), ignore = NULL,
                          sur.type = "ply", tol = 1e-10, deselect = FALSE, inc.check = FALSE,
                          recursive = TRUE, iterations = 3, initproc = TRUE,
                          pairedLM = 0, mc.cores = 1, bending=TRUE,
                          fixRepro = F, stepsize=0.2,
                          missingList = misslist_thoracic)

#Name specimens slided data
dimnames(slided4.all_thoracic[["dataslide"]])[3] <- dimnames(newpts2_thoracic$out)[3]

#Re-estimate missing post sliding
slided4.all_thoracic$dataslide[which(is.na(newpts_thoracic))] <- NA

#Fix missing landmarks
slid.lms_thoracic <- fixLMtps(slided4.all_thoracic$dataslide)

#Extract landmark array from object
slidedlms_thoracic <- slid.lms_thoracic$out

#Check to see how they look (plys must be ASCII format)
specs_tofix_slid_thoracic <- checkLM(slidedlms_thoracic, path="", pt.size = 5, 
                            suffix=".ply", render = "s", begin = 1, point = "s")

##CHANGE WD to source folder in console --->

#Save slided LMs as R data file
save(slidedlms_thoracic, file = 'Output/slidedlms_thoracic.Rdata')

#List of points and curves for different bones - useful for plots
centrum_thoracic <- c(LM_thoracic_table$lm[which(LM_thoracic_table$bone%in%c("centra"))], 
                      my_curves_thoracic$Curves[which(curve_thoracic_table$bone%in%c("centra"))]) %>% unlist(.)%>%unique(.)%>%sort(.)
process_thoracic <- c(LM_thoracic_table$lm[which(LM_thoracic_table$bone%in%c("process"))], 
                      my_curves_thoracic$Curves[which(curve_thoracic_table$bone%in%c("process"))]) %>% unlist(.)%>%unique(.)%>%sort(.)

#MIRROR SYMMETRICAL TRAITS ---- 

#LMs and curves that are on the midline 
midline_thoracic <- as.integer(c(LM_thoracic_table$lm[which(LM_thoracic_table$side%in%c("M"))], 
                                 my_curves_thoracic$Curves[which(curve_thoracic_table$side%in%c("M"))]) %>% unlist(.)%>%unique(.)%>%sort(.))


#Generate length with the both sides - midline
mirrormatrix_thoracic <- dim(slidedlms_thoracic)[1]*2-length(midline_thoracic)

#create new matrix with size of full dataset
m.matrix_thoracic <- array(NA,dim = c(mirrormatrix_thoracic,3,dim(slidedlms_thoracic)[3])) #n landmarks, 3 dimesions, n of specimens
m.matrix_thoracic[1:dim(slidedlms_thoracic)[1],,] <- slidedlms_thoracic

#Left side
left.lm_thoracic <- c(1:dim(slidedlms_thoracic)[1])[-midline_thoracic] # exclude midline points

#Right side
right.lm_thoracic <- c(c(dim(slidedlms_thoracic)[1]+1):mirrormatrix_thoracic) # space for the LM on the left side = last point on the right side + # of point needed on each side

#Landmark pairs
bilat.landmarks_thoracic <- cbind(left.lm_thoracic,right.lm_thoracic) #one column is the rows of the right side LM and the other column the rows of the left side

#Mirror landmarks
mirrored_thoracic <- mirrorfill(A = m.matrix_thoracic, l1 = midline_thoracic, l2 = bilat.landmarks_thoracic) # the NA rows are now filled
mirrored_thoracic

#Plot mirrored points
spheres3d(mirrored_thoracic[,,5],col=2,radius=4)

spheres3d(mirrored_thoracic[left.lm_thoracic,,5],col=40,radius=4)
spheres3d(mirrored_thoracic[midline_thoracic,,5],col=205,radius=4)
spheres3d(mirrored_thoracic[right.lm_thoracic,,5],col=254,radius=4)

#Set specimens names
dimnames(mirrored_thoracic)[3] <- dimnames(slidedlms_thoracic)[3]

#Create new object for analyses with all mirrored data, include only shape data
final_dataset_thoracic <- mirrored_thoracic




