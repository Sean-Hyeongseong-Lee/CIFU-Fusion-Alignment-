# Module A: Curve fitting with p basis functions 
# i) obtain the basis coefficients for each OCT curve
# ii) Compute the FVE
#
###############################################################################
###############################################################################
# Prelimiaries
######################################
#
# Naming Convention: CurveType_MachineType_DiseaseStatus_Original/Normalized_AdditionalInfo
#   CurveType           : NRR (there is no RNFL curves in this study)
#   MachineType         : OCT or Fundus
#   DiseaseStatus       :  Normal or Glau (Glau = Glaucoma) 
#   Original/Normalized : Original = OG, Normalized = Norm
#   AdditionalInfo      : what specifically we are dealing with when it's not just the curves; 
#                         e.g. TV = total variation, FVE = Fraction of Variance Explained, etc
#   Eg: NRR_OCT_Normal_OG_p05 = NRR curve, OCT recorded data, Normal/Non-Glaucoma, 5 basis functions
#
#####################################
#####     Set the libraries     #####
#####################################
# These were the libraries used for the CIFU project
library.set <- c("readr","openxlsx","readxl","plyr","dplyr","reshape2","fda","fdasrvf","funFEM",
                 "car","knitr","kableExtra","xtable","pander","png","grid","gridExtra","pracma",
                 "fda.usc", "ggplot2", "ggpubr", "scales","grid","gridExtra","StatMatch",
                 "slam","movMF", "gplots", "vegan", "RColorBrewer", "writexl") 
# The "Heatplus" package isn't avail for this version of R

# This should install any missing packages
new.packages <- library.set[!(library.set %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0) install.packages(new.packages, dependencies = T)

lapply(library.set, library, character.only = TRUE)

###################################################################
# Read in the data                                                #
# You will need to change the file paths to wherever the files    #
# are located on your machine                                     #
###################################################################
localPath_normal_NRR <- "~/Desktop/CIFU_OCT_vs_Fundus/CIFU_oct_vs_fundus/OCT_Fundus_NRR.csv"

localPath_glaucoma_NRR <- "~/Desktop/CIFU_OCT_vs_Fundus/CIFU_oct_vs_fundus/NRR_OCT_Fundus_Glaucoma_N105.csv"

# Set figurePath
figurePath <- "~/Desktop/CIFU_OCT_vs_Fundus/CIFU_oct_vs_fundus/figures/"

# Set filePath
filePath <- "~/Desktop/CIFU_OCT_vs_Fundus/CIFU_oct_vs_fundus"

# Set the path for tables/spreadsheets
localPath_tables <- "~/Desktop/CIFU_OCT_vs_Fundus/CIFU_oct_vs_fundus/tables/"

# load data
fullOCTdata_normal_NRR <- read.csv(localPath_normal_NRR, stringsAsFactors=TRUE)
fullOCTdata_glaucoma_NRR <- read.csv(localPath_glaucoma_NRR, stringsAsFactors=TRUE)

#################################
# remove duplicates
#################################

fullOCTdata_normal_NRR <- unique(fullOCTdata_normal_NRR)
fullOCTdata_normal_NRR <- na.omit(fullOCTdata_normal_NRR)

fullOCTdata_glaucoma_NRR <- unique(fullOCTdata_glaucoma_NRR)
fullOCTdata_glaucoma_NRR <- na.omit(fullOCTdata_glaucoma_NRR)  # row 79 has na values

###################################################################
# We are keeping the entire data set, so we will not              #
# remove any 'outliers'                                           #
###################################################################
# outlier removal
outlier.identifier <- function(curves){
  #####
  # input: 
  #         curves:   n X d matrix of either NRR or RNFL curves, where 
  #                   n is the number of obs. and d is the number of data points
  # output: 
  #         tracker:  index of outlier
  #####
  n <- dim(curves)[1]
  curves_max_mean <- max(apply(abs(curves),2,mean))
  curves_max_std <- max(apply(abs(curves),2,std))
  curves_outlier_threshold <- curves_max_mean + 3.5*curves_max_std
  tracker <- rep(NA, n)  # outlier tracker
  for( i in 1:n){
    if(max(abs(curves[i,])) > curves_outlier_threshold){
      tracker[i] <- i
    }
  }
  tracker <- na.omit(tracker)
  return(tracker)
}

# # Find and remove any outliers
# E.G.:
# NRR_Normal_OG_tracker <- outlier.identifier(NRR_Normal_OG[,-c(1:8)])
# NRR_Normal_OG_NoOut <- NRR_Normal_OG[-NRR_Normal_OG_tracker,]

# This function corrects for the error we might get if we try to 
# remove observations based on an empty logical vector
check_tracker_length_0 <- function(x,y){
  if(length(x)){  # the conditional will only execute if the logical is non-empty
    z <- y[-x,]  # if the conditional is not empty, remove the outliers
  } else{z <- y}  # else, remove nothing
  return(z)
}


###################################################################
#  Parse the data into OCT and Fundus and remove any NA values    #
###################################################################

NRR_OCT_Normal_OG <- fullOCTdata_normal_NRR[,c(c(1:8),grep("OCT",x=names(fullOCTdata_normal_NRR)))]

NRR_Fundus_Normal_OG <- fullOCTdata_normal_NRR[,c(c(1:8),grep("Fundus",x=names(fullOCTdata_normal_NRR)))]

NRR_OCT_Glau_OG <- fullOCTdata_glaucoma_NRR[,c(c(1:8),grep("OCT",x=names(fullOCTdata_glaucoma_NRR)))]

NRR_Fundus_Glau_OG <- fullOCTdata_glaucoma_NRR[,c(c(1:8),grep("Fundus",x=names(fullOCTdata_glaucoma_NRR)))]


###################################################################
# Normalize the curves

#################################################################################
### Input:  x:        numeric vector of nrr or rnfl datapoints for a single   ###
###                   observation/curve                                       ###
###                                                                           ###
###         method:   numerical integration method from the pracma            ### 
###                   package c("Kronrod", "Clenshaw","Simpson")              ###
###                                                                           ###
###               https://cran.r-project.org/web/packages/pracma/pracma.pdf   ###
###                                                                           ###
### Output:                                                                   ###
###         X_target:     normalized curve                                    ###
###         X_og_area:    area under the original curve                       ###
#################################################################################

normalize.curve = function(curves, method="Kron"){
  m = length(curves)  # number of data points along each curve
  curve_fn = function(t) curves[t] 
  X_og_area = integral(curve_fn,xmin=1,xmax=m,method = method,
                       no_intervals=m, random=FALSE, reltol=1e-8, abstol=0)
  X_target = curves/X_og_area
  return(X_target)
}

stratified.curve.norm <- function(temp_OG_curves){
  nrow_temp_curves <- dim(temp_OG_curves)[1]
  ncol_temp_curves <- dim(temp_OG_curves)[2]
  
  temp_Norm_curves <- matrix(NA, nrow = nrow_temp_curves, ncol = ncol_temp_curves, byrow = TRUE)
  
  for(i in 1:nrow_temp_curves){
    temp_Norm_curves_i <- temp_OG_curves[i,]
    temp_Norm_curves[i,] <- normalize.curve(t(temp_Norm_curves_i))
  }
  temp_Norm_curves <- apply(temp_Norm_curves, 2, as.numeric)
  return(temp_Norm_curves)
}


# 
# NRR_OCT_Normal_Norm <- matrix(NA, nrow = nrow(NRR_OCT_Normal_OG), ncol = ncol(NRR_OCT_Normal_OG), byrow = TRUE)
# NRR_Fundus_Normal_Norm <- matrix(NA, nrow = nrow(NRR_Fundus_Normal_OG), ncol = ncol(NRR_Fundus_Normal_OG), byrow = TRUE)
# NRR_OCT_Glau_Norm <- matrix(NA, nrow = nrow(NRR_OCT_Glau_OG), ncol = ncol(NRR_OCT_Glau_OG), byrow = TRUE)
# NRR_Fundus_Glau_Norm <- matrix(NA, nrow = nrow(NRR_Fundus_Glau_OG), ncol = ncol(NRR_Fundus_Glau_OG), byrow = TRUE)
# 
# # The first 8 columns are non-numeric variables; Col's 9:188 are the 180 numeric NRR/Fundus data points
# NRR_OCT_Normal_Norm[,9:ncol(NRR_OCT_Normal_OG)] <- stratified.curve.norm(NRR_OCT_Normal_OG[,-c(1:8)])
# NRR_Fundus_Normal_Norm[,9:ncol(NRR_Fundus_Normal_OG)] <- stratified.curve.norm(NRR_Fundus_Normal_OG[,-c(1:8)])
# NRR_OCT_Glau_Norm[,9:ncol(NRR_OCT_Glau_OG)] <- stratified.curve.norm(NRR_OCT_Glau_OG[,-c(1:8)])
# NRR_Fundus_Glau_Norm[,9:ncol(NRR_Fundus_Glau_OG)] <- stratified.curve.norm(NRR_Fundus_Glau_OG[,-c(1:8)])
# 
# # Convert to data frames to accommodate multiple data types
# NRR_OCT_Normal_Norm <- as.data.frame(NRR_OCT_Normal_Norm)
# NRR_Fundus_Normal_Norm <- as.data.frame(NRR_Fundus_Normal_Norm)
# NRR_OCT_Glau_Norm <- as.data.frame(NRR_OCT_Glau_Norm)
# NRR_Fundus_Glau_Norm <- as.data.frame(NRR_Fundus_Glau_Norm)
# 
# 
# NRR_OCT_Normal_Norm[,1:8] <- fullOCTdata_normal_NRR[,1:8]
# NRR_Fundus_Normal_Norm[,1:8] <- fullOCTdata_normal_NRR[,1:8]
# 
# NRR_OCT_Glau_Norm[,1:8] <- fullOCTdata_glaucoma_NRR[,1:8]
# NRR_Fundus_Glau_Norm[,1:8] <- fullOCTdata_glaucoma_NRR[,1:8]
# 
# 
# NRR_OCT_Normal_OG_path <- paste(filePath,"/rdata_files/NRR_OCT_Normal_OG.RData",sep="")
# save(NRR_OCT_Normal_OG, file = NRR_OCT_Normal_OG_path)
# load(NRR_OCT_Normal_OG_path)
# 
# NRR_OCT_Normal_Norm_path <- paste(filePath,"/rdata_files/NRR_OCT_Normal_Norm.RData",sep="")
# save(NRR_OCT_Normal_Norm, file = NRR_OCT_Normal_Norm_path)
# load(NRR_OCT_Normal_Norm_path)
# 
# NRR_Fundus_Normal_Norm_path <- paste(filePath,"/rdata_files/NRR_Fundus_Normal_Norm.RData",sep="")
# save(NRR_Fundus_Normal_Norm, file = NRR_Fundus_Normal_Norm_path)
# load(NRR_Fundus_Normal_Norm_path)
# 
# NRR_OCT_Glau_OG_path <- paste(filePath,"/rdata_files/NRR_OCT_Glau_OG.RData",sep="")
# save(NRR_OCT_Glau_OG, file = NRR_OCT_Glau_OG_path)
# load(NRR_OCT_Glau_OG_path)
# 
# NRR_OCT_Glau_Norm_path <- paste(filePath,"/rdata_files/NRR_OCT_Glau_Norm.RData",sep="")
# save(NRR_OCT_Glau_Norm, file = NRR_OCT_Glau_Norm_path)
# load(NRR_OCT_Glau_Norm_path)
# 
# NRR_Fundus_Glau_Norm_path <- paste(filePath,"/rdata_files/NRR_Fundus_Glau_Norm.RData",sep="")
# save(NRR_Fundus_Glau_Norm, file = NRR_Fundus_Glau_Norm_path)
# load(NRR_Fundus_Glau_Norm_path)
# 
# NRR_Fundus_Glau_OG_path <- paste(filePath,"/rdata_files/NRR_Fundus_Glau_OG.RData",sep="")
# save(NRR_Fundus_Glau_OG, file = NRR_Fundus_Glau_OG_path)
# load(NRR_Fundus_Glau_OG_path)



################# BEGIN FUNCTIONAL DATA ANALYSIS ###################################

load(paste(filePath,'/rdata_files/NRR_OCT_Normal_Norm.RData',sep=""))
load(paste(filePath,'/rdata_files/NRR_Fundus_Normal_Norm.RData',sep=""))
load(paste(filePath,'/rdata_files/NRR_OCT_Glau_Norm.RData',sep=""))
load(paste(filePath,'/rdata_files/NRR_Fundus_Glau_Norm.RData',sep=""))

###################################################################
# Generate bases
# sequence of points for for plotting
d_NRR <- seq(from=1, to=180)  # 180 data points for RNFL curves (every other deg.)


# Generate Fourier bases for NRR curves
NRR_basis_03 <- create.fourier.basis(c(1, 180), nbasis=3)
NRR_basis_05 <- create.fourier.basis(c(1, 180), nbasis=5)
NRR_basis_07 <- create.fourier.basis(c(1, 180), nbasis=7)
NRR_basis_09 <- create.fourier.basis(c(1, 180), nbasis=9)
NRR_basis_11 <- create.fourier.basis(c(1, 180), nbasis=11)
NRR_basis_13 <- create.fourier.basis(c(1, 180), nbasis=13)
NRR_basis_15 <- create.fourier.basis(c(1, 180), nbasis=15)
NRR_basis_17 <- create.fourier.basis(c(1, 180), nbasis=17)
NRR_basis_19 <- create.fourier.basis(c(1, 180), nbasis=19)
NRR_basis_21 <- create.fourier.basis(c(1, 180), nbasis=21)

NRR_basis_p <- list(NRR_basis_03, NRR_basis_05,
                    NRR_basis_07, NRR_basis_09,
                    NRR_basis_11, NRR_basis_13,
                    NRR_basis_15, NRR_basis_17,
                    NRR_basis_19, NRR_basis_21)

###################################################################
# Generate fd Objects


# NRR_OCT_Normal_OG
fdNRR_OCT_Normal_OG_p03 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Normal_OG[,-c(1:8)]), basisobj=NRR_basis_03)
fdNRR_OCT_Normal_OG_p05 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Normal_OG[,-c(1:8)]), basisobj=NRR_basis_05)
fdNRR_OCT_Normal_OG_p07 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Normal_OG[,-c(1:8)]), basisobj=NRR_basis_07)
fdNRR_OCT_Normal_OG_p09 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Normal_OG[,-c(1:8)]), basisobj=NRR_basis_09)
fdNRR_OCT_Normal_OG_p11 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Normal_OG[,-c(1:8)]), basisobj=NRR_basis_11)
fdNRR_OCT_Normal_OG_p13 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Normal_OG[,-c(1:8)]), basisobj=NRR_basis_13)
fdNRR_OCT_Normal_OG_p15 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Normal_OG[,-c(1:8)]), basisobj=NRR_basis_15)
fdNRR_OCT_Normal_OG_p17 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Normal_OG[,-c(1:8)]), basisobj=NRR_basis_17)
fdNRR_OCT_Normal_OG_p19 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Normal_OG[,-c(1:8)]), basisobj=NRR_basis_19)
fdNRR_OCT_Normal_OG_p21 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Normal_OG[,-c(1:8)]), basisobj=NRR_basis_21)

fdNRR_OCT_Normal_OG_p <- list(fdNRR_OCT_Normal_OG_p03, fdNRR_OCT_Normal_OG_p05, fdNRR_OCT_Normal_OG_p07, fdNRR_OCT_Normal_OG_p09, 
                        fdNRR_OCT_Normal_OG_p11, fdNRR_OCT_Normal_OG_p13, fdNRR_OCT_Normal_OG_p15, fdNRR_OCT_Normal_OG_p17, 
                        fdNRR_OCT_Normal_OG_p19, fdNRR_OCT_Normal_OG_p21)


# NRR_Fundus_Normal_OG
fdNRR_Fundus_Normal_OG_p03 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Normal_OG[,-c(1:8)]), basisobj=NRR_basis_03)
fdNRR_Fundus_Normal_OG_p05 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Normal_OG[,-c(1:8)]), basisobj=NRR_basis_05)
fdNRR_Fundus_Normal_OG_p07 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Normal_OG[,-c(1:8)]), basisobj=NRR_basis_07)
fdNRR_Fundus_Normal_OG_p09 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Normal_OG[,-c(1:8)]), basisobj=NRR_basis_09)
fdNRR_Fundus_Normal_OG_p11 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Normal_OG[,-c(1:8)]), basisobj=NRR_basis_11)
fdNRR_Fundus_Normal_OG_p13 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Normal_OG[,-c(1:8)]), basisobj=NRR_basis_13)
fdNRR_Fundus_Normal_OG_p15 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Normal_OG[,-c(1:8)]), basisobj=NRR_basis_15)
fdNRR_Fundus_Normal_OG_p17 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Normal_OG[,-c(1:8)]), basisobj=NRR_basis_17)
fdNRR_Fundus_Normal_OG_p19 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Normal_OG[,-c(1:8)]), basisobj=NRR_basis_19)
fdNRR_Fundus_Normal_OG_p21 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Normal_OG[,-c(1:8)]), basisobj=NRR_basis_21)

fdNRR_Fundus_Normal_OG_p <- list(fdNRR_Fundus_Normal_OG_p03, fdNRR_Fundus_Normal_OG_p05, fdNRR_Fundus_Normal_OG_p07, fdNRR_Fundus_Normal_OG_p09, 
                              fdNRR_Fundus_Normal_OG_p11, fdNRR_Fundus_Normal_OG_p13, fdNRR_Fundus_Normal_OG_p15, fdNRR_Fundus_Normal_OG_p17, 
                              fdNRR_Fundus_Normal_OG_p19, fdNRR_Fundus_Normal_OG_p21)


# NRR_OCT_Glau_OG
fdNRR_OCT_Glau_OG_p03 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Glau_OG[,-c(1:8)]), basisobj=NRR_basis_03)
fdNRR_OCT_Glau_OG_p05 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Glau_OG[,-c(1:8)]), basisobj=NRR_basis_05)
fdNRR_OCT_Glau_OG_p07 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Glau_OG[,-c(1:8)]), basisobj=NRR_basis_07)
fdNRR_OCT_Glau_OG_p09 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Glau_OG[,-c(1:8)]), basisobj=NRR_basis_09)
fdNRR_OCT_Glau_OG_p11 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Glau_OG[,-c(1:8)]), basisobj=NRR_basis_11)
fdNRR_OCT_Glau_OG_p13 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Glau_OG[,-c(1:8)]), basisobj=NRR_basis_13)
fdNRR_OCT_Glau_OG_p15 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Glau_OG[,-c(1:8)]), basisobj=NRR_basis_15)
fdNRR_OCT_Glau_OG_p17 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Glau_OG[,-c(1:8)]), basisobj=NRR_basis_17)
fdNRR_OCT_Glau_OG_p19 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Glau_OG[,-c(1:8)]), basisobj=NRR_basis_19)
fdNRR_OCT_Glau_OG_p21 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Glau_OG[,-c(1:8)]), basisobj=NRR_basis_21)

fdNRR_OCT_Glau_OG_p <- list(fdNRR_OCT_Glau_OG_p03, fdNRR_OCT_Glau_OG_p05, fdNRR_OCT_Glau_OG_p07, fdNRR_OCT_Glau_OG_p09, 
                              fdNRR_OCT_Glau_OG_p11, fdNRR_OCT_Glau_OG_p13, fdNRR_OCT_Glau_OG_p15, fdNRR_OCT_Glau_OG_p17, 
                              fdNRR_OCT_Glau_OG_p19, fdNRR_OCT_Glau_OG_p21)


# NRR_Fundus_Glau_OG
fdNRR_Fundus_Glau_OG_p03 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Glau_OG[,-c(1:8)]), basisobj=NRR_basis_03)
fdNRR_Fundus_Glau_OG_p05 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Glau_OG[,-c(1:8)]), basisobj=NRR_basis_05)
fdNRR_Fundus_Glau_OG_p07 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Glau_OG[,-c(1:8)]), basisobj=NRR_basis_07)
fdNRR_Fundus_Glau_OG_p09 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Glau_OG[,-c(1:8)]), basisobj=NRR_basis_09)
fdNRR_Fundus_Glau_OG_p11 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Glau_OG[,-c(1:8)]), basisobj=NRR_basis_11)
fdNRR_Fundus_Glau_OG_p13 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Glau_OG[,-c(1:8)]), basisobj=NRR_basis_13)
fdNRR_Fundus_Glau_OG_p15 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Glau_OG[,-c(1:8)]), basisobj=NRR_basis_15)
fdNRR_Fundus_Glau_OG_p17 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Glau_OG[,-c(1:8)]), basisobj=NRR_basis_17)
fdNRR_Fundus_Glau_OG_p19 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Glau_OG[,-c(1:8)]), basisobj=NRR_basis_19)
fdNRR_Fundus_Glau_OG_p21 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Glau_OG[,-c(1:8)]), basisobj=NRR_basis_21)

fdNRR_Fundus_Glau_OG_p <- list(fdNRR_Fundus_Glau_OG_p03, fdNRR_Fundus_Glau_OG_p05, fdNRR_Fundus_Glau_OG_p07, fdNRR_Fundus_Glau_OG_p09, 
                                 fdNRR_Fundus_Glau_OG_p11, fdNRR_Fundus_Glau_OG_p13, fdNRR_Fundus_Glau_OG_p15, fdNRR_Fundus_Glau_OG_p17, 
                                 fdNRR_Fundus_Glau_OG_p19, fdNRR_Fundus_Glau_OG_p21)



# NRR_OCT_Normal_Norm
fdNRR_OCT_Normal_Norm_p03 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Normal_Norm[,-c(1:8)]), basisobj=NRR_basis_03)
fdNRR_OCT_Normal_Norm_p05 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Normal_Norm[,-c(1:8)]), basisobj=NRR_basis_05)
fdNRR_OCT_Normal_Norm_p07 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Normal_Norm[,-c(1:8)]), basisobj=NRR_basis_07)
fdNRR_OCT_Normal_Norm_p09 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Normal_Norm[,-c(1:8)]), basisobj=NRR_basis_09)
fdNRR_OCT_Normal_Norm_p11 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Normal_Norm[,-c(1:8)]), basisobj=NRR_basis_11)
fdNRR_OCT_Normal_Norm_p13 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Normal_Norm[,-c(1:8)]), basisobj=NRR_basis_13)
fdNRR_OCT_Normal_Norm_p15 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Normal_Norm[,-c(1:8)]), basisobj=NRR_basis_15)
fdNRR_OCT_Normal_Norm_p17 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Normal_Norm[,-c(1:8)]), basisobj=NRR_basis_17)
fdNRR_OCT_Normal_Norm_p19 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Normal_Norm[,-c(1:8)]), basisobj=NRR_basis_19)
fdNRR_OCT_Normal_Norm_p21 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Normal_Norm[,-c(1:8)]), basisobj=NRR_basis_21)

fdNRR_OCT_Normal_Norm_p <- list(fdNRR_OCT_Normal_Norm_p03, fdNRR_OCT_Normal_Norm_p05, fdNRR_OCT_Normal_Norm_p07, fdNRR_OCT_Normal_Norm_p09, 
                              fdNRR_OCT_Normal_Norm_p11, fdNRR_OCT_Normal_Norm_p13, fdNRR_OCT_Normal_Norm_p15, fdNRR_OCT_Normal_Norm_p17, 
                              fdNRR_OCT_Normal_Norm_p19, fdNRR_OCT_Normal_Norm_p21)


# NRR_Fundus_Normal_Norm
fdNRR_Fundus_Normal_Norm_p03 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Normal_Norm[,-c(1:8)]), basisobj=NRR_basis_03)
fdNRR_Fundus_Normal_Norm_p05 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Normal_Norm[,-c(1:8)]), basisobj=NRR_basis_05)
fdNRR_Fundus_Normal_Norm_p07 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Normal_Norm[,-c(1:8)]), basisobj=NRR_basis_07)
fdNRR_Fundus_Normal_Norm_p09 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Normal_Norm[,-c(1:8)]), basisobj=NRR_basis_09)
fdNRR_Fundus_Normal_Norm_p11 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Normal_Norm[,-c(1:8)]), basisobj=NRR_basis_11)
fdNRR_Fundus_Normal_Norm_p13 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Normal_Norm[,-c(1:8)]), basisobj=NRR_basis_13)
fdNRR_Fundus_Normal_Norm_p15 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Normal_Norm[,-c(1:8)]), basisobj=NRR_basis_15)
fdNRR_Fundus_Normal_Norm_p17 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Normal_Norm[,-c(1:8)]), basisobj=NRR_basis_17)
fdNRR_Fundus_Normal_Norm_p19 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Normal_Norm[,-c(1:8)]), basisobj=NRR_basis_19)
fdNRR_Fundus_Normal_Norm_p21 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Normal_Norm[,-c(1:8)]), basisobj=NRR_basis_21)

fdNRR_Fundus_Normal_Norm_p <- list(fdNRR_Fundus_Normal_Norm_p03, fdNRR_Fundus_Normal_Norm_p05, fdNRR_Fundus_Normal_Norm_p07, fdNRR_Fundus_Normal_Norm_p09, 
                                 fdNRR_Fundus_Normal_Norm_p11, fdNRR_Fundus_Normal_Norm_p13, fdNRR_Fundus_Normal_Norm_p15, fdNRR_Fundus_Normal_Norm_p17, 
                                 fdNRR_Fundus_Normal_Norm_p19, fdNRR_Fundus_Normal_Norm_p21)



# NRR_OCT_Glau_Norm
fdNRR_OCT_Glau_Norm_p03 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Glau_Norm[,-c(1:8)]), basisobj=NRR_basis_03)
fdNRR_OCT_Glau_Norm_p05 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Glau_Norm[,-c(1:8)]), basisobj=NRR_basis_05)
fdNRR_OCT_Glau_Norm_p07 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Glau_Norm[,-c(1:8)]), basisobj=NRR_basis_07)
fdNRR_OCT_Glau_Norm_p09 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Glau_Norm[,-c(1:8)]), basisobj=NRR_basis_09)
fdNRR_OCT_Glau_Norm_p11 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Glau_Norm[,-c(1:8)]), basisobj=NRR_basis_11)
fdNRR_OCT_Glau_Norm_p13 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Glau_Norm[,-c(1:8)]), basisobj=NRR_basis_13)
fdNRR_OCT_Glau_Norm_p15 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Glau_Norm[,-c(1:8)]), basisobj=NRR_basis_15)
fdNRR_OCT_Glau_Norm_p17 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Glau_Norm[,-c(1:8)]), basisobj=NRR_basis_17)
fdNRR_OCT_Glau_Norm_p19 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Glau_Norm[,-c(1:8)]), basisobj=NRR_basis_19)
fdNRR_OCT_Glau_Norm_p21 <- Data2fd(argvals=d_NRR, y=t(NRR_OCT_Glau_Norm[,-c(1:8)]), basisobj=NRR_basis_21)

fdNRR_OCT_Glau_Norm_p <- list(fdNRR_OCT_Glau_Norm_p03, fdNRR_OCT_Glau_Norm_p05, fdNRR_OCT_Glau_Norm_p07, fdNRR_OCT_Glau_Norm_p09, 
                            fdNRR_OCT_Glau_Norm_p11, fdNRR_OCT_Glau_Norm_p13, fdNRR_OCT_Glau_Norm_p15, fdNRR_OCT_Glau_Norm_p17, 
                            fdNRR_OCT_Glau_Norm_p19, fdNRR_OCT_Glau_Norm_p21)


# NRR_Fundus_Glau_Norm
fdNRR_Fundus_Glau_Norm_p03 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Glau_Norm[,-c(1:8)]), basisobj=NRR_basis_03)
fdNRR_Fundus_Glau_Norm_p05 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Glau_Norm[,-c(1:8)]), basisobj=NRR_basis_05)
fdNRR_Fundus_Glau_Norm_p07 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Glau_Norm[,-c(1:8)]), basisobj=NRR_basis_07)
fdNRR_Fundus_Glau_Norm_p09 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Glau_Norm[,-c(1:8)]), basisobj=NRR_basis_09)
fdNRR_Fundus_Glau_Norm_p11 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Glau_Norm[,-c(1:8)]), basisobj=NRR_basis_11)
fdNRR_Fundus_Glau_Norm_p13 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Glau_Norm[,-c(1:8)]), basisobj=NRR_basis_13)
fdNRR_Fundus_Glau_Norm_p15 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Glau_Norm[,-c(1:8)]), basisobj=NRR_basis_15)
fdNRR_Fundus_Glau_Norm_p17 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Glau_Norm[,-c(1:8)]), basisobj=NRR_basis_17)
fdNRR_Fundus_Glau_Norm_p19 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Glau_Norm[,-c(1:8)]), basisobj=NRR_basis_19)
fdNRR_Fundus_Glau_Norm_p21 <- Data2fd(argvals=d_NRR, y=t(NRR_Fundus_Glau_Norm[,-c(1:8)]), basisobj=NRR_basis_21)

fdNRR_Fundus_Glau_Norm_p <- list(fdNRR_Fundus_Glau_Norm_p03, fdNRR_Fundus_Glau_Norm_p05, fdNRR_Fundus_Glau_Norm_p07, fdNRR_Fundus_Glau_Norm_p09, 
                               fdNRR_Fundus_Glau_Norm_p11, fdNRR_Fundus_Glau_Norm_p13, fdNRR_Fundus_Glau_Norm_p15, fdNRR_Fundus_Glau_Norm_p17, 
                               fdNRR_Fundus_Glau_Norm_p19, fdNRR_Fundus_Glau_Norm_p21)



###################################################################
###################################################################
### save results for p = ?? for use in ModuleB.R

######################################
# Calculate FVE

sqr.diff.calc <- function(curves){
  ## Takes either NRR or RNFL curves (without PatientID_EYE)
  n <- dim(curves)[1]
  nmb_pts <- dim(curves)[2]
  # average curve values pointwise
  curves_X_Bar = apply(curves, 2, sum)/n
  # integrand
  curves_sqrErr = function(t) (X_i[t] - curves_X_Bar[t])^2
  # output matrix
  curves_sqr_diff <- rep(NA,n)
  for(i in 1:n){
    X_i <- curves[i,]
    curves_sqr_diff[i] <- integral(curves_sqrErr, xmin=1, xmax=nmb_pts, method = "Kron",
                                   no_intervals=nmb_pts, random=FALSE, reltol=1e-8, abstol=0)
    if(i%%100==0) print(i)  # progress tracker
  }
  return(curves_sqr_diff)
}

###################################################################
#####                   Calculate the square diff             #####
###################################################################
# This step takes a significant amount of time, so we will save 
# once the calculations have been processed.
###################################################################
# NRR_OCT_Normal_OG_sqr_diff <- sqr.diff.calc(NRR_OCT_Normal_OG[,-c(1:8)])
NRR_OCT_Normal_OG_sqr_diff_path <- paste(filePath,"/rdata_files/NRR_OCT_Normal_OG_sqr_diff.RData",sep="")
# save(NRR_OCT_Normal_OG_sqr_diff, file = NRR_OCT_Normal_OG_sqr_diff_path)
load(NRR_OCT_Normal_OG_sqr_diff_path)

# NRR_Fundus_Normal_OG_sqr_diff <- sqr.diff.calc(NRR_Fundus_Normal_OG[,-c(1:8)])
NRR_Fundus_Normal_OG_sqr_diff_path <- paste(filePath,"/rdata_files/NRR_Fundus_Normal_OG_sqr_diff.RData",sep="")
# save(NRR_Fundus_Normal_OG_sqr_diff, file = NRR_Fundus_Normal_OG_sqr_diff_path)
load(NRR_Fundus_Normal_OG_sqr_diff_path)

# NRR_OCT_Glau_OG_sqr_diff <- sqr.diff.calc(NRR_OCT_Glau_OG[-79,-c(1:8)])
NRR_OCT_Glau_OG_sqr_diff_path <- paste(filePath,"/rdata_files/NRR_OCT_Glau_OG_sqr_diff.RData",sep="")
# save(NRR_OCT_Glau_OG_sqr_diff, file = NRR_OCT_Glau_OG_sqr_diff_path)
load(NRR_OCT_Glau_OG_sqr_diff_path)

# NRR_Fundus_Glau_OG_sqr_diff <- sqr.diff.calc(NRR_Fundus_Glau_OG[-79,-c(1:8)])
NRR_Fundus_Glau_OG_sqr_diff_path <- paste(filePath,"/rdata_files/NRR_Fundus_Glau_OG_sqr_diff.RData",sep="")
# save(NRR_Fundus_Glau_OG_sqr_diff, file = NRR_Fundus_Glau_OG_sqr_diff_path)
load(NRR_Fundus_Glau_OG_sqr_diff_path)

# NRR_OCT_Normal_Norm_sqr_diff <- sqr.diff.calc(NRR_OCT_Normal_Norm[,-c(1:8)])
NRR_OCT_Normal_Norm_sqr_diff_path <- paste(filePath,"/rdata_files/NRR_OCT_Normal_Norm_sqr_diff.RData",sep="")
# save(NRR_OCT_Normal_Norm_sqr_diff, file = NRR_OCT_Normal_Norm_sqr_diff_path)
load(NRR_OCT_Normal_Norm_sqr_diff_path)

# NRR_Fundus_Normal_Norm_sqr_diff <- sqr.diff.calc(NRR_Fundus_Normal_Norm[,-c(1:8)])
NRR_Fundus_Normal_Norm_sqr_diff_path <- paste(filePath,"/rdata_files/NRR_Fundus_Normal_Norm_sqr_diff.RData",sep="")
# save(NRR_Fundus_Normal_Norm_sqr_diff, file = NRR_Fundus_Normal_Norm_sqr_diff_path)
load(NRR_Fundus_Normal_Norm_sqr_diff_path)

# NRR_OCT_Glau_Norm_sqr_diff <- sqr.diff.calc(NRR_OCT_Glau_Norm[-79,-c(1:8)])
NRR_OCT_Glau_Norm_sqr_diff_path <- paste(filePath,"/rdata_files/NRR_OCT_Glau_Norm_sqr_diff.RData",sep="")
# save(NRR_OCT_Glau_Norm_sqr_diff, file = NRR_OCT_Glau_Norm_sqr_diff_path)
load(NRR_OCT_Glau_Norm_sqr_diff_path)

# NRR_Fundus_Glau_Norm_sqr_diff <- sqr.diff.calc(NRR_Fundus_Glau_Norm[-79,-c(1:8)])
NRR_Fundus_Glau_Norm_sqr_diff_path <- paste(filePath,"/rdata_files/NRR_Fundus_Glau_Norm_sqr_diff.RData",sep="")
# save(NRR_Fundus_Glau_Norm_sqr_diff, file = NRR_Fundus_Glau_Norm_sqr_diff_path)
load(NRR_Fundus_Glau_Norm_sqr_diff_path)


FVE.samp.var.p = function(X_i_p_t, numEyDt, method="Kron"){
  
  #################################################################################
  ###         Function for calculating the Fraction of Variance Explained       ###
  ###                                                                           ###
  ### Input:  numEyDt:  numeric matrix, populated by the retinal eye data       ###
  ###                                                                           ###
  ###         method:   numerical integration function from the pracma          ### 
  ###                   package with methods c("Kronrod", "Clenshaw","Simpson") ###
  ###                                                                           ###
  ###               https://cran.r-project.org/web/packages/pracma/pracma.pdf   ###
  ###                                                                           ###
  ###         X_i_p_t:  numeric matrix, populated by approximated curve values  ###
  ###                   for a range of basis functions.                         ###
  ###                                                                           ###
  ### Output: Sp2:    vector containing sample variance for each value of p     ###
  ###                                                                           ###
  #################################################################################
  
  n = dim(X_i_p_t)[1]  # number of observations
  m = dim(X_i_p_t)[2]  # number of data points along each curve
  
  # integrand
  sqrErr2 = function(t) (X_i_p[t] - X_i[t])^2  # square error function
  
  sqr_diff <- matrix(data = NA, nrow = n, byrow = TRUE)  # (n X 1) vector for collecting values
  for(i in 1:n){
    # integrate the square error function over the range for each curve
    X_i <- numEyDt[i,]
    X_i_p <- X_i_p_t[i,]
    sqr_diff[i] <- integral(sqrErr2, xmin = 1, xmax = m, method = method,
                            no_intervals=m, random=FALSE, reltol=1e-8, abstol=0)
  }
  # sample curve variation for p basis functions
  Sp2 = sum(sqr_diff)/(n-1)  
  return(Sp2)
}



# total variation

# 
# 
# NRR_OCT_Normal_OG_TV <- sum(NRR_OCT_Normal_OG_sqr_diff)/(length(NRR_OCT_Normal_OG_sqr_diff)-1)
# NRR_Fundus_Normal_OG_TV <- sum(NRR_Fundus_Normal_OG_sqr_diff)/(length(NRR_Fundus_Normal_OG_sqr_diff)-1)
# NRR_OCT_Glau_OG_TV <- sum(NRR_OCT_Glau_OG_sqr_diff)/(length(NRR_OCT_Glau_OG_sqr_diff)-1) 
# NRR_Fundus_Glau_OG_TV <- sum(NRR_Fundus_Glau_OG_sqr_diff)/(length(NRR_Fundus_Glau_OG_sqr_diff)-1)
# NRR_OCT_Normal_Norm_TV <- sum(NRR_OCT_Normal_Norm_sqr_diff)/(length(NRR_OCT_Normal_Norm_sqr_diff)-1)
# NRR_Fundus_Normal_Norm_TV <- sum(NRR_Fundus_Normal_Norm_sqr_diff)/(length(NRR_Fundus_Normal_Norm_sqr_diff)-1)
# NRR_OCT_Glau_Norm_TV <- sum(NRR_OCT_Glau_Norm_sqr_diff)/(length(NRR_OCT_Glau_Norm_sqr_diff)-1)
# NRR_Fundus_Glau_Norm_TV <- sum(NRR_Fundus_Glau_Norm_sqr_diff)/(length(NRR_Fundus_Glau_Norm_sqr_diff)-1)
# 
# 
# p_seq_len = length(NRR_basis_p)
# 
# NRR_OCT_Normal_OG_p_basis_var <- rep(NA,p_seq_len)
# NRR_OCT_Normal_OG_FVE <- rep(NA,p_seq_len)
# NRR_Fundus_Normal_OG_p_basis_var <- rep(NA,p_seq_len)
# NRR_Fundus_Normal_OG_FVE <- rep(NA,p_seq_len)
# NRR_OCT_Glau_OG_p_basis_var <- rep(NA,p_seq_len)
# NRR_OCT_Glau_OG_FVE <- rep(NA,p_seq_len)
# NRR_Fundus_Glau_OG_p_basis_var <- rep(NA,p_seq_len)
# NRR_Fundus_Glau_OG_FVE <- rep(NA,p_seq_len)
# 
# NRR_OCT_Normal_Norm_p_basis_var <- rep(NA,p_seq_len)
# NRR_OCT_Normal_Norm_FVE <- rep(NA,p_seq_len)
# NRR_Fundus_Normal_Norm_p_basis_var <- rep(NA,p_seq_len)
# NRR_Fundus_Normal_Norm_FVE <- rep(NA,p_seq_len)
# NRR_OCT_Glau_Norm_p_basis_var <- rep(NA,p_seq_len)
# NRR_OCT_Glau_Norm_FVE <- rep(NA,p_seq_len)
# NRR_Fundus_Glau_Norm_p_basis_var <- rep(NA,p_seq_len)
# NRR_Fundus_Glau_Norm_FVE <- rep(NA,p_seq_len)


###################################################################
# These each take 9-13 hours per iteration of each loop, and
# there are 10 iterations in each sequence
###################################################################
# for(p in 1:p_seq_len){
#   startTime <- Sys.time()
#   print(p)
#   X_ij_p_temp = t(eval.fd(d_NRR, fdNRR_OCT_Normal_OG_p[[p]]))
#   NRR_OCT_Normal_OG_p_basis_var[p] <- FVE.samp.var.p(X_ij_p_temp, NRR_OCT_Normal_OG[,-c(1:8)], method="Kron")
#   NRR_OCT_Normal_OG_FVE[p] <- (NRR_OCT_Normal_OG_TV-NRR_OCT_Normal_OG_p_basis_var[p])/NRR_OCT_Normal_OG_TV
#   endTime <- Sys.time()
#   print(endTime - startTime)
# }

NRR_OCT_Normal_OG_p_basis_var_path <- paste(filePath,"/rdata_files/NRR_OCT_Normal_OG_p_basis_var.RData",sep="")
# save(NRR_OCT_Normal_OG_p_basis_var, file = NRR_OCT_Normal_OG_p_basis_var_path)
load(NRR_OCT_Normal_OG_p_basis_var_path)

NRR_OCT_Normal_OG_FVE_path <- paste(filePath,"/rdata_files/NRR_OCT_Normal_OG_FVE.RData",sep="")
# save(NRR_OCT_Normal_OG_FVE, file = NRR_OCT_Normal_OG_FVE_path)
load(NRR_OCT_Normal_OG_FVE_path)

# for(p in 1:p_seq_len){
#   startTime <- Sys.time()
#   print(p)
#   X_ij_p_temp = t(eval.fd(d_NRR, fdNRR_Fundus_Normal_OG_p[[p]]))
#   NRR_Fundus_Normal_OG_p_basis_var[p] <- FVE.samp.var.p(X_ij_p_temp, NRR_Fundus_Normal_OG[,-c(1:8)], method="Kron")
#   NRR_Fundus_Normal_OG_FVE[p] <- (NRR_Fundus_Normal_OG_TV-NRR_Fundus_Normal_OG_p_basis_var[p])/NRR_Fundus_Normal_OG_TV
#   endTime <- Sys.time()
#   print(endTime - startTime)
# }

NRR_Fundus_Normal_OG_p_basis_var_path <- paste(filePath,"/rdata_files/NRR_Fundus_Normal_OG_p_basis_var.RData",sep="")
# save(NRR_Fundus_Normal_OG_p_basis_var, file = NRR_Fundus_Normal_OG_p_basis_var_path)
load(NRR_Fundus_Normal_OG_p_basis_var_path)

NRR_Fundus_Normal_OG_FVE_path <- paste(filePath,"/rdata_files/NRR_Fundus_Normal_OG_FVE.RData",sep="")
# save(NRR_Fundus_Normal_OG_FVE, file = NRR_Fundus_Normal_OG_FVE_path)
load(NRR_Fundus_Normal_OG_FVE_path)


# for(p in 1:p_seq_len){
#   startTime <- Sys.time()
#   print(p)
#   X_ij_p_temp = t(eval.fd(d_NRR, fdNRR_OCT_Glau_OG_p[[p]]))
#   NRR_OCT_Glau_OG_p_basis_var[p] <- FVE.samp.var.p(X_ij_p_temp, NRR_OCT_Glau_OG[,-c(1:8)], method="Kron")
#   NRR_OCT_Glau_OG_FVE[p] <- (NRR_OCT_Glau_OG_TV-NRR_OCT_Glau_OG_p_basis_var[p])/NRR_OCT_Glau_OG_TV
#   endTime <- Sys.time()
#   print(endTime - startTime)
# }

NRR_OCT_Glau_OG_p_basis_var_path <- paste(filePath,"/rdata_files/NRR_OCT_Glau_OG_p_basis_var.RData",sep="")
# save(NRR_OCT_Glau_OG_p_basis_var, file = NRR_OCT_Glau_OG_p_basis_var_path)
load(NRR_OCT_Glau_OG_p_basis_var_path)

NRR_OCT_Glau_OG_FVE_path <- paste(filePath,"/rdata_files/NRR_OCT_Glau_OG_FVE.RData",sep="")
# save(NRR_OCT_Glau_OG_FVE, file = NRR_OCT_Glau_OG_FVE_path)
load(NRR_OCT_Glau_OG_FVE_path)

# for(p in 1:p_seq_len){
#   startTime <- Sys.time()
#   print(p)
#   X_ij_p_temp = t(eval.fd(d_NRR, fdNRR_Fundus_Glau_OG_p[[p]]))
#   NRR_Fundus_Glau_OG_p_basis_var[p] <- FVE.samp.var.p(X_ij_p_temp, NRR_Fundus_Glau_OG[,-c(1:8)], method="Kron")
#   NRR_Fundus_Glau_OG_FVE[p] <- (NRR_Fundus_Glau_OG_TV-NRR_Fundus_Glau_OG_p_basis_var[p])/NRR_Fundus_Glau_OG_TV
#   endTime <- Sys.time()
#   print(endTime - startTime)
# }

NRR_Fundus_Glau_OG_p_basis_var_path <- paste(filePath,"/rdata_files/NRR_Fundus_Glau_OG_p_basis_var.RData",sep="")
# save(NRR_Fundus_Glau_OG_p_basis_var, file = NRR_Fundus_Glau_OG_p_basis_var_path)
load(NRR_Fundus_Glau_OG_p_basis_var_path)

NRR_Fundus_Glau_OG_FVE_path <- paste(filePath,"/rdata_files/NRR_Fundus_Glau_OG_FVE.RData",sep="")
# save(NRR_Fundus_Glau_OG_FVE, file = NRR_Fundus_Glau_OG_FVE_path)
load(NRR_Fundus_Glau_OG_FVE_path)

###################################################################

# for(p in 1:p_seq_len){
#   print(p)
#   X_ij_p_temp = t(eval.fd(d_NRR, fdNRR_OCT_Normal_Norm_p[[p]]))
#   NRR_OCT_Normal_Norm_p_basis_var[p] <- FVE.samp.var.p(X_ij_p_temp, NRR_OCT_Normal_Norm[,-c(1:8)], method="Kron")
#   NRR_OCT_Normal_Norm_FVE[p] <- (NRR_OCT_Normal_Norm_TV-NRR_OCT_Normal_Norm_p_basis_var[p])/NRR_OCT_Normal_Norm_TV
# }

NRR_OCT_Normal_Norm_p_basis_var_path <- paste(filePath,"/rdata_files/NRR_OCT_Normal_Norm_p_basis_var.RData",sep="")
# save(NRR_OCT_Normal_Norm_p_basis_var, file = NRR_OCT_Normal_Norm_p_basis_var_path)
load(NRR_OCT_Normal_Norm_p_basis_var_path)

NRR_OCT_Normal_Norm_FVE_path <- paste(filePath,"/rdata_files/NRR_OCT_Normal_Norm_FVE.RData",sep="")
# save(NRR_OCT_Normal_Norm_FVE, file = NRR_OCT_Normal_Norm_FVE_path)
load(NRR_OCT_Normal_Norm_FVE_path)

# for(p in 1:p_seq_len){
#   print(p)
#   X_ij_p_temp = t(eval.fd(d_NRR, fdNRR_Fundus_Normal_Norm_p[[p]]))
#   NRR_Fundus_Normal_Norm_p_basis_var[p] <- FVE.samp.var.p(X_ij_p_temp, NRR_Fundus_Normal_Norm[,-c(1:8)], method="Kron")
#   NRR_Fundus_Normal_Norm_FVE[p] <- (NRR_Fundus_Normal_Norm_TV-NRR_Fundus_Normal_Norm_p_basis_var[p])/NRR_Fundus_Normal_Norm_TV
# }

NRR_Fundus_Normal_Norm_p_basis_var_path <- paste(filePath,"/rdata_files/NRR_Fundus_Normal_Norm_p_basis_var.RData",sep="")
# save(NRR_Fundus_Normal_Norm_p_basis_var, file = NRR_Fundus_Normal_Norm_p_basis_var_path)
load(NRR_Fundus_Normal_Norm_p_basis_var_path)

NRR_Fundus_Normal_Norm_FVE_path <- paste(filePath,"/rdata_files/NRR_Fundus_Normal_Norm_FVE.RData",sep="")
# save(NRR_Fundus_Normal_Norm_FVE, file = NRR_Fundus_Normal_Norm_FVE_path)
load(NRR_Fundus_Normal_Norm_FVE_path)

# for(p in 1:p_seq_len){
#   print(p)
#   X_ij_p_temp = t(eval.fd(d_NRR, fdNRR_OCT_Glau_Norm_p[[p]]))
#   NRR_OCT_Glau_Norm_p_basis_var[p] <- FVE.samp.var.p(X_ij_p_temp, NRR_OCT_Glau_Norm[,-c(1:8)], method="Kron")
#   NRR_OCT_Glau_Norm_FVE[p] <- (NRR_OCT_Glau_Norm_TV-NRR_OCT_Glau_Norm_p_basis_var[p])/NRR_OCT_Glau_Norm_TV
# }

NRR_OCT_Glau_Norm_p_basis_var_path <- paste(filePath,"/rdata_files/NRR_OCT_Glau_Norm_p_basis_var.RData",sep="")
# save(NRR_OCT_Glau_Norm_p_basis_var, file = NRR_OCT_Glau_Norm_p_basis_var_path)
load(NRR_OCT_Glau_Norm_p_basis_var_path)

NRR_OCT_Glau_Norm_FVE_path <- paste(filePath,"/rdata_files/NRR_OCT_Glau_Norm_FVE.RData",sep="")
# save(NRR_OCT_Glau_Norm_FVE, file = NRR_OCT_Glau_Norm_FVE_path)
load(NRR_OCT_Glau_Norm_FVE_path)

# for(p in 1:p_seq_len){
#   print(p)
#   X_ij_p_temp = t(eval.fd(d_NRR, fdNRR_Fundus_Glau_Norm_p[[p]]))
#   NRR_Fundus_Glau_Norm_p_basis_var[p] <- FVE.samp.var.p(X_ij_p_temp, NRR_Fundus_Glau_Norm[,-c(1:8)], method="Kron")
#   NRR_Fundus_Glau_Norm_FVE[p] <- (NRR_Fundus_Glau_Norm_TV-NRR_Fundus_Glau_Norm_p_basis_var[p])/NRR_Fundus_Glau_Norm_TV
# }

NRR_Fundus_Glau_Norm_p_basis_var_path <- paste(filePath,"/rdata_files/NRR_Fundus_Glau_Norm_p_basis_var.RData",sep="")
# save(NRR_Fundus_Glau_Norm_p_basis_var, file = NRR_Fundus_Glau_Norm_p_basis_var_path)
load(NRR_Fundus_Glau_Norm_p_basis_var_path)

NRR_Fundus_Glau_Norm_FVE_path <- paste(filePath,"/rdata_files/NRR_Fundus_Glau_Norm_FVE.RData",sep="")
# save(NRR_Fundus_Glau_Norm_FVE, file = NRR_Fundus_Glau_Norm_FVE_path)
load(NRR_Fundus_Glau_Norm_FVE_path)

# Determine p
fve_colors <- c("dodgerBlue3", "firebrick1", "chartreuse", "cornflowerblue", "darkgreen", "darkorchid", "darkred", "goldenrod1")
p_seq <- seq(from=3, to=21, by=2)

par(mfrow = c(2,2))
plot(p_seq, NRR_OCT_Normal_OG_FVE, type="l", pch="*", main="NRR OCT Normal OG FVE", 
     xlab="Nmb of Basis Func.", ylab="FVE", col=fve_colors[1])
abline(v=7, col="red")
plot(p_seq, NRR_Fundus_Normal_OG_FVE, type="l", pch="*", main="NRR Fundus Normal OG FVE", 
     xlab="Nmb of Basis Func.", ylab="FVE", col=fve_colors[3])
abline(v=7, col="red")
plot(p_seq, NRR_OCT_Glau_OG_FVE, type="l", pch="*", main="NRR OCT Glaucoma OG FVE", 
     xlab="Nmb of Basis Func.", ylab="FVE", col=fve_colors[5])
abline(v=7, col="red")
plot(p_seq, NRR_Fundus_Glau_OG_FVE, type="l", pch="*", main="NRR Fundus Glaucoma OG FVE", 
     xlab="Nmb of Basis Func.", ylab="FVE", col=fve_colors[7])
abline(v=7, col="red")

par(mfrow = c(2,2))
plot(p_seq, NRR_OCT_Normal_Norm_FVE, type="l", pch="*", main="NRR OCT Normal Normalized FVE", 
     xlab="Nmb of Basis Func.", ylab="FVE", col=fve_colors[1])
abline(v=15, col="red")
plot(p_seq, NRR_Fundus_Normal_Norm_FVE, type="l", pch="*", main="NRR Fundus Normal Normalized FVE", 
     xlab="Nmb of Basis Func.", ylab="FVE", col=fve_colors[3])
abline(v=15, col="red")
plot(p_seq, NRR_OCT_Glau_Norm_FVE, type="l", pch="*", main="NRR OCT Glaucoma Normalized FVE", 
     xlab="Nmb of Basis Func.", ylab="FVE", col=fve_colors[6])
abline(v=15, col="red")
plot(p_seq, NRR_Fundus_Glau_Norm_FVE, type="l", pch="*", main="NRR Fundus Glaucoma Normalized FVE", 
     xlab="Nmb of Basis Func.", ylab="FVE", col=fve_colors[7])
abline(v=15, col="red")
par(mfrow=c(1,1))

###################################################################
###################################################################

# names(NRR_OCT_Normal_OG_FVE) <- c("p=3","p=5","p=7","p=9","p=11","p=13","p=15","p=17","p=19","p=21")
# names(NRR_Fundus_Normal_OG_FVE) <- c("p=3","p=5","p=7","p=9","p=11","p=13","p=15","p=17","p=19","p=21")
# names(NRR_OCT_Glau_OG_FVE) <- c("p=3","p=5","p=7","p=9","p=11","p=13","p=15","p=17","p=19","p=21")
# names(NRR_Fundus_Glau_OG_FVE) <- c("p=3","p=5","p=7","p=9","p=11","p=13","p=15","p=17","p=19","p=21")
# 
# NRR_OG_FVE_list <- list(NRR_OCT_Normal_OG_FVE,NRR_Fundus_Normal_OG_FVE,NRR_OCT_Glau_OG_FVE,NRR_Fundus_Glau_OG_FVE)
nrr_OG_fve_path <- paste(filePath, "/rdata_files/NRR_OG_FVE_list.RData", sep = "")
# save(NRR_OG_FVE_list, file=nrr_OG_fve_path)
load(nrr_OG_fve_path)


nrr_OG_all_schemes_names <- c("NRR_OCT_Normal_OG", "NRR_Fundus_Normal_OG", 
                              "NRR_OCT_Glau_OG", "NRR_Fundus_Glau_OG")

nrr_OG_all_fve <- cbind(NRR_OCT_Normal_OG_FVE, NRR_Fundus_Normal_OG_FVE, 
                        NRR_OCT_Glau_OG_FVE, NRR_Fundus_Glau_OG_FVE)
nrr_OG_all_fve_path <- paste(filePath, "/rdata_files/nrr_OG_all_fve.RData", sep = "")
# save(nrr_OG_all_fve, file = nrr_OG_all_fve_path)
load(nrr_OG_all_fve_path)

NRR_OCT_Normal_OG_FVE <- NRR_OG_FVE_list[[1]]
NRR_Fundus_Normal_OG_FVE <- NRR_OG_FVE_list[[2]]
NRR_OCT_Glau_OG_FVE <- NRR_OG_FVE_list[[3]]
NRR_Fundus_Glau_OG_FVE <- NRR_OG_FVE_list[[4]]

# Display the FVE values for each age group in the console.
print(nrr_OG_all_fve)

#           NRR_OCT_Normal_OG_FVE   NRR_Fundus_Normal_OG_FVE  NRR_OCT_Glau_OG_FVE     NRR_Fundus_Glau_OG_FVE
# p=3              0.8385771                0.7834860           0.7737223              0.8160698
# p=5              0.9440606                0.9388255           0.9167159              0.9324435
# p=7              0.9729094                0.9596066           0.9619225              0.9668737
# p=9              0.9855419                0.9688663           0.9771509              0.9785326
# p=11             0.9915207                0.9739556           0.9874621              0.9848932
# p=13             0.9952556                0.9769731           0.9927810              0.9886571
# p=15             0.9973323                0.9787198           0.9958425              0.9912233
# p=17             0.9984311                0.9803543           0.9975320              0.9928758
# p=19             0.9990648                0.9815490           0.9985622              0.9942325
# p=21             0.9994191                0.9823769           0.9991777              0.9950911
#
#  p = 7 gets to 95% FVE for all of the schemes

###################################################################
###################################################################

# NRR_OCT_Normal_Norm_FVE
# NRR_Fundus_Normal_Norm_FVE
# NRR_OCT_Glau_Norm_FVE
# NRR_Fundus_Glau_Norm_FVE
# 
# names(NRR_OCT_Normal_Norm_FVE) <- c("p=3","p=5","p=7","p=9","p=11","p=13","p=15","p=17","p=19","p=21")
# names(NRR_Fundus_Normal_Norm_FVE) <- c("p=3","p=5","p=7","p=9","p=11","p=13","p=15","p=17","p=19","p=21")
# names(NRR_OCT_Glau_Norm_FVE) <- c("p=3","p=5","p=7","p=9","p=11","p=13","p=15","p=17","p=19","p=21")
# names(NRR_Fundus_Glau_Norm_FVE) <- c("p=3","p=5","p=7","p=9","p=11","p=13","p=15","p=17","p=19","p=21")
# 
# NRR_FVE_list <- list(NRR_OCT_Normal_Norm_FVE,NRR_Fundus_Normal_Norm_FVE,NRR_OCT_Glau_Norm_FVE,NRR_Fundus_Glau_Norm_FVE)
nrr_fve_path <- paste(filePath, "/rdata_files/NRR_FVE_list.RData", sep = "")
# save(NRR_FVE_list, file=nrr_fve_path)
load(nrr_fve_path)

NRR_OCT_Normal_Norm_FVE <- NRR_FVE_list[[1]]
NRR_Fundus_Normal_Norm_FVE <- NRR_FVE_list[[2]]
NRR_OCT_Glau_Norm_FVE <- NRR_FVE_list[[3]]
NRR_Fundus_Glau_Norm_FVE <- NRR_FVE_list[[4]]


###################################################################
#####  Give all the FVE values for the different schemes and  #####
#####        number of basis functions                        #####
###################################################################
nrr_all_schemes_names <- c("NRR_OCT_Normal_Norm", "NRR_Fundus_Normal_Norm", 
                           "NRR_OCT_Glau_Norm", "NRR_Fundus_Glau_Norm")

# nrr_all_fve <- cbind(NRR_OCT_Normal_Norm_FVE, NRR_Fundus_Normal_Norm_FVE, 
#                      NRR_OCT_Glau_Norm_FVE, NRR_Fundus_Glau_Norm_FVE)
nrr_all_fve_path <- paste(filePath, "/rdata_files/nrr_all_fve.RData", sep = "")
# save(nrr_all_fve, file = nrr_all_fve_path)
load(nrr_all_fve_path)

# Record the FVE values for each age group in a table (pdf file)
# this function doesn't output the entire table, so something 
# should be corrected if we want to use it.
###################################################################
# kable_as_image(kable(nrr_all_fve, col.names = nrr_all_schemes_names,
#                      caption = "FVE Values for OCT and Fundus NRR Normalized Curves",
#                      format = "latex", booktabs = F),
#                filename = paste(filePath,"/NRR_FVE_Table",sep=""), file_format = "pdf" )
###################################################################

# Display the FVE values for each age group in the console.
print(nrr_all_fve)

###################################################################
#####            FVE Output for Normalized Curves             #####
###################################################################

# NRR_OCT_Normal_Norm_FVE NRR_Fundus_Normal_Norm_FVE NRR_OCT_Glau_Norm_FVE NRR_Fundus_Glau_Norm_FVE
# p=3                0.2936284                  0.3327743             0.2816982                0.5280061
# p=5                0.7663243                  0.8257836             0.7085072                0.8070289
# p=7                0.8879440                  0.8937945             0.8574586                0.9027854
# p=9                0.9399380                  0.9241096             0.9164391                0.9356802
# p=11               0.9647831                  0.9390393             0.9516135                0.9539397
# p=13               0.9808578                  0.9480361             0.9713287                0.9660700
# p=15               0.9891840                  0.9529929             0.9829383                0.9735362
# p=17               0.9937817                  0.9574773             0.9901461                0.9785490
# p=19               0.9963349                  0.9605712             0.9942151                0.9823238
# p=21               0.9977617                  0.9629177             0.9966617                0.9847068

# p = 15 achieves over 95% fve for all of the schemes

#######################################################
# Write the functional representation to spreadsheets #
#######################################################
fn_rep_coefs_NRR_OCT_Normal_Norm_p15 <- data.frame(NRR_OCT_Normal_Norm[,1:8], t(fdNRR_OCT_Normal_Norm_p15[[1]]))
names(fn_rep_coefs_NRR_OCT_Normal_Norm_p15)[1:8] <- names(fullOCTdata_normal_NRR)[1:8]

fn_rep_coefs_NRR_Fundus_Normal_Norm_p15 <- data.frame(NRR_Fundus_Normal_Norm[,1:8], t(fdNRR_Fundus_Normal_Norm_p15[[1]]))
names(fn_rep_coefs_NRR_Fundus_Normal_Norm_p15)[1:8] <- names(fullOCTdata_normal_NRR)[1:8]

fn_rep_coefs_NRR_OCT_Glau_Norm_p15 <- data.frame(NRR_OCT_Glau_Norm[,1:8], t(fdNRR_OCT_Glau_Norm_p15[[1]]))
names(fn_rep_coefs_NRR_OCT_Glau_Norm_p15)[1:8] <- names(fullOCTdata_glaucoma_NRR)[1:8]

fn_rep_coefs_NRR_Fundus_Glau_Norm_p15 <- data.frame(NRR_Fundus_Glau_Norm[,1:8], t(fdNRR_Fundus_Glau_Norm_p15[[1]]))
names(fn_rep_coefs_NRR_Fundus_Glau_Norm_p15)[1:8] <- names(fullOCTdata_glaucoma_NRR)[1:8]

fn_rep_coefs_list <- list(fn_rep_coefs_NRR_OCT_Normal_Norm_p15,
                          fn_rep_coefs_NRR_Fundus_Normal_Norm_p15,
                          fn_rep_coefs_NRR_OCT_Glau_Norm_p15,
                          fn_rep_coefs_NRR_Fundus_Glau_Norm_p15)

# write_xlsx(fn_rep_coefs_list, path = paste(localPath_tables,"Coefs_p15.xlsx"))

###################################################################
# Evaluated Functional Representation; also to be used in our 
# elastic distance calculation (15 basis functions)
###################################################################

NRR_OCT_Normal_Norm_p15_evaluated <- eval.fd(d_NRR, fdNRR_OCT_Normal_Norm_p15)  # dim = 180 668
NRR_Fundus_Normal_Norm_p15_evaluated <- eval.fd(d_NRR, fdNRR_Fundus_Normal_Norm_p15)  # dim = 180 668

NRR_OCT_Glau_Norm_p15_evaluated <- eval.fd(d_NRR, fdNRR_OCT_Glau_Norm_p15)  # dim = 180 104
NRR_Fundus_Glau_Norm_p15_evaluated <- eval.fd(d_NRR, fdNRR_Fundus_Glau_Norm_p15)  # dim = 180 104

NRR_OCT_Normal_Norm_p15_eval_full <- data.frame(NRR_OCT_Normal_Norm[,1:8], t(NRR_OCT_Normal_Norm_p15_evaluated))
names(NRR_OCT_Normal_Norm_p15_eval_full)[1:8] <- names(fullOCTdata_normal_NRR)[1:8]

NRR_Fundus_Normal_Norm_p15_eval_full <- data.frame(NRR_Fundus_Normal_Norm[,1:8], t(NRR_Fundus_Normal_Norm_p15_evaluated))
names(NRR_Fundus_Normal_Norm_p15_eval_full)[1:8] <- names(fullOCTdata_normal_NRR)[1:8]

NRR_OCT_Glau_Norm_p15_eval_full <- data.frame(NRR_OCT_Glau_Norm[,1:8], t(NRR_OCT_Glau_Norm_p15_evaluated))
names(NRR_OCT_Glau_Norm_p15_eval_full)[1:8] <- names(fullOCTdata_glaucoma_NRR)[1:8]

NRR_Fundus_Glau_Norm_p15_eval_full <- data.frame(NRR_Fundus_Glau_Norm[,1:8], t(NRR_Fundus_Glau_Norm_p15_evaluated))
names(NRR_Fundus_Glau_Norm_p15_eval_full)[1:8] <- names(fullOCTdata_glaucoma_NRR)[1:8]

fn_rep_eval_list <- list(NRR_OCT_Normal_Norm_p15_eval_full,
                         NRR_Fundus_Normal_Norm_p15_eval_full,
                         NRR_OCT_Glau_Norm_p15_eval_full,
                         NRR_Fundus_Glau_Norm_p15_eval_full)

# write_xlsx(fn_rep_eval_list, path = paste(localPath_tables,"Evaluated_p15.xlsx"))


###################################################################
###################################################################
#####    Evaluate the p=7 values for non-normalized curves    #####
###################################################################
#
# Evaluated Functional Representation; also to be used in our 
# elastic distance calculation (7 basis functions)
#
NRR_OCT_Normal_OG_p07_evaluated <- eval.fd(d_NRR, fdNRR_OCT_Normal_OG_p07)  # dim = 180 668
NRR_Fundus_Normal_OG_p07_evaluated <- eval.fd(d_NRR, fdNRR_Fundus_Normal_OG_p07)  # dim = 180 668

NRR_OCT_Glau_OG_p07_evaluated <- eval.fd(d_NRR, fdNRR_OCT_Glau_OG_p07)  # dim = 180 104
NRR_Fundus_Glau_OG_p07_evaluated <- eval.fd(d_NRR, fdNRR_Fundus_Glau_OG_p07)  # dim = 180 104

NRR_OCT_Normal_OG_p07_eval_full <- data.frame(NRR_OCT_Normal_OG[,1:8], t(NRR_OCT_Normal_OG_p07_evaluated))
names(NRR_OCT_Normal_OG_p07_eval_full)[1:8] <- names(fullOCTdata_normal_NRR)[1:8]

NRR_Fundus_Normal_OG_p07_eval_full <- data.frame(NRR_Fundus_Normal_OG[,1:8], t(NRR_Fundus_Normal_OG_p07_evaluated))
names(NRR_Fundus_Normal_OG_p07_eval_full)[1:8] <- names(fullOCTdata_normal_NRR)[1:8]

NRR_OCT_Glau_OG_p07_eval_full <- data.frame(NRR_OCT_Glau_OG[,1:8], t(NRR_OCT_Glau_OG_p07_evaluated))
names(NRR_OCT_Glau_OG_p07_eval_full)[1:8] <- names(fullOCTdata_glaucoma_NRR)[1:8]

NRR_Fundus_Glau_OG_p07_eval_full <- data.frame(NRR_Fundus_Glau_OG[,1:8], t(NRR_Fundus_Glau_OG_p07_evaluated))
names(NRR_Fundus_Glau_OG_p07_eval_full)[1:8] <- names(fullOCTdata_glaucoma_NRR)[1:8]

fn_rep_OG_eval_list <- list(NRR_OCT_Normal_OG_p07_eval_full,
                            NRR_Fundus_Normal_OG_p07_eval_full,
                            NRR_OCT_Glau_OG_p07_eval_full,
                            NRR_Fundus_Glau_OG_p07_eval_full)

# write_xlsx(fn_rep_OG_eval_list, path = paste(localPath_tables,"Evaluated_OG_p07.xlsx"))

###################################################################
###########             Elastic Distances            ##############
###################################################################
# Functions to be used in our elastic distance calculation of 
# normalized curves (15 basis functions)
#
# NRR_OCT_Normal_Norm_p15_evaluated
# NRR_Fundus_Normal_Norm_p15_evaluated
# NRR_OCT_Glau_Norm_p15_evaluated
# NRR_Fundus_Glau_Norm_p15_evaluated

####################################
#####  Normal Normalized eyes  #####
####################################
#  elastic.distance returns a list, so we create a list of lists to store all of the distances
NRR_Normal_elast_dist <- rep(list(list()), ncol(NRR_OCT_Normal_Norm_p15_evaluated))

for(c in 1:length(NRR_Normal_elast_dist)){
  temp_elast_dist <- elastic.distance(
    f1 = NRR_OCT_Normal_Norm_p15_evaluated[,c],
    f2 = NRR_Fundus_Normal_Norm_p15_evaluated[,c],
    time = d_NRR,
    pen = "geodesic"
  )
  NRR_Normal_elast_dist[[c]] <- temp_elast_dist
}

NRR_Normal_Norm_Dx <- rep(NA, length(NRR_Normal_elast_dist))
NRR_Normal_Norm_Dy <- rep(NA, length(NRR_Normal_elast_dist))

for(i in 1:length(NRR_Normal_elast_dist)){
  NRR_Normal_Norm_Dx[i] <- NRR_Normal_elast_dist[[i]]$Dx
  NRR_Normal_Norm_Dy[i] <- NRR_Normal_elast_dist[[i]]$Dy
}

# Dx: phase distance
# Dy: amplitude distance
scatter_Normal <- qplot(NRR_Normal_Norm_Dx,NRR_Normal_Norm_Dy,
                        xlab = "Phase Distance", 
                        ylab = "Amplitude Distance", 
                        main = "Elastic Distance btwn OCT and Fundus, Normal Eyes (Normalized NRR Curves)")  + 
  scale_x_continuous(limits=c(min(NRR_Normal_Norm_Dx),max(NRR_Normal_Norm_Dx))) + 
  scale_y_continuous(limits=c(min(NRR_Normal_Norm_Dy),max(NRR_Normal_Norm_Dy))) + 
  geom_rug(col=rgb(.5,0,0,alpha=.2))
scatter_Normal




######################################
#####  Glaucoma Normalized eyes  #####
######################################
NRR_Glau_elast_dist <- rep(list(list()), ncol(NRR_OCT_Glau_Norm_p15_evaluated))

for(c in 1:length(NRR_Glau_elast_dist)){
  temp_elast_dist <- elastic.distance(
    f1 = NRR_OCT_Glau_Norm_p15_evaluated[,c],
    f2 = NRR_Fundus_Glau_Norm_p15_evaluated[,c],
    time = d_NRR,
    pen = "geodesic"
  )
  NRR_Glau_elast_dist[[c]] <- temp_elast_dist
}

NRR_Glau_Norm_Dx <- rep(NA, length(NRR_Glau_elast_dist))
NRR_Glau_Norm_Dy <- rep(NA, length(NRR_Glau_elast_dist))

for(i in 1:length(NRR_Glau_elast_dist)){
  NRR_Glau_Norm_Dx[i] <- NRR_Glau_elast_dist[[i]]$Dx
  NRR_Glau_Norm_Dy[i] <- NRR_Glau_elast_dist[[i]]$Dy
}


scatter_Glau <- qplot(NRR_Glau_Norm_Dx,NRR_Glau_Norm_Dy,
                      xlab = "Phase Distance", 
                      ylab = "Amplitude Distance", 
                      main = "Elastic Distance btwn OCT and Fundus, Glaucoma Eyes (Normalized NRR Curves)")  + 
  scale_x_continuous(limits=c(min(NRR_Normal_Norm_Dx),max(NRR_Normal_Norm_Dx))) + 
  scale_y_continuous(limits=c(min(NRR_Normal_Norm_Dy),max(NRR_Normal_Norm_Dy))) + 
  geom_rug(col=rgb(.5,0,0,alpha=.2))
scatter_Glau

########################################
#####  Normal Original Curve eyes  #####
########################################
# Functions to be used in our elastic distance calculation of 
# non-normalized (original) curves (7 basis functions)
#
# NRR_OCT_Normal_OG_p07_evaluated
# NRR_Fundus_Normal_OG_p07_evaluated
# NRR_OCT_Glau_OG_p07_evaluated
# NRR_Fundus_Glau_OG_p07_evaluated

#  elastic.distance returns a list, so we create a list of lists to store all of the distances
NRR_Normal_OG_elast_dist <- rep(list(list()), ncol(NRR_OCT_Normal_OG_p07_evaluated))

for(c in 1:length(NRR_Normal_OG_elast_dist)){
  temp_elast_dist <- elastic.distance(
    f1 = NRR_OCT_Normal_OG_p07_evaluated[,c],
    f2 = NRR_Fundus_Normal_OG_p07_evaluated[,c],
    time = d_NRR,
    pen = "geodesic"
  )
  NRR_Normal_OG_elast_dist[[c]] <- temp_elast_dist
}

NRR_Normal_OG_Dx <- rep(NA, length(NRR_Normal_OG_elast_dist))
NRR_Normal_OG_Dy <- rep(NA, length(NRR_Normal_OG_elast_dist))

for(i in 1:length(NRR_Normal_OG_elast_dist)){
  NRR_Normal_OG_Dx[i] <- NRR_Normal_OG_elast_dist[[i]]$Dx
  NRR_Normal_OG_Dy[i] <- NRR_Normal_OG_elast_dist[[i]]$Dy
}

# Dx: phase distance
# Dy: amplitude distance
scatter_Normal_OG <- qplot(NRR_Normal_OG_Dx,NRR_Normal_OG_Dy,
                           xlab = "Phase Distance", 
                           ylab = "Amplitude Distance", 
                           main = "Elastic Distance btwn OCT and Fundus, Normal Eyes (Original NRR Curves)")  + 
  scale_x_continuous(limits=c(min(NRR_Normal_OG_Dx),max(NRR_Normal_OG_Dx))) + 
  scale_y_continuous(limits=c(min(NRR_Normal_OG_Dy),max(NRR_Normal_OG_Dy))) + 
  geom_rug(col=rgb(.5,0,0,alpha=.2))
scatter_Normal_OG

NRR_Normal_OG_Elastic_Amp_Phase_Dist <- data.frame(NRR_Normal_OG_Dx=NRR_Normal_OG_Dx, 
                                                   NRR_Normal_OG_Dy=NRR_Normal_OG_Dy)

# write_xlsx(NRR_Normal_OG_Elastic_Amp_Phase_Dist, 
#            path = paste(localPath_tables,"NRR_Normal_OG_Elastic_Amp_Phase_Dist.xlsx"))

##########################################
#####  Glaucoma Original Curve eyes  #####
##########################################
NRR_Glau_OG_elast_dist <- rep(list(list()), ncol(NRR_OCT_Glau_OG_p07_evaluated))

for(c in 1:length(NRR_Glau_OG_elast_dist)){
  temp_elast_dist <- elastic.distance(
    f1 = NRR_OCT_Glau_OG_p07_evaluated[,c],
    f2 = NRR_Fundus_Glau_OG_p07_evaluated[,c],
    time = d_NRR,
    pen = "geodesic"
  )
  NRR_Glau_OG_elast_dist[[c]] <- temp_elast_dist
}

NRR_Glau_OG_Dx <- rep(NA, length(NRR_Glau_OG_elast_dist))
NRR_Glau_OG_Dy <- rep(NA, length(NRR_Glau_OG_elast_dist))

for(i in 1:length(NRR_Glau_OG_elast_dist)){
  NRR_Glau_OG_Dx[i] <- NRR_Glau_OG_elast_dist[[i]]$Dx
  NRR_Glau_OG_Dy[i] <- NRR_Glau_OG_elast_dist[[i]]$Dy
}

scatter_Glau_OG <- qplot(NRR_Glau_OG_Dx,NRR_Glau_OG_Dy,
                         xlab = "Phase Distance", 
                         ylab = "Amplitude Distance", 
                         main = "Elastic Distance btwn OCT and Fundus, Glaucoma Eyes (Original NRR Curves)")  + 
  scale_x_continuous(limits=c(min(NRR_Normal_OG_Dx),max(NRR_Normal_OG_Dx))) + 
  scale_y_continuous(limits=c(min(NRR_Normal_OG_Dy),max(NRR_Normal_OG_Dy))) + 
  geom_rug(col=rgb(.5,0,0,alpha=.2))
scatter_Glau_OG

NRR_Glau_OG_Elastic_Amp_Phase_Dist <- data.frame(NRR_Glau_OG_Dx=NRR_Glau_OG_Dx, 
                                                   NRR_Glau_OG_Dy=NRR_Glau_OG_Dy)

# write_xlsx(NRR_Glau_OG_Elastic_Amp_Phase_Dist, 
#            path = paste(localPath_tables,"NRR_Glau_OG_Elastic_Amp_Phase_Dist.xlsx"))

###################################################################
###################################################################
###################################################################

###########################################
#####          Time Warping           #####
###########################################
# An object of class fdawarp which is a list with the following components:
#   
# time: a numeric vector of length M storing the original grid;
# f0: a numeric matrix of shape M \times N storing the original sample of N functions observed on a grid of size M;
# q0: a numeric matrix of the same shape as f0 storing the original SRSFs;
# fn: a numeric matrix of the same shape as f0 storing the aligned functions;
# qn: a numeric matrix of the same shape as f0 storing the aligned SRSFs;
# fmean: a numeric vector of length M storing the mean or median curve;
# mqn: a numeric vector of length M storing the mean or median SRSF;
# warping_functions: a numeric matrix of the same shape as f0 storing the estimated warping functions;
# original_variance: a numeric value storing the variance of the original sample;
# amplitude_variance: a numeric value storing the variance in amplitude of the aligned sample;
# phase_variance: a numeric value storing the variance in phase of the aligned sample;
# qun: a numeric vector of maximum length max_iter + 2 storing the values of the cost function after each iteration;
# lambda: the input parameter lambda which specifies the elasticity;
# centroid_type: the input centroid type;
# optim_method: the input optimization method;
# inverse_average_warping_function: the inverse of the mean estimated warping function;
# rsamps: TO DO.
# 
#
################################################################### 
# Functions to be used in our time warping aligned function 'fmean'
# calculation (15 basis functions)
########################################
##### Using the normalized curves  #####
########################################
# NRR_OCT_Normal_Norm_p15_evaluated -- dim = 180 X 668
# NRR_Fundus_Normal_Norm_p15_evaluated -- dim = 180 X 668
# NRR_OCT_Glau_Norm_p15_evaluated -- dim = 180 104
# NRR_Fundus_Glau_Norm_p15_evaluated -- dim = 180 104
#
# NRR_OCT_Normal_OG_p07_evaluated -- dim = 180 X 668
# NRR_Fundus_Normal_OG_p07_evaluated -- dim = 180 X 668
# NRR_OCT_Glau_OG_p07_evaluated -- dim = 180 104
# NRR_Fundus_Glau_OG_p07_evaluated -- dim = 180 104

 # elastic.distance returns a list, so we create a list of lists to store all of the distances
# NRR_Normal_twarp_dist <- rep(list(list()), ncol(NRR_OCT_Normal_Norm_p15_evaluated))
# 
# for(i in 1:ncol(NRR_OCT_Normal_Norm_p15_evaluated)){
#   temp_twarp_NRR_Norm_geod <- time_warping(cbind(NRR_OCT_Normal_Norm_p15_evaluated[,i],
#                                                  NRR_Fundus_Normal_Norm_p15_evaluated[,i]),
#                                            d_NRR,
#                                            penalty_method = "geodesic")
#   NRR_Normal_twarp_dist[[i]] <- temp_twarp_NRR_Norm_geod
# }
# 
NRR_Normal_twarp_dist_path <- paste(filePath,"/rdata_files/NRR_Normal_twarp_dist.RData",sep="")
# save(NRR_Normal_twarp_dist, file = NRR_Normal_twarp_dist_path)
load(NRR_Normal_twarp_dist_path)

# keeping the dimensions consistent here with rows = observations, cols = curves
NRR_Normal_twarp_dist_fmean <- matrix(data = NA, 
                                      nrow = 180, 
                                      ncol = ncol(NRR_OCT_Normal_Norm_p15_evaluated),
                                      byrow = F)

for(i in 1:length(NRR_Normal_twarp_dist)){
  NRR_Normal_twarp_dist_fmean[,i] <- NRR_Normal_twarp_dist[[i]]$fmean
}

NRR_Normal_Norm_twarp_dist_fmean_p15_full <- data.frame(NRR_OCT_Normal_Norm[,1:8], t(NRR_Normal_twarp_dist_fmean))
names(NRR_Normal_Norm_twarp_dist_fmean_p15_full)[1:8] <- names(fullOCTdata_normal_NRR)[1:8]

# write_xlsx(NRR_Normal_Norm_twarp_dist_fmean_p15_full, path = paste(localPath_tables,"NRR_Normal_Norm_timewarp_aligned_fmean.xlsx"))

###################################################################
##############   For the Glaucoma Eyes now   ######################
###################################################################
# 
#  elastic.distance returns a list, so we create a list of lists to store all of the distances
# NRR_Glau_twarp_dist <- rep(list(list()), ncol(NRR_OCT_Glau_Norm_p15_evaluated))
# 
# for(i in 1:ncol(NRR_OCT_Glau_Norm_p15_evaluated)){
#   temp_twarp_NRR_Norm_geod <- time_warping(cbind(NRR_OCT_Glau_Norm_p15_evaluated[,i],
#                                                  NRR_Fundus_Glau_Norm_p15_evaluated[,i]),
#                                            d_NRR,
#                                            penalty_method = "geodesic")
#   NRR_Glau_twarp_dist[[i]] <- temp_twarp_NRR_Norm_geod
# }

NRR_Glau_twarp_dist_path <- paste(filePath,"/rdata_files/NRR_Glau_twarp_dist.RData",sep="")
# save(NRR_Glau_twarp_dist, file = NRR_Glau_twarp_dist_path)
load(NRR_Glau_twarp_dist_path)

# keeping the dimensions consistent here with rows = observations, cols = curves
NRR_Glau_twarp_dist_fmean <- matrix(data = NA, 
                                      nrow = 180, 
                                      ncol = ncol(NRR_OCT_Glau_Norm_p15_evaluated),
                                      byrow = F)

for(i in 1:length(NRR_Glau_twarp_dist)){
  NRR_Glau_twarp_dist_fmean[,i] <- NRR_Glau_twarp_dist[[i]]$fmean
}

NRR_Glau_twarp_dist_fmean_p15_full <- data.frame(NRR_OCT_Glau_Norm[,1:8], t(NRR_Glau_twarp_dist_fmean))
names(NRR_Glau_twarp_dist_fmean_p15_full)[1:8] <- names(fullOCTdata_normal_NRR)[1:8]

# write_xlsx(NRR_Glau_twarp_dist_fmean_p15_full, path = paste(localPath_tables,"NRR_Glau_Norm_timewarp_aligned_fmean.xlsx"))


###################################################################
###################################################################
#####     For non-normalized curves
#
# NRR_OCT_Normal_OG_p07_evaluated -- dim = 180 X 668
# NRR_Fundus_Normal_OG_p07_evaluated -- dim = 180 X 668
# NRR_OCT_Glau_OG_p07_evaluated -- dim = 180 104
# NRR_Fundus_Glau_OG_p07_evaluated -- dim = 180 104

#  elastic.distance returns a list, so we create a list of lists to store all of the distances
# NRR_Normal_OG_twarp_dist <- rep(list(list()), ncol(NRR_OCT_Normal_OG_p07_evaluated))
# 
# for(i in 1:ncol(NRR_OCT_Normal_OG_p07_evaluated)){
#   temp_twarp_NRR_OG_geod <- time_warping(cbind(NRR_OCT_Normal_OG_p07_evaluated[,i],
#                                                NRR_Fundus_Normal_OG_p07_evaluated[,i]),
#                                          d_NRR,
#                                          penalty_method = "geodesic")
#   NRR_Normal_OG_twarp_dist[[i]] <- temp_twarp_NRR_OG_geod
# }

NRR_Normal_OG_twarp_dist_path <- paste(filePath,"/rdata_files/NRR_Normal_OG_twarp_dist.RData",sep="")
# save(NRR_Normal_OG_twarp_dist, file = NRR_Normal_OG_twarp_dist_path)
load(NRR_Normal_OG_twarp_dist_path)

# keeping the dimensions consistent here with rows = observations, cols = curves
NRR_Normal_OG_twarp_dist_fmean <- matrix(data = NA, 
                                         nrow = 180, 
                                         ncol = ncol(NRR_OCT_Normal_OG_p07_evaluated),
                                         byrow = F)

for(i in 1:length(NRR_Normal_OG_twarp_dist)){
  NRR_Normal_OG_twarp_dist_fmean[,i] <- NRR_Normal_OG_twarp_dist[[i]]$fmean
}

NRR_Normal_OG_twarp_dist_fmean_p07_full <- data.frame(NRR_OCT_Normal_OG[,1:8], t(NRR_Normal_OG_twarp_dist_fmean))
names(NRR_Normal_OG_twarp_dist_fmean_p07_full)[1:8] <- names(fullOCTdata_normal_NRR)[1:8]

# write_xlsx(NRR_Normal_OG_twarp_dist_fmean_p07_full, path = paste(localPath_tables,"NRR_Normal_OG_timewarp_aligned_fmean.xlsx"))



###########  Elastic Distances  ##############

# Example from the package
# distances <- elastic.distance(
#   f1 = simu_data$f[, 1],
#   f2 = simu_data$f[, 2],
#   time = simu_data$time
# )

#  elastic.distance returns a list, so we create a list of lists to store all of the distances
NRR_OG_elast_dist <- rep(list(list()), ncol(NRR_OCT_Normal_OG_p07_evaluated))

for(c in 1:length(NRR_OG_elast_dist)){
  temp_elast_dist <- elastic.distance(
    f1 = NRR_OCT_Normal_OG_p07_evaluated[,c],
    f2 = NRR_Fundus_Normal_OG_p07_evaluated[,c],
    time = d_NRR,
    pen = "geodesic"
  )
  NRR_OG_elast_dist[[c]] <- temp_elast_dist
}



NRR_OG_Dx <- rep(NA, length(NRR_OG_elast_dist))
NRR_OG_Dy <- rep(NA, length(NRR_OG_elast_dist))


for(i in 1:length(NRR_OG_elast_dist)){
  NRR_OG_Dx[i] <- NRR_OG_elast_dist[[i]]$Dx
  NRR_OG_Dy[i] <- NRR_OG_elast_dist[[i]]$Dy
}

NRR_OG_DxDy <- data.frame(PhaseDist=as.numeric(NRR_OG_Dx), Amplitude=as.numeric(NRR_OG_Dy))
# Dx: phase distance
# Dy: amplitude distance

NRR_Normal_OG_p07_ElastDist_DxDy_plot_path <- paste(figurePath,
                                                    "/NRR_Normal_OG_p07_ElastDist_DxDy.pdf",
                                                    sep="")
# pdf(NRR_Normal_OG_p07_ElastDist_DxDy_plot_path, width=6, height=6)
# 
# p <- ggplot(NRR_OG_DxDy, aes(x=PhaseDist, y=Amplitude)) + 
#   geom_point(colour="#000099") + theme_classic() + 
#   geom_vline(xintercept = 0.355) +
#   ggtitle("Elastic Distance OCT/Fundus: Normal Eyes (OG Curves)")
# 
# ggExtra::ggMarginal(p, type = "histogram")
# dev.off()


# Subset for Original Curves eyes - 
# Elastic Distance OCT/Fundus: Normal Eyes (Normalized Curves) phase distance <= 0.355
# NRR_OG_DxDy$PhaseDist

NRR_OG_Fused_p07_DxFiltered <- (NRR_OG_DxDy$PhaseDist <= 0.355)
summary(NRR_OG_DxDy)



##################################
NRR_Norm_elast_dist <- rep(list(list()), ncol(NRR_OCT_Normal_Norm_p15_evaluated))

for(c in 1:length(NRR_Norm_elast_dist)){
  temp_elast_dist <- elastic.distance(
    f1 = NRR_OCT_Normal_Norm_p15_evaluated[,c],
    f2 = NRR_Fundus_Normal_Norm_p15_evaluated[,c],
    time = d_NRR,
    pen = "geodesic"
  )
  NRR_Norm_elast_dist[[c]] <- temp_elast_dist
}



NRR_Norm_Dx <- rep(NA, length(NRR_Norm_elast_dist))
NRR_Norm_Dy <- rep(NA, length(NRR_Norm_elast_dist))


for(i in 1:length(NRR_Norm_elast_dist)){
  NRR_Norm_Dx[i] <- NRR_Norm_elast_dist[[i]]$Dx
  NRR_Norm_Dy[i] <- NRR_Norm_elast_dist[[i]]$Dy
}

NRR_Norm_DxDy <- data.frame(PhaseDist=as.numeric(NRR_Norm_Dx), Amplitude=as.numeric(NRR_Norm_Dy))
# Dx: phase distance
# Dy: amplitude distance

NRR_Normal_Norm_p15_ElastDist_DxDy_plot_path <- paste(figurePath,
                                                    "/NRR_Normal_Norm_p15_ElastDist_DxDy.pdf",
                                                    sep="")
# pdf(NRR_Normal_Norm_p15_ElastDist_DxDy_plot_path, width=6, height=6)
# 
# p <- ggplot(NRR_Norm_DxDy, aes(x=PhaseDist, y=Amplitude)) + 
#   geom_point(colour="#CC0000") + theme_classic() + 
#   geom_vline(xintercept = 0.44) +
#   ggtitle("Elastic Distance OCT/Fundus: Normal Eyes (Normalized Curves)")
# 
# ggExtra::ggMarginal(p, type = "histogram")
# 
# dev.off()

# Subset for Normalized eyes - 
# Elastic Distance OCT/Fundus: Normal Eyes (Normalized Curves) phase distance <= 0.44
# NRR_Norm_DxDy$PhaseDist

NRR_Norm_Fused_p15_DxFiltered <- (NRR_Norm_DxDy$PhaseDist <= 0.44)
View(NRR_Norm_DxDy)
summary(NRR_Norm_DxDy)
write_xlsx(NRR_Norm_DxDy, path = paste(localPath_tables,"NRR_Normal_Norm_Amp_Phase_Dist.xlsx"))
###################################################################
##############   For the Glaucoma Eyes now   ######################
###################################################################
# 
# NRR_Glau_OG_twarp_dist <- rep(list(list()), ncol(NRR_OCT_Glau_OG_p07_evaluated))
# 
# for(i in 1:ncol(NRR_OCT_Glau_OG_p07_evaluated)){
#   temp_twarp_NRR_OG_geod <- time_warping(cbind(NRR_OCT_Glau_OG_p07_evaluated[,i],
#                                                NRR_Fundus_Glau_OG_p07_evaluated[,i]),
#                                          d_NRR,
#                                          penalty_method = "geodesic")
#   NRR_Glau_OG_twarp_dist[[i]] <- temp_twarp_NRR_OG_geod
# }

NRR_Glau_OG_twarp_dist_path <- paste(filePath,"/rdata_files/NRR_Glau_OG_twarp_dist.RData",sep="")
# save(NRR_Glau_OG_twarp_dist, file = NRR_Glau_OG_twarp_dist_path)
load(NRR_Glau_OG_twarp_dist_path)

# keeping the dimensions consistent here with rows = observations, cols = curves
NRR_Glau_OG_twarp_dist_fmean <- matrix(data = NA, 
                                       nrow = 180, 
                                       ncol = ncol(NRR_OCT_Glau_OG_p07_evaluated),
                                       byrow = F)

for(i in 1:length(NRR_Glau_OG_twarp_dist)){
  NRR_Glau_OG_twarp_dist_fmean[,i] <- NRR_Glau_OG_twarp_dist[[i]]$fmean
}

NRR_Glau_OG_twarp_dist_fmean_p07_full <- data.frame(NRR_OCT_Glau_OG[,1:8], t(NRR_Glau_OG_twarp_dist_fmean))
names(NRR_Glau_OG_twarp_dist_fmean_p07_full)[1:8] <- names(fullOCTdata_glaucoma_NRR)[1:8]

# write_xlsx(NRR_Glau_OG_twarp_dist_fmean_p07_full, path = paste(localPath_tables,"NRR_Glau_OG_timewarp_aligned_fmean.xlsx"))
# fdNRR_Fused_Glau_OG_p07 <- Data2fd(argvals=d_NRR, y=t(NRR_Glau_OG_twarp_dist_fmean_p07_full[,-c(1:8)]), basisobj=NRR_basis_07)
fdNRR_Fused_Glau_OG_p07_path <- paste(filePath,"/rdata_files/fdNRR_Fused_Glau_OG_p07.RData",sep="")
# save(fdNRR_Fused_Glau_OG_p07, file = fdNRR_Fused_Glau_OG_p07_path)
load(fdNRR_Fused_Glau_OG_p07_path)



###################################################################
###################################################################



#####  Circular plots
# NRR_OCT_Normal_Norm_p15_evaluated  # dim = 180 668
# NRR_Fundus_Normal_Norm_p15_evaluated  # dim = 180 668
# 
# NRR_OCT_Glau_Norm_p15_evaluated  # dim = 180 104
# NRR_Fundus_Glau_Norm_p15_evaluated  # dim = 180 104

# Potential sample curves to plot
# exmpl_Normal_Norm_curves <- sample(1:668, size=4, replace = F)
# 52 413 581 292

exmpl_dat <- cbind(NRR_OCT_Normal_Norm_p15_evaluated[,52],
                   NRR_Fundus_Normal_Norm_p15_evaluated[,52],
                   NRR_Normal_twarp_dist_fmean[,52])

exmpl_dat <- cbind(NRR_OCT_Normal_Norm_p15_evaluated[,413],
                   NRR_Fundus_Normal_Norm_p15_evaluated[,413],
                   NRR_Normal_twarp_dist_fmean[,413])                   

exmpl_dat <- cbind(NRR_OCT_Normal_Norm_p15_evaluated[,581],
                   NRR_Fundus_Normal_Norm_p15_evaluated[,581],
                   NRR_Normal_twarp_dist_fmean[,581])                   

exmpl_dat <- cbind(NRR_OCT_Normal_Norm_p15_evaluated[,292],
                   NRR_Fundus_Normal_Norm_p15_evaluated[,292],
                   NRR_Normal_twarp_dist_fmean[,292]) 

# plot.curves.on.circle = function(curves, main = NA, col = 1, xy.lim = 1010){
nmb_pts_per_curve <- 180
d = 1:nmb_pts_per_curve

dat_1 <- exmpl_dat[,1]  # OCT curve
dat_2 <- exmpl_dat[,2]  # Fundus curve
dat_3 <- exmpl_dat[,3]  # fmean curve
xc = cos((d/nmb_pts_per_curve)*2*pi)
ys = sin((d/nmb_pts_per_curve)*2*pi)
Xc.mat_1 = t(dat_1) * xc
Ys.mat_1 = t(dat_1) * ys
Xc.mat_2 = t(dat_2) * xc
Ys.mat_2 = t(dat_2) * ys

xy.lim = .01  # max(abs(Xc.mat), abs(Ys.mat)) + 10

# fmean curve
X_bar = t(exmpl_dat[,3])  # need this to be 1 X 180
Xc.bar = X_bar * xc
Ys.bar = X_bar * ys

# plot the average curve last
par(mfrow = c(1,1))
plot(Xc.mat_1, Ys.mat_1, type = "l", col = "red",
     xlim = c(-xy.lim, xy.lim), ylim = c(-xy.lim, xy.lim), 
     xlab = NA, ylab = NA, main="Sample Normal Normalized \n OCT (red) and Fundus (green), and fmean (blue)")


lines(Xc.mat_2, Ys.mat_2, col = "green")  # (i-1)%%10+2)

lines(Xc.bar, Ys.bar, col = "blue", lwd=2)

abline(v = 0, lty = "dashed")
abline(h = 0, lty = "dashed")
abline(0, 1, lty = "dashed")
abline(0, -1, lty = "dashed")

###################################################################
exmpl_Glau_Norm_curves <- sample(1:104, size=4, replace = F)
#  50  7 55 70

exmpl_dat <- cbind(NRR_OCT_Glau_Norm_p15_evaluated[,50],
                   NRR_Fundus_Glau_Norm_p15_evaluated[,50],
                   NRR_Glau_twarp_dist_fmean[,50])

exmpl_dat <- cbind(NRR_OCT_Glau_Norm_p15_evaluated[,7],
                   NRR_Fundus_Glau_Norm_p15_evaluated[,7],
                   NRR_Glau_twarp_dist_fmean[,7])                   

exmpl_dat <- cbind(NRR_OCT_Glau_Norm_p15_evaluated[,55],
                   NRR_Fundus_Glau_Norm_p15_evaluated[,55],
                   NRR_Glau_twarp_dist_fmean[,55]) 

exmpl_dat <- cbind(NRR_OCT_Glau_Norm_p15_evaluated[,70],
                   NRR_Fundus_Glau_Norm_p15_evaluated[,70],
                   NRR_Glau_twarp_dist_fmean[,70]) 


# plot.curves.on.circle = function(curves, main = NA, col = 1, xy.lim = 1010){
nmb_pts_per_curve <- 180
d = 1:nmb_pts_per_curve

dat_1 <- exmpl_dat[,1]  # OCT curve
dat_2 <- exmpl_dat[,2]  # Fundus curve
dat_3 <- exmpl_dat[,3]  # fmean curve
xc = cos((d/nmb_pts_per_curve)*2*pi)
ys = sin((d/nmb_pts_per_curve)*2*pi)
Xc.mat_1 = t(dat_1) * xc
Ys.mat_1 = t(dat_1) * ys
Xc.mat_2 = t(dat_2) * xc
Ys.mat_2 = t(dat_2) * ys

xy.lim = .01  # max(abs(Xc.mat), abs(Ys.mat)) + 10
  
# fmean curve
X_bar = t(exmpl_dat[,3])  # need this to be 1 X 180
Xc.bar = X_bar * xc
Ys.bar = X_bar * ys
  
# plot the average curve last
par(mfrow = c(1,1))
plot(Xc.mat_1, Ys.mat_1, type = "l", col = "red",
     xlim = c(-xy.lim, xy.lim), ylim = c(-xy.lim, xy.lim), 
     xlab = NA, ylab = NA, main="Sample Glaucoma Normalized \n OCT (red) and Fundus (green), and fmean (blue)")


lines(Xc.mat_2, Ys.mat_2, col = "green")
  
lines(Xc.bar, Ys.bar, col = "blue", lwd=2)

abline(v = 0, lty = "dashed")
abline(h = 0, lty = "dashed")
abline(0, 1, lty = "dashed")
abline(0, -1, lty = "dashed")



#################################################################################
#################################################################################
#####                                                                       #####
#####  CIFU Clustering for the NRR Normal Fused Data - OG and Norm          #####
#####                                                                       #####
#################################################################################
#################################################################################

##### Original Curves (OG) #####
# We are using the fused data points from NRR_Normal_OG_twarp_dist_fmean_p07_full
# the first 8 columns are identifier information, followed by 180 points for each
# of 668 observations

View(NRR_Normal_OG_twarp_dist_fmean_p07_full)

# Creating a functional data object, using 7 basis functions (just like the two
# fused curves each used)
# fdNRR_Fused_Normal_OG_p07 <- Data2fd(argvals=d_NRR, 
#                                      y=t(NRR_Normal_OG_twarp_dist_fmean_p07_full[,-c(1:8)]), 
#                                      basisobj=NRR_basis_07)
# 
# NRR_Normal_OG_Fused_p07_k02 = funFEM(fdNRR_Fused_Normal_OG_p07, K=2,  model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_OG_Fused_p07_k03 = funFEM(fdNRR_Fused_Normal_OG_p07, K=3,  model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_OG_Fused_p07_k04 = funFEM(fdNRR_Fused_Normal_OG_p07, K=4,  model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_OG_Fused_p07_k05 = funFEM(fdNRR_Fused_Normal_OG_p07, K=5,  model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_OG_Fused_p07_k06 = funFEM(fdNRR_Fused_Normal_OG_p07, K=6,  model="DkBk", crit='bic', init='kmeans')

# Error in D[k, ((d + 1):p), ((d + 1):p)] <- diag(rep(bk, p - d)) : 
#   number of items to replace is not a multiple of replacement length
# Error in funFEM(fdNRR_Fused_Normal_OG_p07, K = 7, model = "DkBk", crit = "bic",  : 
#                   No reliable results to return (empty clusters in all partitions)!
# Must have fewer clusters than the number of basis functions
# NRR_Normal_OG_Fused_p07_k07 = funFEM(fdNRR_Fused_Normal_OG_p07, K=7,  model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_OG_Fused_p07_k08 = funFEM(fdNRR_Fused_Normal_OG_p07, K=8,  model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_OG_Fused_p07_k09 = funFEM(fdNRR_Fused_Normal_OG_p07, K=9,  model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_OG_Fused_p07_k10 = funFEM(fdNRR_Fused_Normal_OG_p07, K=10, model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_OG_Fused_p07_k11 = funFEM(fdNRR_Fused_Normal_OG_p07, K=11, model="DkBk", crit='bic', init='kmeans')
#####

# Putting all of the output into a list to save as an .RData file
# NRR_Normal_OG_Fused_fem_list <- list(NRR_Normal_OG_Fused_p07_k02,NRR_Normal_OG_Fused_p07_k03,
#                                      NRR_Normal_OG_Fused_p07_k04,NRR_Normal_OG_Fused_p07_k05,
#                                      NRR_Normal_OG_Fused_p07_k06)

NRR_Normal_OG_Fused_fem_list_path <- paste(filePath,"/rdata_files/NRR_Normal_OG_Fused_fem_list.RData",sep = "")

# save(NRR_Normal_OG_Fused_fem_list, file = NRR_Normal_OG_Fused_fem_list_path)

load(NRR_Normal_OG_Fused_fem_list_path)



##### Normalized Curves (Norm) #####
# We are using the fused data points from NRR_Normal_Norm_twarp_dist_fmean_p15_full
# the first 8 columns are identifier information, followed by 180 points for each
# of 668 observations
View(NRR_Normal_Norm_twarp_dist_fmean_p15_full)
# Creating a functional data object, using 15 basis functions (just like the two
# fused curves each used)

# fdNRR_Fused_Normal_Norm_p15 <- Data2fd(argvals=d_NRR,
#                                        y=t(NRR_Normal_Norm_twarp_dist_fmean_p15_full[,-c(1:8)]),
#                                        basisobj=NRR_basis_15)
# 
# NRR_Normal_Norm_Fused_p15_k02 = funFEM(fdNRR_Fused_Normal_Norm_p15, K=2,  model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_Norm_Fused_p15_k03 = funFEM(fdNRR_Fused_Normal_Norm_p15, K=3,  model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_Norm_Fused_p15_k04 = funFEM(fdNRR_Fused_Normal_Norm_p15, K=4,  model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_Norm_Fused_p15_k05 = funFEM(fdNRR_Fused_Normal_Norm_p15, K=5,  model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_Norm_Fused_p15_k06 = funFEM(fdNRR_Fused_Normal_Norm_p15, K=6,  model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_Norm_Fused_p15_k07 = funFEM(fdNRR_Fused_Normal_Norm_p15, K=7,  model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_Norm_Fused_p15_k08 = funFEM(fdNRR_Fused_Normal_Norm_p15, K=8,  model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_Norm_Fused_p15_k09 = funFEM(fdNRR_Fused_Normal_Norm_p15, K=9,  model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_Norm_Fused_p15_k10 = funFEM(fdNRR_Fused_Normal_Norm_p15, K=10, model="DkBk", crit='bic', init='kmeans')
# # Error in funFEM(fdNRR_Fused_Normal_Norm_p15, K = 10, model = "DkBk", crit = "bic",  : 
# #                   No reliable results to return (empty clusters in all partitions)!
# #                   In addition: There were 50 or more warnings (use warnings() to see the first 50)
# NRR_Normal_Norm_Fused_p15_k11 = funFEM(fdNRR_Fused_Normal_Norm_p15, K=11, model="DkBk", crit='bic', init='kmeans')
# 12-14 clusters all result in the same 'empty clusters' error message
# NRR_Normal_Norm_Fused_p15_k12 = funFEM(fdNRR_Fused_Normal_Norm_p15, K=12, model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_Norm_Fused_p15_k13 = funFEM(fdNRR_Fused_Normal_Norm_p15, K=13, model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_Norm_Fused_p15_k14 = funFEM(fdNRR_Fused_Normal_Norm_p15, K=14, model="DkBk", crit='bic', init='kmeans')
#####

# Putting all of the output into a list to save as an .RData file
# NRR_Normal_Norm_Fused_fem_list <- list(NRR_Normal_Norm_Fused_p15_k02,NRR_Normal_Norm_Fused_p15_k03,
#                                      NRR_Normal_Norm_Fused_p15_k04,NRR_Normal_Norm_Fused_p15_k05,
#                                      NRR_Normal_Norm_Fused_p15_k06,NRR_Normal_Norm_Fused_p15_k07,
#                                      NRR_Normal_Norm_Fused_p15_k08,NRR_Normal_Norm_Fused_p15_k09,
#                                      # NRR_Normal_Norm_Fused_p15_k10,  # error
#                                      NRR_Normal_Norm_Fused_p15_k11)

NRR_Normal_Norm_Fused_fem_list_path <- paste(filePath,"/rdata_files/NRR_Normal_Norm_Fused_fem_list.RData",sep = "")

# save(NRR_Normal_Norm_Fused_fem_list, file = NRR_Normal_Norm_Fused_fem_list_path)

load(NRR_Normal_Norm_Fused_fem_list_path)


##### Determine optimal k for Normal OG Fused Curves #####

NRR_Normal_OG_Fused_p07_k02_to_k06_bic <- matrix(data = NA, nrow=5, ncol=1)
NRR_Normal_OG_Fused_p07_k02_to_k06_aic <- matrix(data = NA, nrow=5, ncol=1)
NRR_Normal_OG_Fused_p07_k02_to_k06_icl <- matrix(data = NA, nrow=5, ncol=1)

##### Determine optimal k for Normal Norm Fused Curves  #####

NRR_Normal_Norm_Fused_p15_k02_to_k09_bic <- matrix(data = NA, nrow=8, ncol=1)
NRR_Normal_Norm_Fused_p15_k02_to_k09_aic <- matrix(data = NA, nrow=8, ncol=1)
NRR_Normal_Norm_Fused_p15_k02_to_k09_icl <- matrix(data = NA, nrow=8, ncol=1)

for(i in 1:5){
  NRR_Normal_OG_Fused_p07_k02_to_k06_bic[i] <- NRR_Normal_OG_Fused_fem_list[[i]]$bic
  NRR_Normal_OG_Fused_p07_k02_to_k06_aic[i] <- NRR_Normal_OG_Fused_fem_list[[i]]$aic
  NRR_Normal_OG_Fused_p07_k02_to_k06_icl[i] <- NRR_Normal_OG_Fused_fem_list[[i]]$icl
}
for(i in 1:8){
  NRR_Normal_Norm_Fused_p15_k02_to_k09_bic[i] <- NRR_Normal_Norm_Fused_fem_list[[i]]$bic
  NRR_Normal_Norm_Fused_p15_k02_to_k09_aic[i] <- NRR_Normal_Norm_Fused_fem_list[[i]]$aic
  NRR_Normal_Norm_Fused_p15_k02_to_k09_icl[i] <- NRR_Normal_Norm_Fused_fem_list[[i]]$icl
}

###############################################################################
###                Generate NRR K selection criteria plots                  ###
###############################################################################

# Note: The critical values for k, indicated by the vertical lines in the plots 
# were determined manually, and as such, are hard-coded into the plots in the 
# vertical abline() functions below.
best_K_NRR_Normal_OG_Fused_p07_plot_path <- paste(figurePath,
                                                  "/best_K_NRR_Normal_OG_Fused_p07_plot.pdf",
                                                  sep="")
# pdf(best_K_NRR_Normal_OG_Fused_p07_plot_path, width=6, height=6)

ymin <- min(NRR_Normal_OG_Fused_p07_k02_to_k06_bic[1:5,],
            NRR_Normal_OG_Fused_p07_k02_to_k06_aic[1:5,],
            NRR_Normal_OG_Fused_p07_k02_to_k06_icl[1:5,])
ymax <- max(NRR_Normal_OG_Fused_p07_k02_to_k06_bic[1:5,],
            NRR_Normal_OG_Fused_p07_k02_to_k06_aic[1:5,],
            NRR_Normal_OG_Fused_p07_k02_to_k06_icl[1:5,])


plot(seq(from=2,to=6,by=1), NRR_Normal_OG_Fused_p07_k02_to_k06_bic[1:5,], type='l', col=1,
     ylab="Selection Criteria", xlab="Nmb Clusters: K", ylim=c(ymin,ymax),
     main="NRR Normal (Original)Fused, p = 07; Model: DkBk")
lines(seq(from=2,to=6,by=1), NRR_Normal_OG_Fused_p07_k02_to_k06_aic[1:5,],col=2)
lines(seq(from=2,to=6,by=1), NRR_Normal_OG_Fused_p07_k02_to_k06_icl[1:5,],col=3)
# abline(v=6, lty=2, col="goldenrod")
legend(2,(ymax-500), legend = c("BIC","AIC","ICL"), col = c(1,2,3),
       lty = 1, cex = .8, y.intersp = 1)
# Looks like the maximum number of clusters (6) works best
# dev.off()



best_K_NRR_Normal_Norm_Fused_p15_plot_path <- paste(figurePath,
                                                  "/best_K_NRR_Normal_Norm_Fused_p15_plot.pdf",
                                                  sep="")
# pdf(best_K_NRR_Normal_Norm_Fused_p15_plot_path, width=6, height=6)

ymin <- min(NRR_Normal_Norm_Fused_p15_k02_to_k09_bic[1:8,],
            NRR_Normal_Norm_Fused_p15_k02_to_k09_aic[1:8,],
            NRR_Normal_Norm_Fused_p15_k02_to_k09_icl[1:8,])
ymax <- max(NRR_Normal_Norm_Fused_p15_k02_to_k09_bic[1:8,],
            NRR_Normal_Norm_Fused_p15_k02_to_k09_aic[1:8,],
            NRR_Normal_Norm_Fused_p15_k02_to_k09_icl[1:8,])


plot(seq(from=2,to=9,by=1), NRR_Normal_Norm_Fused_p15_k02_to_k09_bic[1:8,], type='l', col=1,
     ylab="Selection Criteria", xlab="Nmb Clusters: K", ylim=c(ymin,ymax),
     main="NRR Normal (Normalized) Fused, p = 15; Model: DkBk")
lines(seq(from=2,to=9,by=1), NRR_Normal_Norm_Fused_p15_k02_to_k09_aic[1:8,],col=2)
lines(seq(from=2,to=9,by=1), NRR_Normal_Norm_Fused_p15_k02_to_k09_icl[1:8,],col=3)
abline(v=8, lty=2, col="goldenrod")
legend(2,(ymax-500), legend = c("BIC","AIC","ICL"), col = c(1,2,3),
       lty = 1, cex = .8, y.intersp = 1)
# Looks like either 8 or the max works best
# dev.off()




#################################################################################
#################################################################################
#####                                                                       #####
#####  CIFU Clustering for the Phase Distance Filtered                      #####
#####   NRR Normal Fused Data - OG and Norm                                 #####
#####                                                                       #####
#################################################################################
#################################################################################

##### Original Curves (OG) #####
# We now filter the fused data points based on phase amplitude threshold
# the first 8 columns are identifier information, followed by 180 points for each
# of 668 observations
# Creating a functional data object, using 7 basis functions (just like the two
# fused curves each used)
# NRR_OG_Fused_p07_DxFiltered

# 272 obs
fdNRR_Fused_Normal_OG_p07_DxFiltered <- NRR_Normal_OG_twarp_dist_fmean_p07_full[NRR_OG_Fused_p07_DxFiltered,]
###########Save this for ALEX!!!!!  

fdNRR_Fused_Normal_OG_p07_Filtered <- Data2fd(argvals=d_NRR,
                                       y=t(fdNRR_Fused_Normal_OG_p07_DxFiltered[,-c(1:8)]),
                                       basisobj=NRR_basis_07)

# NRR_Normal_OG_Fused_p07_k02_Filtered = funFEM(fdNRR_Fused_Normal_OG_p07_Filtered, K=2,  model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_OG_Fused_p07_k03_Filtered = funFEM(fdNRR_Fused_Normal_OG_p07_Filtered, K=3,  model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_OG_Fused_p07_k04_Filtered = funFEM(fdNRR_Fused_Normal_OG_p07_Filtered, K=4,  model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_OG_Fused_p07_k05_Filtered = funFEM(fdNRR_Fused_Normal_OG_p07_Filtered, K=5,  model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_OG_Fused_p07_k06_Filtered = funFEM(fdNRR_Fused_Normal_OG_p07_Filtered, K=6,  model="DkBk", crit='bic', init='kmeans')
# #####
# 
# # Putting all of the output into a list to save as an .RData file
# NRR_Normal_OG_Fused_fem_list_Filtered <- list(NRR_Normal_OG_Fused_p07_k02_Filtered,
#                                               NRR_Normal_OG_Fused_p07_k03_Filtered,
#                                               NRR_Normal_OG_Fused_p07_k04_Filtered,
#                                               NRR_Normal_OG_Fused_p07_k05_Filtered,
#                                               NRR_Normal_OG_Fused_p07_k06_Filtered)

NRR_Normal_OG_Fused_fem_list_path_Filtered <- paste(filePath,"/rdata_files/NRR_Normal_OG_Fused_fem_list_Filtered.RData",sep = "")

# save(NRR_Normal_OG_Fused_fem_list_Filtered, file = NRR_Normal_OG_Fused_fem_list_path_Filtered)

load(NRR_Normal_OG_Fused_fem_list_path_Filtered)



##### Normalized Curves (Norm) #####
# View(NRR_Normal_Norm_twarp_dist_fmean_p15_full)
# Creating a functional data object, using 15 basis functions (just like the two
# fused curves each used) from the filtered data (filtered by phase distance)

# 305 obs
fdNRR_Fused_Normal_Norm_p15_DxFiltered <- NRR_Normal_Norm_twarp_dist_fmean_p15_full[NRR_Norm_Fused_p15_DxFiltered,]
###########Save this for ALEX!!!!!  

fdNRR_Fused_Normal_Norm_p15_Filtered <- Data2fd(argvals=d_NRR,
                                              y=t(fdNRR_Fused_Normal_Norm_p15_DxFiltered[,-c(1:8)]),
                                              basisobj=NRR_basis_15)

# NRR_Normal_Norm_Fused_p15_k02_Filtered = funFEM(fdNRR_Fused_Normal_Norm_p15_Filtered, K=2,  model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_Norm_Fused_p15_k03_Filtered = funFEM(fdNRR_Fused_Normal_Norm_p15_Filtered, K=3,  model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_Norm_Fused_p15_k04_Filtered = funFEM(fdNRR_Fused_Normal_Norm_p15_Filtered, K=4,  model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_Norm_Fused_p15_k05_Filtered = funFEM(fdNRR_Fused_Normal_Norm_p15_Filtered, K=5,  model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_Norm_Fused_p15_k06_Filtered = funFEM(fdNRR_Fused_Normal_Norm_p15_Filtered, K=6,  model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_Norm_Fused_p15_k07_Filtered = funFEM(fdNRR_Fused_Normal_Norm_p15_Filtered, K=7,  model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_Norm_Fused_p15_k08_Filtered = funFEM(fdNRR_Fused_Normal_Norm_p15_Filtered, K=8,  model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_Norm_Fused_p15_k09_Filtered = funFEM(fdNRR_Fused_Normal_Norm_p15_Filtered, K=9,  model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_Norm_Fused_p15_k10_Filtered = funFEM(fdNRR_Fused_Normal_Norm_p15_Filtered, K=10, model="DkBk", crit='bic', init='kmeans')
# NRR_Normal_Norm_Fused_p15_k11_Filtered = funFEM(fdNRR_Fused_Normal_Norm_p15_Filtered, K=11, model="DkBk", crit='bic', init='kmeans')
# #####
# 
# # Putting all of the output into a list to save as an .RData file
# NRR_Normal_Norm_Fused_fem_list_Filtered <- list(NRR_Normal_Norm_Fused_p15_k02_Filtered,
#                                                 NRR_Normal_Norm_Fused_p15_k03_Filtered,
#                                                 NRR_Normal_Norm_Fused_p15_k04_Filtered,
#                                                 NRR_Normal_Norm_Fused_p15_k05_Filtered,
#                                                 NRR_Normal_Norm_Fused_p15_k06_Filtered,
#                                                 NRR_Normal_Norm_Fused_p15_k07_Filtered,
#                                                 NRR_Normal_Norm_Fused_p15_k08_Filtered,
#                                                 NRR_Normal_Norm_Fused_p15_k09_Filtered,
#                                                 NRR_Normal_Norm_Fused_p15_k10_Filtered,
#                                                 NRR_Normal_Norm_Fused_p15_k11_Filtered)

NRR_Normal_Norm_Fused_fem_list_path_Filtered <- paste(filePath,"/rdata_files/NRR_Normal_Norm_Fused_fem_list_Filtered.RData",sep = "")

# save(NRR_Normal_Norm_Fused_fem_list_Filtered, file = NRR_Normal_Norm_Fused_fem_list_path_Filtered)

load(NRR_Normal_Norm_Fused_fem_list_path_Filtered)


##### Determine optimal k for Normal OG Fused Curves #####

NRR_Normal_OG_Fused_p07_k02_to_k06_bic_Filtered <- matrix(data = NA, nrow=5, ncol=1)
NRR_Normal_OG_Fused_p07_k02_to_k06_aic_Filtered <- matrix(data = NA, nrow=5, ncol=1)
NRR_Normal_OG_Fused_p07_k02_to_k06_icl_Filtered <- matrix(data = NA, nrow=5, ncol=1)

##### Determine optimal k for Normal Norm Fused Curves  #####

NRR_Normal_Norm_Fused_p15_k02_to_k11_bic_Filtered <- matrix(data = NA, nrow=10, ncol=1)
NRR_Normal_Norm_Fused_p15_k02_to_k11_aic_Filtered <- matrix(data = NA, nrow=10, ncol=1)
NRR_Normal_Norm_Fused_p15_k02_to_k11_icl_Filtered <- matrix(data = NA, nrow=10, ncol=1)

for(i in 1:5){
  NRR_Normal_OG_Fused_p07_k02_to_k06_bic_Filtered[i] <- NRR_Normal_OG_Fused_fem_list_Filtered[[i]]$bic
  NRR_Normal_OG_Fused_p07_k02_to_k06_aic_Filtered[i] <- NRR_Normal_OG_Fused_fem_list_Filtered[[i]]$aic
  NRR_Normal_OG_Fused_p07_k02_to_k06_icl_Filtered[i] <- NRR_Normal_OG_Fused_fem_list_Filtered[[i]]$icl
}
for(i in 1:10){
  NRR_Normal_Norm_Fused_p15_k02_to_k11_bic_Filtered[i] <- NRR_Normal_Norm_Fused_fem_list_Filtered[[i]]$bic
  NRR_Normal_Norm_Fused_p15_k02_to_k11_aic_Filtered[i] <- NRR_Normal_Norm_Fused_fem_list_Filtered[[i]]$aic
  NRR_Normal_Norm_Fused_p15_k02_to_k11_icl_Filtered[i] <- NRR_Normal_Norm_Fused_fem_list_Filtered[[i]]$icl
}

###############################################################################
###                Generate NRR K selection criteria plots                  ###
###############################################################################

# Note: The critical values for k, indicated by the vertical lines in the plots 
# were determined manually, and as such, are hard-coded into the plots in the 
# vertical abline() functions below.
best_K_NRR_Normal_OG_Fused_p07_plot_path_Filtered <- paste(figurePath,
                                                           "/best_K_NRR_Normal_OG_Fused_p07_plot_Filtered.pdf",
                                                           sep="")
# pdf(best_K_NRR_Normal_OG_Fused_p07_plot_path_Filtered, width=7, height=6)

ymin <- min(NRR_Normal_OG_Fused_p07_k02_to_k06_bic_Filtered[1:5,],
            NRR_Normal_OG_Fused_p07_k02_to_k06_aic_Filtered[1:5,],
            NRR_Normal_OG_Fused_p07_k02_to_k06_icl_Filtered[1:5,])
ymax <- max(NRR_Normal_OG_Fused_p07_k02_to_k06_bic_Filtered[1:5,],
            NRR_Normal_OG_Fused_p07_k02_to_k06_aic_Filtered[1:5,],
            NRR_Normal_OG_Fused_p07_k02_to_k06_icl_Filtered[1:5,])


plot(seq(from=2,to=6,by=1), NRR_Normal_OG_Fused_p07_k02_to_k06_bic_Filtered[1:5,], type='l', col=1,
     ylab="Selection Criteria", xlab="Nmb Clusters: K", ylim=c(ymin,ymax),
     main="NRR Normal (Original) Fused (Dx Filtered), p = 07; Model: DkBk")
lines(seq(from=2,to=6,by=1), NRR_Normal_OG_Fused_p07_k02_to_k06_aic_Filtered[1:5,],col=2)
lines(seq(from=2,to=6,by=1), NRR_Normal_OG_Fused_p07_k02_to_k06_icl_Filtered[1:5,],col=3)
abline(v=5, lty=2, col="goldenrod")
legend(2,(ymax-100), legend = c("BIC","AIC","ICL"), col = c(1,2,3),
       lty = 1, cex = .8, y.intersp = 1)
# Looks like the maximum number of clusters (6) works best

# dev.off()



best_K_NRR_Normal_Norm_Fused_p15_plot_path_Filtered <- paste(figurePath,
                                                             "/best_K_NRR_Normal_Norm_Fused_p15_plot_Filtered.pdf",
                                                             sep="")
# pdf(best_K_NRR_Normal_Norm_Fused_p15_plot_path_Filtered, width=7, height=6)

ymin <- min(NRR_Normal_Norm_Fused_p15_k02_to_k11_bic_Filtered[1:10,],
            NRR_Normal_Norm_Fused_p15_k02_to_k11_aic_Filtered[1:10,],
            NRR_Normal_Norm_Fused_p15_k02_to_k11_icl_Filtered[1:10,])
ymax <- max(NRR_Normal_Norm_Fused_p15_k02_to_k11_bic_Filtered[1:10,],
            NRR_Normal_Norm_Fused_p15_k02_to_k11_aic_Filtered[1:10,],
            NRR_Normal_Norm_Fused_p15_k02_to_k11_icl_Filtered[1:10,])


plot(seq(from=2,to=11,by=1), NRR_Normal_Norm_Fused_p15_k02_to_k11_bic_Filtered[1:10,], type='l', col=1,
     ylab="Selection Criteria", xlab="Nmb Clusters: K", ylim=c(ymin,ymax),
     main="NRR Normal (Normalized) Fused (Dx Filtered), p = 15; Model: DkBk")
lines(seq(from=2,to=11,by=1), NRR_Normal_Norm_Fused_p15_k02_to_k11_aic_Filtered[1:10,],col=2)
lines(seq(from=2,to=11,by=1), NRR_Normal_Norm_Fused_p15_k02_to_k11_icl_Filtered[1:10,],col=3)
abline(v=4, lty=2, col="goldenrod")
legend(2,(ymax-100), legend = c("BIC","AIC","ICL"), col = c(1,2,3),
       lty = 1, cex = .8, y.intersp = 1)
# Looks like 4 works best
# dev.off()


#################################################################################
#################################################################################
#################################################################################
#This is the winner. We want to analyze the output of this model


NRR_Normal_Norm_Fused_p15_k04_Filtered <- NRR_Normal_Norm_Fused_fem_list_Filtered[[3]]

# NRR_Normal_Norm_Fused_p15_k04_Filtered$cls - cluster assignments
# NRR_Norm_Fused_p15_DxFiltered
NRR_Fused_Normal_Norm_p15_DxFiltered <- NRR_Normal_Norm_twarp_dist_fmean_p15_full[NRR_Norm_Fused_p15_DxFiltered,]

NRR_OCT_Normal_Norm_p15_eval_full_filtered <- NRR_OCT_Normal_Norm_p15_eval_full[NRR_Norm_Fused_p15_DxFiltered,]


#################################################################################
#################################################################################
#################################################################################



# Subset for Original Curves eyes - 
# Elastic Distance OCT/Fundus: Normal Eyes (Original Curves) phase distance <= 0.355
# NRR_OG_Fused_p07_DxFiltered <- (NRR_OG_DxDy$PhaseDist <= 0.355)
# 272 obs
NRR_Fused_Normal_OG_p07_fullData_DxFiltered <- NRR_Normal_OG_twarp_dist_fmean_p07_full[NRR_OG_Fused_p07_DxFiltered,]

# Subset for Normalized eyes - 
# Elastic Distance OCT/Fundus: Normal Eyes (Normalized Curves) phase distance <= 0.44
# NRR_Norm_Fused_p15_DxFiltered <- (NRR_Norm_DxDy$PhaseDist <= 0.44)

# 305 obs
NRR_Fused_Normal_Norm_p15_fullData_DxFiltered <- NRR_Normal_Norm_twarp_dist_fmean_p15_full[NRR_Norm_Fused_p15_DxFiltered,]


#################################################################################
#################################################################################
# plot this by the cluster
NRR_Fused_Normal_Norm_p15_fullData_DxFiltered_clust_k04 <- cbind(NRR_Fused_Normal_Norm_p15_fullData_DxFiltered[,1:8], 
                                                                 cls=NRR_Normal_Norm_Fused_p15_k04_Filtered$cls, 
                                                                 NRR_Fused_Normal_Norm_p15_fullData_DxFiltered[,9:188])

# View(NRR_Fused_Normal_Norm_p15_fullData_DxFiltered_clust_k04)
# write_xlsx(NRR_Fused_Normal_Norm_p15_fullData_DxFiltered_clust_k04, 
#            path = paste(localPath_tables,"NRR_Fused_Normal_Norm_p15_fullData_DxFiltered_clust_k04.xlsx"))


#################################################################################
#################################################################################



# Phase_dist_filtered_data <- list(NRR_Fused_Normal_OG_p07_fullData_DxFiltered,
#                                  NRR_Fused_Normal_Norm_p15_fullData_DxFiltered)

Phase_dist_filtered_data_path <- paste(filePath,"/rdata_files/Phase_dist_filtered_data.RData",sep = "")
# save(Phase_dist_filtered_data, file = Phase_dist_filtered_data_path)
load(Phase_dist_filtered_data_path)
############################################################


# # Skin study 
# MMa_DxDy <- data.frame(PhaseDist=as.numeric(MMa_Dx), Amplitude=as.numeric(MMa_Dy))
# # Dx: phase distance
# # Dy: amplitude distance
# # qplot is deprecated. Updated version below.
# # scatter <- qplot(MMa_Dx,MMa_Dy, 
# #                  xlab = "Phase Distance", 
# #                  ylab = "Amplitude Distance", 
# #                  main = "Elastic Distances (Normalized Curves)")  + 
# #   scale_x_continuous(limits=c(min(MMa_Dx),max(MMa_Dx))) + 
# #   scale_y_continuous(limits=c(min(MMa_Dy),max(MMa_Dy))) + 
# #   geom_rug(col=rgb(.5,0,0,alpha=.2))
# # scatter
# 
# p <- ggplot(MMa_DxDy, aes(x=PhaseDist, y=Amplitude)) + geom_point() + theme_classic()
# ggExtra::ggMarginal(p, type = "histogram")
# 

###########
# try to replicate the plots from previous clustering paper
# but with NRR_Normal_Norm_Fused_p15_k04_Filtered$cls
# the Filtered, Normal eyes, normalized curves, clustered into 
# 4 clusters by the funFem algorithm, utilizing the DkBk model


plot.mean.curves.on.circle = function(curves, main = NA, col = 1){
  nmb_pts_per_curve <- dim(curves)[2]
  d = 1:nmb_pts_per_curve
  dat = curves #curves  # eval.fd(d, fd)    # Convert fd to data matrix (181 by N)
  xc = cos((d/nmb_pts_per_curve)*2*pi)
  ys = sin((d/nmb_pts_per_curve)*2*pi)
  Xc.mat = t(dat) * xc
  # Xc.mat <- c(Xc.mat, Xc.mat[1])
  Ys.mat = t(dat) * ys
  # Ys.mat <- c(Ys.mat, Ys.mat[1])
  xy.lim = 0.010  # hard coded to better fit the plots
  # xy.lim = max(abs(Xc.mat), abs(Ys.mat)) + 0.0025  #original
  print(xy.lim)
  # Find the average curve
  X_bar = apply(curves, 2, sum)/dim(curves)[1]
  Xc.bar = X_bar * xc
  Ys.bar = X_bar * ys
  
  # plot the average curve last
  plot(c(Xc.mat[,1],Xc.mat[1,1]), c(Ys.mat[,1],Ys.mat[1,1]), type = "l", main = main, col = col,
       xlim = c(-xy.lim, xy.lim), ylim = c(-xy.lim, xy.lim), 
       xlab = NA, ylab = NA)
  
  for(i in 2:nrow(dat)){
    lines(c(Xc.mat[,i],Xc.mat[1,i]), c(Ys.mat[,i],Ys.mat[1,i]), col = col)  # (i-1)%%10+2)
  }
  
  lines(c(Xc.bar,Xc.bar[1]), c(Ys.bar,Ys.bar[1]), col = 1, lwd=2)
  
  # # plot the average curve
  # plot(Xc.bar, Ys.bar, col = col, lwd=3, type = "l", main = main,
  #      xlim = c(-xy.lim, xy.lim), ylim = c(-xy.lim, xy.lim),
  #      xlab = NA, ylab = NA)
  
  mtext("T", side = 4, cex = 1.25, col = 2, font = 2, las = 2)
  mtext("I", side = 1, cex = 1.25, col = 2, font = 2)
  mtext("N", side = 2, cex = 1.25, col = 2, font = 2, las = 2)
  mtext("S", side = 3, cex = 1.25, col = 2, font = 2)
  
  abline(v = 0, lty = "dashed")
  abline(h = 0, lty = "dashed")
  abline(0, 1, lty = "dashed")
  abline(0, -1, lty = "dashed")
}

##### NRR Curves by Cluster  #####
nrr_cluster_colors <- c("coral2","darkgoldenrod","chartreuse3","deeppink1","forestgreen",
                        "cyan3","deepskyblue1","dodgerblue2","darkorchid1")


##### NRR Curves by Cluster  #####
nrr_normal_norm_fused_p15_k04_cls_1 <- NRR_Fused_Normal_Norm_p15_fullData_DxFiltered_clust_k04[NRR_Fused_Normal_Norm_p15_fullData_DxFiltered_clust_k04$cls==1, 10:189]
nrr_normal_norm_fused_p15_k04_cls_2 <- NRR_Fused_Normal_Norm_p15_fullData_DxFiltered_clust_k04[NRR_Fused_Normal_Norm_p15_fullData_DxFiltered_clust_k04$cls==2, 10:189]
nrr_normal_norm_fused_p15_k04_cls_3 <- NRR_Fused_Normal_Norm_p15_fullData_DxFiltered_clust_k04[NRR_Fused_Normal_Norm_p15_fullData_DxFiltered_clust_k04$cls==3, 10:189]
nrr_normal_norm_fused_p15_k04_cls_4 <- NRR_Fused_Normal_Norm_p15_fullData_DxFiltered_clust_k04[NRR_Fused_Normal_Norm_p15_fullData_DxFiltered_clust_k04$cls==4, 10:189]

NRR_Normal_Norm_fused_p15_k04_cls_1_plot_path <- paste(figurePath,"/NRR_Normal_Norm_fused_p15_k04_cls_1_plot.pdf",sep="")
# pdf(NRR_Normal_Norm_fused_p15_k04_cls_1_plot_path, width=6, height=6)
plot.mean.curves.on.circle(nrr_normal_norm_fused_p15_k04_cls_1, main = "NRR Normal Normalized Fused Cluster 1", col = nrr_cluster_colors[1])
# dev.off()

NRR_Normal_Norm_fused_p15_k04_cls_2_plot_path <- paste(figurePath,"/NRR_Normal_Norm_fused_p15_k04_cls_2_plot.pdf",sep="")
# pdf(NRR_Normal_Norm_fused_p15_k04_cls_2_plot_path, width=6, height=6)
plot.mean.curves.on.circle(nrr_normal_norm_fused_p15_k04_cls_2, main = "NRR Normal Normalized Fused Cluster 2", col = nrr_cluster_colors[2])
# dev.off()

NRR_Normal_Norm_fused_p15_k04_cls_3_plot_path <- paste(figurePath,"/NRR_Normal_Norm_fused_p15_k04_cls_3_plot.pdf",sep="")
# pdf(NRR_Normal_Norm_fused_p15_k04_cls_3_plot_path, width=6, height=6)
plot.mean.curves.on.circle(nrr_normal_norm_fused_p15_k04_cls_3, main = "NRR Normal Normalized Fused Cluster 3", col = nrr_cluster_colors[3])
# dev.off()

NRR_Normal_Norm_fused_p15_k04_cls_4_plot_path <- paste(figurePath,"/NRR_Normal_Norm_fused_p15_k04_cls_4_plot.pdf",sep="")
# pdf(NRR_Normal_Norm_fused_p15_k04_cls_4_plot_path, width=6, height=6)
plot.mean.curves.on.circle(nrr_normal_norm_fused_p15_k04_cls_4, main = "NRR Normal Normalized Fused Cluster 4", col = nrr_cluster_colors[4])
# dev.off()

NRR_Normal_Norm_fused_p15_k04_cls_plots_path <- paste(figurePath,"NRR_Normal_Norm_fused_p15_k04_cls_plots.pdf",sep="")
# pdf(NRR_Normal_Norm_fused_p15_k04_cls_plots_path, width=8, height=8)
par(mfrow=c(2,2))
plot.mean.curves.on.circle(nrr_normal_norm_fused_p15_k04_cls_1, main = "Cluster 1", col = nrr_cluster_colors[1])
plot.mean.curves.on.circle(nrr_normal_norm_fused_p15_k04_cls_2, main = "Cluster 2", col = nrr_cluster_colors[2])
plot.mean.curves.on.circle(nrr_normal_norm_fused_p15_k04_cls_3, main = "Cluster 3", col = nrr_cluster_colors[3])
plot.mean.curves.on.circle(nrr_normal_norm_fused_p15_k04_cls_4, main = "Cluster 4", col = nrr_cluster_colors[4])
par(mfrow=c(1,1))
# dev.off()



#################################################################################
#################################################################################
# Double filtered - DxDy - sample curves
#################################################################################
#################################################################################

plot.example.curves.on.circle = function(exmpl_dat, main = NA, col = 1){
  nmb_pts_per_curve <- dim(exmpl_dat)[2]
  d = 1:nmb_pts_per_curve 
  dat_1 <- exmpl_dat[1,]  # OCT curve
  dat_2 <- exmpl_dat[2,]  # Fundus curve
  dat_3 <- exmpl_dat[3,]  # fmean curve
  xc = cos((d/nmb_pts_per_curve)*2*pi)
  ys = sin((d/nmb_pts_per_curve)*2*pi)
  Xc.mat_1 = t(dat_1) * xc
  Xc.mat_1 <- c(Xc.mat_1, Xc.mat_1[1])
  Ys.mat_1 = t(dat_1) * ys
  Ys.mat_1 <- c(Ys.mat_1, Ys.mat_1[1])
  Xc.mat_2 = t(dat_2) * xc
  Xc.mat_2 <- c(Xc.mat_2, Xc.mat_2[1])
  Ys.mat_2 = t(dat_2) * ys
  Ys.mat_2 <- c(Ys.mat_2, Ys.mat_2[1])
  xy.lim = 0.010
  
  
  # fmean curve
  X_bar = t(dat_3)  # need this to be 1 X 180
  Xc.bar = X_bar * xc
  Xc.bar <- c(Xc.bar, Xc.bar[1])
  Ys.bar = X_bar * ys
  Ys.bar <- c(Ys.bar, Ys.bar[1])
  
  # plot the average curve last
  par(mfrow = c(1,1))
  plot(Xc.mat_1, Ys.mat_1, type = "l", col = "blue",
       xlim = c(-xy.lim, xy.lim), ylim = c(-xy.lim, xy.lim), 
       xlab = NA, ylab = NA, main=paste("OCT (blue) and Fundus (red), and Fused (magenta) \n", main, sep = ""))
  
  
  lines(Xc.mat_2, Ys.mat_2, col = "red")
  
  lines(Xc.bar, Ys.bar, col = "magenta", lwd=2)
  
  abline(v = 0, lty = "dashed")
  abline(h = 0, lty = "dashed")
  abline(0, 1, lty = "dashed")
  abline(0, -1, lty = "dashed")
}

###############################################################################
###############################################################################
###############################################################################













# NRR_OCT_Normal_Norm_p15_evaluated  # dim = 180 668
# NRR_Fundus_Normal_Norm_p15_evaluated  # dim = 180 668

plot(NRR_Norm_DxDy$Amplitude,NRR_Norm_DxDy$PhaseDist)
abline(v=median(NRR_Norm_DxDy$Amplitude))
abline(h=median(NRR_Norm_DxDy$PhaseDist))

# median(NRR_Norm_DxDy$Amplitude)
# 0.05387742

# median(NRR_Norm_DxDy$PhaseDist)
# 0.4469626

dxdy_Filtered <- (NRR_Norm_DxDy$Amplitude < median(NRR_Norm_DxDy$Amplitude) & NRR_Norm_DxDy$PhaseDist < median(NRR_Norm_DxDy$PhaseDist))
NRR_Norm_Fused_p15_dxdy_Filtered <- NRR_Normal_Norm_twarp_dist_fmean_p15_full[dxdy_Filtered, ]
# View(NRR_Norm_Fused_p15_dxdy_Filtered)

Dxdy_Filtered <- (NRR_Norm_DxDy$Amplitude < median(NRR_Norm_DxDy$Amplitude) & NRR_Norm_DxDy$PhaseDist > median(NRR_Norm_DxDy$PhaseDist))
NRR_Norm_Fused_p15_Dxdy_Filtered <- NRR_Normal_Norm_twarp_dist_fmean_p15_full[Dxdy_Filtered, ]
# View(NRR_Norm_Fused_p15_Dxdy_Filtered)

dxDy_Filtered <- (NRR_Norm_DxDy$Amplitude > median(NRR_Norm_DxDy$Amplitude) & NRR_Norm_DxDy$PhaseDist < median(NRR_Norm_DxDy$PhaseDist))
NRR_Norm_Fused_p15_dxDy_Filtered <- NRR_Normal_Norm_twarp_dist_fmean_p15_full[dxDy_Filtered, ]
# View(NRR_Norm_Fused_p15_dxDy_Filtered)

DxDy_Filtered <- (NRR_Norm_DxDy$Amplitude > median(NRR_Norm_DxDy$Amplitude) & NRR_Norm_DxDy$PhaseDist > median(NRR_Norm_DxDy$PhaseDist))
NRR_Norm_Fused_p15_DxDy_Filtered <- NRR_Normal_Norm_twarp_dist_fmean_p15_full[DxDy_Filtered, ]
# View(NRR_Norm_Fused_p15_DxDy_Filtered)

#################################################################################
#####      Troubleshooting

# 
# dx_Filtered <- (NRR_Norm_DxDy$Amplitude < median(NRR_Norm_DxDy$Amplitude))
# sum(dx_Filtered)
# 
# dy_Filtered <- (NRR_Norm_DxDy$PhaseDist < median(NRR_Norm_DxDy$PhaseDist))
# sum(dy_Filtered)
# 
# dxdy_Filtered <- (NRR_Norm_DxDy$Amplitude < median(NRR_Norm_DxDy$Amplitude) & NRR_Norm_DxDy$PhaseDist < median(NRR_Norm_DxDy$PhaseDist))
# sum(dxdy_Filtered)
# 
# Dxdy_Filtered <- (NRR_Norm_DxDy$Amplitude > median(NRR_Norm_DxDy$Amplitude) & NRR_Norm_DxDy$PhaseDist < median(NRR_Norm_DxDy$PhaseDist))
# sum(Dxdy_Filtered)
# 
# dxDy_Filtered <- (NRR_Norm_DxDy$Amplitude < median(NRR_Norm_DxDy$Amplitude) & NRR_Norm_DxDy$PhaseDist > median(NRR_Norm_DxDy$PhaseDist))
# sum(dxDy_Filtered)
# 
# DxDy_Filtered <- (NRR_Norm_DxDy$Amplitude > median(NRR_Norm_DxDy$Amplitude) & NRR_Norm_DxDy$PhaseDist > median(NRR_Norm_DxDy$PhaseDist))
# sum(DxDy_Filtered)
# 
# 
# NRR_Normal_Norm_twarp_dist_fmean_p15_full[(NRR_Norm_DxDy$Amplitude < median(NRR_Norm_DxDy$Amplitude) & NRR_Norm_DxDy$PhaseDist < median(NRR_Norm_DxDy$PhaseDist)), ]
# (NRR_Norm_DxDy$Amplitude < median(NRR_Norm_DxDy$Amplitude) & NRR_Norm_DxDy$PhaseDist < median(NRR_Norm_DxDy$PhaseDist))
# 
# NRR_Norm_DxDy[1:5,]
# 

# sample(NRR_Norm_Fused_p15_dxdy_Filtered$PatientID_EYE, size = 1)




#  Randomly sample eyes from the different combinations of above and below 
#  the mean of the amplitude and phase distances
# (egEyeCurves_Fundus_OCT_Fused_dxdy <- sample(NRR_Norm_Fused_p15_dxdy_Filtered$PatientID_EYE, size = 1))  # 18845OS
# (egEyeCurves_Fundus_OCT_Fused_Dxdy <- sample(NRR_Norm_Fused_p15_Dxdy_Filtered$PatientID_EYE, size = 1))  # 29102OD
# (egEyeCurves_Fundus_OCT_Fused_dxDy <- sample(NRR_Norm_Fused_p15_dxDy_Filtered$PatientID_EYE, size = 1))  # 14026OS
# (egEyeCurves_Fundus_OCT_Fused_DxDy <- sample(NRR_Norm_Fused_p15_DxDy_Filtered$PatientID_EYE, size = 1))  # 16423OS

egEyeCurves_Fundus_OCT_Fused_dxdy <- '18845OS'
egEyeCurves_Fundus_OCT_Fused_Dxdy <- '29102OD'
egEyeCurves_Fundus_OCT_Fused_dxDy <- '14026OS'
egEyeCurves_Fundus_OCT_Fused_DxDy <- '16423OS'


# dim = 188 variables & 668 observations

eg_OCT_dxdy <- NRR_Fundus_Normal_Norm_p15_eval_full[NRR_Fundus_Normal_Norm_p15_eval_full$PatientID_EYE == egEyeCurves_Fundus_OCT_Fused_dxdy, ]  
eg_Fundus_dxdy <- NRR_OCT_Normal_Norm_p15_eval_full[NRR_OCT_Normal_Norm_p15_eval_full$PatientID_EYE == egEyeCurves_Fundus_OCT_Fused_dxdy, ]  
eg_Fused_dxdy <- NRR_Normal_Norm_twarp_dist_fmean_p15_full[NRR_Normal_Norm_twarp_dist_fmean_p15_full$PatientID_EYE == egEyeCurves_Fundus_OCT_Fused_dxdy, ]  

eg_dxdy <- rbind(eg_OCT_dxdy[9:188], eg_Fundus_dxdy[9:188], eg_Fused_dxdy[9:188])
eg_dxdy_amplitude <- round(NRR_Norm_DxDy$Amplitude[which(NRR_Fundus_Normal_Norm_p15_eval_full$PatientID_EYE == egEyeCurves_Fundus_OCT_Fused_dxdy)],digits = 4)
eg_dxdy_phaseDist <- round(NRR_Norm_DxDy$PhaseDist[which(NRR_Fundus_Normal_Norm_p15_eval_full$PatientID_EYE == egEyeCurves_Fundus_OCT_Fused_dxdy)],digits = 4)

plot.example.curves.on.circle(eg_dxdy, main = paste("Phase Distance < 0.447 (med): ", eg_dxdy_phaseDist, " ; Amplitude < 0.0539 (med): ", eg_dxdy_amplitude))

# median(NRR_Norm_DxDy$Amplitude)
# 0.05387742

# median(NRR_Norm_DxDy$PhaseDist)
# 0.4469626
#####

eg_Fundus_Dxdy <- NRR_OCT_Normal_Norm_p15_eval_full[NRR_OCT_Normal_Norm_p15_eval_full$PatientID_EYE == egEyeCurves_Fundus_OCT_Fused_Dxdy, ]  
eg_OCT_Dxdy <- NRR_Fundus_Normal_Norm_p15_eval_full[NRR_Fundus_Normal_Norm_p15_eval_full$PatientID_EYE == egEyeCurves_Fundus_OCT_Fused_Dxdy, ]  
eg_Fused_Dxdy <- NRR_Normal_Norm_twarp_dist_fmean_p15_full[NRR_Normal_Norm_twarp_dist_fmean_p15_full$PatientID_EYE == egEyeCurves_Fundus_OCT_Fused_Dxdy, ]  

eg_Dxdy <- rbind(eg_OCT_Dxdy[9:188], eg_Fundus_Dxdy[9:188], eg_Fused_Dxdy[9:188])
eg_Dxdy_amplitude <- round(NRR_Norm_DxDy$Amplitude[which(NRR_Fundus_Normal_Norm_p15_eval_full$PatientID_EYE == egEyeCurves_Fundus_OCT_Fused_Dxdy)],digits = 4)
eg_Dxdy_phaseDist <- round(NRR_Norm_DxDy$PhaseDist[which(NRR_Fundus_Normal_Norm_p15_eval_full$PatientID_EYE == egEyeCurves_Fundus_OCT_Fused_Dxdy)],digits = 4)

plot.example.curves.on.circle(eg_Dxdy, main = paste("Phase Distance > 0.447 (med): ", eg_Dxdy_phaseDist, " ; Amplitude < 0.0539 (med): ", eg_Dxdy_amplitude))

#####

eg_Fundus_dxDy <- NRR_OCT_Normal_Norm_p15_eval_full[NRR_OCT_Normal_Norm_p15_eval_full$PatientID_EYE == egEyeCurves_Fundus_OCT_Fused_dxDy, ]  
eg_OCT_dxDy <- NRR_Fundus_Normal_Norm_p15_eval_full[NRR_Fundus_Normal_Norm_p15_eval_full$PatientID_EYE == egEyeCurves_Fundus_OCT_Fused_dxDy, ]  
eg_Fused_dxDy <- NRR_Normal_Norm_twarp_dist_fmean_p15_full[NRR_Normal_Norm_twarp_dist_fmean_p15_full$PatientID_EYE == egEyeCurves_Fundus_OCT_Fused_dxDy, ]  


eg_dxDy <- rbind(eg_OCT_dxDy[9:188], eg_Fundus_dxDy[9:188], eg_Fused_dxDy[9:188])
eg_dxDy_amplitude <- round(NRR_Norm_DxDy$Amplitude[which(NRR_Fundus_Normal_Norm_p15_eval_full$PatientID_EYE == egEyeCurves_Fundus_OCT_Fused_dxDy)],digits = 4)
eg_dxDy_phaseDist <- round(NRR_Norm_DxDy$PhaseDist[which(NRR_Fundus_Normal_Norm_p15_eval_full$PatientID_EYE == egEyeCurves_Fundus_OCT_Fused_dxDy)],digits = 4)

plot.example.curves.on.circle(eg_dxDy, main = paste("Phase Distance < 0.447 (med): ", eg_dxDy_phaseDist, " ; Amplitude > 0.0539 (med): ", eg_dxDy_amplitude))


#####

eg_Fundus_DxDy <- NRR_OCT_Normal_Norm_p15_eval_full[NRR_OCT_Normal_Norm_p15_eval_full$PatientID_EYE == egEyeCurves_Fundus_OCT_Fused_DxDy, ]  
eg_OCT_DxDy <- NRR_Fundus_Normal_Norm_p15_eval_full[NRR_Fundus_Normal_Norm_p15_eval_full$PatientID_EYE == egEyeCurves_Fundus_OCT_Fused_DxDy, ]  
eg_Fused_DxDy <- NRR_Normal_Norm_twarp_dist_fmean_p15_full[NRR_Normal_Norm_twarp_dist_fmean_p15_full$PatientID_EYE == egEyeCurves_Fundus_OCT_Fused_DxDy, ]  

eg_DxDy <- rbind(eg_OCT_DxDy[9:188], eg_Fundus_DxDy[9:188], eg_Fused_DxDy[9:188])
eg_DxDy_amplitude <- round(NRR_Norm_DxDy$Amplitude[which(NRR_Fundus_Normal_Norm_p15_eval_full$PatientID_EYE == egEyeCurves_Fundus_OCT_Fused_DxDy)],digits = 4)
eg_DxDy_phaseDist <- round(NRR_Norm_DxDy$PhaseDist[which(NRR_Fundus_Normal_Norm_p15_eval_full$PatientID_EYE == egEyeCurves_Fundus_OCT_Fused_DxDy)],digits = 4)

plot.example.curves.on.circle(eg_DxDy, main = paste("Phase Distance > 0.447 (med): ", eg_DxDy_phaseDist, " ; Amplitude > 0.0539 (med): ", eg_DxDy_amplitude))


NRR_Normal_Norm_OCT_Fundus_Fused_p15_plots_path <- paste(figurePath,"NRR_Normal_Norm_OCT_Fundus_Fused_p15_plots.pdf",sep="")
# pdf(NRR_Normal_Norm_OCT_Fundus_Fused_p15_plots_path, width=6, height=12)

par(mfrow = c(2,2))

plot.example.curves.on.circle(eg_dxdy, main = paste("Phase Distance < 0.447 (med): ", eg_dxdy_phaseDist, " ; Amplitude < 0.0539 (med): ", eg_dxdy_amplitude))
plot.example.curves.on.circle(eg_Dxdy, main = paste("Phase Distance > 0.447 (med): ", eg_Dxdy_phaseDist, " ; Amplitude < 0.0539 (med): ", eg_Dxdy_amplitude))
plot.example.curves.on.circle(eg_dxDy, main = paste("Phase Distance < 0.447 (med): ", eg_dxDy_phaseDist, " ; Amplitude > 0.0539 (med): ", eg_dxDy_amplitude))
plot.example.curves.on.circle(eg_DxDy, main = paste("Phase Distance > 0.447 (med): ", eg_DxDy_phaseDist, " ; Amplitude > 0.0539 (med): ", eg_DxDy_amplitude))
par(mfrow=c(1,1))
# dev.off()


NRR_Normal_Norm_OCT_Fundus_Fused_p15_dxdy_plot_path <- paste(figurePath,"NRR_Normal_Norm_OCT_Fundus_Fused_p15_dxdy_plot1.pdf",sep="")
# pdf(NRR_Normal_Norm_OCT_Fundus_Fused_p15_dxdy_plot_path, width=8, height=6)
plot.example.curves.on.circle(eg_dxdy, main = paste("Phase Distance < 0.447 (med): ", eg_dxdy_phaseDist, " ; Amplitude < 0.0539 (med): ", eg_dxdy_amplitude))
# dev.off()

NRR_Normal_Norm_OCT_Fundus_Fused_p15_Dxdy_plot_path <- paste(figurePath,"NRR_Normal_Norm_OCT_Fundus_Fused_p15_Dxdy_plot2.pdf",sep="")
# pdf(NRR_Normal_Norm_OCT_Fundus_Fused_p15_Dxdy_plot_path, width=8, height=6)
plot.example.curves.on.circle(eg_Dxdy, main = paste("Phase Distance > 0.447 (med): ", eg_Dxdy_phaseDist, " ; Amplitude < 0.0539 (med): ", eg_Dxdy_amplitude))
# dev.off()

NRR_Normal_Norm_OCT_Fundus_Fused_p15_dxDy_plot_path <- paste(figurePath,"NRR_Normal_Norm_OCT_Fundus_Fused_p15_dxDy_plot3.pdf",sep="")
# pdf(NRR_Normal_Norm_OCT_Fundus_Fused_p15_dxDy_plot_path, width=8, height=6)
plot.example.curves.on.circle(eg_dxDy, main = paste("Phase Distance < 0.447 (med): ", eg_dxDy_phaseDist, " ; Amplitude > 0.0539 (med): ", eg_dxDy_amplitude))
# dev.off()

NRR_Normal_Norm_OCT_Fundus_Fused_p15_DxDy_plot_path <- paste(figurePath,"NRR_Normal_Norm_OCT_Fundus_Fused_p15_DxDy_plot4.pdf",sep="")
# pdf(NRR_Normal_Norm_OCT_Fundus_Fused_p15_DxDy_plot_path, width=8, height=6)
plot.example.curves.on.circle(eg_DxDy, main = paste("Phase Distance > 0.447 (med): ", eg_DxDy_phaseDist, " ; Amplitude > 0.0539 (med): ", eg_DxDy_amplitude))
# dev.off()


par(mfrow=c(2,2))
plot.mean.curves.on.circle(nrr_normal_norm_fused_p15_k04_cls_1, main = "Cluster 1", col = nrr_cluster_colors[1])
plot.mean.curves.on.circle(nrr_normal_norm_fused_p15_k04_cls_2, main = "Cluster 2", col = nrr_cluster_colors[2])
plot.mean.curves.on.circle(nrr_normal_norm_fused_p15_k04_cls_3, main = "Cluster 3", col = nrr_cluster_colors[3])
plot.mean.curves.on.circle(nrr_normal_norm_fused_p15_k04_cls_4, main = "Cluster 4", col = nrr_cluster_colors[4])
par(mfrow=c(1,1))







