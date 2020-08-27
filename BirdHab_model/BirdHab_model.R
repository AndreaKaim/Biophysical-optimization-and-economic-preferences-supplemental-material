##############################################################################################################
#                                                                                                            #
# BIRD HABITAT MODEL                                                                                         #
#                                                                                                            #
# ---------------------------------------------------------------------------------------------------------- #
#                                                                                                            #
# Date: 13-07-2020                                                                                           #
#                                                                                                            #
# Author: Anne Jungandreas                                                                                   #
# Edits (+Willingness to pay function): Andrea Kaim                                                          #
#                                                                                                            #
# Description: Model as used in Kaim et al. "Combining biophysical optimization with economic preference     #
#              analysis for agricultural land use allocation". For a land use map, the model calculates an   #
#              indicator that describes the percental change in suitable habitat for 9 different bird        #
#              species compared to status quo.                                                               #
#                                                                                                            #
# Input data:                                                                                                #
# - Optimization map                                                                                         #
# - See folder "source"                                                                                      #
#                                                                                                            #
# ---------------------------------------------------------------------------------------------------------- #
#                                                                                                            #
# 0. Load Packages                                                                                           #
# 1. Data Preparation                                                                                        #
#    1.1 Read input data                                                                                     #
#    1.2 Convert optimization land use into bird model format                                                #
#    1.3 Identify/Load predictors                                                                            #
#    1.4 Select random raster points                                                                         #
# 2. Species Prediction                                                                                      #
# 3. Output                                                                                                  #
#    3.1 Bird indicator                                                                                      #
#    3.2 Willingness to pay function                                                                         #
#                                                                                                            #
##############################################################################################################


setwd("C:/.../BirdHab_model")
sink("C:/.../BirdHab/console.txt", append=FALSE)


##############################################################################################################
# 0. LOAD PACKAGES
##############################################################################################################

library(sp)
library(raster)
library(randomForest)
library(foreign)
library(rgdal)

##############################################################################################################
# 1. DATA PREPARATION
##############################################################################################################

## 1.1 Read input data ---------------------------------------------------------------------------------------

# read optimization land use map and remove blank lines in ascii header (which cause problems under linux)
ascfile <- file("map.asc","r+")
lines <- readLines(ascfile)
close(ascfile)
lines <- lines[which(lines!="")]
ascfile <- file("map.asc","w+")
writeLines(lines,ascfile)
close(ascfile)

lu_opt <- read.asciigrid("map.asc", as.image = FALSE, plot.image = FALSE, colname = "lu",
                         proj4string = CRS(as.character(NA)))

# read status quo land use and optimization patch ID map
lu_sq <- read.asciigrid("source/maps_michaS/bird_lu_sq.asc", as.image = FALSE, plot.image = FALSE, colname = "bird_lu",
                        proj4string = CRS(as.character(NA)))

patch <- read.asciigrid("source/maps_michaS/hru_fullpatchID_map.asc", as.image = FALSE, plot.image = FALSE, colname = "ID",
                        proj4string = CRS(as.character(NA)))

linE_A_mio <- read.asciigrid("source/maps_michaS/line_a_mio.asc", as.image = FALSE, plot.image = FALSE, colname = "line_a",
                             proj4string = CRS(as.character(NA)))

linE_P_mio <- read.asciigrid("source/maps_michaS/line_p_mio.asc", as.image = FALSE, plot.image = FALSE, colname = "line_p",
                             proj4string = CRS(as.character(NA)))

# read land use translation key and linear element values for each scenario
lu_key <- read.table("source/maps_michaS/scen_lu.txt", sep=";", h=T, as.is=T)


## 1.2 Convert optimization land use into bird model format --------------------------------------------------

lu <- lu_sq
linE_A <- linE_A_mio
linE_P <- linE_P_mio

for(i in 1:dim(lu_key)[1]){
  k <- unique(lu_opt$lu[which(patch$ID==lu_key$OBJECTID_1[i])])
  lu$bird_lu[which(patch$ID==lu_key$OBJECTID_1[i])] <- lu_key[i,k+1]
  #  linE_A$line_a[which(patch$ID==lu_key$OBJECTID_1[i])] <- linE_A_key[i,k+1]
  #  linE_P$line_p[which(patch$ID==lu_key$OBJECTID_1[i])] <- linE_P_key[i,k+1]
}

# save translated land use map
# write.asciigrid(lu, "lu_map_transl.asc", attr = 1, na.value = -2)

## 1.3 Identify/Load predictors ------------------------------------------------------------------------------

r_LU <- raster(lu)
linE_A$line_a <- linE_A$line_a/1000000
linE_P$line_p <- linE_P$line_p/1000000

linE_AntFlHRU <- raster(linE_A)
linE_AntUmfHRU <- raster(linE_P)

# Laubwald
LFor <- r_LU
LFor [!LFor %in% c(NA,1) ] <- 0
LFor <- focal(LFor, w=ceiling(focalWeight(LFor,250,'circle')),sum, na.rm=TRUE) /441/1*100

# Mischwald
MFor <- r_LU
MFor [!MFor %in% c(NA,2) ] <- 0
MFor <- focal(MFor, w=ceiling(focalWeight(MFor,250,'circle')),sum, na.rm=TRUE) /441/2*100

# Nadelwald
NFor <- r_LU
NFor [!NFor %in% c(NA,3) ] <- 0
NFor <- focal(NFor, w=ceiling(focalWeight(NFor,250,'circle')),sum, na.rm=TRUE) /441/3*100

# Forest (all)
For <- r_LU
For [!For %in% c(NA,1,2,3) ] <- 0
For [For %in% c(1,2,3) ] <- 1
For <- focal(For, w=ceiling(focalWeight(For,250,'circle')),sum, na.rm=TRUE) /441/1*100

# Settlement
Set <- r_LU
Set [!Set %in% c(NA,4) ] <- 0
Set <- focal(Set, w=ceiling(focalWeight(Set,250,'circle')),sum, na.rm=TRUE) /441/4*100

# Intensive grassland
IntGL <- r_LU
IntGL [!IntGL %in% c(NA,5) ] <- 0
IntGL <- focal(IntGL, w=ceiling(focalWeight(IntGL,250,'circle')),sum, na.rm=TRUE) /441/5*100

# Extensive grassland
ExtGL <- r_LU
ExtGL [!ExtGL %in% c(NA,8) ] <- 0
ExtGL <- focal(ExtGL, w=ceiling(focalWeight(ExtGL,250,'circle')),sum, na.rm=TRUE) /441/8*100

# Grassland (all)
GL <- r_LU
GL [!GL %in% c(NA,5,8) ] <- 0
GL [GL %in% c(5,8) ] <- 1
GL <- focal(GL, w=ceiling(focalWeight(GL,250,'circle')),sum, na.rm=TRUE) /441/1*100

# Cropland
Ag <- r_LU
Ag [!Ag %in% c(NA,6) ] <- 0
Ag <- focal(Ag, w=ceiling(focalWeight(Ag,250,'circle')),sum, na.rm=TRUE) /441/6*100

# Water bodies
Wat <- r_LU
Wat [!Wat %in% c(NA,7) ] <- 0
Wat <- focal(Wat, w=ceiling(focalWeight(Wat,250,'circle')),sum, na.rm=TRUE) /441/7*100

# Traffic
Ver <- r_LU
Ver [!Ver %in% c(NA,9) ] <- 0
Ver <- focal(Ver, w=ceiling(focalWeight(Ver,250,'circle')),sum, na.rm=TRUE) /441/9*100

# Feuchtgebiete
FeuchtG <- r_LU
FeuchtG [!FeuchtG %in% c(NA,10) ] <- 0
FeuchtG <- focal(FeuchtG, w=ceiling(focalWeight(FeuchtG,250,'circle')),sum, na.rm=TRUE) /441/10*100

# No vegetation
NoVeg <- r_LU
NoVeg [!NoVeg %in% c(NA,11) ] <- 0
NoVeg <- focal(NoVeg, w=ceiling(focalWeight(NoVeg,250,'circle')),sum, na.rm=TRUE) /441/11*100

# Permanent crops
DauerK <- r_LU
DauerK [!DauerK %in% c(NA,12) ] <- 0
DauerK <- focal(DauerK, w=ceiling(focalWeight(DauerK,250,'circle')),sum, na.rm=TRUE) /441/12*100

# Forest edge count
For_edge <- as.matrix(r_LU)
For_edge[For_edge %in% c(4,7,9,11)] <- NA
For_edge[For_edge %in% c(5,6,8,10,12)] <- 0
For_edge[For_edge %in% c(1,2,3)] <- 1
r<-nrow(For_edge)
c<-ncol(For_edge)
y <- matrix(ncol=ncol(For_edge), nrow=nrow(For_edge))
z <- y
for (i in 1:(r-1)){
  for (j in 1:(c)){
    y[i,j]<-(as.numeric((For_edge[i,j]==For_edge[i+1,j])==FALSE))
  }}
for (j in 1:(c-1)){
  for (i in 1:(r)){
    z[i,j]<-(as.numeric((For_edge[i,j]==For_edge[i,j+1])==FALSE))
  }}
For_edge <- y + z
For_edge <- raster (For_edge, template=r_LU)
For_edge <- focal(For_edge, w=ceiling(focalWeight(For_edge,250,'circle')),sum, na.rm=TRUE)
For_edge [For_edge %in% NA] <- 0

# Load and stack remaining predictors
cdist_wasser <- crop(raster ("source/Prediktoren_tifs/cdist_wasser.tif",as.is=TRUE),extent(r_LU))
CDist_Strassen <- crop(raster ("source/Prediktoren_tifs/CDist_Strassen.tif",as.is=TRUE),extent(r_LU))
Boden_AWC <- crop(raster ("source/Prediktoren_tifs/Boden_AWC.tif",as.is=TRUE),extent(r_LU))
Boden_BD <- crop(raster ("source/Prediktoren_tifs/Boden_BD.tif",as.is=TRUE),extent(r_LU))
Boden_CBN <- crop(raster ("source/Prediktoren_tifs/Boden_CBN.tif",as.is=TRUE),extent(r_LU))
Boden_SolK <- crop(raster ("source/Prediktoren_tifs/Boden_SolK.tif",as.is=TRUE),extent(r_LU))
Klima_PCP <- crop(raster ("source/Prediktoren_tifs/Klima_PCP.tif",as.is=TRUE),extent(r_LU))
Klima_Tmean <- crop(raster ("source/Prediktoren_tifs/Klima_Tmean.tif",as.is=TRUE),extent(r_LU))
Klima_Trange <- crop(raster ("source/Prediktoren_tifs/Klima_Trange.tif",as.is=TRUE),extent(r_LU))

predictors <- stack (cdist_wasser,CDist_Strassen,linE_AntFlHRU,linE_AntUmfHRU,
                     Boden_AWC,Boden_BD,Boden_CBN,Boden_SolK,Klima_PCP,Klima_Tmean,Klima_Trange,
                     LFor, MFor, NFor, For, Set, IntGL, Ag, Wat, ExtGL, GL, Ver, FeuchtG, NoVeg, 
                     DauerK,For_edge)

names (predictors) = c("cdist_wasser", "CDist_Strassen", 
                       "linE_AntFlHRU", "linE_AntUmfHRU", "Boden_AWC", "Boden_BD",
                       "Boden_CBN", "Boden_SolK", "Klima_PCP", "Klima_Tmean", "Klima_Trange", "LFor", 
                       "MFor", "NFor", "For", "Set", "IntGL", "Ag", "Wat", "ExtGL", "GL", "Ver", "FeuchtG",
                       "NoVeg", "DauerK", "For_edge")

rm (cdist_wasser,CDist_Strassen,linE_AntFlHRU,linE_AntUmfHRU,
    Boden_AWC,Boden_BD,Boden_CBN,Boden_SolK,Klima_PCP,Klima_Tmean,Klima_Trange,
    LFor, MFor, NFor, For, Set, IntGL, Ag, Wat, ExtGL, GL, Ver, FeuchtG, NoVeg, 
    DauerK,For_edge)


## 1.4 Select random raster points ---------------------------------------------------------------------------

# 40000 raster points (minimum: 1 per HRU)
samp_size = 40000

set.seed(25000)
temp_ras <- sampleStratified (raster(patch),1,xy=TRUE, sp=TRUE, na.rm=TRUE)
coord4Bio <- coordinates(temp_ras)
temp_ras <- sampleRandom (raster(patch),samp_size,xy=TRUE, sp=TRUE, na.rm=TRUE)
coord4Bio <- rbind (coord4Bio,coordinates(temp_ras))
coord4Bio <- coord4Bio [-which(duplicated(coord4Bio)==TRUE),]
colnames (coord4Bio) <- c("X","Y")
pred4Bio <- extract (predictors, coord4Bio)
for (i in 1:ncol(pred4Bio)) { # i=24
  pred4Bio <- pred4Bio[!is.na(pred4Bio[,i]),]
}

##############################################################################################################
# 2. SPECIES PREDICTION
##############################################################################################################

# Load reference values
load("source/refvalues2.RData")


# Whinchat (Saxicola rubetra)

# load thresholds
load ("source/Braunkehlchen_BGV_UK250/tresholds.RData")


load ("source/Braunkehlchen_BGV_UK250/RF1.RData") # load Random Forest model results
pred_RF1 <- predict (theRF,newdata=pred4Bio)
pred_RF1 [pred_RF1 > tresholds[1]] <- 1; pred_RF1 [pred_RF1 < 1] <- 0
load ("source/Braunkehlchen_BGV_UK250/RF2.RData")
pred_RF2 <- predict (theRF,newdata=pred4Bio)
pred_RF2 [pred_RF2 > tresholds[2]] <- 1; pred_RF2 [pred_RF2 < 1] <- 0
load ("source/Braunkehlchen_BGV_UK250/RF3.RData")
pred_RF3 <- predict (theRF,newdata=pred4Bio)
pred_RF3 [pred_RF3 > tresholds[3]] <- 1; pred_RF3 [pred_RF3 < 1] <- 0
load ("source/Braunkehlchen_BGV_UK250/RF4.RData")
pred_RF4 <- predict (theRF,newdata=pred4Bio)
pred_RF4 [pred_RF4 > tresholds[4]] <- 1; pred_RF4 [pred_RF4 < 1] <- 0
load ("source/Braunkehlchen_BGV_UK250/RF5.RData")
pred_RF5 <- predict (theRF,newdata=pred4Bio)
pred_RF5 [pred_RF5 > tresholds[5]] <- 1; pred_RF5 [pred_RF5 < 1] <- 0
load ("source/Braunkehlchen_BGV_UK250/RF6.RData")
pred_RF6 <- predict (theRF,newdata=pred4Bio)
pred_RF6 [pred_RF6 > tresholds[6]] <- 1; pred_RF6 [pred_RF6 < 1] <- 0
load ("source/Braunkehlchen_BGV_UK250/RF7.RData")
pred_RF7 <- predict (theRF,newdata=pred4Bio)
pred_RF7 [pred_RF7 > tresholds[7]] <- 1; pred_RF7 [pred_RF7 < 1] <- 0
load ("source/Braunkehlchen_BGV_UK250/RF8.RData")
pred_RF8 <- predict (theRF,newdata=pred4Bio)
pred_RF8 [pred_RF8 > tresholds[8]] <- 1; pred_RF8 [pred_RF8 < 1] <- 0
load ("source/Braunkehlchen_BGV_UK250/RF9.RData")
pred_RF9 <- predict (theRF,newdata=pred4Bio)
pred_RF9 [pred_RF9 > tresholds[9]] <- 1; pred_RF9 [pred_RF9 < 1] <- 0
load ("source/Braunkehlchen_BGV_UK250/RF10.RData")
pred_RF10 <- predict (theRF,newdata=pred4Bio)
pred_RF10 [pred_RF10 > tresholds[10]] <- 1; pred_RF10 [pred_RF10 < 1] <- 0
pred_sumRF <- pred_RF1+pred_RF2+pred_RF3+pred_RF4+pred_RF5+pred_RF6+pred_RF7+pred_RF8+pred_RF9+pred_RF10
pred_sumRF [pred_sumRF < 5] <- 0; pred_sumRF [pred_sumRF > 4] <- 1
occ_Braunkehlchen <- (sum (pred_sumRF) / length (pred_sumRF)) / ref_Braunkehlchen

# Western Jackdaw (Coloeus monedula)

load ("source/Dohle_BGV_UK250/tresholds.RData")

load ("source/Dohle_BGV_UK250/RF1.RData")
pred_RF1 <- predict (theRF,newdata=pred4Bio)
pred_RF1 [pred_RF1 > tresholds[1]] <- 1; pred_RF1 [pred_RF1 < 1] <- 0
load ("source/Dohle_BGV_UK250/RF2.RData")
pred_RF2 <- predict (theRF,newdata=pred4Bio)
pred_RF2 [pred_RF2 > tresholds[2]] <- 1; pred_RF2 [pred_RF2 < 1] <- 0
load ("source/Dohle_BGV_UK250/RF3.RData")
pred_RF3 <- predict (theRF,newdata=pred4Bio)
pred_RF3 [pred_RF3 > tresholds[3]] <- 1; pred_RF3 [pred_RF3 < 1] <- 0
load ("source/Dohle_BGV_UK250/RF4.RData")
pred_RF4 <- predict (theRF,newdata=pred4Bio)
pred_RF4 [pred_RF4 > tresholds[4]] <- 1; pred_RF4 [pred_RF4 < 1] <- 0
load ("source/Dohle_BGV_UK250/RF5.RData")
pred_RF5 <- predict (theRF,newdata=pred4Bio)
pred_RF5 [pred_RF5 > tresholds[5]] <- 1; pred_RF5 [pred_RF5 < 1] <- 0
load ("source/Dohle_BGV_UK250/RF6.RData")
pred_RF6 <- predict (theRF,newdata=pred4Bio)
pred_RF6 [pred_RF6 > tresholds[6]] <- 1; pred_RF6 [pred_RF6 < 1] <- 0
load ("source/Dohle_BGV_UK250/RF7.RData")
pred_RF7 <- predict (theRF,newdata=pred4Bio)
pred_RF7 [pred_RF7 > tresholds[7]] <- 1; pred_RF7 [pred_RF7 < 1] <- 0
load ("source/Dohle_BGV_UK250/RF8.RData")
pred_RF8 <- predict (theRF,newdata=pred4Bio)
pred_RF8 [pred_RF8 > tresholds[8]] <- 1; pred_RF8 [pred_RF8 < 1] <- 0
load ("source/Dohle_BGV_UK250/RF9.RData")
pred_RF9 <- predict (theRF,newdata=pred4Bio)
pred_RF9 [pred_RF9 > tresholds[9]] <- 1; pred_RF9 [pred_RF9 < 1] <- 0
load ("source/Dohle_BGV_UK250/RF10.RData")
pred_RF10 <- predict (theRF,newdata=pred4Bio)
pred_RF10 [pred_RF10 > tresholds[10]] <- 1; pred_RF10 [pred_RF10 < 1] <- 0
pred_sumRF <- pred_RF1+pred_RF2+pred_RF3+pred_RF4+pred_RF5+pred_RF6+pred_RF7+pred_RF8+pred_RF9+pred_RF10
pred_sumRF [pred_sumRF < 5] <- 0; pred_sumRF [pred_sumRF > 4] <- 1

occ_Dohle <- (sum (pred_sumRF) / length (pred_sumRF)) / ref_Dohle

# Common Kingfisher (Alcedo atthis)

load ("source/Eisvogel_BGV_UK250/tresholds.RData")

load ("source/Eisvogel_BGV_UK250/RF1.RData")
pred_RF1 <- predict (theRF,newdata=pred4Bio)
pred_RF1 [pred_RF1 > tresholds[1]] <- 1; pred_RF1 [pred_RF1 < 1] <- 0
load ("source/Eisvogel_BGV_UK250/RF2.RData")
pred_RF2 <- predict (theRF,newdata=pred4Bio)
pred_RF2 [pred_RF2 > tresholds[2]] <- 1; pred_RF2 [pred_RF2 < 1] <- 0
load ("source/Eisvogel_BGV_UK250/RF3.RData")
pred_RF3 <- predict (theRF,newdata=pred4Bio)
pred_RF3 [pred_RF3 > tresholds[3]] <- 1; pred_RF3 [pred_RF3 < 1] <- 0
load ("source/Eisvogel_BGV_UK250/RF4.RData")
pred_RF4 <- predict (theRF,newdata=pred4Bio)
pred_RF4 [pred_RF4 > tresholds[4]] <- 1; pred_RF4 [pred_RF4 < 1] <- 0
load ("source/Eisvogel_BGV_UK250/RF5.RData")
pred_RF5 <- predict (theRF,newdata=pred4Bio)
pred_RF5 [pred_RF5 > tresholds[5]] <- 1; pred_RF5 [pred_RF5 < 1] <- 0
load ("source/Eisvogel_BGV_UK250/RF6.RData")
pred_RF6 <- predict (theRF,newdata=pred4Bio)
pred_RF6 [pred_RF6 > tresholds[6]] <- 1; pred_RF6 [pred_RF6 < 1] <- 0
load ("source/Eisvogel_BGV_UK250/RF7.RData")
pred_RF7 <- predict (theRF,newdata=pred4Bio)
pred_RF7 [pred_RF7 > tresholds[7]] <- 1; pred_RF7 [pred_RF7 < 1] <- 0
load ("source/Eisvogel_BGV_UK250/RF8.RData")
pred_RF8 <- predict (theRF,newdata=pred4Bio)
pred_RF8 [pred_RF8 > tresholds[8]] <- 1; pred_RF8 [pred_RF8 < 1] <- 0
load ("source/Eisvogel_BGV_UK250/RF9.RData")
pred_RF9 <- predict (theRF,newdata=pred4Bio)
pred_RF9 [pred_RF9 > tresholds[9]] <- 1; pred_RF9 [pred_RF9 < 1] <- 0
load ("source/Eisvogel_BGV_UK250/RF10.RData")
pred_RF10 <- predict (theRF,newdata=pred4Bio)
pred_RF10 [pred_RF10 > tresholds[10]] <- 1; pred_RF10 [pred_RF10 < 1] <- 0
pred_sumRF <- pred_RF1+pred_RF2+pred_RF3+pred_RF4+pred_RF5+pred_RF6+pred_RF7+pred_RF8+pred_RF9+pred_RF10
pred_sumRF [pred_sumRF < 5] <- 0; pred_sumRF [pred_sumRF > 4] <- 1

occ_Eisvogel <- (sum (pred_sumRF) / length (pred_sumRF)) / ref_Eisvogel

# Common Redstart (Phoenicurus phoenicurus)

load ("source/Gartenrotschwanz_BGV_UK250/tresholds.RData")

load ("source/Gartenrotschwanz_BGV_UK250/RF1.RData")
pred_RF1 <- predict (theRF,newdata=pred4Bio)
pred_RF1 [pred_RF1 > tresholds[1]] <- 1; pred_RF1 [pred_RF1 < 1] <- 0
load ("source/Gartenrotschwanz_BGV_UK250/RF2.RData")
pred_RF2 <- predict (theRF,newdata=pred4Bio)
pred_RF2 [pred_RF2 > tresholds[2]] <- 1; pred_RF2 [pred_RF2 < 1] <- 0
load ("source/Gartenrotschwanz_BGV_UK250/RF3.RData")
pred_RF3 <- predict (theRF,newdata=pred4Bio)
pred_RF3 [pred_RF3 > tresholds[3]] <- 1; pred_RF3 [pred_RF3 < 1] <- 0
load ("source/Gartenrotschwanz_BGV_UK250/RF4.RData")
pred_RF4 <- predict (theRF,newdata=pred4Bio)
pred_RF4 [pred_RF4 > tresholds[4]] <- 1; pred_RF4 [pred_RF4 < 1] <- 0
load ("source/Gartenrotschwanz_BGV_UK250/RF5.RData")
pred_RF5 <- predict (theRF,newdata=pred4Bio)
pred_RF5 [pred_RF5 > tresholds[5]] <- 1; pred_RF5 [pred_RF5 < 1] <- 0
load ("source/Gartenrotschwanz_BGV_UK250/RF6.RData")
pred_RF6 <- predict (theRF,newdata=pred4Bio)
pred_RF6 [pred_RF6 > tresholds[6]] <- 1; pred_RF6 [pred_RF6 < 1] <- 0
load ("source/Gartenrotschwanz_BGV_UK250/RF7.RData")
pred_RF7 <- predict (theRF,newdata=pred4Bio)
pred_RF7 [pred_RF7 > tresholds[7]] <- 1; pred_RF7 [pred_RF7 < 1] <- 0
load ("source/Gartenrotschwanz_BGV_UK250/RF8.RData")
pred_RF8 <- predict (theRF,newdata=pred4Bio)
pred_RF8 [pred_RF8 > tresholds[8]] <- 1; pred_RF8 [pred_RF8 < 1] <- 0
load ("source/Gartenrotschwanz_BGV_UK250/RF9.RData")
pred_RF9 <- predict (theRF,newdata=pred4Bio)
pred_RF9 [pred_RF9 > tresholds[9]] <- 1; pred_RF9 [pred_RF9 < 1] <- 0
load ("source/Gartenrotschwanz_BGV_UK250/RF10.RData")
pred_RF10 <- predict (theRF,newdata=pred4Bio)
pred_RF10 [pred_RF10 > tresholds[10]] <- 1; pred_RF10 [pred_RF10 < 1] <- 0
pred_sumRF <- pred_RF1+pred_RF2+pred_RF3+pred_RF4+pred_RF5+pred_RF6+pred_RF7+pred_RF8+pred_RF9+pred_RF10
pred_sumRF [pred_sumRF < 5] <- 0; pred_sumRF [pred_sumRF > 4] <- 1

occ_Gartenrotschwanz <- (sum (pred_sumRF) / length (pred_sumRF)) / ref_Gartenrotschwanz

# Wood lark (Lullula arborea)

load ("source/Haubenlerche_BGV_UK250/tresholds.RData")

load ("source/Haubenlerche_BGV_UK250/RF1.RData")
pred_RF1 <- predict (theRF,newdata=pred4Bio)
pred_RF1 [pred_RF1 > tresholds[1]] <- 1; pred_RF1 [pred_RF1 < 1] <- 0
load ("source/Haubenlerche_BGV_UK250/RF2.RData")
pred_RF2 <- predict (theRF,newdata=pred4Bio)
pred_RF2 [pred_RF2 > tresholds[2]] <- 1; pred_RF2 [pred_RF2 < 1] <- 0
load ("source/Haubenlerche_BGV_UK250/RF3.RData")
pred_RF3 <- predict (theRF,newdata=pred4Bio)
pred_RF3 [pred_RF3 > tresholds[3]] <- 1; pred_RF3 [pred_RF3 < 1] <- 0
load ("source/Haubenlerche_BGV_UK250/RF4.RData")
pred_RF4 <- predict (theRF,newdata=pred4Bio)
pred_RF4 [pred_RF4 > tresholds[4]] <- 1; pred_RF4 [pred_RF4 < 1] <- 0
load ("source/Haubenlerche_BGV_UK250/RF5.RData")
pred_RF5 <- predict (theRF,newdata=pred4Bio)
pred_RF5 [pred_RF5 > tresholds[5]] <- 1; pred_RF5 [pred_RF5 < 1] <- 0
load ("source/Haubenlerche_BGV_UK250/RF6.RData")
pred_RF6 <- predict (theRF,newdata=pred4Bio)
pred_RF6 [pred_RF6 > tresholds[6]] <- 1; pred_RF6 [pred_RF6 < 1] <- 0
load ("source/Haubenlerche_BGV_UK250/RF7.RData")
pred_RF7 <- predict (theRF,newdata=pred4Bio)
pred_RF7 [pred_RF7 > tresholds[7]] <- 1; pred_RF7 [pred_RF7 < 1] <- 0
load ("source/Haubenlerche_BGV_UK250/RF8.RData")
pred_RF8 <- predict (theRF,newdata=pred4Bio)
pred_RF8 [pred_RF8 > tresholds[8]] <- 1; pred_RF8 [pred_RF8 < 1] <- 0
load ("source/Haubenlerche_BGV_UK250/RF9.RData")
pred_RF9 <- predict (theRF,newdata=pred4Bio)
pred_RF9 [pred_RF9 > tresholds[9]] <- 1; pred_RF9 [pred_RF9 < 1] <- 0
load ("source/Haubenlerche_BGV_UK250/RF10.RData")
pred_RF10 <- predict (theRF,newdata=pred4Bio)
pred_RF10 [pred_RF10 > tresholds[10]] <- 1; pred_RF10 [pred_RF10 < 1] <- 0
pred_sumRF <- pred_RF1+pred_RF2+pred_RF3+pred_RF4+pred_RF5+pred_RF6+pred_RF7+pred_RF8+pred_RF9+pred_RF10
pred_sumRF [pred_sumRF < 5] <- 0; pred_sumRF [pred_sumRF > 4] <- 1

occ_Haubenlerche <- (sum (pred_sumRF) / length (pred_sumRF)) / ref_Haubenlerche

# Crested Lark (Galerida cristata)

load ("source/Heidelerche_BGV_UK250/tresholds.RData")

load ("source/Heidelerche_BGV_UK250/RF1.RData")
pred_RF1 <- predict (theRF,newdata=pred4Bio)
pred_RF1 [pred_RF1 > tresholds[1]] <- 1; pred_RF1 [pred_RF1 < 1] <- 0
load ("source/Heidelerche_BGV_UK250/RF2.RData")
pred_RF2 <- predict (theRF,newdata=pred4Bio)
pred_RF2 [pred_RF2 > tresholds[2]] <- 1; pred_RF2 [pred_RF2 < 1] <- 0
load ("source/Heidelerche_BGV_UK250/RF3.RData")
pred_RF3 <- predict (theRF,newdata=pred4Bio)
pred_RF3 [pred_RF3 > tresholds[3]] <- 1; pred_RF3 [pred_RF3 < 1] <- 0
load ("source/Heidelerche_BGV_UK250/RF4.RData")
pred_RF4 <- predict (theRF,newdata=pred4Bio)
pred_RF4 [pred_RF4 > tresholds[4]] <- 1; pred_RF4 [pred_RF4 < 1] <- 0
load ("source/Heidelerche_BGV_UK250/RF5.RData")
pred_RF5 <- predict (theRF,newdata=pred4Bio)
pred_RF5 [pred_RF5 > tresholds[5]] <- 1; pred_RF5 [pred_RF5 < 1] <- 0
load ("source/Heidelerche_BGV_UK250/RF6.RData")
pred_RF6 <- predict (theRF,newdata=pred4Bio)
pred_RF6 [pred_RF6 > tresholds[6]] <- 1; pred_RF6 [pred_RF6 < 1] <- 0
load ("source/Heidelerche_BGV_UK250/RF7.RData")
pred_RF7 <- predict (theRF,newdata=pred4Bio)
pred_RF7 [pred_RF7 > tresholds[7]] <- 1; pred_RF7 [pred_RF7 < 1] <- 0
load ("source/Heidelerche_BGV_UK250/RF8.RData")
pred_RF8 <- predict (theRF,newdata=pred4Bio)
pred_RF8 [pred_RF8 > tresholds[8]] <- 1; pred_RF8 [pred_RF8 < 1] <- 0
load ("source/Heidelerche_BGV_UK250/RF9.RData")
pred_RF9 <- predict (theRF,newdata=pred4Bio)
pred_RF9 [pred_RF9 > tresholds[9]] <- 1; pred_RF9 [pred_RF9 < 1] <- 0
load ("source/Heidelerche_BGV_UK250/RF10.RData")
pred_RF10 <- predict (theRF,newdata=pred4Bio)
pred_RF10 [pred_RF10 > tresholds[10]] <- 1; pred_RF10 [pred_RF10 < 1] <- 0
pred_sumRF <- pred_RF1+pred_RF2+pred_RF3+pred_RF4+pred_RF5+pred_RF6+pred_RF7+pred_RF8+pred_RF9+pred_RF10
pred_sumRF [pred_sumRF < 5] <- 0; pred_sumRF [pred_sumRF > 4] <- 1

occ_Heidelerche <- (sum (pred_sumRF) / length (pred_sumRF)) / ref_Heidelerche

# Northern Lapwing (Venellus vanellus)

load ("source/Kiebitz_BGV_UK250/tresholds.RData")

load ("source/Kiebitz_BGV_UK250/RF1.RData")
pred_RF1 <- predict (theRF,newdata=pred4Bio)
pred_RF1 [pred_RF1 > tresholds[1]] <- 1; pred_RF1 [pred_RF1 < 1] <- 0
load ("source/Kiebitz_BGV_UK250/RF2.RData")
pred_RF2 <- predict (theRF,newdata=pred4Bio)
pred_RF2 [pred_RF2 > tresholds[2]] <- 1; pred_RF2 [pred_RF2 < 1] <- 0
load ("source/Kiebitz_BGV_UK250/RF3.RData")
pred_RF3 <- predict (theRF,newdata=pred4Bio)
pred_RF3 [pred_RF3 > tresholds[3]] <- 1; pred_RF3 [pred_RF3 < 1] <- 0
load ("source/Kiebitz_BGV_UK250/RF4.RData")
pred_RF4 <- predict (theRF,newdata=pred4Bio)
pred_RF4 [pred_RF4 > tresholds[4]] <- 1; pred_RF4 [pred_RF4 < 1] <- 0
load ("source/Kiebitz_BGV_UK250/RF5.RData")
pred_RF5 <- predict (theRF,newdata=pred4Bio)
pred_RF5 [pred_RF5 > tresholds[5]] <- 1; pred_RF5 [pred_RF5 < 1] <- 0
load ("source/Kiebitz_BGV_UK250/RF6.RData")
pred_RF6 <- predict (theRF,newdata=pred4Bio)
pred_RF6 [pred_RF6 > tresholds[6]] <- 1; pred_RF6 [pred_RF6 < 1] <- 0
load ("source/Kiebitz_BGV_UK250/RF7.RData")
pred_RF7 <- predict (theRF,newdata=pred4Bio)
pred_RF7 [pred_RF7 > tresholds[7]] <- 1; pred_RF7 [pred_RF7 < 1] <- 0
load ("source/Kiebitz_BGV_UK250/RF8.RData")
pred_RF8 <- predict (theRF,newdata=pred4Bio)
pred_RF8 [pred_RF8 > tresholds[8]] <- 1; pred_RF8 [pred_RF8 < 1] <- 0
load ("source/Kiebitz_BGV_UK250/RF9.RData")
pred_RF9 <- predict (theRF,newdata=pred4Bio)
pred_RF9 [pred_RF9 > tresholds[9]] <- 1; pred_RF9 [pred_RF9 < 1] <- 0
load ("source/Kiebitz_BGV_UK250/RF10.RData")
pred_RF10 <- predict (theRF,newdata=pred4Bio)
pred_RF10 [pred_RF10 > tresholds[10]] <- 1; pred_RF10 [pred_RF10 < 1] <- 0
pred_sumRF <- pred_RF1+pred_RF2+pred_RF3+pred_RF4+pred_RF5+pred_RF6+pred_RF7+pred_RF8+pred_RF9+pred_RF10
pred_sumRF [pred_sumRF < 5] <- 0; pred_sumRF [pred_sumRF > 4] <- 1

occ_Kiebitz <- (sum (pred_sumRF) / length (pred_sumRF)) / ref_Kiebitz

# Barn Owl (Tyto alba)

load ("source/Schleiereule_BGV_UK250/tresholds.RData")

load ("source/Schleiereule_BGV_UK250/RF1.RData")
pred_RF1 <- predict (theRF,newdata=pred4Bio)
pred_RF1 [pred_RF1 > tresholds[1]] <- 1; pred_RF1 [pred_RF1 < 1] <- 0
load ("source/Schleiereule_BGV_UK250/RF2.RData")
pred_RF2 <- predict (theRF,newdata=pred4Bio)
pred_RF2 [pred_RF2 > tresholds[2]] <- 1; pred_RF2 [pred_RF2 < 1] <- 0
load ("source/Schleiereule_BGV_UK250/RF3.RData")
pred_RF3 <- predict (theRF,newdata=pred4Bio)
pred_RF3 [pred_RF3 > tresholds[3]] <- 1; pred_RF3 [pred_RF3 < 1] <- 0
load ("source/Schleiereule_BGV_UK250/RF4.RData")
pred_RF4 <- predict (theRF,newdata=pred4Bio)
pred_RF4 [pred_RF4 > tresholds[4]] <- 1; pred_RF4 [pred_RF4 < 1] <- 0
load ("source/Schleiereule_BGV_UK250/RF5.RData")
pred_RF5 <- predict (theRF,newdata=pred4Bio)
pred_RF5 [pred_RF5 > tresholds[5]] <- 1; pred_RF5 [pred_RF5 < 1] <- 0
load ("source/Schleiereule_BGV_UK250/RF6.RData")
pred_RF6 <- predict (theRF,newdata=pred4Bio)
pred_RF6 [pred_RF6 > tresholds[6]] <- 1; pred_RF6 [pred_RF6 < 1] <- 0
load ("source/Schleiereule_BGV_UK250/RF7.RData")
pred_RF7 <- predict (theRF,newdata=pred4Bio)
pred_RF7 [pred_RF7 > tresholds[7]] <- 1; pred_RF7 [pred_RF7 < 1] <- 0
load ("source/Schleiereule_BGV_UK250/RF8.RData")
pred_RF8 <- predict (theRF,newdata=pred4Bio)
pred_RF8 [pred_RF8 > tresholds[8]] <- 1; pred_RF8 [pred_RF8 < 1] <- 0
load ("source/Schleiereule_BGV_UK250/RF9.RData")
pred_RF9 <- predict (theRF,newdata=pred4Bio)
pred_RF9 [pred_RF9 > tresholds[9]] <- 1; pred_RF9 [pred_RF9 < 1] <- 0
load ("source/Schleiereule_BGV_UK250/RF10.RData")
pred_RF10 <- predict (theRF,newdata=pred4Bio)
pred_RF10 [pred_RF10 > tresholds[10]] <- 1; pred_RF10 [pred_RF10 < 1] <- 0
pred_sumRF <- pred_RF1+pred_RF2+pred_RF3+pred_RF4+pred_RF5+pred_RF6+pred_RF7+pred_RF8+pred_RF9+pred_RF10
pred_sumRF [pred_sumRF < 5] <- 0; pred_sumRF [pred_sumRF > 4] <- 1

occ_Schleiereule <- (sum (pred_sumRF) / length (pred_sumRF)) / ref_Schleiereule

# Northern Wheatear (Oenanthe oenanthe)

load ("source/Steinschmaetzer_BGV_UK250/tresholds.RData")

load ("source/Steinschmaetzer_BGV_UK250/RF1.RData")
pred_RF1 <- predict (theRF,newdata=pred4Bio)
pred_RF1 [pred_RF1 > tresholds[1]] <- 1; pred_RF1 [pred_RF1 < 1] <- 0
load ("source/Steinschmaetzer_BGV_UK250/RF2.RData")
pred_RF2 <- predict (theRF,newdata=pred4Bio)
pred_RF2 [pred_RF2 > tresholds[2]] <- 1; pred_RF2 [pred_RF2 < 1] <- 0
load ("source/Steinschmaetzer_BGV_UK250/RF3.RData")
pred_RF3 <- predict (theRF,newdata=pred4Bio)
pred_RF3 [pred_RF3 > tresholds[3]] <- 1; pred_RF3 [pred_RF3 < 1] <- 0
load ("source/Steinschmaetzer_BGV_UK250/RF4.RData")
pred_RF4 <- predict (theRF,newdata=pred4Bio)
pred_RF4 [pred_RF4 > tresholds[4]] <- 1; pred_RF4 [pred_RF4 < 1] <- 0
load ("source/Steinschmaetzer_BGV_UK250/RF5.RData")
pred_RF5 <- predict (theRF,newdata=pred4Bio)
pred_RF5 [pred_RF5 > tresholds[5]] <- 1; pred_RF5 [pred_RF5 < 1] <- 0
load ("source/Steinschmaetzer_BGV_UK250/RF6.RData")
pred_RF6 <- predict (theRF,newdata=pred4Bio)
pred_RF6 [pred_RF6 > tresholds[6]] <- 1; pred_RF6 [pred_RF6 < 1] <- 0
load ("source/Steinschmaetzer_BGV_UK250/RF7.RData")
pred_RF7 <- predict (theRF,newdata=pred4Bio)
pred_RF7 [pred_RF7 > tresholds[7]] <- 1; pred_RF7 [pred_RF7 < 1] <- 0
load ("source/Steinschmaetzer_BGV_UK250/RF8.RData")
pred_RF8 <- predict (theRF,newdata=pred4Bio)
pred_RF8 [pred_RF8 > tresholds[8]] <- 1; pred_RF8 [pred_RF8 < 1] <- 0
load ("source/Steinschmaetzer_BGV_UK250/RF9.RData")
pred_RF9 <- predict (theRF,newdata=pred4Bio)
pred_RF9 [pred_RF9 > tresholds[9]] <- 1; pred_RF9 [pred_RF9 < 1] <- 0
load ("source/Steinschmaetzer_BGV_UK250/RF10.RData")
pred_RF10 <- predict (theRF,newdata=pred4Bio)
pred_RF10 [pred_RF10 > tresholds[10]] <- 1; pred_RF10 [pred_RF10 < 1] <- 0
pred_sumRF <- pred_RF1+pred_RF2+pred_RF3+pred_RF4+pred_RF5+pred_RF6+pred_RF7+pred_RF8+pred_RF9+pred_RF10
pred_sumRF [pred_sumRF < 5] <- 0; pred_sumRF [pred_sumRF > 4] <- 1

occ_Steinschmaetzer <- (sum (pred_sumRF) / length (pred_sumRF)) / ref_Steinschmaetzer


##############################################################################################################
# 3.  Output
##############################################################################################################

## 3.1 Bird indicator ----------------------------------------------------------------------------------------

# bird indicator (percental change from status quo)
occ_mean <- mean(c(occ_Braunkehlchen,occ_Dohle,occ_Eisvogel,occ_Gartenrotschwanz,occ_Haubenlerche,
                   occ_Heidelerche,occ_Kiebitz,occ_Schleiereule,occ_Steinschmaetzer))


## 3.2 Willingness to pay function ----------------------------------------------------------------------------

# piecewise: for occ_mean < 1 => wtp = -1, for occ_mean < 1.5999999 => wtp = 213.36 (max wtp), else wtp follows polynom of 6th degree

if (occ_mean < 1){
  wtp <- -1
} else if(occ_mean > 1.5999999){
  wtp <- 213.36
} else {
  wtp <- -5417*occ_mean^6 + 51621*occ_mean^5 - 199662*occ_mean^4 + 400856*occ_mean^3 - 440436*occ_mean^2 + 251474*occ_mean - 58436
}

write.table(round(wtp, 2), "BirdHab_output.csv",append=FALSE ,sep =";",col.names=FALSE ,row.names=FALSE)



sink()
