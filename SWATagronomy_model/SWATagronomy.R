##############################################################################################################
#                                                                                                            #
# SWAT AGRONOMY MODEL                                                                                        #
#                                                                                                            #
# ---------------------------------------------------------------------------------------------------------- #
#                                                                                                            #
# Date: 15-07-2020                                                                                           #
#                                                                                                            #
# Author: Michael Strauch                                                                                    #
# Edits: Andrea Kaim                                                                                         #
#                                                                                                            #
# Description: Script as used in Kaim et al. "Combining biophysical optimization with economic preference    #
#              analysis for agricultural land use allocation". The model calculates the agricultural gross   #
#              margin based on SWAT simulation results. SWAT results must be given on HRU-level for each     #
#              possible land use/cover class.                                                                #
#              Considered here: (1) cropland represented by representative crop rotations                    #
#                               (2) intensively used grassland (more fertilizer and cuttings)                #
#                               (3) extensively used grassland (less fertilizer and cuttings)                #
#              SWAT simulation results refer to average crop yields simulated for the period from 1996 to    #
#              2008 combined with crop-specific prices and costs from "Leistungsrechnung der KTBL"           #                                                    #
#                                                                                                            #
# Input data:                                                                                                #
# - HRU list                                                                                                 #
# - Genome of optimized land use map                                                                         #
# - Gross margins                                                                                            #
#                                                                                                            #
# ---------------------------------------------------------------------------------------------------------- #
#                                                                                                            #
# 1. Data Preparation                                                                                        #
#    1.1 Read input data                                                                                     #
#    1.2 Combine data                                                                                        #
# 2. Calculate gross margins                                                                                 #
# 3. Output                                                                                                  #
#                                                                                                            #
##############################################################################################################


setwd("C:/.../SWATagronomy_model")
sink("C:/.../SWATagronomy/console.txt", append=FALSE)


##############################################################################################################
# 1. DATA PREPARATION
##############################################################################################################

## 1.1 Read input data ---------------------------------------------------------------------------------------

HRUs <- read.table("hru_list.txt")
genome <- read.csv("genom.csv")
margin <- read.table("margin_lossa.txt", sep=";",h=T)


## 1.2 Combine data ------------------------------------------------------------------------------------------

# 1 = status quo (cropland or grassland),
# 2 = intensive grassland,
# 3 = extensive grassland,
# else = forest

# Combine HRUs and genome
HRUs.genome <- cbind.data.frame(HRUs,genome)

# Merge margin and genome by HRU
margin.optim <- merge.data.frame(HRUs.genome, margin, by.x="V1", by.y="HRUGIS")


##############################################################################################################
# 2. CALCULATE GROSS MARGINS
##############################################################################################################

margin.sum = 0

# Calculate margin for first land cover in genome.list
genome.list <- unique(genome$genom)
subset1 <- subset(margin.optim, margin.optim$genom==genome.list[1])
if(genome.list[1]==1){
  margin.sum <- sum(subset1$margintotal_SQ)
} else{
  if(genome.list[1]==2){
    margin.sum <- sum(subset1$margintotal_int)
  } else{
    if(genome.list[1]==3){
    margin.sum <- sum(subset1$margintotal_ext)
    } else{
    margin.sum <- sum(subset1$margintotal_forest)
    }
  }
}

# Calculate margin for all other land covers in genome.list and derive total sum
if(length(genome.list)>1){
  for(i in 2:length(genome.list)){
    subsetx <- subset(margin.optim, margin.optim$genom==genome.list[i])
    if(genome.list[i]==1){
      margin.sumx <- sum(subsetx$margintotal_SQ)
    } else{
      if(genome.list[i]==2){
        margin.sumx <- sum(subsetx$margintotal_int)
      } else{
        if(genome.list[i]==3){
          margin.sumx <- sum(subsetx$margintotal_ext)
        } else{
          margin.sumx <- sum(subsetx$margintotal_forest)
        }
      }
    }
    margin.sum <- sum(margin.sum, margin.sumx)
  }
}

# Add margin from HRUs not considered in optimization
margin.sum <- margin.sum + 159136.3

margin.sum.ha <- margin.sum/14076.3125 # Total area of Lossa River Basin: 14076.3125 ha


##############################################################################################################
# 3. OUTPUT
##############################################################################################################

# Write model output
write.table(round(margin.sum.ha,2) , "SWATagronomy_output.csv",append=FALSE ,sep =";",col.names=FALSE ,row.names=FALSE)



sink()
