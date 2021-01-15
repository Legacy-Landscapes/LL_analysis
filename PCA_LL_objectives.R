#-#-# Clear memory #-#-#
rm(list=ls())

library(devtools) 
library(ggbiplot) 
library(factoextra)
library(easystats)
library(tidyverse)
library(ggplot2)
library(ggfortify)


#-#-# Get the final data frame (gridded PAs with values for all indicators) #-#-#
setwd("/Users/alkevoskamp/Documents/Legacy Landscapes/Legacy_landscapes_analysis/LL_analysis/Processed_data/")
PAdata <- read.csv("Final_dataset_IUCN_WHS_KBA_18_12.csv")


#-#-# Subset to columns that are included in the weighting #-#-#
PAdataS <- PAdata[c("SP_ID","Int_Name","PA_type","PA_DEF","IUCN_CAT","RealmNr","Area_Rep","Area_GIS","Prop_SR_B","Prop_SR_M","Prop_SR_A","Prop_SR_R",
                    "RSR_B","RSR_M","RSR_A","RSR_R","PE_WMean_B","PE_WMean_M","PE_WMean_A","PE_WMean_R","TO_PA_weightbird","TO_PA_weightmammal",
                    "TO_PA_weightamphibian","TO_PA_weightreptile","HFP_WMean","BaseCarbon_WMean","VulCarbon_WMean","IrrCarbon_WMean","BII_WMean",
                    "BioAnth_WMean","BioAnth_Max","TreeCover_WMean","BioCropRisk","CropRisk","PastureRisk")]
head(PAdataS)
ncol(PAdataS)

#-#-# Add Realms #-#-#
PAdataS <- subset(PAdataS,RealmNr=="1"|RealmNr=="3"|RealmNr=="4"|RealmNr=="5"|RealmNr=="6"|RealmNr=="8")
as.numeric(as.character(PAdataS$RealmNr))
PAdataS$RealmName <- 0
PAdataS$RealmName[PAdataS$RealmNr == 1] <- "Australasian"
PAdataS$RealmName[PAdataS$RealmNr == 3] <- "Afrotropic"
PAdataS$RealmName[PAdataS$RealmNr == 4] <- "Indomalaya"
PAdataS$RealmName[PAdataS$RealmNr == 5] <- "Nearctic"
PAdataS$RealmName[PAdataS$RealmNr == 6] <- "Neotropic"
PAdataS$RealmName[PAdataS$RealmNr == 8] <- "Palearctic"

#-#-# Set the base data #-#-#
Base <- PAdataS[c(1:6,36)]
head(Base)

#-#-# Set the variable data #-#-#
VariableData <- PAdataS[c(8:35)]
head(VariableData)
VariableData[] <- lapply(VariableData[c(1:ncol(VariableData))], function(x) as.numeric(as.character(x))) # Change columns to numeric

#-#-# Format PCA input data #-#-#
PCAdata <- cbind(Base,VariableData)
head(PCAdata)
ncol(PCAdata)
nrow(PCAdata)

#-#-# Combine variables #-#-#
head(PCAdata)
PCAdata$Landuse_Stability <- rowMeans(PCAdata[c("BioCropRisk","CropRisk","PastureRisk")],na.rm=T)
head(PCAdata)

PCAdata$Species_Richness  <- rowMeans(PCAdata[c("Prop_SR_B","Prop_SR_M","Prop_SR_A","Prop_SR_R")],na.rm=T)
head(PCAdata)

PCAdata$Species_Endemism  <- rowMeans(PCAdata[c("RSR_B","RSR_M","RSR_A","RSR_R")],na.rm=T)
head(PCAdata)

PCAdata$Evolutionary_Diversity  <- rowMeans(PCAdata[c("PE_WMean_B","PE_WMean_M","PE_WMean_A","PE_WMean_R")],na.rm=T)
head(PCAdata)

PCAdata$Climatic_Stability  <- rowMeans(PCAdata[c("TO_PA_weightbird","TO_PA_weightmammal","TO_PA_weightamphibian","TO_PA_weightreptile")],na.rm=T)
head(PCAdata)


#-#-# Select variables used in PCA and change names for plot #-#-#
PCAdata <- PCAdata[c("SP_ID","Int_Name","RealmName","Species_Richness","Species_Endemism","Evolutionary_Diversity","BII_WMean",
                     "HFP_WMean","BioAnth_WMean","Climatic_Stability","TreeCover_WMean","Landuse_Stability","BaseCarbon_WMean","VulCarbon_WMean","IrrCarbon_WMean","Area_GIS")]


colnames(PCAdata) <- c("SP_ID","Int_Name","Realm","High_species_richness","High_species_endemism","High_evolutionary_diversity","High_biodiversity_intactness",
                       "Low_human_footprint","Low_anthrome_increase","High_climatic_stability","Low_forestcover_change","High_landuse_stability","High_carbon_storage",
                       "High_vulnerable_carbon","High_irreplaceable_carbon","Area")

#-#-# Remove NAs #-#-#
PCAdata <- na.omit(PCAdata)
nrow(PCAdata)

## Turn CS, LUS, HFP and Anthrome around!
PCAdata$Low_human_footprint <- (PCAdata$Low_human_footprint)* -1
PCAdata$High_landuse_stability <- (PCAdata$High_landuse_stability)* -1
PCAdata$High_climatic_stability <- (PCAdata$High_climatic_stability)* -1
PCAdata$Low_anthrome_increase <- (PCAdata$Low_anthrome_increase)* -1
PCAdata$Low_forestcover_change <- (PCAdata$Low_forestcover_change)* -1

#-#-# Scale variables and compute PCA #-#-#
PA.pca.cl <- prcomp(PCAdata[c(4:16)],center=T,scale=T) 
summary(PA.pca.cl) 
str(PA.pca.cl)

#-#-# Explained variance of components #-#-#
fviz_eig(PA.pca.cl)


#-#-# Axis plot (dimension 1 and 2) #-#-#
fviz_pca_var(PA.pca.cl, 
             col.var = "contrib", # Color by contributions to the PC 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping 
)

fviz_pca_ind(PA.pca.cl, label="none", habillage=PCAdata$Realm,
             palette="Dark2")


