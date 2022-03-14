#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#          Supplementary analysis "LL Site Selection Tool"          # 
#   Checking for correlation among indicator and objective ranks    #
#            Figure S7 and S8 supplementary material                #
#                       Alke December 2021                          #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#


#-#-# Clear memory #-#-#
rm(list=ls())


#-#-# Load libraries #-#-#
library(plyr)
library(corrplot)
library(cowplot)
library(ggplot2)
library(grid)
library(gridExtra)
library(matrixStats)


#-#-# Get the final data frame (gridded PAs with values for all indicators) #-#-#
setwd("/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/LL_analysis/Processed_data/")
PAdata <- read.csv("Final_dataset_IUCN_WHS_KBA_country_correct_names.csv")
Outpath <- "/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/LL_analysis/Processed_data/Ranks_for_correlation_matrix/"


#-#-# Subset to columns that are included in the weighting #-#-#
PAdataS <- PAdata[c("SP_ID","Int_Name","PA_type","PA_DEF","IUCN_CAT","RealmNr","Area_Rep","Area_GIS","Prop_SR_B","Prop_SR_M","Prop_SR_A","Prop_SR_R",
                    "RSR_B","RSR_M","RSR_A","RSR_R","PE_WMean_B","PE_WMean_M","PE_WMean_A","PE_WMean_R","TO_PA_weightbird","TO_PA_weightmammal",
                    "TO_PA_weightamphibian","TO_PA_weightreptile","HFP_WMean","BaseCarbon_WMean","VulCarbon_WMean","IrrCarbon_WMean","BII_WMean",
                    "BioAnth_WMean","BioAnth_Max","TreeCover_WMean","BioCropRisk","CropRisk","PastureRisk")]
head(PAdataS)
ncol(PAdataS)


#-#-# List realm numbers for individual ranking tables #-#-#
RealmNr <- c("1","2","3","4","5","6","8","All")


#-#-# Scaling function #-#-#
scaleDF <- function(x){
  Scaling <- lapply(1:ncol(x),function(s){
    print(s)
    oneCol <- x[s]
    name <- colnames(oneCol)
    oneColS <- (oneCol - min(oneCol,na.rm=T)) / (max(oneCol,na.rm=T) - min(oneCol,na.rm=T))
    colnames(oneColS) <- name
    return(oneColS)
  })
  ScaledData <- do.call(cbind,Scaling)
  return(ScaledData)
}

#-#-# Ranking function #-#-#
rankDF <- function(x){
  ## Run through variables in two v
  Ranking <- lapply(2:ncol(x),function(r){
    print(r)
    oneCol <- x[c(1,r)]
    name <- colnames(oneCol)
    colnames(oneCol) <- c("ID","Rank")
    oneColR <- as.data.frame(oneCol[order(oneCol$Rank, decreasing = TRUE), ])
    oneColR$Rank <- seq(nrow(oneColR))
    oneColR[] <- lapply(oneColR, function(x) as.numeric(as.character(x)))
    colnames(oneColR) <- name
    return(oneColR)
  })
  
  RankedData <- Reduce(function(...) merge(..., all=T, by="ID"), Ranking)

  return(RankedData)
}


#--------# Scale the data and format for correlation matrix analysis #--------#
ScaleData <- lapply(RealmNr,function(x){
  
  #-#-# Subset the dataframe to one realm #-#-#
  if(x =="All"){
    PAdataR <- PAdataS #Global dataset  
  }else{
    PAdataR <- subset(PAdataS,RealmNr == x)
  }
  
  #-#-# Set realm name #-#-#
  if(x == "All"){Rname <- "Global"}
  if(x == 1){Rname <- "Australasian"}
  if(x == 2){Rname <- NA}
  if(x == 3){Rname <- "Afrotropical"}
  if(x == 4){Rname <- "Indomalayan"}
  if(x == 5){Rname <- "Nearctic"}
  if(x == 6){Rname <- "Neotropical"}
  if(x == 8){Rname <- "Palearctic"}
  
  print(Rname)
  
  #-#-# Summarize the variables into different objectives #-#-#
  BaseData <- PAdataR[1:6]
  Variables <- PAdataR[7:35]
  
  ## Change variable columns to numeric
  Variables[] <- lapply(Variables, function(x) as.numeric(as.character(x)))
  
  #-#-# Scale individual variables #-#-#
  VariablesScaled <- scaleDF(Variables[c(2:24,26:ncol(Variables))])
  ## Rename here already for variable names used in manuscript
  colnames(VariablesScaled) <- c("Large_size", "Hig_SR_B", "High_SR_M", "High_SR_A", "High_SR_R", 
                                 "High_SE_B", "High_SE_M", "High_SE_A", "High_SE_R", 
                                 "High_ED_B", "High_ED_M", "High_ED_A", "High_ED_R", 
                                 "High_CS_B", "High_CS_M", "High_CS_A", "High_CS_R",
                                 "Low_HFP", "High_CaS", "High_VCaS", "High_ICaS",
                                 "High_BII", "Low_BTAC", "Low_FCC", "Low_BCrop", "Low_Crop", "Low_Pasture")
  
  
  ## Scale individual variables between 0 and 1 and take mean across biodiversity, climatic stability land-use, 
  ## climate protection, and wilderness before final scaling
  
  #---# Biodiversity variables #---#
  ## Biodiversity has only variables were high values are good
  Biodiv <- Variables[3:14]
  Biodiv <- scaleDF(Biodiv)
  Biodiv <- as.data.frame(rowWeightedMeans(as.matrix(Biodiv), na.rm = T)) # Mean across all biodiversity variables - everyting counts equal 
  colnames(Biodiv) <- "Biodiversity"
  
  #---# Calculate wilderness measure #---#
  ## Wilderness has one variable were high values are good (BII) and two were high values are bad (HFP and biome to anthrome)
  Wild <- Variables[c(19,23,24)]
  Wild <- scaleDF(Wild)
  
  ## Turn HFP and Biome to anthrome around first so that 1 is the best site and 0 the worst
  Wild$HFP_WMean <- 1 - (Wild$HFP_WMean)
  Wild$BioAnth_WMean <- 1 - (Wild$BioAnth_WMean)
  
  ## Take mean accross HFP, BII and Biome to anthrome to get the combined wilderness measure
  Wild <- as.data.frame(rowMeans(Wild,na.rm=T))
  colnames(Wild) <- "Wilderness"
  
  #---# Climatic stability #---#
  ## For all climatic stability variables high values are bad
  ## Change tree cover to absolute numbers - high tree cover change either way is regarded as negative change
  Tree <- Variables[26]
  Tree <- abs(Tree)
  Tree <- scaleDF(Tree)
  
  ## Climatic stability of species communities 
  Clim <- Variables[15:18]
  Clim <- scaleDF(Clim)
  Clim <- as.data.frame(rowMeans(cbind(Clim,Tree),na.rm=T))
  colnames(Clim) <- "ClimateStability"
  Clim$ClimateStability <- 1 - (Clim$ClimateStability)
  
  #---# Land-use stability #---#
  ## For all land-use stability variables high values are bad
  Land <- Variables[27:29]
  Land <- scaleDF(Land)
  Land <- as.data.frame(rowMeans(Land,na.rm=T))
  colnames(Land) <- "LandUseStability"
  Land$LandUseStability <- 1 - (Land$LandUseStability)
  
  #---# Climate protection #---#
  Carbon <- Variables[20:22]
  Carbon <- scaleDF(Carbon) 
  Carbon <- as.data.frame(rowMeans(Carbon,na.rm=T))
  colnames(Carbon) <- "ClimateProtection"
  
  
  #-#-# Combine variables that are averaged and the ones that stayed the same in one dataframe #-#-#
  #---# Size #---#
  VariablesSub <- Variables[c("Area_GIS")]
  colnames(VariablesSub) <- "Size"
  ## Log size because there are a handful of sites that are huge and skew the distribution
  VariablesSub$Size <- log(VariablesSub$Size)
  
  
  ## Combined objectives
  VariablesSub$Biodiversity <- Biodiv$Biodiversity
  VariablesSub$Wilderness <- Wild$Wilderness
  VariablesSub$ClimateStability <- Clim$ClimateStability
  VariablesSub$LandUseStability <- Land$LandUseStability
  VariablesSub$ClimateProtection <- Carbon$ClimateProtection
  
  
  #-#-# Rescale the Objectives so that each spreads from 0 to 1 #-#-#
  Scaling <- lapply(1:ncol(VariablesSub),function(s){
    print(s)
    oneCol <- VariablesSub[s]
    name <- colnames(oneCol)
    oneColS <- (oneCol - min(oneCol,na.rm=T)) / (max(oneCol,na.rm=T) - min(oneCol,na.rm=T))
    colnames(oneColS) <- name
    return(oneColS)
  })
  
  ## Bind the scaled data 
  ScaledData <- do.call(cbind,Scaling)
  
  #-#-# Put final dataframe together #-#-#
  ## Rename base columns and new columns
  colnames(BaseData) <- c("SP_ID","Int_Name","PA_type","PA_DEF","IUCN_CAT","RealmNr")
  
  AggData <- cbind(BaseData,VariablesScaled,ScaledData[c(2,3,4,5,6,1)])
  
  #-#-# Write final dataframe #-#-#
  write.csv(AggData,paste0(Outpath,"Final_dataset_IUCN_WHS_KBA_for_correlation_matrix",Rname,".csv"))
})


#-#-# Rank the indicators and objectives individually and explore correlations #-#-#
Global <- read.csv(paste0(Outpath, "Final_dataset_IUCN_WHS_KBA_for_correlation_matrixGlobal.csv"))
head(Global)
ID <- 1:nrow(Global)
Global <- cbind(ID,Global)

## Select columns for analysis
Base <- Global[c(1,4:8)]
Data <- Global[c(1,9:ncol(Global))]

head(Data)

## Reorder variables with opposite impact
#The individual variables in the dataframe are scaled but climatic stability of species communities
#and land-use change are opposite of the other variables (high values are bad) and for tree cover
#change high values in both directions are bad

Data$Low_FCC <- 1-(Data$Low_FCC)
Data$High_CS_B <- 1-(Data$High_CS_B) #(i.e. low species turnover naming is confusing here)
Data$High_CS_M <- 1-(Data$High_CS_M)
Data$High_CS_A <- 1-(Data$High_CS_A)
Data$High_CS_R <- 1-(Data$High_CS_R)
Data$Low_BCrop <- 1-(Data$Low_BCrop)
Data$Low_Crop <- 1-(Data$Low_Crop)
Data$Low_Pasture <- 1-(Data$Low_Pasture)
Data$Low_BTAC <- 1-(Data$Low_BTAC)
Data$Low_HFP <- 1-(Data$Low_HFP)
# 
# plot(RankedData$Low_FCC~RankedData$High_CS_B)
# reg <- lm(RankedData$Low_FCC~RankedData$High_CS_B)
# coeff=coefficients(reg)
# # equation of the line : 
# eq = paste0("y = ", round(coeff[2],1), "*x ", round(coeff[1],1))
# # plot
# abline(reg, col="blue")
# 
# plot(Data$Low_FCC~Data$High_CS_M)
# reg <- lm(Data$Low_FCC~Data$High_CS_M)
# coeff=coefficients(reg)
# # equation of the line : 
# eq = paste0("y = ", round(coeff[2],1), "*x ", round(coeff[1],1))
# # plot
# abline(reg, col="blue")
# 
# plot(Data$Low_FCC~Data$High_CS_A)
# reg <- lm(Data$Low_FCC~Data$High_CS_A)
# coeff=coefficients(reg)
# # equation of the line : 
# eq = paste0("y = ", round(coeff[2],1), "*x ", round(coeff[1],1))
# # plot
# abline(reg, col="blue")
# 
# plot(Data$Low_FCC~Data$High_CS_R)
# reg <- lm(Data$Low_FCC~Data$High_CS_R)
# coeff=coefficients(reg)
# # equation of the line : 
# eq = paste0("y = ", round(coeff[2],1), "*x ", round(coeff[1],1))
# # plot
# abline(reg, col="blue")
# # hist(Data$High_CS_B)
# # min(na.omit(Data$Low_FCC))
# # max(na.omit(Data$Low_FCC))

#-# Calculate correlations #-#
CorD <- cor(Data[8:ncol(Data)], use = "pairwise.complete.obs") ## Calculate correlation between variables (Pearson by default)
CorForP <- abs(CorD) 

## Plot correlations
par(mfrow=c(1,1))
corrplot(CorD, method="circle", bg = "white",addgrid.col = "gray10", 
         tl.col = "black",tl.cex = 0.8, p.mat = CorForP, sig.level = 0.7)


#-# Look at correlation of ranks across objectives #-#
ObjectivesData <- Data[c("Biodiversity", "Wilderness", "ClimateStability", "LandUseStability", "ClimateProtection", "Size")]
ObjectivesDataForPlot <- ObjectivesData
colnames(ObjectivesDataForPlot) <- c("Biodiversity", "Ecosystem integrity", "Climatic stability", "Land-use stability", "Climate protection", "Size")
head(ObjectivesDataForPlot)
ObjectivesDataForPlot <- na.omit(ObjectivesDataForPlot)

## Calculate correlations
Cor0 <- cor(ObjectivesDataForPlot, use = "pairwise.complete.obs") ## Calculate correlation between variables (Pearson by default)
Cor0_M  <- Cor0
Cor0_M[Cor0_M == 1] <- NA
mean(as.matrix(Cor0_M), na.rm=TRUE)

corrplot(Cor0, type = 'lower', bg = "white", addCoef.col = 'black',
         addgrid.col = "gray10", tl.col = "black", tl.cex = 1, 
         title = "Conservation objectives",
         mar=c(0,0,1.5,0))


#-# Correlations cited in paper #-#
cor.test(ObjectivesData$Biodiversity, ObjectivesData$ClimateProtection,  method = c("pearson"))
## between ranks (r(1346) = .59, p < 0.01) # t= 26.61 / between ranks (r(1346) = .64, p < 0.01) # t= 30.74

cor.test(ObjectivesData$Biodiversity, ObjectivesData$LandUseStability,  method = c("pearson"))
## between ranks(r(1346) = -.35, p < 0.01) # t = -13.65 / between ranks(r(1346) = -.44, p < 0.01) # t = -18.37

cor.test(ObjectivesData$Biodiversity, ObjectivesData$ClimateStability,  method = c("pearson"))
## between ranks(r(1346) = .42, p < 0.01) # t = 16.77 / between ranks(r(1346) = .51, p < 0.01) # t = 21.98


#-# Calculate correlations for indicators within objectives #-#
## Biodiversity
BiodiversityDataSR <- as.data.frame(rowMeans(Data[c("Hig_SR_B", "High_SR_M", "High_SR_A", "High_SR_R")]))
colnames(BiodiversityDataSR) <- "Species_richness"
head(BiodiversityDataSR)

BiodiversityDataRSR <- as.data.frame(rowMeans(Data[c("High_SE_B", "High_SE_M", "High_SE_A", "High_SE_R")]))
colnames(BiodiversityDataRSR) <- "Species_endemism"
head(BiodiversityDataRSR)

BiodiversityDataED <- as.data.frame(rowMeans(Data[c("High_SE_R", "High_ED_B", "High_ED_M", "High_ED_A")]))
colnames(BiodiversityDataED) <- "Evolutionary_diversity"
head(BiodiversityDataED)

BiodiversityData <- cbind(BiodiversityDataSR,BiodiversityDataRSR,BiodiversityDataED)
CorB <- cor(BiodiversityData, use = "pairwise.complete.obs") ## Calculate correlation between variables (Pearson by default)
CorB_M  <- CorB
CorB_M[CorB_M == 1] <- NA
mean(as.matrix(CorB_M), na.rm=TRUE)
corrplot(CorB, addCoef.col = 'black', type = 'lower', bg = "white",
         addgrid.col = "gray10", tl.col = "black", tl.cex = 1, 
         title = "Biodiversity indicators",
         mar=c(0,0,1.5,0))
 
## Ecosystem integrity
EcosysIntegrData <- as.data.frame(Data[c("Low_HFP", "High_BII", "Low_BTAC")])
colnames(EcosysIntegrData) <- c("Human_footprint", "Biodiversity_intactness_index", "Recent_LU_change")
head(EcosysIntegrData)

plot(Data$High_BII,Data$Low_BTAC)
reg <- lm(Data$Low_HFP~Data$Low_BTAC)
coeff=coefficients(reg)
# equation of the line :
eq = paste0("y = ", round(coeff[2],1), "*x ", round(coeff[1],1))
# plot
abline(reg, col="blue")

CorE <- cor(EcosysIntegrData, use = "pairwise.complete.obs") ## Calculate correlation between variables (Pearson by default)
CorE_M  <- CorE
CorE_M[CorE_M == 1] <- NA
mean(as.matrix(CorE_M), na.rm=TRUE)
corrplot(CorE, addCoef.col = 'black', type = 'lower', bg = "white",
         addgrid.col = "gray10", tl.col = "black", tl.cex = 1, 
         title = "Ecosystem integrity indicators",
         mar=c(0,0,4,0))

## Climate protection
ClimProtData <- as.data.frame(Data[c("High_CaS", "High_VCaS", "High_ICaS")])
colnames(ClimProtData) <- c("Manageable_carbon", "Vulnerable_carbon", "Irreplaceable_carbon")
head(ClimProtData)

CorCP <- cor(ClimProtData, use = "pairwise.complete.obs") ## Calculate correlation between variables (Pearson by default)
CorCP_M  <- CorCP
CorCP_M[CorCP_M == 1] <- NA
mean(as.matrix(CorCP_M), na.rm=TRUE)
corrplot(CorCP, addCoef.col = 'black', type = 'lower', bg = "white",
         addgrid.col = "gray10", tl.col = "black", tl.cex = 1, 
         title = "Climate protection indicators",
         mar=c(0,0,4,0))


## Climatic stability
CSData <- Data[c("High_CS_B", "High_CS_M", "High_CS_A", "High_CS_R", "Low_FCC")]
colnames(CSData) <- c("Turnover_birds", "Turnover_mammals", "Turnover_Amphibians", "Turnover_Reptiles", "Tree_cover_change")
head(CSData)

## Calculate correlations
CorCS <- cor(CSData, use = "pairwise.complete.obs") ## Calculate correlation between variables (Pearson by default)
CorCS_M <- CorCS
CorCS_M[CorCS_M == 1] <- NA
mean(as.matrix(CorCS), na.rm=TRUE)

corrplot(CorCS, addCoef.col = 'black', type = 'lower', bg = "white",
         addgrid.col = "gray10", tl.col = "black",tl.cex = 1, 
         title = "Climatic stability",
         mar=c(0,0,1.5,0))


#-#-# Investigate impact of indicators on rank changes #-#-#
## Rank sites for all indicators and objectives
RankedData <- rankDF(Data)
head(RankedData)
nrow(RankedData)

## Combine site info and ranked data
DataMatrix <- merge(Base, RankedData, by="ID")
head(DataMatrix)
nrow(DataMatrix)

#-# Function too loop through all columns and get rank changes with every other column
RankChanges <- function(x){
  ## Run through all columns of the DF
  Data <- x
  Changes <- lapply(1:(ncol(x)-1), function(c){
    colNr <- c
    print(c)
    colList <- (colNr+1):ncol(x)
    Sub <- lapply(colList, function(rc){
      Rchange <- Data[colNr] - Data[rc]
      allChanges <- abs(Rchange)
      meanChange <- as.data.frame(colMeans(allChanges))
      Name <- paste0(colnames(Data[colNr]),"_",colnames(Data[rc]))
      #colnames(meanChange) <- Name
      rownames(meanChange) <- Name
      print(meanChange)
      return(meanChange)
    })
    OneObjective <- do.call(rbind,Sub)
    return(OneObjective)
  })
  return(Changes)
}

ReturnAllChanges <- function(x){
  ## Run through all columns of the DF
  Data <- x
  Changes <- lapply(1:(ncol(x)-1), function(c){
    colNr <- c
    print(c)
    colList <- (colNr+1):ncol(x)
    Sub <- lapply(colList, function(rc){
      Rchange <- Data[colNr] - Data[rc]
      #allChanges <- abs(Rchange)
      #meanChange <- as.data.frame(colMeans(allChanges))
      #Name <- paste0(colnames(Data[colNr]),"_",colnames(Data[rc]))
      #colnames(meanChange) <- Name
      colnames(Rchange) <- "allranks"
      print(Rchange)
      return(Rchange)
    })
    OneObjective <- do.call(rbind,Sub)
    return(OneObjective)
  })
  return(Changes)
}

head(DataMatrix)

## All objectives
AllObs <- DataMatrix[c(34:39)]
head(AllObs)
AllObs_Rank <- RankChanges(AllObs)
AllObs_Rank <- do.call(rbind,AllObs_Rank)
colnames(AllObs_Rank) <- "Mean_rank_change"
AllObs_Rank_Mean <- mean(AllObs_Rank$Mean_rank_change)

AllObs_Rank_all_ind <- ReturnAllChanges(AllObs)
AllObs_Rank_all_ind <- as.data.frame(do.call(rbind,AllObs_Rank_all_ind))
AllObsN <- nrow(AllObs_Rank_all_ind)
hist(AllObs_Rank_all_ind$allranks)
AllObsM <- round(mean(abs(AllObs_Rank_all_ind$allranks)),0)

#-# All objectives rank plot
AllAllObsM <- paste0("change = ",AllObsM)
AllAllObs <- ggplot(AllObs_Rank_all_ind, aes(x=allranks)) + 
  geom_histogram(color="black", fill="grey") +
  labs(title="All conservation objectives", x= "", y = "") +
  scale_x_continuous(limits = c(-1500, 1500)) +
  annotate(geom="text", x=680, y=2300, label="mean", color="black", size=12) +
  annotate(geom="text", x=800, y=2200, label=AllAllObsM, color="black", size=12) +
  theme(axis.text.y = element_text(size=24,colour="black"),
        axis.text.x = element_text(size=24,colour="black"),
        plot.title = element_text(size=32,colour="black",face="bold"),
        panel.background=element_rect(fill='white',colour="white"))
plot(AllAllObs)

#-#-# Save the final histogram plot #-#-#
setwd("/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/Result_plots/")
ggsave("Rank_changes_Histogram_AllObjectives_submission.tiff",AllAllObs,width=25, height=15, unit="in", dpi=300, bg="white")



## Biodiversity
Biodiv <- DataMatrix[c(8:19)]
head(Biodiv)
Biodiv_Rank <- RankChanges(Biodiv)
Biodiv_Rank <- do.call(rbind,Biodiv_Rank)
colnames(Biodiv_Rank) <- "Mean_rank_change"
Biodiv_Rank_Mean <- mean(Biodiv_Rank$Mean_rank_change)

Biodiv_Rank_all_ind <- ReturnAllChanges(Biodiv)
Biodiv_Rank_all_ind <- as.data.frame(do.call(rbind,Biodiv_Rank_all_ind))
BiodN <- nrow(Biodiv_Rank_all_ind)
hist(Biodiv_Rank_all_ind$allranks)
BiodM <- round(mean(abs(Biodiv_Rank_all_ind$allranks)),0)

## Split by taxa
Bird <- Biodiv[c("Hig_SR_B","High_SE_B","High_ED_B")]
head(Bird)
Bird_Rank_all_ind <- ReturnAllChanges(Bird)
Bird_Rank_all_ind <- as.data.frame(do.call(rbind,Bird_Rank_all_ind))
nrow(Bird_Rank_all_ind)
hist(Bird_Rank_all_ind$allranks)
BirdM <- round(mean(abs(Bird_Rank_all_ind$allranks)),0)

Mam <- Biodiv[c("High_SR_M","High_SE_M","High_ED_M")]
head(Mam)
Mam_Rank_all_ind <- ReturnAllChanges(Mam)
Mam_Rank_all_ind <- as.data.frame(do.call(rbind,Mam_Rank_all_ind))
nrow(Mam_Rank_all_ind)
hist(Mam_Rank_all_ind$allranks)
MammalM <- round(mean(abs(Mam_Rank_all_ind$allranks)),0)

Amp <- Biodiv[c("High_SR_A","High_SE_A","High_ED_A")]
head(Amp)
Amp_Rank_all_ind <- ReturnAllChanges(Amp)
Amp_Rank_all_ind <- as.data.frame(do.call(rbind,Amp_Rank_all_ind))
nrow(Amp_Rank_all_ind)
hist(Amp_Rank_all_ind$allranks)
AmphibianM <- round(mean(abs(Amp_Rank_all_ind$allranks)),0)

Rep <- Biodiv[c("High_SR_R","High_SE_R","High_ED_R")]
head(Rep)
Rep_Rank_all_ind <- ReturnAllChanges(Rep)
Rep_Rank_all_ind <- as.data.frame(do.call(rbind,Rep_Rank_all_ind))
nrow(Rep_Rank_all_ind)
hist(Rep_Rank_all_ind$allranks)
ReptileM <- round(mean(abs(Rep_Rank_all_ind$allranks)),0)

## Split by measure
SR <- Biodiv[c("Hig_SR_B","High_SR_M","High_SR_A","High_SR_R")]
head(SR)
SR_Rank_all_ind <- ReturnAllChanges(SR)
SR_Rank_all_ind <- as.data.frame(do.call(rbind,SR_Rank_all_ind))
nrow(SR_Rank_all_ind)
hist(SR_Rank_all_ind$allranks)
SRM <- round(mean(abs(SR_Rank_all_ind$allranks)),0)

SE <- Biodiv[c("High_SE_B","High_SE_M","High_SE_A","High_SE_R")]
head(SE)
SE_Rank_all_ind <- ReturnAllChanges(SE)
SE_Rank_all_ind <- as.data.frame(do.call(rbind,SE_Rank_all_ind))
nrow(SE_Rank_all_ind)
hist(SE_Rank_all_ind$allranks)
SEM <- round(mean(abs(SE_Rank_all_ind$allranks)),0)

ED <- Biodiv[c("High_ED_B","High_ED_M","High_ED_A","High_ED_R")]
head(ED)
ED_Rank_all_ind <- ReturnAllChanges(ED)
ED_Rank_all_ind <- as.data.frame(do.call(rbind,ED_Rank_all_ind))
nrow(ED_Rank_all_ind)
hist(ED_Rank_all_ind$allranks)
EDM <- round(mean(abs(ED_Rank_all_ind$allranks)),0)

#-# Biodiversity rank plot
AllBioM <- paste0("change = ",BiodM)
AllBio <- ggplot(Biodiv_Rank_all_ind, aes(x=allranks)) + 
            geom_histogram(color="black", fill="darkgreen") +
            labs(title="All biodiversity indicators", x= "", y = "") +
            scale_x_continuous(limits = c(-1500, 1500)) +
            annotate(geom="text", x=460, y=13800, label="mean", color="black", size=8) +
            annotate(geom="text", x=800, y=13000, label=AllBioM, color="black", size=8) +
          theme(axis.text.y = element_text(size=18,colour="black"),
            axis.text.x = element_text(size=18,colour="black"),
            plot.title = element_text(size=24,colour="black",face="bold"),
            panel.background=element_rect(fill='white',colour="white"))
plot(AllBio)

AllBirdM <- paste0("change = ",BirdM)
Birds <- ggplot(Bird_Rank_all_ind, aes(x=allranks)) + 
  geom_histogram(color="black", fill="olivedrab3", binwidth=50) +
  labs(title="Birds", x="", y = "") +
  scale_x_continuous(limits = c(-1500, 1500)) +
  scale_y_continuous(limits = c(0, 700)) +
  annotate(geom="text", x=460, y=530, label="mean", color="black", size=8) +
  annotate(geom="text", x=800, y=500, label=AllBirdM, color="black", size=8) +
  theme(axis.text.y = element_text(size=18,colour="black"),
        axis.text.x = element_text(size=18,colour="black"),
        plot.title = element_text(size=24,colour="black",face="bold"),
        panel.background=element_rect(fill='white',colour="white"))
plot(Birds)

AllMammalM <- paste0("change = ",MammalM)
Mammals <- ggplot(Mam_Rank_all_ind, aes(x=allranks)) + 
  geom_histogram(color="black", fill="olivedrab3", binwidth=50) +
  labs(title="Mammals", x="", y = "") +
  scale_x_continuous(limits = c(-1500, 1500)) +
  scale_y_continuous(limits = c(0, 700)) +
  annotate(geom="text", x=460, y=530, label="mean", color="black", size=8) +
  annotate(geom="text", x=800, y=500, label=AllMammalM, color="black", size=8) +
  theme(axis.text.y = element_text(size=18,colour="black"),
        axis.text.x = element_text(size=18,colour="black"),
        plot.title = element_text(size=24,colour="black",face="bold"),
        panel.background=element_rect(fill='white',colour="white"))
plot(Mammals)

AllAmphibianM <- paste0("change = ",AmphibianM)
Amphibians <- ggplot(Amp_Rank_all_ind, aes(x=allranks)) + 
  geom_histogram(color="black", fill="olivedrab3", binwidth=50) +
  labs(title="Ambhibians", x="", y = "") +
  scale_x_continuous(limits = c(-1500, 1500)) +
  scale_y_continuous(limits = c(0, 700)) +
  annotate(geom="text", x=460, y=530, label="mean", color="black", size=8) +
  annotate(geom="text", x=800, y=500, label=AllAmphibianM, color="black", size=8) +
  theme(axis.text.y = element_text(size=18,colour="black"),
        axis.text.x = element_text(size=18,colour="black"),
        plot.title = element_text(size=24,colour="black",face="bold"),
        panel.background=element_rect(fill='white',colour="white"))
plot(Amphibians)

AllReptilesM <- paste0("change = ",ReptileM)
Reptiles <- ggplot(Rep_Rank_all_ind, aes(x=allranks)) + 
  geom_histogram(color="black", fill="olivedrab3", binwidth=50) +
  labs(title="Reptiles", x="", y = "") +
  scale_x_continuous(limits = c(-1500, 1500)) +
  scale_y_continuous(limits = c(0, 1000)) +
  annotate(geom="text", x=460, y=750, label="mean", color="black", size=8) +
  annotate(geom="text", x=800, y=700, label=AllReptilesM, color="black", size=8) +
  theme(axis.text.y = element_text(size=18,colour="black"),
        axis.text.x = element_text(size=18,colour="black"),
        plot.title = element_text(size=24,colour="black",face="bold"),
        panel.background=element_rect(fill='white',colour="white"))
plot(Reptiles)

AllSRM <- paste0("change = ",SRM)
SR <- ggplot(SR_Rank_all_ind, aes(x=allranks)) + 
  geom_histogram(color="black", fill="darkolivegreen4", binwidth=50) +
  labs(title="Species richness", x="", y = "") +
  scale_x_continuous(limits = c(-1500, 1500)) +
  scale_y_continuous(limits = c(0, 1200)) +
  annotate(geom="text", x=540, y=1060, label="mean", color="black", size=8) +
  annotate(geom="text", x=900, y=1000, label=AllSRM, color="black", size=8) +
  theme(axis.text.y = element_text(size=18,colour="black"),
        axis.text.x = element_text(size=18,colour="black"),
        plot.title = element_text(size=24,colour="black",face="bold"),
        panel.background=element_rect(fill='white',colour="white"))
plot(SR)

AllSEM <- paste0("change = ",SEM)
SE <- ggplot(SE_Rank_all_ind, aes(x=allranks)) + 
  geom_histogram(color="black", fill="darkolivegreen4", binwidth=50) +
  labs(title="Species endemism", x="", y = "") +
  scale_x_continuous(limits = c(-1500, 1500)) +
  scale_y_continuous(limits = c(0, 1200)) +
  annotate(geom="text", x=540, y=1060, label="mean", color="black", size=8) +
  annotate(geom="text", x=900, y=1000, label=AllSEM, color="black", size=8) +
  theme(axis.text.y = element_text(size=18,colour="black"),
        axis.text.x = element_text(size=18,colour="black"),
        plot.title = element_text(size=24,colour="black",face="bold"),
        panel.background=element_rect(fill='white',colour="white"))
plot(SE)

AllEDM <- paste0("change = ",EDM)
ED <- ggplot(ED_Rank_all_ind, aes(x=allranks)) + 
  geom_histogram(color="black", fill="darkolivegreen4", binwidth=50) +
  labs(title="Genetic diversity", x="", y = "") +
  scale_x_continuous(limits = c(-1500, 1500)) +
  scale_y_continuous(limits = c(0, 1200)) +
  annotate(geom="text", x=540, y=1060, label="mean", color="black", size=8) +
  annotate(geom="text", x=900, y=1000, label=AllEDM, color="black", size=8) +
  theme(axis.text.y = element_text(size=18,colour="black"),
        axis.text.x = element_text(size=18,colour="black"),
        plot.title = element_text(size=24,colour="black",face="bold"),
        panel.background=element_rect(fill='white',colour="white"))
plot(ED)

CombHist <- arrangeGrob(AllBio,SR,SE,ED,Birds,Mammals,Amphibians,Reptiles,
                          widths = c(2,2,2,2),
                          heights = c(1,1),
                          ncol = 4,
                          nrow = 2)
plot(CombHist)

#-#-# Save the final histogram plot #-#-#
setwd("/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/Result_plots/")
ggsave("Rank_changes_Histogram_Biodiversity_submission.tiff",CombHist,width=25, height=15, unit="in", dpi=300, bg="white")



## Ecosystem Integrity
Wild <- DataMatrix[c("Low_HFP", "High_BII", "Low_BTAC")]
head(Wild)
Wild_Rank <- RankChanges(Wild)
Wild_Rank <- do.call(rbind,Wild_Rank)
colnames(Wild_Rank) <- "Mean_rank_change"
Wild_Rank_Mean <- mean(Wild_Rank$Mean_rank_change)

Wild_Rank_all_ind <- ReturnAllChanges(Wild)
Wild_Rank_all_ind <- as.data.frame(do.call(rbind,Wild_Rank_all_ind))
WildN <- nrow(Wild_Rank_all_ind)
hist(Wild_Rank_all_ind$allranks)
WildM <- round(mean(abs(Wild_Rank_all_ind$allranks)),0)

#-# Combined rank plot Wild
AllWildM <- paste0("change = ",WildM)
AllWild <- ggplot(Wild_Rank_all_ind, aes(x=allranks)) + 
  geom_histogram(color="black", fill="orangered3") +
  labs(title="All ecosystem integrity indicators", x= "", y = "") +
  scale_x_continuous(limits = c(-1500, 1500)) +
  scale_y_continuous(limits = c(0, 450)) +
  annotate(geom="text", x=940, y=320, label="mean", color="black", size=8) +
  annotate(geom="text", x=1100, y=300, label=AllBioM, color="black", size=8) +
  theme(axis.text.y = element_text(size=18,colour="black"),
        axis.text.x = element_text(size=18,colour="black"),
        plot.title = element_text(size=24,colour="black",face="bold"),
        panel.background=element_rect(fill='white',colour="white"))
plot(AllWild)


## Climate protection
ClimP <- DataMatrix[c("High_CaS", "High_VCaS", "High_ICaS")]
head(ClimP)
ClimP_Rank <- RankChanges(ClimP)
ClimP_Rank <- do.call(rbind,ClimP_Rank)
colnames(ClimP_Rank) <- "Mean_rank_change"
ClimP_Rank_Mean <- mean(ClimP_Rank$Mean_rank_change)

ClimP_Rank_all_ind <- ReturnAllChanges(ClimP)
ClimP_Rank_all_ind <- as.data.frame(do.call(rbind,ClimP_Rank_all_ind))
ClimPN <- nrow(ClimP_Rank_all_ind)
hist(ClimP_Rank_all_ind$allranks)
ClimPM <- round(mean(abs(ClimP_Rank_all_ind$allranks)),0)

#-# Rank plot Carbon
AllClimPM <- paste0("change = ",ClimPM)
AllClimP <- ggplot(ClimP_Rank_all_ind, aes(x=allranks)) + 
  geom_histogram(color="black", fill="royalblue4") +
  labs(title="All climate protection indicators", x= "", y = "") +
  scale_x_continuous(limits = c(-1500, 1500)) +
  annotate(geom="text", x=360, y=1600, label="mean", color="black", size=8) +
  annotate(geom="text", x=500, y=1500, label=AllClimPM, color="black", size=8) +
  theme(axis.text.y = element_text(size=18,colour="black"),
        axis.text.x = element_text(size=18,colour="black"),
        plot.title = element_text(size=24,colour="black",face="bold"),
        panel.background=element_rect(fill='white',colour="white"))
plot(AllClimP)


## Climatic stability
ClimS <- DataMatrix[c("High_CS_B", "High_CS_M", "High_CS_A", "High_CS_R", "Low_FCC")]
## Combine climatic stability of communities across taxa
# ClimSSumTax <- round(as.data.frame(rowMeans(DataMatrix[c("High_CS_B", "High_CS_M", "High_CS_A", "High_CS_R")])),0)
# colnames(ClimSSumTax) <- "High_CS_All"
# ClimS <- cbind(ClimS[5],ClimSSumTax)
head(ClimS)
ClimS_Rank <- RankChanges(ClimS)
ClimS_Rank <- do.call(rbind,ClimS_Rank)
colnames(ClimS_Rank) <- "Mean_rank_change"
ClimS_Rank_Mean <- mean(ClimS_Rank$Mean_rank_change)

ClimS_Rank_all_ind <- ReturnAllChanges(ClimS)
ClimS_Rank_all_ind <- as.data.frame(do.call(rbind,ClimS_Rank_all_ind))
ClimSN <- nrow(ClimS_Rank_all_ind)
hist(ClimS_Rank_all_ind$allranks)
ClimSM <- round(mean(abs(ClimS_Rank_all_ind$allranks)),0)


ClimSTaxa <- DataMatrix[c("High_CS_B", "High_CS_M", "High_CS_A", "High_CS_R")]
head(ClimSTaxa)
ClimSTaxa_Rank <- RankChanges(ClimSTaxa)
ClimSTaxa_Rank <- do.call(rbind,ClimSTaxa_Rank)
colnames(ClimSTaxa_Rank) <- "Mean_rank_change"
ClimSTaxa_Rank_Mean <- mean(ClimSTaxa_Rank$Mean_rank_change)

ClimSTaxa_Rank_all_ind <- ReturnAllChanges(ClimSTaxa)
ClimSTaxa_Rank_all_ind <- as.data.frame(do.call(rbind,ClimSTaxa_Rank_all_ind))
ClimSTaxaN <- nrow(ClimSTaxa_Rank_all_ind)
hist(ClimSTaxa_Rank_all_ind$allranks)
ClimSTaxaM <- round(mean(abs(ClimSTaxa_Rank_all_ind$allranks)),0)

#-# Combined rank plot Climatic stability
AllClimSM <- paste0("change = ",ClimSM)
AllClimS <- ggplot(ClimS_Rank_all_ind, aes(x=allranks)) + 
  geom_histogram(color="black", fill="gold1") +
  labs(title="All climatic stability indicators", x= "", y = "") +
  scale_x_continuous(limits = c(-1500, 1500)) +
  scale_y_continuous(limits = c(0, 2500)) +
  annotate(geom="text", x=490, y=2000, label="mean", color="black", size=8) +
  annotate(geom="text", x=650, y=1800, label=AllClimSM, color="black", size=8) +
  theme(axis.text.y = element_text(size=18,colour="black"),
        axis.text.x = element_text(size=18,colour="black"),
        plot.title = element_text(size=24,colour="black",face="bold"),
        panel.background=element_rect(fill='white',colour="white"))
plot(AllClimS)

#-# Combined rank plot Climatic stability taxa
AllClimSTaxaM <- paste0("change = ",ClimSTaxaM)
AllClimSTaxa <- ggplot(ClimSTaxa_Rank_all_ind, aes(x=allranks)) + 
  geom_histogram(color="black", fill="orange") +
  labs(title="Climatic stability of communities", x= "", y = "") +
  scale_x_continuous(limits = c(-1500, 1500)) +
  scale_y_continuous(limits = c(0, 2000)) +
  annotate(geom="text", x=490, y=1100, label="mean", color="black", size=8) +
  annotate(geom="text", x=650, y=1000, label=AllClimSTaxaM, color="black", size=8) +
  theme(axis.text.y = element_text(size=18,colour="black"),
        axis.text.x = element_text(size=18,colour="black"),
        plot.title = element_text(size=24,colour="black",face="bold"),
        panel.background=element_rect(fill='white',colour="white"))
plot(AllClimSTaxa)

CombHistOther <- arrangeGrob(AllWild,AllClimP,AllClimS,AllClimSTaxa,
                        widths = c(2,2),
                        heights = c(1,1),
                        ncol = 2,
                        nrow = 2)
plot(CombHistOther)

#-#-# Save the final histogram plot #-#-#
setwd("/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/Result_plots/")
ggsave("Rank_changes_Histogram_other_objectives_submission.tiff",CombHistOther,width=25, height=15, unit="in", dpi=300, bg="white")


