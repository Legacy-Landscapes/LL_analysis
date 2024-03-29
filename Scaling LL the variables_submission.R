#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#     Ranking the potential LL sites based on conservation objectives     #
#           Biodiversity - Wilderness - Climatic stability Area           #
#                Land-use stability - Climate protection                  #
#     Variables are scaled 0 to 1 and and ranking order is set to 1       #
#                 1 = best and 0 worst for each variable                  #
#                           Alke - June 2020                              #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

#-#-# Clear memory #-#-#
rm(list=ls())


#-#-# Load libraries #-#-#
library(matrixStats)
library(e1071)


#-#-# Get the final data frame (gridded PAs with values for all indicators) #-#-#
setwd("/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/LL_analysis/Processed_data/")
PAdata <- read.csv("Final_dataset_IUCN_WHS_KBA_country_oda.csv")
Outpath <- "/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/LL_analysis/Processed_data/Scaled_variables_country/"


#-#-# Subset to columns that are included in the weighting #-#-#
PAdataS <- PAdata[c("SP_ID","Int_Name","PA_type","IUCN_CAT","COUNTRY","ISO3","ODA_status","ODA","RealmNr","Area_Rep","Area_GIS","Prop_SR_B","Prop_SR_M","Prop_SR_A","Prop_SR_R",
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


#--------# Scale the data and format as input for the shiny app #--------#

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
    BaseData <- PAdataR[1:9] #basic information about the site
    Variables <- PAdataR[10:38] #site data
    
    ## Change variable columns to numeric
    Variables[] <- lapply(Variables, function(x) as.numeric(as.character(x)))

    ## Scale individual variables between 0 and 1 and take mean across biodiversity, climatic stability land-use, 
    ## climate protection, and ecosystem integrity before final scaling
  
    #---# Biodiversity variables #---#
     ## Biodiversity has only variables were high values are good
      Biodiv <- Variables[3:14]
      Biodiv <- scaleDF(Biodiv)
      Biodiv <- as.data.frame(rowWeightedMeans(as.matrix(Biodiv), na.rm = T)) # Mean across all biodiversity variables - everything counts equal 
      colnames(Biodiv) <- "Biodiversity"
    
    #---# Calculate wilderness measure #---#
      ## Ecosystem integrity has one variable were high values are good (BII) and two were high values are bad (HFP and biome to anthrome)
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
    AggData <- cbind(BaseData,ScaledData$Biodiversity,ScaledData$Wilderness,ScaledData$ClimateStability,
                 ScaledData$LandUseStability,ScaledData$Size,ScaledData$ClimateProtection)

    colnames(AggData) <- c("SP_ID","Int_Name","PA_type","IUCN_CAT","COUNTRY","ISO3","ODA_status","ODA","RealmNr","Biodiversity","Wilderness","ClimateStability","LandUseStability","Size","ClimateProtection")
    
    
  #-#-# Write final dataframe #-#-#
  write.csv(AggData,paste0(Outpath,"Final_dataset_IUCN_WHS_KBA_for_weighting_",Rname,".csv"))
})

