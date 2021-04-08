#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#   Spatial distribution of top sites based on individual objectives  #
#        Map per realm top 5 sites per conservation objective         #
#                       Figure 3 main manuscript                      #
#                           Alke March 2021                           #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#


#-#-# Clear memory #-#-#
rm(list=ls())


#-#-# Load libraries #-#-#
library(ggplot2)


#-#-# Load the world map #-#-#
load_worldmap <- function(filename) {
  worldmap <- sf::st_read(filename, layer = "ne_50m_admin_0_countries")
  worldmap <- simplify_polygons(worldmap)
  return(worldmap)
}

  
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# Weightig the different objectives #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## Top sites Biodiversity, Wilderness Combined (Bio25%, Wild 25%, Size25%, LandUseStab10%, ClimateStab10%, CarbonStorage5%)
Filepath <- "/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/LL_analysis/Processed_data/Scaled_variables/"
AllDF <- list.files(Filepath)

#-#-# Set the weights for the scenarioes that should be plottet #-#-#
## Weights need to sum up to 1

#---# Biodiversity only #---#
Bio_WBio <- 1
Bio_WWild <- 0
Bio_WClimStab <- 0
Bio_WLUStab <- 0
Bio_WArea <- 0
Bio_WClimPro <- 0

Bio_weights <- as.data.frame(rbind(Bio_WBio,Bio_WClimStab,Bio_WLUStab,Bio_WWild,Bio_WArea,Bio_WClimPro))

#---# Wilderness only #---#
Wild_WBio <- 0
Wild_WWild <- 1
Wild_WClimStab <- 0
Wild_WLUStab <- 0
Wild_WArea <- 0
Wild_WClimPro <- 0

Wild_weights <- as.data.frame(rbind(Wild_WBio,Wild_WClimStab,Wild_WLUStab,Wild_WWild,Wild_WArea,Wild_WClimPro))

#---# Combination only #---#
Comb_WBio <- 0
Comb_WWild <- 1
Comb_WClimStab <- 0
Comb_WLUStab <- 0
Comb_WArea <- 0
Comb_WClimPro <- 0

Comb_weights <- as.data.frame(rbind(Comb_WBio,Comb_WClimStab,Comb_WLUStab,Comb_WWild,Comb_WArea,Comb_WClimPro))

Scenarioes <- cbind(Bio_weights,Wild_weights,Comb_weights)
colnames(Scenarioes) <- c("Biodiversity", "Wilderness", "Combination")


Ranks <- lapply(AllDF, function(d){
  
  RealmName <- strsplit(strsplit(d, split= "_")[[1]][8], split = ".csv")[[1]][1]
  print(RealmName)
  
  File <- read.csv(paste0(Filepath,d))
  
  Individual_weighting <- lapply(c(1:3), function(s){

    Scenario <- colnames(Scenarioes)[s]

    Weights <- as.vector(as.list(Scenarioes[s]))[[1]]

    WeightTable <- File[c("Int_Name","Biodiversity","ClimateStability","LandUseStability","Wilderness","Size","ClimateProtection")]

    ## Check min and max again
    for(i in 2:ncol(WeightTable)){
      print(max(na.omit(WeightTable[,i])))
    }

    for(i in 2:ncol(WeightTable)){
      print(min(na.omit(WeightTable[,i])))
    }

    ## Apply weighting to the different objectives
    Data <- as.matrix(WeightTable[2:7])
    Rank <- rowWeightedMeans(Data, w = Weights, na.rm = TRUE)
    WeightTable$Rank <- Rank

    ## Order and select the top 5 sites
    WeightTable <- WeightTable[order(WeightTable$Rank, decreasing = TRUE),]
    WeightTable05 <- WeightTable[1:5,]
    max(na.omit(WeightTable$Rank))
    min(na.omit(WeightTable$Rank))
    head(WeightTable)
    head(WeightTable05)

    WeightTable05$Scenario <- Scenario

    return(WeightTable05)

  })

  Result_table <-  as.data.frame(do.call(rbind,Individual_weighting))
  Result_table$Realm <- RealmName
  return(Result_table)
   
})

Map_table <- do.call(rbind,Ranks)
head(Map_table)

## Save the ranking 
write.csv(WeightTable,"Ranking_indomalaya_4Bio_2Wild_05ClimSt_05LUStab_2Size_1ClimP.csv")
#write.csv(WeightTable20,"Ranking_palearctic_top20_100_ClimProt.csv")








plot <- ggplot(selected_sites) +
    geom_sf(data = worldmap, fill = "black", color = "grey60", size = 0.2) +
    # Four types of points to show all and selected sites
    geom_point(data = pa_centroids,
               aes(x = x, y = y),
               shape = 16,
               size = 0.2,
               colour = "grey90") +
    geom_point(data = selected_sites,
               aes(x = x, y = y, color = suitability),
               shape = 18,
               size = 2.5) +
    scale_color_manual(name = "Suitability top sites:",
                       values = c("Very high" = "red",
                                  "High" = "orange",
                                  "Good" = "gold")) +
    coord_sf(xlim = c(-170, 180), ylim = c(-60, 90), expand = FALSE) +
    # Add shortened equator line
    geom_segment(
      aes(x = -180, xend = 180, y = 0, yend = 0),
      colour = "black",
      linetype = "dashed"
    ) +
    theme(legend.key = element_blank(), legend.position = "bottom",
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12)) +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    theme(axis.title = element_text(size = 16)) +
    theme(panel.background = element_rect(fill = "white", colour = "white")) +
    labs(x = "", y = "", title = "") +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()) +
    ggtitle(paste("Location top", n_sites, "sites", selection)) +
    theme(plot.title = element_text(size = 21, face = "bold", hjust = 0))
  return(plot)
}
