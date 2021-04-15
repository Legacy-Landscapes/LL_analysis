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
library(matrixStats)
library(scales)

  
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
Comb_WBio <- 0.25
Comb_WWild <- 0.25
Comb_WClimStab <- 0.1
Comb_WLUStab <- 0.1
Comb_WArea <- 0.25
Comb_WClimPro <- 0.05

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
Plot_map <- Map_table[c("Int_Name", "Scenario", "Realm")]
head(Plot_map)

## Save the ranking 
#write.csv(WeightTable,"Ranking_indomalaya_4Bio_2Wild_05ClimSt_05LUStab_2Size_1ClimP.csv")

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# Plot the hotspots map #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

#-#-# Set the file paths #-#-#
world <- "/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/Large_data_files/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp"
centroids <- "/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/LL_analysis/Processed_data/Global_IUCN_WHS_KBA_centroids.csv"
  
  
#-#-# Load the world map #-#-#
worldmap <- load_worldmap(world)
#africa <- subset(worldmap, CONTINENT == "Africa")
#plot(africa)

#-#-# Get the coordinates for the selected sites #-#-#
sites <- read.csv(centroids)[2:4]
colnames(sites) <- c("Int_Name", "x", "y")
#selected <- Plot_map$Int_Name
#selected_sites <- filter(sites, PAname %in% selected)
selected_sites <- merge(Plot_map, sites, by = "Int_Name", all.x = T)
head(selected_sites)


#-#-# Nearctic #-#-# 
selected_sites_n <- subset(selected_sites, Realm == "Nearctic")

nearctic <- ggplot(selected_sites_n) +
  geom_sf(data = worldmap, fill = "black", color = "grey60", size = 0.2) +
  geom_point(data = selected_sites,
             aes(x = x, y = y, color = Scenario),
             shape = 18,
             size = 2.5) +
  scale_color_manual(name = "Objective:",
                     values = c("Biodiversity" = "forestgreen",
                                "Wilderness" = "red",
                                "Combination" = "grey")) +
  coord_sf(xlim = c(-168, -17), ylim = c(12, 90), expand = FALSE) +
  ## Add shortened equator line
  geom_segment(
    aes(x = -180, xend = -17, y = 0, yend = 0),
    colour = "darkred",
    linetype = "dashed") +
  theme(legend.position = "none")+
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme(axis.title = element_text(size = 16)) +
  theme(panel.background = element_rect(fill = "white", colour = "white")) +
  labs(x = "", y = "", title = "") +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()) +
  ggtitle("Nearctic") +
  theme(plot.title = element_text(size = 21, face = "bold", hjust = 0))
plot(nearctic)


#-#-# Neotropical #-#-# 
neotropic <- ggplot(selected_sites) +
  geom_sf(data = worldmap, fill = "black", color = "grey60", size = 0.2) +
  geom_point(data = selected_sites,
             aes(x = x, y = y, color = Scenario),
             shape = 18,
             size = 2.5) +
  scale_color_manual(name = "Objective:",
                     values = c("Biodiversity" = "forestgreen",
                                "Wilderness" = "red",
                                "Combination" = "grey")) +
  coord_sf(xlim = c(-90, -30), ylim = c(-55, 15), expand = FALSE) +
  ## Add shortened equator line
  geom_segment(
    aes(x = -90, xend = -30, y = 0, yend = 0),
    colour = "darkred",
    linetype = "dashed") +
  theme(legend.position = "none")+
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme(axis.title = element_text(size = 16)) +
  theme(panel.background = element_rect(fill = "white", colour = "white")) +
  labs(x = "", y = "", title = "") +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()) +
  ggtitle("Neotropic") +
  theme(plot.title = element_text(size = 21, face = "bold", hjust = 0))
plot(neotropic)


#-#-# Palearctic #-#-# 
palearctic <- ggplot(selected_sites) +
  geom_sf(data = worldmap, fill = "black", color = "grey60", size = 0.2) +
  geom_point(data = selected_sites,
             aes(x = x, y = y, color = Scenario),
             shape = 18,
             size = 2.5) +
  scale_color_manual(name = "Objective:",
                     values = c("Biodiversity" = "forestgreen",
                                "Wilderness" = "red",
                                "Combination" = "grey")) +
  coord_sf(xlim = c(-12, 180), ylim = c(25, 83), expand = FALSE) +
  theme(legend.position = "none")+
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme(axis.title = element_text(size = 16)) +
  theme(panel.background = element_rect(fill = "white", colour = "white")) +
  labs(x = "", y = "", title = "") +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()) +
  ggtitle("Palearctic") +
  theme(plot.title = element_text(size = 21, face = "bold", hjust = 0))
plot(palearctic)


#-#-# Afrotropical #-#-# 
afrotropic <- ggplot(selected_sites) +
  geom_sf(data = worldmap, fill = "black", color = "grey60", size = 0.2) +
  geom_point(data = selected_sites,
             aes(x = x, y = y, color = Scenario),
             shape = 18,
             size = 2.5) +
  scale_color_manual(name = "Objective:",
                     values = c("Biodiversity" = "forestgreen",
                                "Wilderness" = "red",
                                "Combination" = "grey")) +
  coord_sf(xlim = c(-20, 60), ylim = c(-40, 20), expand = FALSE) +
  ## Add shortened equator line
  geom_segment(
    aes(x = -20, xend = 60, y = 0, yend = 0),
    colour = "darkred",
    linetype = "dashed") +
  theme(legend.position = "none")+
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme(axis.title = element_text(size = 16)) +
  theme(panel.background = element_rect(fill = "white", colour = "white")) +
  labs(x = "", y = "", title = "") +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()) +
  ggtitle("Afrotropic") +
  theme(plot.title = element_text(size = 21, face = "bold", hjust = 0))
plot(afrotropic)


#-#-# Indomalayan #-#-# 
indomalaya<- ggplot(selected_sites) +
  geom_sf(data = worldmap, fill = "black", color = "grey60", size = 0.2) +
  geom_point(data = selected_sites,
             aes(x = x, y = y, color = Scenario),
             shape = 18,
             size = 2.5) +
  scale_color_manual(name = "Objective:",
                     values = c("Biodiversity" = "forestgreen",
                                "Wilderness" = "red",
                                "Combination" = "grey")) +
  coord_sf(xlim = c(67, 130), ylim = c(25, -9), expand = FALSE) +
  ## Add shortened equator line
  geom_segment(
    aes(x = 67, xend = 130, y = 0, yend = 0),
    colour = "darkred",
    linetype = "dashed") +
  theme(legend.position = "none")+
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme(axis.title = element_text(size = 16)) +
  theme(panel.background = element_rect(fill = "white", colour = "white")) +
  labs(x = "", y = "", title = "") +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()) +
  ggtitle("Indomalaya") +
  theme(plot.title = element_text(size = 21, face = "bold", hjust = 0))
plot(indomalaya)


#-#-# Australasia #-#-# 
#austrasia <- subset(worldmap, CONTINENT == "Oceania" | CONTINENT == "Asia")

australasia <- ggplot(selected_sites) +
  geom_sf(data = worldmap, fill = "black", color = "grey60", size = 0.2) +
  coord_sf(xlim = c(110, 160), ylim = c(1, -50), expand = FALSE) +
  geom_point(data = selected_sites,
             aes(x = x, y = y, color = Scenario),
             shape = 18,
             size = 2.5) +
  scale_color_manual(name = "Objective:",
                     values = c("Biodiversity" = "forestgreen",
                                "Wilderness" = "red",
                                "Combination" = "grey")) +
  theme(legend.position = "none")+
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme(axis.title = element_text(size = 16)) +
  theme(panel.background = element_rect(fill = "white", colour = "white")) +
  labs(x = "", y = "", title = "") +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()) +
  ggtitle("Australasia") +
  theme(plot.title = element_text(size = 21, face = "bold", hjust = 0))
plot(australasia)


## Create legend for the PCA plot
## Function to extract plot legends
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

## ggplot to extract legend from for final plot
legend_map <- ggplot(selected_sites) +
  coord_sf(xlim = c(110, 160), ylim = c(1, -50), expand = FALSE) +
  geom_point(data = selected_sites,
             aes(x = x, y = y, color = Scenario),
             shape = 18,
             size = 2.5) +
  scale_color_manual(name = "Objective:",
                     values = c("Biodiversity" = "forestgreen",
                                "Wilderness" = "red",
                                "Combination" = "grey")) +
      theme(legend.key = element_blank(), legend.position = "bottom",
            legend.title = element_text(size = 14, face = "bold"),
            legend.text = element_text(size = 12)) +
      guides(colour = guide_legend(override.aes = list(size = 5)))

plot(legend_map)
ScatterLegend <- g_legend(legend_map)
plot(ScatterLegend)


#-#-# Combine the PCA plots and the legend #-#-#  
CombGlobal <- arrangeGrob(nearctic,palearctic,neotropic,afrotropic,indomalaya,australasia,ScatterLegend,ScatterLegend,
                          widths = c(0.5,0.5),
                          heights = c(1,1,1,0.2),
                          ncol = 2,
                          nrow = 4)
plot(CombGlobal)

#-#-# Save the final map #-#-#
setwd("/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/Result_plots/")
ggsave("Map Legacy Landscapes Conservation Objectives_separated realm plots.tiff",CombGlobal,width=25, height=15, unit="in", dpi=300, bg="white")


#---#---#---# Other plot otion #---#---#---#
#-#-# Set the file paths #-#-#
realmpath <- "/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Site selection analysis/Other files/Realms and Regions/CMEC regions & realms/WWF_Realms_gridded_05.Rdata"
centroids <- "/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/LL_analysis/Processed_data/Global_IUCN_WHS_KBA_centroids.csv"
world <- "/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/Large_data_files/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp"
centroids <- "/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/LL_analysis/Processed_data/Global_IUCN_WHS_KBA_centroids.csv"


#-#-# Get the coordinates for the selected sites #-#-#
sites <- read.csv(centroids)[2:4]
colnames(sites) <- c("Int_Name", "x", "y")
#selected <- Plot_map$Int_Name
#selected_sites <- filter(sites, PAname %in% selected)
selected_sites <- merge(Plot_map, sites, by = "Int_Name", all.x = T)
head(selected_sites)


#-#-# Load the world map #-#-#
worldmap <- load_worldmap(world)


#-#-# Load the realm map #-#-#
realmmap <- get(load(realmpath))
head(realmmap)
realmmap <- na.omit(realmmap)
realmmap <- subset(realmmap, !(RealmWWF == 2))
realmmap <- subset(realmmap, !(RealmWWF == 7))

realmmap$RealmName <- 0
realmmap$RealmName[realmmap$RealmWWF == 1] <- "Australasia"
realmmap$RealmName[realmmap$RealmWWF == 3] <- "Afrotropic"
realmmap$RealmName[realmmap$RealmWWF == 4] <- "Indomalaya"
realmmap$RealmName[realmmap$RealmWWF == 5] <- "Nearctic"
realmmap$RealmName[realmmap$RealmWWF == 6] <- "Neotropic"
realmmap$RealmName[realmmap$RealmWWF == 8] <- "Palearctic"

colours <- c("steelblue3","orange3","navajowhite3","darkolivegreen3","plum4","brown")
 
map <- ggplot(selected_sites) +
  geom_raster(data = realmmap, aes( x = x, y = y, fill = as.factor(RealmName), alpha = 0.2)) +
  scale_fill_manual("Realm",values = colours) +
  geom_sf(data = worldmap, fill = NA, color = "grey60", size = 0.2) +
  coord_sf(xlim = c(-170, 180), ylim = c(90, -60), expand = FALSE) +
  geom_point(data = selected_sites,
             aes(x = x, y = y, color = Scenario),
             shape = 18,
             size = 8,
             alpha = 0.7) +
  scale_color_manual(name = "Objective:",
                     values = c("Biodiversity" = "darkgreen",
                                "Wilderness" = "red",
                                "Combination" = "black")) +
  geom_segment(
    aes(x = -180, xend = 180, y = 0, yend = 0), colour = "black", linetype = "dashed") +
  theme(legend.key = element_blank(), legend.position = "bottom",
        legend.title = element_text(size = 22, face = "bold"),
        legend.text = element_text(size = 20)) +
  theme(axis.text=element_text(size=20))+
  theme(panel.background=element_rect(fill='white',colour="white"))+ # Remove the background
  labs(x="", y="", title="") + # Remove axis titles
  guides(colour = guide_legend(override.aes = list(size = 8)))

plot(map) 

#-#-# Save the final map #-#-#
setwd("/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/Result_plots/")
ggsave("Map Legacy Landscapes Conservation Objectives.tiff",map,width=25, height=15, unit="in", dpi=300, bg="white")



