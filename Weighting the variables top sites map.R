#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#   Spatial distribution of top sites based on individual objectives    #
#           flexible for diffrent conservation scenarioes               #
#        Map per realm top 5 sites per conservation objective           #
#                       Figure 3 main manuscript                        #
#                           Alke March 2021                             #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#


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

#-#-# Set the weights for the scenarioes that should be plotted #-#-#
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

#---# Combination only (LL example) #---#
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


#---#---#---# Other plot option #---#---#---#
# #-#-# Set the file paths #-#-#
# realmpath <- "/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Site selection analysis/Other files/Realms and Regions/CMEC regions & realms/WWF_Realms_gridded_05.Rdata"
# centroids <- "/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/LL_analysis/Processed_data/Global_IUCN_WHS_KBA_centroids.csv"
# world <- "/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/Large_data_files/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp"
# 
# 
# #-#-# Get the coordinates for the selected sites #-#-#
# sites <- read.csv(centroids)[2:4]
# colnames(sites) <- c("Int_Name", "x", "y")
# #selected <- Plot_map$Int_Name
# #selected_sites <- filter(sites, PAname %in% selected)
# selected_sites <- merge(Plot_map, sites, by = "Int_Name", all.x = T)
# head(selected_sites)
# 
# 
# #-#-# Load the world map #-#-#
# worldmap <- load_worldmap(world)
# 
# 
# #-#-# Load the realm map #-#-#
# realmmap <- get(load(realmpath))
# head(realmmap)
# realmmap <- na.omit(realmmap)
# realmmap <- subset(realmmap, !(RealmWWF == 2))
# realmmap <- subset(realmmap, !(RealmWWF == 7))
# 
# realmmap$RealmName <- 0
# realmmap$RealmName[realmmap$RealmWWF == 1] <- "Australasia"
# realmmap$RealmName[realmmap$RealmWWF == 3] <- "Afrotropic"
# realmmap$RealmName[realmmap$RealmWWF == 4] <- "Indomalaya"
# realmmap$RealmName[realmmap$RealmWWF == 5] <- "Nearctic"
# realmmap$RealmName[realmmap$RealmWWF == 6] <- "Neotropic"
# realmmap$RealmName[realmmap$RealmWWF == 8] <- "Palearctic"
# 
# 
# colours <- c("steelblue3","orange3","navajowhite3","darkolivegreen3","plum4","brown")
#  
# map <- ggplot(selected_sites) +
#   geom_raster(data = realmmap, aes( x = x, y = y, fill = as.factor(RealmName), alpha = 0.2)) +
#   scale_fill_manual("Realm",values = colours) +
#   geom_sf(data = worldmap, fill = NA, color = "grey40", size = 0.2) +
#   coord_sf(xlim = c(-170, 180), ylim = c(90, -60), expand = FALSE) +
#   geom_point(data = selected_sites,
#              aes(x = x, y = y, color = Scenario),
#              shape = 18,
#              size = 8,
#              alpha = 0.7) +
#   scale_color_manual(name = "Objective:",
#                      values = c("Biodiversity" = "darkgreen",
#                                 "Wilderness" = "red",
#                                 "Combination" = "black")) +
#   geom_segment(
#     aes(x = -180, xend = 180, y = 0, yend = 0), colour = "black", linetype = "dashed") +
#   theme(legend.key = element_blank(), legend.position = "bottom",
#         legend.title = element_text(size = 22, face = "bold"),
#         legend.text = element_text(size = 20)) +
#   theme(axis.text=element_text(size=20))+
#   theme(panel.background=element_rect(fill='white',colour="white"))+ # Remove the background
#   labs(x="", y="", title="") + # Remove axis titles
#   guides(colour = guide_legend(override.aes = list(size = 8)))
# 
# plot(map) 
# 
# #-#-# Save the final map #-#-#
# setwd("/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/Result_plots/")
# ggsave("Map Legacy Landscapes Conservation Objectives.tiff",map,width=25, height=15, unit="in", dpi=300, bg="white")
# 


#-#-# Map global but in black (after first feedback) #-#-#
#library(tidyverse)
library(sf)
library(rgdal)
library(packcircles)
#library(VennDiagram)
#library(ggVennDiagram)
#library(ggvenn)
library(colorfulVennPlot)
library(gridExtra)
library(png)
library(cowplot)
library(dplyr)

#-#-# Set the file paths #-#-#
realmpath <- "/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Site selection analysis/Other files/Realms and Regions/CMEC regions & realms/WWF_Realms_gridded_05.Rdata"
centroids <- "/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/LL_analysis/Processed_data/Global_IUCN_WHS_KBA_centroids.csv"
world <- "/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/Large_data_files/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp"


#-#-# Get the coordinates for the selected sites #-#-#
sites <- read.csv(centroids)[2:4]
colnames(sites) <- c("Int_Name", "x", "y")
selected_sites <- merge(Plot_map, sites, by = "Int_Name", all.x = T)

sites$sites <- "Included site"
sites$realm <- "Realm boarder"

#-#-# Summarize categories #-#-#
selected_sites <- subset(selected_sites, !Realm == "Global") # Remove globally best sites
head(selected_sites)

selected_sites <- selected_sites[order(factor(selected_sites$Scenario, levels = c("Biodiversity", "Combination", "Wilderness"))),]

dublicates <- selected_sites %>% 
  group_by(Int_Name) %>% 
  filter(n()>1)

sum_dub_names <- unique(dublicates$Int_Name)

double_selections <- lapply(sum_dub_names, function(d){
  one_PA <- subset(dublicates, Int_Name == d)
  Int_name <- as.data.frame(one_PA$Int_Name[[1]][1])
  Realm <- one_PA$Realm[1]
  x <- one_PA$x[1]
  y <- one_PA$y[1]
  
  scenario_names <- one_PA$Scenario
  
  combi_length <- nrow(one_PA)
  
  if(combi_length == 2){
  Scenario <- paste(scenario_names[1],"&",scenario_names[2])
  }
  
  if(combi_length == 3){
    Scenario <- paste(scenario_names[1],"&",scenario_names[2],"&",scenario_names[3])
  }
  
  sum_PA <- cbind(Int_name[1],Realm,Scenario,x,y)
  colnames(sum_PA) <- c("Int_Name","Realm","Scenario","x","y")
  
  return(sum_PA)
  
})

combined_hotspots <- do.call(rbind,double_selections)

singles <- setdiff(selected_sites$Int_Name,dublicates$Int_Name)
single_selections <- lapply(singles, function(s){
  print(s)
  single <- subset(selected_sites, Int_Name == s)
  return(single)
})
single_hotspots <- do.call(rbind,single_selections)

selected_sites_new <- rbind(single_hotspots[1:5],combined_hotspots)

#-#-# Space the points out to avoid overlay #-#-#
selected_sites_new$radius <- 5
head(selected_sites_new)

selected_sites_repel <- circleRepelLayout(
  selected_sites_new,
  xysizecols = c(4, 5, 6),
  sizetype = c("area"),
  maxiter = 1000,
  wrap = TRUE,
  weights = 1
)

selected_sites_repel_df <- as.data.frame(selected_sites_repel[[1]])
selected_sites_repel_df <- cbind(selected_sites_repel_df,selected_sites_new)
selected_sites_repel_df <- selected_sites_repel_df[1:6]

#-#-# Load the world map #-#-#
worldmap <- load_worldmap(world)


#-#-# Load the realm map #-#-#
setwd("/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Site selection analysis/Other files/Generalised_Biogeographic_Realms_2004/")
realm_poly <- readOGR("brpol_fin_dd.shp")
realm_st <- st_read("brpol_fin_dd.shp")
plot(realm_poly)
plot(realm_st)

realmmap <- get(load(realmpath))
realmmap <- na.omit(realmmap)
realmmap <- subset(realmmap, y > -60)

map <- ggplot() + 
  geom_raster(data = realmmap, aes( x = x, y = y), fill = "gray60") +
  geom_sf(data = realm_st, fill = NA, size = 0.6, color = "red") +
  geom_polygon(colour="black", fill='white') +
  geom_sf(data = worldmap, fill = NA, color = "black", size = 0.3) +
  coord_sf(xlim = c(-170, 179), ylim = c(-60, 89), expand = FALSE) +
  geom_point(data = sites,
             aes(x = x, y = y),
             shape = 16,
             size = 0.8,
             colour = "black") +
  geom_point(data = selected_sites_repel_df,
             aes(x = x, y = y, color = Scenario),
             shape = 18,
             size = 8,
             alpha = 1) +
  scale_color_manual(name = "Objective:",
                     values = c("Biodiversity" = "darkgreen",
                                "Wilderness" = "red",
                                "Combination" = "gold1",
                                "Biodiversity & Combination" = "sienna",
                                "Biodiversity & Wilderness" = "orange2",
                                "Combination & Wilderness" = "olivedrab3",
                                "Biodiversity & Combination & Wilderness" = "royalblue2")) +
  geom_segment(
  aes(x = -180, xend = 180, y = 0, yend = 0), colour = "black", linetype = "dashed") +
  theme(legend.key = element_blank(), legend.position = "none",
        legend.title = element_text(size = 22, face = "bold"),
        legend.text = element_text(size = 20)) +
  theme(axis.text=element_text(size=20))+
  theme(panel.background=element_rect(fill='white',colour="white"))+ # Remove the background
  labs(x="", y="", title="") + # Remove axis titles
  guides(colour = guide_legend(override.aes = list(size = 8))) 

plot(map) 


#-#-# Create the legends #-#-#
## Function to extract plot legends
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

## Included sites
points <- ggplot() + 
  geom_point(data = sites,
             aes(x = x, y = y, color = sites),
             shape = 16,
             size = 0.8) +
  scale_color_manual(name = "", values = c("Included site" = 'black')) +
  theme(legend.key = element_blank(),
        legend.title = element_text(size = 22, face = "bold"),
        legend.text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 8)))
plot(points) 

Site_leg <- g_legend(points)

## Realm boarder
realm <- ggplot() + 
  geom_point(data = sites,
             aes(x = x, y = y, color = realm),
             shape = 12,
             size = 0.8) +
  scale_color_manual(name = "", values = c("Realm boarder" = 'red')) +
  theme(legend.key = element_blank(),
        legend.title = element_text(size = 22, face = "bold"),
        legend.text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 8)))
plot(realm) 

Realm_leg <- g_legend(realm)

## Combine map and legend
Comb_leg <- arrangeGrob(Site_leg,Realm_leg,
                        widths = c(5,5),
                        heights = c(0.5),
                        ncol = 2,
                        nrow = 1)
plot(Comb_leg)

Comb_map <- arrangeGrob(map,Comb_leg,
                           widths = c(5),
                           heights = c(5,0.5),
                           ncol = 1,
                           nrow = 2)
plot(Comb_map)

## Venn diagram
Scenarioes <-  unique(selected_sites_new$Scenario)

Venn_data <- lapply(Scenarioes, function(v){
  All <- nrow(selected_sites_new)
  Subset <- as.data.frame(nrow(subset(selected_sites_new, Scenario == v)))
  #SubsetP <- round(Subset/(All/100),0) # if you need percentages
  #colnames(SubsetP) <- v
  return(Subset)
})

Venn_df <- do.call(cbind,Venn_data)

if(ncol(Venn_df) == 7){
  colnames(Venn_df) <- c("Biodiversity","Combination","Wilderness","Biodiversity & Combination",           
                         "Biodiversity & Wilderness","Biodiversity & Combination & Wilderness",
                         "Combination & Wilderness")
}

if(ncol(Venn_df) == 6){
  Venn_df$add <- 0
  colnames(Venn_df) <- c("Biodiversity", "Combination", "Wilderness", "Biodiversity & Combination", 
                         "Biodiversity & Wilderness", "Combination & Wilderness", "Biodiversity & Wilderness & Combination")
}

## Reorder DF in case there were only 6 combinations in the data
Venn_df <- Venn_df[c("Biodiversity","Combination","Wilderness","Biodiversity & Combination",           
                   "Biodiversity & Wilderness","Biodiversity & Combination & Wilderness",
                   "Combination & Wilderness")]

#colnames(Venn_df) <- c("111","011","101","001","110","010","100") ##original order Venn package
 
## Rename columns to their position in the venn diagram                      
colnames(Venn_df) <- c("100","001","010","101","110","111","011") #Shape order
Labels <- c("Biodiversity","Combination","Wilderness","Biodiversity & Combination",           
            "Biodiversity & Wilderness","Biodiversity & Combination & Wilderness",
            "Combination & Wilderness")

Venn_input <- list(Venn_df,Labels)

#library(patchwork)
Comb_p <- {map + plot_spacer() + plot_layout(ncol=2, nrow=1, heights=c(2)) + theme_classic()}
#Comb_p <- Comb_map + plot_spacer() + plot_layout(ncol=1, heights=c(1,1))

plot(Comb_p)

vp1 <- viewport(width=0.1, height=0.3, x = 0, y = 0, just=c("left","bottom"))
vp2 <- viewport(width=0.7, height=0.7, x = 0.5, y = 0.65, just=c("left", "bottom"))
png(file=paste0("/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/Result_plots/Venn_diagram.png"), 
    width=6, height=2, bg="transparent", units="in", res = 900)


plot.new()
  print(Comb_p,vp=vp1)
  # pushViewport(vp2)
  # plotVenn(Venn_input[[1]], Venn_input[[2]],
  #        Colors = c("darkgreen","gold1","red","olivedrab3","sienna","royalblue2","orange2"))
dev.off()






#-#-# Blank plot #-#-#  
blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank() +
    annotation_raster(img, ymin = 0,ymax= 5,xmin = 0,xmax = 5) 
  )
plot(blankPlot)
plot(img)

 


#-#-# Save the final map #-#-#
setwd("/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/Result_plots/")
ggsave("Map Legacy Landscapes Conservation Objectives_black.tiff",map,width=25, height=15, unit="in", dpi=300, bg="white")
