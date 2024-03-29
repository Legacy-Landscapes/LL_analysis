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


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# Weighting the different objectives #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## Top sites Biodiversity, Ecosystem integrity Combined (Bio25%, Wild 25%, Size25%, LandUseStab10%, ClimateStab10%, CarbonStorage5%)
Filepath <- "/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/LL_analysis/Processed_data/Scaled_variables_country/"
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
    WeightTable05 <- WeightTable[1:5,] ## Change here for the number of sites included
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

## Bind top sites into one dataframe to plot the hotspot map
Map_table <- do.call(rbind,Ranks)
Plot_map <- Map_table[c("Int_Name", "Scenario", "Realm")]
head(Plot_map)

## Save the ranking 
#write.csv(WeightTable,"Ranking_indomalaya_4Bio_2Wild_05ClimSt_05LUStab_2Size_1ClimP.csv")


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# Plot the hotspot map #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

#-#-# Load the plot libraries #-#-#
library(sf)
library(rgdal)
library(packcircles)
library(colorfulVennPlot)
library(gridExtra)
#library(png)
library(cowplot)
library(dplyr)
library(raster)

# Additional packages
library(ggpubr)
library(geosphere)

#-#-# Set the file paths #-#-#
countrypath <- "/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/Large_data_files/ne_50m_admin_0_countries/"
realmpath <- "/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/Large_data_files/Olson_realms/"
realmgrid <- "/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Site selection analysis/Other files/Realms and Regions/CMEC regions & realms/"
centroids <- "/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/LL_analysis/Processed_data/Global_IUCN_WHS_KBA_centroids.csv"
world <- "/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/Large_data_files/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp"


# #-#-# Load shapefile function #-#-#
# load_worldmap <- function(filename) {
#   worldmap <- sf::st_read(filename, layer = "ne_50m_admin_0_countries")
#   #worldmap <- simplify_polygons(worldmap)
#   return(worldmap)
# }

#-#-# Get the coordinates for the selected sites #-#-#
sites <- read.csv(centroids)
sites <- sites[c("Int_Name", "x", "y")]
selected_sites <- merge(Plot_map, sites, by = "Int_Name", all.x = T)
# setwd("/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/Test")
# write.csv(Plot_map,"Plot_map.csv")
# write.csv(selected_sites,"Selected_sites.csv")

## Add columns for legend
sites$sites <- "Included site"
sites$realm <- "Realm boarder"


#-#-# Summarize categories #-#-#
selected_sites <- subset(selected_sites, !Realm == "Global") # Remove globally best sites
head(selected_sites)

selected_sites <- selected_sites[order(factor(selected_sites$Scenario, levels = c("Biodiversity", "Combination", "Wilderness"))),]

## Select sites that are hotspots for more scenarioes
dublicates <- selected_sites %>% 
  group_by(Int_Name) %>% 
  filter(n()>1)

sum_dub_names <- unique(dublicates$Int_Name)

## Summarize that are selected under multiple scenarioes
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

## Select the sites that are only hotspots under one scenario
singles <- setdiff(selected_sites$Int_Name,dublicates$Int_Name)
single_selections <- lapply(singles, function(s){
  print(s)
  single <- subset(selected_sites, Int_Name == s)
  return(single)
})
single_hotspots <- do.call(rbind,single_selections)

## New dataframe with the single and the summarized combined hotspots 
selected_sites_new <- rbind(single_hotspots[1:5],combined_hotspots)


#-#-# Space the points out to avoid overlay #-#-#
selected_sites_new$radius <- 15 #5
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

#-#-# Transfer the points to Robinson projection #-#-#
selected_sites_repel_df <- st_as_sf(selected_sites_repel_df, coords = c("x", "y"), crs = 4326)
selected_sites_repel_df <- st_transform(selected_sites_repel_df, crs = "+proj=robin +over")

sites <- st_as_sf(sites, coords = c("x", "y"), crs = 4326)
sites <- st_transform(sites, crs = "+proj=robin +over")


#-#-# Load the world map #-#-#
setwd(countrypath)
worldmap <- st_read("ne_50m_admin_0_countries.shp")
worldmap <- st_transform(worldmap, "ESRI:54030")


#-#-# Load the realm map #-#-#
# setwd("/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Site selection analysis/Other files/Generalised_Biogeographic_Realms_2004/")
# #realm_poly <- readOGR("brpol_fin_dd.shp")
# realm_st <- st_read("brpol_fin_dd.shp")
# #shp <- readOGR(dsn = file.path(tempdir(), "GBR_adm1.shp"), stringsAsFactors = F)
# plot(realm_st)
# #plot(realm_st)

setwd(realmgrid)
realmmap <- get(load("WWF_Realms_gridded_05.Rdata"))
#realmmap <- na.omit(realmmap)
realmmap <- subset(realmmap, y > -60)
realmmap <- na.omit(realmmap)
realmmap <- subset(realmmap, !(RealmWWF == 2))
realmmap <- subset(realmmap, !(RealmWWF == 7))

# realmmapN <- realmmap
# realmmapN$RealmName <- 0
# realmmapN$RealmName[realmmap$RealmWWF == 1] <- "Australasia"
# realmmapN$RealmName[realmmap$RealmWWF == 3] <- "Afrotropic"
# realmmapN$RealmName[realmmap$RealmWWF == 4] <- "Indomalaya"
# realmmapN$RealmName[realmmap$RealmWWF == 5] <- "Nearctic"
# realmmapN$RealmName[realmmap$RealmWWF == 6] <- "Neotropic"
# realmmapN$RealmName[realmmap$RealmWWF == 8] <- "Palearctic"

realmraster<- rasterFromXYZ(realmmap)
crs(realmraster) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
realmraster <- projectRaster(realmraster, crs = "+proj=robin +over")
plot(realmraster)

realmraster_df <- as.data.frame(realmraster, xy = TRUE) 
realmraster_df$RealmWWF <- round(realmraster_df$RealmWWF,0)
realmraster_df[is.na(realmraster_df)] <- 0

realmraster_df$RealmName <- 0
realmraster_df$RealmName[realmraster_df$RealmWWF == 1] <- "Australasia"
realmraster_df$RealmName[realmraster_df$RealmWWF == 3] <- "Afrotropic"
realmraster_df$RealmName[realmraster_df$RealmWWF == 4] <- "Indomalaya"
realmraster_df$RealmName[realmraster_df$RealmWWF == 5] <- "Nearctic"
realmraster_df$RealmName[realmraster_df$RealmWWF == 6] <- "Neotropic"
realmraster_df$RealmName[realmraster_df$RealmWWF == 8] <- "Palearctic"

setwd(realmpath)
realm_WWF_st <- st_read("wwf_terr_ecos.shp") 



realm_colours <- c("steelblue3","orange3","navajowhite3","darkolivegreen3","plum4","brown")


#-#-# Plot the global hotspots map #-#-#
map <- ggplot() + 
  geom_raster(data = realmraster_df, aes( x = x, y = y, fill = as.factor(RealmName))) +
  scale_fill_manual(name = "Realm",
                     values = alpha(c("Nearctic" = "darkolivegreen3",
                                "Palearctic" = "brown",
                                "Indomalaya" = "navajowhite3",
                                "Neotropic" = "plum4",
                                "Afrotropic" = "orange3",
                                "Australasia" = "steelblue3",
                                "0" = "white"), 0.2)) +
  geom_sf(data = worldmap, fill = NA, color = "black", size = 0.2) +
  geom_sf(data = sites,
             shape = 16,
             size = 0.6,
             colour = "slategrey") +
  geom_sf(data = selected_sites_repel_df,
             aes(colour = Scenario),
             shape = 17,
             size = 4,
             alpha = 1) +
  coord_sf(xlim = c(-17262571, 17275829),ylim = c(-6137346, 8575154), datum = sf::st_crs("ESRI:54030")) +
  scale_color_manual(name = "Objective:",
                     values = c("Biodiversity" = "darkgreen",
                                "Wilderness" = "red",
                                "Combination" = "royalblue2",
                                "Biodiversity & Combination" = "sienna",
                                "Biodiversity & Wilderness" = "orange2",
                                "Combination & Wilderness" = "olivedrab3",
                                "Biodiversity & Combination & Wilderness" = "gold1")) +
  geom_segment(aes(x = -170, xend = 180, y = 0, yend = 0), colour = "gray60", linetype = "dashed") +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) 

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
  geom_sf(data = sites,
             aes(color = sites),
             shape = 16,
             size = 0.8) +
  scale_color_manual(name = "", values = c("Included site" = 'slategrey')) +
  theme(legend.key = element_blank(),
        legend.title = element_text(size = 22, face = "bold"),
        legend.text = element_text(size = 25)) +
  guides(colour = guide_legend(override.aes = list(size = 8)))
plot(points) 

Site_leg <- g_legend(points)

Site_leg2 <- ggpubr::as_ggplot(ggpubr::get_legend(points))
Site_leg2

Comb_map <- arrangeGrob(map,Site_leg2,
                        widths = c(5),
                        heights = c(5,0.5),
                        ncol = 1,
                        nrow = 2)
plot(Comb_map)


#-#-# Venn diagram for the hotspot legend #-#-#
Scenarioes <-  unique(selected_sites_new$Scenario)

## Format the input data for the Venn diagram
Venn_data <- lapply(Scenarioes, function(v){
  All <- nrow(selected_sites_new)
  Subset <- as.data.frame(nrow(subset(selected_sites_new, Scenario == v)))
  #SubsetP <- round(Subset/(All/100),0) # if you need percentages
  #colnames(SubsetP) <- v
  return(Subset)
})

Venn_df <- do.call(cbind,Venn_data)

## Rename the columns 
if(ncol(Venn_df) == 7){
  colnames(Venn_df) <- Scenarioes
} else if(ncol(Venn_df) == 6){
  Venn_df$add <- 0
  colnames(Venn_df) <- c("Biodiversity", "Combination", "Wilderness", "Biodiversity & Combination", 
                         "Biodiversity & Wilderness", "Combination & Wilderness", 
                         "Biodiversity & Wilderness & Combination")
}
## Reorder DF
Venn_df <- Venn_df[c("Biodiversity", "Wilderness", "Combination", 
                     "Biodiversity & Wilderness", "Biodiversity & Combination",
                     "Combination & Wilderness", "Biodiversity & Combination & Wilderness")]
Venn_df
Labels <- colnames(Venn_df)
Labels

# Now, the values show up correctly

ggvenn_dat <- list("Biodiversity"= c(rep(Labels[1], Venn_df[1,1]), 
                                     rep(Labels[4], Venn_df[1,4]),
                                     rep(Labels[5],Venn_df[1,5]),
                                     rep(Labels[7],Venn_df[1,7])),
                   "Ecosystem integrity"= c(rep(Labels[2], Venn_df[1,2]), 
                                    rep(Labels[4], Venn_df[1,4]),
                                    rep(Labels[6], Venn_df[1,6]),
                                    rep(Labels[7],Venn_df[1,7])),
                   "LLF combination" = c(rep(Labels[3], Venn_df[1,3]), 
                                     rep(Labels[5], Venn_df[1,5]),
                                     rep(Labels[6], Venn_df[1,6]),
                                     rep(Labels[7], Venn_df[1,7])))
ggvenn_dat

library(ggVennDiagram)
venn <- Venn(ggvenn_dat)
ggvenn_dat <- process_data(venn)
ggvenn_dat@region$name
ggvenn_dat@region$count <- as.numeric(Venn_df[1,])
venn_plot <- ggplot() +
  geom_sf(aes(fill=name), data = venn_region(ggvenn_dat)) +
  geom_sf(size = 0.25, color = "black", data = venn_setedge(ggvenn_dat), show.legend = F) +
  geom_sf_text(aes(label = name), data = venn_setlabel(ggvenn_dat), size=8) +
  geom_sf_label(aes(label = count), fontface = "bold", data = venn_region(ggvenn_dat), size=7) +
  scale_fill_manual(values=c("darkgreen", "orange2", "gold1", "sienna", "red", "olivedrab3", "royalblue2")) + 
  theme_void() + theme(legend.position="none") +
  xlim(c(-8,20)) #+ylim(c(-20,20))
venn_plot

# Always need to re-check colours are correctly ordered, if code changes again!!!

# Combine Figures

## Combine map and legend
library(patchwork)
mapLg <- map + theme(legend.position = "right")
FinalPlot <- (map + venn_plot + plot_layout(widths=c(3,2))) + #/ (Site_leg2) + 
  plot_layout(heights=c(3))

plot(FinalPlot)

#-#-# Save the final map #-#-#
setwd("/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/Result_plots/")
ggsave("Map Legacy Landscapes Conservation Objectives_robinson_country.tiff",FinalPlot,width=20, height=6, unit="in", dpi=300, bg="white")
