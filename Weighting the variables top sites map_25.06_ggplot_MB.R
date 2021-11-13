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
library(png)
library(cowplot)
library(dplyr)

# Additional packages
library(ggpubr)

#-#-# Set the file paths #-#-#
realmpath <- "/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Site selection analysis/Other files/Realms and Regions/CMEC regions & realms/WWF_Realms_gridded_05.Rdata"
centroids <- "/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/LL_analysis/Processed_data/Global_IUCN_WHS_KBA_centroids.csv"
world <- "/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/Large_data_files/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp"


#-#-# Get the coordinates for the selected sites #-#-#
sites <- read.csv(centroids)[2:4]
colnames(sites) <- c("Int_Name", "x", "y")
selected_sites <- merge(Plot_map, sites, by = "Int_Name", all.x = T)

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


#-#-# Plot the global hotspots map #-#-#
map <- ggplot() + 
  geom_raster(data = realmmap, aes( x = x, y = y), fill = "black") +
  geom_sf(data = realm_st, fill = NA, size = 0.6, color = "red") +
  geom_polygon(colour="black", fill='white') +
  geom_sf(data = worldmap, fill = NA, color = "gray60", size = 0.3) +
  coord_sf(xlim = c(-170, 179), ylim = c(-60, 89), expand = FALSE) +
  geom_point(data = sites,
             aes(x = x, y = y),
             shape = 16,
             size = 0.6,
             colour = "white") +
  geom_point(data = selected_sites_repel_df,
             aes(x = x, y = y, color = Scenario),
             shape = 18,
             size = 4,
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
  aes(x = -180, xend = 180, y = 0, yend = 0), colour = "gray60", linetype = "dashed") +
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
        legend.text = element_text(size = 25)) +
  guides(colour = guide_legend(override.aes = list(size = 8)))
plot(points) 

Site_leg <- g_legend(points)

Site_leg2 <- ggpubr::as_ggplot(ggpubr::get_legend(points))
Site_leg2

## Realm boarder
realm <- ggplot() + 
  geom_point(data = sites,
             aes(x = x, y = y, color = realm),
             shape = 12,
             size = 0.8) +
  scale_color_manual(name = "", values = c("Realm boarder" = 'red')) +
  theme(legend.key = element_blank(),
        legend.title = element_text(size = 22, face = "bold"),
        legend.text = element_text(size = 25)) +
  guides(colour = guide_legend(override.aes = list(size = 8)))
plot(realm) 

Realm_leg <- g_legend(realm)

Realm_leg2 <- ggpubr::as_ggplot(ggpubr::get_legend(realm))
Realm_leg2

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
                   "Wilderness" = c(rep(Labels[2], Venn_df[1,2]), 
                                    rep(Labels[4], Venn_df[1,4]),
                                    rep(Labels[6], Venn_df[1,6]),
                                    rep(Labels[7],Venn_df[1,7])),
                   "Combination" = c(rep(Labels[3], Venn_df[1,3]), 
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
  geom_sf_text(aes(label = name), data = venn_setlabel(ggvenn_dat), size=10) +
  geom_sf_label(aes(label = count), fontface = "bold", data = venn_region(ggvenn_dat), size=8) +
  scale_fill_manual(values=c("darkgreen", "sienna", "orange2", "royalblue2", "gold1", "red", "olivedrab3")) + 
  theme_void() + theme(legend.position="none") +
  xlim(c(-8,20)) #+ylim(c(-20,20))
venn_plot

# Always need to re-check colours are correctly ordered, if code changes again!!!

# Combine Figures

## Combine map and legend
library(patchwork)
mapLg <- map + theme(legend.position = "right")
FinalPlot <- (map + venn_plot + plot_layout(widths=c(3,2))) / (Site_leg2 + Realm_leg2) + 
              plot_layout(heights=c(8,1))

#-#-# Save the final map #-#-#
setwd("/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/Result_plots/")
ggsave("Map Legacy Landscapes Conservation Objectives_black.tiff",FinalPlot,width=20, height=6, unit="in", dpi=300, bg="white")

