#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#    Principal component analysis of the conservation variables     # 
#               included in the "LL Site Selection Tool"            #
#                       Figure 2 Main manuscript                    #
#                         Alke November 2020                        #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#


#-#-# Clear memory #-#-#
rm(list=ls())


#-#-# Load libraries #-#-#
library(devtools) 
library(factoextra)
library(cowplot)
library(tidyverse)
library(ggplot2)
library(ggfortify)
library(grid)
library(gridExtra)


#-#-# Load the data frame containing values for a sites included in the analysis #-#-#
setwd("/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/LL_analysis/Processed_data/")
PAdata <- read.csv("Final_dataset_IUCN_WHS_KBA_country_correct_names.csv")


#-#-# Format the dataframe for the analysis #-#-#
## Subset to columns that are included in the LL site selection tool
PAdataS <- PAdata[c("SP_ID", "Int_Name", "PA_type", "IUCN_CAT", "COUNTRY", "ISO3", "RealmNr", "Area_Rep", "Area_GIS", "Prop_SR_B",
                    "Prop_SR_M", "Prop_SR_A", "Prop_SR_R", "RSR_B","RSR_M", "RSR_A", "RSR_R", "PE_WMean_B", "PE_WMean_M",
                    "PE_WMean_A", "PE_WMean_R", "TO_PA_weightbird", "TO_PA_weightmammal", "TO_PA_weightamphibian",
                    "TO_PA_weightreptile", "HFP_WMean", "BaseCarbon_WMean", "VulCarbon_WMean", "IrrCarbon_WMean", "BII_WMean",
                    "BioAnth_WMean", "BioAnth_Max", "TreeCover_WMean", "BioCropRisk", "CropRisk", "PastureRisk")]
head(PAdataS)
ncol(PAdataS)

## Add the realm names for analysis and plotting 
PAdataS <- subset(PAdataS, RealmNr == "1" | RealmNr == "3" | RealmNr == "4" | RealmNr == "5" | RealmNr == "6" | RealmNr == "8")
as.numeric(as.character(PAdataS$RealmNr))
PAdataS$RealmName <- 0
PAdataS$RealmName[PAdataS$RealmNr == 1] <- "Australasia"
PAdataS$RealmName[PAdataS$RealmNr == 3] <- "Afrotropic"
PAdataS$RealmName[PAdataS$RealmNr == 4] <- "Indomalaya"
PAdataS$RealmName[PAdataS$RealmNr == 5] <- "Nearctic"
PAdataS$RealmName[PAdataS$RealmNr == 6] <- "Neotropic"
PAdataS$RealmName[PAdataS$RealmNr == 8] <- "Palearctic"

## Select the base data that remains the same throughout the analysis (site info)
Base <- PAdataS[c(1:7,37)]
head(Base)

## Select the variables and ensure all columns are numeric
VariableData <- PAdataS[c(8:36)]
head(VariableData)
VariableData[] <- lapply(VariableData[c(1:ncol(VariableData))], function(x) as.numeric(as.character(x))) # Change columns to numeric

## Combine the two dataframes containing all formated variables needed for the analysis
PCAdata <- cbind(Base,VariableData)
head(PCAdata)
ncol(PCAdata)
nrow(PCAdata)


#-#-# Combine the sub variables for the analysis#-#-#
## Three land-use change variables as summed land-use change
PCAdata$Landuse_Stability <- rowSums(PCAdata[c("BioCropRisk","CropRisk","PastureRisk")],na.rm=T)

## Species richness of the individual taxa into as mean species richness 
PCAdata$Species_Richness  <- rowMeans(PCAdata[c("Prop_SR_B","Prop_SR_M","Prop_SR_A","Prop_SR_R")],na.rm=T)
head(PCAdata)

## Species endemism of the individual taxa as mean endemism
PCAdata$Species_Endemism  <- rowMeans(PCAdata[c("RSR_B","RSR_M","RSR_A","RSR_R")],na.rm=T)
head(PCAdata)

## Evolutionary diversity of the individual taxa as mean evolutionary diversity
PCAdata$Evolutionary_Diversity  <- rowMeans(PCAdata[c("PE_WMean_B","PE_WMean_M","PE_WMean_A","PE_WMean_R")],na.rm=T)
head(PCAdata)

## Species turnover of the individual taxa as mean turnover
PCAdata$Climatic_Stability  <- rowMeans(PCAdata[c("TO_PA_weightbird","TO_PA_weightmammal","TO_PA_weightamphibian","TO_PA_weightreptile")],na.rm=T)
head(PCAdata)


#-#-# Select variables used in the PCA rename and format for analysis #-#-#
## Select PCA variables
PCAdata <- PCAdata[c("SP_ID","Int_Name","RealmName","Species_Richness","Species_Endemism","Evolutionary_Diversity","BII_WMean",
                     "HFP_WMean","BioAnth_WMean","Climatic_Stability","TreeCover_WMean","Landuse_Stability","BaseCarbon_WMean","VulCarbon_WMean","IrrCarbon_WMean","Area_GIS")]


## Remove NAs for analysis 
PCAdata <- na.omit(PCAdata)
nrow(PCAdata)

## Turn climatic stability, land use stability, human footprint and biome to anthrome around!
## For these variables high values are negative for all other variables they are positive
PCAdata$HFP_WMean <- (PCAdata$HFP_WMean)* -1
PCAdata$Landuse_Stability <- (PCAdata$Landuse_Stability)* -1
PCAdata$Climatic_Stability <- (PCAdata$Climatic_Stability)* -1
PCAdata$BioAnth_WMean <- (PCAdata$BioAnth_WMean)* -1
PCAdata$TreeCover_WMean <- abs(PCAdata$TreeCover_WMean)* -1 #change to absolute values because both high values below and above 0 are negative (high risk) 

## Change the column names for the PCA plot 
colnames(PCAdata) <- c("SP_ID","Int_Name","Realm","High SR","High SE","High ED","High BII",
                       "Low HFP","Low RLC","High CS","Low TCC","High LUS","High CaS",
                       "High VCaS","High ICaS","Large size")



#-#-# Run the PCA globally and for all realms #-#-#
## Center and scale variables set to TRUE for PCA

#---# Global PCA #---#
## Scale variables and compute PCA
PA.pca.cl <- prcomp(PCAdata[c(4:16)],center=T,scale=T) 
summary(PA.pca.cl) 
str(PA.pca.cl)

## Explained variance of components 
fviz_eig(PA.pca.cl)

#---# Global PCA scatterplot coloured by realm for manuscript (Fig. 2) #---#
pal <- colorRampPalette(c("steelblue3","orange3","navajowhite3","darkolivegreen3","plum4","brown"))
colour <- pal(6)


# fviz_pca_var(Neotropic.pca, labelsize = 8, repel = TRUE, 
#              col.var=List, palette = colour) +


GlobalScatter <- fviz_pca_biplot(PA.pca.cl, pointsize = 4,
                              label="var", 
                              repel = TRUE, 
                              labelsize = 8,
                              habillage=PCAdata$Realm,
                              palette=colour,
                              col.var="black",
                              alpha = 0.6) +
  scale_shape_manual(values=c(18,18,18,18,18,18)) +
  coord_equal() +
  theme(legend.position = "none") +
  theme(text = element_text(size = 26),
        axis.title = element_text(size = 26),
        axis.text = element_text(size = 26)) +
  scale_x_continuous(breaks = c(-2,0,2,4)) +
  scale_y_continuous(breaks = c(-4,-2,0,4)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm")) +
  ggtitle("")
plot(GlobalScatter)


#-#_# Simple realm map #-#-#
## Set the file paths 
realmpath <- "/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Site selection analysis/Other files/Realms and Regions/CMEC regions & realms/WWF_Realms_gridded_05.Rdata"

realmmap <- get(load(realmpath))
head(realmmap)
realmmap <- na.omit(realmmap)
realmmap <- subset(realmmap, !(RealmWWF == 2))
realmmap <- subset(realmmap, !(RealmWWF == 7))

colours <- c("steelblue3","orange3","navajowhite3","darkolivegreen3","plum4","brown")

realms <- ggplot(realmmap) +
  geom_raster(aes( x = x, y = y, fill = as.factor(RealmWWF), alpha = 0.2)) +
  scale_fill_manual("Realm",values = colours) +
  coord_sf(xlim = c(-170, 180), ylim = c(90, -60), expand = FALSE) +
  theme(legend.position = "none") +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        plot.title = element_text(size = 26, hjust = 0.5))+
  theme(panel.background=element_rect(fill='white',colour="white")) + # Remove the background
  labs(x="", y="", title="") + # Remove axis titles
  ggtitle("Biogeographic realms")
plot(realms) 


#-#-# Arrange plot parts #-#-#
## Blank plot 
blankPlot <- ggplot()+
             geom_blank(aes(1,1))+
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
             axis.line = element_blank()) 
plot(blankPlot)

## Function to extract legend from plot 
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

## Extract legend from scatter plot (fist draw plot with legend)
GlobalScatterLeg <- fviz_pca_ind(PA.pca.cl, pointsize = 4,
                              label="none", 
                              habillage=PCAdata$Realm,
                              palette=colour,
                              alpha = 0.6,
                              legend.title = "") +
  scale_shape_manual(values=c(18,18,18,18,18,18)) +
  theme(legend.position = "right",
        legend.text=element_text(size=26)) 
plot(GlobalScatterLeg)

legScat <- g_legend(GlobalScatterLeg)
plot(legScat)


#-#-# Combine the PCA scatter plot and the realmmap #-#-#  
CombMap <- plot_grid(blankPlot,realms,legScat,blankPlot,ncol = 1)
CombScatter <- plot_grid(GlobalScatter,CombMap, rel_widths = c(2, 1))

plot(CombScatter)

#-#-# Save the final PCA scatter plot #-#-#
setwd("/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/Result_plots/")
ggsave("PCA Legacy Landscapes scatter realms resize rescale_FCCcorrected_country.tiff",CombScatter,width=20, height=15, unit="in", dpi=600, bg="white")



#---# Nearctic PCA #---#
Nearctic_data <- subset(PCAdata,Realm =="Nearctic")

## Scale variables and compute PCA
Nearctic.pca <- prcomp(Nearctic_data[c(4:16)],center=T,scale=T) 
summary(Nearctic.pca) 
str(Nearctic.pca)

#---# Palearctic PCA #---#
Palearctic_data <- subset(PCAdata,Realm =="Palearctic")

## Scale variables and compute PCA
Palearctic.pca <- prcomp(Palearctic_data[c(4:16)],center=T,scale=T) 
summary(Palearctic.pca) 
str(Palearctic.pca)

#---# Neotropic PCA #---#
Neotropic_data <- subset(PCAdata,Realm =="Neotropic")

## Scale variables and compute PCA
Neotropic.pca <- prcomp(Neotropic_data[c(4:16)],center=T,scale=T) 
summary(Neotropic.pca) 
str(Neotropic.pca)

#---# Afrotropic PCA #---#
Afrotropic_data <- subset(PCAdata,Realm =="Afrotropic")

## Scale variables and compute PCA
Afrotropic.pca <- prcomp(Afrotropic_data[c(4:16)],center=T,scale=T) 
summary(Afrotropic.pca) 
str(Afrotropic.pca)

#---# Australasia PCA #---#
Australasia_data <- subset(PCAdata,Realm =="Australasia")

## Scale variables and compute PCA
Australasia.pca <- prcomp(Australasia_data[c(4:16)],center=T,scale=T) 
summary(Australasia.pca) 
str(Australasia.pca)

#---# Indomalaya PCA #---#
Indomalaya_data <- subset(PCAdata,Realm =="Indomalaya")

## Scale variables and compute PCA
Indomalaya.pca <- prcomp(Indomalaya_data[c(4:16)],center=T,scale=T) 
summary(Indomalaya.pca) 
str(Indomalaya.pca)


#-#-# Combined figure for all PCA results #-#-#
## PCA axis plots (dimension 1 and 2) 

## Set colours for the "variables" (colour loadings of the PCA) by variables
palVar <- colorRampPalette(c("tomato2", "dodgerblue4", "orange", "chartreuse3", "skyblue1",
                          "peru", "forestgreen", "darkgreen", "dodgerblue2",
                          "saddlebrown", "gold", "hotpink1", "lightpink1"))
colourVar <- palVar(13)

## List variable names to colour plot by
List <- c("High SR","High SE","High ED","High BII",
          "Low HFP","Low RLC","High CS","Low FCC","High LUS","High CaS",
          "High VCaS","High ICaS","Large size")


## Global PCA plot with coloured loadings
Global <- fviz_pca_var(PA.pca.cl, labelsize = 8, repel = TRUE,
                       col.var=List, palette = colourVar) +
          theme(axis.title = element_text(size = 24),
          axis.text = element_text(size = 20),
          plot.title = element_text(size = 26),
          legend.position = "none") +
          ggtitle("a) Global")
plot(Global)

## Nearctic PCA plot with coloured loadings
Nearctic <- fviz_pca_var(Nearctic.pca, labelsize = 8, repel = TRUE,
                         col.var=List, palette = colourVar) +
        theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 26),
        legend.position = "none") +
  ggtitle("b) Nearctic")
plot(Nearctic)

## Palearctic PCA plot with coloured loadings
Palearctic <- fviz_pca_var(Palearctic.pca, labelsize = 8, repel = TRUE,
                           col.var=List, palette = colourVar) +
        theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 26),
        legend.position = "none") +
  ggtitle("c) Palearctic")
plot(Palearctic)

## Neotropic PCA plot with coloured loadings
Neotropic <- fviz_pca_var(Neotropic.pca, labelsize = 8, repel = TRUE, 
                          col.var=List, palette = colourVar) +
        theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 26),
        legend.position = "none") +
  ggtitle("e) Neotropic")
plot(Neotropic)

## Afrotropic PCA plot with coloured loadings
Afrotropic <- fviz_pca_var(Afrotropic.pca, labelsize = 8, repel = TRUE,
                           col.var=List, palette = colourVar) +
        theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 26),
        legend.position = "none") +
  ggtitle("f) Afrotropic")
plot(Afrotropic)

## Australasia PCA plot with coloured loadings
Australasia <- fviz_pca_var(Australasia.pca, labelsize = 8, repel = TRUE,
                            col.var=List, palette = colourVar) +
        theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 26),
        legend.position = "none") +
  ggtitle("g) Australasia")
plot(Australasia)

## Indomalaya PCA plot with coloured loadings
Indomalaya <- fviz_pca_var(Indomalaya.pca, labelsize = 8, repel = TRUE,
                           col.var=List, palette = colourVar) +
        theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 26),
        legend.position = "none") +
  ggtitle("d) Indomalaya")
plot(Indomalaya)


## Create legend for the PCA plot
## Function to extract plot legends
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

## Create simple dataframe containing the list of variable names to be displayed in legend
#c("High SR","High SE","High ED","High BII","Low HFP","Low BTAC","High CS",
#  "Low FCC","High LUS","High CaS","High VCaS","High ICaS","Large size"))
Variables <- c("High species richness (SR)","High species endemism (SE)","High evolutionary diversity (ED)","High biodiversity intactness (BII)",
          "Low human footprint (HFP)","Low recent land-use change (RLC)","High climatic stability (CS)","Low tree cover change (TCC)",
          "High land-use stability (LUS)","High manageable carbon (CaS)","High vulnerable carbon (VCaS)","High irreplaceable carbon (ICaS)","Large size")
Vals <- c(rep(1,13))
ValsII <- c(rep(1,13))
Fake_data <- as.data.frame(cbind(Variables,Vals,ValsII))
Fake_data$Variables <- factor(Fake_data$Variables, levels = c("High species richness (SR)","High species endemism (SE)","High evolutionary diversity (ED)","High biodiversity intactness (BII)",
                                                              "Low human footprint (HFP)","Low recent land-use change (RLC)","High climatic stability (CS)","Low tree cover change (TCC)","High land-use stability (LUS)",
                                                              "High manageable carbon (CaS)","High vulnerable carbon (VCaS)","High irreplaceable carbon (ICaS)","Large size"))

## ggplot to extract legend from for final plot
In <- ggplot(data = Fake_data, aes(Vals, ValsII, color = Variables)) +
  geom_point() +
  scale_color_manual(values = c("High species richness (SR)" = "darkgreen",
                                "High species endemism (SE)" = "forestgreen",
                                "High evolutionary diversity (ED)" = "chartreuse3",
                                "High biodiversity intactness (BII)" = "tomato2",
                                "Low human footprint (HFP)" = "hotpink1",
                                "Low recent land-use change (RLC)" = "lightpink1",
                                "High manageable carbon (CaS)" = "dodgerblue4",
                                "High vulnerable carbon (VCaS)" = "dodgerblue2",
                                "High irreplaceable carbon (ICaS)" = "skyblue1",
                                "Large size" = "saddlebrown",
                                "High land-use stability (LUS)" = "peru",
                                "High climatic stability (CS)" = "orange",
                                "Low tree cover change (TCC)" = "gold")) +
  guides(color = guide_legend(override.aes = list(linetype = 0, size=8))) +
  theme(legend.key.size = unit(1,"point"), #change legend key size
        #legend.key.height = unit(1, 'cm'), #change legend key height
        #legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=24), #change legend title font size
        legend.text = element_text(size=22)) #change legend text font size
  
ScatterLegend <- g_legend(In)
plot(ScatterLegend)


#-#-# Combine the PCA plots and the legend #-#-#  
CombGlobal <- arrangeGrob(Global,Nearctic,Palearctic,Indomalaya,ScatterLegend,Neotropic,Afrotropic,Australasia,
                        widths = c(1,1,1,1),
                        heights = c(1,1),
                        ncol = 4,
                        nrow = 2)
plot(CombGlobal)

CombGlobal <- plot_grid(Global,Nearctic,Palearctic,Indomalaya,ScatterLegend,Neotropic,Afrotropic,Australasia,
                          ncol = 4, nrow = 2)
plot(CombGlobal)


#-#-# Save the final PCA plot #-#-#
setwd("/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/Result_plots/")
ggsave("PCA Legacy Landscapes Conservation Objectives_full_names_Abrev_labelled_FCCcorrected_final_country.tiff",CombGlobal,width=25, height=15, unit="in", dpi=300, bg="white")


#--------------------     # Supplementary plot - Scree plots #---------------------------------#
#-#_# Simple realm map #-#-#
## Set the file paths 
realmpath <- "/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Site selection analysis/Other files/Realms and Regions/CMEC regions & realms/WWF_Realms_gridded_05.Rdata"

realmmap <- get(load(realmpath))
head(realmmap)
realmmap <- na.omit(realmmap)
realmmap <- subset(realmmap, !(RealmWWF == 2))
realmmap <- subset(realmmap, !(RealmWWF == 7))

colours <- c("steelblue3","orange3","navajowhite3","darkolivegreen3","plum4","brown")

realms <- ggplot(realmmap) +
  geom_raster(aes( x = x, y = y, fill = as.factor(RealmWWF), alpha = 0.2)) +
  scale_fill_manual("Realm",values = colours) +
  coord_sf(xlim = c(-170, 180), ylim = c(90, -60), expand = FALSE) +
  theme(legend.position = "none") +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        plot.title = element_text(size = 25))+
  theme(panel.background=element_rect(fill='white',colour="white"), # Remove the background
        plot.title = element_text(size = 22, hjust = 0.5)) + 
  labs(x="", y="", title="") + # Remove axis titles
  guides(colour = guide_legend(override.aes = list(size = 8))) +
  ggtitle("Biogeographic realms")
plot(realms) 

Glob_scree <- fviz_eig(PA.pca.cl, choice = "variance",  #"eigenvalue"
                       geom = c("line"),
                       ggtheme = theme_minimal()) +
                       theme(axis.title = element_blank(),
                             axis.text = element_text(size = 16),
                             plot.title = element_text(size = 22, face="bold")) +
                         ggtitle("Global")
plot(Glob_scree)

Nearctic_scree <- fviz_eig(Nearctic.pca, choice = "variance",
                       geom = c("line"),
                       linecolor = "darkolivegreen3",
                       ggtheme = theme_minimal()) +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 22, face="bold")) +
  ggtitle("Nearctic")
plot(Nearctic_scree)

Palearctic_scree <- fviz_eig(Palearctic.pca, choice = "variance",  
                           geom = c("line"),
                           linecolor = "brown",
                           ggtheme = theme_minimal()) +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 22, face="bold")) +
  ggtitle("Palearctic")
plot(Palearctic_scree)

Indomalaya_scree <- fviz_eig(Indomalaya.pca, choice = "variance", 
                             geom = c("line"),
                             linecolor = "navajowhite3",
                             ggtheme = theme_minimal()) +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 22, face="bold")) +
  ggtitle("Indomalaya")
plot(Indomalaya_scree)

Neotropic_scree <- fviz_eig(Neotropic.pca, choice = "variance", 
                             geom = c("line"),
                             linecolor = "plum4",
                             ggtheme = theme_minimal()) +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 22, face="bold")) +
  ggtitle("Neotropic")
plot(Neotropic_scree)

Afrotropic_scree <- fviz_eig(Afrotropic.pca, choice = "variance", 
                            geom = c("line"),
                            linecolor = "orange3",
                            ggtheme = theme_minimal()) +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 22, face="bold")) +
  ggtitle("Afrotropic")
plot(Afrotropic_scree)

Australasia_scree <- fviz_eig(Australasia.pca, choice = "variance",
                             geom = c("line"),
                             linecolor = "steelblue3",
                             ggtheme = theme_minimal()) +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 16),
        plot.title = element_text(size = 22, face="bold")) +
  ggtitle("Australasia")
plot(Australasia_scree)

#Blank <- ggplot() + theme_void() # Use blank plot instead of realm map

#-#-# Combine the PCA plots and the legend #-#-#  
CombScree <- arrangeGrob(Glob_scree,Nearctic_scree,Palearctic_scree,Indomalaya_scree,realms,Neotropic_scree,Afrotropic_scree,Australasia_scree,
                          widths = c(1,1,1,1),
                          heights = c(1,1),
                          ncol = 4,
                          nrow = 2,
                          bottom = textGrob("Dimensions", hjust = 1, gp = gpar(fontface = "bold", cex = 2)),
                          left = textGrob("Percentage of explained variances", rot = 90, vjust = 1, gp = gpar(fontface = "bold", cex = 2)))
plot(CombScree)


#-#-# Save the final PCA plot #-#-#
setwd("/Users/alkevoskamp/AG BGFM Dropbox/Voskamp/Legacy Landscapes/Legacy_landscapes_analysis/Result_plots/")
ggsave("Spp Legacy Landscapes PCA scree plot colour_country.tiff",CombScree,width=25, height=12, unit="in", dpi=300, bg="white")
