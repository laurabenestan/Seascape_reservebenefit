### Download libraries
library(mapproj)
library(raster)
library(gpclib)
library(ggplot2)
library(shapefiles)
library(ggmap)
library(sp)
library(RgoogleMaps)
library(dplyr)
library(cowplot)
library(rgdal)
library(mapdata)
library(wesanderson)

############### DOWNLOAD GEOGRAPHIC AND MPA INFORMATION #######################

# Import geographic coordinates file
data<-read.table("../00-data/distances_to_mpa_no_coteblue.txt",header=TRUE, dec=".",sep="\t",na.strings="NA",strip.white=T)
summary(data) 

# Add species information
data$POP <- substr(data$IND,1,3)
unique(data$POP)
data$SPECIES <- ifelse(data$POP =="Mul"| data$POP =="mul","Mullus surmuletus", ifelse(data$POP=="dip","Diplodus sargus", ifelse(data$POP=="ser","Serranus cabrilla","Palinurus elephas")))

### Explore how many individuals per MPA
table(data$MPA.NAME_EN,data$SPECIES) 

### Remove Cote Bleue
data2 <- data[data$MPA.NAME_EN!="Cote Bleue",]

### Remove palinurus
data3 <- data2[data2$SPECIES!="Palinurus elephas",]

### Create three subsets for each species
diplodus <- subset(data3, data3$SPECIES=="Diplodus sargus")
mullus <- subset(data3, data3$SPECIES=="Mullus surmuletus")
serranus <- subset(data3, data3$SPECIES=="Serranus cabrilla")

### Download the map for the Mediterranean Sea
wH <- map_data("worldHires",  xlim=c(-8,37), ylim=c(29.5,47)) # subset polygons surrounding med sea

### Change order of MPA for latitude graidnet on the map
data3$MPA.NAME_EN <- factor(data3$MPA.NAME_EN, levels = c("Cabo de Gata Níjar", "Cabo de Palos", "Illa De Tabarca","Illes Columbretes","Llevant De Mallorca","Norte de Menorca","Cap de Creus","Cerbere-Banyuls"))

### Map with MPA coloured
graph1 <- ggplot() +
  geom_polygon(data = wH, aes(x=long, y = lat, group = group), fill = "gray80", color = NA) +
  coord_fixed(xlim = c(-6, 6), ylim = c(35, 45), ratio=1.2)+
  geom_point(aes(x = LON, y = LAT, fill=MPA.NAME_EN,shape=SPECIES), size=2,data =data3, size=1.5)+
  scale_shape_manual(values=c(21, 23, 25))+
  scale_fill_brewer(palette="RdYlBu", name="Marine reserves") +
  xlab("Longitude") + 
  ylab("Latitude") +
  theme_classic()+
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "black")))
graph1
ggsave("Fig1.pdf", width=8, height=8)
