###########################################################################################################
# Distance-based Redundancy Analysis (db-RDA)
#
# Author : Laura Benestan
# Date : 27-07-2022
# 
#Input:
#Euclidean distances or genepop file = response variable
#Environmental table = explanatory variables
#Order of the samples 

#Libraries that we will need
library(adespatial)
library(ade4)
library(vcfR)
library(codep)
library(adespatial)
library(adegraphics)
library(vegan)
library(ape)
library(car)
library(adegenet)
library(dplyr)
library(data.table)
library(RColorBrewer)
library(tidyr)
library(fastDummies)
library(ggplot2)
library(factoextra)
library(PCPS)
library(spdep)

##### Y = GENOMIC DATA ####

### Genomic data to be used

#Import vcf file
vcf <- read.vcfR("../../../00-Data/00-All/13101snps_468ind.recode.vcf")

#Transform vcf to genind file
genind <- vcfR2genind(vcf)

### Calculate Euclidean distances
distgenEUCL <- dist(genind, method = 
                      "euclidean", diag = FALSE, upper = FALSE, p = 2)

##### X = PCA ENV ####

### Download the environmental matrix
env = read.table('../00-data/all_env_data_4species.txt', row.names=1, header=T)
EnvMatrix <- env[2:ncol(env)]

### Scale before PCA
sEnvMatrix = scale(EnvMatrix, center = T, scale = T)

# We can use the prcomp() function to perform the principal component analysis.
ePCA = prcomp(sEnvMatrix)

# The principal components can be found in the $x matrix:
pca_axis <- as.data.frame(ePCA$x)
fviz_eig(ePCA)

#Subset the variables Env file
env_ser <- pca_axis[588:1055,]

#### X = HABITATS ####

### Import habitats variables
habitat <- read.table("../../../00-Data/05-Habitats/1299ind_habitats.txt", header=TRUE, sep="\t",dec=".")
habitat_select <- habitat[588:1055,]

### Produce dummy variables
habitats_dummy <- fastDummies::dummy_cols(habitat_select, select_columns = "Habitats",remove_first_dummy = TRUE)

##### ORDINATION ON Y ####
Pcoa=pcoa(distgenEUCL)
Pcoa

### Check the percent of variation explained by each axis
Pcoa$values$Cumul_eig

### Keep only significant Pcoa principal components, which will be the response variable in the db-RDA
X=as.data.frame(Pcoa$vectors)

#Look at genotypes distribution in relation to the first 2 Pcoa axes
plot(X[,1], X[,2])

############### SELECT ENV ####

#OrdiR2step will start working from an empty model without explanatory variables, just the intercept
rda0<-rda(X ~ 1, env_ser[,1:3])

#OrdiR2step will move towards the global model with all explanatory variables
rdaG<- rda(X ~ ., env_ser[,1:3])
anova(rdaG, perm=10000)

### Select the variables
Sel <- ordiR2step(rda0, scope = formula(rdaG), direction="both") 

### Summary table with selected variables    
Sel$anova

############### SELECT MEM ####

### Download geographic info
mpa <- read.table("../../../00-Data/04-Geo/distance_to_mpa_serran_468ind.txt", sep="\t", header=TRUE)
Coor=select(mpa, LON, LAT)
Coorxy=select(mpa, LAT, LON)

### Check the distribution of longitude and lattitude
summary(Coor)

### Get bathymetric dataset
bathydata <- marmap::getNOAA.bathy(lon1= 10,
                                   lon2= -4,
                                   lat1= 34,
                                   lat2= 45,
                                   resolution = 1)

### Create a map to ensure you get the right geographic positions
plot(bathydata, lwd = c(0.3, 1), lty = c(1, 1),
     deep = c(-4500, 0), shallow = c(-50, 0), 
     step = c(500, 0),
     col = c("grey", "black"), drawlabels = c(FALSE, FALSE))
points(mpa$LON, mpa$LAT, pch = 21, col = "black", bg = "grey", cex = 1)

### Calculate the shortest by-water path
#first, define the constraint for calculating a least-cost-path. make transition object.
t <- marmap::trans.mat(bathydata, min.depth = -0.1) #travel cant be shallower than 0.1 metres depth. It can take a while to calculate. This is a "transition object"

### Get km distance matrix
leastDist.km <- marmap::lc.dist(t, Coor, res = "dist") #use "dist" instead of path to get a kilometres distance matrix between all sites.

#Compute in-water MEM
dbMEM_inwater <- adespatial::dbmem(leastDist.km, MEM.autocor = "non-null", store.listw = TRUE)

### Transform into a dataframe
dbMEM.vectors.inwater <- as.data.frame(dbMEM_inwater)

### Transform spatial coordinates into cartesian coordinates 
geo <- SoDA::geoXY(mpa$LAT, mpa$LON, unit=1000)

### Create euclidean distance matrix
euclidean_distances <- dist(geo, method="euclidean") 

#### Perform dbMEM on euclidean distances
dbMEM <- adespatial::dbmem(euclidean_distances, MEM.autocor = "non-null", store.listw = TRUE)

### Transform into a dataframe
dbMEM.vectors.euclidian <- as.data.frame(dbMEM)

### Observe the correlation between dbMEM inwater and dbMEM from euclidian distances
cor.test(dbMEM$MEM1,dbMEM.vectors.inwater$MEM1)

# Keeping default setting for truncation (length of the longest edge of the minimum spanning tree will be used as the threshold) and just positive MEM
dbmem = dbmem(Coor) 

#Look at general output
summary(dbmem)
dbmem <- as.data.frame(dbmem)

###Compute other MEM relatively to spatial weighting scheme
# Gabriel graph
coor_matrix <- as.matrix(Coor)
DistSpatial=gcd.hf(Coor) 

### Distance relative to Gabriel graph
nbnear1 <- dnearneigh(coor_matrix, 0, 0.2)
gabriel <- graph2nb(relativeneigh(coor_matrix, nnmult=10))
list_gabriel <- nb2listw(gabriel,zero.policy = TRUE)
dbmem_gab = mem(list_gabriel) 

### Distance relative to Delaunay triangulation
coor_matrix_jitter <- jitter(coor_matrix)
nbtri <- tri2nb(coor_matrix_jitter)
list_delaunay <- nb2listw(nbtri,zero.policy = TRUE)
dbmem_delaunay = mem(list_delaunay) 

###### Selection for PCNM function 

#OrdiR2step will start working from an empty model without explanatory variables, just the intercept
rda0<-rda(X ~ 1, dbmem)

#OrdiR2step will move towards the global model with all explanatory variables
rdaG<- rda(X ~ ., dbmem)
anova(rdaG, perm=10000)

### Select variables
Sel <- ordiR2step(rda0, scope = formula(rdaG), direction="forward") 
Sel$anova

###### Selection for Gabriel graph
#OrdiR2step will start working from an empty model without explanatory variables, just the intercept
rda0<-rda(X ~ 1, dbmem_gab[,1:300])

#OrdiR2step will move towards the global model with all explanatory variables
rdaG<- rda(X ~ ., dbmem_gab[,1:300])
anova(rdaG, perm=10000)

### Select variables
Sel <- ordiR2step(rda0, scope = formula(rdaG), direction="forward") 
Sel$anova

###### Selection for Delaunay triangulation
#OrdiR2step will start working from an empty model without explanatory variables, just the intercept
rda0<-rda(X ~ 1, dbmem_delaunay[,1:300])

#OrdiR2step will move towards the global model with all explanatory variables
rdaG<- rda(X ~ ., dbmem_delaunay[,1:300])
anova(rdaG, perm=10000)

### Select the variables
Sel <- ordiR2step(rda0, scope = formula(rdaG), direction="forward") 

### Summary table with selected variables    
Sel$anova
rownames(Sel$anova)

dbmem_gab_df <- as.data.frame(dbmem_gab)
#mem_select <- dbmem_gab_df%>% select(MEM22, MEM32,MEM334,MEM127,MEM307, MEM245,MEM157,MEM69,MEM36,MEM398,MEM112,MEM134,MEM211, MEM296) 
mem_select <- dbmem_gab_df%>% select("MEM32","MEM7","MEM196","MEM293","MEM31","MEM263","MEM254","MEM210")   # mullus     

dbmem_delaunay_df <- as.data.frame(dbmem_delaunay)
#mem_select <- dbmem_delaunay_df%>% select(MEM6,MEM246,MEM20,MEM198,MEM4,MEM138,MEM36,MEM41,MEM27,MEM33) # serranus

### Select 
#no mems selected for diplodus
mem_select <- cbind(dbmem$MEM1,dbmem$MEM2, dbmem$MEM5, dbmem$MEM6,dbmem$MEM22)# mullus
mem_select <- cbind(dbmem$MEM4) # serranus
  
### Check the model ith only dbmem
rdamem <- rda(X ,mem_select)
summary(rdamem, scaling=1)  

# Check the RDA summary
RsquareAdj(rdamem)
anova(rdamem, perm=10000)

############### SELECT INSIDE/OUTSIDE ####

### Add MPA info
mpa_nocoteblue <- filter(mpa, MPA!="Cote Bleue")

### Produce dummy variables
inside_dummy <- fastDummies::dummy_cols(mpa_nocoteblue, select_columns = "CATEGORY",remove_first_dummy = TRUE)

#OrdiR2step will start working from an empty model without explanatory variables, just the intercept
rda_category <-rda(X ~ inside_dummy$CATEGORY)

# Check the RDA summary
RsquareAdj(rda_category)
anova(rda_category, perm=10000)

############### SELECT MPA ####

### Create dummy variables
results <- fastDummies::dummy_cols(mpa_nocoteblue, select_columns = "MPA",remove_first_dummy = TRUE)
mpa_dummy <- results[6:ncol(results)] 

#OrdiR2step will start working from an empty model without explanatory variables, just the intercept
rda0<-rda(X ~ 1, mpa_dummy)

#OrdiR2step will move towards the global model with all explanatory variables
rdaG<- rda(X ~ ., mpa_dummy)
anova(rdaG, perm=10000, by="margin")

### Select the variables
Sel <- ordiR2step(rda0, scope = formula(rdaG), direction="forward") 

### Summary table with selected variables    
Sel$anova

### Create a final object with env, mem and mpa variables
var_selected <- as.data.frame(cbind(env_ser$PC2,mpa_dummy$`MPA_Illes Columbretes`)) # serranus

### Select mem variables
### Build the model for diplodus
rdafinal <- rda(X ~ ., var_selected)
anova(rdafinal, perm=10000) 
RsquareAdj(rdafinal)


### Partial RDA
mpa_selected <- var_selected[1:2]
rda_partial <- rda(X, mpa_selected$V2, var_selected$V1)
anova(rda_partial, perm=10000) 
RsquareAdj(rda_partial)

### Build the model
rdaS<- rdafinal

### Check the RDA summary
RsquareAdj(rdaS)

### Test the significance
anova(rdaS, perm=999)
anova(rdaS, perm=999, by="margin")
anova(rdaS, perm=999, by="axis")

################# VISUALISE DB-RDA results ##################

### Getting the scores for plotting.
scrs <- scores(rdaS, display = c("species", "sites", "bp"), choices = 1:4, scaling = 2)
scrs

### Collect information on the pcao axes
species_centroids <- data.frame(scrs$species)
species_centroids
species_centroids$PC_names <- rownames(species_centroids)

### Add the PC components names.
species_centroids$PC_names <- rownames(species_centroids) 

### Add information of arrows
continuous_arrows <- data.frame(scrs$biplot)
continuous_arrows
rownames(continuous_arrows) <- c("MEM4","Columbretes", "Cabo de Palos") # diplodus
rownames(continuous_arrows) <- c("PC1 env","MEM1","MEM2","MEM5","MEM16","MEM20","MEM24","Cerbere-Banyuls","Cap de Creus")
rownames(continuous_arrows) <- c("PC2 env","Illes de Columbretes")
rownames(continuous_arrows) <- c("Cerbere-Banyuls") # palinurus

### Add names of variables
continuous_arrows$number <- rownames(continuous_arrows)
continuous_arrows
baseplot<-plot(rdaS, scaling = 2)
mult <- attributes(baseplot$biplot)$arrow.mul

### Add mpa names
sitenames<- mpa$MPA
sites_centroids <- data.frame(scrs$sites)
sites_centroids$SITE <- sitenames
head(sites_centroids)

### Add protection status
sitestatus <- mpa$CATEGORY
sites_centroids$STATUS <- sitestatus
  
### Specify the color
br_pal <- brewer.pal(8,"RdYlBu") 

### Have the right order of labels
sites_centroids$SITE <- factor(sites_centroids$SITE, levels = c("Cabo de Gata Níjar", "Cabo de Palos", "Illa De Tabarca","Illes Columbretes","Llevant De Mallorca","Norte de Menorca","Cap de Creus","Cerbere-Banyuls"))

### Make a ggplot
RDA_plot <- ggplot(data = sites_centroids, aes(x = RDA1, y = RDA2))+
  geom_point(data = sites_centroids, size = 2, aes(fill = SITE, shape = STATUS))+
  scale_shape_manual(values=c(21,24))+
  scale_fill_manual(values = br_pal, name="Marine reserves")+
  geom_text(data = continuous_arrows,
            aes(x= (mult + mult/10) * RDA1, y = (mult + mult/10) * RDA2,
                label = number), size = 4, hjust = 0.5)+
  geom_segment(data = continuous_arrows,
               aes(x = 0, xend = mult * RDA1, y = 0, yend = mult * RDA2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, lty = "dotted") + geom_vline(xintercept = 0, lty = "dotted") +
  labs(x = paste("RDA1 (", round(summary(rdaS)$cont$importance[2,1]*100, 2), "%)", sep = ""), y = paste("RDA2 (", round(summary(rdaS)$cont$importance[2,2]*100, 2), "%)", sep = ""))+
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "black")))
RDA_plot

### For palinurus
RDA_plot <- ggplot(data = sites_centroids, aes(x = RDA1, y = PC1))+
  geom_point(data = sites_centroids, size = 2, aes(fill = SITE, shape = STATUS))+
  scale_shape_manual(values=c(21,24))+
  scale_fill_manual(values = br_pal, name="Marine reserves")+
  geom_text(data = continuous_arrows,
            aes(x= (mult + mult/10) * RDA1, y = (mult + mult/10) * PC1,
                label = number), size = 4, hjust = 0.5)+
  geom_segment(data = continuous_arrows,
               aes(x = 0, xend = mult * RDA1, y = 0, yend = mult * PC1),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, lty = "dotted") + geom_vline(xintercept = 0, lty = "dotted") +
  labs(x = paste("RDA1 (", round(summary(rdaS)$cont$importance[2,1]*100, 2), "%)", sep = ""), y = paste("PC1 (", round(summary(rdaS)$cont$importance[2,2]*100, 2), "%)", sep = ""))+
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "black")))
RDA_plot


ggsave("RDA_serranus_april.pdf", height=7, width=7)


