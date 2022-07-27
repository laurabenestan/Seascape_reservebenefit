---
title: "Guidelines for the dbRDA"
output:
  html_document:
    toc: true
    theme: united
---

# Is space, environment or protection that explain the most the genetic variation observed among species ?
Individual-based distance)based Redundancy Analysis 

# 1. Y = GENOMIC DATA 

## Import your vcf file using the function `read.vcfR`
```{r include = FALSE}
vcf <- read.vcfR("yourvcffile")
```

## Transform vcf to genind file
```{r include = FALSE}
genind <- vcfR2genind(vcf)
```

## Calculate Euclidean distances
```{r include = FALSE}
distgenEUCL <- dist(genind, method = 
                      "euclidean", diag = FALSE, upper = FALSE, p = 2)
                      ```
                      
## Ordination on Y
```{r include = FALSE}
Pcoa=pcoa(distgenEUCL)
Pcoa
 ```
 
## Check the percent of variation explained by each axis
```{r}
Pcoa$values$Cumul_eig*100
```

## One of the rules is based on the cumulative percentage explained, i.e. retain the components which capture, say, 70% or 90% of the variation. Keep the number of Pcoa principal components that explained up to 90% of the genomic variation, which will be the response variable in the db-RDA.
```{r}
X=as.data.frame(Pcoa$vectors[,1:390])
```

## Look at genotypes distribution in relation to the first 2 Pcoa axes
```{r}
plot(X[,1], X[,2])
```

# 2.  X = PCA ENV

## Download the environmental matrix
```{r}
env = read.table('../00-data/all_env_data_4species.txt', row.names=1, header=T)
EnvMatrix <- env[2:ncol(env)]
```

## Scale before PCA
```{r}
sEnvMatrix = scale(EnvMatrix, center = T, scale = T)
```

## We can use the prcomp() function to perform the principal component analysis.
```{r}
ePCA = prcomp(sEnvMatrix)
```

## The principal components can be found in the matrix
```{r}
pca_axis <- as.data.frame(ePCA$x)
fviz_eig(ePCA)
```

## Subset the individuals to keep in the Env file containing the four species information
```{r}
env_ser <- pca_axis[588:1055,]
```

## OrdiR2step will start working from an empty model without explanatory variables, just the intercept
```{r}
rda0<-rda(X ~ 1, env_ser[,1:3])
```

## OrdiR2step will move towards the global model with all explanatory variables
```{r}
rdaG<- rda(X ~ ., env_ser[,1:3])
anova(rdaG, perm=10000)
```

### Select the variables
```{r}
Sel <- ordiR2step(rda0, scope = formula(rdaG), direction="both") 
```

### Summary table with selected variables   
```{r}
Sel$anova
```

# 3. X = HABITATS

## Import habitats variables
```{r}
habitat <- read.table("../../../00-Data/05-Habitats/1299ind_habitats.txt", header=TRUE, sep="\t",dec=".")
habitat_select <- habitat[588:1055,]
```

## Produce dummy variables
```{r}
habitats_dummy <- fastDummies::dummy_cols(habitat_select, select_columns = "Habitats",remove_first_dummy = TRUE)
```

# 4. SPACE : MEM

## Download geographic info
```{r}
mpa <- read.table("../../../00-Data/04-Geo/distance_to_mpa_serran_468ind.txt", sep="\t", header=TRUE)
Coor=select(mpa, LON, LAT)
Coorxy=select(mpa, LAT, LON)
```

## Check the distribution of longitude and lattitude
```{r}
summary(Coor)
```

## Get bathymetric dataset
```{r}
bathydata <- marmap::getNOAA.bathy(lon1= 10,
                                   lon2= -4,
                                   lat1= 34,
                                   lat2= 45,
                                   resolution = 1)
```

## Create a map to ensure you get the right geographic positions
```{r}
plot(bathydata, lwd = c(0.3, 1), lty = c(1, 1),
     deep = c(-4500, 0), shallow = c(-50, 0), 
     step = c(500, 0),
     col = c("grey", "black"), drawlabels = c(FALSE, FALSE))
points(mpa$LON, mpa$LAT, pch = 21, col = "black", bg = "grey", cex = 1)
```

## Calculate the shortest by-water path. First, define the constraint for calculating a least-cost-path. make transition object.
```{r}
t <- marmap::trans.mat(bathydata, min.depth = -0.1) #travel cant be shallower than 0.1 metres depth. It can take a while to calculate. This is a "transition object"
```

## Get km distance matrix
```{r}
leastDist.km <- marmap::lc.dist(t, Coor, res = "dist") #use "dist" instead of path to get a kilometres distance matrix between all sites.
```

## Compute in-water MEM
```{r}
dbMEM_inwater <- adespatial::dbmem(leastDist.km, MEM.autocor = "non-null", store.listw = TRUE)
```

## Transform into a dataframe
```{r}
dbMEM.vectors.inwater <- as.data.frame(dbMEM_inwater)
```

## Transform spatial coordinates into cartesian coordinates
```{r}
geo <- SoDA::geoXY(mpa$LAT, mpa$LON, unit=1000)
```

## Create euclidean distance matrix
```{r}
euclidean_distances <- dist(geo, method="euclidean") 
```

## Perform dbMEM on euclidean distances
```{r}
dbMEM <- adespatial::dbmem(euclidean_distances, MEM.autocor = "non-null", store.listw = TRUE)
```

## Transform into a dataframe
```{r}
dbMEM.vectors.euclidian <- as.data.frame(dbMEM)
```

## Observe the correlation between dbMEM inwater and dbMEM from euclidian distances
```{r}
cor.test(dbmem$MEM1,dbMEM.vectors.inwater$MEM1)
```

## Keeping default setting for truncation (length of the longest edge of the minimum spanning tree will be used as the threshold) and just positive MEM
```{r}
dbmem = dbmem(Coor) 
```

## Look at general output
```{r}
summary(dbmem)
dbmem <- as.data.frame(dbmem)
```

## Select the MEMs variables using `ordistep` function
```{r}
rda0<-rda(X ~ 1, dbmem)
```

## OrdiR2step will move towards the global model with all explanatory variables
```{r}
rdaG<- rda(X ~ ., dbmem)
anova(rdaG, perm=10000)
```

## Select variables
```{r}
Sel <- ordiR2step(rda0, scope = formula(rdaG), direction="forward") 
Sel$anova
```

# 5. PROTECTION LEVEL : INSIDE/OUTSIDE ####

### Remove Cote bleue samples
```{r}
mpa_nocoteblue <- filter(mpa, MPA!="Cote Bleue")
```

## Produce dummy variables for inside/outside level
```{r}
inside_dummy <- fastDummies::dummy_cols(mpa_nocoteblue, select_columns = "CATEGORY",remove_first_dummy = TRUE)
```

## Select significant variables
```{r}
rda_category <-rda(X ~ inside_dummy$CATEGORY)
```

## Check the RDA summary
```{r}
RsquareAdj(rda_category)
anova(rda_category, perm=10000)
```

### Create dummy variables for marine reserves locations
```{r}
results <- fastDummies::dummy_cols(mpa_nocoteblue, select_columns = "MPA",remove_first_dummy = TRUE)
mpa_dummy <- results[6:ncol(results)] 
```

## Select significant variables
```{r}
rda0<-rda(X ~ 1, mpa_dummy)
rdaG<- rda(X ~ ., mpa_dummy)
anova(rdaG, perm=10000, by="margin")
Sel <- ordiR2step(rda0, scope = formula(rdaG), direction="forward") 
```

## Summary table with selected variables    
```{r}
Sel$anova
```

## Create a final object with env, mem and mpa variables
```{r}
var_selected <- as.data.frame(cbind(env_ser$PC2,dbmem$MEM4,dbmem$MEM9,dbmem$MEM17,mpa_dummy$`MPA_Illes Columbretes`))
```

# 6. FINAL MODEL AND RESULTS

## Build the global model
```{r}
rdafinal <- rda(X ~ ., var_selected)
anova(rdafinal, perm=10000) 
RsquareAdj(rdafinal)
```

## Getting the scores for plotting.
```{r}
scrs <- scores(rdafinal, display = c("species", "sites", "bp"), choices = 1:4, scaling = 2)
scrs
```

## Collect information on the pcao axes
```{r}
species_centroids <- data.frame(scrs$species)
species_centroids
species_centroids$PC_names <- rownames(species_centroids)
```

### Add the PC components names.
```{r}
species_centroids$PC_names <- rownames(species_centroids) 
```

### Add information of arrows
```{r}
rownames(continuous_arrows) <- c("PC2 env","MEM4","MEM9","MEM17","Illes de Columbretes")
```

### Add names of variables
```{r}
continuous_arrows$number <- rownames(continuous_arrows)
continuous_arrows
baseplot<-plot(rdafinal, scaling = 2)
mult <- attributes(baseplot$biplot)$arrow.mul
```

### Add mpa names
```{r}
sitenames<- mpa$MPA
sites_centroids <- data.frame(scrs$sites)
sites_centroids$SITE <- sitenames
head(sites_centroids)
```

### Add protection status
```{r}
sitestatus <- mpa$CATEGORY
sites_centroids$STATUS <- sitestatus
```

### Specify the color
```{r}
br_pal <- brewer.pal(8,"RdYlBu") 
```

### Have the right order of labels
```{r}
sites_centroids$SITE <- factor(sites_centroids$SITE, levels = c("Cabo de Gata Níjar", "Cabo de Palos", "Illa De Tabarca","Illes Columbretes","Llevant De Mallorca","Norte de Menorca","Cap de Creus","Cerbere-Banyuls"))
```

### Make a ggplot
```{r}
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
  labs(x = paste("RDA1 (", round(summary(rdafinal)$cont$importance[2,1]*100, 2), "%)", sep = ""), y = paste("RDA2 (", round(summary(rdafinal)$cont$importance[2,2]*100, 2), "%)", sep = ""))+
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "black")))
RDA_plot
```

