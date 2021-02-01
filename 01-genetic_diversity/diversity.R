### Download libraries
pkgs <- c('tidyr','ggplot2','cowplot','broom','dplyr','ggpubr','factoextra','here')
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))


### Download genetic diversity, distances and environmental information
distances_mpa <- read.table("distances_to_mpa_all.txt",header=TRUE,sep="\t")
het <- read.table("genetic_diversity_all.txt",header=TRUE,sep="\t", dec=".")
env_habitats <- read.table("env_data_landscape.txt", header=TRUE)

################### TEST FOR INSIDE/OUTSIDE RESERVE ####

### Merge morpho and distances to mpa information with category information
morpho_distances_all <- merge(x=het,y=distances_mpa,by.x=c("INDV"),by.y=c("IND"))
morpho_distances_all$SPECIES <- substr(morpho_distances_all$INDV, 1,3)
morpho_distances_all$SPECIES <- ifelse(morpho_distances_all$SPECIES =="mul"| morpho_distances_all$SPECIES =="Mul"| morpho_distances_all$SPECIES =="mul","Mullus surmuletus", ifelse(morpho_distances_all$SPECIES=="dip","Diplodus sargus", ifelse(morpho_distances_all$SPECIES=="ser","Serranus cabrilla","Palinurus elephas")))

### Create a column representing the buffer at 5km (if higher, unbalanced number of samples)
morpho_distances_all$BUFFER5 <- ifelse(morpho_distances_all$distance < 5000, "Inside","Outside")

### Check how is distributed the dataset
aggregate(NEUTRAL_HET ~ BUFFER5+SPECIES+MPA.NAME_EN, data = morpho_distances_all, FUN = length)

### Subset by species
diplodus_morpho_fish_het <- subset(morpho_distances_all, subset=morpho_distances_all$SPECIES=="Diplodus sargus")
mullus_morpho_fish_het <- subset(morpho_distances_all, subset=morpho_distances_all$SPECIES=="Mullus surmuletus")
serranus_morpho_fish_het <- subset(morpho_distances_all, subset=morpho_distances_all$SPECIES=="Serranus cabrilla")

### Check the distribution
ggplot(morpho_distances_all, aes(x=NEUTRAL_HET))+
  geom_density()+
  facet_grid(~SPECIES)+
  theme_classic()

# the distribution is highly skew to positive values. It seems like a beta distribution

### Normality test
morpho_distances_all %>% 
     group_by(SPECIES)  %>% 
    do(tidy(shapiro.test(.$NEUTRAL_HET))) %>% 
    ungroup() %>% 
    select(-method)

### Check linear regression 
morpho_distances_all %>% group_by(SPECIES) %>%
  do(fitGen = tidy(lm(NEUTRAL_HET ~ distance, data = .))) %>% unnest(fitGen)


### Non parametric tests
morpho_distances_all %>% 
  group_by(SPECIES)  %>% 
  do(tidy(wilcox.test(ADAPTIVE_HET~ BUFFER5, data = .))) %>% 
  ungroup()

### Remove south effect
`%notin%` <- Negate(`%in%`)
serranus_morpho_fish_het_north <-subset(serranus_morpho_fish_het, MPA.NAME_EN %notin% c('Cabo de Gata Nijar', 'Cabo de Palos'))
wilcox.test(NEUTRAL_HET ~ BUFFER5, data = serranus_morpho_fish_het_north)

### Make a graph
g1 <- ggline(serranus_morpho_fish_het, x = "BUFFER5", y = "NEUTRAL_HET", 
       add = c("mean_se"), 
       order = c("Inside", "Outside"),
       ylab = "Neutral heterozygosity", xlab = "Position to marine reserve")
g2 <- ggline(serranus_morpho_fish_het, x = "BUFFER5", y = "ADAPTIVE_HET", 
             add = c("mean_se"), 
             order = c("Inside", "Outside"),
             ylab = "Adaptive heterozygosity", xlab = "Position to marine reserve")
g2

### Save the graph
pdf("Fig1.pdf", height=5, width=10)
plot_grid(g1,g2)
dev.off()

################### PCA FOR DATA ####

### Merge morpho_distances_all information with env information
env <- read.table("env_data_landscape.txt",header=TRUE)
EnvMatrix <- env[2:24]

### Scale before PCA
sEnvMatrix = scale(EnvMatrix, center = T, scale = T)

# We can use the prcomp() function to perform the principal component analysis.
ePCA = prcomp(sEnvMatrix)

# The principal components can be found in the $x matrix:
pca_axis <- as.data.frame(ePCA$x)
pca_axis$IND <- env$labels

### PCA on environmental variables
morpho_distances_env <- merge(x=morpho_distances_all,y=pca_axis,by.x=c("INDV"),by.y=c("IND"))

### Check the number of axis to keep
fviz_eig(ePCA)


################### DISTRIBUTION OF HET ####

### Download libraries
library("fitdistrplus")
library("logspline")
library("betareg")
library("grid")

### Find the distribution of genetic diversity
descdist(as.vector(diplodus_morpho_fish_het$NEUTRAL_HET), discrete = FALSE)

#The kurtosis and squared skewness of your sample is plotted as a blue point named "Observation". 
#It seems that possible distributions include beta distribution.
fit.beta <- fitdist(as.vector(diplodus_morpho_fish_het$NEUTRAL_HET), "beta")
fit.norm <- fitdist(as.vector(diplodus_morpho_fish_het$NEUTRAL_HET), "norm")

### Compare the residuals plots for the two models
plot(fit.norm)
plot(fit.beta)

### Compare the Akaike Criterion of the two models
fit.beta$aic
fit.norm$aic

################### BUILD MODELS BETWEEN GENETIC DIVERSITY AND ENVIRONMENTAL DATA ####

### test several models
  morpho_distances_env %>% 
  group_by(SPECIES)  %>% 
  do(tidy(betareg(NEUTRAL_HET~ PC2, data = .))) %>% 
  ungroup()
# the only significant variable is PC2

### Check the correlation between these PC2 and distance
cor.test(morpho_distances_env$distance, morpho_distances_env$PC2)

### Build a model 
grob <- grid::grobTree(textGrob("Pseudo R-squared = 0.0331", x=0.4,  y=0.95, hjust=0,
                          gp=gpar(col="black", fontsize=10, fontface="bold")))

diplodus <- ggplot(data = subset(morpho_distances_env, subset=morpho_distances_env$SPECIES=="Diplodus sargus"),aes(x=PC2,y=NEUTRAL_HET)) + 
  geom_point() +
  geom_smooth(method="loess",color="magenta")+
  theme_classic()+
  ylab("Genetic diversity")+
  xlab("")+
  annotation_custom(grob)
diplodus

### Remove the outliers individuals from Serranus cabrilla
no_outliers <- subset(morpho_distances_env, subset=morpho_distances_env$PC2 > -4)

### beta model
mod <- betareg(NEUTRAL_HET ~ PC2, data = subset(no_outliers, subset=no_outliers$SPECIES=="Serranus cabrilla"))
summary(mod)

# mullus no link with env
grob <- grobTree(textGrob("Pseudo R-squared = 0.03343", x=0.4,  y=0.95, hjust=0,
                          gp=gpar(col="black", fontsize=10, fontface="bold")))
serranus <- ggplot(data = subset(no_outliers, subset=no_outliers$SPECIES=="Serranus cabrilla"),aes(x=PC2,y=NEUTRAL_HET)) + 
  geom_point() +
  geom_smooth(method="loess",color="red")+
  theme_classic()+
  ylab("Genetic diversity")+
  xlab("Principal Component 2 of the environmental variation")+
  annotation_custom(grob)
serranus

pdf("Fig2.pdf",width = 4, height=5)
plot_grid(diplodus, serranus, names="auto", ncol=1, nrow=2)
dev.off()
save.image("Fig2.RData")

