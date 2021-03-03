### Download librairies
library(ggplot2)
library(cowplot)
library(RColorBrewer)

### Download the environmental matrix
env = read.table('env_data_landscape.txt', row.names=1, header=T)
EnvMatrix <- env[2:24]

### Scale before PCA
sEnvMatrix = scale(EnvMatrix, center = T, scale = T)

### The columns correspond to different environmental variables, the rows to the different sampling locations.
heatmap(sEnvMatrix) 

### We can use the prcomp() function to perform the principal component analysis.
ePCA = prcomp(sEnvMatrix)
ePCA

### Check the percent of variation
S=summary(ePCA)
S

### The principal components can be found in the $x matrix:
pca_axis <- as.data.frame(ePCA$x)
pca_axis$IND <- row.names(pca_axis)

### Note that the principal components are non-correlated  (orthogonal) between each other, and therefore describe non overlapping patterns of environmental variations. 
cor(ePCA$x)

### The question now is, what does this main axis of variation mean? 
### We can have the answer by looking at another matrix stored in the ePCA object: the rotation matrix 
ePCA$rotation

### This matrix informs us on which environmental variables participate to the composition of each principal components. The contribution of each variable is called "loading". 
### We can visualize this information for PC1:
pc1 <- as.data.frame(ePCA$rotation[,'PC1'])
pc1$variables <- row.names(pc1) 
colnames(pc1) <- c("contribution", "variables")

### We can visualize this information for PC2:
pc2 <- as.data.frame(ePCA$rotation[,'PC2'])
pc2$variables <- row.names(pc2) 
colnames(pc2) <- c("contribution", "variables")

### Make a ggplot for PC1
g1 <- ggplot(pc1, aes(x = variables, y = contribution)) +
  geom_bar(
    stat = "identity", position = position_stack(),
    color = "black", fill = "grey"
  ) +
  theme_bw()+
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),  
        axis.title.x = element_text(color = "black", size = 14),
        axis.title.y = element_text(color = "black", size = 14))+
  xlab("Environmental variables")+
  ylab("Variable loading")+
  coord_flip()
g1

### Make a ggplot for PC2
g2 <- ggplot(pc2, aes(x = variables, y = contribution)) +
  geom_bar(
    stat = "identity", position = position_stack(),
    color = "black", fill = "grey"
  ) +
  xlab("Environmental variables")+
  ylab("Variable loading")+
  coord_flip()+
  theme_bw()+
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),  
        axis.title.x = element_text(color = "black", size = 14),
        axis.title.y = element_text(color = "black", size = 14))
g2

### Download grouping information
group <- read.table("../00-data/distances_to_mpa_no_coteblue.txt",header=T,sep="\t", dec=".")

### Check the PCA results
pca_env_mpa <- merge(pca_axis, group, by="IND")

# Add MPA information
pca_env_mpa$POP <- substr(pca_env_mpa$IND,1,3)
unique(pca_env_mpa$POP)
pca_env_mpa$SPECIES <- ifelse(pca_env_mpa$POP =="Mul"| pca_env_mpa$POP =="mul","Mullus surmuletus", ifelse(pca_env_mpa$POP=="dip","Diplodus sargus", "Serranus cabrilla"))

### Remove Cote Bleue sample
pca_env_mpa <- subset(pca_env_mpa, subset=pca_env_mpa$MPA.NAME_EN!="Cote Bleue")

### Have the right order of labels relatively to latitude gradient
pca_env_mpa$MPA.NAME_EN <- factor(pca_env_mpa$MPA.NAME_EN, levels = c("Cabo de Gata Níjar", "Cabo de Palos", "Illa De Tabarca","Illes Columbretes","Llevant De Mallorca","Norte de Menorca","Cap de Creus","Cerbere-Banyuls"))

### Make a ggplot
g3 <- ggplot(pca_env_mpa, aes(x=PC1,y=PC2, shape=SPECIES, fill=MPA.NAME_EN))+
  geom_point(size=3,aes(fill=MPA.NAME_EN))+
  scale_shape_manual(values=c(21,23,25),name="Species")+
  scale_fill_brewer(palette="RdYlBu", name="Marine reserves")+
  theme_bw()+
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),  
        axis.title.x = element_text(color = "black", size = 14),
        axis.title.y = element_text(color = "black", size = 14))+
  xlab("Axis 1 (58.2%)")+
  ylab("Axis 2 (15.6%)")+
  guides(fill=guide_legend(override.aes=list(shape=21)))
g3

### Make a combined graph
pdf("Figure2.pdf", width=10, height=15)
first_row <- plot_grid(g3, labels = c('A'))
second_row <- plot_grid(g1,g2, labels = c('B', 'C'), nrow = 1)
plot_grid(first_row,second_row, ncol=1,nrow=2)
dev.off()

