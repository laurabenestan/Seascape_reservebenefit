# seascape_reservebenefit

The objective of this project is to dissect the link between genomics and environmental gradients of differentiation in a context of conservation of marine reserves.
Using a state-of-the-art approach based on seascape genomics, we aim to:

01. Estimating genetic diversity within and outside a network of marine reserve, delineating the influence of spatial and environmental gradient on the genetic diversity observed
02. Comparing historical and contemporary population structure across species occupying a similar seascape
03. Characterising the seascape differentiation experienced by the species
04. Determining the influence of the marine reserve position, seascape features and spatial distribution on the population structure 
05. Bringing useful recommendations for marine reserve design


# 01-Genetic diversity

## Genetic diversity estimation

Genetic diversity is estimated with [Plink1.9](https://www.cog-genomics.org/plink/1.9/basic_stats), though the use of the --het function that adjust the calculation with the sample size. 
We transformed the output obtained from PLINK to Observed heterozygosity values for each SNP datasets (neutral and adaptive).

## Testing the difference between genetic diversity of individuals sampled inside and outside one of the 8 marine reserves sampled for the project

We want to know the differences of the genetic diversity (hereafter observed heterozygosity) for each species regarding the category "Inside/outside".
First, we check the distribution of observed heterozygosity and we test the normality of the distribution using Shapiro test.

The observed heterozygosity does not follow normality.
As is highly skew to positive values. It seems like a beta distribution but we will check it later.

First, we considered a Wilcoxon test to compare observed heterozygosity between samples inside/outside a marine reserve.

#### Neutral genetic diversity inside/outside for the three species

| Species | statistic | p.value |
|--------|--------------------------------------------------|-------------|
| Diplodus sargus | 8817 | 6.01e-1 | 
| Mullus surmuletus | 11424 | 7.42e-1 |
| Serranus cabrilla | 20314 | 1.71e-6 |

#### Adaptive genetic diversity inside/outside for the three species

| Species | statistic | p.value |
|--------|--------------------------------------------------|-------------|
| Diplodus sargus | 8984 | 7.90e-1 | 
| Mullus surmuletus | 11175 | 9.99e-1 |
| Serranus cabrilla | 22082 | 3.50e-4 |

Wilcoxon Signed Rank test revealed that both neutral and adaptive observed heterozygosity were significantly different for Serranus cabrilla but not for Diplodus sargus and Mullus surmuletus.


<img align="center" height="240" src="01-genetic_diversity/FigS2.pdf"></img>

According to our expectation, genetic diversity would be higher in non impacted areas such as a marine reserve.
Yet, here we observed the opposite trend. 
So, we want to explore more on the environmental determinant of genetic diversity and identify which determinant may influence genetic diversity in the three species.

##  Principal components on 24 seascape features that includes salinity, chlorophyll and temperature

We perform a Principal Component Analysis (PCA) on 24 seascape features to avoid issues regarding to collinearity.
We keep 3 axes of the PCA based on the Kaiser criterion: keeping component that accounts for at least 5% of the total variance.

<img align="center" height="240" src="01-genetic_diversity/pca_axes.png"></img>

## Exploring the distribution of genetic diversity

We use the function `descdist`available in the package `fitdistrplus` to compute descriptive parameters of our empirical distribution.
The kurtosis and squared skewness of our sample is plotted as a blue point named "Observation".

<img align="center" height="240" src="01-genetic_diversity/observation.png"></img>

It seems that possible distributions include beta distribution.
We then compare the The Akaike information criterion (AIC) for the two models.

| Model | AIC|
|--------|-------------|
| Beta distribution | -1909.775 | 
| normal distribution | -1929.488 |

Here, AIC provides a means for model selection and the preferred model is the one built using beta distribution.