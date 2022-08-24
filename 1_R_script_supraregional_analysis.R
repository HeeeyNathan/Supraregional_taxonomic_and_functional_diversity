###################################################################################################
# Supplementary Data 1                                                                            #
#                                                                                                 #
# R script to create results presented in:                                                        # 
#                                                                                                 #
# Seasonal and spatial variation of stream macroinvertebrate taxonomic and functional diversity   #
# across three boreal regions                                                                     #
#                                                                                                 #
# Nathan J. Baker, Ellen A. R. Welti, Francesca Pilotto, Jonas Jourdan, Burkhard Beudert,         #
# Kaisa-Leena Huttunen, Timo Muotka, Riku Paavola, Emma Göthe, Peter Haase                        #
#                                                                                                 #
# Code written by Nathan J. Baker and Francesca Pilotto                                           #
# Code used for functional analyses adapted from Cayetano Gutiérrez-Cánovas & Gabone Iturralde    #
# (https://github.com/tanogc/overarching_functional_space)                                        #
#                                                                                                 #
# Email for queries: Nathan93Baker@gmail.com                                                      #
# Code written using R version 4.2.1 "Funny-Looking Kid"                                          #
# Code written in Rstudio version 2022.07.1+554 "Spotted Wakerobin"                               #
#                                                                                                 #
###################################################################################################

##### CHANGE YOUR WORKING DIRECTORY --------------------
setwd() #set your own WD

##### LOAD LIBRARIES --------------------
# Load the required packages
library(ade4)
library(adespatial) # beta diversity
library(FD) # functional diversity 
library(geodist) # geographic distance matrix
library(ggplot2) # plotting
library(ks) # kernal density estimation
library(pacman) # package manager
library(plyr)
library(RColorBrewer) #colours
library(vegan) # community ecology, taxonomic diversity

##### LOAD 3RD-PARTY ADD-ONS & ADDITIONAL FUNCTIONS --------------------
# Load add-ons from various sources
source("triplot.rda.R") # triplots as per Boccard et al. (2018)
source("HighstatLibV10.R") # pair plots& FUNCTIO

# Additional functions
source("0_FD_functions.R")
source("0_quality_funct_space_fromdist.R")
source("0_vegetarianRpackage_functions.R")

##### LOAD THE DATA --------------------
#### TAXANOMIC DATA
### Full site-by-species matrix
comm <- read.csv(file = "full_site_taxa_mat.csv", h = TRUE, sep = ",", row.names = c(1) ,stringsAsFactors = FALSE)

#### FUNCTIONAL DATA
### Species-by-trait matrix
comm_traits_raw <- read.table("func_taxa_trait_mat.csv", h = T, sep = ",", row.names = 1, stringsAsFactors = FALSE)
comm_traits_raw <- comm_traits_raw[-13,] # remove 'chelifsp' because it has no trait information
comm_traits <- subset(comm_traits_raw, select = size1:fhab8)
rownames(comm_traits) <- comm_traits_raw$Code
comm_traits <- comm_traits[, -38] #remove row 38 (resp5) because rowSums = 0 - trait is redundant

## Changes trait data to proportions (between 0 & 1)
# Trait category blocks
comm_traits_fuzzy <- prep.fuzzy.var(comm_traits, c(7, 2, 3, 4, 8, 4, 5, 4, 8, 9, 8)) # These numbers show how the 11 traits are grouped
comm_traits <- comm_traits_fuzzy
rowSums(comm_traits)

### functional site-by-species matrix
func_comm <- read.table("func_site_taxa_mat.csv", header = TRUE, sep = ",", row.names = 1) # Loading site x taxa matrix
func_comm <- func_comm[,-13] # remove 'chelifsp' because it has no trait information
reg.mat <- read.table("func_reg_taxa_mat.csv", header = TRUE, sep = ",", row.names = 1)
reg.mat <- reg.mat[,-13] # remove 'chelifsp' because it has no trait information

#### ENVIRONMENTAL DATA
env <- read.csv(file = "site_env_mat.csv", header = TRUE, sep = ",", row.names = c(1), na.strings = "NA")

#### SPATIAL DATA
geo <- read.csv(file = "site_geo_mat.csv", header = TRUE, sep = ",", row.names = c(1), stringsAsFactors = FALSE)
dist <- geodist(geo, measure = "geodesic") # A matrix representing physical distances between sites
dist <- as.dist(dist) # distances are in metres (m)

#### SITE GROUPING FACTORS
group <- read.csv("site_grouping_mat.csv", header = TRUE, sep = ",", row.names = c(1), na.strings = "NA")

#### HOW DO THE FUNCITONAL TAXA-COMPLEXES IMPACT REGIONAL SPECIES RICHNESS --------------------
### Determine the number of species lost due to functional taxonomic adjustments
## All regions
d(comm, lev = "gamma", q = 0)
d(func_comm, lev = "gamma", q = 0)
## Germany
d(comm[1:16,], lev = "gamma", q = 0)
d(func_comm[1:16,], lev = "gamma", q = 0)
## Sweden
d(comm[17:120,], lev = "gamma", q = 0)
d(func_comm[17:120,], lev = "gamma", q = 0)
## Finland
d(comm[121:140,], lev = "gamma", q = 0)
d(func_comm[121:140,], lev = "gamma", q = 0)



#### TAXONOMIC ALPHA DIVERSITY --------------------
### Compute alpha diversity indices of the communities
A <- rowSums(comm)               # Overall abundance
N0 <- rowSums(comm > 0)          # Species richness
H <- diversity(comm)             # Shannon entropy (base e)
N1 <- exp(H)                     # Shannon diversity (base e) i.e., the number of abundant species
E10 <- N1 / N0                   # Shannon evenness (Hill's ratio)

div <- data.frame(A, N0, N1, E10)

#### TAXONOMIC BETA DIVERSITY --------------------
# Beta div = The variation in species composition among sites within a geographic area of interest

### Partitioning Beta Diversity into Replacement, Richness Difference and Nestedness 
# Species replacement = Turnover = The change in community composition between two gradients
# Richness difference = species gain and/or loss
# Nestedness = the species at poorer sites are a strict subset of the species at richer sites

## Jaccard-based Podani indices (abundance data)
comm_pod_pd <- beta.div.comp(comm, coef = "S", quant = TRUE)
# What is in the output object?
summary(comm_pod_pd)
# Display summary statistics:
comm_pod_pd$part
comm_pod_pd$Note
# boxplot
boxplot(comm_pod_pd$D)

# Germany
ger_pod_pd <- beta.div.comp(comm[1:16,], coef = "S", quant = TRUE)
boxplot(ger_pod_pd$D)
ger_pod_pd_spr <- beta.div.comp(comm[1:8,], coef = "S", quant = TRUE)
boxplot(ger_pod_pd_spr$D)
ger_pod_pd_aut <- beta.div.comp(comm[9:16,], coef = "S", quant = TRUE)
boxplot(ger_pod_pd_aut$D)

# Sweden
swe_pod_pd <- beta.div.comp(comm[17:120,], coef = "S", quant = TRUE)
boxplot(swe_pod_pd$D)
swe_pod_pd_spr <- beta.div.comp(comm[17:68,], coef = "S", quant = TRUE)
boxplot(swe_pod_pd_spr$D)
swe_pod_pd_aut <- beta.div.comp(comm[69:120,], coef = "S", quant = TRUE)
boxplot(swe_pod_pd_aut$D)

# Finland
fin_pod_pd <- beta.div.comp(comm[121:140,], coef = "S", quant = TRUE)
boxplot(fin_pod_pd$D)
fin_pod_pd_spr <- beta.div.comp(comm[121:130,], coef = "S", quant = TRUE)
boxplot(fin_pod_pd_spr$D)
fin_pod_pd_aut <- beta.div.comp(comm[131:140,], coef = "S", quant = TRUE)
boxplot(fin_pod_pd_aut$D)

# The sum of the second and third values (rounded) = 0.14 + 0.26 = 0.40, 
# i.e., that total replacement diversity and total richness diversity sum to BDTotal. 

#### TAXONOMIC TURNOVER --------------------
### Podani family, percentage difference (Sorensen dissimilarity on quantitative data)
## Germany
FB_beta_div <- beta.div.comp(comm[c(1, 9), ], coef = "S", quant = TRUE)
GO2_beta_div <- beta.div.comp(comm[c(2, 10), ], coef = "S", quant = TRUE)
GO4_beta_div <- beta.div.comp(comm[c(3, 11), ], coef = "S", quant = TRUE)
KR1_beta_div <- beta.div.comp(comm[c(4, 12), ], coef = "S", quant = TRUE)
KR4_beta_div <- beta.div.comp(comm[c(5, 13), ], coef = "S", quant = TRUE)
KB1_beta_div <- beta.div.comp(comm[c(6, 14), ], coef = "S", quant = TRUE)
KB4_beta_div <- beta.div.comp(comm[c(7, 15), ], coef = "S", quant = TRUE)
SWR2_beta_div <- beta.div.comp(comm[c(8, 16), ], coef = "S", quant = TRUE)

(ger_beta_div <- rbind(FB_beta_div, 
                       GO2_beta_div, GO4_beta_div, 
                       KR1_beta_div, KR4_beta_div, 
                       KB1_beta_div, KB4_beta_div, 
                       SWR2_beta_div))

## Sweden
A2_beta_div <- beta.div.comp(comm[c(17, 69), ], coef = "S", quant = TRUE)
A3_beta_div <- beta.div.comp(comm[c(18, 70), ], coef = "S", quant = TRUE)
A4_beta_div <- beta.div.comp(comm[c(19, 71), ], coef = "S", quant = TRUE)
A5_beta_div <- beta.div.comp(comm[c(20, 72), ], coef = "S", quant = TRUE)
A6_beta_div <- beta.div.comp(comm[c(21, 73), ], coef = "S", quant = TRUE)
A7_beta_div <- beta.div.comp(comm[c(22, 74), ], coef = "S", quant = TRUE)
A8_beta_div <- beta.div.comp(comm[c(23, 75), ], coef = "S", quant = TRUE)
A9_beta_div <- beta.div.comp(comm[c(24, 76), ], coef = "S", quant = TRUE)
ANN12_beta_div <- beta.div.comp(comm[c(25, 77), ], coef = "S", quant = TRUE)
ANN15_beta_div <- beta.div.comp(comm[c(26, 78), ], coef = "S", quant = TRUE)
ANN18_beta_div <- beta.div.comp(comm[c(27, 79), ], coef = "S", quant = TRUE)
ANN21_beta_div <- beta.div.comp(comm[c(28, 80), ], coef = "S", quant = TRUE)
ANN24_beta_div <- beta.div.comp(comm[c(29, 81), ], coef = "S", quant = TRUE)
ANN6_beta_div <- beta.div.comp(comm[c(30, 82), ], coef = "S", quant = TRUE)
ANN9_beta_div <- beta.div.comp(comm[c(31, 83), ], coef = "S", quant = TRUE)
ANO10_beta_div <- beta.div.comp(comm[c(32, 84), ], coef = "S", quant = TRUE)
ANO11_beta_div <- beta.div.comp(comm[c(33, 85), ], coef = "S", quant = TRUE)
ANO13_beta_div <- beta.div.comp(comm[c(34, 86), ], coef = "S", quant = TRUE)
ANO14_beta_div <- beta.div.comp(comm[c(35, 87), ], coef = "S", quant = TRUE)
ANO16_beta_div <- beta.div.comp(comm[c(36, 88), ], coef = "S", quant = TRUE)
ANO17_beta_div <- beta.div.comp(comm[c(37, 89), ], coef = "S", quant = TRUE)
ANO19_beta_div <- beta.div.comp(comm[c(38, 90), ], coef = "S", quant = TRUE)
ANO20_beta_div <- beta.div.comp(comm[c(39, 91), ], coef = "S", quant = TRUE)
ANO22_beta_div <- beta.div.comp(comm[c(40, 92), ], coef = "S", quant = TRUE)
ANO23_beta_div <- beta.div.comp(comm[c(41, 93), ], coef = "S", quant = TRUE)
ANO4XTR_beta_div <- beta.div.comp(comm[c(42, 94), ], coef = "S", quant = TRUE)
ANO5XTR_beta_div <- beta.div.comp(comm[c(43, 95), ], coef = "S", quant = TRUE)
ANO7_beta_div <- beta.div.comp(comm[c(44, 96), ], coef = "S", quant = TRUE)
ANO8_beta_div <- beta.div.comp(comm[c(45, 97), ], coef = "S", quant = TRUE)
B1_beta_div <- beta.div.comp(comm[c(46, 98), ], coef = "S", quant = TRUE)
B2_beta_div <- beta.div.comp(comm[c(47, 99), ], coef = "S", quant = TRUE)
B2XTR_beta_div <- beta.div.comp(comm[c(48, 100), ], coef = "S", quant = TRUE)
B3_beta_div <- beta.div.comp(comm[c(49, 101), ], coef = "S", quant = TRUE)
B3XTR_beta_div <- beta.div.comp(comm[c(50, 102), ], coef = "S", quant = TRUE)
B4_beta_div <- beta.div.comp(comm[c(51, 103), ], coef = "S", quant = TRUE)
B5_beta_div <- beta.div.comp(comm[c(52, 104), ], coef = "S", quant = TRUE)
BNN15_beta_div <- beta.div.comp(comm[c(53, 105), ], coef = "S", quant = TRUE)
BNN3_beta_div <- beta.div.comp(comm[c(54, 106), ], coef = "S", quant = TRUE)
BNN6_beta_div <- beta.div.comp(comm[c(55, 107), ], coef = "S", quant = TRUE)
BNN9_beta_div <- beta.div.comp(comm[c(56, 108), ], coef = "S", quant = TRUE)
BNO1_beta_div <- beta.div.comp(comm[c(57, 109), ], coef = "S", quant = TRUE)
BNO13_beta_div <- beta.div.comp(comm[c(58, 110), ], coef = "S", quant = TRUE)
BNO14_beta_div <- beta.div.comp(comm[c(59, 111), ], coef = "S", quant = TRUE)
BNO2_beta_div <- beta.div.comp(comm[c(60, 112), ], coef = "S", quant = TRUE)
BNO4_beta_div <- beta.div.comp(comm[c(61, 113), ], coef = "S", quant = TRUE)
BNO5_beta_div <- beta.div.comp(comm[c(62, 114), ], coef = "S", quant = TRUE)
BNO7_beta_div <- beta.div.comp(comm[c(63, 115), ], coef = "S", quant = TRUE)
BNO8_beta_div <- beta.div.comp(comm[c(64, 116), ], coef = "S", quant = TRUE)
C1_beta_div <- beta.div.comp(comm[c(65, 117), ], coef = "S", quant = TRUE)
CNN3_beta_div <- beta.div.comp(comm[c(66, 118), ], coef = "S", quant = TRUE)
CNO1_beta_div <- beta.div.comp(comm[c(67, 119), ], coef = "S", quant = TRUE)
CNO2_beta_div <- beta.div.comp(comm[c(68, 120), ], coef = "S", quant = TRUE)

(swe_beta_div <- rbind(A2_beta_div, A3_beta_div, A4_beta_div, A5_beta_div, A6_beta_div, A7_beta_div, A8_beta_div, A9_beta_div, 
                       ANN12_beta_div, ANN15_beta_div, ANN18_beta_div, ANN21_beta_div, ANN24_beta_div, ANN6_beta_div, ANN9_beta_div, 
                       ANO10_beta_div, ANO11_beta_div, ANO13_beta_div, ANO14_beta_div, ANO16_beta_div, ANO17_beta_div, ANO19_beta_div,
                       ANO20_beta_div, ANO22_beta_div, ANO23_beta_div, ANO4XTR_beta_div, ANO5XTR_beta_div, ANO7_beta_div, ANO8_beta_div,
                       B1_beta_div, B2_beta_div, B2XTR_beta_div, B3_beta_div, B3XTR_beta_div, B4_beta_div, B5_beta_div,
                       BNN15_beta_div, BNN3_beta_div, BNN6_beta_div, BNN9_beta_div, 
                       BNO1_beta_div, BNO13_beta_div, BNO14_beta_div, BNO2_beta_div, BNO4_beta_div, BNO5_beta_div, BNO7_beta_div, BNO8_beta_div, 
                       C1_beta_div, CNN3_beta_div, CNO1_beta_div, CNO2_beta_div))

## Finland
(Hang_beta_div <- beta.div.comp(comm[c(121, 131), ], coef = "S", quant = TRUE))
(Kant_beta_div <- beta.div.comp(comm[c(122, 132), ], coef = "S", quant = TRUE))
(Koti_beta_div <- beta.div.comp(comm[c(123, 133), ], coef = "S", quant = TRUE))
(Mati_beta_div <- beta.div.comp(comm[c(124, 134), ], coef = "S", quant = TRUE))
(Pess_beta_div <- beta.div.comp(comm[c(125, 135), ], coef = "S", quant = TRUE))
(Poro_beta_div <- beta.div.comp(comm[c(126, 136), ], coef = "S", quant = TRUE))
(Puta_beta_div <- beta.div.comp(comm[c(127, 137), ], coef = "S", quant = TRUE))
(Salm_beta_div <- beta.div.comp(comm[c(128, 138), ], coef = "S", quant = TRUE))
(Uopa_beta_div <- beta.div.comp(comm[c(129, 139), ], coef = "S", quant = TRUE))
(Vans_beta_div <- beta.div.comp(comm[c(130, 140), ], coef = "S", quant = TRUE))

(fin_beta_div <- rbind(Hang_beta_div, Kant_beta_div, Koti_beta_div, 
                       Mati_beta_div, Pess_beta_div, Poro_beta_div, 
                       Puta_beta_div, Salm_beta_div, Uopa_beta_div, 
                       Vans_beta_div))

(comm_beta_div <- as.data.frame(rbind(ger_beta_div, swe_beta_div, fin_beta_div)))
(rownames(comm_beta_div) <- rownames(comm[1:70, ]))
(comm_beta_div <- comm_beta_div[, c(1:3)])
colnames(comm_beta_div) <- c("repl", "rich_diff", "total_bd")
comm_beta_div <- apply(comm_beta_div, 2, as.character) # Used to unlist the variables & convert to dataframe
comm_beta_div <- as.data.frame(comm_beta_div)

#### FUNCTIONAL ALPHA DIVERSITY --------------------
### CREATING THE SUPRAREGIONAL FS (PCoA) #####
comm_traits <- comm_traits[(intersect(rownames(comm_traits), colnames(reg.mat))),]

# Gower dissimilarity
trait_df <- ktab.list.df(list(comm_traits[which(rowSums(comm_traits)==11),]))
tr.dist <- dist.ktab(trait_df, type= c("F"))

# Estimating the optimum number of dimensions
qual_fs <- quality_funct_space_fromdist(tr.dist, nbdim=15) 
qual_fs$meanSD # 0.006 at 9D

# Supraregional FS (PCoA)
supreg.pco <- dudi.pco(tr.dist, scan = F, nf = 9)

cumsum(supreg.pco$eig)[2]/sum(supreg.pco$eig)*100 # Variance explained by FS (PCoA); 68.1%

# Spearman rank correlations between original trait categories and functional space axes
cor.res <- round(cor(comm_traits[which(rowSums(comm_traits) == 11),], supreg.pco$li, method = "spearman"), 2) #11 trait groups used
cor.res

## FIGURE 3A: PROBABILITY DENSITY USING KERNEL DENSITY ESTIMATION #####
# Spearman rank correlations between original trait categories and functional space axes, which provides the contribution of each trait on the PCoA axes
cor.res
supreg.pco$li

# functional space
# Use dudi.pco for PCA for being consistent with scaling of scores in analysis
pco12 <- supreg.pco$li[,1:2]

# Kernel Density Estimation
Ho12 <- Hpi(x = pco12) # optimal bandwidth estimation
esto12 <- kde(x = pco12, H = Ho12, compute.cont = TRUE)     # kernel density estimation
plot(esto12)

# Set contour probabilities for drawing contour levels
clo12 <- contourLevels(esto12, prob = c(0.5, 0.05, 0.001), approx = TRUE)

# Setting arrows based on correlations
# pc1:All trait categories selected on axis1 -> life1 (0.56), life2 (-0.56), cycl2 (0.54), food4 (0.56), food8 (-0.64), fhab4 (0.59), fhab7 (-0.5)
# pc2:All trait categories selected on axis2-> cycl3 (0.52), repr3 (-0.57), loco4 (-0.51), loco7 (0.57), food3 (-0.50), fhab3 (-0.58)
cor.res # select traits with highest correlation (e.g. >= 0.5, <= -0.5)
comm_traits_sel <- subset(comm_traits, select = c(life1, life2, cycl2, food4, food8, fhab4, fhab7,
                                                  cycl3, repr3, loco4, loco7, food3, fhab3))
fit12 <- envfit(pco12, comm_traits_sel)  # use envfit to draw arrows, can be also done using trait loadings
fit122 <- fit12$vectors$arrows*-1 # drawing line segments in arrow opposites direction for pretty layout

# Plot Kernel Density Estimations
# pc1 and pc2
plot(esto12, cont = seq(1, 100, by = 1.79), display  = "filled.contour", add = FALSE, ylab = "PC2", xlab = "PC1", 
     cex.axis = 0.75, ylim = c(-1, 1), xlim = c(-1, 1) , las = 1)

pdf(file="FIGURE 3A - TRAIT PROBABILITY DENSITY PLOT.pdf",onefile = T, width = 6, height = 6)
# svg(file="Supraregional_FS_density_vectors.svg",onefile = T, width = 6, height = 6)
par(cex.lab = 1.25, cex.axis = 1.25, mfrow = c(1,1))
plot(esto12, cont = seq(1, 100, by = 1.79), display = "filled.contour", add = FALSE, 
     ylab = paste("PCo-2"), 
     xlab = paste("PCo-1"),
     main = paste("Supraregional FS"),
     #    ylab = paste("PC2","(",round(supreg.pco$eig[2]/sum(supreg.pco$eig), 3),"%)"), 
     #    xlab = paste("PC1","(",round(supreg.pco$eig[1]/sum(supreg.pco$eig), 3),"%)"), 
     cex.axis = 0.75, ylim=c(-1, 1), xlim = c(-1, 1), las = 1)
abline(h = 0, lty = 2, col = grey(0.5, alpha=0.2))
abline(v = 0, lty = 2, col = grey(0.5, alpha=0.2))
plot(esto12, abs.cont = clo12[1], labels = c(0.5), labcex = 0.75, add = TRUE, lwd = 0.75, col = "grey30")
plot(esto12, abs.cont = clo12[2], labels = c(0.95), labcex = 0.75, add = TRUE, lwd = 0.5, col = "grey60")
plot(esto12, abs.cont = clo12[3], labels = c(0.99), labcex = 0.75, add = TRUE, lwd = 0.5, col = "grey60")
points(pco12[,], pch = 16, cex = 0.25, col = "black")
plot(fit12, cex = 0.90, col = 1)
dev.off()

## FIGURE 3B: POSITION OF TAXONOMIC GROUPS WITHIN THE FS #####
library(adegraphics)
supra.gr <- data.frame(group2 = comm_traits_raw$Group2, comm_traits)[(intersect(rownames(comm_traits), colnames(reg.mat))),]

pdf(file = "FIGURE 3B - TAXA GROUPS PLOT.pdf", onefile = T, width = 6, height = 6)
# svg(file = "Supraregional_FS_taxa_grouping.svg", onefile = T, width = 4, height = 4)
par(mfrow = c(1,1), cex.axis = 1.85, cex.lab = 2, cex.main = 2, mar = c(5,5,4,1))
s.class(supreg.pco$li, fac = as.factor(supra.gr$group2[which(rowSums(comm_traits)==11)]), plines.col = 1:18, col = T)
dev.off()

## FIGURE 3C-D: REGIONAL TAXON POOL REPRESENTATION WITHIN THE SUPRAREGIONAL FSs #####
pdf(file = "FIGURE 3 - SUPRAREGIONAL FS (WITHOUT TAXA GROUPS PLOT).pdf",onefile = T, width = 12, height = 16) 
# svg(file = "Supraregional_and_regional_FS_convexhull.svg",onefile = T, width = 12, height = 16) 
par(mfrow = c(4,3), cex.axis = 1.85, cex.lab = 2, cex.main = 2, mar = c(5,5,3,1))

c("gold", "#228b22","royalblue4", 
  "gold", "#228b22","royalblue4", 
  "gold", "#228b22","royalblue4") -> col.reg
c("Germany - Overall FS", "Germany - Spring FS", "Germany - Autumn FS", 
  "Sweden - Overall FS", "Sweden - Spring FS", "Sweden - Autumn FS", 
  "Finland - Overall FS", "Finland - Spring FS", "Finland - Autumn FS") -> lab.reg

# Loop to plot the 6 FSs

plot(esto12, cont = seq(1, 100, by = 1.79), display = "filled.contour", add = FALSE, 
     ylab = paste("PCo-2"), 
     xlab = paste("PCo-1"),
     main = paste("Supraregional FS"),
     #    ylab = paste("PC2","(",round(supreg.pco$eig[2]/sum(supreg.pco$eig), 3),"%)"), 
     #    xlab = paste("PC1","(",round(supreg.pco$eig[1]/sum(supreg.pco$eig), 3),"%)"), 
     cex.axis = 1, ylim=c(-1, 1), xlim = c(-1, 1), las = 1)
abline(h = 0, lty = 2, col = grey(0.5, alpha=0.2))
abline(v = 0, lty = 2, col = grey(0.5, alpha=0.2))
plot(esto12, abs.cont = clo12[1], labels = c(0.5), labcex = 0.60, add = TRUE, lwd = 0.75, col = "grey30")
plot(esto12, abs.cont = clo12[2], labels = c(0.95), labcex = 0.60, add = TRUE, lwd = 0.5, col = "grey60")
plot(esto12, abs.cont = clo12[3], labels = c(0.99), labcex = 0.60, add = TRUE, lwd = 0.5, col = "grey60")
points(pco12[,], pch = 16, cex = 0.25, col = "black")
plot(fit12, cex = 0.90, col = 1)

plot(range(supreg.pco$li[1]), range(supreg.pco$li[2]),type = "n", main="Supraregional FS",xlab = "PCo-1", cex.axis = 1, ylab = "PCo-1")
points(supreg.pco$li[,c(1,2)],col="#4D4D4D",pch=".",cex=7)
plot_chull2D(supreg.pco$li[,c(1,2)],col_ch="#1E90FF50",border_ch="#4D4D4D")
# plot_chull2D(supreg.pco$li[(intersect(rownames(supreg.pco$li), colnames(reg.mat[which(colSums(reg.mat)>0)]))),c(1,2)],col_ch=adjustcolor("blue", alpha.f = 0.4),border_ch=adjustcolor("blue", alpha.f = 0.4))
colMeans(supreg.pco$li[,c(1,2)])->cent_ov
points(cent_ov[1],cent_ov[2],col="black",pch="+",cex=4)
colMeans(supreg.pco$li[(intersect(rownames(supreg.pco$li), colnames(reg.mat))),c(1,2)])->cent_r
points(cent_r[1],cent_r[2],col="red",pch="+",cex=2.5)

plot.new()

for (i in 1:9)
  
{
  plot(range(supreg.pco$li[1]), range(supreg.pco$li[2]),type = "n", main=lab.reg[i], cex.axis = 1, xlab= "PCo-1", ylab="PCo-1")
  points(supreg.pco$li[,c(1,2)],col="#4D4D4D",pch=".",cex=7)
  plot_chull2D(supreg.pco$li[,c(1,2)],col_ch="#1E90FF50",border_ch="#4D4D4D")
  plot_chull2D(supreg.pco$li[(intersect(rownames(supreg.pco$li), colnames(reg.mat[which(reg.mat[i,]>0)]))),c(1,2)],col_ch=adjustcolor(col.reg[i], alpha.f = 0.35),border_ch=adjustcolor(col.reg[i], alpha.f = 0.4))
  colMeans(supreg.pco$li[,c(1,2)])-> cent_ov
  points(cent_ov[1],cent_ov[2],col="black",pch="+",cex=4)
  colMeans(supreg.pco$li[(intersect(rownames(supreg.pco$li), colnames(reg.mat[which(reg.mat[i,]>0)]))),c(1,2)])-> cent_r
  points(cent_r[1],cent_r[2],col="red",pch="+",cex=2.5)
}

dev.off()

### COMPUTING FD METRICS (FRIC, FDIS, FEVE)#####
# Using supraregional-scale FS
chull_3d(supreg.pco$li, m = 6, prec = "Qt") -> ch_st # estimates the convex hull of a Functional Space
fric_3d(func_comm, supreg.pco$li, m = 6, prec = "QJ", fric.3d.max = ch_st) -> FRic_co
feve_k(supreg.pco$li, func_comm, m = 6) -> FEve_co
fdisp_k_sub(tr.dist, func_comm, tax_sub = colnames(func_comm), m = 6)$FDis -> FDis_co

# Substituting NAs by 0 (only if necessary)
FRic_co[which(is.na(FRic_co) == T)] <- 0
FEve_co[which(is.na(FEve_co) == T)] <- 0

### Add these indices to the div data frame
div$FRic <- FRic_co
div$FEve <- FEve_co
div$FDis <- FDis_co
div

## Computing FRed
func_comm
tr.dist
library(adiv)
# Calculated from the adiv package (Redundancy = 1-Uniqueness)
Redundancy <- uniqueness(func_comm, tr.dist, tol = 1e-08, abundance = TRUE) # creates 3 dataframs: kbar, V, and Red
# Isolate redundancy metrics
div$FRed <- Redundancy$red$R

### NULL MODELS OF FD METRICS #####
### shuffling null model: keeping the functional space the same, but shuffling taxa within the FS
# FRic
## inputs
func_comm
supreg.pco$li
ch_st # the maximum size of the convex hull

## model
FRic.shuff <- function(x){
  rownames(x) <- sample(rownames(x), length(rownames(x)), replace = F)
  x <- x[order(rownames(x)),]
  fric_3d(func_comm, x, m = 6, prec = "QJ", fric.3d.max = ch_st)
}
set.seed(1) # make results repeatable
FRic.obs.null.output.shuff <- cbind(fric_3d(func_comm, supreg.pco$li, m = 6, prec = "QJ", fric.3d.max = ch_st), 
                                    replicate(3, FRic.shuff(supreg.pco$li))) # change this number to 999 after you have determined the code runs

FRic.ses.value <- (FRic.obs.null.output.shuff[,1] - 
                     apply(FRic.obs.null.output.shuff, 1, mean))/
  apply(FRic.obs.null.output.shuff,1,sd)

qFRic <- NaN*FRic.obs.null.output.shuff[,1]

for(i in seq(qFRic)){
  qFRic[i] <- sum(FRic.obs.null.output.shuff[,1][i] > FRic.obs.null.output.shuff[i,]) / length(FRic.obs.null.output.shuff[i,])
}

sigFRic <- qFRic < 0.05 | qFRic > 0.95 # test if outside distribution

FRic_output <- as.data.frame(cbind(FRic.ses.value, sigFRic))
colnames(FRic_output) <- c("FRic.SES", "FRic.SES.sig")
boxplot(FRic.obs.null.output.shuff[,1], FRic_output[,1])#obsevred vs standardised FRic

null_outputs <- FRic_output

# write.csv(FRic_output, "FRic_null_output.csv")

# FEve
## inputs
func_comm
supreg.pco$li

## model
FEve.shuff <- function(x){
  rownames(x) <- sample(rownames(x), length(rownames(x)), replace = F)
  x <- x[order(rownames(x)),]
  feve_k(x, func_comm, m = 6)
}
set.seed(1) # make results repeatable
FEve.obs.null.output.shuff <- cbind(feve_k(supreg.pco$li, func_comm, m = 6), 
                                    replicate(3, FEve.shuff(supreg.pco$li))) # change this number to 999 after you have determined the code runs

FEve.ses.value <- (FEve.obs.null.output.shuff[,1] - 
                     apply(FEve.obs.null.output.shuff, 1, mean))/
  apply(FEve.obs.null.output.shuff,1,sd)

qFEve <- NaN*FEve.obs.null.output.shuff[,1]

for(i in seq(qFEve)){
  qFEve[i] <- sum(FEve.obs.null.output.shuff[,1][i] > FEve.obs.null.output.shuff[i,]) / length(FEve.obs.null.output.shuff[i,])
}

sigFEve <- qFEve < 0.05 | qFEve > 0.95 # test if outside distribution

FEve_output <- as.data.frame(cbind(FEve.ses.value, sigFEve))
colnames(FEve_output) <- c("FEve.SES", "FEve.SES.sig")
boxplot(FEve.obs.null.output.shuff[,1], FEve_output[,1])#obsevred vs standardised FRic

null_outputs <- cbind(FRic_output, FEve_output)

# write.csv(FEve_output, "FEve_null_output.csv")

# FDisp
## inputs
func_comm
tr.dist

## model - changed species identities by shuffling the site-by-taxa matrix
FDis.shuff <- function(x){
  colnames(x) <- sample(colnames(x), length(colnames(x)), replace = F)
  x <- x[, order(names(x))]
  colnames(x) <- colnames(func_comm)
  fdisp_k(tr.dist, x, m = 6)$FDis
}
set.seed(1) # make results repeatable
FDis.obs.null.output.shuff <- cbind(fdisp_k(tr.dist, func_comm, m = 6)$FDis, 
                                    replicate(3, FDis.shuff(func_comm))) # change this number to 999 after you have determined the code runs

FDis.ses.value <- (FDis.obs.null.output.shuff[,1] - 
                     apply(FDis.obs.null.output.shuff, 1, mean))/
  apply(FDis.obs.null.output.shuff,1,sd)

qFDis <- NaN*FDis.obs.null.output.shuff[,1]

for(i in seq(qFDis)){
  qFDis[i] <- sum(FDis.obs.null.output.shuff[,1][i] > FDis.obs.null.output.shuff[i,]) / length(FDis.obs.null.output.shuff[i,])
}

sigFDis <- qFDis < 0.05 | qFDis > 0.95 # test if outside distribution

FDis_output <- as.data.frame(cbind(FDis.ses.value, sigFDis))
colnames(FDis_output) <- c("FDis.SES", "FDis.SES.sig")
boxplot(FDis.obs.null.output.shuff[,1], FDis_output[,1])# obsevred vs standardised FDis

null_outputs <- cbind(FRic_output, FEve_output, FDis_output)

# Substituting NAs by 0 (only if necessary)
null_outputs$FRic.SES[which(is.na(null_outputs$FRic.SES) == T)] <- 0
null_outputs$FEve.SES[which(is.na(null_outputs$FEve.SES) == T)] <- 0

# write.csv(FDis_output, "FDis_null_output.csv")
# write.csv(null_outputs, "null_model_outputs.csv")

# Read in computed null model values
null.outputs.combined <- read.csv("func_null_model_outputs.csv", h = TRUE, sep = ",", row.names = c(1) ,stringsAsFactors = FALSE)

# Substituting NAs by 0 (only if necessary)
null.outputs.combined$FRic.SES[which(is.na(null.outputs.combined$FRic.SES) == T)] <- 0
null.outputs.combined$FEve.SES[which(is.na(null.outputs.combined$FEve.SES) == T)] <- 0

#### FUNCTIONAL BETA DIVERSITY --------------------
library(BAT)
### Functional dissimilarity matrix (Gower dissimilarity)
(comm_traits_gowd <- gowdis(comm_traits))
(comm_traits_gd_mat <- as.matrix(gowdis(comm_traits)))
## Plot the tree and the dendrogram
comm_traits_tree <- hclust(comm_traits_gowd, "ward.D2")
plot(comm_traits_tree, hang = -1, main = "")
## Compute functional dissimilarity matrix with abundance informtation
(func_dissim <- beta(func_comm, comm_traits_tree, func = "Soerensen", abund = TRUE))
(func_dissim_mat <- func_dissim$Btotal)

#### Separate matrices for each region
### Germany 
## Combined seasons
ger_func_comm_full <- func_comm[1:16, ]
ger_func_comm_full <- ger_func_comm_full[, colSums(ger_func_comm_full) > 0] # only keeps columns with column sums > 0
# functional trait data
spp.ger.full = (row.names(comm_traits) %in% colnames(ger_func_comm_full)) # make a list of T or F is there is a match in spp
ger_traits.full = comm_traits[spp.ger.full, ] # subset by this list
# Check data
identical(row.names(ger_traits.full), colnames(ger_func_comm_full)) # check if names are the same
# Compute functional dissimilarity matrix
ger_gower_full <- gowdis(ger_traits.full) # gower disimilarity matrix
plot(ger_clust_full <- hclust(ger_gower_full, method = "ward.D"))
ger_func_beta_full <- beta(ger_func_comm_full, ger_clust_full, func = "Soerensen", abund = TRUE) # percentage difference
(ger_func_mat_full <- ger_func_beta_full$Btotal)
## Spring
ger_func_comm_spr <- func_comm[1:8, ]
ger_func_comm_spr <- ger_func_comm_spr[, colSums(ger_func_comm_spr) > 0] # only keeps columns with column sums > 0
# functional trait data
spp.ger.spr = (row.names(comm_traits) %in% colnames(ger_func_comm_spr)) # make a list of T or F is there is a match in spp
ger_traits.spr = comm_traits[spp.ger.spr, ] # subset by this list
# Check data
identical(row.names(ger_traits.spr), colnames(ger_func_comm_spr)) # check if names are the same
# Compute functional dissimilarity matrix
ger_gower_spr <- gowdis(ger_traits.spr) # gower disimilarity matrix
plot(ger_clust_spr <- hclust(ger_gower_spr, method = "ward.D"))
ger_func_beta_spr <- beta(ger_func_comm_spr, ger_clust_spr, func = "Soerensen", abund = TRUE) # percentage difference
(ger_func_mat_spr <- ger_func_beta_spr$Btotal)
## Autumn
ger_func_comm_aut <- func_comm[9:16, ]
ger_func_comm_aut <- ger_func_comm_aut[, colSums(ger_func_comm_aut) > 0] # only keeps columns with column sums > 0
# functional trait data
spp.ger.aut = (row.names(comm_traits) %in% colnames(ger_func_comm_aut)) # make a list of T or F is there is a match in spp
ger_traits.aut = comm_traits[spp.ger.aut, ] # subset by this list
# Check data
identical(row.names(ger_traits.aut), colnames(ger_func_comm_aut)) # check if names are the same
# Compute functional dissimilarity matrix
ger_gower_aut <- gowdis(ger_traits.aut) # gower disimilarity matrix
plot(ger_clust_aut <- hclust(ger_gower_aut, method = "ward.D"))
ger_func_beta_aut <- beta(ger_func_comm_aut, ger_clust_aut, func = "Soerensen", abund = TRUE) # percentage difference
(ger_func_mat_aut <- ger_func_beta_aut$Btotal)

### Sweden 
## Combined seasons
swe_func_comm_full <- func_comm[17:120, ]
swe_func_comm_full <- swe_func_comm_full[, colSums(swe_func_comm_full) > 0] # only keeps columns with column sums > 0
# functional trait data
spp.swe.full = (row.names(comm_traits) %in% colnames(swe_func_comm_full)) # make a list of T or F is there is a match in spp
swe_traits.full = comm_traits[spp.swe.full, ] # subset by this list
# Check data
identical(row.names(swe_traits.full), colnames(swe_func_comm_full)) # check if names are the same
# Compute functional dissimilarity matrix
swe_gower_full <- gowdis(swe_traits.full) # gower disimilarity matrix
plot(swe_clust_full <- hclust(swe_gower_full, method = "ward.D"))
swe_func_beta_full <- beta(swe_func_comm_full, swe_clust_full, func = "Soerensen", abund = TRUE) # percentage difference
(swe_func_mat_full <- swe_func_beta_full$Btotal)
## Spring
swe_func_comm_spr <- func_comm[17:68, ]
swe_func_comm_spr <- swe_func_comm_spr[, colSums(swe_func_comm_spr) > 0] # only keeps columns with column sums > 0
# functional trait data
spp.swe.spr = (row.names(comm_traits) %in% colnames(swe_func_comm_spr)) # make a list of T or F is there is a match in spp
swe_traits.spr = comm_traits[spp.swe.spr, ] # subset by this list
# Check data
identical(row.names(swe_traits.spr), colnames(swe_func_comm_spr)) # check if names are the same
# Compute functional dissimilarity matrix
swe_gower_spr <- gowdis(swe_traits.spr) # gower disimilarity matrix
plot(swe_clust_spr <- hclust(swe_gower_spr, method = "ward.D"))
swe_func_beta_spr <- beta(swe_func_comm_spr, swe_clust_spr, func = "Soerensen", abund = TRUE) # percentage difference
(swe_func_mat_spr <- swe_func_beta_spr$Btotal)
## Autumn
swe_func_comm_aut <- func_comm[69:120, ]
swe_func_comm_aut <- swe_func_comm_aut[, colSums(swe_func_comm_aut) > 0] # only keeps columns with column sums > 0
# functional trait data
spp.swe.aut = (row.names(comm_traits) %in% colnames(swe_func_comm_aut)) # make a list of T or F is there is a match in spp
swe_traits.aut = comm_traits[spp.swe.aut, ] # subset by this list
# Check data
identical(row.names(swe_traits.aut), colnames(swe_func_comm_aut)) # check if names are the same
# Compute functional dissimilarity matrix
swe_gower_aut <- gowdis(swe_traits.aut) # gower disimilarity matrix
plot(swe_clust_aut <- hclust(swe_gower_aut, method = "ward.D"))
swe_func_beta_aut <- beta(swe_func_comm_aut, swe_clust_aut, func = "Soerensen", abund = TRUE) # percentage difference
(swe_func_mat_aut <- swe_func_beta_aut$Btotal)

### Finland 
## Combined seasons
fin_func_comm_full <- func_comm[121:140, ]
fin_func_comm_full <- fin_func_comm_full[, colSums(fin_func_comm_full) > 0] # only keeps columns with column sums > 0
# functional trait data
spp.fin.full = (row.names(comm_traits) %in% colnames(fin_func_comm_full)) # make a list of T or F is there is a match in spp
fin_traits.full = comm_traits[spp.fin.full, ] # subset by this list
# Check data
identical(row.names(fin_traits.full), colnames(fin_func_comm_full)) # check if names are the same
# Compute functional dissimilarity matrix
fin_gower_full <- gowdis(fin_traits.full) # gower disimilarity matrix
plot(fin_clust_full <- hclust(fin_gower_full, method = "ward.D"))
fin_func_beta_full <- beta(fin_func_comm_full, fin_clust_full, func = "Soerensen", abund = TRUE) # percentage difference
(fin_func_mat_full <- fin_func_beta_full$Btotal)
## Spring
fin_func_comm_spr <- func_comm[121:130, ]
fin_func_comm_spr <- fin_func_comm_spr[, colSums(fin_func_comm_spr) > 0] # only keeps columns with column sums > 0
# functional trait data
spp.fin.spr = (row.names(comm_traits) %in% colnames(fin_func_comm_spr)) # make a list of T or F is there is a match in spp
fin_traits.spr = comm_traits[spp.fin.spr, ] # subset by this list
# Check data
identical(row.names(fin_traits.spr), colnames(fin_func_comm_spr)) # check if names are the same
# Compute functional dissimilarity matrix
fin_gower_spr <- gowdis(fin_traits.spr) # gower disimilarity matrix
plot(fin_clust_spr <- hclust(fin_gower_spr, method = "ward.D"))
fin_func_beta_spr <- beta(fin_func_comm_spr, fin_clust_spr, func = "Soerensen", abund = TRUE) # percentage difference
(fin_func_mat_spr <- fin_func_beta_spr$Btotal)
## Autumn
fin_func_comm_aut <- func_comm[131:140, ]
fin_func_comm_aut <- fin_func_comm_aut[, colSums(fin_func_comm_aut) > 0] # only keeps columns with column sums > 0
# functional trait data
spp.fin.aut = (row.names(comm_traits) %in% colnames(fin_func_comm_aut)) # make a list of T or F is there is a match in spp
fin_traits.aut = comm_traits[spp.fin.aut, ] # subset by this list
# Check data
identical(row.names(fin_traits.aut), colnames(fin_func_comm_aut)) # check if names are the same
# Compute functional dissimilarity matrix
fin_gower_aut <- gowdis(fin_traits.aut) # gower disimilarity matrix
plot(fin_clust_aut <- hclust(fin_gower_aut, method = "ward.D"))
fin_func_beta_aut <- beta(fin_func_comm_aut, fin_clust_aut, func = "Soerensen", abund = TRUE) # percentage difference
(fin_func_mat_aut <- fin_func_beta_aut$Btotal)

#### FUNCTIONAL TURNOVER --------------------
comm_CWM <- dbFD(comm_traits, as.matrix(func_comm), calc.CWM = TRUE, CWM.type = "all")$CWM
### Podani family, percentage difference (Sorensen dissimilarity on quantitative data)
## Germany
FB_func_beta_div <- beta.div.comp(comm_CWM[c(1, 9), ], coef = "S", quant = TRUE)
GO2_func_beta_div <- beta.div.comp(comm_CWM[c(2, 10), ], coef = "S", quant = TRUE)
GO4_func_beta_div <- beta.div.comp(comm_CWM[c(3, 11), ], coef = "S", quant = TRUE)
KR1_func_beta_div <- beta.div.comp(comm_CWM[c(4, 12), ], coef = "S", quant = TRUE)
KR4_func_beta_div <- beta.div.comp(comm_CWM[c(5, 13), ], coef = "S", quant = TRUE)
KB1_func_beta_div <- beta.div.comp(comm_CWM[c(6, 14), ], coef = "S", quant = TRUE)
KB4_func_beta_div <- beta.div.comp(comm_CWM[c(7, 15), ], coef = "S", quant = TRUE)
SWR2_func_beta_div <- beta.div.comp(comm_CWM[c(8, 16), ], coef = "S", quant = TRUE)

(ger_func_beta_div <- rbind(FB_func_beta_div, 
                            GO2_func_beta_div, GO4_func_beta_div, 
                            KR1_func_beta_div, KR4_func_beta_div, 
                            KB1_func_beta_div, KB4_func_beta_div, 
                            SWR2_func_beta_div))

## Sweden
A2_func_beta_div <- beta.div.comp(comm_CWM[c(17, 69), ], coef = "S", quant = TRUE)
A3_func_beta_div <- beta.div.comp(comm_CWM[c(18, 70), ], coef = "S", quant = TRUE)
A4_func_beta_div <- beta.div.comp(comm_CWM[c(19, 71), ], coef = "S", quant = TRUE)
A5_func_beta_div <- beta.div.comp(comm_CWM[c(20, 72), ], coef = "S", quant = TRUE)
A6_func_beta_div <- beta.div.comp(comm_CWM[c(21, 73), ], coef = "S", quant = TRUE)
A7_func_beta_div <- beta.div.comp(comm_CWM[c(22, 74), ], coef = "S", quant = TRUE)
A8_func_beta_div <- beta.div.comp(comm_CWM[c(23, 75), ], coef = "S", quant = TRUE)
A9_func_beta_div <- beta.div.comp(comm_CWM[c(24, 76), ], coef = "S", quant = TRUE)
ANN12_func_beta_div <- beta.div.comp(comm_CWM[c(25, 77), ], coef = "S", quant = TRUE)
ANN15_func_beta_div <- beta.div.comp(comm_CWM[c(26, 78), ], coef = "S", quant = TRUE)
ANN18_func_beta_div <- beta.div.comp(comm_CWM[c(27, 79), ], coef = "S", quant = TRUE)
ANN21_func_beta_div <- beta.div.comp(comm_CWM[c(28, 80), ], coef = "S", quant = TRUE)
ANN24_func_beta_div <- beta.div.comp(comm_CWM[c(29, 81), ], coef = "S", quant = TRUE)
ANN6_func_beta_div <- beta.div.comp(comm_CWM[c(30, 82), ], coef = "S", quant = TRUE)
ANN9_func_beta_div <- beta.div.comp(comm_CWM[c(31, 83), ], coef = "S", quant = TRUE)
ANO10_func_beta_div <- beta.div.comp(comm_CWM[c(32, 84), ], coef = "S", quant = TRUE)
ANO11_func_beta_div <- beta.div.comp(comm_CWM[c(33, 85), ], coef = "S", quant = TRUE)
ANO13_func_beta_div <- beta.div.comp(comm_CWM[c(34, 86), ], coef = "S", quant = TRUE)
ANO14_func_beta_div <- beta.div.comp(comm_CWM[c(35, 87), ], coef = "S", quant = TRUE)
ANO16_func_beta_div <- beta.div.comp(comm_CWM[c(36, 88), ], coef = "S", quant = TRUE)
ANO17_func_beta_div <- beta.div.comp(comm_CWM[c(37, 89), ], coef = "S", quant = TRUE)
ANO19_func_beta_div <- beta.div.comp(comm_CWM[c(38, 90), ], coef = "S", quant = TRUE)
ANO20_func_beta_div <- beta.div.comp(comm_CWM[c(39, 91), ], coef = "S", quant = TRUE)
ANO22_func_beta_div <- beta.div.comp(comm_CWM[c(40, 92), ], coef = "S", quant = TRUE)
ANO23_func_beta_div <- beta.div.comp(comm_CWM[c(41, 93), ], coef = "S", quant = TRUE)
ANO4XTR_func_beta_div <- beta.div.comp(comm_CWM[c(42, 94), ], coef = "S", quant = TRUE)
ANO5XTR_func_beta_div <- beta.div.comp(comm_CWM[c(43, 95), ], coef = "S", quant = TRUE)
ANO7_func_beta_div <- beta.div.comp(comm_CWM[c(44, 96), ], coef = "S", quant = TRUE)
ANO8_func_beta_div <- beta.div.comp(comm_CWM[c(45, 97), ], coef = "S", quant = TRUE)
B1_func_beta_div <- beta.div.comp(comm_CWM[c(46, 98), ], coef = "S", quant = TRUE)
B2_func_beta_div <- beta.div.comp(comm_CWM[c(47, 99), ], coef = "S", quant = TRUE)
B2XTR_func_beta_div <- beta.div.comp(comm_CWM[c(48, 100), ], coef = "S", quant = TRUE)
B3_func_beta_div <- beta.div.comp(comm_CWM[c(49, 101), ], coef = "S", quant = TRUE)
B3XTR_func_beta_div <- beta.div.comp(comm_CWM[c(50, 102), ], coef = "S", quant = TRUE)
B4_func_beta_div <- beta.div.comp(comm_CWM[c(51, 103), ], coef = "S", quant = TRUE)
B5_func_beta_div <- beta.div.comp(comm_CWM[c(52, 104), ], coef = "S", quant = TRUE)
BNN15_func_beta_div <- beta.div.comp(comm_CWM[c(53, 105), ], coef = "S", quant = TRUE)
BNN3_func_beta_div <- beta.div.comp(comm_CWM[c(54, 106), ], coef = "S", quant = TRUE)
BNN6_func_beta_div <- beta.div.comp(comm_CWM[c(55, 107), ], coef = "S", quant = TRUE)
BNN9_func_beta_div <- beta.div.comp(comm_CWM[c(56, 108), ], coef = "S", quant = TRUE)
BNO1_func_beta_div <- beta.div.comp(comm_CWM[c(57, 109), ], coef = "S", quant = TRUE)
BNO13_func_beta_div <- beta.div.comp(comm_CWM[c(58, 110), ], coef = "S", quant = TRUE)
BNO14_func_beta_div <- beta.div.comp(comm_CWM[c(59, 111), ], coef = "S", quant = TRUE)
BNO2_func_beta_div <- beta.div.comp(comm_CWM[c(60, 112), ], coef = "S", quant = TRUE)
BNO4_func_beta_div <- beta.div.comp(comm_CWM[c(61, 113), ], coef = "S", quant = TRUE)
BNO5_func_beta_div <- beta.div.comp(comm_CWM[c(62, 114), ], coef = "S", quant = TRUE)
BNO7_func_beta_div <- beta.div.comp(comm_CWM[c(63, 115), ], coef = "S", quant = TRUE)
BNO8_func_beta_div <- beta.div.comp(comm_CWM[c(64, 116), ], coef = "S", quant = TRUE)
C1_func_beta_div <- beta.div.comp(comm_CWM[c(65, 117), ], coef = "S", quant = TRUE)
CNN3_func_beta_div <- beta.div.comp(comm_CWM[c(66, 118), ], coef = "S", quant = TRUE)
CNO1_func_beta_div <- beta.div.comp(comm_CWM[c(67, 119), ], coef = "S", quant = TRUE)
CNO2_func_beta_div <- beta.div.comp(comm_CWM[c(68, 120), ], coef = "S", quant = TRUE)

(swe_func_beta_div <- rbind(A2_func_beta_div, A3_func_beta_div, A4_func_beta_div, A5_func_beta_div, A6_func_beta_div, A7_func_beta_div, A8_func_beta_div, A9_func_beta_div, 
                            ANN12_func_beta_div, ANN15_func_beta_div, ANN18_func_beta_div, ANN21_func_beta_div, ANN24_func_beta_div, ANN6_func_beta_div, ANN9_func_beta_div, 
                            ANO10_func_beta_div, ANO11_func_beta_div, ANO13_func_beta_div, ANO14_func_beta_div, ANO16_func_beta_div, ANO17_func_beta_div, ANO19_func_beta_div,
                            ANO20_func_beta_div, ANO22_func_beta_div, ANO23_func_beta_div, ANO4XTR_func_beta_div, ANO5XTR_func_beta_div, ANO7_func_beta_div, ANO8_func_beta_div,
                            B1_func_beta_div, B2_func_beta_div, B2XTR_func_beta_div, B3_func_beta_div, B3XTR_func_beta_div, B4_func_beta_div, B5_func_beta_div,
                            BNN15_func_beta_div, BNN3_func_beta_div, BNN6_func_beta_div, BNN9_func_beta_div, 
                            BNO1_func_beta_div, BNO13_func_beta_div, BNO14_func_beta_div, BNO2_func_beta_div, BNO4_func_beta_div, BNO5_func_beta_div, BNO7_func_beta_div, BNO8_func_beta_div, 
                            C1_func_beta_div, CNN3_func_beta_div, CNO1_func_beta_div, CNO2_func_beta_div))

## Finland
(Hang_func_beta_div <- beta.div.comp(comm_CWM[c(121, 131), ], coef = "S", quant = TRUE))
(Kant_func_beta_div <- beta.div.comp(comm_CWM[c(122, 132), ], coef = "S", quant = TRUE))
(Koti_func_beta_div <- beta.div.comp(comm_CWM[c(123, 133), ], coef = "S", quant = TRUE))
(Mati_func_beta_div <- beta.div.comp(comm_CWM[c(124, 134), ], coef = "S", quant = TRUE))
(Pess_func_beta_div <- beta.div.comp(comm_CWM[c(125, 135), ], coef = "S", quant = TRUE))
(Poro_func_beta_div <- beta.div.comp(comm_CWM[c(126, 136), ], coef = "S", quant = TRUE))
(Puta_func_beta_div <- beta.div.comp(comm_CWM[c(127, 137), ], coef = "S", quant = TRUE))
(Salm_func_beta_div <- beta.div.comp(comm_CWM[c(128, 138), ], coef = "S", quant = TRUE))
(Uopa_func_beta_div <- beta.div.comp(comm_CWM[c(129, 139), ], coef = "S", quant = TRUE))
(Vans_func_beta_div <- beta.div.comp(comm_CWM[c(130, 140), ], coef = "S", quant = TRUE))

(fin_func_beta_div <- rbind(Hang_func_beta_div, Kant_func_beta_div, Koti_func_beta_div, 
                            Mati_func_beta_div, Pess_func_beta_div, Poro_func_beta_div, 
                            Puta_func_beta_div, Salm_func_beta_div, Uopa_func_beta_div, 
                            Vans_func_beta_div))

## Combine the datasets
(comm_func_beta_div <- as.data.frame(rbind(ger_func_beta_div, swe_func_beta_div, fin_func_beta_div)))
(rownames(comm_func_beta_div) <- rownames(comm_CWM[1:70, ]))
(comm_func_beta_div <- comm_func_beta_div[, c(1:3)])
colnames(comm_func_beta_div) <- c("repl", "rich_diff", "total_bd")
comm_func_beta_div <- apply(comm_func_beta_div, 2, as.character) # Used to unlist the variables & convert to dataframe
comm_func_beta_div <- as.data.frame(comm_func_beta_div)

# Create turnover database
# change values to numeric
comm_beta_div$repl <- as.numeric(comm_beta_div$repl)
comm_beta_div$rich_diff <- as.numeric(comm_beta_div$rich_diff)
comm_beta_div$total_bd <- as.numeric(comm_beta_div$total_bd)

comm_func_beta_div$repl <- as.numeric(comm_func_beta_div$repl)
comm_func_beta_div$rich_diff <- as.numeric(comm_func_beta_div$rich_diff)
comm_func_beta_div$total_bd <- as.numeric(comm_func_beta_div$total_bd)

# create database
turnover <- cbind(comm_beta_div$repl, comm_func_beta_div$repl)
colnames(turnover) <- c("Tturn", "Fturn")
turnover <- as.data.frame(turnover)

#### FIGURE 2: BOXPLOTS OF DIVERSITY INDICES ----------------
### Datasets used
div
null <- null.outputs.combined[, c(1, 3, 5)] # Regional null models are used as they are computed on the regional species pools

### Panel parametres
pdf("FIGURE 2 - BOXPLOTS.pdf", width = 8, height = 16, pointsize = 12, onefile = T)
# svg("Boxplots - All diversity metrics - vertical.svg", width = 8, height = 16, pointsize = 12, onefile = T)
par(mfrow = c(5, 2), cex.axis = 1.85, cex.lab = 2, cex.main = 2, mar = c(5,5,2,1))

### Taxonomic richness
# boxplot(div$N0[1:8], div$N0[9:16], div$N0[17:68], div$N0[69:120], div$N0[121:130], div$N0[131:140]) # test to check axes
(TRic_plot <- boxplot(div$N0[1:8], div$N0[9:16], div$N0[17:68], div$N0[69:120], div$N0[121:130], div$N0[131:140], 
                      yaxt = "n", xaxt = "n",  # should the axes be displayed
                      ylim = c(0, 60), # limits of the y-axis
                      ylab = "TRic", xlab = "", # names of the axis labels
                      main = "A",
                      cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, # size of labels and axes
                      horizontal = FALSE, #this argument reverses the axis orientation
                      varwidth = FALSE, # varwidth makes the boxplot widths proportional to the sample size
                      notch = FALSE, # this makes a notch in the boxplot. In the notched boxplot, if two boxes' notches do not overlap this is 'strong evidence' their medians differ (Chambers et al., 1983, p. 62).
                      border = TRUE,
                      outline = FALSE,
                      cex = 1, boxlty = 1, col = c("#638a634d", "#00688c4d", "#638a634d", "#00688c4d", "#638a634d", "#00688c4d")))
axis(1, col = "black", col.axis = "black", at = seq(1, 6, 1), cex.axis = 1.5, labels = FALSE) # X-axis
axis(2, col = "black", col.axis = "black", at = seq(0, 60, 10), cex.axis = 1.5) # Y-axis
x <- c(1:6)
labs <- c("Ger/Spr", "Ger/Aut", 
          "Swe/Spr", "Swe/Aut", 
          "Fin/Spr", "Fin/Aut")
# text(cex = 1.5, x = x, y = -10, labs, xpd = TRUE, srt = 45, font = 1)

### Taxonomic evenness
#boxplot(div$E10[1:8], div$E10[9:16], div$E10[17:68], div$E10[69:120], div$E10[121:130], div$E10[131:140]) # test to check axes
(TEve_plot <- boxplot(div$E10[1:8], div$E10[9:16], div$E10[17:68], div$E10[69:120], div$E10[121:130], div$E10[131:140], 
                      yaxt = "n", xaxt = "n",  # should the axes be displayed
                      ylim = c(0.0, 1.0), # limits of the y-axis
                      ylab = "TEve", xlab = "", # names of the axis labels
                      main = "B",
                      cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, # size of labels and axes
                      horizontal = FALSE, #this argument reverses the axis orientation
                      varwidth = FALSE, # varwidth makes the boxplot widths proportional to the sample size
                      notch = FALSE, # this makes a notch in the boxplot. In the notched boxplot, if two boxes' notches do not overlap this is 'strong evidence' their medians differ (Chambers et al., 1983, p. 62).
                      border = TRUE,
                      outline = FALSE,
                      cex = 1, boxlty = 1, col = c("#638a634d", "#00688c4d", "#638a634d", "#00688c4d", "#638a634d", "#00688c4d")))
axis(1, col = "black", col.axis = "black", at = seq(1, 6, 1), cex.axis = 1.5, labels = FALSE) # X-axis
axis(2, col = "black", col.axis = "black", at = seq(0.0, 1.0, 0.2), cex.axis = 1.5) # Y-axis
x <- c(1:6)
labs <- c("Ger/Spr", "Ger/Aut", 
          "Swe/Spr", "Swe/Aut", 
          "Fin/Spr", "Fin/Aut")
# text(cex = 1.5, x = x, y = -10, labs, xpd = TRUE, srt = 45, font = 1)

### Shannon diversity
#boxplot(div$N1[1:8], div$N1[9:16], div$N1[17:68], div$N1[69:120], div$N1[121:130], div$N1[131:140]) # test to check axes
(Shan_plot <- boxplot(div$N1[1:8], div$N1[9:16], div$N1[17:68], div$N1[69:120], div$N1[121:130], div$N1[131:140], 
                      yaxt = "n", xaxt = "n",  # should the axes be displayed
                      ylim = c(0, 20), # limits of the y-axis
                      ylab = "Shan", xlab = "", # names of the axis labels
                      main = "C",
                      cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, # size of labels and axes
                      horizontal = FALSE, #this argument reverses the axis orientation
                      varwidth = FALSE, # varwidth makes the boxplot widths proportional to the sample size
                      notch = FALSE, # this makes a notch in the boxplot. In the notched boxplot, if two boxes' notches do not overlap this is 'strong evidence' their medians differ (Chambers et al., 1983, p. 62).
                      border = TRUE,
                      outline = FALSE,
                      cex = 1, boxlty = 1, col = c("#638a634d", "#00688c4d", "#638a634d", "#00688c4d", "#638a634d", "#00688c4d")))
axis(1, col = "black", col.axis = "black", at = seq(1, 6, 1), cex.axis = 1.5, labels = FALSE) # X-axis
axis(2, col = "black", col.axis = "black", at = seq(0, 20, 4), cex.axis = 1.5) # Y-axis
x <- c(1:6)
labs <- c("Ger/Spr", "Ger/Aut", 
          "Swe/Spr", "Swe/Aut", 
          "Fin/Spr", "Fin/Aut")
# text(cex = 1.5, x = x, y = -10, labs, xpd = TRUE, srt = 45, font = 1)

### Overall abundance (log(Abund))
# boxplot(log(div$A[1:8]), log(div$A[9:16]), log(div$A[17:68]), log(div$A[69:120]), log(div$A[121:130]), log(div$A[131:140])) # test to check axes
(Abund_plot <- boxplot(log(div$A[1:8]), log(div$A[9:16]), log(div$A[17:68]), log(div$A[69:120]), log(div$A[121:130]), log(div$A[131:140]), 
                       yaxt = "n", xaxt = "n",  # should the axes be displayed
                       ylim = c(3, 10), # limits of the y-axis
                       ylab = "log(Abund)", xlab = "", # names of the axis labels
                       main = "D",
                       cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, # size of labels and axes
                       horizontal = FALSE, #this argument reverses the axis orientation
                       varwidth = FALSE, # varwidth makes the boxplot widths proportional to the sample size
                       notch = FALSE, # this makes a notch in the boxplot. In the notched boxplot, if two boxes' notches do not overlap this is 'strong evidence' their medians differ (Chambers et al., 1983, p. 62).
                       border = TRUE,
                       outline = FALSE,
                       cex = 1, boxlty = 1, col = c("#638a634d", "#00688c4d", "#638a634d", "#00688c4d", "#638a634d", "#00688c4d")))
axis(1, col = "black", col.axis = "black", at = seq(1, 6, 1), cex.axis = 1.5, labels = FALSE) # X-axis
axis(2, col = "black", col.axis = "black", at = seq(3, 10, 1), cex.axis = 1.5) # Y-axis
x <- c(1:6)
labs <- c("Ger/Spr", "Ger/Aut", 
          "Swe/Spr", "Swe/Aut", 
          "Fin/Spr", "Fin/Aut")
# text(cex = 1.5, x = x, y = -10, labs, xpd = TRUE, srt = 45, font = 1)

### Corrected functional richness
#boxplot(null$FRic.SES[1:8], null$FRic.SES[9:16], null$FRic.SES[17:68], null$FRic.SES[69:120], null$FRic.SES[121:130], null$FRic.SES[131:140]) # test to check axes
(FRic.SES_plot <- boxplot(null$FRic.SES[1:8], null$FRic.SES[9:16], null$FRic.SES[17:68], null$FRic.SES[69:120], null$FRic.SES[121:130], null$FRic.SES[131:140], 
                          yaxt = "n", xaxt = "n",  # should the axes be displayed
                          ylim = c(-2.7, 0.3), # limits of the y-axis
                          ylab = "FRic S.E.S.", xlab = "", # names of the axis labels
                          main = "E",
                          cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, # size of labels and axes
                          horizontal = FALSE, #this argument reverses the axis orientation
                          varwidth = FALSE, # varwidth makes the boxplot widths proportional to the sample size
                          notch = FALSE, # this makes a notch in the boxplot. In the notched boxplot, if two boxes' notches do not overlap this is 'strong evidence' their medians differ (Chambers et al., 1983, p. 62).
                          border = TRUE,
                          outline = FALSE,
                          cex = 1, boxlty = 1, col = c("#638a634d", "#00688c4d", "#638a634d", "#00688c4d", "#638a634d", "#00688c4d")))
axis(1, col = "black", col.axis = "black", at = seq(1, 6, 1), cex.axis = 1.5, labels = FALSE)
axis(2, col = "black", col.axis = "black", at = seq(-2.7, 0.3, 0.5), cex.axis = 1.5)
x <- c(1:6)
labs <- c("Ger/Spr", "Ger/Aut", 
          "Swe/Spr", "Swe/Aut", 
          "Fin/Spr", "Fin/Aut")
# text(cex = 1.5, x = x, y = -3.2, labs, xpd = TRUE, srt = 45, font = 1)

### Corrected functional evenness
# boxplot(null$FEve.SES[1:8], null$FEve.SES[9:16], null$FEve.SES[17:68], null$FEve.SES[69:120], null$FEve.SES[121:130], null$FEve.SES[131:140]) # test to check axes
(FEve.SES_plot <- boxplot(null$FEve.SES[1:8], null$FEve.SES[9:16], null$FEve.SES[17:68], null$FEve.SES[69:120], null$FEve.SES[121:130], null$FEve.SES[131:140], 
                          yaxt = "n", xaxt = "n",  # should the axes be displayed
                          ylim = c(-2.0, 3.0), # limits of the y-axis
                          ylab = "FEve S.E.S.", xlab = "", # names of the axis labels
                          main = "F",
                          cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, # size of labels and axes
                          horizontal = FALSE, #this argument reverses the axis orientation
                          varwidth = FALSE, # varwidth makes the boxplot widths proportional to the sample size
                          notch = FALSE, # this makes a notch in the boxplot. In the notched boxplot, if two boxes' notches do not overlap this is 'strong evidence' their medians differ (Chambers et al., 1983, p. 62).
                          border = TRUE,
                          outline = FALSE,
                          cex = 1, boxlty = 1, col = c("#638a634d", "#00688c4d", "#638a634d", "#00688c4d", "#638a634d", "#00688c4d")))
axis(1, col = "black", col.axis = "black", at = seq(1, 6, 1), cex.axis = 1.5, labels = FALSE)
axis(2, col = "black", col.axis = "black", at = seq(-2.0, 3.0, 1), cex.axis = 1.5)
x <- c(1:6)
labs <- c("Ger/Spr", "Ger/Aut", 
          "Swe/Spr", "Swe/Aut", 
          "Fin/Spr", "Fin/Aut")
# text(cex = 1.5, x = x, y = -2.8, labs, xpd = TRUE, srt = 45, font = 1)

### Corrected functional dispersion
#boxplot(null$FDis.SES[1:8], null$FDis.SES[9:16], null$FDis.SES[17:68], null$FDis.SES[69:120], null$FDis.SES[121:130], null$FDis.SES[131:140]) # test to check axes
(FDis.SES_plot <- boxplot(null$FDis.SES[1:8], null$FDis.SES[9:16], null$FDis.SES[17:68], null$FDis.SES[69:120], null$FDis.SES[121:130], null$FDis.SES[131:140], 
                          yaxt = "n", xaxt = "n",  # should the axes be displayed
                          ylim = c(-2.0, 2.0), # limits of the y-axis
                          ylab = "FDis S.E.S.", xlab = "", # names of the axis labels
                          main = "G",
                          cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, # size of labels and axes
                          horizontal = FALSE, #this argument reverses the axis orientation
                          varwidth = FALSE, # varwidth makes the boxplot widths proportional to the sample size
                          notch = FALSE, # this makes a notch in the boxplot. In the notched boxplot, if two boxes' notches do not overlap this is 'strong evidence' their medians differ (Chambers et al., 1983, p. 62).
                          border = TRUE,
                          outline = FALSE,
                          cex = 1, boxlty = 1, col = c("#638a634d", "#00688c4d", "#638a634d", "#00688c4d", "#638a634d", "#00688c4d")))
axis(1, col = "black", col.axis = "black", at = seq(1, 6, 1), cex.axis = 1.5, labels = FALSE)
axis(2, col = "black", col.axis = "black", at = seq(-2.0, 2.0, 1), cex.axis = 1.5)
x <- c(1:6)
labs <- c("Ger/Spr", "Ger/Aut", 
          "Swe/Spr", "Swe/Aut", 
          "Fin/Spr", "Fin/Aut")
text(cex = 1.5, x = x, y = -2.75, labs, xpd = TRUE, srt = 45, font = 1)

# Functional redundancy
# boxplot(div$FRed[1:8], div$FRed[9:16], div$FRed[17:68], div$FRed[69:120], div$FRed[121:130], div$FRed[131:140]) # test to check axes
(func_redund_plot <- boxplot(div$FRed[1:8], div$FRed[9:16], div$FRed[17:68], div$FRed[69:120], div$FRed[121:130], div$FRed[131:140], 
                             yaxt = "n", xaxt = "n",  # should the axes be displayed
                             ylim = c(0.30, 0.60), # limits of the y-axis
                             ylab = "FRed", xlab = "", # names of the axis labels
                             main = "H",
                             cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, # size of labels and axes
                             horizontal = FALSE, #this argument reverses the axis orientation
                             varwidth = FALSE, # varwidth makes the boxplot widths proportional to the sample size
                             notch = FALSE, # this makes a notch in the boxplot. In the notched boxplot, if two boxes' notches do not overlap this is 'strong evidence' their medians differ (Chambers et al., 1983, p. 62).
                             border = TRUE,
                             outline = FALSE,
                             cex = 1, boxlty = 1, col = c("#638a634d", "#00688c4d", "#638a634d", "#00688c4d", "#638a634d", "#00688c4d")))
axis(1, col = "black", col.axis = "black", at = seq(1, 6, 1), cex.axis = 1.5, labels = FALSE)
axis(2, col = "black", col.axis = "black", at = seq(0.30, 0.60, 0.05), cex.axis = 1.5)
x <- c(1:6)
labs <- c("Ger/Spr", "Ger/Aut", 
          "Swe/Spr", "Swe/Aut", 
          "Fin/Spr", "Fin/Aut")
text(cex = 1.5, x = x, y = 0.240, labs, xpd = TRUE, srt = 45, font = 1)

### Tturn
# boxplot(turnover$Tturn[1:8], turnover$Tturn[9:60], turnover$Tturn[61:70]) # test to check axes
(Tturn_plot <- boxplot(turnover$Tturn[1:8], turnover$Tturn[9:60], turnover$Tturn[61:70], 
                       yaxt = "n", xaxt = "n",  # should the axes be displayed
                       ylim = c(0.0, 0.7), # limits of the y-axis
                       ylab = "Tturn", xlab = "", # names of the axis labels
                       main = "I",
                       cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, # size of labels and axes
                       horizontal = FALSE, #this argument reverses the axis orientation
                       varwidth = FALSE, # varwidth makes the boxplot widths proportional to the sample size
                       notch = FALSE, # this makes a notch in the boxplot. In the notched boxplot, if two boxes' notches do not overlap this is 'strong evidence' their medians differ (Chambers et al., 1983, p. 62).
                       border = TRUE,
                       outline = FALSE,
                       cex = 1, boxlty = 1, col = c("#638a634d", "#00688c4d", "#638a634d", "#00688c4d", "#638a634d", "#00688c4d")))
axis(1, col = "black", col.axis = "black", at = seq(1, 6, 1), cex.axis = 1.5, labels = FALSE) # X-axis
axis(2, col = "black", col.axis = "black", at = seq(0.0, 0.7, 0.1), cex.axis = 1.5) # Y-axis
x <- c(1:3)
labs <- c("Germany",
          "Sweden",
          "Finland")
text(cex = 1.5, x = x, y = -0.15, labs, xpd = TRUE, srt = 45, font = 1)

### Fturn
# boxplot(turnover$Fturn[1:8], turnover$Fturn[9:60], turnover$Fturn[61:70]) # test to check axes
(Tturn_plot <- boxplot(turnover$Fturn[1:8], turnover$Fturn[9:60], turnover$Fturn[61:70], 
                       yaxt = "n", xaxt = "n",  # should the axes be displayed
                       ylim = c(0.0, 0.4), # limits of the y-axis
                       ylab = "Fturn", xlab = "", # names of the axis labels
                       main = "J",
                       cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, # size of labels and axes
                       horizontal = FALSE, #this argument reverses the axis orientation
                       varwidth = FALSE, # varwidth makes the boxplot widths proportional to the sample size
                       notch = FALSE, # this makes a notch in the boxplot. In the notched boxplot, if two boxes' notches do not overlap this is 'strong evidence' their medians differ (Chambers et al., 1983, p. 62).
                       border = TRUE,
                       outline = FALSE,
                       cex = 1, boxlty = 1, col = c("#638a634d", "#00688c4d", "#638a634d", "#00688c4d", "#638a634d", "#00688c4d")))
axis(1, col = "black", col.axis = "black", at = seq(1, 6, 1), cex.axis = 1.5, labels = FALSE) # X-axis
axis(2, col = "black", col.axis = "black", at = seq(0.0, 0.4, 0.05), cex.axis = 1.5) # Y-axis
x <- c(1:3)
labs <- c("Germany",
          "Sweden",
          "Finland")
text(cex = 1.5, x = x, y = -0.09, labs, xpd = TRUE, srt = 45, font = 1)

dev.off()
par(mfrow = c(1,1))

#### DESCRIPTIVE STATISTICS OF DIVERSITY METRICS -----------------
div_all <- cbind(div, null)
turnover
library(EnvStats)

### Beta diversity
## Taxonomic beta diversity
# Germany
(ger_taxa_bd_stat <- summaryFull(ger_pod_pd$D))
(ger_taxa_bd_stat_spr <- summaryFull(ger_pod_pd_spr$D))
(ger_taxa_bd_stat_aut <- summaryFull(ger_pod_pd_aut$D))
(colnames(ger_taxa_bd_describe <- rbind(t(ger_taxa_bd_stat), t(ger_taxa_bd_stat_spr), t(ger_taxa_bd_stat_aut))))
(ger_taxa_bd_describe <- as.data.frame(ger_taxa_bd_describe[, c(8:10, 3, 2, 13)]))
rownames(ger_taxa_bd_describe) <- c("Bt_ger", "Bt_ger_spr", "Bt_ger_aut")
ger_descriptives <- rbind(ger_taxa_bd_describe)
# Sweden
(swe_taxa_bd_stat <- summaryFull(swe_pod_pd$D))
(swe_taxa_bd_stat_spr <- summaryFull(swe_pod_pd_spr$D))
(swe_taxa_bd_stat_aut <- summaryFull(swe_pod_pd_aut$D))
(colnames(swe_taxa_bd_describe <- rbind(t(swe_taxa_bd_stat), t(swe_taxa_bd_stat_spr), t(swe_taxa_bd_stat_aut))))
(swe_taxa_bd_describe <- as.data.frame(swe_taxa_bd_describe[, c(8:10, 3, 2, 13)]))
rownames(swe_taxa_bd_describe) <- c("Bt_swe", "Bt_swe_spr", "Bt_swe_aut")
swe_descriptives <- rbind(swe_taxa_bd_describe)
# Finland
(fin_taxa_bd_stat <- summaryFull(fin_pod_pd$D))
(fin_taxa_bd_stat_spr <- summaryFull(fin_pod_pd_spr$D))
(fin_taxa_bd_stat_aut <- summaryFull(fin_pod_pd_aut$D))
(colnames(fin_taxa_bd_describe <- rbind(t(fin_taxa_bd_stat), t(fin_taxa_bd_stat_spr), t(fin_taxa_bd_stat_aut))))
(fin_taxa_bd_describe <- as.data.frame(fin_taxa_bd_describe[, c(8:10, 3, 2, 13)]))
rownames(fin_taxa_bd_describe) <- c("Bt_fin", "Bt_fin_spr", "Bt_fin_aut")
fin_descriptives <- rbind(fin_taxa_bd_describe)

## Functional beta diversity
# Germany
(ger_func_bd_stat <- summaryFull(ger_func_mat_full))
(ger_func_bd_stat_spr <- summaryFull(ger_func_mat_spr))
(ger_func_bd_stat_aut <- summaryFull(ger_func_mat_aut))
(colnames(ger_func_bd_describe <- rbind(t(ger_func_bd_stat), t(ger_func_bd_stat_spr), t(ger_func_bd_stat_aut))))
(ger_func_bd_describe <- as.data.frame(ger_func_bd_describe[, c(8:10, 3, 2, 13)]))
rownames(ger_func_bd_describe) <- c("Bf_ger", "Bf_ger_spr", "Bf_ger_aut")
ger_descriptives <- rbind(ger_descriptives, ger_func_bd_describe)
# Sweden
(swe_func_bd_stat <- summaryFull(swe_func_mat_full))
(swe_func_bd_stat_spr <- summaryFull(swe_func_mat_spr))
(swe_func_bd_stat_aut <- summaryFull(swe_func_mat_aut))
(colnames(swe_func_bd_describe <- rbind(t(swe_func_bd_stat), t(swe_func_bd_stat_spr), t(swe_func_bd_stat_aut))))
(swe_func_bd_describe <- as.data.frame(swe_func_bd_describe[, c(8:10, 3, 2, 13)]))
rownames(swe_func_bd_describe) <- c("Bf_swe", "Bf_swe_spr", "Bf_swe_aut")
swe_descriptives <- rbind(swe_descriptives, swe_func_bd_describe)
# Finland
(fin_func_bd_stat <- summaryFull(fin_func_mat_full))
(fin_func_bd_stat_spr <- summaryFull(fin_func_mat_spr))
(fin_func_bd_stat_aut <- summaryFull(fin_func_mat_aut))
(colnames(fin_func_bd_describe <- rbind(t(fin_func_bd_stat), t(fin_func_bd_stat_spr), t(fin_func_bd_stat_aut))))
(fin_func_bd_describe <- as.data.frame(fin_func_bd_describe[, c(8:10, 3, 2, 13)]))
rownames(fin_func_bd_describe) <- c("Bf_fin", "Bf_fin_spr", "Bf_fin_aut")
fin_descriptives <- rbind(fin_descriptives, fin_func_bd_describe)

### Alpha diversity
## Germany
(ger_stat <- summaryFull(div[1:16, ]))
(ger_stat_spr <- summaryFull(div[1:8, ]))
(ger_stat_aut <- summaryFull(div[9:16, ]))
(colnames(ger_describe <- rbind(t(ger_stat), t(ger_stat_spr), t(ger_stat_aut))))
(ger_describe <- ger_describe[, c(8:10, 3, 2, 13)])
rownames(ger_describe) <- c("A_ger", "E10_ger", "FDis_ger", "FEve_ger", "FRed_ger", "FRic_ger", "N0_ger", "N1_ger",
                            "A_ger_spr", "E10_ger_spr", "FDis_ger_spr", "FEve_ger_spr", "FRed_ger_spr", "FRic_ger_spr", "N0_ger_spr", "N1_ger_spr",
                            "A_ger_aut", "E10_ger_aut", "FDis_ger_aut", "FEve_ger_aut", "FRed_ger_aut", "FRic_ger_aut", "N0_ger_aut", "N1_ger_aut")
ger_descriptives <- rbind(ger_descriptives, ger_describe)
# Corrected functional metrics
(ger_stat1 <- summaryFull(div_all[1:16, 9:11]))
(ger_stat_spr1 <- summaryFull(div_all[1:8, 9:11]))
(ger_stat_aut1 <- summaryFull(div_all[9:16, 9:11]))
(colnames(ger_describe1 <- rbind(t(ger_stat1), t(ger_stat_spr1), t(ger_stat_aut1))))
(ger_describe1 <- as.data.frame(ger_describe1[, c(7:9, 3, 2, 12)]))
rownames(ger_describe1) <- c("FDis.SES_ger", "FEve.SES_ger", "FRic.SES_ger",
                             "FDis.SES_ger_spr", "FEve.SES_ger_spr", "FRic.SES_ger_spr",
                             "FDis.SES_ger_aut", "FEve.SES_ger_aut", "FRic.SES_ger_aut")
ger_descriptives <- rbind(ger_descriptives, ger_describe1)
# Turnover
(ger_turn_stat <- summaryFull(turnover[1:8,]))
(colnames(ger_turn_describe <- rbind(t(ger_turn_stat))))
(ger_turn_describe <- ger_turn_describe[, c(8:10, 3, 2, 13)])
rownames(ger_turn_describe) <- c("Fturn_ger", "Tturn_ger")
ger_descriptives <- rbind(ger_descriptives, ger_turn_describe)

## Sweden
(swe_stat <- summaryFull(div[17:120, ]))
(swe_stat_spr <- summaryFull(div[17:68, ]))
(swe_stat_aut <- summaryFull(div[69:120, ]))
(colnames(swe_describe <- rbind(t(swe_stat), t(swe_stat_spr), t(swe_stat_aut))))
(swe_describe <- swe_describe[, c(7:9, 3, 2, 12)])
rownames(swe_describe) <- c("A_swe", "E10_swe", "FDis_swe", "FEve_swe", "FRed_swe", "FRic_swe", "N0_swe", "N1_swe",
                            "A_swe_spr", "E10_swe_spr", "FDis_swe_spr", "FEve_swe_spr", "FRed_swe_spr", "FRic_swe_spr", "N0_swe_spr", "N1_swe_spr",
                            "A_swe_aut", "E10_swe_aut", "FDis_swe_aut", "FEve_swe_aut", "FRed_swe_aut", "FRic_swe_aut", "N0_swe_aut", "N1_swe_aut")
swe_descriptives <- rbind(swe_descriptives, swe_describe)
# Corrected functional metrics
(swe_stat1 <- summaryFull(div_all[17:120, 9:11]))
(swe_stat_spr1 <- summaryFull(div_all[17:68, 9:11]))
(swe_stat_aut1 <- summaryFull(div_all[69:120, 9:11]))
(colnames(swe_describe1 <- rbind(t(swe_stat1), t(swe_stat_spr1), t(swe_stat_aut1))))
(swe_describe1 <- as.data.frame(swe_describe1[, c(7:9, 3, 2, 12)]))
rownames(swe_describe1) <- c("FDis.SES_swe", "FEve.SES_swe", "FRic.SES_swe",
                             "FDis.SES_swe_spr", "FEve.SES_swe_spr", "FRic.SES_swe_spr",
                             "FDis.SES_swe_aut", "FEve.SES_swe_aut", "FRic.SES_swe_aut")
swe_descriptives <- rbind(swe_descriptives, swe_describe1)
# Turnover
(swe_turn_stat <- summaryFull(turnover[9:60,]))
(colnames(swe_turn_describe <- rbind(t(swe_turn_stat))))
(swe_turn_describe <- swe_turn_describe[, c(8:10, 3, 2, 13)])
rownames(swe_turn_describe) <- c("Fturn_swe", "Tturn_swe")
swe_descriptives <- rbind(swe_descriptives, swe_turn_describe)

## Finland
(fin_stat <- summaryFull(div[121:140, ]))
(fin_stat_spr <- summaryFull(div[121:130, ]))
(fin_stat_aut <- summaryFull(div[131:140, ]))
(colnames(fin_describe <- rbind(t(fin_stat), t(fin_stat_spr), t(fin_stat_aut))))
(fin_describe <- fin_describe[, c(8:10, 3, 2, 13)])
rownames(fin_describe) <- c("A_fin", "E10_fin", "FDis_fin", "FEve_fin", "FRed_fin", "FRic_fin", "N0_fin", "N1_fin",
                            "A_fin_spr", "E10_fin_spr", "FDis_fin_spr", "FEve_fin_spr", "FRed_fin_spr", "FRic_fin_spr", "N0_fin_spr", "N1_fin_spr",
                            "A_fin_aut", "E10_fin_aut", "FDis_fin_aut", "FEve_fin_aut", "FRed_fin_aut", "FRic_fin_aut", "N0_fin_aut", "N1_fin_aut")
fin_descriptives <- rbind(fin_descriptives, fin_describe)
# Corrected functional metrics
(fin_stat1 <- summaryFull(div_all[121:140, 9:11]))
(fin_stat_spr1 <- summaryFull(div_all[121:130, 9:11]))
(fin_stat_aut1 <- summaryFull(div_all[131:140, 9:11]))
(colnames(fin_describe1 <- rbind(t(fin_stat1), t(fin_stat_spr1), t(fin_stat_aut1))))
(fin_describe1 <- as.data.frame(fin_describe1[, c(7:9, 3, 2, 12)]))
rownames(fin_describe1) <- c("FDis.SES_fin", "FEve.SES_fin", "FRic.SES_fin",
                             "FDis.SES_fin_spr", "FEve.SES_fin_spr", "FRic.SES_fin_spr",
                             "FDis.SES_fin_aut", "FEve.SES_fin_aut", "FRic.SES_fin_aut")
fin_descriptives <- rbind(fin_descriptives, fin_describe1)
# Turnover
(fin_turn_stat <- summaryFull(turnover[61:70, ]))
(colnames(fin_turn_describe <- rbind(t(fin_turn_stat))))
(fin_turn_describe <- fin_turn_describe[, c(8:10, 3, 2, 13)])
rownames(fin_turn_describe) <- c("Fturn_fin", "Tturn_fin")
fin_descriptives <- rbind(fin_descriptives, fin_turn_describe)

all_desciptives <- rbind(ger_descriptives, swe_descriptives, fin_descriptives)
rownames(all_desciptives)
all_desciptives_sel <- all_desciptives[c(1, 4, 7:14, 31:33, 40:41,
                                         42, 45, 48:55, 72:74, 81:82,
                                         83, 86, 89:96, 113:115, 122:123),]
all_desciptives_sel <- all_desciptives_sel[c(3, 1, 9, 4, 10, 15, 2, 8, 6, 5, 13, 12, 11, 14, 7,
                                             18, 16, 24, 19, 25, 30, 17, 23, 21, 20, 28, 27, 26, 29, 22,
                                             33, 31, 39, 34, 40, 45, 32, 38, 36, 35, 43, 42, 41, 44, 37),]
all_desciptives_sel <- all_desciptives_sel[-c(2, 7, 17, 22, 32, 37),]

boxplot(log(div_all$A[1:8]), log(div_all$A[9:16]), log(div_all$A[17:68]), log(div_all$A[69:120]), log(div_all$A[121:130]), log(div_all$A[131:140]))
boxplot(div_all$N0[1:8], div_all$N0[9:16], div_all$N0[17:68], div_all$N0[69:120], div_all$N0[121:130], div_all$N0[131:140])
boxplot(div_all$E10[1:8], div_all$E10[9:16], div_all$E10[17:68], div_all$E10[69:120], div_all$E10[121:130], div_all$E10[131:140])
boxplot(div_all$N1[1:8], div_all$N1[9:16], div_all$N1[17:68], div_all$N1[69:120], div_all$N1[121:130], div_all$N1[131:140])
boxplot(div_all$FRic[1:8], div_all$FRic[9:16], div_all$FRic[17:68], div_all$FRic[69:120], div_all$FRic[121:130], div_all$FRic[131:140])
boxplot(div_all$FEve[1:8], div_all$FEve[9:16], div_all$FEve[17:68], div_all$FEve[69:120], div_all$FEve[121:130], div_all$FEve[131:140])
boxplot(div_all$FDis[1:8], div_all$FDis[9:16], div_all$FDis[17:68], div_all$FDis[69:120], div_all$FDis[121:130], div_all$FDis[131:140])
boxplot(div_all$FRic.SES[1:8], div_all$FRic.SES[9:16], div_all$FRic.SES[17:68], div_all$FRic.SES[69:120], div_all$FRic.SES[121:130], div_all$FRic.SES[131:140])
boxplot(div_all$FEve.SES[1:8], div_all$FEve.SES[9:16], div_all$FEve.SES[17:68], div_all$FEve.SES[69:120], div_all$FEve.SES[121:130], div_all$FEve.SES[131:140])
boxplot(div_all$FDis.SES[1:8], div_all$FDis.SES[9:16], div_all$FDis.SES[17:68], div_all$FDis.SES[69:120], div_all$FDis.SES[121:130], div_all$FDis.SES[131:140])
boxplot(div_all$FRed[1:8], div_all$FRed[9:16], div_all$FRed[17:68], div_all$FRed[69:120], div_all$FRed[121:130], div_all$FRed[131:140])

### Paired t-tests
## Germany
(stat_A_ger <- t.test(div_all$A[1:8], div_all$A[9:16], paired = T))
(stat_TRic_ger <- t.test(div_all$N0[1:8], div_all$N0[9:16], paired = T))
(stat_TEve_ger <- t.test(div_all$E10[1:8], div_all$E10[9:16], paired = T))
(stat_Shan_ger <- t.test(div_all$N1[1:8], div_all$N1[9:16], paired = T))
(stat_FRic_ger <- t.test(div_all$FRic[1:8], div_all$FRic[9:16], paired = T))
(stat_FEve_ger <- t.test(div_all$FEve[1:8], div_all$FEve[9:16], paired = T))
(stat_FDis_ger <- t.test(div_all$FDis[1:8], div_all$FDis[9:16], paired = T))
(stat_FRic.SES_ger <- t.test(div_all$FRic.SES[1:8], div_all$FRic.SES[9:16], paired = T))
(stat_FEve.SES_ger <- t.test(div_all$FEve.SES[1:8], div_all$FEve.SES[9:16], paired = T))
(stat_FDis.SES_ger <- t.test(div_all$FDis.SES[1:8], div_all$FDis.SES[9:16], paired = T))
(stat_FRed_ger <- t.test(div_all$FRed[1:8], div_all$FRed[9:16], paired = T))

ger_stat <- as.data.frame(rbind(cbind(stat_A_ger$statistic, stat_A_ger$p.value),
                                cbind(stat_TRic_ger$statistic, stat_TRic_ger$p.value),
                                cbind(stat_TEve_ger$statistic, stat_TEve_ger$p.value),
                                cbind(stat_Shan_ger$statistic, stat_Shan_ger$p.value),
                                cbind(stat_FRic_ger$statistic, stat_FRic_ger$p.value),
                                cbind(stat_FEve_ger$statistic, stat_FEve_ger$p.value),
                                cbind(stat_FDis_ger$statistic, stat_FDis_ger$p.value),
                                cbind(stat_FRic.SES_ger$statistic, stat_FRic.SES_ger$p.value),
                                cbind(stat_FEve.SES_ger$statistic, stat_FEve.SES_ger$p.value),
                                cbind(stat_FDis.SES_ger$statistic, stat_FDis.SES_ger$p.value),
                                cbind(stat_FRed_ger$statistic, stat_FRed_ger$p.value)))
colnames(ger_stat) <- c("t-stat", "p-val")
rownames(ger_stat) <- c("A_ger", "N0_ger", "E10_ger", "N1_ger", 
                        "FRic_ger", "FEve_ger", "FDis_ger",
                        "FRic.SES_ger", "FEve.SES_ger", "FDis.SES_ger",
                        "FRed_ger")
ger_stat

## Sweden
(stat_A_swe <- t.test(div_all$A[17:68], div_all$A[69:120], paired = T))
(stat_TRic_swe <- t.test(div_all$N0[17:68], div_all$N0[69:120], paired = T))
(stat_TEve_swe <- t.test(div_all$E10[17:68], div_all$E10[69:120], paired = T))
(stat_Shan_swe <- t.test(div_all$N1[17:68], div_all$N1[69:120], paired = T))
(stat_FRic_swe <- t.test(div_all$FRic[17:68], div_all$FRic[69:120], paired = T))
(stat_FEve_swe <- t.test(div_all$FEve[17:68], div_all$FEve[69:120], paired = T))
(stat_FDis_swe <- t.test(div_all$FDis[17:68], div_all$FDis[69:120], paired = T))
(stat_FRic.SES_swe <- t.test(div_all$FRic.SES[17:68], div_all$FRic.SES[69:120], paired = T))
(stat_FEve.SES_swe <- t.test(div_all$FEve.SES[17:68], div_all$FEve.SES[69:120], paired = T))
(stat_FDis.SES_swe <- t.test(div_all$FDis.SES[17:68], div_all$FDis.SES[69:120], paired = T))
(stat_FRed_swe <- t.test(div_all$FRed[17:68], div_all$FRed[69:120], paired = T))

swe_stat <- as.data.frame(rbind(cbind(stat_A_swe$statistic, stat_A_swe$p.value),
                                cbind(stat_TRic_swe$statistic, stat_TRic_swe$p.value),
                                cbind(stat_TEve_swe$statistic, stat_TEve_swe$p.value),
                                cbind(stat_Shan_swe$statistic, stat_Shan_swe$p.value),
                                cbind(stat_FRic_swe$statistic, stat_FRic_swe$p.value),
                                cbind(stat_FEve_swe$statistic, stat_FEve_swe$p.value),
                                cbind(stat_FDis_swe$statistic, stat_FDis_swe$p.value),
                                cbind(stat_FRic.SES_swe$statistic, stat_FRic.SES_swe$p.value),
                                cbind(stat_FEve.SES_swe$statistic, stat_FEve.SES_swe$p.value),
                                cbind(stat_FDis.SES_swe$statistic, stat_FDis.SES_swe$p.value),
                                cbind(stat_FRed_swe$statistic, stat_FRed_swe$p.value)))
colnames(swe_stat) <- c("t-stat", "p-val")
rownames(swe_stat) <- c("A_swe", "N0_swe", "E10_swe", "N1_swe", 
                        "FRic_swe", "FEve_swe", "FDis_swe",
                        "FRic.SES_swe", "FEve.SES_swe", "FDis.SES_swe",
                        "FRed_swe")
swe_stat

## Finland
(stat_A_fin <- t.test(div_all$A[121:130], div_all$A[131:140], paired = T))
(stat_TRic_fin <- t.test(div_all$N0[121:130], div_all$N0[131:140], paired = T))
(stat_TEve_fin <- t.test(div_all$E10[121:130], div_all$E10[131:140], paired = T))
(stat_Shan_fin <- t.test(div_all$N1[121:130], div_all$N1[131:140], paired = T))
(stat_FRic_fin <- t.test(div_all$FRic[121:130], div_all$FRic[131:140], paired = T))
(stat_FEve_fin <- t.test(div_all$FEve[121:130], div_all$FEve[131:140], paired = T))
(stat_FDis_fin <- t.test(div_all$FDis[121:130], div_all$FDis[131:140], paired = T))
(stat_FRic.SES_fin <- t.test(div_all$FRic.SES[121:130], div_all$FRic.SES[131:140], paired = T))
(stat_FEve.SES_fin <- t.test(div_all$FEve.SES[121:130], div_all$FEve.SES[131:140], paired = T))
(stat_FDis.SES_fin <- t.test(div_all$FDis.SES[121:130], div_all$FDis.SES[131:140], paired = T))
(stat_FRed_fin <- t.test(div_all$FRed[121:130], div_all$FRed[131:140], paired = T))

fin_stat <- as.data.frame(rbind(cbind(stat_A_fin$statistic, stat_A_fin$p.value),
                                cbind(stat_TRic_fin$statistic, stat_TRic_fin$p.value),
                                cbind(stat_TEve_fin$statistic, stat_TEve_fin$p.value),
                                cbind(stat_Shan_fin$statistic, stat_Shan_fin$p.value),
                                cbind(stat_FRic_fin$statistic, stat_FRic_fin$p.value),
                                cbind(stat_FEve_fin$statistic, stat_FEve_fin$p.value),
                                cbind(stat_FDis_fin$statistic, stat_FDis_fin$p.value),
                                cbind(stat_FRic.SES_fin$statistic, stat_FRic.SES_fin$p.value),
                                cbind(stat_FEve.SES_fin$statistic, stat_FEve.SES_fin$p.value),
                                cbind(stat_FDis.SES_fin$statistic, stat_FDis.SES_fin$p.value),
                                cbind(stat_FRed_fin$statistic, stat_FRed_fin$p.value)))
colnames(fin_stat) <- c("t-stat", "p-val")
rownames(fin_stat) <- c("A_fin", "N0_fin", "E10_fin", "N1_fin", 
                        "FRic_fin", "FEve_fin", "FDis_fin",
                        "FRic.SES_fin", "FEve.SES_fin", "FDis.SES_fin",
                        "FRed_fin")
fin_stat

all_stats <- rbind(ger_stat, swe_stat, fin_stat)

##### DB-RDAs OF COMMUNITY COMPOSITION --------------------
library(labdsv)
#### GERMANY --------------------
### TAXONOMIC COMPOSITION --------------------
## Organize data
# Spatial data
ger_geo <- geo[1:16,]
ger_dist <- geodist(ger_geo, measure = "geodesic") # A matrix representing physical distances between sites
ger_dist <- as.dist(ger_dist) # distances are in metres (m)
ger_dist_pco <- as.matrix(decostand(pco(ger_dist)$points, "standardize"))  # pco of geographic distances

# Environmental data
ger_env <- env[c(1:16), c(1:3, 6, 9, 14:20, 23:27, 41, 46:48, 50:52)] # extract german env variables
ger_env_sel <- decostand(ger_env[, c(1, 2, 5, 15:17)], "standardize") # select variables
ger_env_sel$Season <- group[1:16, ]$Season
rownames(ger_env_sel) = rownames(div[1:16, ])
ger_env_sel

# Check env data for collinearity
pairs(ger_env_sel[ , c(1:6)], lower.panel = panel.smooth, upper.panel = panel.cor, 
      diag.panel = panel.hist, main = "Pearson Correlation Matrix")
# Check for big correlations
ger_Corrs <- cor(ger_env_sel[, c(1:6)])
ger_BigCorrs <- which(ger_Corrs > 0.7 & ger_Corrs < 1, arr.ind = TRUE)
# pairs(ger_env_sel[ , unique(rownames(ger_BigCorrs))], lower.panel = panel.smooth, upper.panel = panel.cor, diag.panel = panel.hist, main = "Pearson Correlation Matrix")

# Community data (response variables)
# Remove rare species
# Look at the data, eliminate species with less than 2 total abundance
(ger_num_occur_species <- apply(comm[1:16, ] > 0,2, sum)) 
(ger_sorted_abundance <- sort(ger_num_occur_species))
# keep only spp that occur in more than 2 sites
(ger_comm <- comm[c(1:16), ger_num_occur_species >= 1])
# Compute the dissimilarity response matrix with vegan's vegdist()
(ger_comm_bray <- vegdist(ger_comm, "bray")) # Percentage difference dissimilarity matrix

## PCOA of raw community data (unconstrained)
# Run a PCoA - distance based ordination with no environmental constraints (i.e. unconstrained analysis)
is.euclid(sqrt(ger_comm_bray))
ger_bray.pcoa <- cmdscale(sqrt(ger_comm_bray), k = (nrow(ger_comm) - 1), eig = TRUE)
# Plot of the sites
ordiplot(scores(ger_bray.pcoa, choices = c(1, 2)),
         type = "t",
         main = "PCoA with species weighted averages")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
# Add weighted average projection of species
ger_comm.wa <- wascores(ger_bray.pcoa$points[, 1:2], ger_comm)
text(ger_comm.wa, rownames(ger_comm.wa), cex = 0.7, col = "red")
# A posteriori projection of environmental variables
(ger_bray.pcoa.env <- envfit(ger_bray.pcoa, ger_env_sel))
# Plot significant variables with a user-selected colour
plot(ger_bray.pcoa.env, p.max = 0.05, col = 3)

## dbRDA (square-root dissimilarity matrix)
# Function dbrda() applied to a square-rooted dissimilarity matrix will provide the test
# with a correct type I error rate.
# dbrda() on the square-rooted dissimilarity matrix
ger_bray.env.dbrda <- dbrda(sqrt(ger_comm_bray) ~ 
                              Altitude + 
                              pH + 
                              DOC + 
                              SO4 + 
                              NH4.N + 
                              NO3 + 
                              Season + 
                              Condition(ger_dist_pco),
                            data = as.data.frame(ger_env_sel),
                            add = FALSE)
# Model output
summary(ger_bray.env.dbrda, axes = FALSE)
# Global test of the RDA result
anova(ger_bray.env.dbrda, permutations = how(nperm = 999))
# Tests of all canonical axes
anova(ger_bray.env.dbrda, by = "axis", permutations = how(nperm = 999))
# Tests of all explanatory terms
anova(ger_bray.env.dbrda, by = "terms", permutations = how(nperm = 999))

## Capscale (raw species data)
# capscale() with raw (site by species) data
ger_bray.env.cap <- capscale(ger_comm ~ 
                               Altitude + 
                               pH + 
                               DOC + 
                               SO4 + 
                               NH4.N + 
                               NO3 + 
                               Season + 
                               Condition(ger_dist_pco),
                             data = as.data.frame(ger_env_sel),
                             distance = "bray",
                             add = "lingoes",
                             comm = ger_comm)
# Model output
summary(ger_bray.env.cap, axes = FALSE)
# Global test of the RDA result
anova(ger_bray.env.cap, permutations = how(nperm = 999))
# Tests of all canonical axes
anova(ger_bray.env.cap, by = "axis", permutations = how(nperm = 999))
# Tests of all explanatory terms
anova(ger_bray.env.cap, by = "terms", permutations = how(nperm = 999))
# Canonical coefficients from the rda object
coef(ger_bray.env.cap)
# Unadjusted R^2 retrieved from the rda object
(ger_bray.env.cap_R2 <- RsquareAdj(ger_bray.env.cap)$r.squared)
# Adjusted R^2 retrieved from the rda object
(ger_bray.env.cap_R2adj <- RsquareAdj(ger_bray.env.cap)$adj.r.squared)

## Plot with "wa" scores to see dispersion of sites around the factor levels
triplot.rda(ger_bray.env.cap, site.sc = "wa", scaling = 1)
triplot.rda(ger_bray.env.cap, site.sc = "wa", scaling = 2)

## Look at the variance inflation factors, values > 10 indicate there is redundancy in the model.
vif.cca(ger_bray.env.cap)

# Remove variables with high variance inflation
ger_env_sel <- ger_env_sel[, -6] # Removed: NO3

## Rerun dbRDA & Capscale() with removed variables
# dbrda() on the square-rooted dissimilarity matrix
ger_bray.env.dbrda.base <- dbrda(sqrt(ger_comm_bray) ~ 
                                   Altitude + 
                                   pH + 
                                   DOC + 
                                   SO4 + 
                                   NH4.N + 
                                   Season + 
                                   Condition(ger_dist_pco),
                                 data = as.data.frame(ger_env_sel),
                                 add = FALSE)
# Model output
summary(ger_bray.env.dbrda.base, axes = FALSE)
# Global test of the RDA result
anova(ger_bray.env.dbrda.base, permutations = how(nperm = 999))
# Tests of all canonical axes
anova(ger_bray.env.dbrda.base, by = "axis", permutations = how(nperm = 999))
# Tests of all explanatory terms
anova(ger_bray.env.dbrda.base, by = "terms", permutations = how(nperm = 999))

# capscale() with raw (site by species) data
ger_bray.env.cap.base <- capscale(ger_comm ~ 
                                    Altitude + 
                                    pH + 
                                    DOC + 
                                    SO4 + 
                                    NH4.N + 
                                    Season + 
                                    Condition(ger_dist_pco),
                                  data = as.data.frame(ger_env_sel),
                                  distance = "bray",
                                  add = "lingoes",
                                  comm = ger_comm)
# Model output
summary(ger_bray.env.cap.base, axes = FALSE)
# Global test of the RDA result
anova(ger_bray.env.cap.base, permutations = how(nperm = 999))
# Tests of all canonical axes
anova(ger_bray.env.cap.base, by = "axis", permutations = how(nperm = 999))
# Tests of all explanatory terms
anova(ger_bray.env.cap.base, by = "terms", permutations = how(nperm = 999))
# Look at the variance inflation factors, values > 10 indicate there is redundancy in the model.
vif.cca(ger_bray.env.cap.base)
# Canonical coefficients from the rda object
coef(ger_bray.env.cap.base)
# Unadjusted R^2 retrieved from the rda object
(ger_bray.env.cap.base_R2 <- RsquareAdj(ger_bray.env.cap.base)$r.squared)
# Adjusted R^2 retrieved from the rda object
(ger_bray.env.cap.base_R2adj <- RsquareAdj(ger_bray.env.cap.base)$adj.r.squared)
# Plot with "wa" scores to see dispersion of sites around the factor levels
triplot.rda(ger_bray.env.cap.base, site.sc = "wa", scaling = 1)
triplot.rda(ger_bray.env.cap.base, site.sc = "wa", scaling = 2)

## Variable selection
## Forward selection using vegan's ordiR2step()
# Ordination with all explanatory variables (full model)
ger_cap_par_full <- capscale(ger_comm ~ . + Condition(ger_dist_pco),
                             data = as.data.frame(ger_env_sel),
                             distance = "bray",
                             add = "lingoes",
                             comm = ger_comm)
(ger_cap_par_full_R2a <- RsquareAdj(ger_cap_par_full)$adj.r.squared)
# Ordination with no explanatory variables except conditional variables (null model)
ger_cap_par_null <- capscale(ger_comm ~ 1 + Condition(ger_dist_pco),
                             data = as.data.frame(ger_env_sel),
                             distance = "bray",
                             add = "lingoes",
                             comm = ger_comm)
# Forward selection
(ger_cap_par.step2.forward <- ordiR2step(ger_cap_par_null,
                                         scope = formula(ger_cap_par_full),
                                         R2scope = ger_cap_par_full_R2a,
                                         direction = "forward",
                                         permutations = how(nperm = 999)))
summary(ger_cap_par.step2.forward, axes = FALSE)
RsquareAdj(ger_cap_par.step2.forward) # Condition(ger_dist_pco) + Season + pH

## Check of explanatory variables against full model
vif.cca(ger_bray.env.cap.base)
vif.cca(ger_cap_par.step2.forward)

## ANOVA: test top models
ger_cap_par_red <- capscale(ger_comm ~ 
                              Season + 
#                             pH +
                              Condition(ger_dist_pco),
                            data = as.data.frame(ger_env_sel),
                            distance = "bray",
                            add = "lingoes",
                            comm = ger_comm)
# Model output
summary(ger_cap_par_red, axes = FALSE)
# Global test of the RDA result
anova(ger_cap_par_red, permutations = how(nperm = 999))
# Tests of all canonical axes
anova(ger_cap_par_red,permutations = how(nperm = 999), by = "axis")
# Tests of all explanatory terms
anova(ger_cap_par_red, permut = 999, by = "terms")
# Unadjusted R^2 retrieved from the rda object
(ger_cap_par_red_R2 <- RsquareAdj(ger_cap_par_red)$r.squared)
# Adjusted R^2 retrieved from the rda object
(ger_cap_par_red_R2a <- RsquareAdj(ger_cap_par_red)$adj.r.squared)
# Plot with "wa" scores to see dispersion of sites around the factor levels
# triplot.rda(ger_cap_par_red, site.sc = "wa", scaling = 1)
# triplot.rda(ger_cap_par_red, site.sc = "wa", scaling = 2)

## Parsimonious dbRDA
## dbRDA (provides test)
ger_dbrda_pars <- dbrda(sqrt(ger_comm_bray) ~ 
                          Season + 
                         pH +
                          Condition(ger_dist_pco),
                        data = as.data.frame(ger_env_sel),
                        add = FALSE)
summary(ger_dbrda_pars, axes = FALSE)
set.seed(1)
anova(ger_dbrda_pars, permutations = how(nperm = 999))
set.seed(1)
anova(ger_dbrda_pars, permutations = how(nperm = 999), by = "axis")
set.seed(1)
anova(ger_dbrda_pars, permutations = how(nperm = 999), by = "terms")
(ger_dbrda_pars_R2 <- RsquareAdj(ger_dbrda_pars)$r.squared)
(ger_dbrda_pars_R2a <- RsquareAdj(ger_dbrda_pars)$adj.r.squared)

## capscale (provides material for plotting)
ger_cap_pars <- capscale(ger_comm ~ 
                           Season + 
                           pH +
                           Condition(ger_dist_pco), 
                         data = as.data.frame(ger_env_sel),
                         distance = "bray",
                         add = "lingoes",
                         comm = ger_comm)
summary(ger_cap_pars, axes = FALSE)
anova(ger_cap_pars, permutations = how(nperm = 999))
anova(ger_cap_pars, permutations = how(nperm = 999), by = "axis")
anova(ger_cap_pars, permutations = how(nperm = 999), by = "terms")
(ger_cap_pars_R2 <- RsquareAdj(ger_cap_pars)$r.squared)
(ger_cap_pars_R2a <- RsquareAdj(ger_cap_pars)$adj.r.squared)

op <- par(mfrow = c(1, 2))
# Scaling 1
triplot.rda(ger_cap_pars,
            site.sc = "wa",
            scaling = 1,
            cex.char2 = 0.8,
            pos.env = 2,
            mult.spe = 0.9,
            mult.arrow = 0.92,
            mar.percent = 0.01)
# Scaling 2
triplot.rda(ger_cap_pars,
            site.sc = "wa",
            scaling = 2,
            cex.char2 = 0.8,
            pos.env = 1,
            pos.spe = 2,
            mult.spe = 0.9,
            mult.arrow = 0.92,
            mar.percent = 0.01)
par(op)

op <- par(mfrow = c(1, 2))

#Scaling 1
ger_fig1 <- ordiplot(ger_cap_pars, type = "none", scaling = 1, main = "RDA triplot - Scaling 1 - wa", 
                     xlim = c(-1.5, 2.5), ylim = c(-5.0, 1.0), cex.main = 1, xlab = "RDA 1", ylab = "RDA 2")
ordihull(ger_cap_pars, ger_env_sel$Season, draw = "polygon", lty = 1, col = c("#01698e64", "#658c6564"), scaling = 1, label = F)
# point symbols according to a specific factor
n <- points(ger_fig1, "sites", type = "n")
scores <- n$sites
s1 <-scores[ ,1]
s2 <-scores[ ,2]
listy1 = ger_env_sel$Season == "Spring" # make a list of T or F, for you this will be e.g. "spring" instead of "1"
s1_1 = s1[listy1] # subset by this list
s2_1 = s2[listy1] # subset by this list
listy2 = ger_env_sel$Season == "Autumn" # e.g. Autumn instead of 2
s1_2 = s1[listy2] # subset by this list
s2_2 = s2[listy2] # subset by this list
points(x = s1_1, y = s2_1, display = "wa", pch = 19, col = "black", cex = 1.0, scaling = 1, pos = 4)
points(x = s1_2, y = s2_2, display = "wa", pch = 17, col = "black", cex = 1.0, scaling = 1, pos = 4)
text(ger_fig1, "species", cex = 0.8, col = "red", pos = 4, scaling = 1)
spe.sc <- scores(ger_cap_pars, choices = 1:2, display = "sp", scaling = 1)
arrows(0,0, spe.sc[,1], spe.sc[,2], length = 0.05, angle = 30, lty = 1, col = 'red', scaling = 1)
text(ger_cap_pars, display = "wa", pos = 2, cex = 0.7, scaling = 1, col = "black")
text(ger_cap_pars, dis = "cn", col = "blue", cex = 1.0)
legend("bottomright", legend = c("Spring","Autumn"), bty = "n", cex = 1.0, pch = c(19, 17), col =  c("black", "black"))

#Scaling 2
ger_fig2 <- ordiplot(ger_cap_pars, type = "none", scaling = 2, main = "RDA triplot - Scaling 2 - wa", 
                     xlim = c(-1.5, 1.5), ylim = c(-3.0, 2.0), cex.main = 1, xlab = "RDA 1", ylab = "RDA 2")
ordihull(ger_cap_pars, ger_env_sel$Season, draw = "polygon", lty = 1, col = c("#01698e64", "#658c6564"), scaling = 2, label = F)
n <- points(ger_fig2, "sites", type = "n")
scores <- n$sites
s1 <-scores[ ,1]
s2 <-scores[ ,2]
listy1 = ger_env_sel$Season == "Spring" #make a list of T or F, for you this will be e.g. "spring" instead of "1"
s1_1 = s1[listy1] #subset by this list
s2_1 = s2[listy1] #subset by this list
listy2 = ger_env_sel$Season == "Autumn" # e.g. Fall instead of 2
s1_2 = s1[listy2] #subset by this list
s2_2 = s2[listy2] #subset by this list
points(x = s1_1, y = s2_1, display = "wa", pch = 19, col = "black", cex = 1.0, scaling = 2, pos = 4)
points(x = s1_2, y = s2_2, display = "wa", pch = 17, col = "black", cex = 1.0, scaling = 2, pos = 4)
text(ger_fig2, "species", cex = 0.8, col = "red", pos = 4, scaling = 2)
spe.sc <- scores(ger_cap_pars, choices = 1:2, display = "sp", scaling = 2)
arrows(0,0, spe.sc[,1], spe.sc[,2], length = 0.05, angle = 30, lty = 1, col = 'red', scaling = 2)
text(ger_cap_pars, display = "wa", pos = 2, cex = 0.7, scaling = 2, col = "black")
text(ger_cap_pars, dis = "cn", col = "blue", cex = 1.0)
legend("bottomright", legend = c("Spring","Autumn"), bty = "n", cex = 1.0, pch = c(19, 17), col =  c("black", "black"))
par(op)

### FUNCTIONAL COMPOSITION --------------------
## Organize data
# Spatial data
ger_dist_pco
# Environmental data
ger_env_sel

# Functional data (response variables)
# Create regional trait and abundance datasets
# funcitonal abundance data
ger_func_comm_full <- func_comm[1:16, ]
ger_func_comm_full <- ger_func_comm_full[, colSums(ger_func_comm_full) > 0] # only keeps columns with column sums > 0
# functional trait data
spp.ger.full = (row.names(comm_traits) %in% colnames(ger_func_comm_full)) # make a list of T or F is there is a match in spp
ger_traits.full = comm_traits[spp.ger.full, ] # subset by this list
# Check data
identical(row.names(ger_traits.full), colnames(ger_func_comm_full)) # check if names are the same
# Remove rare species
# Look at the data, eliminate species with less than 2 total abundance
ger_num_occur_species_func <- apply(ger_func_comm_full > 0,2, sum) 
ger_sorted_abundance_func <- sort(ger_num_occur_species_func)
# keep only spp that occur in more than 2 sites
ger_func_comm <- ger_func_comm_full[, ger_num_occur_species_func >= 1]
# create trait subset from dataset with rare taxa removed
spp.ger <- (row.names(ger_traits.full) %in% colnames(ger_func_comm))
ger_traits = ger_traits.full[spp.ger, ]
#Check if data is correct & that row names match
identical(row.names(ger_traits), colnames(ger_func_comm)) # check if names are the same
# Compute the functional dissimilarity response matrix
ger_gower <- gowdis(ger_traits) # gower disimilarity matrix
plot(ger_clust <- hclust(ger_gower, method = "ward.D"))
ger_func_beta <- beta(ger_func_comm, ger_clust, func = "Soerensen", abund = TRUE) # percentage difference
(ger_func_mat <- ger_func_beta$Btotal)

## PCOA of raw community data (unconstrained)
# Run a PCoA - distance based ordination with no environmental constraints (i.e. unconstrained analysis)
is.euclid(sqrt(ger_func_mat))
ger_soren.pcoa <- cmdscale(sqrt(ger_func_mat), k = (nrow(ger_func_comm) - 1), eig = TRUE)
## Plot of the sites
ordiplot(scores(ger_soren.pcoa, choices = c(1, 2)),
         type = "t",
         main = "PCoA with species weighted averages")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
## Add weighted average projection of species
ger_func.wa <- wascores(ger_soren.pcoa$points[, 1:2], ger_func_comm)
text(ger_func.wa, rownames(ger_func.wa), cex = 0.7, col = "red")
## A posteriori projection of environmental variables
(ger_soren.pcoa.env <- envfit(ger_soren.pcoa, ger_env_sel))
## Plot significant variables with a user-selected colour
plot(ger_soren.pcoa.env, p.max = 0.05, col = 3)

## 6. run dbRDA & Capscale
## dbrda() on the square-rooted dissimilarity matrix
ger_soren.env.dbrda.base <- dbrda(sqrt(ger_func_mat) ~ 
                                    Altitude + 
                                    pH + 
                                    DOC + 
                                    SO4 + 
                                    NH4.N + 
                                    Season + 
                                    Condition(ger_dist_pco),
                                  data = as.data.frame(ger_env_sel),
                                  add = FALSE)
# Model output
summary(ger_soren.env.dbrda.base, axes = FALSE)
# Global test of the RDA result
anova(ger_soren.env.dbrda.base, permutations = how(nperm = 999))
# Tests of all canonical axes
anova(ger_soren.env.dbrda.base, by = "axis", permutations = how(nperm = 999))
# Tests of all explanatory terms
anova(ger_soren.env.dbrda.base, by = "terms", permutations = how(nperm = 999))
# Unadjusted R^2 retrieved from the rda object
(ger_soren.env.dbrda.base_R2 <- RsquareAdj(ger_soren.env.dbrda.base)$r.squared)
# Adjusted R^2 retrieved from the rda object
(ger_soren.env.dbrda.base_R2adj <- RsquareAdj(ger_soren.env.dbrda.base)$adj.r.squared)

## capscale() with raw (site by species) data
ger_soren.env.cap.base <- capscale(ger_func_mat ~ 
                                     Altitude + 
                                     pH + 
                                     DOC + 
                                     SO4 + 
                                     NH4.N + 
                                     Season + 
                                     Condition(ger_dist_pco),
                                   data = as.data.frame(ger_env_sel),
                                   add = "lingoes",
                                   comm = ger_func_comm)
# Model output
summary(ger_soren.env.cap.base, axes = FALSE)
# Global test of the RDA result
anova(ger_soren.env.cap.base, permutations = how(nperm = 999))
# Tests of all canonical axes
anova(ger_soren.env.cap.base, by = "axis", permutations = how(nperm = 999))
# Tests of all explanatory terms
anova(ger_soren.env.cap.base, by = "terms", permutations = how(nperm = 999))
# Look at the variance inflation factors, values > 10 indicate there is redundancy in the model.
vif.cca(ger_soren.env.cap.base)
# Canonical coefficients from the rda object
coef(ger_soren.env.cap.base)
# Unadjusted R^2 retrieved from the rda object
(ger_soren.env.cap.base_R2 <- RsquareAdj(ger_soren.env.cap.base)$r.squared)
# Adjusted R^2 retrieved from the rda object
(ger_soren.env.cap.base_R2adj <- RsquareAdj(ger_soren.env.cap.base)$adj.r.squared)
# Plot with "wa" scores to see dispersion of sites around the factor levels
triplot.rda(ger_soren.env.cap.base, site.sc = "wa", scaling = 1)
triplot.rda(ger_soren.env.cap.base, site.sc = "wa", scaling = 2)

## Variable selection
## Forward selection using vegan's ordiR2step()
# Ordination with all explanatory variables (full model)
ger_soren_cap_par_full <- capscale(ger_func_mat ~ . + Condition(ger_dist_pco),
                                   data = as.data.frame(ger_env_sel),
                                   add = "lingoes",
                                   comm = ger_func_comm)
(ger_soren_cap_par_full_R2a <- RsquareAdj(ger_soren_cap_par_full)$adj.r.squared)
# Ordination with no explanatory variables except conditional variables (null model)
ger_soren_cap_par_null <- capscale(ger_func_mat ~ 1 + Condition(ger_dist_pco),
                                   data = as.data.frame(ger_env_sel),
                                   add = "lingoes",
                                   comm = ger_func_comm)
# Forward selection
(ger_soren_cap_par.step2.forward <- ordiR2step(ger_soren_cap_par_null,
                                               scope = formula(ger_soren_cap_par_full),
                                               R2scope = ger_soren_cap_par_full_R2a,
                                               direction = "forward",
                                               permutations = how(nperm = 999)))
summary(ger_soren_cap_par.step2.forward, axes = FALSE)
RsquareAdj(ger_soren_cap_par.step2.forward) # Condition(ger_pco_axis1 + ger_pco_axis2) + EC

## Check of explanatory variables against full model
vif.cca(ger_soren.env.cap.base)
#vif.cca(ger_soren_cap_par.step2.forward)

## ANOVA: test top models
ger_soren_cap_par_red <- capscale(ger_func_mat ~ 
                                    1 + Condition(ger_dist_pco),
                                  data = as.data.frame(ger_env_sel),
                                  add = "lingoes",
                                  comm = ger_func_comm)
# Model output
summary(ger_soren_cap_par_red, axes = FALSE)
# Global test of the RDA result
anova(ger_soren_cap_par_red, permutations = how(nperm = 999))
# Tests of all canonical axes
anova(ger_soren_cap_par_red,permutations = how(nperm = 999), by = "axis")
# Tests of all explanatory terms
anova(ger_soren_cap_par_red, permut = 999, by = "terms")
# Unadjusted R^2 retrieved from the rda object
(ger_soren_cap_par_red_R2 <- RsquareAdj(ger_soren_cap_par_red)$r.squared)
# Adjusted R^2 retrieved from the rda object
(ger_soren_cap_par_red_R2a <- RsquareAdj(ger_soren_cap_par_red)$adj.r.squared)
# Plot with "wa" scores to see dispersion of sites around the factor levels
#triplot.rda(ger_soren_cap_par_red, site.sc = "wa", scaling = 1)
#triplot.rda(ger_soren_cap_par_red, site.sc = "wa", scaling = 2)

## Parsimonious RDA
# dbRDA (provides test)
ger_func_dbrda_pars <- dbrda(sqrt(ger_func_mat) ~ 
                               1 + Condition(ger_dist_pco),
                             data = as.data.frame(ger_env_sel),
                             add = FALSE)
summary(ger_func_dbrda_pars, axes = FALSE)
set.seed(1)
anova(ger_func_dbrda_pars, permutations = how(nperm = 999))
set.seed(1)
anova(ger_func_dbrda_pars, permutations = how(nperm = 999), by = "axis")
set.seed(1)
anova(ger_func_dbrda_pars, permutations = how(nperm = 999), by = "terms")
(ger_func_dbrda_pars_R2 <- RsquareAdj(ger_func_dbrda_pars)$r.squared)
(ger_func_dbrda_pars_R2a <- RsquareAdj(ger_func_dbrda_pars)$adj.r.squared)

## capscale (provides material for plotting)
ger_func_cap_pars <- capscale(ger_func_mat ~ 
                                1 + Condition(ger_dist_pco), 
                              data = as.data.frame(ger_env_sel),
                              add = "lingoes",
                              comm = ger_func_comm)
summary(ger_func_cap_pars, axes = FALSE)
anova(ger_func_cap_pars, permutations = how(nperm = 999))
anova(ger_func_cap_pars, permutations = how(nperm = 999), by = "axis")
anova(ger_func_cap_pars, permutations = how(nperm = 999), by = "terms")
(ger_func_cap_pars_R2 <- RsquareAdj(ger_func_cap_pars)$r.squared)
(ger_func_cap_pars_R2a <- RsquareAdj(ger_func_cap_pars)$adj.r.squared)

# Plot with "wa" scores to see dispersion of sites around the factor levels
plot(ger_func_cap_pars, site.sc = "wa", scaling = 1)
plot(ger_func_cap_pars, site.sc = "wa", scaling = 2)

op <- par(mfrow = c(1, 2))
#Scaling 1
ger_fig1 <- ordiplot(ger_func_cap_pars, type = "none", scaling = 1, main = "RDA triplot - Scaling 1 - wa", 
                     xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), cex.main = 1, xlab = "RDA 1", ylab = "RDA 2")
ordihull(ger_func_cap_pars, ger_env_sel$Season, draw = "polygon", lty = 1, col = c("#01698e64", "#658c6564"), scaling = 1, label = F)
n <- points(ger_fig1, "sites", type = "n")
scores <- n$sites
s1 <-scores[ ,1]
s2 <-scores[ ,2]
listy1 = ger_env_sel$Season == "Spring" #make a list of T or F, for you this will be e.g. "spring" instead of "1"
s1_1 = s1[listy1] #subset by this list
s2_1 = s2[listy1] #subset by this list
listy2 = ger_env_sel$Season == "Autumn" # e.g. Fall instead of 2
s1_2 = s1[listy2] #subset by this list
s2_2 = s2[listy2] #subset by this list
points(x = s1_1, y = s2_1, display = "wa", pch = 19, col = "black", cex = 1.0, scaling = 1, pos = 4)
points(x = s1_2, y = s2_2, display = "wa", pch = 17, col = "black", cex = 1.0, scaling = 1, pos = 4)
text(ger_fig1, "species", cex = 0.8, col = "red", pos = 4, scaling = 1)
spe.sc <- scores(ger_func_cap_pars, choices = 1:2, display = "sp", scaling = 1)
arrows(0,0, spe.sc[,1], spe.sc[,2], length = 0.05, angle = 30, lty = 1, col = 'red', scaling = 1)
text(ger_func_cap_pars, display = "wa", pos = 2, cex = 0.7, scaling = 1, col = "black")
#(fit <- envfit(ger_func_cap_pars, ger_env_sel[, 2], perm = 999, display = "lc", scaling = 1))
#scores(fit, "vectors", scaling = 1)
#plot(fit, p.max = 0.05, col = "blue", cex = 1.0, lab = "pH")
legend("bottomright", legend = c("Spring","Autumn"), bty = "n", cex = 1.0, pch = c(19, 17), col =  c("black", "black"))

#Scaling 2
ger_fig2 <- ordiplot(ger_func_cap_pars, type = "none", scaling = 2, main = "RDA triplot - Scaling 2 - wa", 
                     xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), cex.main = 1, xlab = "RDA 1", ylab = "RDA 2")
ordihull(ger_func_cap_pars, ger_env_sel$Season, draw = "polygon", lty = 1, col = c("#01698e64", "#658c6564"), scaling = 2, label = F)
n <- points(ger_fig2, "sites", type = "n")
scores <- n$sites
s1 <-scores[ ,1]
s2 <-scores[ ,2]
listy1 = ger_env_sel$Season == "Spring" #make a list of T or F, for you this will be e.g. "spring" instead of "1"
s1_1 = s1[listy1] #subset by this list
s2_1 = s2[listy1] #subset by this list
listy2 = ger_env_sel$Season == "Autumn" # e.g. Fall instead of 2
s1_2 = s1[listy2] #subset by this list
s2_2 = s2[listy2] #subset by this list
points(x = s1_1, y = s2_1, display = "wa", pch = 19, col = "black", cex = 1.0, scaling = 2, pos = 4)
points(x = s1_2, y = s2_2, display = "wa", pch = 17, col = "black", cex = 1.0, scaling = 2, pos = 4)
text(ger_fig2, "species", cex = 0.8, col = "red", pos = 4, scaling = 2)
spe.sc <- scores(ger_func_cap_pars, choices = 1:2, display = "sp", scaling = 2)
arrows(0,0, spe.sc[,1], spe.sc[,2], length = 0.05, angle = 30, lty = 1, col = 'red', scaling = 2)
text(ger_func_cap_pars, display = "wa", pos = 2, cex = 0.7, scaling = 2, col = "black")
#(fit <- envfit(ger_func_cap_pars, ger_env_sel[, 2], perm = 999, display = "lc", scaling = 2))
#scores(fit, "vectors", scaling = 2)
#plot(fit, p.max = 0.05, col = "blue", cex = 1.0, lab = "pH")
legend("bottomright", legend = c("Spring","Autumn"), bty = "n", cex = 1.0, pch = c(19, 17), col =  c("black", "black"))
par(op)

#### SWEDEN ------------------
### TAXONOMIC STRUCTURE --------------------
## Organise data
## Geographic / Spatial data (covariables)
swe_geo <- geo[17:120,]
swe_dist <- geodist(swe_geo, measure = "geodesic") # A matrix representing physical distances between sites
swe_dist <- as.dist(swe_dist) # distances are in metres (m)
swe_dist_pco <- as.matrix(decostand(pco(swe_dist)$points, "standardize"))  # pco of geographic distances

# Environmental data
(swe_env <- env[c(17:120), c(1:2, 7:8, 10, 12:22, 24:25, 28:31, 37:51)]) # extract sweden env variables
(swe_env_sel <- decostand(swe_env[, c(1, 2, 5, 18:21)], "standardize"))
swe_env_sel$Season <- group[17:120, ]$Season
swe_env_sel
rownames(swe_env_sel) = rownames(div[17:120, ])
# Check env data for collinearity
pairs(swe_env_sel[ , c(1:7)], lower.panel = panel.smooth, upper.panel = panel.cor, 
      diag.panel = panel.hist, main = "Pearson Correlation Matrix")
# Check for big correlations
swe_Corrs <- cor(swe_env_sel[, c(1:7)])
swe_BigCorrs <- which(swe_Corrs > 0.7 & swe_Corrs < 1, arr.ind = TRUE)
pairs(swe_env_sel[ , unique(rownames(swe_BigCorrs))], lower.panel = panel.smooth, 
      upper.panel = panel.cor, diag.panel = panel.hist, main = "Pearson Correlation Matrix")
# Remove collinear variables
swe_env_sel <- swe_env_sel[, -6] # Remove Tot.N

## Community data (response variables)
# Remove rare species
# Look at the data, eliminate species with less than 2 total abundance
(swe_num_occur_species <- apply(comm[17:120, ] > 0,2, sum)) 
(swe_sorted_abundance <- sort(swe_num_occur_species))
# keep only spp that occur in more than 2 sites
(swe_comm <- comm[c(17:120), swe_num_occur_species >= 3])
# Compute the dissimilarity response matrix with vegan's vegdist()
(swe_comm_bray <- vegdist(swe_comm, "bray")) # Percentage difference dissimilarity matrix

## PCOA of raw community data (unconstrained)
## Run a PCoA - distance based ordination with no environmental constraints (i.e. unconstrained analysis)
is.euclid(sqrt(swe_comm_bray))
swe_bray.pcoa <- cmdscale(sqrt(swe_comm_bray), k = (nrow(swe_comm) - 1), eig = TRUE)
## Plot of the sites
ordiplot(scores(swe_bray.pcoa, choices = c(1, 2)),
         type = "t",
         main = "PCoA with species weighted averages")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
## Add weighted average projection of species
swe_comm.wa <- wascores(swe_bray.pcoa$points[, 1:2], swe_comm)
text(swe_comm.wa, rownames(swe_comm.wa), cex = 0.7, col = "red")
## A posteriori projection of environmental variables
(swe_bray.pcoa.env <- envfit(swe_bray.pcoa, swe_env_sel))
## Plot significant variables with a user-selected colour
plot(swe_bray.pcoa.env, p.max = 0.05, col = 3)

## dbRDA (square root dissimilarity matrix)
# Function dbrda() applied to a square-rooted dissimilarity matrix will provide the test
# with a correct type I error rate.
swe_bray.env.dbrda <- dbrda(sqrt(swe_comm_bray) ~ 
                              Altitude + 
                              pH + 
                              TOC + 
                              SO4 + 
                              NO2.NO3_N +
                              Tot.P +
                              Season +
                              Condition(swe_dist_pco),
                            data = as.data.frame(swe_env_sel),
                            add = FALSE)
# Model output
summary(swe_bray.env.dbrda, axes = FALSE)
# Global test of the RDA result
anova(swe_bray.env.dbrda, permutations = how(nperm = 999))
# Tests of all canonical axes
anova(swe_bray.env.dbrda, by = "axis", permutations = how(nperm = 999))
# Tests of all explanatory terms
anova(swe_bray.env.dbrda, by = "terms", permutations = how(nperm = 999))

## Capscale (raw species data)
# capscale() with raw (site by species) data
swe_bray.env.cap <- capscale(swe_comm ~ 
                               Altitude + 
                               pH + 
                               TOC + 
                               SO4 + 
                               NO2.NO3_N +
                               Tot.P +
                               Season +
                               Condition(swe_dist_pco),
                             data = as.data.frame(swe_env_sel),
                             distance = "bray",
                             add = "lingoes",
                             comm = swe_comm)
# Model output
summary(swe_bray.env.cap, axes = FALSE)
# Global test of the RDA result
anova(swe_bray.env.cap, permutations = how(nperm = 999))
# Tests of all canonical axes
anova(swe_bray.env.cap, by = "axis", permutations = how(nperm = 999))
# Tests of all explanatory terms
anova(swe_bray.env.cap, by = "terms", permutations = how(nperm = 999))
# Canonical coefficients from the rda object
coef(swe_bray.env.cap)
# Unadjusted R^2 retrieved from the rda object
(swe_bray.env.cap_R2 <- RsquareAdj(swe_bray.env.cap)$r.squared)
# Adjusted R^2 retrieved from the rda object
(swe_bray.env.cap_R2adj <- RsquareAdj(swe_bray.env.cap)$adj.r.squared)

## Plot with "wa" scores to see dispersion of sites around the factor levels
triplot.rda(swe_bray.env.cap, site.sc = "wa", scaling = 1)
triplot.rda(swe_bray.env.cap, site.sc = "wa", scaling = 2)

## Look at the variance inflation factors, values > 10 indicate there is redundancy in the model.
vif.cca(swe_bray.env.cap)

## Rerun dbRDA & capscale with removed variables
# dbrda() on the square-rooted dissimilarity matrix
swe_bray.env.dbrda.base <- dbrda(sqrt(swe_comm_bray) ~ 
                                   Altitude + 
                                   pH + 
                                   TOC + 
                                   SO4 + 
                                   NO2.NO3_N +
                                   Tot.P +
                                   Season +
                                   Condition(swe_dist_pco),
                                 data = as.data.frame(swe_env_sel),
                                 add = FALSE)
# Model output
summary(swe_bray.env.dbrda.base, axes = FALSE)
# Global test of the RDA result
anova(swe_bray.env.dbrda.base, permutations = how(nperm = 999))
# Tests of all canonical axes
anova(swe_bray.env.dbrda.base, by = "axis", permutations = how(nperm = 999))
# Tests of all explanatory terms
anova(swe_bray.env.dbrda.base, by = "terms", permutations = how(nperm = 999))

## capscale() with raw (site by species) data
swe_bray.env.cap.base <- capscale(swe_comm ~ 
                                    Altitude + 
                                    pH + 
                                    TOC + 
                                    SO4 + 
                                    NO2.NO3_N +
                                    Tot.P +
                                    Season +
                                    Condition(swe_dist_pco),
                                  data = as.data.frame(swe_env_sel),
                                  distance = "bray",
                                  add = "lingoes",
                                  comm = swe_comm)
# Model output
summary(swe_bray.env.cap.base, axes = FALSE)
# Global test of the RDA result
anova(swe_bray.env.cap.base, permutations = how(nperm = 999))
# Tests of all canonical axes
anova(swe_bray.env.cap.base, by = "axis", permutations = how(nperm = 999))
# Tests of all explanatory terms
anova(swe_bray.env.cap.base, by = "terms", permutations = how(nperm = 999))
# Look at the variance inflation factors, values > 10 indicate there is redundancy in the model.
vif.cca(swe_bray.env.cap.base)
# Canonical coefficients from the rda object
coef(swe_bray.env.cap.base)
# Unadjusted R^2 retrieved from the rda object
(swe_bray.env.cap.base_R2 <- RsquareAdj(swe_bray.env.cap.base)$r.squared)
# Adjusted R^2 retrieved from the rda object
(swe_bray.env.cap.base_R2adj <- RsquareAdj(swe_bray.env.cap.base)$adj.r.squared)
# Plot with "wa" scores to see dispersion of sites around the factor levels
triplot.rda(swe_bray.env.cap.base, site.sc = "wa", scaling = 1)
triplot.rda(swe_bray.env.cap.base, site.sc = "wa", scaling = 2)

## Variable selection
# Forward selection using vegan's ordiR2step()
# Ordination with all explanatory variables (full model)
swe_cap_par_full <- capscale(swe_comm ~ . + Condition(swe_dist_pco),
                             data = as.data.frame(swe_env_sel),
                             distance = "bray",
                             add = "lingoes",
                             comm = swe_comm)
(swe_cap_par_full_R2a <- RsquareAdj(swe_cap_par_full)$adj.r.squared)
# Ordination with no explanatory variables except conditional variables (null model)
swe_cap_par_null <- capscale(swe_comm ~ 1 + Condition(swe_dist_pco),
                             data = as.data.frame(swe_env_sel),
                             distance = "bray",
                             add = "lingoes",
                             comm = swe_comm)
# Forward selection
(swe_cap_par.step2.forward <- ordiR2step(swe_cap_par_null,
                                         scope = formula(swe_cap_par_full),
                                         R2scope = swe_cap_par_full_R2a,
                                         direction = "forward",
                                         permutations = how(nperm = 999)))
summary(swe_cap_par.step2.forward, axes = FALSE)
RsquareAdj(swe_cap_par.step2.forward) # Condition(swe_dist_pco) + Season + Altitude + TOC + SO4

## Check of explanatory variables against full model
vif.cca(swe_bray.env.cap.base)
vif.cca(swe_cap_par.step2.forward)

## ANOVA: test top models
swe_cap_par_red <- capscale(swe_comm ~ 
                              Season + 
                              Altitude + 
                              TOC + 
                              SO4 +
                              Condition(swe_dist_pco),
                            data = as.data.frame(swe_env_sel),
                            distance = "bray",
                            add = "lingoes",
                            comm = swe_comm)
# Model output
summary(swe_cap_par_red, axes = FALSE)
# Global test of the RDA result
anova(swe_cap_par_red, permutations = how(nperm = 999))
# Tests of all canonical axes
anova(swe_cap_par_red,permutations = how(nperm = 999), by = "axis")
# Tests of all explanatory terms
anova(swe_cap_par_red, permut = 999, by = "terms")
# Unadjusted R^2 retrieved from the rda object
(swe_cap_par_red_R2 <- RsquareAdj(swe_cap_par_red)$r.squared)
# Adjusted R^2 retrieved from the rda object
(swe_cap_par_red_R2a <- RsquareAdj(swe_cap_par_red)$adj.r.squared)
# Plot with "wa" scores to see dispersion of sites around the factor levels
triplot.rda(swe_cap_par_red, site.sc = "wa", scaling = 1)
triplot.rda(swe_cap_par_red, site.sc = "wa", scaling = 2)

## Parsimonious rda
# dbRDA (provides test)
swe_dbrda_pars <- dbrda(sqrt(swe_comm_bray) ~ 
                          Season + 
                          Altitude + 
                          TOC + 
                          SO4 +
                          Condition(swe_dist_pco),
                        data = as.data.frame(swe_env_sel),
                        add = FALSE)
summary(swe_dbrda_pars, axes = FALSE)
set.seed(1)
anova(swe_dbrda_pars, permutations = how(nperm = 999))
set.seed(1)
anova(swe_dbrda_pars, permutations = how(nperm = 999), by = "axis")
set.seed(1)
anova(swe_dbrda_pars, permutations = how(nperm = 999), by = "terms")
(swe_dbrda_pars_R2 <- RsquareAdj(swe_dbrda_pars)$r.squared)
(swe_dbrda_pars_R2a <- RsquareAdj(swe_dbrda_pars)$adj.r.squared)

## capscale (provides material for plotting)
swe_cap_pars <- capscale(swe_comm ~ 
                           Season + 
                           Altitude + 
                           TOC + 
                           SO4 +
                           Condition(swe_dist_pco),
                         data = as.data.frame(swe_env_sel),
                         distance = "bray",
                         add = "lingoes",
                         comm = swe_comm)
summary(swe_cap_pars, axes = FALSE)
anova(swe_cap_pars, permutations = how(nperm = 999))
anova(swe_cap_pars, permutations = how(nperm = 999), by = "axis")
anova(swe_cap_pars, permutations = how(nperm = 999), by = "terms")
(swe_cap_pars_R2 <- RsquareAdj(swe_cap_pars)$r.squared)
(swe_cap_pars_R2a <- RsquareAdj(swe_cap_pars)$adj.r.squared)

# Plot with "wa" scores to see dispersion of sites around the factor levels
triplot.rda(swe_cap_pars, site.sc = "wa", scaling = 1)
triplot.rda(swe_cap_pars, site.sc = "wa", scaling = 2)

op <- par(mfrow = c(1, 2))
# Scaling 1
triplot.rda(swe_cap_pars,
            site.sc = "wa",
            scaling = 1,
            cex.char2 = 0.8,
            pos.env = 2,
            mult.spe = 0.9,
            mult.arrow = 0.92,
            mar.percent = 0.01)
# Scaling 2
triplot.rda(swe_cap_pars,
            site.sc = "wa",
            scaling = 2,
            cex.char2 = 0.8,
            pos.env = 2,
            mult.spe = 0.9,
            mult.arrow = 0.92,
            mar.percent = 0.01)
par(op)

op <- par(mfrow = c(1, 2))
# Scaling 1
swe_fig1 <- ordiplot(swe_cap_pars, type = "none", scaling = 1, main = "RDA triplot - Scaling 1 - wa", 
                     xlim = c(-5.0, 8.0), ylim = c(-8.0, 2.0), cex.main = 1, xlab = "RDA 1", ylab = "RDA 2")
ordihull(swe_cap_pars, swe_env_sel$Season, draw = "polygon", lty = 1, col = c("#01698e64", "#658c6564"), scaling = 1, label = F)
# point symbols according to a specific factor
n <- points(swe_fig1, "sites", type = "n")
scores <- n$sites
s1 <-scores[ ,1]
s2 <-scores[ ,2]
listy1 = swe_env_sel$Season == "Spring" # make a list of T or F, for you this will be e.g. "spring" instead of "1"
s1_1 = s1[listy1] # subset by this list
s2_1 = s2[listy1] # subset by this list
listy2 = swe_env_sel$Season == "Autumn" # e.g. Autumn instead of 2
s1_2 = s1[listy2] # subset by this list
s2_2 = s2[listy2] # subset by this list
points(x = s1_1, y = s2_1, display = "wa", pch = 19, col = "black", cex = 1.0, scaling = 1, pos = 4)
points(x = s1_2, y = s2_2, display = "wa", pch = 17, col = "black", cex = 1.0, scaling = 1, pos = 4)
text(swe_fig1, "species", cex = 0.8, col = "red", pos = 4, scaling = 1)
spe.sc <- scores(swe_cap_pars, choices = 1:2, display = "sp", scaling = 1)
arrows(0,0, spe.sc[,1], spe.sc[,2], length = 0.05, angle = 30, lty = 1, col = 'red', scaling = 1)
text(swe_cap_pars, display = "wa", pos = 2, cex = 0.7, scaling = 1, col = "black")
text(swe_cap_pars, dis = "cn", col = "blue", cex = 1.0)
legend("bottomright", legend = c("Spring","Autumn"), bty = "n", cex = 1.0, pch = c(19, 17), col =  c("black", "black"))

# Scaling 2
swe_fig2 <- ordiplot(swe_cap_pars, type = "none", scaling = 2, main = "RDA triplot - Scaling 2 - wa", 
                     xlim = c(-3.0, 4.0), ylim = c(-5.0, 4.0), cex.main = 1, xlab = "RDA 1", ylab = "RDA 2")
ordihull(swe_cap_pars, swe_env_sel$Season, draw = "polygon", lty = 1, col = c("#01698e64", "#658c6564"), scaling = 2, label = F)
n <- points(swe_fig2, "sites", type = "n")
scores <- n$sites
s1 <-scores[ ,1]
s2 <-scores[ ,2]
listy1 = swe_env_sel$Season == "Spring" #make a list of T or F, for you this will be e.g. "spring" instead of "1"
s1_1 = s1[listy1] #subset by this list
s2_1 = s2[listy1] #subset by this list
listy2 = swe_env_sel$Season == "Autumn" # e.g. Fall instead of 2
s1_2 = s1[listy2] #subset by this list
s2_2 = s2[listy2] #subset by this list
points(x = s1_1, y = s2_1, display = "wa", pch = 19, col = "black", cex = 1.0, scaling = 2, pos = 4)
points(x = s1_2, y = s2_2, display = "wa", pch = 17, col = "black", cex = 1.0, scaling = 2, pos = 4)
text(swe_fig2, "species", cex = 0.8, col = "red", pos = 4, scaling = 2)
spe.sc <- scores(swe_cap_pars, choices = 1:2, display = "sp", scaling = 2)
arrows(0,0, spe.sc[,1], spe.sc[,2], length = 0.05, angle = 30, lty = 1, col = 'red', scaling = 2)
text(swe_cap_pars, display = "wa", pos = 2, cex = 0.7, scaling = 2, col = "black")
text(swe_cap_pars, dis = "cn", col = "blue", cex = 1)
legend("bottomright", legend = c("Spring","Autumn"), bty = "n", cex = 1.0, pch = c(19, 17), col =  c("black", "black"))
par(op)

### FUNCTIONAL STRUCTURE --------------------
## Organise data
# Distance matrix for covariable (i.e. the variable used as a condition factor)
swe_dist_pco

# Environmental data (predictor variables)
swe_env_sel

# Functional data (response variables)
# Create regional trait and abundance datasets
# funcitonal abundance data
swe_func_comm_full <- func_comm[17:120, ]
swe_func_comm_full <- swe_func_comm_full[, colSums(swe_func_comm_full) > 0] # only keeps columns with column sums > 0
# functional trait data
spp.swe.full = (row.names(comm_traits) %in% colnames(swe_func_comm_full)) # make a list of T or F is there is a match in spp
swe_traits.full = comm_traits[spp.swe.full, ] # subset by this list
# Check data
identical(rownames(swe_traits.full), colnames(swe_func_comm_full)) # check if names are the same
# Remove rare species
# Look at the data, eliminate species with less than 2 total abundance
swe_num_occur_species_func <- apply(swe_func_comm_full > 0,2, sum) 
swe_sorted_abundance_func <- sort(swe_num_occur_species_func)
# keep only spp that occur in more than 2 sites
swe_func_comm <- swe_func_comm_full[, swe_num_occur_species_func >= 3]
# create trait subset from dataset with rare taxa removed
spp.swe <- (row.names(swe_traits.full) %in% colnames(swe_func_comm))
swe_traits = swe_traits.full[spp.swe, ]
#Check if data is correct & that row names match
identical(row.names(swe_traits), colnames(swe_func_comm)) # check if names are the same
# Compute the functional dissimilarity response matrix
swe_gower <- gowdis(swe_traits) # gower disimilarity matrix
plot(swe_clust <- hclust(swe_gower, method = "ward.D"))
swe_func_beta <- beta(swe_func_comm, swe_clust, func = "Soerensen", abund = TRUE) # percentage difference
(swe_func_mat <- swe_func_beta$Btotal)

## PCOA of raw community data (unconstrained)
# Run a PCoA - distance based ordination with no environmental constraints (i.e. unconstrained analysis)
is.euclid(sqrt(swe_func_mat))
swe_soren.pcoa <- cmdscale(sqrt(swe_func_mat), k = (nrow(swe_func_comm) - 1), eig = TRUE)
## Plot of the sites
ordiplot(scores(swe_soren.pcoa, choices = c(1, 2)),
         type = "t",
         main = "PCoA with species weighted averages")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
## Add weighted average projection of species
swe_func.wa <- wascores(swe_soren.pcoa$points[, 1:2], swe_func_comm)
text(swe_func.wa, rownames(swe_func.wa), cex = 0.7, col = "red")
## A posteriori projection of environmental variables
(swe_soren.pcoa.env <- envfit(swe_soren.pcoa, swe_env_sel))
## Plot significant variables with a user-selected colour
plot(swe_soren.pcoa.env, p.max = 0.05, col = 3)

## Run dbRDA & capscale with removed variables
## dbrda() on the square-rooted dissimilarity matrix
swe_soren.env.dbrda.base <- dbrda(sqrt(swe_func_mat) ~ 
                                    Altitude + 
                                    pH + 
                                    TOC + 
                                    SO4 + 
                                    NO2.NO3_N +
                                    Tot.P +
                                    Season +
                                    Condition(swe_dist_pco),
                                  data = as.data.frame(swe_env_sel),
                                  add = FALSE)
# Model output
summary(swe_soren.env.dbrda.base, axes = FALSE)
# Global test of the RDA result
anova(swe_soren.env.dbrda.base, permutations = how(nperm = 999))
# Tests of all canonical axes
anova(swe_soren.env.dbrda.base, by = "axis", permutations = how(nperm = 999))
# Tests of all explanatory terms
anova(swe_soren.env.dbrda.base, by = "terms", permutations = how(nperm = 999))
# Unadjusted R^2 retrieved from the rda object
(swe_soren.env.dbrda.base_R2 <- RsquareAdj(swe_soren.env.dbrda.base)$r.squared)
# Adjusted R^2 retrieved from the rda object
(swe_soren.env.dbrda.base_R2adj <- RsquareAdj(swe_soren.env.dbrda.base)$adj.r.squared)

## capscale() with raw (site by species) data
swe_soren.env.cap.base <- capscale(swe_func_mat ~ 
                                     Altitude + 
                                     pH + 
                                     TOC + 
                                     SO4 + 
                                     NO2.NO3_N +
                                     Tot.P +
                                     Season +
                                     Condition(swe_dist_pco),
                                   data = as.data.frame(swe_env_sel),
                                   add = "lingoes",
                                   comm = swe_func_comm)
# Model output
summary(swe_soren.env.cap.base, axes = FALSE)
# Global test of the RDA result
anova(swe_soren.env.cap.base, permutations = how(nperm = 999))
# Tests of all canonical axes
anova(swe_soren.env.cap.base, by = "axis", permutations = how(nperm = 999))
# Tests of all explanatory terms
anova(swe_soren.env.cap.base, by = "terms", permutations = how(nperm = 999))
# Look at the variance inflation factors, values > 10 indicate there is redundancy in the model.
vif.cca(swe_soren.env.cap.base)
# Canonical coefficients from the rda object
coef(swe_soren.env.cap.base)
# Unadjusted R^2 retrieved from the rda object
(swe_soren.env.cap.base_R2 <- RsquareAdj(swe_soren.env.cap.base)$r.squared)
# Adjusted R^2 retrieved from the rda object
(swe_soren.env.cap.base_R2adj <- RsquareAdj(swe_soren.env.cap.base)$adj.r.squared)
# Plot with "wa" scores to see dispersion of sites around the factor levels
triplot.rda(swe_soren.env.cap.base, site.sc = "wa", scaling = 1)
triplot.rda(swe_soren.env.cap.base, site.sc = "wa", scaling = 2)

## Variable selection
# Forward selection using vegan's ordiR2step()
# Ordination with all explanatory variables (full model)
swe_soren_cap_par_full <- capscale(swe_func_mat ~ . + Condition(swe_dist_pco),
                                   data = as.data.frame(swe_env_sel),
                                   add = "lingoes",
                                   comm = swe_func_comm)
(swe_soren_cap_par_full_R2a <- RsquareAdj(swe_soren_cap_par_full)$adj.r.squared)
# Ordination with no explanatory variables except conditional variables (null model)
swe_soren_cap_par_null <- capscale(swe_func_mat ~ 1 + Condition(swe_dist_pco),
                                   data = as.data.frame(swe_env_sel),
                                   add = "lingoes",
                                   comm = swe_func_comm)
# Forward selection
(swe_soren_cap_par.step2.forward <- ordiR2step(swe_soren_cap_par_null,
                                               scope = formula(swe_soren_cap_par_full),
                                               R2scope = swe_soren_cap_par_full_R2a,
                                               direction = "forward",
                                               permutations = how(nperm = 999)))
summary(swe_soren_cap_par.step2.forward, axes = FALSE)
RsquareAdj(swe_soren_cap_par.step2.forward) # Condition(swe_dist_pco) + TOC + Season + Altitude + SO4

## Check of explanatory variables against full model
vif.cca(swe_soren.env.cap.base)
vif.cca(swe_soren_cap_par.step2.forward)

## ANOVA: test top models
swe_soren_cap_par_red <- capscale(swe_func_mat ~ 
                                    TOC +
                                    Season +
                                    Altitude +
                                    SO4 +
                                    Condition(swe_dist_pco),
                                  data = as.data.frame(swe_env_sel),
                                  add = "lingoes",
                                  comm = swe_func_comm)
# Model output
summary(swe_soren_cap_par_red, axes = FALSE)
# Global test of the RDA result
anova(swe_soren_cap_par_red, permutations = how(nperm = 999))
# Tests of all canonical axes
anova(swe_soren_cap_par_red,permutations = how(nperm = 999), by = "axis")
# Tests of all explanatory terms
anova(swe_soren_cap_par_red, permut = 999, by = "terms")
# Unadjusted R^2 retrieved from the rda object
(swe_soren_cap_par_red_R2 <- RsquareAdj(swe_soren_cap_par_red)$r.squared)
# Adjusted R^2 retrieved from the rda object
(swe_soren_cap_par_red_R2a <- RsquareAdj(swe_soren_cap_par_red)$adj.r.squared)
# Plot with "wa" scores to see dispersion of sites around the factor levels
triplot.rda(swe_soren_cap_par_red, site.sc = "wa", scaling = 1)
triplot.rda(swe_soren_cap_par_red, site.sc = "wa", scaling = 2)

## Parsimonious rda
# dbRDA (provides test)
swe_func_dbrda_pars <- dbrda(sqrt(swe_func_mat) ~ 
                               TOC +
                               Season +
                               Altitude +
                               SO4 +
                               Condition(swe_dist_pco),
                             data = as.data.frame(swe_env_sel),
                             add = FALSE)
summary(swe_func_dbrda_pars, axes = FALSE)
set.seed(1)
anova(swe_func_dbrda_pars, permutations = how(nperm = 999))
set.seed(1)
anova(swe_func_dbrda_pars, permutations = how(nperm = 999), by = "axis")
set.seed(1)
anova(swe_func_dbrda_pars, permutations = how(nperm = 999), by = "terms")
(swe_func_dbrda_pars_R2 <- RsquareAdj(swe_func_dbrda_pars)$r.squared)
(swe_func_dbrda_pars_R2a <- RsquareAdj(swe_func_dbrda_pars)$adj.r.squared)

## capscale (provides material for plotting)
swe_func_cap_pars <- capscale(swe_func_mat ~ 
                                TOC +
                                Season +
                                Altitude +
                                SO4 +
                                Condition(swe_dist_pco), 
                              data = as.data.frame(swe_env_sel),
                              add = "lingoes",
                              comm = swe_func_comm)
summary(swe_func_cap_pars, axes = FALSE)
anova(swe_func_cap_pars, permutations = how(nperm = 999))
anova(swe_func_cap_pars, permutations = how(nperm = 999), by = "axis")
anova(swe_func_cap_pars, permutations = how(nperm = 999), by = "terms")
(swe_func_cap_pars_R2 <- RsquareAdj(swe_func_cap_pars)$r.squared)
(swe_func_cap_pars_R2a <- RsquareAdj(swe_func_cap_pars)$adj.r.squared)

op <- par(mfrow = c(1, 2))
# Scaling 1
triplot.rda(swe_func_cap_pars,
            site.sc = "wa",
            scaling = 1,
            cex.char2 = 0.8,
            pos.env = 2,
            mult.spe = 0.9,
            mult.arrow = 0.92,
            mar.percent = 0.01)
# Scaling 2
triplot.rda(swe_func_cap_pars,
            site.sc = "wa",
            scaling = 2,
            cex.char2 = 0.8,
            pos.env = 2,
            mult.spe = 0.9,
            mult.arrow = 0.92,
            mar.percent = 0.01)
par(op)

op <- par(mfrow = c(1, 2))
#Scaling 1
swe_fig1 <- ordiplot(swe_func_cap_pars, type = "none", scaling = 1, main = "RDA triplot - Scaling 1 - wa", 
                     xlim = c(-2.0, 9.0), ylim = c(-6.0, 8.0), cex.main = 1, xlab = "RDA 1", ylab = "RDA 2")
ordihull(swe_func_cap_pars, swe_env_sel$Season, draw = "polygon", lty = 1, col = c("#01698e64", "#658c6564"), scaling = 1, label = F)
# point symbols according to a specific factor
n <- points(swe_fig1, "sites", type = "n")
scores <- n$sites
s1 <-scores[ ,1]
s2 <-scores[ ,2]
listy1 = swe_env_sel$Season == "Spring" # make a list of T or F, for you this will be e.g. "spring" instead of "1"
s1_1 = s1[listy1] # subset by this list
s2_1 = s2[listy1] # subset by this list
listy2 = swe_env_sel$Season == "Autumn" # e.g. Autumn instead of 2
s1_2 = s1[listy2] # subset by this list
s2_2 = s2[listy2] # subset by this list
points(x = s1_1, y = s2_1, display = "wa", pch = 19, col = "black", cex = 1.0, scaling = 1, pos = 4)
points(x = s1_2, y = s2_2, display = "wa", pch = 17, col = "black", cex = 1.0, scaling = 1, pos = 4)
text(swe_fig1, "species", cex = 0.8, col = "red", pos = 4, scaling = 1)
spe.sc <- scores(swe_func_cap_pars, choices = 1:2, display = "sp", scaling = 1)
arrows(0,0, spe.sc[,1], spe.sc[,2], length = 0.05, angle = 30, lty = 1, col = 'red', scaling = 1)
text(swe_func_cap_pars, display = "wa", pos = 2, cex = 0.7, scaling = 1, col = "black")
text(swe_func_cap_pars, dis = "cn", col = "blue", cex = 1.0)
legend("bottomright", legend = c("Spring","Autumn"), bty = "n", cex = 1.0, pch = c(19, 17), col =  c("black", "black"))

#Scaling 2
swe_fig2 <- ordiplot(swe_func_cap_pars, type = "none", scaling = 2, main = "RDA triplot - Scaling 1 - wa", 
                     xlim = c(-4.0, 8.0), ylim = c(-5.0, 7.0), cex.main = 1, xlab = "RDA 1", ylab = "RDA 2")
ordihull(swe_func_cap_pars, swe_env_sel$Season, draw = "polygon", lty = 1, col = c("#01698e64", "#658c6564"), scaling = 2, label = F)
# point symbols according to a specific factor
n <- points(swe_fig2, "sites", type = "n")
scores <- n$sites
s1 <-scores[ ,1]
s2 <-scores[ ,2]
listy1 = swe_env_sel$Season == "Spring" # make a list of T or F, for you this will be e.g. "spring" instead of "1"
s1_1 = s1[listy1] # subset by this list
s2_1 = s2[listy1] # subset by this list
listy2 = swe_env_sel$Season == "Autumn" # e.g. Autumn instead of 2
s1_2 = s1[listy2] # subset by this list
s2_2 = s2[listy2] # subset by this list
points(x = s1_1, y = s2_1, display = "wa", pch = 19, col = "black", cex = 1.0, scaling = 2, pos = 4)
points(x = s1_2, y = s2_2, display = "wa", pch = 17, col = "black", cex = 1.0, scaling = 2, pos = 4)
text(swe_fig2, "species", cex = 0.8, col = "red", pos = 4, scaling = 2)
spe.sc <- scores(swe_func_cap_pars, choices = 1:2, display = "sp", scaling = 2)
arrows(0,0, spe.sc[,1], spe.sc[,2], length = 0.05, angle = 30, lty = 1, col = 'red', scaling = 2)
text(swe_func_cap_pars, display = "wa", pos = 2, cex = 0.7, scaling = 2, col = "black")
text(swe_func_cap_pars, dis = "cn", col = "blue", cex = 1.0)
legend("bottomright", legend = c("Spring","Autumn"), bty = "n", cex = 1.0, pch = c(19, 17), col =  c("black", "black"))
par(op)

#### FINLAND --------------------
# Our aim is double: we want a test with a correct type I error rate, and we want to
# plot the results with the species scores as well. Function dbrda() applied to a
# square-rooted dissimilarity matrix will provide the test, and function capscale()
# with the Lingoes correction will provide the material for the plot.

### TAXONOMIC STRUCTURE --------------------
## Organize data
# Spatial data
fin_geo <- geo[121:140,]
fin_dist <- geodist(fin_geo, measure = "geodesic") # A matrix representing physical distances between sites
fin_dist <- as.dist(fin_dist) # distances are in metres (m)
fin_dist_pco <- as.matrix(decostand(pco(fin_dist)$points, "standardize"))  # pco of geographic distances

# Environmental data
fin_env <- env[c(121:140), c(1:6, 9, 11, 29:37)] # extract finland env variables
fin_env_sel <- decostand(fin_env[, c(1, 2, 7, 9:10)], "standardize") # select variables
fin_env_sel$Season <- group[121:140, ]$Season
rownames(fin_env_sel) = rownames(div[121:140, ])
fin_env_sel
# Check env data for collinearity
pairs(fin_env_sel[ , c(1:5)], lower.panel = panel.smooth, upper.panel = panel.cor, 
      diag.panel = panel.hist, main = "Pearson Correlation Matrix")
# Check for big correlations
fin_Corrs <- cor(fin_env_sel[, c(1:5)])
fin_BigCorrs <- which(fin_Corrs > 0.7 & fin_Corrs < 1, arr.ind = TRUE)
#pairs(fin_env_sel[ , unique(rownames(fin_BigCorrs))], lower.panel = panel.smooth, upper.panel = panel.cor, diag.panel = panel.hist, main = "Pearson Correlation Matrix")

## Community data (response variables)
# Remove rare species
# Look at the data, eliminate species with less than 2 total abundance
(fin_num_occur_species <- apply(comm[121:140, ] > 0,2, sum)) 
(fin_sorted_abundance <- sort(fin_num_occur_species))
# keep only spp that occur in more than 2 sites
(fin_comm <- comm[c(121:140), fin_num_occur_species >= 1])
# Compute the dissimilarity response matrix with vegan's vegdist()
(fin_comm_bray <- vegdist(fin_comm, "bray")) # Percentage difference dissimilarity matrix

## PCOA of raw community data (unconstrained)
## Run a PCoA - distance based ordination with no environmental constraints (i.e. unconstrained analysis)
is.euclid(sqrt(fin_comm_bray))
fin_bray.pcoa <- cmdscale(sqrt(fin_comm_bray), k = (nrow(fin_comm) - 1), eig = TRUE)
## Plot of the sites
ordiplot(scores(fin_bray.pcoa, choices = c(1, 2)),
         type = "t",
         main = "PCoA with species weighted averages")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
## Add weighted average projection of species
fin_comm.wa <- wascores(fin_bray.pcoa$points[, 1:2], fin_comm)
text(fin_comm.wa, rownames(fin_comm.wa), cex = 0.7, col = "red")
## A posteriori projection of environmental variables
(fin_bray.pcoa.env <- envfit(fin_bray.pcoa, fin_env_sel))
## Plot significant variables with a user-selected colour
plot(fin_bray.pcoa.env, p.max = 0.05, col = 3)

## dbRDA (quare-rooted dissimilarity matrix)
# Function dbrda() applied to a square-rooted dissimilarity matrix will provide the test
# with a correct type I error rate.
# dbrda() on the square-rooted dissimilarity matrix
fin_bray.env.dbrda <- dbrda(sqrt(fin_comm_bray) ~ 
                              Altitude + 
                              pH + 
                              DOC + 
                              Tot.N +
                              Tot.P + 
                              Season + 
                              Condition(fin_dist_pco),
                            data = as.data.frame(fin_env_sel),
                            add = FALSE)
# Model output
summary(fin_bray.env.dbrda, axes = FALSE)
# Global test of the RDA result
anova(fin_bray.env.dbrda, permutations = how(nperm = 999))
# Tests of all canonical axes
anova(fin_bray.env.dbrda, by = "axis", permutations = how(nperm = 999))
# Tests of all explanatory terms
anova(fin_bray.env.dbrda, by = "terms", permutations = how(nperm = 999))

## Capscale  (raw species data)
## capscale() with raw (site by species) data
fin_bray.env.cap <- capscale(fin_comm ~ 
                               Altitude + 
                               pH + 
                               DOC + 
                               Tot.N +
                               Tot.P + 
                               Season + 
                               Condition(fin_dist_pco),
                             data = as.data.frame(fin_env_sel),
                             distance = "bray",
                             add = "lingoes",
                             comm = fin_comm)
# Model output
summary(fin_bray.env.cap, axes = FALSE)
# Global test of the RDA result
anova(fin_bray.env.cap, permutations = how(nperm = 999))
# Tests of all canonical axes
anova(fin_bray.env.cap, by = "axis", permutations = how(nperm = 999))
# Tests of all explanatory terms
anova(fin_bray.env.cap, by = "terms", permutations = how(nperm = 999))
# Canonical coefficients from the rda object
coef(fin_bray.env.cap)
# Unadjusted R^2 retrieved from the rda object
(fin_bray.env.cap_R2 <- RsquareAdj(fin_bray.env.cap)$r.squared)
# Adjusted R^2 retrieved from the rda object
(fin_bray.env.cap_R2adj <- RsquareAdj(fin_bray.env.cap)$adj.r.squared)

## Plot with "wa" scores to see dispersion of sites around the factor levels
triplot.rda(fin_bray.env.cap, site.sc = "wa", scaling = 1)
triplot.rda(fin_bray.env.cap, site.sc = "wa", scaling = 2)

## Look at the variance inflation factors, values > 10 indicate there is redundancy in the model.
vif.cca(fin_bray.env.cap)

## Rerun dbRDA & Capscale with removed variables
# dbrda() on the square-rooted dissimilarity matrix
fin_bray.env.dbrda.base <- dbrda(sqrt(fin_comm_bray) ~ 
                                   Altitude + 
                                   pH + 
                                   DOC + 
                                   Tot.N +
                                   Tot.P + 
                                   Season + 
                                   Condition(fin_dist_pco),
                                 data = as.data.frame(fin_env_sel),
                                 add = FALSE)
# Model output
summary(fin_bray.env.dbrda.base, axes = FALSE)
# Global test of the RDA result
anova(fin_bray.env.dbrda.base, permutations = how(nperm = 999))
# Tests of all canonical axes
anova(fin_bray.env.dbrda.base, by = "axis", permutations = how(nperm = 999))
# Tests of all explanatory terms
anova(fin_bray.env.dbrda.base, by = "terms", permutations = how(nperm = 999))

## capscale() with raw (site by species) data
fin_bray.env.cap.base <- capscale(fin_comm ~ 
                                    Altitude + 
                                    pH + 
                                    DOC + 
                                    Tot.N +
                                    Tot.P + 
                                    Season + 
                                    Condition(fin_dist_pco),
                                  data = as.data.frame(fin_env_sel),
                                  distance = "bray",
                                  add = "lingoes",
                                  comm = fin_comm)
# Model output
summary(fin_bray.env.cap.base, axes = FALSE)
# Global test of the RDA result
anova(fin_bray.env.cap.base, permutations = how(nperm = 999))
# Tests of all canonical axes
anova(fin_bray.env.cap.base, by = "axis", permutations = how(nperm = 999))
# Tests of all explanatory terms
anova(fin_bray.env.cap.base, by = "terms", permutations = how(nperm = 999))
# Look at the variance inflation factors, values > 10 indicate there is redundancy in the model.
vif.cca(fin_bray.env.cap.base)
# Canonical coefficients from the rda object
coef(fin_bray.env.cap.base)
# Unadjusted R^2 retrieved from the rda object
(fin_bray.env.cap.base_R2 <- RsquareAdj(fin_bray.env.cap.base)$r.squared)
# Adjusted R^2 retrieved from the rda object
(fin_bray.env.cap.base_R2adj <- RsquareAdj(fin_bray.env.cap.base)$adj.r.squared)
# Plot with "wa" scores to see dispersion of sites around the factor levels
triplot.rda(fin_bray.env.cap.base, site.sc = "wa", scaling = 1)
triplot.rda(fin_bray.env.cap.base, site.sc = "wa", scaling = 2)

## Variable selection
## Forward selection using vegan's ordiR2step()
# Ordination with all explanatory variables (full model)
fin_cap_par_full <- capscale(fin_comm ~ . + Condition(fin_dist_pco),
                             data = as.data.frame(fin_env_sel),
                             distance = "bray",
                             add = "lingoes",
                             comm = fin_comm)
(fin_cap_par_full_R2a <- RsquareAdj(fin_cap_par_full)$adj.r.squared)
# Ordination with no explanatory variables except conditional variables (null model)
fin_cap_par_null <- capscale(fin_comm ~ 1 + Condition(fin_dist_pco),
                             data = as.data.frame(fin_env_sel),
                             distance = "bray",
                             add = "lingoes",
                             comm = fin_comm)
# Forward selection
(fin_cap_par.step2.forward <- ordiR2step(fin_cap_par_null,
                                         scope = formula(fin_cap_par_full),
                                         R2scope = fin_cap_par_full_R2a,
                                         direction = "forward",
                                         permutations = how(nperm = 999)))
summary(fin_cap_par.step2.forward, axes = FALSE)
RsquareAdj(fin_cap_par.step2.forward) # Condition(fin_dist_pco) + pH

## Check of explanatory variables against full model
vif.cca(fin_bray.env.cap.base)
vif.cca(fin_cap_par.step2.forward)

## ANOVA: test top models
fin_cap_par_red <- capscale(fin_comm ~ 
                              pH +
                              Condition(fin_dist_pco),
                            data = as.data.frame(fin_env_sel),
                            distance = "bray",
                            add = "lingoes",
                            comm = fin_comm)
# Model output
summary(fin_cap_par_red, axes = FALSE)
# Global test of the RDA result
anova(fin_cap_par_red, permutations = how(nperm = 999))
# Tests of all canonical axes
anova(fin_cap_par_red,permutations = how(nperm = 999), by = "axis")
# Tests of all explanatory terms
anova(fin_cap_par_red, permut = 999, by = "terms")
# Unadjusted R^2 retrieved from the rda object
(fin_cap_par_red_R2 <- RsquareAdj(fin_cap_par_red)$r.squared)
# Adjusted R^2 retrieved from the rda object
(fin_cap_par_red_R2a <- RsquareAdj(fin_cap_par_red)$adj.r.squared)
# Plot with "wa" scores to see dispersion of sites around the factor levels
# Will not run with only one variable - see other plotting code below
plot(fin_cap_par_red, site.sc = "wa", scaling = 1)
plot(fin_cap_par_red, site.sc = "wa", scaling = 2)

## Parsimonious RDA
## dbRDA (provides test)
fin_dbrda_pars <- dbrda(sqrt(fin_comm_bray) ~ 
                          pH + 
                          Condition(fin_dist_pco),
                        data = as.data.frame(fin_env_sel),
                        add = FALSE)
summary(fin_dbrda_pars, axes = FALSE)
set.seed(1)
anova(fin_dbrda_pars, permutations = how(nperm = 999))
set.seed(1)
anova(fin_dbrda_pars, permutations = how(nperm = 999), by = "axis")
set.seed(1)
anova(fin_dbrda_pars, permutations = how(nperm = 999), by = "terms")
(fin_dbrda_pars_R2 <- RsquareAdj(fin_dbrda_pars)$r.squared)
(fin_dbrda_pars_R2a <- RsquareAdj(fin_dbrda_pars)$adj.r.squared)

## capscale (provides material for plotting)
fin_cap_pars <- capscale(fin_comm ~ 
                           pH +
                           Condition(fin_dist_pco), 
                         data = as.data.frame(fin_env_sel),
                         distance = "bray",
                         add = "lingoes",
                         comm = fin_comm)
summary(fin_cap_pars, axes = FALSE)
anova(fin_cap_pars, permutations = how(nperm = 999))
anova(fin_cap_pars, permutations = how(nperm = 999), by = "axis")
anova(fin_cap_pars, permutations = how(nperm = 999), by = "terms")
(fin_cap_pars_R2 <- RsquareAdj(fin_cap_pars)$r.squared)
(fin_cap_pars_R2a <- RsquareAdj(fin_cap_pars)$adj.r.squared)

# Plot with "wa" scores to see dispersion of sites around the factor levels
plot(fin_cap_pars, site.sc = "wa", scaling = 1)
plot(fin_cap_pars, site.sc = "wa", scaling = 2)

op <- par(mfrow = c(1, 2))
# Scaling 1
fin_fig1 <- ordiplot(fin_cap_pars, type = "none", scaling = 1, main = "RDA triplot - Scaling 1 - wa", 
                     xlim = c(-1.0, 3.0), ylim = c(-3.0, 1.0), cex.main = 1, xlab = "RDA 1", ylab = "RDA 2")
ordihull(fin_cap_pars, fin_env_sel$Season, draw = "polygon", lty = 1, col = c("#01698e64", "#658c6564"), scaling = 1, label = F)
# point symbols according to a specific factor
n <- points(fin_fig1, "sites", type = "n")
scores <- n$sites
s1 <-scores[ ,1]
s2 <-scores[ ,2]
listy1 = fin_env_sel$Season == "Spring" # make a list of T or F, for you this will be e.g. "spring" instead of "1"
s1_1 = s1[listy1] # subset by this list
s2_1 = s2[listy1] # subset by this list
listy2 = fin_env_sel$Season == "Autumn" # e.g. Autumn instead of 2
s1_2 = s1[listy2] # subset by this list
s2_2 = s2[listy2] # subset by this list
points(x = s1_1, y = s2_1, display = "wa", pch = 19, col = "black", cex = 1.0, scaling = 1, pos = 4)
points(x = s1_2, y = s2_2, display = "wa", pch = 17, col = "black", cex = 1.0, scaling = 1, pos = 4)
text(fin_fig1, "species", cex = 0.8, col = "red", pos = 4, scaling = 1)
spe.sc <- scores(fin_cap_pars, choices = 1:2, display = "sp", scaling = 1)
arrows(0,0, spe.sc[,1], spe.sc[,2], length = 0.05, angle = 30, lty = 1, col = 'red', scaling = 1)
text(fin_cap_pars, display = "wa", pos = 2, cex = 0.7, scaling = 1, col = "black")
(fit <- envfit(fin_cap_pars, fin_env_sel[, 2], perm = 999, display = "lc", scaling = 1))
scores(fit, "vectors", scaling = 1)
plot(fit, p.max = 0.05, col = "blue", cex = 1.0, lab = "pH")
legend("bottomright", legend = c("Spring","Autumn"), bty = "n", cex = 1.0, pch = c(19, 17), col =  c("black", "black"))

# Scaling 2
fin_fig1 <- ordiplot(fin_cap_pars, type = "none", scaling = 2, main = "RDA triplot - Scaling 2 - wa", 
                     xlim = c(-2.0, 2.0), ylim = c(-1.5, 1.5), cex.main = 1, xlab = "RDA 1", ylab = "RDA 2")
ordihull(fin_cap_pars, fin_env_sel$Season, draw = "polygon", lty = 1, col = c("#01698e64", "#658c6564"), scaling = 2, label = F)
n <- points(fin_fig1, "sites", type = "n")
scores <- n$sites
s1 <-scores[ ,1]
s2 <-scores[ ,2]
listy1 = fin_env_sel$Season == "Spring" #make a list of T or F, for you this will be e.g. "spring" instead of "1"
s1_1 = s1[listy1] #subset by this list
s2_1 = s2[listy1] #subset by this list
listy2 = fin_env_sel$Season == "Autumn" # e.g. Fall instead of 2
s1_2 = s1[listy2] #subset by this list
s2_2 = s2[listy2] #subset by this list
points(x = s1_1, y = s2_1, display = "wa", pch = 19, col = "black", cex = 1.0, scaling = 2, pos = 4)
points(x = s1_2, y = s2_2, display = "wa", pch = 17, col = "black", cex = 1.0, scaling = 2, pos = 4)
text(fin_fig1, "species", cex = 0.8, col = "red", pos = 4, scaling = 2)
spe.sc <- scores(fin_cap_pars, choices = 1:2, display = "sp", scaling = 2)
arrows(0,0, spe.sc[,1], spe.sc[,2], length = 0.05, angle = 30, lty = 1, col = 'red', scaling = 2)
text(fin_cap_pars, display = "wa", pos = 2, cex = 0.7, scaling = 2, col = "black")
(fit <- envfit(fin_cap_pars, fin_env_sel[, 2], perm = 999, display = "lc", scaling = 2))
scores(fit, "vectors", scaling = 2)
plot(fit, p.max = 0.05, col = "blue", cex = 1.0, lab = "pH")
legend("bottomright", legend = c("Spring","Autumn"), bty = "n", cex = 1.0, pch = c(19, 17), col =  c("black", "black"))
par(op)

### FUNCTIONAL STRUCTURE --------------------
## Organize data
## Distance matrix for covariable (i.e. the variable used as a condition factor)
fin_dist_pco

## Environmental data (predictor variables)
fin_env_sel

## Functional data (response variables)
# Create regional trait and abundance datasets
# funcitonal abundance data
fin_func_comm_full <- func_comm[121:140, ]
fin_func_comm_full <- fin_func_comm_full[, colSums(fin_func_comm_full) > 0] # only keeps columns with column sums > 0
# functional trait data
spp.fin.full = (row.names(comm_traits) %in% colnames(fin_func_comm_full)) # make a list of T or F is there is a match in spp
fin_traits.full = comm_traits[spp.fin.full, ] # subset by this list
# Check data
identical(row.names(fin_traits.full), colnames(fin_func_comm_full)) # check if names are the same
# Remove rare species
# Look at the data, eliminate species with less than 2 total abundance
fin_num_occur_species_func <- apply(fin_func_comm_full > 0,2, sum) 
fin_sorted_abundance_func <- sort(fin_num_occur_species_func)
# keep only spp that occur in more than 2 sites
fin_func_comm <- fin_func_comm_full[, fin_num_occur_species_func >= 1]
# create trait subset from dataset with rare taxa removed
spp.fin <- (row.names(fin_traits.full) %in% colnames(fin_func_comm))
fin_traits = fin_traits.full[spp.fin, ]
#Check if data is correct & that row names match
identical(row.names(fin_traits), colnames(fin_func_comm)) # check if names are the same
# Compute the functional dissimilarity response matrix
fin_gower <- gowdis(fin_traits) # gower disimilarity matrix
plot(fin_clust <- hclust(fin_gower, method = "ward.D"))
fin_func_beta <- beta(fin_func_comm, fin_clust, func = "Soerensen", abund = TRUE) # percentage difference
(fin_func_mat <- fin_func_beta$Btotal)

## PCOA of raw community data (unconstrained)
# Run a PCoA - distance based ordination with no environmental constraints (i.e. unconstrained analysis)
is.euclid(sqrt(fin_func_mat))
fin_soren.pcoa <- cmdscale(sqrt(fin_func_mat), k = (nrow(fin_func_comm) - 1), eig = TRUE)
## Plot of the sites
ordiplot(scores(fin_soren.pcoa, choices = c(1, 2)),
         type = "t",
         main = "PCoA with species weighted averages")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
## Add weighted average projection of species
fin_func.wa <- wascores(fin_soren.pcoa$points[, 1:2], fin_func_comm)
text(fin_func.wa, rownames(fin_func.wa), cex = 0.7, col = "red")
## A posteriori projection of environmental variables
(fin_soren.pcoa.env <- envfit(fin_soren.pcoa, fin_env_sel))
## Plot significant variables with a user-selected colour
plot(fin_soren.pcoa.env, p.max = 0.05, col = 3)

## Run dbRDA & Capscale
## dbrda() on the square-rooted dissimilarity matrix
fin_soren.env.dbrda.base <- dbrda(sqrt(fin_func_mat) ~ 
                                    Altitude + 
                                    pH + 
                                    DOC + 
                                    Tot.N +
                                    Tot.P +
                                    Season +
                                    Condition(fin_dist_pco),
                                  data = as.data.frame(fin_env_sel),
                                  add = FALSE)
# Model output
summary(fin_soren.env.dbrda.base, axes = FALSE)
# Global test of the RDA result
anova(fin_soren.env.dbrda.base, permutations = how(nperm = 999))
# Tests of all canonical axes
anova(fin_soren.env.dbrda.base, by = "axis", permutations = how(nperm = 999))
# Tests of all explanatory terms
anova(fin_soren.env.dbrda.base, by = "terms", permutations = how(nperm = 999))
# Unadjusted R^2 retrieved from the rda object
(fin_soren.env.dbrda.base_R2 <- RsquareAdj(fin_soren.env.dbrda.base)$r.squared)
# Adjusted R^2 retrieved from the rda object
(fin_soren.env.dbrda.base_R2adj <- RsquareAdj(fin_soren.env.dbrda.base)$adj.r.squared)

## capscale() with raw (site by species) data
fin_soren.env.cap.base <- capscale(fin_func_mat ~ 
                                     Altitude + 
                                     pH + 
                                     DOC + 
                                     Tot.N +
                                     Tot.P +
                                     Season +
                                     Condition(fin_dist_pco),
                                   data = as.data.frame(fin_env_sel),
                                   add = "lingoes",
                                   comm = fin_func_comm)
# Model output
summary(fin_soren.env.cap.base, axes = FALSE)
# Global test of the RDA result
anova(fin_soren.env.cap.base, permutations = how(nperm = 999))
# Tests of all canonical axes
anova(fin_soren.env.cap.base, by = "axis", permutations = how(nperm = 999))
# Tests of all explanatory terms
anova(fin_soren.env.cap.base, by = "terms", permutations = how(nperm = 999))
# Look at the variance inflation factors, values > 10 indicate there is redundancy in the model.
vif.cca(fin_soren.env.cap.base)
# Canonical coefficients from the rda object
coef(fin_soren.env.cap.base)
# Unadjusted R^2 retrieved from the rda object
(fin_soren.env.cap.base_R2 <- RsquareAdj(fin_soren.env.cap.base)$r.squared)
# Adjusted R^2 retrieved from the rda object
(fin_soren.env.cap.base_R2adj <- RsquareAdj(fin_soren.env.cap.base)$adj.r.squared)
# Plot with "wa" scores to see dispersion of sites around the factor levels
triplot.rda(fin_soren.env.cap.base, site.sc = "wa", scaling = 1)
triplot.rda(fin_soren.env.cap.base, site.sc = "wa", scaling = 2)

## Variable selection
## Forward selection using vegan's ordiR2step()
# Ordination with all explanatory variables (full model)
fin_soren_cap_par_full <- capscale(fin_func_mat ~ . + Condition(fin_dist_pco),
                                   data = as.data.frame(fin_env_sel),
                                   add = "lingoes",
                                   comm = fin_func_comm)
(fin_soren_cap_par_full_R2a <- RsquareAdj(fin_soren_cap_par_full)$adj.r.squared)
# Ordination with no explanatory variables except conditional variables (null model)
fin_soren_cap_par_null <- capscale(fin_func_mat ~ 1 + Condition(fin_dist_pco),
                                   data = as.data.frame(fin_env_sel),
                                   add = "lingoes",
                                   comm = fin_func_comm)
# Forward selection
(fin_soren_cap_par.step2.forward <- ordiR2step(fin_soren_cap_par_null,
                                               scope = formula(fin_soren_cap_par_full),
                                               R2scope = fin_soren_cap_par_full_R2a,
                                               direction = "forward",
                                               permutations = how(nperm = 999)))
summary(fin_soren_cap_par.step2.forward, axes = FALSE)
RsquareAdj(fin_soren_cap_par.step2.forward) # Condition(fin_pco_axis1 + fin_pco_axis2)

## Check of explanatory variables against full model
vif.cca(fin_soren.env.cap.base)
#vif.cca(fin_soren_cap_par.step2.forward)

## ANOVA: test top models
fin_soren_cap_par_red <- capscale(fin_func_mat ~ 
                                    1 + Condition(fin_dist_pco),
                                  data = as.data.frame(fin_env_sel),
                                  add = "lingoes",
                                  comm = fin_func_comm)
# Model output
summary(fin_soren_cap_par_red, axes = FALSE)
# Global test of the RDA result
anova(fin_soren_cap_par_red, permutations = how(nperm = 999))
# Tests of all canonical axes
anova(fin_soren_cap_par_red,permutations = how(nperm = 999), by = "axis")
# Tests of all explanatory terms
anova(fin_soren_cap_par_red, permut = 999, by = "terms")
# Unadjusted R^2 retrieved from the rda object
(fin_soren_cap_par_red_R2 <- RsquareAdj(fin_soren_cap_par_red)$r.squared)
# Adjusted R^2 retrieved from the rda object
(fin_soren_cap_par_red_R2a <- RsquareAdj(fin_soren_cap_par_red)$adj.r.squared)
# Plot with "wa" scores to see dispersion of sites around the factor levels
plot(fin_soren_cap_par_red, site.sc = "wa", scaling = 1)
plot(fin_soren_cap_par_red, site.sc = "wa", scaling = 2)

## Parsimonious RDA
## dbRDA (provides test)
fin_func_dbrda_pars <- dbrda(sqrt(fin_func_mat) ~ 
                               1 + Condition(fin_dist_pco),
                             data = as.data.frame(fin_env_sel),
                             add = FALSE)
summary(fin_func_dbrda_pars, axes = FALSE)
set.seed(1)
anova(fin_func_dbrda_pars, permutations = how(nperm = 999))
set.seed(1)
anova(fin_func_dbrda_pars, permutations = how(nperm = 999), by = "axis")
set.seed(1)
anova(fin_func_dbrda_pars, permutations = how(nperm = 999), by = "terms")
(fin_func_dbrda_pars_R2 <- RsquareAdj(fin_func_dbrda_pars)$r.squared)
(fin_func_dbrda_pars_R2a <- RsquareAdj(fin_func_dbrda_pars)$adj.r.squared)

## capscale (provides material for plotting)
fin_func_cap_pars <- capscale(fin_func_mat ~ 
                                1 + Condition(fin_dist_pco), 
                              data = as.data.frame(fin_env_sel),
                              add = "lingoes",
                              comm = fin_func_comm)
summary(fin_func_cap_pars, axes = FALSE)
anova(fin_func_cap_pars, permutations = how(nperm = 999))
anova(fin_func_cap_pars, permutations = how(nperm = 999), by = "axis")
anova(fin_func_cap_pars, permutations = how(nperm = 999), by = "terms")
(fin_func_cap_pars_R2 <- RsquareAdj(fin_func_cap_pars)$r.squared)
(fin_func_cap_pars_R2a <- RsquareAdj(fin_func_cap_pars)$adj.r.squared)

# Plot with "wa" scores to see dispersion of sites around the factor levels
plot(fin_func_cap_pars, site.sc = "wa", scaling = 1)
plot(fin_func_cap_pars, site.sc = "wa", scaling = 2)

op <- par(mfrow = c(1, 2))
#Scaling 1
fin_fig1 <- ordiplot(fin_func_cap_pars, type = "none", scaling = 1, main = "RDA triplot - Scaling 1 - wa", 
                     xlim = c(-3.0, 1.0), ylim = c(-0.5, 3.0), cex.main = 1, xlab = "RDA 1", ylab = "RDA 2")
ordihull(fin_func_cap_pars, fin_env_sel$Season, draw = "polygon", lty = 1, col = c("#01698e64", "#658c6564"), scaling = 1, label = F)
# point symbols according to a specific factor
n <- points(fin_fig1, "sites", type = "n")
scores <- n$sites
s1 <-scores[ ,1]
s2 <-scores[ ,2]
listy1 = fin_env_sel$Season == "Spring" # make a list of T or F, for you this will be e.g. "spring" instead of "1"
s1_1 = s1[listy1] # subset by this list
s2_1 = s2[listy1] # subset by this list
listy2 = fin_env_sel$Season == "Autumn" # e.g. Autumn instead of 2
s1_2 = s1[listy2] # subset by this list
s2_2 = s2[listy2] # subset by this list
points(x = s1_1, y = s2_1, display = "wa", pch = 19, col = "black", cex = 1.0, scaling = 1, pos = 4)
points(x = s1_2, y = s2_2, display = "wa", pch = 17, col = "black", cex = 1.0, scaling = 1, pos = 4)
text(fin_fig1, "species", cex = 0.8, col = "red", pos = 4, scaling = 1)
spe.sc <- scores(fin_func_cap_pars, choices = 1:2, display = "sp", scaling = 1)
arrows(0,0, spe.sc[,1], spe.sc[,2], length = 0.05, angle = 30, lty = 1, col = 'red', scaling = 1)
text(fin_func_cap_pars, display = "wa", pos = 2, cex = 0.7, scaling = 1, col = "black")
#(fit <- envfit(fin_func_cap_pars, fin_env_sel[, 2], perm = 999, display = "lc", scaling = 1))
#scores(fit, "vectors", scaling = 1)
#plot(fit, p.max = 0.05, col = "blue", cex = 1.0, lab = "pH")
legend("bottomright", legend = c("Spring","Autumn"), bty = "n", cex = 1.0, pch = c(19, 17), col =  c("black", "black"))

#Scaling 2
fin_fig2 <- ordiplot(fin_func_cap_pars, type = "none", scaling = 2, main = "RDA triplot - Scaling 2 - wa", 
                     xlim = c(-2.0, 2.0), ylim = c(-1.0, 2.5), cex.main = 1, xlab = "RDA 1", ylab = "RDA 2")
ordihull(fin_func_cap_pars, fin_env_sel$Season, draw = "polygon", lty = 1, col = c("#01698e64", "#658c6564"), scaling = 2, label = F)
n <- points(fin_fig2, "sites", type = "n")
scores <- n$sites
s1 <-scores[ ,1]
s2 <-scores[ ,2]
listy1 = fin_env_sel$Season == "Spring" #make a list of T or F, for you this will be e.g. "spring" instead of "1"
s1_1 = s1[listy1] #subset by this list
s2_1 = s2[listy1] #subset by this list
listy2 = fin_env_sel$Season == "Autumn" # e.g. Fall instead of 2
s1_2 = s1[listy2] #subset by this list
s2_2 = s2[listy2] #subset by this list
points(x = s1_1, y = s2_1, display = "wa", pch = 19, col = "black", cex = 1.0, scaling = 2, pos = 4)
points(x = s1_2, y = s2_2, display = "wa", pch = 17, col = "black", cex = 1.0, scaling = 2, pos = 4)
text(fin_fig2, "species", cex = 0.8, col = "red", pos = 4, scaling = 2)
spe.sc <- scores(fin_func_cap_pars, choices = 1:2, display = "sp", scaling = 2)
arrows(0,0, spe.sc[,1], spe.sc[,2], length = 0.05, angle = 30, lty = 1, col = 'red', scaling = 2)
text(fin_func_cap_pars, display = "wa", pos = 2, cex = 0.7, scaling = 2, col = "black")
#(fit <- envfit(fin_func_cap_pars, fin_env_sel[, 2], perm = 999, display = "lc", scaling = 2))
#scores(fit, "vectors", scaling = 2)
#plot(fit, p.max = 0.05, col = "blue", cex = 1.0, lab = "pH")
legend("bottomright", legend = c("Spring","Autumn"), bty = "n", cex = 1.0, pch = c(19, 17), col =  c("black", "black"))
par(op)

#### FIGURE 4A-F: DB-RDA PLOTS OF COMPOSITION DATA -----------------
### Scaling 2
## Open new plot window
# dev.new(title = "Capscale plot with all species.rda",
#        width = 12,
#        height = 18,
#        noRStudioGD = TRUE)
## Set image properties
pdf(file = "FIGURE 4 - dbRDA PLOTS.pdf", width = 12, height = 16) 
# svg(file = "All_dbRDAs - Community composition - R1_1.svg", width = 12, height = 16) 
## Set layout properties
op <- par(mfrow = c(3, 2))

## Germany
# Taxonomic structure
# Scaling 2
ger_fig1 <- ordiplot(ger_cap_pars, type = "none", scaling = 2, main = "RDA triplot - Scaling 2 - wa", 
                     xlim = c(-1.5, 1.5), ylim = c(-2.0, 2.0), cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, xlab = "RDA 1", ylab = "RDA 2")
ordihull(ger_cap_pars, ger_env_sel$Season, draw = "polygon", lty = 0, col = c("#00688c4d", "#638a634d"), scaling = 2, label = F)
n <- points(ger_fig1, "sites", type = "n")
scores <- n$sites
s1 <-scores[ ,1]
s2 <-scores[ ,2]
listy1 = ger_env_sel$Season == "Spring" #make a list of T or F, for you this will be e.g. "spring" instead of "1"
s1_1 = s1[listy1] #subset by this list
s2_1 = s2[listy1] #subset by this list
listy2 = ger_env_sel$Season == "Autumn" # e.g. Fall instead of 2
s1_2 = s1[listy2] #subset by this list
s2_2 = s2[listy2] #subset by this list
points(x = s1_1, y = s2_1, display = "wa", pch = 19, col = "black", cex = 1.5, scaling = 2, pos = 4)
points(x = s1_2, y = s2_2, display = "wa", pch = 17, col = "black", cex = 1.5, scaling = 2, pos = 4)
text(ger_cap_pars, display = "wa", pos = 2, cex = 0.8, scaling = 2, col = "black")
text(ger_cap_pars, dis = "cn", cex = 1.2, font = 2, col = "blue")
#text(ger_fig1, "species", cex = 1.2, col = "red", pos = 4, scaling = 2)
orditorp(ger_fig1, display = "species", priority = colSums(ger_comm), cex = 1, pcex = 0.5, 
         col = "red", pcol="grey40", pch = 3, air = 1, scaling = 2)
spe.sc <- scores(ger_cap_pars, choices = 1:2, display = "sp", scaling = 2)
arrows(0,0, spe.sc[,1], spe.sc[,2], length = 0.05, angle = 30, lty = 1, col = 'red', scaling = 2)
legend("bottomright", legend = c("Spring","Autumn"), bty = "n", cex = 1.5, pch = c(19, 17), col =  c("black", "black"))

# Functional structure
ger_fig2 <- ordiplot(ger_func_cap_pars, type = "none", scaling = 2, main = "RDA triplot - Scaling 2 - wa", 
                     xlim = c(-1.5, 1.5), ylim = c(-2.0, 2.0), cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, xlab = "RDA 1", ylab = "RDA 2")
ordihull(ger_func_cap_pars, ger_env_sel$Season, draw = "polygon", lty = 0, col = c("#00688c4d", "#638a634d"), scaling = 2, label = F)
n <- points(ger_fig2, "sites", type = "n")
scores <- n$sites
s1 <-scores[ ,1]
s2 <-scores[ ,2]
listy1 = ger_env_sel$Season == "Spring" #make a list of T or F, for you this will be e.g. "spring" instead of "1"
s1_1 = s1[listy1] #subset by this list
s2_1 = s2[listy1] #subset by this list
listy2 = ger_env_sel$Season == "Autumn" # e.g. Fall instead of 2
s1_2 = s1[listy2] #subset by this list
s2_2 = s2[listy2] #subset by this list
points(x = s1_1, y = s2_1, display = "wa", pch = 19, col = "black", cex = 1.5, scaling = 2, pos = 4)
points(x = s1_2, y = s2_2, display = "wa", pch = 17, col = "black", cex = 1.5, scaling = 2, pos = 4)
#text(ger_fig2, "species", cex = 1.2, col = "red", pos = 4, scaling = 2)
orditorp(ger_fig2, display = "species", priority = colSums(ger_func_comm), cex = 1, pcex = 0.5, 
         col = "red", pcol="grey40", pch = 3, air = 1, scaling = 2)
spe.sc <- scores(ger_func_cap_pars, choices = 1:2, display = "sp", scaling = 2)
arrows(0,0, spe.sc[,1], spe.sc[,2], length = 0.05, angle = 30, lty = 1, col = 'red', scaling = 2)
text(ger_func_cap_pars, display = "wa", pos = 2, cex = 0.8, scaling = 2, col = "black")
legend("bottomright", legend = c("Spring","Autumn"), bty = "n", cex = 1.5, pch = c(19, 17), col =  c("black", "black"))

## Sweden
# Taxonomic structure
swe_fig1 <- ordiplot(swe_cap_pars, type = "none", scaling = 2, main = "RDA triplot - Scaling 2 - wa", 
                     xlim = c(-3.0, 3.0), ylim = c(-3.0, 3.0), cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, xlab = "RDA 1", ylab = "RDA 2")
ordihull(swe_cap_pars, swe_env_sel$Season, draw = "polygon", lty = 0, col = c("#00688c4d", "#638a634d"), scaling = 2, label = F)
n <- points(swe_fig1, "sites", type = "n")
scores <- n$sites
s1 <-scores[ ,1]
s2 <-scores[ ,2]
listy1 = swe_env_sel$Season == "Spring" #make a list of T or F, for you this will be e.g. "spring" instead of "1"
s1_1 = s1[listy1] #subset by this list
s2_1 = s2[listy1] #subset by this list
listy2 = swe_env_sel$Season == "Autumn" # e.g. Fall instead of 2
s1_2 = s1[listy2] #subset by this list
s2_2 = s2[listy2] #subset by this list
points(x = s1_1, y = s2_1, display = "wa", pch = 19, col = "black", cex = 1.5, scaling = 2, pos = 4)
points(x = s1_2, y = s2_2, display = "wa", pch = 17, col = "black", cex = 1.5, scaling = 2, pos = 4)
text(swe_cap_pars, display = "wa", pos = 2, cex = 0.8, scaling = 2, col = "black")
text(swe_cap_pars, dis = "cn", col = "blue", cex = 1.2, font = 2)
#text(swe_fig1, "species", cex = 1.2, col = "red", pos = 4, scaling = 2)
orditorp(swe_fig1, display = "species", priority = colSums(swe_comm), cex = 1, pcex = 0.5, 
         col = "red", pcol="grey40", pch = 3, air = 1, scaling = 2)
spe.sc <- scores(swe_cap_pars, choices = 1:2, display = "sp", scaling = 2)
arrows(0,0, spe.sc[,1], spe.sc[,2], length = 0.05, angle = 30, lty = 1, col = 'red', scaling = 2)
legend("bottomright", legend = c("Spring","Autumn"), bty = "n", cex = 1.5, pch = c(19, 17), col =  c("black", "black"))

# functional structure
swe_fig1 <- ordiplot(swe_func_cap_pars, type = "none", scaling = 2, main = "RDA triplot - Scaling 2 - wa", 
                     xlim = c(-4.0, 4.0), ylim = c(-5.0, 5.0), cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, xlab = "RDA 1", ylab = "RDA 2")
ordihull(swe_func_cap_pars, swe_env_sel$Season, draw = "polygon", lty = 0, col = c("#00688c4d", "#638a634d"), scaling = 2, label = F)
n <- points(swe_fig1, "sites", type = "n")
scores <- n$sites
s1 <-scores[ ,1]
s2 <-scores[ ,2]
listy1 = swe_env_sel$Season == "Spring" #make a list of T or F, for you this will be e.g. "spring" instead of "1"
s1_1 = s1[listy1] #subset by this list
s2_1 = s2[listy1] #subset by this list
listy2 = swe_env_sel$Season == "Autumn" # e.g. Fall instead of 2
s1_2 = s1[listy2] #subset by this list
s2_2 = s2[listy2] #subset by this list
points(x = s1_1, y = s2_1, display = "wa", pch = 19, col = "black", cex = 1.5, scaling = 2, pos = 4)
points(x = s1_2, y = s2_2, display = "wa", pch = 17, col = "black", cex = 1.5, scaling = 2, pos = 4)
text(swe_func_cap_pars, display = "wa", pos = 2, cex = 0.8, scaling = 2, col = "black")
text(swe_func_cap_pars, dis = "cn", col = "blue", cex = 1.2, font = 2)
#text(swe_fig1, "species", cex = 1.2, col = "red", pos = 4, scaling = 2)
orditorp(swe_fig1, display = "species", priority = colSums(swe_func_comm), cex = 1, pcex = 0.5, 
         col = "red", pcol="grey40", pch = 3, air = 1, scaling = 2)
spe.sc <- scores(swe_func_cap_pars, choices = 1:2, display = "sp", scaling = 2)
arrows(0,0, spe.sc[,1], spe.sc[,2], length = 0.05, angle = 30, lty = 1, col = 'red', scaling = 2)
legend("bottomright", legend = c("Spring","Autumn"), bty = "n", cex = 1.5, pch = c(19, 17), col =  c("black", "black"))

## Finland
# Taxonomic structure
fin_fig1 <- ordiplot(fin_cap_pars, type = "none", scaling = 2, main = "RDA triplot - Scaling 2 - wa", 
                     xlim = c(-2.0, 2.0), ylim = c(-2.0, 2.0), cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, xlab = "RDA 1", ylab = "RDA 2")
ordihull(fin_cap_pars, fin_env_sel$Season, draw = "polygon", lty = 0, col = c("#00688c4d", "#638a634d"), scaling = 2, label = F)
n <- points(fin_fig1, "sites", type = "n")
scores <- n$sites
s1 <-scores[ ,1]
s2 <-scores[ ,2]
listy1 = fin_env_sel$Season == "Spring" #make a list of T or F, for you this will be e.g. "spring" instead of "1"
s1_1 = s1[listy1] #subset by this list
s2_1 = s2[listy1] #subset by this list
listy2 = fin_env_sel$Season == "Autumn" # e.g. Fall instead of 2
s1_2 = s1[listy2] #subset by this list
s2_2 = s2[listy2] #subset by this list
points(x = s1_1, y = s2_1, display = "wa", pch = 19, col = "black", cex = 1.5, scaling = 2, pos = 4)
points(x = s1_2, y = s2_2, display = "wa", pch = 17, col = "black", cex = 1.5, scaling = 2, pos = 4)
text(fin_cap_pars, display = "wa", pos = 2, cex = 0.8, scaling = 2, col = "black")
(fit <- envfit(fin_cap_pars, fin_env_sel[, 2], perm = 999, display = "lc", scaling = 2))
scores(fit, "vectors", scaling = 2)
plot(fit, p.max = 0.05, col = "blue", cex = 1.2, font = 2, lab = "pH")
#text(fin_fig1, "species", cex = 1.2, col = "red", pos = 4, scaling = 2)
orditorp(fin_fig1, display = "species", priority = colSums(fin_comm), cex = 1, pcex = 0.5, 
         col = "red", pcol="grey40", pch = 3, air = 1, scaling = 2)
spe.sc <- scores(fin_cap_pars, choices = 1:2, display = "sp", scaling = 2)
arrows(0,0, spe.sc[,1], spe.sc[,2], length = 0.05, angle = 30, lty = 1, col = 'red', scaling = 2)
legend("bottomright", legend = c("Spring","Autumn"), bty = "n", cex = 1.5, pch = c(19, 17), col =  c("black", "black"))

# Functional structure
fin_fig2 <- ordiplot(fin_func_cap_pars, type = "none", scaling = 2, main = "RDA triplot - Scaling 2 - wa", 
                     xlim = c(-2.0, 2.0), ylim = c(-1.0, 3.0), cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, xlab = "RDA 1", ylab = "RDA 2")
ordihull(fin_func_cap_pars, fin_env_sel$Season, draw = "polygon", lty = 0, col = c("#00688c4d", "#638a634d"), scaling = 2, label = F)
n <- points(fin_fig2, "sites", type = "n")
scores <- n$sites
s1 <-scores[ ,1]
s2 <-scores[ ,2]
listy1 = fin_env_sel$Season == "Spring" #make a list of T or F, for you this will be e.g. "spring" instead of "1"
s1_1 = s1[listy1] #subset by this list
s2_1 = s2[listy1] #subset by this list
listy2 = fin_env_sel$Season == "Autumn" # e.g. Fall instead of 2
s1_2 = s1[listy2] #subset by this list
s2_2 = s2[listy2] #subset by this list
points(x = s1_1, y = s2_1, display = "wa", pch = 19, col = "black", cex = 1.5, scaling = 2, pos = 4)
points(x = s1_2, y = s2_2, display = "wa", pch = 17, col = "black", cex = 1.5, scaling = 2, pos = 4)
#text(fin_fig2, "species", cex = 1.2, col = "red", pos = 4, scaling = 2)
orditorp(fin_fig2, display = "species", priority = colSums(fin_func_comm), cex = 1, pcex = 0.5, 
         col = "red", pcol="grey40", pch = 3, air = 1, scaling = 2)
spe.sc <- scores(fin_func_cap_pars, choices = 1:2, display = "sp", scaling = 2)
arrows(0,0, spe.sc[,1], spe.sc[,2], length = 0.05, angle = 30, lty = 1, col = 'red', scaling = 2)
text(fin_func_cap_pars, display = "wa", pos = 2, cex = 0.8, scaling = 2, col = "black")
#(fit <- envfit(fin_func_cap_pars, fin_env_sel[, c(4, 8)], perm = 999, display = "lc", scaling = 1))
#scores(fit, "vectors", scaling = 1)
#plot(fit, p.max = 0.05, col = "blue", cex = 1.2, font = 2, lab = c("TOC", "Si"))
legend("bottomright", legend = c("Spring","Autumn"), bty = "n", cex = 1.5, pch = c(19, 17), col =  c("black", "black"))

par(op)
dev.off()

#### DESCRIPTIVE STATISTICS OF ENVIRONMENTAL VARIABLES -----------------
## Germany
ger_env
(ger_env_stat <- summaryFull(ger_env[, c(1, 2, 5, 15, 16)]))
(ger_env_stat_spr <- summaryFull(ger_env[1:8, c(1, 2, 5, 15, 16)]))
(ger_env_stat_aut <- summaryFull(ger_env[9:16, c(1, 2, 5, 15, 16)]))
(colnames(ger_env_describe <- rbind(t(ger_env_stat), t(ger_env_stat_spr), t(ger_env_stat_aut))))
(ger_env_describe <- ger_env_describe[, c(8, 9, 10, 2, 3, 13)])

## Sweden
swe_env
(swe_env_stat <- summaryFull(swe_env[, c(1, 2, 5, 18, 19, 21)]))
(swe_env_stat_spr <- summaryFull(swe_env[1:52, c(1, 2, 5, 18, 19, 21)]))
(swe_env_stat_aut <- summaryFull(swe_env[53:104, c(1, 2, 5, 18, 19, 21)]))
(colnames(swe_env_describe <- rbind(t(swe_env_stat), t(swe_env_stat_spr), t(swe_env_stat_aut))))
(swe_env_describe <- swe_env_describe[, c(7, 8, 9, 2, 3, 12)])

## Finland
fin_env
(fin_env_stat <- summaryFull(fin_env[, c(1, 2, 7, 9, 10)]))
(fin_env_stat_spr <- summaryFull(fin_env[1:10, c(1, 2, 7, 9, 10)]))
(fin_env_stat_aut <- summaryFull(fin_env[11:20, c(1, 2, 7, 9, 10)]))
(colnames(fin_env_describe <- rbind(t(fin_env_stat), t(fin_env_stat_spr), t(fin_env_stat_aut))))
(fin_env_describe <- fin_env_describe[, c(7, 8, 9, 2, 3, 12)])

#### VARIANCE PARTITIONING --------------------
### GERMANY --------------------
## ORGANIZE DATA
# Spatial data
ger_dist_pco <- as.data.frame(ger_dist_pco)  # pco of geographic distances
colnames(ger_dist_pco) = c("PC1", "PC2")
rownames(ger_dist_pco) = rownames(div[1:16, ])

# Environmental data
ger_env_sel$Season <- as.numeric(factor(ger_env_sel$Season)) # spring is 2, autumn is 1

## diversity data
ger_div_vp <- div_all[1:16, 1:11]

## Create variable groups
# Phenological
ger_env_sel$Season
# Geographical
ger_geographic <- cbind(ger_env_sel$Altitude, ger_dist_pco)
colnames(ger_geographic) <- c("Alt", "PC1", "PC2")
ger_geographic
# pH
ger_env_sel$pH
# Nutrients
(ger_nutrients <- ger_env_sel[, c("DOC", "SO4", "NH4.N")])

# Check env data for collinearity
ger_env_varpar <- cbind(ger_env_sel, ger_dist_pco)
ger_env_varpar <- ger_env_varpar[, c(6, 1, 7:8, 2:5)]
pairs(ger_env_varpar, lower.panel = panel.smooth, upper.panel = panel.cor, 
      diag.panel = panel.hist, main = "Pearson Correlation Matrix")
# Check for big correlations
ger_Corrs <- cor(ger_env_varpar)
ger_BigCorrs <- which(ger_Corrs > 0.7 & ger_Corrs < 1, arr.ind = TRUE)
# pairs(ger_env_sel[ , unique(rownames(ger_BigCorrs))], lower.panel = panel.smooth,upper.panel = panel.cor, diag.panel = panel.hist, main = "Pearson Correlation Matrix")
set.seed(1)
### ANALYZE DATA
## A
hist(log(ger_div_vp$A))
ger_varpart_A <- varpart(decostand(log(ger_div_vp$A), "standardize"), 
                         ger_env_sel$Season,
                         ger_geographic,
                         ger_env_sel$pH,
                         ger_nutrients)
plot(ger_varpart_A, bg = 2:5)
plot(ger_varpart_A, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
ger_A_varpart_fract <- ger_varpart_A$part$fract[1:4,]
ger_A_varpart_fract$Metric <- "A"
ger_A_varpart_fract
set.seed(1)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure seasonal
anova(rda(decostand(log(ger_div_vp$A), "standardize"), 
          ger_env_sel$Season, cbind(ger_geographic, ger_env_sel$pH, ger_nutrients)))
# Fraction [b], pure geography
anova(rda(decostand(log(ger_div_vp$A), "standardize"), 
          ger_geographic, cbind(ger_env_sel$Season, ger_env_sel$pH, ger_nutrients)))
# Fraction [c], pure acidity
anova(rda(decostand(log(ger_div_vp$A), "standardize"), 
          ger_env_sel$pH, cbind(ger_env_sel$Season, ger_geographic, ger_nutrients)))
# Fraction [d], pure nutrients
anova(rda(decostand(log(ger_div_vp$A), "standardize"), 
          ger_nutrients, cbind(ger_env_sel$Season, ger_geographic, ger_env_sel$pH)))

## N0
hist(ger_div_vp$N0)
ger_varpart_N0 <- varpart(decostand(ger_div_vp$N0, "standardize"), 
                          as.numeric(factor(ger_env_sel$Season)),
                          ger_geographic,
                          ger_env_sel$pH,
                          ger_nutrients)
plot(ger_varpart_N0, bg = 2:5)
plot(ger_varpart_N0, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
ger_N0_varpart_fract <- ger_varpart_N0$part$fract[1:4,]
ger_N0_varpart_fract$Metric <- "N0"
ger_N0_varpart_fract
set.seed(1)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure seasonal
anova(rda(decostand(ger_div_vp$N0, "standardize"), 
          ger_env_sel$Season, cbind(ger_geographic, ger_env_sel$pH, ger_nutrients)))
# Fraction [b], pure geographic
anova(rda(decostand(ger_div_vp$N0, "standardize"), 
          ger_geographic, cbind(ger_env_sel$Season, ger_env_sel$pH, ger_nutrients)))
# Fraction [c], pure acidity
anova(rda(decostand(ger_div_vp$N0, "standardize"), 
          ger_env_sel$pH, cbind(ger_env_sel$Season, ger_geographic, ger_nutrients)))
# Fraction [d], pure nutrients
anova(rda(decostand(ger_div_vp$N0, "standardize"), 
          ger_nutrients, cbind(ger_env_sel$Season, ger_geographic, ger_env_sel$pH)))

## E10
hist(ger_div_vp$E10)
ger_varpart_E10 <- varpart(decostand(ger_div_vp$E10, "standardize"), 
                           as.numeric(factor(ger_env_sel$Season)),
                           ger_geographic,
                           ger_env_sel$pH,
                           ger_nutrients)
plot(ger_varpart_E10, bg = 2:5)
plot(ger_varpart_E10, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
ger_E10_varpart_fract <- ger_varpart_E10$part$fract[1:4,]
ger_E10_varpart_fract$Metric <- "E10"
ger_E10_varpart_fract
set.seed(1)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure season
anova(rda(decostand(ger_div_vp$E10, "standardize"), 
          ger_env_sel$Season, cbind(ger_geographic, ger_env_sel$pH, ger_nutrients)))
# Fraction [b], pure geography
anova(rda(decostand(ger_div_vp$E10, "standardize"), 
          ger_geographic, cbind(ger_env_sel$Season, ger_env_sel$pH, ger_nutrients)))
# Fraction [c], pure acidity
anova(rda(decostand(ger_div_vp$E10, "standardize"), 
          ger_env_sel$pH, cbind(ger_env_sel$Season, ger_geographic, ger_nutrients)))
# Fraction [d], pure mutrients
anova(rda(decostand(ger_div_vp$E10, "standardize"), 
          ger_nutrients, cbind(ger_env_sel$Season, ger_geographic, ger_env_sel$pH)))

## N1
hist(ger_div_vp$N1)
ger_varpart_N1 <- varpart(decostand(ger_div_vp$N1, "standardize"), 
                          as.numeric(factor(ger_env_sel$Season)),
                          ger_geographic,
                          ger_env_sel$pH,
                          ger_nutrients)
plot(ger_varpart_N1, bg = 2:5)
plot(ger_varpart_N1, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
ger_N1_varpart_fract <- ger_varpart_N1$part$fract[1:4,]
ger_N1_varpart_fract$Metric <- "N1"
ger_N1_varpart_fract
set.seed(1)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(decostand(ger_div_vp$N1, "standardize"), 
          ger_env_sel$Season, cbind(ger_geographic, ger_env_sel$pH, ger_nutrients)))
# Fraction [b], pure trend
anova(rda(decostand(ger_div_vp$N1, "standardize"), 
          ger_geographic, cbind(ger_env_sel$Season, ger_env_sel$pH, ger_nutrients)))
# Fraction [c], pure broad scale spatial
anova(rda(decostand(ger_div_vp$N1, "standardize"), 
          ger_env_sel$pH, cbind(ger_env_sel$Season, ger_geographic, ger_nutrients)))
# Fraction [d], pure fine scale spatial
anova(rda(decostand(ger_div_vp$N1, "standardize"), 
          ger_nutrients, cbind(ger_env_sel$Season, ger_geographic, ger_env_sel$pH)))

## FRic
hist(log(ger_div_vp$FRic))
ger_varpart_FRic <- varpart(decostand(log(ger_div_vp$FRic), "standardize"), 
                            as.numeric(factor(ger_env_sel$Season)),
                            ger_geographic,
                            ger_env_sel$pH,
                            ger_nutrients)
plot(ger_varpart_FRic, bg = 2:5)
plot(ger_varpart_FRic, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
ger_FRic_varpart_fract <- ger_varpart_FRic$part$fract[1:4,]
ger_FRic_varpart_fract$Metric <- "FRic"
ger_FRic_varpart_fract
set.seed(1)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(decostand(ger_div_vp$FRic, "standardize"), 
          ger_env_sel$Season, cbind(ger_geographic, ger_env_sel$pH, ger_nutrients)))
# Fraction [b], pure trend
anova(rda(decostand(ger_div_vp$FRic, "standardize"), 
          ger_geographic, cbind(ger_env_sel$Season, ger_env_sel$pH, ger_nutrients)))
# Fraction [c], pure broad scale spatial
anova(rda(decostand(ger_div_vp$FRic, "standardize"), 
          ger_env_sel$pH, cbind(ger_env_sel$Season, ger_geographic, ger_nutrients)))
# Fraction [d], pure fine scale spatial
anova(rda(decostand(ger_div_vp$FRic, "standardize"), 
          ger_nutrients, cbind(ger_env_sel$Season, ger_geographic, ger_env_sel$pH)))

## FEve
hist(ger_div_vp$FEve)
ger_varpart_FEve <- varpart(decostand(ger_div_vp$FEve, "standardize"), 
                            as.numeric(factor(ger_env_sel$Season)),
                            ger_geographic,
                            ger_env_sel$pH,
                            ger_nutrients)
plot(ger_varpart_FEve, bg = 2:5)
plot(ger_varpart_FEve, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
ger_FEve_varpart_fract <- ger_varpart_FEve$part$fract[1:4,]
ger_FEve_varpart_fract$Metric <- "FEve"
ger_FEve_varpart_fract
set.seed(1)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(decostand(ger_div_vp$FEve, "standardize"), 
          ger_env_sel$Season, cbind(ger_geographic, ger_env_sel$pH, ger_nutrients)))
# Fraction [b], pure trend
anova(rda(decostand(ger_div_vp$FEve, "standardize"), 
          ger_geographic, cbind(ger_env_sel$Season, ger_env_sel$pH, ger_nutrients)))
# Fraction [c], pure broad scale spatial
anova(rda(decostand(ger_div_vp$FEve, "standardize"), 
          ger_env_sel$pH, cbind(ger_env_sel$Season, ger_geographic, ger_nutrients)))
# Fraction [d], pure fine scale spatial
anova(rda(decostand(ger_div_vp$FEve, "standardize"), 
          ger_nutrients, cbind(ger_env_sel$Season, ger_geographic, ger_env_sel$pH)))

## FDis
hist(ger_div_vp$FDis)
ger_varpart_FDis <- varpart(decostand(ger_div_vp$FDis, "standardize"), 
                            as.numeric(factor(ger_env_sel$Season)),
                            ger_geographic,
                            ger_env_sel$pH,
                            ger_nutrients)
plot(ger_varpart_FDis)
plot(ger_varpart_FDis, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
ger_FDis_varpart_fract <- ger_varpart_FDis$part$fract[1:4,]
ger_FDis_varpart_fract$Metric <- "FDis"
ger_FDis_varpart_fract
set.seed(1)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(decostand(ger_div_vp$FDis, "standardize"), 
          ger_env_sel$Season, cbind(ger_geographic, ger_env_sel$pH, ger_nutrients)))
# Fraction [b], pure trend
anova(rda(decostand(ger_div_vp$FDis, "standardize"), 
          ger_geographic, cbind(ger_env_sel$Season, ger_env_sel$pH, ger_nutrients)))
# Fraction [c], pure broad scale spatial
anova(rda(decostand(ger_div_vp$FDis, "standardize"), 
          ger_env_sel$pH, cbind(ger_env_sel$Season, ger_geographic, ger_nutrients)))
# Fraction [d], pure fine scale spatial
anova(rda(decostand(ger_div_vp$FDis, "standardize"), 
          ger_nutrients, cbind(ger_env_sel$Season, ger_geographic, ger_env_sel$pH)))

## FRic.SES
hist(ger_div_vp$FRic.SES)
ger_varpart_FRic.SES <- varpart(decostand(ger_div_vp$FRic.SES, "standardize"), 
                                as.numeric(factor(ger_env_sel$Season)),
                                ger_geographic,
                                ger_env_sel$pH,
                                ger_nutrients)
plot(ger_varpart_FRic.SES, bg = 2:5)
plot(ger_varpart_FRic.SES, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
ger_FRic.SES_varpart_fract <- ger_varpart_FRic.SES$part$fract[1:4,]
ger_FRic.SES_varpart_fract$Metric <- "FRic.SES"
ger_FRic.SES_varpart_fract
set.seed(1)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(decostand(ger_div_vp$FRic.SES, "standardize"), 
          ger_env_sel$Season, cbind(ger_geographic, ger_env_sel$pH, ger_nutrients)))
# Fraction [b], pure trend
anova(rda(decostand(ger_div_vp$FRic.SES, "standardize"), 
          ger_geographic, cbind(ger_env_sel$Season, ger_env_sel$pH, ger_nutrients)))
# Fraction [c], pure broad scale spatial
anova(rda(decostand(ger_div_vp$FRic.SES, "standardize"), 
          ger_env_sel$pH, cbind(ger_env_sel$Season, ger_geographic, ger_nutrients)))
# Fraction [d], pure fine scale spatial
anova(rda(decostand(ger_div_vp$FRic.SES, "standardize"), 
          ger_nutrients, cbind(ger_env_sel$Season, ger_geographic, ger_env_sel$pH)))

## FEve.SES
hist(ger_div_vp$FEve.SES)
ger_varpart_FEve.SES <- varpart(decostand(ger_div_vp$FEve.SES, "standardize"), 
                                as.numeric(factor(ger_env_sel$Season)),
                                ger_geographic,
                                ger_env_sel$pH,
                                ger_nutrients)
plot(ger_varpart_FEve.SES, bg = 2:5)
plot(ger_varpart_FEve.SES, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
ger_FEve.SES_varpart_fract <- ger_varpart_FEve.SES$part$fract[1:4,]
ger_FEve.SES_varpart_fract$Metric <- "FEve.SES"
ger_FEve.SES_varpart_fract
set.seed(1)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(decostand(ger_div_vp$FEve.SES, "standardize"), 
          ger_env_sel$Season, cbind(ger_geographic, ger_env_sel$pH, ger_nutrients)))
# Fraction [b], pure trend
anova(rda(decostand(ger_div_vp$FEve.SES, "standardize"), 
          ger_geographic, cbind(ger_env_sel$Season, ger_env_sel$pH, ger_nutrients)))
# Fraction [c], pure broad scale spatial
anova(rda(decostand(ger_div_vp$FEve.SES, "standardize"), 
          ger_env_sel$pH, cbind(ger_env_sel$Season, ger_geographic, ger_nutrients)))
# Fraction [d], pure fine scale spatial
anova(rda(decostand(ger_div_vp$FEve.SES, "standardize"), 
          ger_nutrients, cbind(ger_env_sel$Season, ger_geographic, ger_env_sel$pH)))

## FDis.SES
hist(ger_div_vp$FDis.SES)
ger_varpart_FDis.SES <- varpart(decostand(ger_div_vp$FDis.SES, "standardize"), 
                                as.numeric(factor(ger_env_sel$Season)),
                                ger_geographic,
                                ger_env_sel$pH,
                                ger_nutrients)
plot(ger_varpart_FDis.SES, bg = 2:5)
plot(ger_varpart_FDis.SES, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
ger_FDis.SES_varpart_fract <- ger_varpart_FDis.SES$part$fract[1:4,]
ger_FDis.SES_varpart_fract$Metric <- "FDis.SES"
ger_FDis.SES_varpart_fract
set.seed(1)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(decostand(ger_div_vp$FDis.SES, "standardize"), 
          ger_env_sel$Season, cbind(ger_geographic, ger_env_sel$pH, ger_nutrients)))
# Fraction [b], pure trend
anova(rda(decostand(ger_div_vp$FDis.SES, "standardize"), 
          ger_geographic, cbind(ger_env_sel$Season, ger_env_sel$pH, ger_nutrients)))
# Fraction [c], pure broad scale spatial
anova(rda(decostand(ger_div_vp$FDis.SES, "standardize"), 
          ger_env_sel$pH, cbind(ger_env_sel$Season, ger_geographic, ger_nutrients)))
# Fraction [d], pure fine scale spatial
anova(rda(decostand(ger_div_vp$FDis.SES, "standardize"), 
          ger_nutrients, cbind(ger_env_sel$Season, ger_geographic, ger_env_sel$pH)))

## FRed
hist(ger_div_vp$FRed)
ger_varpart_FRed <- varpart(decostand(ger_div_vp$FRed, "standardize"), 
                            as.numeric(factor(ger_env_sel$Season)),
                            ger_geographic,
                            ger_env_sel$pH,
                            ger_nutrients)
plot(ger_varpart_FRed, bg = 2:5)
plot(ger_varpart_FRed, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
ger_FRed_varpart_fract <- ger_varpart_FRed$part$fract[1:4,]
ger_FRed_varpart_fract$Metric <- "FRed"
ger_FRed_varpart_fract
set.seed(1)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(decostand(ger_div_vp$FRed, "standardize"), 
          ger_env_sel$Season, cbind(ger_geographic, ger_env_sel$pH, ger_nutrients)))
# Fraction [b], pure trend
anova(rda(decostand(ger_div_vp$FRed, "standardize"), 
          ger_geographic, cbind(ger_env_sel$Season, ger_env_sel$pH, ger_nutrients)))
# Fraction [c], pure broad scale spatial
anova(rda(decostand(ger_div_vp$FRed, "standardize"), 
          ger_env_sel$pH, cbind(ger_env_sel$Season, ger_geographic, ger_nutrients)))
# Fraction [d], pure fine scale spatial
anova(rda(decostand(ger_div_vp$FRed, "standardize"), 
          ger_nutrients, cbind(ger_env_sel$Season, ger_geographic, ger_env_sel$pH)))

ger_varpart_fract <- rbind(ger_A_varpart_fract, ger_N0_varpart_fract, ger_E10_varpart_fract, ger_N1_varpart_fract, 
                           ger_FRic_varpart_fract, ger_FEve_varpart_fract, ger_FDis_varpart_fract,
                           ger_FRic.SES_varpart_fract, ger_FEve.SES_varpart_fract, ger_FDis.SES_varpart_fract,
                           ger_FRed_varpart_fract)
(ger_varpart_fract <- ger_varpart_fract[, c(3, 5)])
rownames(ger_varpart_fract) <- c("Season", "Geography", "Acidity", "Nutrients", 
                                 "Season0", "Geography0", "Acidity0", "Nutrients0", 
                                 "Season1", "Geography1", "Acidity1", "Nutrients1",
                                 "Season2", "Geography2", "Acidity2", "Nutrients2", 
                                 "Season3", "Geography3", "Acidity3", "Nutrients3", 
                                 "Season4", "Geography4", "Acidity4", "Nutrients4", 
                                 "Season5", "Geography5", "Acidity5", "Nutrients5",
                                 "Season6", "Geography6", "Acidity6", "Nutrients6", 
                                 "Season7", "Geography7", "Acidity7", "Nutrients7", 
                                 "Season8", "Geography8", "Acidity8", "Nutrients8", 
                                 "Season9", "Geography9", "Acidity9", "Nutrients9")

### SWEDEN --------------------
## ORGANIZE DATA
# Spatial data
swe_dist_pco <- as.data.frame(swe_dist_pco)  # pco of geographic distances
colnames(swe_dist_pco) = c("PC1", "PC2")
rownames(swe_dist_pco) = rownames(div[17:120, ])

# Environmental data
swe_env_sel$Season <- as.numeric(factor(swe_env_sel$Season)) # spring is 2, autumn is 1

## diversity data
swe_div_vp <- div_all[17:120, 1:11]

## Create variable groups
# Phenological
swe_env_sel$Season
# Geographical
swe_geographic <- cbind(swe_env_sel$Altitude, swe_dist_pco)
colnames(swe_geographic) <- c("Alt", "PC1", "PC2")
swe_geographic
# pH
swe_env_sel$pH
# Nutrients
(swe_nutrients <- swe_env_sel[, c("TOC", "SO4", "NO2.NO3_N", "Tot.P")])

### ANALYZE DATA
## A
hist(log(swe_div_vp$A))
swe_varpart_A <- varpart(decostand(log(swe_div_vp$A), "standardize"), 
                         swe_env_sel$Season,
                         swe_geographic,
                         swe_env_sel$pH,
                         swe_nutrients)
plot(swe_varpart_A, bg = 2:5)
plot(swe_varpart_A, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
swe_A_varpart_fract <- swe_varpart_A$part$fract[1:4,]
swe_A_varpart_fract$Metric <- "A"
swe_A_varpart_fract
set.seed(1)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(decostand(log(swe_div_vp$A), "standardize"), 
          swe_env_sel$Season, cbind(swe_geographic, swe_env_sel$pH, swe_nutrients)))
# Fraction [b], pure trend
anova(rda(decostand(log(swe_div_vp$A), "standardize"), 
          swe_geographic, cbind(swe_env_sel$Season, swe_env_sel$pH, swe_nutrients)))
# Fraction [c], pure broad scale spatial
anova(rda(decostand(log(swe_div_vp$A), "standardize"), 
          swe_env_sel$pH, cbind(swe_env_sel$Season, swe_geographic, swe_nutrients)))
# Fraction [d], pure fine scale spatial
anova(rda(decostand(log(swe_div_vp$A), "standardize"), 
          swe_nutrients, cbind(swe_env_sel$Season, swe_geographic, swe_env_sel$pH)))

## N0
hist(swe_div_vp$N0)
swe_varpart_N0 <- varpart(decostand(swe_div_vp$N0, "standardize"), 
                          as.numeric(factor(swe_env_sel$Season)),
                          swe_geographic,
                          swe_env_sel$pH,
                          swe_nutrients)
plot(swe_varpart_N0, bg = 2:5)
plot(swe_varpart_N0, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
swe_N0_varpart_fract <- swe_varpart_N0$part$fract[1:4,]
swe_N0_varpart_fract$Metric <- "N0"
swe_N0_varpart_fract
set.seed(1)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(decostand(swe_div_vp$N0, "standardize"), 
          swe_env_sel$Season, cbind(swe_geographic, swe_env_sel$pH, swe_nutrients)))
# Fraction [b], pure trend
anova(rda(decostand(swe_div_vp$N0, "standardize"), 
          swe_geographic, cbind(swe_env_sel$Season, swe_env_sel$pH, swe_nutrients)))
# Fraction [c], pure broad scale spatial
anova(rda(decostand(swe_div_vp$N0, "standardize"), 
          swe_env_sel$pH, cbind(swe_env_sel$Season, swe_geographic, swe_nutrients)))
# Fraction [d], pure fine scale spatial
anova(rda(decostand(swe_div_vp$N0, "standardize"), 
          swe_nutrients, cbind(swe_env_sel$Season, swe_geographic, swe_env_sel$pH)))

## E10
hist(swe_div_vp$E10)
swe_varpart_E10 <- varpart(decostand(swe_div_vp$E10, "standardize"), 
                           as.numeric(factor(swe_env_sel$Season)),
                           swe_geographic,
                           swe_env_sel$pH,
                           swe_nutrients)
plot(swe_varpart_E10, bg = 2:5)
plot(swe_varpart_E10, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
swe_E10_varpart_fract <- swe_varpart_E10$part$fract[1:4,]
swe_E10_varpart_fract$Metric <- "E10"
swe_E10_varpart_fract
set.seed(1)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(decostand(swe_div_vp$E10, "standardize"), 
          swe_env_sel$Season, cbind(swe_geographic, swe_env_sel$pH, swe_nutrients)))
# Fraction [b], pure trend
anova(rda(decostand(swe_div_vp$E10, "standardize"), 
          swe_geographic, cbind(swe_env_sel$Season, swe_env_sel$pH, swe_nutrients)))
# Fraction [c], pure broad scale spatial
anova(rda(decostand(swe_div_vp$E10, "standardize"), 
          swe_env_sel$pH, cbind(swe_env_sel$Season, swe_geographic, swe_nutrients)))
# Fraction [d], pure fine scale spatial
anova(rda(decostand(swe_div_vp$E10, "standardize"), 
          swe_nutrients, cbind(swe_env_sel$Season, swe_geographic, swe_env_sel$pH)))

## N1
hist(swe_div_vp$N1)
swe_varpart_N1 <- varpart(decostand(swe_div_vp$N1, "standardize"), 
                          as.numeric(factor(swe_env_sel$Season)),
                          swe_geographic,
                          swe_env_sel$pH,
                          swe_nutrients)
plot(swe_varpart_N1, bg = 2:5)
plot(swe_varpart_N1, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
swe_N1_varpart_fract <- swe_varpart_N1$part$fract[1:4,]
swe_N1_varpart_fract$Metric <- "N1"
swe_N1_varpart_fract
set.seed(1)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(decostand(swe_div_vp$N1, "standardize"), 
          swe_env_sel$Season, cbind(swe_geographic, swe_env_sel$pH, swe_nutrients)))
# Fraction [b], pure trend
anova(rda(decostand(swe_div_vp$N1, "standardize"), 
          swe_geographic, cbind(swe_env_sel$Season, swe_env_sel$pH, swe_nutrients)))
# Fraction [c], pure broad scale spatial
anova(rda(decostand(swe_div_vp$N1, "standardize"), 
          swe_env_sel$pH, cbind(swe_env_sel$Season, swe_geographic, swe_nutrients)))
# Fraction [d], pure fine scale spatial
anova(rda(decostand(swe_div_vp$N1, "standardize"), 
          swe_nutrients, cbind(swe_env_sel$Season, swe_geographic, swe_env_sel$pH)))

## FRic
hist(log(swe_div_vp$FRic))
swe_varpart_FRic <- varpart(decostand(swe_div_vp$FRic, "standardize"), 
                            as.numeric(factor(swe_env_sel$Season)),
                            swe_geographic,
                            swe_env_sel$pH,
                            swe_nutrients)
plot(swe_varpart_FRic, bg = 2:5)
plot(swe_varpart_FRic, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
swe_FRic_varpart_fract <- swe_varpart_FRic$part$fract[1:4,]
swe_FRic_varpart_fract$Metric <- "FRic"
swe_FRic_varpart_fract
set.seed(1)
## Test fraction [a, b, c, d]
anova(swe_rda_FRic_a <- rda(decostand(swe_div_vp$FRic, "standardize") ~ swe_env_sel$Season +
                              Condition(swe_geographic$Alt + swe_geographic$PC1 + swe_geographic$PC2) +
                              Condition(swe_env_sel$pH) +
                              Condition(swe_nutrients$TOC + swe_nutrients$SO4 + swe_nutrients$NO2.NO3_N + swe_nutrients$Tot.P)))
anova(swe_rda_FRic_b <- rda(decostand(swe_div_vp$FRic, "standardize") ~ swe_geographic$Alt + swe_geographic$PC1 + swe_geographic$PC2 +
                              Condition(swe_env_sel$Season) +
                              Condition(swe_env_sel$pH) +
                              Condition(swe_nutrients$TOC + swe_nutrients$SO4 + swe_nutrients$NO2.NO3_N + swe_nutrients$Tot.P)))
anova(swe_rda_FRic_c <- rda(decostand(swe_div_vp$FRic, "standardize") ~ swe_env_sel$pH +
                              Condition(swe_env_sel$Season) +
                              Condition(swe_geographic$Alt + swe_geographic$PC1 + swe_geographic$PC2) +
                              Condition(swe_nutrients$TOC + swe_nutrients$SO4 + swe_nutrients$NO2.NO3_N + swe_nutrients$Tot.P)))
anova(swe_rda_FRic_d <- rda(decostand(swe_div_vp$FRic, "standardize") ~ swe_nutrients$TOC + swe_nutrients$SO4 + swe_nutrients$NO2.NO3_N + swe_nutrients$Tot.P +
                              Condition(swe_env_sel$Season) +
                              Condition(swe_geographic$Alt + swe_geographic$PC1 + swe_geographic$PC2) +
                              Condition(swe_env_sel$pH)))

## FEve
hist(swe_div_vp$FEve)
swe_varpart_FEve <- varpart(decostand(swe_div_vp$FEve, "standardize"), 
                            as.numeric(factor(swe_env_sel$Season)),
                            swe_geographic,
                            swe_env_sel$pH,
                            swe_nutrients)
plot(swe_varpart_FEve, bg = 2:5)
plot(swe_varpart_FEve, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
swe_FEve_varpart_fract <- swe_varpart_FEve$part$fract[1:4,]
swe_FEve_varpart_fract$Metric <- "FEve"
swe_FEve_varpart_fract
set.seed(1)
## Test fraction [a, b, c, d]
anova(swe_rda_FEve_a <- rda(decostand(swe_div_vp$FEve, "standardize") ~ swe_env_sel$Season +
                              Condition(swe_geographic$Alt + swe_geographic$PC1 + swe_geographic$PC2) +
                              Condition(swe_env_sel$pH) +
                              Condition(swe_nutrients$TOC + swe_nutrients$SO4 + swe_nutrients$NO2.NO3_N + swe_nutrients$Tot.P)))
anova(swe_rda_FEve_b <- rda(decostand(swe_div_vp$FEve, "standardize") ~ swe_geographic$Alt + swe_geographic$PC1 + swe_geographic$PC2 +
                              Condition(swe_env_sel$Season) +
                              Condition(swe_env_sel$pH) +
                              Condition(swe_nutrients$TOC + swe_nutrients$SO4 + swe_nutrients$NO2.NO3_N + swe_nutrients$Tot.P)))
anova(swe_rda_FEve_c <- rda(decostand(swe_div_vp$FEve, "standardize") ~ swe_env_sel$pH +
                              Condition(swe_env_sel$Season) +
                              Condition(swe_geographic$Alt + swe_geographic$PC1 + swe_geographic$PC2) +
                              Condition(swe_nutrients$TOC + swe_nutrients$SO4 + swe_nutrients$NO2.NO3_N + swe_nutrients$Tot.P)))
anova(swe_rda_FEve_d <- rda(decostand(swe_div_vp$FEve, "standardize") ~ swe_nutrients$TOC + swe_nutrients$SO4 + swe_nutrients$NO2.NO3_N + swe_nutrients$Tot.P +
                              Condition(swe_env_sel$Season) +
                              Condition(swe_geographic$Alt + swe_geographic$PC1 + swe_geographic$PC2) +
                              Condition(swe_env_sel$pH)))

## FDis
hist(swe_div_vp$FDis)
swe_varpart_FDis <- varpart(decostand(swe_div_vp$FDis, "standardize"), 
                            as.numeric(factor(swe_env_sel$Season)),
                            swe_geographic,
                            swe_env_sel$pH,
                            swe_nutrients)
plot(swe_varpart_FDis)
plot(swe_varpart_FDis, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
swe_FDis_varpart_fract <- swe_varpart_FDis$part$fract[1:4,]
swe_FDis_varpart_fract$Metric <- "FDis"
swe_FDis_varpart_fract
set.seed(1)
## Test fraction [a, b, c, d]
anova(swe_rda_FDis_a <- rda(decostand(swe_div_vp$FDis, "standardize") ~ swe_env_sel$Season +
                              Condition(swe_geographic$Alt + swe_geographic$PC1 + swe_geographic$PC2) +
                              Condition(swe_env_sel$pH) +
                              Condition(swe_nutrients$TOC + swe_nutrients$SO4 + swe_nutrients$NO2.NO3_N + swe_nutrients$Tot.P)))
anova(swe_rda_FDis_b <- rda(decostand(swe_div_vp$FDis, "standardize") ~ swe_geographic$Alt + swe_geographic$PC1 + swe_geographic$PC2 +
                              Condition(swe_env_sel$Season) +
                              Condition(swe_env_sel$pH) +
                              Condition(swe_nutrients$TOC + swe_nutrients$SO4 + swe_nutrients$NO2.NO3_N + swe_nutrients$Tot.P)))
anova(swe_rda_FDis_c <- rda(decostand(swe_div_vp$FDis, "standardize") ~ swe_env_sel$pH +
                              Condition(swe_env_sel$Season) +
                              Condition(swe_geographic$Alt + swe_geographic$PC1 + swe_geographic$PC2) +
                              Condition(swe_nutrients$TOC + swe_nutrients$SO4 + swe_nutrients$NO2.NO3_N + swe_nutrients$Tot.P)))
anova(swe_rda_FDis_d <- rda(decostand(swe_div_vp$FDis, "standardize") ~ swe_nutrients$TOC + swe_nutrients$SO4 + swe_nutrients$NO2.NO3_N + swe_nutrients$Tot.P +
                              Condition(swe_env_sel$Season) +
                              Condition(swe_geographic$Alt + swe_geographic$PC1 + swe_geographic$PC2) +
                              Condition(swe_env_sel$pH)))

## FRic.SES
hist(swe_div_vp$FRic.SES)
swe_varpart_FRic.SES <- varpart(decostand(swe_div_vp$FRic.SES, "standardize"), 
                                as.numeric(factor(swe_env_sel$Season)),
                                swe_geographic,
                                swe_env_sel$pH,
                                swe_nutrients)
plot(swe_varpart_FRic.SES, bg = 2:5)
plot(swe_varpart_FRic.SES, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
swe_FRic.SES_varpart_fract <- swe_varpart_FRic.SES$part$fract[1:4,]
swe_FRic.SES_varpart_fract$Metric <- "FRic.SES"
swe_FRic.SES_varpart_fract
set.seed(1)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(decostand(swe_div_vp$FRic.SES, "standardize"), 
          swe_env_sel$Season, cbind(swe_geographic, swe_env_sel$pH, swe_nutrients)))
# Fraction [b], pure trend
anova(rda(decostand(swe_div_vp$FRic.SES, "standardize"), 
          swe_geographic, cbind(swe_env_sel$Season, swe_env_sel$pH, swe_nutrients)))
# Fraction [c], pure broad scale spatial
anova(rda(decostand(swe_div_vp$FRic.SES, "standardize"), 
          swe_env_sel$pH, cbind(swe_env_sel$Season, swe_geographic, swe_nutrients)))
# Fraction [d], pure fine scale spatial
anova(rda(decostand(swe_div_vp$FRic.SES, "standardize"), 
          swe_nutrients, cbind(swe_env_sel$Season, swe_geographic, swe_env_sel$pH)))

## FEve.SES
hist(swe_div_vp$FEve.SES)
swe_varpart_FEve.SES <- varpart(decostand(swe_div_vp$FEve.SES, "standardize"), 
                                as.numeric(factor(swe_env_sel$Season)),
                                swe_geographic,
                                swe_env_sel$pH,
                                swe_nutrients)
plot(swe_varpart_FEve.SES, bg = 2:5)
plot(swe_varpart_FEve.SES, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
swe_FEve.SES_varpart_fract <- swe_varpart_FEve.SES$part$fract[1:4,]
swe_FEve.SES_varpart_fract$Metric <- "FEve.SES"
swe_FEve.SES_varpart_fract
set.seed(1)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(decostand(swe_div_vp$FEve.SES, "standardize"), 
          swe_env_sel$Season, cbind(swe_geographic, swe_env_sel$pH, swe_nutrients)))
# Fraction [b], pure trend
anova(rda(decostand(swe_div_vp$FEve.SES, "standardize"), 
          swe_geographic, cbind(swe_env_sel$Season, swe_env_sel$pH, swe_nutrients)))
# Fraction [c], pure broad scale spatial
anova(rda(decostand(swe_div_vp$FEve.SES, "standardize"), 
          swe_env_sel$pH, cbind(swe_env_sel$Season, swe_geographic, swe_nutrients)))
# Fraction [d], pure fine scale spatial
anova(rda(decostand(swe_div_vp$FEve.SES, "standardize"), 
          swe_nutrients, cbind(swe_env_sel$Season, swe_geographic, swe_env_sel$pH)))

## FDis.SES
hist(swe_div_vp$FDis.SES)
swe_varpart_FDis.SES <- varpart(decostand(swe_div_vp$FDis.SES, "standardize"), 
                                as.numeric(factor(swe_env_sel$Season)),
                                swe_geographic,
                                swe_env_sel$pH,
                                swe_nutrients)
plot(swe_varpart_FDis.SES, bg = 2:5)
plot(swe_varpart_FDis.SES, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
swe_FDis.SES_varpart_fract <- swe_varpart_FDis.SES$part$fract[1:4,]
swe_FDis.SES_varpart_fract$Metric <- "FDis.SES"
swe_FDis.SES_varpart_fract
set.seed(1)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(decostand(swe_div_vp$FDis.SES, "standardize"), 
          swe_env_sel$Season, cbind(swe_geographic, swe_env_sel$pH, swe_nutrients)))
# Fraction [b], pure trend
anova(rda(decostand(swe_div_vp$FDis.SES, "standardize"), 
          swe_geographic, cbind(swe_env_sel$Season, swe_env_sel$pH, swe_nutrients)))
# Fraction [c], pure broad scale spatial
anova(rda(decostand(swe_div_vp$FDis.SES, "standardize"), 
          swe_env_sel$pH, cbind(swe_env_sel$Season, swe_geographic, swe_nutrients)))
# Fraction [d], pure fine scale spatial
anova(rda(decostand(swe_div_vp$FDis.SES, "standardize"), 
          swe_nutrients, cbind(swe_env_sel$Season, swe_geographic, swe_env_sel$pH)))

## FRed
hist(swe_div_vp$FRed)
swe_varpart_FRed <- varpart(decostand(swe_div_vp$FRed, "standardize"), 
                            as.numeric(factor(swe_env_sel$Season)),
                            swe_geographic,
                            swe_env_sel$pH,
                            swe_nutrients)
plot(swe_varpart_FRed, bg = 2:5)
plot(swe_varpart_FRed, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
swe_FRed_varpart_fract <- swe_varpart_FRed$part$fract[1:4,]
swe_FRed_varpart_fract$Metric <- "FRed"
swe_FRed_varpart_fract
set.seed(1)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(decostand(swe_div_vp$FRed, "standardize"), 
          swe_env_sel$Season, cbind(swe_geographic, swe_env_sel$pH, swe_nutrients)))
# Fraction [b], pure trend
anova(rda(decostand(swe_div_vp$FRed, "standardize"), 
          swe_geographic, cbind(swe_env_sel$Season, swe_env_sel$pH, swe_nutrients)))
# Fraction [c], pure broad scale spatial
anova(rda(decostand(swe_div_vp$FRed, "standardize"), 
          swe_env_sel$pH, cbind(swe_env_sel$Season, swe_geographic, swe_nutrients)))
# Fraction [d], pure fine scale spatial
anova(rda(decostand(swe_div_vp$FRed, "standardize"), 
          swe_nutrients, cbind(swe_env_sel$Season, swe_geographic, swe_env_sel$pH)))

swe_varpart_fract <- rbind(swe_A_varpart_fract, swe_N0_varpart_fract, swe_E10_varpart_fract, swe_N1_varpart_fract, 
                           swe_FRic_varpart_fract, swe_FEve_varpart_fract, swe_FDis_varpart_fract,
                           swe_FRic.SES_varpart_fract, swe_FEve.SES_varpart_fract, swe_FDis.SES_varpart_fract,
                           swe_FRed_varpart_fract)
(swe_varpart_fract <- swe_varpart_fract[, c(3, 5)])
rownames(swe_varpart_fract) <- c("Season", "Geography", "Acidity", "Nutrients", 
                                 "Season0", "Geography0", "Acidity0", "Nutrients0", 
                                 "Season1", "Geography1", "Acidity1", "Nutrients1",
                                 "Season2", "Geography2", "Acidity2", "Nutrients2", 
                                 "Season3", "Geography3", "Acidity3", "Nutrients3", 
                                 "Season4", "Geography4", "Acidity4", "Nutrients4", 
                                 "Season5", "Geography5", "Acidity5", "Nutrients5",
                                 "Season6", "Geography6", "Acidity6", "Nutrients6", 
                                 "Season7", "Geography7", "Acidity7", "Nutrients7", 
                                 "Season8", "Geography8", "Acidity8", "Nutrients8", 
                                 "Season9", "Geography9", "Acidity9", "Nutrients9")

### FINLAND --------------------
## ORGANIZE DATA
# Spatial data
fin_dist_pco <- as.data.frame(fin_dist_pco)  # pco of geographic distances
colnames(fin_dist_pco) = c("PC1", "PC2")
rownames(fin_dist_pco) = rownames(div[121:140, ])

# Environmental data
fin_env_sel$Season <- as.numeric(factor(fin_env_sel$Season)) # spring is 2, autumn is 1
fin_env_sel

## diversity data
fin_div_vp <- div_all[121:140, 1:11]

## Create variable groups
# Phenological
fin_env_sel$Season
# Geographical
fin_geographic <- cbind(fin_env_sel$Altitude, fin_dist_pco)
colnames(fin_geographic) <- c("Alt", "PC1", "PC2")
fin_geographic
# pH
fin_env_sel$pH
# Nutrients
(fin_nutrients <- fin_env_sel[, c("DOC", "Tot.N", "Tot.P")])

### ANALYZE DATA
## A
hist(log(fin_div_vp$A))
fin_varpart_A <- varpart(decostand(log(fin_div_vp$A), "standardize"), 
                         fin_env_sel$Season,
                         fin_geographic,
                         fin_env_sel$pH,
                         fin_nutrients)
plot(fin_varpart_A, bg = 2:5)
plot(fin_varpart_A, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
fin_A_varpart_fract <- fin_varpart_A$part$fract[1:4,]
fin_A_varpart_fract$Metric <- "A"
fin_A_varpart_fract
set.seed(1)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(decostand(log(fin_div_vp$A), "standardize"), 
          fin_env_sel$Season, cbind(fin_geographic, fin_env_sel$pH, fin_nutrients)))
# Fraction [b], pure trend
anova(rda(decostand(log(fin_div_vp$A), "standardize"), 
          fin_geographic, cbind(fin_env_sel$Season, fin_env_sel$pH, fin_nutrients)))
# Fraction [c], pure broad scale spatial
anova(rda(decostand(log(fin_div_vp$A), "standardize"), 
          fin_env_sel$pH, cbind(fin_env_sel$Season, fin_geographic, fin_nutrients)))
# Fraction [d], pure fine scale spatial
anova(rda(decostand(log(fin_div_vp$A), "standardize"), 
          fin_nutrients, cbind(fin_env_sel$Season, fin_geographic, fin_env_sel$pH)))

## N0
hist(fin_div_vp$N0)
fin_varpart_N0 <- varpart(decostand(fin_div_vp$N0, "standardize"), 
                          as.numeric(factor(fin_env_sel$Season)),
                          fin_geographic,
                          fin_env_sel$pH,
                          fin_nutrients)
plot(fin_varpart_N0, bg = 2:5)
plot(fin_varpart_N0, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
fin_N0_varpart_fract <- fin_varpart_N0$part$fract[1:4,]
fin_N0_varpart_fract$Metric <- "N0"
fin_N0_varpart_fract
set.seed(1)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(decostand(fin_div_vp$N0, "standardize"), 
          fin_env_sel$Season, cbind(fin_geographic, fin_env_sel$pH, fin_nutrients)))
# Fraction [b], pure trend
anova(rda(decostand(fin_div_vp$N0, "standardize"), 
          fin_geographic, cbind(fin_env_sel$Season, fin_env_sel$pH, fin_nutrients)))
# Fraction [c], pure broad scale spatial
anova(rda(decostand(fin_div_vp$N0, "standardize"), 
          fin_env_sel$pH, cbind(fin_env_sel$Season, fin_geographic, fin_nutrients)))
# Fraction [d], pure fine scale spatial
anova(rda(decostand(fin_div_vp$N0, "standardize"), 
          fin_nutrients, cbind(fin_env_sel$Season, fin_geographic, fin_env_sel$pH)))

## E10
hist(fin_div_vp$E10)
fin_varpart_E10 <- varpart(decostand(fin_div_vp$E10, "standardize"), 
                           as.numeric(factor(fin_env_sel$Season)),
                           fin_geographic,
                           fin_env_sel$pH,
                           fin_nutrients)
plot(fin_varpart_E10, bg = 2:5)
plot(fin_varpart_E10, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
fin_E10_varpart_fract <- fin_varpart_E10$part$fract[1:4,]
fin_E10_varpart_fract$Metric <- "E10"
fin_E10_varpart_fract
set.seed(1)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(decostand(fin_div_vp$E10, "standardize"), 
          fin_env_sel$Season, cbind(fin_geographic, fin_env_sel$pH, fin_nutrients)))
# Fraction [b], pure trend
anova(rda(decostand(fin_div_vp$E10, "standardize"), 
          fin_geographic, cbind(fin_env_sel$Season, fin_env_sel$pH, fin_nutrients)))
# Fraction [c], pure broad scale spatial
anova(rda(decostand(fin_div_vp$E10, "standardize"), 
          fin_env_sel$pH, cbind(fin_env_sel$Season, fin_geographic, fin_nutrients)))
# Fraction [d], pure fine scale spatial
anova(rda(decostand(fin_div_vp$E10, "standardize"), 
          fin_nutrients, cbind(fin_env_sel$Season, fin_geographic, fin_env_sel$pH)))

## N1
hist(fin_div_vp$N1)
fin_varpart_N1 <- varpart(decostand(fin_div_vp$N1, "standardize"), 
                          as.numeric(factor(fin_env_sel$Season)),
                          fin_geographic,
                          fin_env_sel$pH,
                          fin_nutrients)
plot(fin_varpart_N1, bg = 2:5)
plot(fin_varpart_N1, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
fin_N1_varpart_fract <- fin_varpart_N1$part$fract[1:4,]
fin_N1_varpart_fract$Metric <- "N1"
fin_N1_varpart_fract
set.seed(1)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(decostand(fin_div_vp$N1, "standardize"), 
          fin_env_sel$Season, cbind(fin_geographic, fin_env_sel$pH, fin_nutrients)))
# Fraction [b], pure trend
anova(rda(decostand(fin_div_vp$N1, "standardize"), 
          fin_geographic, cbind(fin_env_sel$Season, fin_env_sel$pH, fin_nutrients)))
# Fraction [c], pure broad scale spatial
anova(rda(decostand(fin_div_vp$N1, "standardize"), 
          fin_env_sel$pH, cbind(fin_env_sel$Season, fin_geographic, fin_nutrients)))
# Fraction [d], pure fine scale spatial
anova(rda(decostand(fin_div_vp$N1, "standardize"), 
          fin_nutrients, cbind(fin_env_sel$Season, fin_geographic, fin_env_sel$pH)))

## FRic
hist(log(fin_div_vp$FRic))
fin_varpart_FRic <- varpart(decostand(fin_div_vp$FRic, "standardize"), 
                            as.numeric(factor(fin_env_sel$Season)),
                            fin_geographic,
                            fin_env_sel$pH,
                            fin_nutrients)
plot(fin_varpart_FRic, bg = 2:5)
plot(fin_varpart_FRic, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
fin_FRic_varpart_fract <- fin_varpart_FRic$part$fract[1:4,]
fin_FRic_varpart_fract$Metric <- "FRic"
fin_FRic_varpart_fract
set.seed(1)
## Test fraction [a, b, c, d]
anova(fin_rda_FRic_a <- rda(decostand(fin_div_vp$FRic, "standardize") ~ fin_env_sel$Season +
                              Condition(fin_geographic$Alt + fin_geographic$PC1 + fin_geographic$PC2) +
                              Condition(fin_env_sel$pH) +
                              Condition(fin_nutrients$DOC + fin_nutrients$Tot.N + fin_nutrients$Tot.P)))
anova(fin_rda_FRic_b <- rda(decostand(fin_div_vp$FRic, "standardize") ~ fin_geographic$Alt + fin_geographic$PC1 + fin_geographic$PC2 +
                              Condition(fin_env_sel$Season) +
                              Condition(fin_env_sel$pH) +
                              Condition(fin_nutrients$DOC + fin_nutrients$Tot.N + fin_nutrients$Tot.P)))
anova(fin_rda_FRic_c <- rda(decostand(fin_div_vp$FRic, "standardize") ~ fin_env_sel$pH +
                              Condition(fin_env_sel$Season) +
                              Condition(fin_geographic$Alt + fin_geographic$PC1 + fin_geographic$PC2) +
                              Condition(fin_nutrients$DOC + fin_nutrients$Tot.N + fin_nutrients$Tot.P)))
anova(fin_rda_FRic_d <- rda(decostand(fin_div_vp$FRic, "standardize") ~ fin_nutrients$DOC + fin_nutrients$Tot.N + fin_nutrients$Tot.P +
                              Condition(fin_env_sel$Season) +
                              Condition(fin_geographic$Alt + fin_geographic$PC1 + fin_geographic$PC2) +
                              Condition(fin_env_sel$pH)))

## FEve
hist(fin_div_vp$FEve)
fin_varpart_FEve <- varpart(decostand(fin_div_vp$FEve, "standardize"), 
                            as.numeric(factor(fin_env_sel$Season)),
                            fin_geographic,
                            fin_env_sel$pH,
                            fin_nutrients)
plot(fin_varpart_FEve, bg = 2:5)
plot(fin_varpart_FEve, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
fin_FEve_varpart_fract <- fin_varpart_FEve$part$fract[1:4,]
fin_FEve_varpart_fract$Metric <- "FEve"
fin_FEve_varpart_fract
set.seed(1)
## Test fraction [a, b, c, d]
anova(fin_rda_FEve_a <- rda(decostand(fin_div_vp$FEve, "standardize") ~ fin_env_sel$Season +
                              Condition(fin_geographic$Alt + fin_geographic$PC1 + fin_geographic$PC2) +
                              Condition(fin_env_sel$pH) +
                              Condition(fin_nutrients$DOC + fin_nutrients$Tot.N + fin_nutrients$Tot.P)))
anova(fin_rda_FEve_b <- rda(decostand(fin_div_vp$FEve, "standardize") ~ fin_geographic$Alt + fin_geographic$PC1 + fin_geographic$PC2 +
                              Condition(fin_env_sel$Season) +
                              Condition(fin_env_sel$pH) +
                              Condition(fin_nutrients$DOC + fin_nutrients$Tot.N + fin_nutrients$Tot.P)))
anova(fin_rda_FEve_c <- rda(decostand(fin_div_vp$FEve, "standardize") ~ fin_env_sel$pH +
                              Condition(fin_env_sel$Season) +
                              Condition(fin_geographic$Alt + fin_geographic$PC1 + fin_geographic$PC2) +
                              Condition(fin_nutrients$DOC + fin_nutrients$Tot.N + fin_nutrients$Tot.P)))
anova(fin_rda_FEve_d <- rda(decostand(fin_div_vp$FEve, "standardize") ~ fin_nutrients$DOC + fin_nutrients$Tot.N + fin_nutrients$Tot.P +
                              Condition(fin_env_sel$Season) +
                              Condition(fin_geographic$Alt + fin_geographic$PC1 + fin_geographic$PC2) +
                              Condition(fin_env_sel$pH)))

## FDis
hist(fin_div_vp$FDis)
fin_varpart_FDis <- varpart(decostand(fin_div_vp$FDis, "standardize"), 
                            as.numeric(factor(fin_env_sel$Season)),
                            fin_geographic,
                            fin_env_sel$pH,
                            fin_nutrients)
plot(fin_varpart_FDis)
plot(fin_varpart_FDis, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
fin_FDis_varpart_fract <- fin_varpart_FDis$part$fract[1:4,]
fin_FDis_varpart_fract$Metric <- "FDis"
fin_FDis_varpart_fract
set.seed(1)
## Test fraction [a, b, c, d]
anova(fin_rda_FDis_a <- rda(decostand(fin_div_vp$FDis, "standardize") ~ fin_env_sel$Season +
                              Condition(fin_geographic$Alt + fin_geographic$PC1 + fin_geographic$PC2) +
                              Condition(fin_env_sel$pH) +
                              Condition(fin_nutrients$DOC + fin_nutrients$Tot.N + fin_nutrients$Tot.P)))
anova(fin_rda_FDis_b <- rda(decostand(fin_div_vp$FDis, "standardize") ~ fin_geographic$Alt + fin_geographic$PC1 + fin_geographic$PC2 +
                              Condition(fin_env_sel$Season) +
                              Condition(fin_env_sel$pH) +
                              Condition(fin_nutrients$DOC + fin_nutrients$Tot.N + fin_nutrients$Tot.P)))
anova(fin_rda_FDis_c <- rda(decostand(fin_div_vp$FDis, "standardize") ~ fin_env_sel$pH +
                              Condition(fin_env_sel$Season) +
                              Condition(fin_geographic$Alt + fin_geographic$PC1 + fin_geographic$PC2) +
                              Condition(fin_nutrients$DOC + fin_nutrients$Tot.N + fin_nutrients$Tot.P)))
anova(fin_rda_FDis_d <- rda(decostand(fin_div_vp$FDis, "standardize") ~ fin_nutrients$DOC + fin_nutrients$Tot.N + fin_nutrients$Tot.P +
                              Condition(fin_env_sel$Season) +
                              Condition(fin_geographic$Alt + fin_geographic$PC1 + fin_geographic$PC2) +
                              Condition(fin_env_sel$pH)))

## FRic.SES
hist(fin_div_vp$FRic.SES)
fin_varpart_FRic.SES <- varpart(decostand(fin_div_vp$FRic.SES, "standardize"), 
                                as.numeric(factor(fin_env_sel$Season)),
                                fin_geographic,
                                fin_env_sel$pH,
                                fin_nutrients)
plot(fin_varpart_FRic.SES, bg = 2:5)
plot(fin_varpart_FRic.SES, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
fin_FRic.SES_varpart_fract <- fin_varpart_FRic.SES$part$fract[1:4,]
fin_FRic.SES_varpart_fract$Metric <- "FRic.SES"
fin_FRic.SES_varpart_fract
set.seed(1)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(decostand(fin_div_vp$FRic.SES, "standardize"), 
          fin_env_sel$Season, cbind(fin_geographic, fin_env_sel$pH, fin_nutrients)))
# Fraction [b], pure trend
anova(rda(decostand(fin_div_vp$FRic.SES, "standardize"), 
          fin_geographic, cbind(fin_env_sel$Season, fin_env_sel$pH, fin_nutrients)))
# Fraction [c], pure broad scale spatial
anova(rda(decostand(fin_div_vp$FRic.SES, "standardize"), 
          fin_env_sel$pH, cbind(fin_env_sel$Season, fin_geographic, fin_nutrients)))
# Fraction [d], pure fine scale spatial
anova(rda(decostand(fin_div_vp$FRic.SES, "standardize"), 
          fin_nutrients, cbind(fin_env_sel$Season, fin_geographic, fin_env_sel$pH)))

## FEve.SES
hist(fin_div_vp$FEve.SES)
fin_varpart_FEve.SES <- varpart(decostand(fin_div_vp$FEve.SES, "standardize"), 
                                as.numeric(factor(fin_env_sel$Season)),
                                fin_geographic,
                                fin_env_sel$pH,
                                fin_nutrients)
plot(fin_varpart_FEve.SES, bg = 2:5)
plot(fin_varpart_FEve.SES, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
fin_FEve.SES_varpart_fract <- fin_varpart_FEve.SES$part$fract[1:4,]
fin_FEve.SES_varpart_fract$Metric <- "FEve.SES"
fin_FEve.SES_varpart_fract
set.seed(1)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(decostand(fin_div_vp$FEve.SES, "standardize"), 
          fin_env_sel$Season, cbind(fin_geographic, fin_env_sel$pH, fin_nutrients)))
# Fraction [b], pure trend
anova(rda(decostand(fin_div_vp$FEve.SES, "standardize"), 
          fin_geographic, cbind(fin_env_sel$Season, fin_env_sel$pH, fin_nutrients)))
# Fraction [c], pure broad scale spatial
anova(rda(decostand(fin_div_vp$FEve.SES, "standardize"), 
          fin_env_sel$pH, cbind(fin_env_sel$Season, fin_geographic, fin_nutrients)))
# Fraction [d], pure fine scale spatial
anova(rda(decostand(fin_div_vp$FEve.SES, "standardize"), 
          fin_nutrients, cbind(fin_env_sel$Season, fin_geographic, fin_env_sel$pH)))

## FDis.SES
hist(fin_div_vp$FDis.SES)
fin_varpart_FDis.SES <- varpart(decostand(fin_div_vp$FDis.SES, "standardize"), 
                                as.numeric(factor(fin_env_sel$Season)),
                                fin_geographic,
                                fin_env_sel$pH,
                                fin_nutrients)
plot(fin_varpart_FDis.SES, bg = 2:5)
plot(fin_varpart_FDis.SES, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
fin_FDis.SES_varpart_fract <- fin_varpart_FDis.SES$part$fract[1:4,]
fin_FDis.SES_varpart_fract$Metric <- "FDis.SES"
fin_FDis.SES_varpart_fract
set.seed(1)
## Test fraction [a, b, c, d]
anova(fin_rda_FDis.SES_a <- rda(decostand(fin_div_vp$FDis.SES, "standardize") ~ fin_env_sel$Season +
                                  Condition(fin_geographic$Alt + fin_geographic$PC1 + fin_geographic$PC2) +
                                  Condition(fin_env_sel$pH) +
                                  Condition(fin_nutrients$DOC + fin_nutrients$Tot.N + fin_nutrients$Tot.P)))
anova(fin_rda_FDis.SES_b <- rda(decostand(fin_div_vp$FDis.SES, "standardize") ~ fin_geographic$Alt + fin_geographic$PC1 + fin_geographic$PC2 +
                                  Condition(fin_env_sel$Season) +
                                  Condition(fin_env_sel$pH) +
                                  Condition(fin_nutrients$DOC + fin_nutrients$Tot.N + fin_nutrients$Tot.P)))
anova(fin_rda_FDis.SES_c <- rda(decostand(fin_div_vp$FDis.SES, "standardize") ~ fin_env_sel$pH +
                                  Condition(fin_env_sel$Season) +
                                  Condition(fin_geographic$Alt + fin_geographic$PC1 + fin_geographic$PC2) +
                                  Condition(fin_nutrients$DOC + fin_nutrients$Tot.N + fin_nutrients$Tot.P)))
anova(fin_rda_FDis.SES_d <- rda(decostand(fin_div_vp$FDis.SES, "standardize") ~ fin_nutrients$DOC + fin_nutrients$Tot.N + fin_nutrients$Tot.P +
                                  Condition(fin_env_sel$Season) +
                                  Condition(fin_geographic$Alt + fin_geographic$PC1 + fin_geographic$PC2) +
                                  Condition(fin_env_sel$pH)))
## FRed
hist(fin_div_vp$FRed)
fin_varpart_FRed <- varpart(decostand(fin_div_vp$FRed, "standardize"), 
                            as.numeric(factor(fin_env_sel$Season)),
                            fin_geographic,
                            fin_env_sel$pH,
                            fin_nutrients)
plot(fin_varpart_FRed, bg = 2:5)
plot(fin_varpart_FRed, cutoff = -Inf, cex = 0.7, bg = 2:5) # Show values for all partitions by putting 'cutoff' low enough:
fin_FRed_varpart_fract <- fin_varpart_FRed$part$fract[1:4,]
fin_FRed_varpart_fract$Metric <- "FRed"
fin_FRed_varpart_fract
set.seed(1)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(decostand(fin_div_vp$FRed, "standardize"), 
          fin_env_sel$Season, cbind(fin_geographic, fin_env_sel$pH, fin_nutrients)))
# Fraction [b], pure trend
anova(rda(decostand(fin_div_vp$FRed, "standardize"), 
          fin_geographic, cbind(fin_env_sel$Season, fin_env_sel$pH, fin_nutrients)))
# Fraction [c], pure broad scale spatial
anova(rda(decostand(fin_div_vp$FRed, "standardize"), 
          fin_env_sel$pH, cbind(fin_env_sel$Season, fin_geographic, fin_nutrients)))
# Fraction [d], pure fine scale spatial
anova(rda(decostand(fin_div_vp$FRed, "standardize"), 
          fin_nutrients, cbind(fin_env_sel$Season, fin_geographic, fin_env_sel$pH)))

fin_varpart_fract <- rbind(fin_A_varpart_fract, fin_N0_varpart_fract, fin_E10_varpart_fract, fin_N1_varpart_fract, 
                           fin_FRic_varpart_fract, fin_FEve_varpart_fract, fin_FDis_varpart_fract,
                           fin_FRic.SES_varpart_fract, fin_FEve.SES_varpart_fract, fin_FDis.SES_varpart_fract,
                           fin_FRed_varpart_fract)
(fin_varpart_fract <- fin_varpart_fract[, c(3, 5)])
rownames(fin_varpart_fract) <- c("Season", "Geography", "Acidity", "Nutrients", 
                                 "Season0", "Geography0", "Acidity0", "Nutrients0", 
                                 "Season1", "Geography1", "Acidity1", "Nutrients1",
                                 "Season2", "Geography2", "Acidity2", "Nutrients2", 
                                 "Season3", "Geography3", "Acidity3", "Nutrients3", 
                                 "Season4", "Geography4", "Acidity4", "Nutrients4", 
                                 "Season5", "Geography5", "Acidity5", "Nutrients5",
                                 "Season6", "Geography6", "Acidity6", "Nutrients6", 
                                 "Season7", "Geography7", "Acidity7", "Nutrients7", 
                                 "Season8", "Geography8", "Acidity8", "Nutrients8", 
                                 "Season9", "Geography9", "Acidity9", "Nutrients9")

#### REGRESSIONS OF VARPART OUTPUTS (TO GET DIRECTIONS OF CHANGE) --------------------
### GERMANY --------------------
ger_env_sel$Season
ger_geographic
ger_env_sel$pH
ger_nutrients

### significant variables as identified by varpart
## log(Abund)
ger_A_season_lm <- lm(log(ger_div_vp$A) ~ ger_env_sel$Season)
summary(ger_A_season_lm)
ger_A_geo_lm <- lm(log(ger_div_vp$A) ~ ger_geographic$Alt)
summary(ger_A_geo_lm)
ger_A_pH_lm <- lm(log(ger_div_vp$A) ~ ger_env_sel$pH)
summary(ger_A_pH_lm)
ger_A_nut_lm <- lm(log(ger_div_vp$A) ~ ger_nutrients$DOC + ger_nutrients$SO4 + ger_nutrients$NH4.N)
summary(ger_A_nut_lm)

## TRic
ger_N0_season_lm <- lm(ger_div_vp$N0 ~ ger_env_sel$Season)
summary(ger_N0_season_lm)
ger_N0_geo_lm <- lm(ger_div_vp$N0 ~ ger_geographic$Alt)
summary(ger_N0_geo_lm)
ger_N0_pH_lm <- lm(ger_div_vp$N0 ~ ger_env_sel$pH)
summary(ger_N0_pH_lm)
ger_N0_nut_lm <- lm(ger_div_vp$N0 ~ ger_nutrients$DOC + ger_nutrients$SO4 + ger_nutrients$NH4.N)
summary(ger_N0_nut_lm)

## TEve
ger_E10_season_lm <- lm(ger_div_vp$E10 ~ ger_env_sel$Season)
summary(ger_E10_season_lm)
ger_E10_geo_lm <- lm(ger_div_vp$E10 ~ ger_geographic$Alt)
summary(ger_E10_geo_lm)
ger_E10_pH_lm <- lm(ger_div_vp$E10 ~ ger_env_sel$pH)
summary(ger_E10_pH_lm)
ger_E10_nut_lm <- lm(ger_div_vp$E10 ~ ger_nutrients$DOC + ger_nutrients$SO4 + ger_nutrients$NH4.N)
summary(ger_E10_nut_lm)

## Shan
ger_N1_season_lm <- lm(ger_div_vp$N1 ~ ger_env_sel$Season)
summary(ger_N1_season_lm)
ger_N1_geo_lm <- lm(ger_div_vp$N1 ~ ger_geographic$Alt)
summary(ger_N1_geo_lm)
ger_N1_pH_lm <- lm(ger_div_vp$N1 ~ ger_env_sel$pH)
summary(ger_N1_pH_lm)
ger_N1_nut_lm <- lm(ger_div_vp$N1 ~ ger_nutrients$DOC + ger_nutrients$SO4 + ger_nutrients$NH4.N)
summary(ger_N1_nut_lm)

## FRic.SES
ger_FRic.SES_season_lm <- lm(ger_div_vp$FRic.SES ~ ger_env_sel$Season)
summary(ger_FRic.SES_season_lm)
ger_FRic.SES_geo_lm <- lm(ger_div_vp$FRic.SES ~ ger_geographic$Alt)
summary(ger_FRic.SES_geo_lm)
ger_FRic.SES_pH_lm <- lm(ger_div_vp$FRic.SES ~ ger_env_sel$pH)
summary(ger_FRic.SES_pH_lm)
ger_FRic.SES_nut_lm <- lm(ger_div_vp$FRic.SES ~ ger_nutrients$DOC + ger_nutrients$SO4 + ger_nutrients$NH4.N)
summary(ger_FRic.SES_nut_lm)

## FEve.SES
ger_FEve.SES_season_lm <- lm(ger_div_vp$FEve.SES ~ ger_env_sel$Season)
summary(ger_FEve.SES_season_lm)
ger_FEve.SES_geo_lm <- lm(ger_div_vp$FEve.SES ~ ger_geographic$Alt)
summary(ger_FEve.SES_geo_lm)
ger_FEve.SES_pH_lm <- lm(ger_div_vp$FEve.SES ~ ger_env_sel$pH)
summary(ger_FEve.SES_pH_lm)
ger_FEve.SES_nut_lm <- lm(ger_div_vp$FEve.SES ~ ger_nutrients$DOC + ger_nutrients$SO4 + ger_nutrients$NH4.N)
summary(ger_FEve.SES_nut_lm)

## FDis.SES
ger_FDis.SES_season_lm <- lm(ger_div_vp$FDis.SES ~ ger_env_sel$Season)
summary(ger_FDis.SES_season_lm)
ger_FDis.SES_geo_lm <- lm(ger_div_vp$FDis.SES ~ ger_geographic$Alt)
summary(ger_FDis.SES_geo_lm)
ger_FDis.SES_pH_lm <- lm(ger_div_vp$FDis.SES ~ ger_env_sel$pH)
summary(ger_FDis.SES_pH_lm)
ger_FDis.SES_nut_lm <- lm(ger_div_vp$FDis.SES ~ ger_nutrients$DOC + ger_nutrients$SO4 + ger_nutrients$NH4.N)
summary(ger_FDis.SES_nut_lm)

## FRed
ger_FRed_season_lm <- lm(ger_div_vp$FRed ~ ger_env_sel$Season)
summary(ger_FRed_season_lm)
ger_FRed_geo_lm <- lm(ger_div_vp$FRed ~ ger_geographic$Alt)
summary(ger_FRed_geo_lm)
ger_FRed_pH_lm <- lm(ger_div_vp$FRed ~ ger_env_sel$pH)
summary(ger_FRed_pH_lm)
ger_FRed_nut_lm <- lm(ger_div_vp$FRed ~ ger_nutrients$DOC + ger_nutrients$SO4 + ger_nutrients$NH4.N)
summary(ger_FRed_nut_lm)

### SWEDEN --------------------
swe_env_sel$Season
swe_geographic
swe_env_sel$pH
swe_nutrients

### significant variables as identified by varpart
## log(Abund)
swe_A_season_lm <- lm(log(swe_div_vp$A) ~ swe_env_sel$Season)
summary(swe_A_season_lm)
swe_A_geo_lm <- lm(log(swe_div_vp$A) ~ swe_geographic$Alt)
summary(swe_A_geo_lm)
swe_A_pH_lm <- lm(log(swe_div_vp$A) ~ swe_env_sel$pH)
summary(swe_A_pH_lm)
swe_A_nut_lm <- lm(log(swe_div_vp$A) ~ swe_nutrients$TOC + swe_nutrients$SO4 + swe_nutrients$NO2.NO3_N + swe_nutrients$Tot.P)
summary(swe_A_nut_lm)

## TRic
swe_N0_season_lm <- lm(swe_div_vp$N0 ~ swe_env_sel$Season)
summary(swe_N0_season_lm)
swe_N0_geo_lm <- lm(swe_div_vp$N0 ~ swe_geographic$Alt)
summary(swe_N0_geo_lm)
swe_N0_pH_lm <- lm(swe_div_vp$N0 ~ swe_env_sel$pH)
summary(swe_N0_pH_lm)
swe_N0_nut_lm <- lm(swe_div_vp$N0 ~ swe_nutrients$TOC + swe_nutrients$SO4 + swe_nutrients$NO2.NO3_N + swe_nutrients$Tot.P)
summary(swe_N0_nut_lm)

## TEve
swe_E10_season_lm <- lm(swe_div_vp$E10 ~ swe_env_sel$Season)
summary(swe_E10_season_lm)
swe_E10_geo_lm <- lm(swe_div_vp$E10 ~ swe_geographic$Alt)
summary(swe_E10_geo_lm)
swe_E10_pH_lm <- lm(swe_div_vp$E10 ~ swe_env_sel$pH)
summary(swe_E10_pH_lm)
swe_E10_nut_lm <- lm(swe_div_vp$E10 ~ swe_nutrients$TOC + swe_nutrients$SO4 + swe_nutrients$NO2.NO3_N + swe_nutrients$Tot.P)
summary(swe_E10_nut_lm)

## Shan
swe_N1_season_lm <- lm(swe_div_vp$N1 ~ swe_env_sel$Season)
summary(swe_N1_season_lm)
swe_N1_geo_lm <- lm(swe_div_vp$N1 ~ swe_geographic$Alt)
summary(swe_N1_geo_lm)
swe_N1_pH_lm <- lm(swe_div_vp$N1 ~ swe_env_sel$pH)
summary(swe_N1_pH_lm)
swe_N1_nut_lm <- lm(swe_div_vp$N1 ~ swe_nutrients$TOC + swe_nutrients$SO4 + swe_nutrients$NO2.NO3_N + swe_nutrients$Tot.P)
summary(swe_N1_nut_lm)

## FRic.SES
swe_FRic.SES_season_lm <- lm(swe_div_vp$FRic.SES ~ swe_env_sel$Season)
summary(swe_FRic.SES_season_lm)
swe_FRic.SES_geo_lm <- lm(swe_div_vp$FRic.SES ~ swe_geographic$Alt)
summary(swe_FRic.SES_geo_lm)
swe_FRic.SES_pH_lm <- lm(swe_div_vp$FRic.SES ~ swe_env_sel$pH)
summary(swe_FRic.SES_pH_lm)
swe_FRic.SES_nut_lm <- lm(swe_div_vp$FRic.SES ~ swe_nutrients$TOC + swe_nutrients$SO4 + swe_nutrients$NO2.NO3_N + swe_nutrients$Tot.P)
summary(swe_FRic.SES_nut_lm)

## FEve.SES
swe_FEve.SES_season_lm <- lm(swe_div_vp$FEve.SES ~ swe_env_sel$Season)
summary(swe_FEve.SES_season_lm)
swe_FEve.SES_geo_lm <- lm(swe_div_vp$FEve.SES ~ swe_geographic$Alt)
summary(swe_FEve.SES_geo_lm)
swe_FEve.SES_pH_lm <- lm(swe_div_vp$FEve.SES ~ swe_env_sel$pH)
summary(swe_FEve.SES_pH_lm)
swe_FEve.SES_nut_lm <- lm(swe_div_vp$FEve.SES ~ swe_nutrients$TOC + swe_nutrients$SO4 + swe_nutrients$NO2.NO3_N + swe_nutrients$Tot.P)
summary(swe_FEve.SES_nut_lm)

## FDis.SES
swe_FDis.SES_season_lm <- lm(swe_div_vp$FDis.SES ~ swe_env_sel$Season)
summary(swe_FDis.SES_season_lm)
swe_FDis.SES_geo_lm <- lm(swe_div_vp$FDis.SES ~ swe_geographic$Alt)
summary(swe_FDis.SES_geo_lm)
swe_FDis.SES_pH_lm <- lm(swe_div_vp$FDis.SES ~ swe_env_sel$pH)
summary(swe_FDis.SES_pH_lm)
swe_FDis.SES_nut_lm <- lm(swe_div_vp$FDis.SES ~ swe_nutrients$TOC + swe_nutrients$SO4 + swe_nutrients$NO2.NO3_N + swe_nutrients$Tot.P)
summary(swe_FDis.SES_nut_lm)

## FRed
swe_FRed_season_lm <- lm(swe_div_vp$FRed ~ swe_env_sel$Season)
summary(swe_FRed_season_lm)
swe_FRed_geo_lm <- lm(swe_div_vp$FRed ~ swe_geographic$Alt)
summary(swe_FRed_geo_lm)
swe_FRed_pH_lm <- lm(swe_div_vp$FRed ~ swe_env_sel$pH)
summary(swe_FRed_pH_lm)
swe_FRed_nut_lm <- lm(swe_div_vp$FRed ~ swe_nutrients$TOC + swe_nutrients$SO4 + swe_nutrients$NO2.NO3_N + swe_nutrients$Tot.P)
summary(swe_FRed_nut_lm)

### FINLAND --------------------
fin_env_sel$Season
fin_geographic
fin_env_sel$pH
fin_nutrients

### significant variables as identified by varpart
## log(Abund)
fin_A_season_lm <- lm(log(fin_div_vp$A) ~ fin_env_sel$Season)
summary(fin_A_season_lm)
fin_A_geo_lm <- lm(log(fin_div_vp$A) ~ fin_geographic$Alt)
summary(fin_A_geo_lm)
fin_A_pH_lm <- lm(log(fin_div_vp$A) ~ fin_env_sel$pH)
summary(fin_A_pH_lm)
fin_A_nut_lm <- lm(log(fin_div_vp$A) ~ fin_nutrients$DOC + fin_nutrients$Tot.P + fin_nutrients$Tot.N)
summary(fin_A_nut_lm)

## TRic
fin_N0_season_lm <- lm(fin_div_vp$N0 ~ fin_env_sel$Season)
summary(fin_N0_season_lm)
fin_N0_geo_lm <- lm(fin_div_vp$N0 ~ fin_geographic$Alt)
summary(fin_N0_geo_lm)
fin_N0_pH_lm <- lm(fin_div_vp$N0 ~ fin_env_sel$pH)
summary(fin_N0_pH_lm)
fin_N0_nut_lm <- lm(fin_div_vp$N0 ~ fin_nutrients$DOC + fin_nutrients$Tot.P + fin_nutrients$Tot.N)
summary(fin_N0_nut_lm)

## TEve
fin_E10_season_lm <- lm(fin_div_vp$E10 ~ fin_env_sel$Season)
summary(fin_E10_season_lm)
fin_E10_geo_lm <- lm(fin_div_vp$E10 ~ fin_geographic$Alt)
summary(fin_E10_geo_lm)
fin_E10_pH_lm <- lm(fin_div_vp$E10 ~ fin_env_sel$pH)
summary(fin_E10_pH_lm)
fin_E10_nut_lm <- lm(fin_div_vp$E10 ~ fin_nutrients$DOC + fin_nutrients$Tot.P + fin_nutrients$Tot.N)
summary(fin_E10_nut_lm)

## Shan
fin_N1_season_lm <- lm(fin_div_vp$N1 ~ fin_env_sel$Season)
summary(fin_N1_season_lm)
fin_N1_geo_lm <- lm(fin_div_vp$N1 ~ fin_geographic$Alt)
summary(fin_N1_geo_lm)
fin_N1_pH_lm <- lm(fin_div_vp$N1 ~ fin_env_sel$pH)
summary(fin_N1_pH_lm)
fin_N1_nut_lm <- lm(fin_div_vp$N1 ~ fin_nutrients$DOC + fin_nutrients$Tot.P + fin_nutrients$Tot.N)
summary(fin_N1_nut_lm)

## FRic.SES
fin_FRic.SES_season_lm <- lm(fin_div_vp$FRic.SES ~ fin_env_sel$Season)
summary(fin_FRic.SES_season_lm)
fin_FRic.SES_geo_lm <- lm(fin_div_vp$FRic.SES ~ fin_geographic$Alt)
summary(fin_FRic.SES_geo_lm)
fin_FRic.SES_pH_lm <- lm(fin_div_vp$FRic.SES ~ fin_env_sel$pH)
summary(fin_FRic.SES_pH_lm)
fin_FRic.SES_nut_lm <- lm(fin_div_vp$FRic.SES ~ fin_nutrients$DOC + fin_nutrients$Tot.P + fin_nutrients$Tot.N)
summary(fin_FRic.SES_nut_lm)

## FEve.SES
fin_FEve.SES_season_lm <- lm(fin_div_vp$FEve.SES ~ fin_env_sel$Season)
summary(fin_FEve.SES_season_lm)
fin_FEve.SES_geo_lm <- lm(fin_div_vp$FEve.SES ~ fin_geographic$Alt)
summary(fin_FEve.SES_geo_lm)
fin_FEve.SES_pH_lm <- lm(fin_div_vp$FEve.SES ~ fin_env_sel$pH)
summary(fin_FEve.SES_pH_lm)
fin_FEve.SES_nut_lm <- lm(fin_div_vp$FEve.SES ~ fin_nutrients$DOC + fin_nutrients$Tot.P + fin_nutrients$Tot.N)
summary(fin_FEve.SES_nut_lm)

## FDis.SES
fin_FDis.SES_season_lm <- lm(fin_div_vp$FDis.SES ~ fin_env_sel$Season)
summary(fin_FDis.SES_season_lm)
fin_FDis.SES_geo_lm <- lm(fin_div_vp$FDis.SES ~ fin_geographic$Alt)
summary(fin_FDis.SES_geo_lm)
fin_FDis.SES_pH_lm <- lm(fin_div_vp$FDis.SES ~ fin_env_sel$pH)
summary(fin_FDis.SES_pH_lm)
fin_FDis.SES_nut_lm <- lm(fin_div_vp$FDis.SES ~ fin_nutrients$DOC + fin_nutrients$Tot.P + fin_nutrients$Tot.N)
summary(fin_FDis.SES_nut_lm)

## FRed
fin_FRed_season_lm <- lm(fin_div_vp$FRed ~ fin_env_sel$Season)
summary(fin_FRed_season_lm)
fin_FRed_geo_lm <- lm(fin_div_vp$FRed ~ fin_geographic$Alt)
summary(fin_FRed_geo_lm)
fin_FRed_pH_lm <- lm(fin_div_vp$FRed ~ fin_env_sel$pH)
summary(fin_FRed_pH_lm)
fin_FRed_nut_lm <- lm(fin_div_vp$FRed ~ fin_nutrients$DOC + fin_nutrients$Tot.P + fin_nutrients$Tot.N)
summary(fin_FRed_nut_lm)

##### CLEAN UP --------------------

# Clear data
rm(list = ls())  # Removes all objects from environment

# Clear packages
p_unload(all)  # Remove all contributed packages

# Clear plots
graphics.off()  # Clears plots, closes all graphics devices

# Clear console
cat("\014")  # Mimics ctrl+L

# Clear mind :)