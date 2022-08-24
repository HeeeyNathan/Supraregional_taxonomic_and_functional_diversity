# Analysing supraregional taxonomic and functional diversity

Code to re-create analysis for "Seasonal and spatial variation of stream macroinvertebrate taxonomic and functional diversity across three boreal regions".

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img src="https://i.creativecommons.org/l/by/4.0/88x31.png" alt="Creative Commons License" style="border-width:0"/></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

## Original article:

Please, use this citation to reference the code:

    Baker et al. 2022. Seasonal and spatial variation of stream macroinvertebrate 
    taxonomic and functional diversity across three boreal regions. 
    Insect Conservation and Diversity xx, xx-xx. DOI:

## Description of R files:

-   1_R_script_supraregional_analysis.R: Main code to reproduce the results presented in the paper

#### Source:

Zuur, A. F., Ieno, E. N. & Smith, G. M. Analyzing Ecological Data. Methods (2007). <doi:10.1016/B978-0-12-387667-6.00013-0>.

-   HighstatLibV10.R - For pairplots of response and explanatory variables

#### Source:

Borcard, D., Gillet, F. & Legendre, P. Numerical Ecology with R. (Springer New York, 2018). <doi:10.1007/978-1-4419-7976-6>.

Additional functions needed for ordinations 
-   hcoplot.R 
-   triplot.rda.R
-   plot.lda.R
-   polyvars.R
-   screestick.R
-   panelutils.R
-   Rao.R

#### Source:

<https://github.com/tanogc/overarching_functional_space>

-   0_FD_functions.R - R script to estimate Functional Diversity (FD) metrics This function is a modification of the functions within the 'FD' R package (Laliberté et al., 2014)
-   0_quality_funct_space_fromdist.R - R function for computing the quality of functional dendrogramm and multidimensional functional spaces. This function is a simplified version of the Appendix S1 associated to Maire et al. 2015 (Global Ecol. and Biogeogr.)

#### Source:

<https://github.com/cran/vegetarian>

-   0_vegetarianRpackage_functions.R - Manual addition of the functions in the vegetarian R package - Jost Diversity Measures for Community Data. This package is no longer hosted on CRAN

## Original data:

-   full_site_taxa_mat.csv: Abundances of taxa across all sites. These data are the full dataset used for taxonomic diversity analyses (221 river macroinvertebrate taxa).
-   func_taxa_trait_mat.csv: Fuzzy coded trait data for 93 river macroinvertebrate taxa.
-   func_site_taxa_mat.csv: Abundances of taxa across all sites. These data are the taxonomically dataset used for functional diversity analyses (93 river macroinvertebrate taxa). 
-   func_reg_taxa_mat.csv: Regional and seasonal abundances of riverine macroinvertebrate taxa used for functional analyses.
-   site_env_mat.csv: Site-based environmental variable measurements. 
-   site_geo_mat.csv: Geographic coordinates of each site (in degree decimal) 
-   site_grouping_mat.csv: Site grouping factors

## Dependencies:

To run the code and functions from this repository, you need to install the following packages: 'ade4', 'adegraphics', 'adespatial', 'adiv', 'BAT', 'corrplot', 'EnvStats', 'FD', 'funrar', 'geodist', 'ggplot2', 'ks', 'labdsv', 'pacman', 'plyr', 'RColorBrewer', 'vegan'. Use this code to install them:

```{r}
install.packages(c("ade4", "adegraphics", "adespatial", "adiv", "BAT", "corrplot", "EnvStats", "FD", "funrar", "geodist", "ggplot2", "ks", "labdsv", "pacman", "plyr", "RColorBrewer", "vegan"))
```

Code written by Nathan J. Baker and Francesca Pilotto                                   Code used for functional analyses adapted from Cayetano Gutiérrez-Cánovas and Gabone Iturralde (<https://github.com/tanogc/overarching_functional_space>)

Email for queries: <Nathan93Baker@gmail.com>
Code written using R version 4.2.1 "Funny-Looking Kid"
Code written in Rstudio version 2022.07.1+554 "Spotted Wakerobin"