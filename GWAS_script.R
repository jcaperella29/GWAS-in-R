library(statgenGWAS)
data(dropsMarkers)
data(dropsMap)
data(dropsPheno)
library(dplyr)

## Add genotypes as row names of dropsMarkers and drop Ind column.
rownames(dropsMarkers) <- dropsMarkers[["Ind"]]
dropsMarkers <- dropsMarkers[colnames(dropsMarkers) != "Ind"]

## Add genotypes as row names of dropsMap.
rownames(dropsMap) <- dropsMap[["SNP.names"]]
## Rename Chomosome and Position columns.
colnames(dropsMap)[match(c("Chromosome", "Position"), colnames(dropsMap))] <- c("chr", "pos")

## Create a gData object containing map and marker information.
gDataDrops <- createGData(geno = dropsMarkers, map = dropsMap)

## Rename Variety_ID to genotype.
colnames(dropsPheno)[colnames(dropsPheno) == "Variety_ID"] <- "genotype"
## Select relevant columns and convert data to a list.
dropsPhenoList <- split(x = dropsPheno[c("genotype", "grain.yield",
                                         "grain.number", "seed.size",
                                         "anthesis", "silking", "plant.height",
                                         "tassel.height", "ear.height")], 
                        f = dropsPheno[["Experiment"]])
## Add phenotypic data to gDataDrops.
gDataDrops <- createGData(gData = gDataDrops, pheno = dropsPhenoList) 

## Summarize gDataDrops.
summary(gDataDrops, trials = "Mur13W")


## Plot genetic map.
plot(gDataDrops)
## Remove duplicate SNPs from gDataDrops.
gDataDropsDedup <- codeMarkers(gDataDrops, impute = FALSE, verbose = TRUE) 


## Run single trait GWAS for traits 'grain.yield' and 'anthesis' for trial Mur13W.
GWASDrops <- runSingleTraitGwas(gData = gDataDropsDedup,
                                trials = "Mur13W",
                                traits = c("grain.yield", "anthesis"))

summary(GWASDrops)

GWAS_DF<-data.frame(GWASDrops$GWAResult$Mur13W)
GWAS_DF_SORTED<- GWAS_DF[order(GWAS_DF$pValue),]

#outputting a CVS file with the SNP sorted by P-value
write.csv(GWAS_DF_SORTED,"path choosen by user")

#checking  for SNPs of Genome Wide Significance in the  data.

GWAS_DF_ONLY_GENOME_WIDE_SIG<- GWAS_DF_SORTED %>% filter(pValue<= 5e-8)
#outputting a text file containing a list of SNPs that have pvalues less than or equal  or to 5e-8

write.table(GWAS_DF_ONLY_GENOME_WIDE_SIG$snp,"path choosen by the user")


## Plot a manhattan plot of GWAS Drops with significance threshold 8.
## Assume PZE-106021410 and PZE-105012420 are SNPs with known effects.
plot(GWASDrops, plotType = "manhattan", trait = "grain.yield", 
     effects = c("PZE-106021410", "PZE-105012420"),lod = 8)





#orignal code is from https://cran.r-project.org/web/packages/statgenGWAS/vignettes/GWAS.html , I made changes to add in reporting by filtering for SNP that meet the cutoff for Genome Wide Sigficance and outputting results as a csv and text file.
