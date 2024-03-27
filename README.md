# GWAS-in-R
An R script for processing GWAS data.This produces a Manhattan plot , a cvs file of the results and a text file containing the SNPs met the p-value cutoff for genome wide significance
After reading in sample and data and needed libraries  we create  add a column to the dropsmap object  cotaining the SNP name and  rename chrosome and positon names as chr and pos respectively.
We then make a gdata object using the map and marker information.
Next we rename the Variety_ID column as genotype.
Then we isolate the phenotype information and make it into a list.
Next we add that information to the to the gdata object.
We then read the summary of the gdata object, before ploting a genetic map and removing duplicates.
We then perform the GWAS looking at our traits of interest.
Results of the GWAS are captured in a dataframe and sorted in asscending order by p-value
that dataframe is then output as csv file
Next, we isolate SNPs that met the p-value cutoff for genenome wide significanace (<5e-80)
Those SNPs are output as a text file.
Finally we produce a mahattan plot.
