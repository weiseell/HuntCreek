# Hunt Creek
All scripts and data analysis associated with the Hunt Creek translocation event and subsequent genetic analysis (Weise et al 2020). 

## Guide to navigating the repository:
The repository structure mirrors the Methods section in Weise et al 2020. There are 4 folders that contain a different type of analysis. Each folder contains an input and output folder, as well as a source folder that contains scripts for required source functions. Finally, there is a script that will run the analysis. If multiple scripts are required in each folder, they are numbered and should be executed in numerical order.
###Analysis overviews:

Measures of genetic diversity:
	-Main objective: Calculate several population genetic summary statistics for each population:
		-average number of alleles
		-Allelic richness
		-Expected and observed heterozygosity
		-Fis
		-Fst
	-Outputs: 
		-a table containing several summary statistics organized by source population
		-a table containing significance tests and adjusted p-values for some of the summary statistics
		-a figure visualizing the results of the Fst calculations

Randomization:
	-Main objective: test the significance of a chi-squared analysis using randomization procedure to generate 1000 randomly mating population. The chi-squared values for the random populations are used as a null distribution to test for significance of the observed chi-squared values.
	-Outputs:

###Function overviews:

