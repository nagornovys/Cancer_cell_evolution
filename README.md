Authors: Iurii Nagornov, Mamoru Kato, Department of Bioinformatics, Research Institute, National Cancer Center, Tokyo, Japan
"Cancer cell evolution"

This is a script in R to simulate the cancer cell evolution in the framework of Hallmarks model proposed by Mamoru Kato (Head of Bioinformatics Department, Research Institute, Nation Cancer Center, Tokyo, JAPAN).

The directory has several files:

	•	Tests.zip - archive with tests for program code. Each test is in the separate directory and has a note with explanation of details of test. Some tests need to change a code a bit, you can find information about change in the note and also in the code, for this, please, find word "test".
    
	•	presentation_of_tests.pdf - presentation file with results of all tests. You can find all tests and results in a brief manner.
    
	•	Final_version with clones.zip - the directory this final version of program code, which has a several files:
    
	    ◦	CancerProgress_Final_Version.R - the R script for simulation (please, see input files in the code and see a format of the input files in the presentation). Other script files are only for analysis the results of simulations.
        
	    ◦	ineq.R - the R script for ineq function with Gini parameter, which is used in the Diversity of clones.R file. First script to run.
        
	    ◦	Analize.R - the R script with functions for Diversity of clones.R script, which allows to calculate time evolutions of number of metastasis and normal cells, Hallmarks variables, probabilities variables and to calculate the order of gene dysfunction. Second script to run.
        
	    ◦	Diversity of clones.R - the R script with additional analytical functions: diversity if clones, number of cells in each clones, inequality coefficient and variant allele frequency. Third script to run.
        
        
        
