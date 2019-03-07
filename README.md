Authors: Iurii Nagornov, Mamoru Kato, Department of Bioinformatics, Research Institute, National Cancer Center, Tokyo, Japan
"Cancer cell evolution"

This is a script in R to simulate the cancer cell evolution in the framework of Hallmarks model proposed by Mamoru Kato (Head of Bioinformatics Department, Research Institute, Nation Cancer Center, Tokyo, JAPAN).

The directory has several files:

	•	Tests.zip - archive with tests for program code. Each test is in the separate directory and has a note with explanation of details of test. Some tests need to change a code a bit, you can find information about change in the note and also in the code, for this, please, find word "test".
    
	•	presentation_of_tests.pdf - presentation file with results of all tests. You can find all tests and results in a brief manner.

	•	CancerProgress_Final_Version.R - the R script for simulation (please, see input files in the code and see a format of the input files in the presentation). Other script files are only for analysis the results of simulations. The script allows to make repeating simulation/replicas to get statistical data. 
	
	•	Report.pdf - the file with a plots and graphs of results of simulations for three cases. 
	
	•	Analyze_of_clones.R - the R script to analize the calculation data. The analyzing script allows to calculate next dependences: an evolution of the cell’s number for the normal and mutated cells, an evolution of the Hallmarks and the probabilities, a number of the cells in clones vs clone’s ID (for drivers and for all genes), a number of clones vs the replica’s ID (for drivers and for all genes and the relations between them) and it’s distribution, a histogram for the inequalities coefficients for drivers and for all genes, the distribution of the variant allele function and the order of the gene’s dysfunction.
	
	•	Functions.R - script to calculate the functions for Analyze_of_clones.R, so please run it first. 
		
	•	The examples of the simulation results for the dummy cases and for a more realistic case as colorectal cancer:
        
	    ◦	single.zip - example of simulation for simplest dependence of Hallmark on genes: one gene - one Hallmarks.
	    ◦	multi.zip - example of simulation for complex dependence of Hallmark on genes.
	    ◦	cancer.zip - example of simulation for colorectal cancer.

