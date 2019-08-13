tugHall
====================
tugHall: a simulator of cancer cell evolution based on the hallmarks of cancer, linked to the mutational states of tumor-related genes

Short description
---
This is a script in R to simulate the cancer cell evolution in the framework of Hallmarks model proposed by Mamoru Kato (Head of Bioinformatics Department, Research Institute, Nation Cancer Center, Tokyo, JAPAN).

Licence Information
------
GNU GENERAL PUBLIC LICENSE Version 3 from 29 June 2007 (READ LICENSE)

Project source can be downloaded from websites  
--- 
https://github.com/nagornovys/Cancer_cell_evolution  -  developed resource
https://www.doi.org/10.5281/zenodo.2667073  -  fixed resource in the form of an archive

Authors and contributor list: 
--- 
*Iurii Nagornov* 

*Mamoru Kato*

Department of Bioinformatics, Research Institute, National Cancer Center Japan, Tokyo, Japan

All other known bugs and fixes can be sent to inagonov@ncc.go.jp

Reported bugs/fixes will be submitted to correction.

File List
---
``` 
	•	Tests.zip - archive with tests for program code. Each test is in the separate directory and has a note with explanation of details of test. Some tests need to change a code a bit, you can find information about change in the note and also in the code, for this, please, find word "test".
    
	•	presentation_of_tests.pdf - presentation file with results of all tests. You can find all tests and results in a brief manner.
	•	Colorectal_cancer.zip - the zip-arhive with R script and input files for simulation (please, see input files in the code and see a format of the input files in the presentation). Other script files are only for analysis the results of simulations. The set of genes is included for colorectal cancer as a example. Please, to change the genes - Hallmarks relations to another simulation. 
```

Functions.R and Analyze_of_clones.R files
---
The zip arhives of simulations include the R scrips to analize the results of the simulations:
```
	•	Functions.R (ineq.R + Analyze.R) - script/scripts to calculate the functions for Analyze_of_clones.R (or Dyversity of clones.R), so please run it first. 
	•	Analyze_of_clones.R (or Dyversity of clones.R) - the R script to analize the calculation data, please, run it after Functions.R (ineq.R + Analyze.R). The analyzing script allows to calculate next dependences: an evolution of the cell’s number for the normal and mutated cells, an evolution of the Hallmarks and the probabilities, a number of the cells in clones vs clone’s ID (for drivers and for all genes), a number of clones vs the replica’s ID (for drivers and for all genes and the relations between them) and it’s distribution, a histogram for the inequalities coefficients for drivers and for all genes, and the order of the gene’s dysfunction.
```

Results of simulations
---
The script allows to make repeating simulation/replicas to get statistical data:
```
	•	REPEAT.zip - zip arhive with the script to calculate			
	
	•	Results of simulation.pdf - the file with a plots and graphs of results of simulations for colorectal cancer, for one simulation and for one hundred simulations. 

	•	last_step.txt.aa.zip and last_step.txt.ab.zip - the two parties (zip arhive) of the file "last_step.txt", which includes of the data for the last step of the one hundred simulations. 
```

