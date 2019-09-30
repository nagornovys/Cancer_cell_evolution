Tests for tugHall
====================

---
## File List

* _Tests.zip_ - archive with tests for program code. Each test is in the separate directory and has a note with explanation of details of test. Some tests need to change a code a bit, you can find information about change in the note and also in the code, for this, please, find word "test".
* _presentation_of_tests.pdf_ - presentation file with results of all tests. You can find all tests and results in a brief manner.

## Input files
* _gene_cds2.txt_ - the information about _Names of genes_, _length of CDS_ (CoDing Sequence - the coding region of a gene), _relation to hallmarks_, _oncogene or suppressor_ (o or s) and _relative weights for hallmarks_.
* _cellinit.txt_ - the list of _initial cells_ with _ID_ with a _list of mutated genes_. 

## Output files in TESTS
* _log.txt_ - the file with information about input parameters and names of input files.
* _geneout.txt_ - the file with information about _Names of genes_, _length of CDS_, _relation to hallmarks_, _oncogene or suppressor_ and _real weights for hallmarks_, which are used in the simulation.
* _Order_of_gene_disfunction.txt_ - the file with information of order of gene disfunction for clones at last time step.
* _cellout.txt_ - the file with **all** output data during the simulation for each time step.


## Script to make 100 trials for same simulation

This script is a tests for repeating simulation and checking that statistical data are around some average value.
The script allows to make repeating simulation/replicas to get statistical data:
* _REPEAT.zip_ - zip arhive with the script to calculate
* _Results of simulation.pdf_ - the file with a plots and graphs of results of simulations for colorectal cancer, for one simulation and for one hundred simulations.
* _ast_step.txt.aa.zip_ and _last_step.txt.ab.zip_ - the two parties (zip arhive) of the file "_last_step.txt_", which includes of them both.

How to RUN the simulation
--
In order to make the simulation, please, follow next procedure:

1. Open the script _CancerProgress_Final_Version.R_ in R Stiduo, for example.
2. To check or/and change the parameters of simulation in the end of file and in the input files.
3. Run the simulation using _"model()"_ function. 
4. To see data, please, see _cellout.txt_ file. 
5. To analize results of simulation and to plot them, please, use _ineq.R_ and _Analyze.R_ scripts.


Funding
---
This work was supported by CREST-JST (14531766); MEXT (15K06916); and AMED (16ck0106013h0003).

Acknowledgements
---
We thank Asmaa Elzawahry, Yusuke Suenaga, Sana Yokoi, Yoshitaka Hippo, Atsushi Niida, and Daniel A. Vasco for useful suggestions.
