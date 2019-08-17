tugHall
====================

[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://github.com/nagornovys/Cancer_cell_evolution/blob/master/LICENSE)


**tugHall** _(**tu**mor **g**ene-**Hall**mark)_ is a simulator of a cancer-cell evolution model, wherein gene mutations are linked to the tumor cell behaviors that are influenced by the hallmarks of cancer.

This is a script in _**R**_ to simulate the cancer cell evolution in the framework of the model proposed by _**prof. Mamoru Kato**_, 
_Head of Bioinformatics Department, Research Institute, Nation Cancer Center, Tokyo, JAPAN_.

Authors and contributor list: 
--- 
_**Iurii Nagornov**_ 

_**Mamoru Kato**_

_Department of Bioinformatics, Research Institute, National Cancer Center Japan, Tokyo, Japan_

All questions and requests can be sent to inagonov@ncc.go.jp

Contents
--- 
* [Project source](#Project-source)
* [Short description](#Short-description)
* [Content of package](#Content-of-package)
  * [File List](#File-List)
  * [Script to simulate](#Script-to-simulate)
  * [Input files](#Input-files)
  * [Output files](#Output-files)
  - [Scripts to analyse the output file](#Scripts-to-analyse-the-output-file)
* [Results of simulations](#Results-of-simulations)
* [How to RUN the simulation](#How-to-RUN-the-simulation)
* [Funding](#Funding)
* [Acknowledgements](#Acknowledgements)

Project source can be downloaded from websites  
--- 
https://github.com/nagornovys/Cancer_cell_evolution  -  developed resource

https://www.doi.org/10.5281/zenodo.2667073  -  fixed resource in the form of an archive

Short description
---
The flood of recent cancer genomic data requires a coherent model that can sort out the findings to systematically explain clonal evolution and the resultant intra-tumor heterogeneity (ITH). Here, we present a new mathematical model designed to computationally simulate the evolution of cancer cells. The model connects the well-known hallmarks of cancer with the specific mutational states of tumor-related genes. The cell behavior phenotypes are stochastically determined and the hallmarks probabilistically interfere with the phenotypic probabilities. In turn, the hallmark variables depend on the mutational states of tumor-related genes. Thus, our software can deepen our understanding of cancer-cell evolution and generation of ITH.

Content of package
---
### File List

* _Tests.zip_ - archive with tests for program code. Each test is in the separate directory and has a note with explanation of details of test. Some tests need to change a code a bit, you can find information about change in the note and also in the code, for this, please, find word "test".
* _presentation_of_tests.pdf_ - presentation file with results of all tests. You can find all tests and results in a brief manner.
* _Colorectal_cancer.zip_ - the zip-arhive with R script and input files for simulation (please, see input files in the code and see a format of the input files in the presentation). Other script files are only for analysis the results of simulations. The set of genes is included for colorectal cancer as a example. Please, to change the genes - Hallmarks relations to another simulation.

#### Script to simulate
* _CancerProgress_Final_Version.R_ - the script to simulate cell evolution.

#### Input files
* _gene_cds2.txt_ - the information about _Names of genes_, _length of CDS_ (CoDing Sequence - the coding region of a gene), _relation to hallmarks_, _oncogene or suppressor_ (o or s) and _relative weights for hallmarks_.
* _cellinit.txt_ - the list of _initial cells_ with _ID_ with a _list of mutated genes_. 

#### Output files
* _log.txt_ - the file with information about input parameters and names of input files.
* _geneout.txt_ - the file with information about _Names of genes_, _length of CDS_, _relation to hallmarks_, _oncogene or suppressor_ and _real weights for hallmarks_, which are used in the simulation.
* _Order_of_gene_disfunction.txt_ - the file with information of order of gene disfunction for clones at last time step.
* _cellout.txt_ - the file with **all** output data during the simulation for each time step.

#### Scripts to analyse the output file
The zip arhives of simulations include the R scrips to analize the results of the simulations:
* Functions.R (_ineq.R_ + _Analyze.R_) - script/scripts to calculate the functions for Analyze_of_clones.R (or Dyversity of clones.R), so please run it first.
* _Analyze_of_clones.R_ (or _Dyversity_of_clones.R_) - the R script to analize the calculation data, please, run it after Functions.R (_ineq.R_ + _Analyze.R_). The analyzing script allows to calculate next dependences: an evolution of the cell’s number for the normal and mutated cells, an evolution of the Hallmarks and the probabilities, a number of the cells in clones vs clone’s ID (for drivers and for all genes), a number of clones vs the replica’s ID (for drivers and for all genes and the relations between them) and it’s distribution, a histogram for the inequalities coefficients for drivers and for all genes, and the order of the gene’s dysfunction.

### Results of simulations
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
