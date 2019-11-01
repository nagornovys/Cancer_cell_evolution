tugHall version 1.1
====================

[![License](https://img.shields.io/badge/License-GPLv3-orange.svg)](https://github.com/nagornovys/Cancer_cell_evolution/blob/master/LICENSE)


**tugHall** _(**tu**mor **g**ene-**Hall**mark)_ is cancer-cell evolution model simulator, wherein gene mutations are linked to the tumor cell behaviors, which are influenced by the hallmarks of cancer.

This is an _**R**_-based script to simulate the cancer cell evolution in the framework of the model proposed by _**Prof. Mamoru Kato**_,
_Head of Bioinformatics Department, Research Institute, National Cancer Center, Tokyo, JAPAN_.

Authors and contributor list:
---
_**Iurii Nagornov**_

_**Mamoru Kato**_

_Department of Bioinformatics, Research Institute, National Cancer Center Japan, Tokyo, Japan_

All questions and requests can be sent to inagonov@ncc.go.jp

Project source can be downloaded from websites  
--- 
https://github.com/nagornovys/Cancer_cell_evolution  -  the developing resource

Short description
---
The wide availability of recent cancer genomic data requires a coherent model that can sort out the relevant findings to systematically explain the clonal evolution and resultant intra-tumor heterogeneity (ITH). Here, we present a new mathematical model designed to computationally simulate the evolution of cancer cells. The model connects well-known cancer hallmarks with the specific mutational states of tumor-related genes. The cell behavior phenotypes are stochastically determined and the hallmarks interfere probabilistically with the phenotypic probabilities. In turn, the hallmark variables depend on the mutational states of tumor-related genes. Thus, it is expected our software can be used to deepen our understanding of cancer-cell evolution and generation of ITH.

Content of package
---

* _User-Guide-tugHall.Rmd_ and _User-Guide-tugHall.html_ are the user guides to install, run, and use tugHall simulator in R Markdown and html formats. 
* _User-Guide-analysis.Rmd_ and _User-Guide-analysis.html_ are the user guide to analyze the results of a single simulation in R Markdown and html formats.
* **/tugHall/** is a directory with the program code, input, and output data.
* **/TESTS/** is a directory with tests, please, see README_TESTS.md file for explanation.
* **/ABC/** is a directory with a simple example of 10000 trials of simulation and Approximate Bayesian computation. For details, please, see README_ABC.md. 
