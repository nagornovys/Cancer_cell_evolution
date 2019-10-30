tugHall
====================

[![License](https://img.shields.io/badge/License-GPLv3-orange.svg)](https://github.com/nagornovys/Cancer_cell_evolution/blob/master/LICENSE)


**tugHall** _(**tu**mor **g**ene-**Hall**mark)_ is a simulator of a cancer-cell evolution model, wherein gene mutations are linked to the tumor cell behaviors that are influenced by the hallmarks of cancer.

This is a script in _**R**_ to simulate the cancer cell evolution in the framework of the model proposed by _**prof. Mamoru Kato**_,
_Head of Bioinformatics Department, Research Institute, Nation Cancer Center, Tokyo, JAPAN_.

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
The flood of recent cancer genomic data requires a coherent model that can sort out the findings to systematically explain clonal evolution and the resultant intra-tumor heterogeneity (ITH). Here, we present a new mathematical model designed to computationally simulate the evolution of cancer cells. The model connects the well-known hallmarks of cancer with the specific mutational states of tumor-related genes. The cell behavior phenotypes are stochastically determined and the hallmarks probabilistically interfere with the phenotypic probabilities. In turn, the hallmark variables depend on the mutational states of tumor-related genes. Thus, our software can deepen our understanding of cancer-cell evolution and generation of ITH.

Content of package
---

* _User-Guide-tugHall.Rmd_ and _User-Guide-tugHall.html_ are the user guide to install, run and use tugHall simulator in R Markdown and html formats. 
* _User-Guide-analysis.Rmd_ and _User-Guide-analysis.html_ are the user guide to analyze results of single simulation in R Markdown and html formats.
* _/tugHall/_ is a directory with program code, input and output data.
* _/TESTS/_ is a directory with tests, please, see README_TESTS.md file for explanation.
* _/ABC/_ is a directory with a simple example of 10000 trials of simulation and Approximate Bayesian computation. For details, please, see README_ABC.md. 

Funding
---
This work was supported by CREST-JST (14531766); MEXT (15K06916); and AMED (16ck0106013h0003).
