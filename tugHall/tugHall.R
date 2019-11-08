
##########################################################################################
###  The simulation uses the functions and classes in the "Code/tugHall_functions.R" 

library(stringr)

source(file = "Code/tugHall_functions.R")


## Create folders:  /Input, /Output and /Figures 

mainDir <- getwd()
subDir <- "Output"
if (! file.exists(subDir)){  dir.create(file.path(mainDir, subDir)) }

subDir <- "Input"
if (! file.exists(subDir)){  dir.create(file.path(mainDir, subDir)) }

subDir <- "Figures"
if (! file.exists(subDir)){  dir.create(file.path(mainDir, subDir)) }


##########################################################################################
### Files to output and input data

genefile <- 'Input/gene_cds2.txt'    # gene file 
cellfile <- 'Input/cellinit.txt'     # initial Cells 

### Output files
geneoutfile <- 'Output/geneout.txt'  # Gene Out file with Hallmarks 
celloutfile <- 'Output/cellout.txt'  # output information of simulation
logoutfile <-  'Output/log.txt'      # log file to save the input information of simulation - "log.txt"
### Output/Weights.txt               # file with gene weights for hallmarks


##########################################################################################
# Probabilities of processes

E0 <<- 0.5E-3       # parameter in the division probability  
F0 <<- 10        # parameter in the division probability  
m <<-  1E-5       # mutation probability  
uo <<- 0.5        # oncogene mutation probability  
us <<- 0.5        # suppressor mutation probability  
s <<-  10         # parameter in the sigmoid function  
k <<-  0.1        # Environmental death probability  

### Additional parameters of simulation
censore_n <<- 30000       # Max cell number where the program forcibly stops
censore_t <<- 100         # Max time where the program forcibly stops

##########################################################################################
# if you have a new format of gene file, please, use change of columns function like: 
# genefile <- changeCol(genefile)

### Making of the input file for initial cells
x <- 1:1000
xz <- data.frame(x,"")
write.table(xz,file = cellfile, col.names = FALSE,sep = "\t",row.names = FALSE)


##########################################################################################
### Simulation of the cancer cell evolution:

model(genefile, cellfile, geneoutfile, celloutfile, logoutfile, E0, F0, m, uo, us, s, k, censore_n, censore_t)



##########################################################################################
#### Analysis of the output data:

# Note: if output files have no data to plot Code/Analysis.R produces errors during plotting

source("Code/Analysis.R")



##########################################################################################

# In order to make report, please, use USER-GUIDE.Rmd to show results of simulation

# In order to improve the output plot, please, use Code/Functions.R and Code/Analysis.R scripts

