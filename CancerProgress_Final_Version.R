# Cell class
cell <- setRefClass(
    # name of the class
    Class = "Cell",
    # Field
    fields = list(
        id = "numeric",          # identificator
        parent = "numeric",      # parent ID (for first - 0)
        c = "numeric",           # split counter
        d = "numeric",           # probability of division
        i = "numeric",           # probability of hayflick limit
        m = "numeric",           # probability that gene normal function is destroyed due to epigenome abnormality
        a = "numeric",           # probability of apoptosis
        s = "numeric",           # coefficient in the of apoptosis
        k = "numeric",           # probability of cell death by environment
        E = "numeric",           # coefficient of friction term against to the split probability,
                                 # coefficient for determination the max number of cells that can exist in the primary tumor (Nmax = 1/E)
        Nmax = "numeric",        # the max number of cells that can exist in the primary tumor (Nmax = 1/E)
        im = "numeric",          # invasion/ metastasis probability
        Ha = "numeric",          # apoptosis probability difference (apoptosis gene x weight) 
        Him = "numeric",         # invasion/ metastasis probability difference (invasion/ metastasis gene x weight) 
        Hi = "numeric",          # mitotic restriction probability difference (immortalized gene x weight)
        Hd = "numeric",          # divide probability difference (cancer gene x weight)
        Hb = "numeric",          # friction coefficient (angiogenic gene x weight)
        gene = "numeric",        # flag for cancer gene function deficit 
        pasgene = "numeric",     # flag for cancer gene as passenger dysfunction 
        posdriver = "character", # position of cancer gene damage (function deficit)
        pospasngr = "character", # position of cancer gene damage (maintenance of function)
        mutden = "numeric",      # gene mutation density
        invasion = "logical",    # Wetting/Displacement flag:    TRUE: Wetting/Displacement      FALSE: Limited pattern
        birthday = "numeric"     # time step of birth of cell
    ),

    # Method
    methods = list(
        # Initialize
        initialize = function(gene_size, id=1, parent=0, c=0, d=0, i=1, m=5*10E-8,
                              mutden=0, a=0, k=runif(1), E=0.1, Nmax=Nmax, gene=NULL, pasgene=NULL,
                              posdriver=NULL, pospasngr=NULL, invasion=FALSE, s=10,birthday=0) {
            id <<- id
            parent <<- parent
            c <<- c
            d <<- d
            i <<- i
            m <<- m
            s <<- s
            birthday <<- birthday
            mutden <<- mutden
            if (is.null(a)) {
                a <<- 1/(1+exp(-s*(mutden - 0.5)))
            } else {
                a <<- a
            }
            k <<- k
            E <<- E
            Nmax <<- 1.0 / E
            im <<- 0
            Ha <<- 0
            Him <<- 0
            Hi <<- 0
            Hd <<- 0
            Hb <<- 0

            if (is.null(gene)) {
                gene <<- rep(0, gene_size)
            } else {
                gene <<- gene
            }
            if (is.null(pasgene)) {
              pasgene <<- rep(0, gene_size)
            } else {
              pasgene <<- pasgene
            }
            if (is.null(posdriver)) {
                posdriver <<- rep("", gene_size)
            } else {
                posdriver <<- posdriver
            }
            if (is.null(pospasngr)) {
                pospasngr <<- rep("", gene_size)
            } else {
                pospasngr <<- pospasngr
            }
            invasion <<- invasion
        },
        # Apoptosis
        calcApoptosis = function() {
#            if (mutden <= sum(gene)/length(gene)) {
#                a1 = 1/(1+exp(s*(mutden - 0.5)))
#                mutden <<- sum(gene)/length(gene)
#                a2 = 1/(1+exp(s*(mutden - 0.5)))
#                a <<- a - (a1 - a2)

          mutden <<- sum(gene)/length(gene)
          a <<- 1/(1+exp(-1*s*(mutden - 0.5)))
                if (a < 0) {
                    a <<- 0
                }

#            }
        },
        # Aggregate
        calcMutden = function() {
            mutden <<- sum(gene)/length(gene)
        }
    )
)

# Environ class
environ <- setRefClass(
    # the class name
    Class = "Environ",

    # Fields
    fields = list(
        T = "numeric",           # time counter
        N = "numeric",           # localized cell number
        M = "numeric",           # number of infiltrting / metastatic cells
        F = "numeric",           # a coeffitient (Nmax = F/E) that determines the maximal number of cells 
                                 # that can exist in the primary tumor when the hallmark is engraved  
        c = "numeric",           # average number of divisions 
        d = "numeric",           # mean value of spliting probability
        i = "numeric",           # average value of spliting probability
        a = "numeric",           # average value of apoptosis probability
        k = "numeric",           # average probability of cell death
        E = "numeric",           # average value of coefficients of friction term proportional to N, for splitting probability
        Nmax = "numeric",        # Maximal number of cells that can exist 
        im = "numeric",          # average value of invasion / metastasis probability
        Ha = "numeric",
        Him = "numeric",
        Hi = "numeric",
        Hd = "numeric",
        Hb = "numeric",
        type = "numeric",        # invasion / metastatic ratio
        gene = "numeric",        # cancer gene damage rate
        posdriver = "character", # cancer gene damage position (function deficit)
        pospasngr = "character", # cancer gene damage position (maintaince of function)
        mutden = "numeric",      # average mutation rate
        last_id = "numeric"
    ),

    # Methods
    methods = list(
        # Initialize
        initialize = function(F0) {
            T <<- 0
            N <<- 0
            M <<- 0
            F <<- F0
        }
    )
)


# Cancer gene
oncogene <- setRefClass(
    # name of the class
    Class = "OncoGene",

    # Fields
    fields = list(
        name = "character",   # Cancer gene name list
        onsp = "character",   # oncogene/suppressor indicator
        cds = "numeric",      # cancer gene CDS base number list
        len = "numeric"       # number of cancer genes
    ),

    # Methods
    methods = list(
        # read the configuration file
        read = function(file) {
            data = read.table(file, sep="\t")
            name0 = NULL
            onsp0 = NULL
            cds0 = NULL
            for (i in 1:nrow(data)) {
                name <<- as.character(data[i, 1])
                if (!is.element(name, name0)) {
                    name0 = c(name0, name)
                    type = as.character(data[i, 4])
                    if (type == "?") {
                        if (runif(1) > 0.5) {
                            type = "o"
                        } else {
                            type = "s"
                        }
                    }
                    onsp0 = c(onsp0, type)
                    cds0 = c(cds0, as.numeric(as.character(data[i, 2])))
                }
            }
            name <<- name0
            onsp <<- onsp0
            cds <<- cds0
            len <<- length(name0)
        }
    )
)


hallmark <- setRefClass(
    # 
    Class = "HallMark",

    # 
    fields = list(
        Ha = "numeric",       # (Evading apoptosis)
        Hi = "numeric",       # (immortalization limit)
        Hd = "numeric",       # (Insensitivity to anti-growth signals || Self-sufficiency in growth signals)
        Hb = "numeric",       # (Sustained angiogenesis)
        Him = "numeric",      # (Tissue invasion & metastasis)
        Ha_w = "numeric",     # 
        Hi_w = "numeric",     # 
        Hd_w = "numeric",     # 
        Hb_w = "numeric",     # 
        Him_w = "numeric",    # 
        notHa = "numeric"
    ),

    # 
    methods = list(
        # 
        read = function(file, names) {
            data <- read.table(file, sep="\t")
            Ha0 = NULL
            Hi0 = NULL
            Hd0 = NULL
            Hb0 = NULL
            Him0 = NULL
            Ha0_w = NULL
            Hi0_w = NULL
            Hd0_w = NULL
            Hb0_w = NULL
            Him0_w = NULL
            w_flg = FALSE
            if (ncol(data) >= 5) {
                w_flg = TRUE
            }
            # Acquire gene name and Hallmark coefficient by function from definition file
            for (i in 1:nrow(data)) {
                if (data[i, 3] == "apoptosis") {
                    Ha0 = c(Ha0, as.character(data[i, 1]))
                    if (w_flg) {
                        Ha0_w = c(Ha0_w, as.numeric(as.character(data[i, 5])))
                    }
                } else if (data[i, 3] == "immortalization") {
                    Hi0 = c(Hi0, as.character(data[i, 1]))
                    if (w_flg) {
                        Hi0_w = c(Hi0_w, as.numeric(as.character(data[i, 5])))
                    }
                } else if (data[i, 3] == "anti-growth" | data[i, 3] == "growth") {
                    Hd0 = c(Hd0, as.character(data[i, 1]))
                    if (w_flg) {
                        Hd0_w = c(Hd0_w, as.numeric(as.character(data[i, 5])))
                    }
                } else if (data[i, 3] == "angiogenesis") {
                    Hb0 = c(Hb0, as.character(data[i, 1]))
                    if (w_flg) {
                        Hb0_w = c(Hb0_w, as.numeric(as.character(data[i, 5])))
                    }
                } else if (data[i, 3] == "invasion") {
                    Him0 = c(Him0, as.character(data[i, 1]))
                    if (w_flg) {
                        Him0_w = c(Him0_w, as.numeric(as.character(data[i, 5])))
                    }
                }
            }
            Ha <<- match(Ha0, names)
            notHa <<- setdiff(seq(1,length(names)),Ha)
            Hi <<- match(Hi0, names)
            Hd <<- match(Hd0, names)
            Hb <<- match(Hb0, names)
            Him <<- match(Him0, names)

            # if there is no Hallmark coefficient then generate a Hallmark coefficient as a random number - beta distribution
            if (!w_flg) {
                if (length(Ha) > 0) {
                    Ha_rnd = 1:length(Ha)
                } else {
                    Ha_rnd = c()
                }
                total0 = length(Ha)
                if (length(Hi) > 0) {
                    Hi_rnd = (total0 + 1):(total0 + length(Hi))
                } else {
                    Hi_rnd = c()
                }
                total0 = total0 + length(Hi)
                if (length(Hd) > 0) {
                    Hd_rnd = (total0 + 1):(total0 + length(Hd))
                } else {
                    Hd_rnd = c()
                }
                total0 = total0 + length(Hd)
                if (length(Hb) > 0) {
                    Hb_rnd = (total0 + 1):(total0 + length(Hb))
                } else {
                    Hb_rnd = c()
                }
                total0 = total0 + length(Hb)
                if (length(Him) > 0) {
                    Him_rnd = (total0 + 1):(total0 + length(Him))
                } else {
                    Him_rnd = c()
                }
                total = total0 + length(Him)
                # random = runif(total)
                random = rbeta(total, 0.01, 1)
                Ha0_w = random[Ha_rnd]
                Hi0_w = random[Hi_rnd]
                Hd0_w = random[Hd_rnd]
                Hb0_w = random[Hb_rnd]
                Him0_w = random[Him_rnd]
            }
            # Total by genetic mode 
            Ha_sum = sum(Ha0_w)
            Hi_sum = sum(Hi0_w)
            Hd_sum = sum(Hd0_w)
            Hb_sum = sum(Hb0_w)
            Him_sum = sum(Him0_w)
            Ha_w <<- Ha0_w/Ha_sum
            Hi_w <<- Hi0_w/Hi_sum
            Hd_w <<- Hd0_w/Hd_sum
            Hb_w <<- Hb0_w/Hb_sum
            Him_w <<- Him0_w/Him_sum
        },
        # Change the cell variables
        # mode = 2 Corresponding (Hallmark) Gene Mode
        updateCell = function(cell1, F) {
          # Apoptosis
            cell1$calcApoptosis()
            cell1$Ha = sum(cell1$gene[Ha]*Ha_w)
            cell1$a = cell1$a - cell1$Ha
            if (cell1$a < 0) {
                cell1$a = 0
            }
            # Not dead - Immortalized
            cell1$Hi = sum(cell1$gene[Hi]*Hi_w)
            cell1$i = 1 - cell1$Hi
            if (cell1$i < 0) {
                cell1$i = 0
            }
            # Angiogenesis
            cell1$Hb = sum(cell1$gene[Hb]*Hb_w)
            
            cell1$E = E0 / (1 + F * cell1$Hb)
            cell1$Nmax = 1.0 / cell1$E
            
            # Cancer gene, tumor suppressor gene
            cell1$Hd = sum(cell1$gene[Hd]*Hd_w)
            cell1$Him = sum(cell1$gene[Him]*Him_w)
            
            cell1$d = cell1$Hd
                if (cell1$d > 1) {cell1$d = 1}
            if (!cell1$invasion) {
              cell1$d = cell1$d - cell1$E * env$N
              if (cell1$d < 0) {cell1$d = 0}
            }

                # Invasion metastasis
                if (!cell1$invasion) {
                cell1$im = cell1$Him
            } else {
                cell1$im = 1
            }
          },
        # Change the environment variables
        updateEnviron = function(env, cells) {
            env$T = env$T + 1
            sum_cell(env, cells)
            env$M = ceiling(length(cells) * env$type)
            env$N = length(cells)  - env$M
        }
    )
)

# Function to update Hallmark and variable after division or under initialization
update_Hallmarks <- function(cell1) {
  # Hallmark
  hall$updateCell(cell1, env$F)
}

# trial function
trial <- function(cell1) {
    gtotal = onco$len       # length of onco 
    # random numbers for gene mutation
    # random numbers for mutation position
    random = runif(2*gtotal + 5)
    # Hallmark
#    hall$updateCell(cell1, env$F)
    # return values
    # 0: death 1:raw 2:split
    # trial for Environmental death of cell 
    if (cell1$k >= random[2*gtotal + 1]) {
        return(0)
    }
    # Apoptosis trial
    if (cell1$a >= random[2*gtotal + 2]) {
        return(0)
    }

    # invasion / metastasis trial 
    if (!cell1$invasion) {
      if (cell1$im > 0) {
        if (cell1$im >= random[gtotal + 5]) {
          cell1$invasion = TRUE
        } else {
          return(0)
        }
      }
    }    
    
    # Fragmentation restriction trial
    divi = TRUE
    if (cell1$c > 50) {
        if (cell1$i >= random[2*gtotal + 3]) {
            divi = FALSE
        }
    }
    # Devide trial 
    status = 1
    if (divi) {
        # if it is NOT invasion / metastatic type 
        d1 = cell1$d
        # dividing 
        if (d1 >= random[2*gtotal + 4]) {
            cell1$c = cell1$c + 1
            status = 2
        }
    }

    return(status)
}

# mutagenesis trial
trial_mutagenesis <- function(cell1) {
  gtotal = onco$len       # length of onco 
  # random numbers for gene mutation
  # random numbers for mutation position
  random = runif(2*gtotal + 5)
#  print(c(cell1$id,m,cell1$m,onco$cds,gtotal))

  
# mut1 = ifelse(cell1$m*onco$cds >= random[1:gtotal], 1, 0)
mut1 = ifelse(m*onco$cds >= random[1:gtotal], 1, 0)
posg = ceiling(random[(gtotal + 1):(2*gtotal)]*onco$cds)[mut1==1]
random2 = runif(length(mut1[mut1==1]))

# print(c(cell1$id,onco$onsp))

mut2 = ifelse((onco$onsp[mut1 == 1] == 'o' & uo >= random2) |
                (onco$onsp[mut1 == 1] == 's' & us >= random2), 1, 0)
cell1$gene[mut1==1] = ifelse(cell1$gene[mut1==1] == 1 | mut2 == 1, 1, 0)
cell1$pasgene[mut1==1] = ifelse(cell1$pasgene[mut1==1] == 1 | mut2 == 0, 1, 0)
cell1$posdriver[mut1==1][mut2==1] = ifelse(cell1$posdriver[mut1==1][mut2==1] == '',
                                           paste(posg[mut2==1],env$T,sep = ":"), paste(cell1$posdriver[mut1==1][mut2==1], paste(posg[mut2==1],env$T,sep = ":"), sep=","))
cell1$pospasngr[mut1==1][mut2==0] = ifelse(cell1$pospasngr[mut1==1][mut2==0] == '',
                                           paste(posg[mut2==0],env$T,sep = ":"), paste(cell1$pospasngr[mut1==1][mut2==0], paste(posg[mut2==0],env$T,sep = ":"), sep=","))

}


# to make one copy for cell1 in cell_init function
cell_copy <- function(cell1) {
  env$last_id = env$last_id + 1
      return(cell$new(id=env$last_id, parent=cell1$id, c=cell1$c, d=cell1$d, i=cell1$i, mutden=cell1$mutden,
                      a=cell1$a, k=cell1$k, E=cell1$E, Nmax=cell1$Nmax, gene=cell1$gene, pasgene=cell1$pasgene, posdriver=cell1$posdriver,
                      pospasngr=cell1$pospasngr, invasion=cell1$invasion, s=cell1$s, birthday=env$T))
}


# aggregate
sum_cell <- function(env, cells) {
    if (length(cells) > 0) {
        avg = apply(matrix(unlist(lapply(cells, sum_mutation)),ncol=length(cells)),1,sum)/length(cells)
        env$c = avg[1]
        env$d = avg[2]
        env$i = avg[3]
        env$a = avg[4]
        env$k = avg[5]
        env$E = avg[6]
        env$Nmax = avg[7]
        env$im = avg[8]
        env$Ha = avg[9]
        env$Him = avg[10]
        env$Hi = avg[11]
        env$Hb = avg[12]
        env$Hd = avg[13]
        env$type = avg[14]
        env$mutden = avg[15]
        env$gene = avg[16:length(avg)]
    } else {
        env$c = 0
        env$d = 0
        env$i = 0
        env$a = 0
        env$k = 0
        env$E = 0
        env$Nmax = 0
        env$im = 0
        env$Ha = 0
        env$Him = 0
        env$Hi = 0
        env$Hb = 0
        env$Hd = 0
        env$type = 0
        env$gene = rep(0, length(onco$name))
    }
    env$posdriver = rep("", length(onco$name))
    env$pospasngr = rep("", length(onco$name))
}


sum_mutation <- function(cell1) {
    return(c(cell1$c, cell1$d, cell1$i, cell1$a, cell1$k, cell1$E,
             cell1$Nmax, cell1$im, cell1$Ha, cell1$Him, cell1$Hi,
             cell1$Hb, cell1$Hd, ifelse(cell1$invasion,1,0), cell1$mutden, cell1$gene))
}

# write log file
write_log <- function(genefile, cellfile, geneoutfile, celloutfile, logoutfile, E0, F0, m, uo, us, s, k, censore_n, censore_t) {
    data <- c("genefile", "cellfile", "geneoutfile", "celloutfile", "logoutfile",
              "E", "F", "m", "uo", "us", "s", "k", "censore_n", "censore_t")
    data <- rbind(data, c(genefile, cellfile, geneoutfile, celloutfile, logoutfile,
                          E0, F0, m, uo, us, s, k, censore_n, censore_t))
    write(data, logoutfile, ncolumn=2, sep="\t")
}
write_geneout <- function(outfile, hall) {
    data <- c(onco$name[hall$Ha], onco$name[hall$Hi], onco$name[hall$Hd], onco$name[hall$Hb], onco$name[hall$Him])
    data <- rbind(data, c(rep("apoptosis", length(onco$name[hall$Ha])),
                          rep("immortalization", length(onco$name[hall$Hi])),
                          rep("growth|anti-growth", length(onco$name[hall$Hd])),
                          rep("angiogenesis", length(onco$name[hall$Hb])),
                          rep("invasion", length(onco$name[hall$Him]))))
    data <- rbind(data, c(hall$Ha_w, hall$Hi_w, hall$Hd_w, hall$Hb_w, hall$Him_w))
    data <- rbind(data, c(onco$onsp[hall$Ha], onco$onsp[hall$Hi], onco$onsp[hall$Hd], onco$onsp[hall$Hb], onco$onsp[hall$Him]))
    write(data, outfile, ncolumn=4, sep="\t")
}

write_header <- function(outfile, env) {
    header <- c('Trial', 'Time', 'AvgOrIndx', 'ID', 'ParentID:Birthday', 'c\'', 'd\'', 'i\'', 'im\'', 'a\'',
                'k\'', 'E\'', 'N', 'Nmax\'', 'M', 'Ha', 'Him', 'Hi', 'Hd', 'Hb', 'type', 'mut_den',
                paste("PosDriver:", onco$name, sep=""), paste("PosPasngr:", onco$name, sep=""), 'Clone number', 'Passengers Clone number', 'Mix Clone number')
    write(header, outfile, append=FALSE, ncolumn=length(header), sep="\t")
}

write_cellout <- function(outfile, env, cells, isFirst) {
    data <- c(env$T, 'avg', '-', '-', env$c, env$d, env$i, env$im, env$a, env$k, env$E, env$N,
              env$Nmax, env$M, env$Ha, env$Him, env$Hi, env$Hd, env$Hb, env$type, env$mutden,
              env$posdriver, env$pospasngr, '-', '-', '-')
    write(data, outfile, append=TRUE, ncolumn=length(data), sep="\t")
    if (length(cells) > 0 & isFirst) {
        for (i in 1:length(cells)) {
            cell1 = cells[[i]]
            vc <- 0
            for (j in 1:length(cell1$gene)) {vc<- vc + 2^(length(cell1$gene)-j) * cell1$gene[j]}
            pasvc <- 0
            for (j in 1:length(cell1$pasgene)) {pasvc<- pasvc + 2^(length(cell1$pasgene)-j) * cell1$pasgene[j]}
            mixvc <- 0
            for (j in 1:length(cell1$pasgene)) {mixvc<- mixvc + 2 * 4^(length(cell1$pasgene)-j) * cell1$pasgene[j] + 4^(length(cell1$gene)-j) * cell1$gene[j]}

            data <- c(env$T, i, cell1$id, paste(cell1$parent,cell1$birthday,sep = ":"), cell1$c, cell1$d, cell1$i, cell1$im, cell1$a,
                      cell1$k, cell1$E, env$N, cell1$Nmax, env$M,
                      cell1$Ha, cell1$Him, cell1$Hi, cell1$Hd, cell1$Hb, ifelse(cell1$invasion,1,0), cell1$mutden,
                      cell1$posdriver, cell1$pospasngr, vc, pasvc, mixvc)
            write(data, outfile, append=TRUE, ncolumn=length(data), sep="\t")
        }
    }
}

# To write last step
write_last <- function(outfile, env, cells, isFirst,rp) {
  if (length(cells) > 0 & isFirst) {
    for (i in 1:length(cells)) {
      cell1 = cells[[i]]
      vc <- 0
      for (j in 1:length(cell1$gene)) {vc<- vc + 2^(length(cell1$gene)-j) * cell1$gene[j]}
      pasvc <- 0
      for (j in 1:length(cell1$pasgene)) {pasvc<- pasvc + 2^(length(cell1$pasgene)-j) * cell1$pasgene[j]}
      mixvc <- 0
      for (j in 1:length(cell1$pasgene)) {mixvc<- mixvc + 2 * 4^(length(cell1$pasgene)-j) * cell1$pasgene[j] + 4^(length(cell1$gene)-j) * cell1$gene[j]}
      
      data <- c(rp, env$T, i, cell1$id, paste(cell1$parent,cell1$birthday,sep = ":"), cell1$c, cell1$d, cell1$i, cell1$im, cell1$a,
                cell1$k, cell1$E, env$N, cell1$Nmax, env$M,
                cell1$Ha, cell1$Him, cell1$Hi, cell1$Hd, cell1$Hb, ifelse(cell1$invasion,1,0), cell1$mutden,
                cell1$posdriver, cell1$pospasngr, vc, pasvc, mixvc)
      write(data, outfile, append=TRUE, ncolumn=length(data), sep="\t")
    }
  }
}

# To write only average values

write_average <- function(outfile, env, cells, isFirst,rp) {
  data <- c(rp,env$T, 'avg', '-', '-', env$c, env$d, env$i, env$im, env$a, env$k, env$E, env$N,
            env$Nmax, env$M, env$Ha, env$Him, env$Hi, env$Hd, env$Hb, env$type, env$mutden,
            env$posdriver, env$pospasngr, '-', '-', '-')
  write(data, outfile, append=TRUE, ncolumn=length(data), sep="\t")
 }




# initial cell setting
init_cells <- function(cellfile, cell1) {
    mpos <- regexpr("\\.", cellfile)[1]
    if (mpos != -1) {
        name <- substr(cellfile, 1, mpos - 1)
    } else {
        name <- cellfile
    }
    cells = NULL
    n <- as.numeric(name)
    if (!is.na(n) && is.numeric(n)) {
        factor = n / sum(cell1$m*onco$cds)
        f2 = 1.0
        while (TRUE) {
            if (sum(floor(cell1$m*onco$cds*factor*f2 + 0.5)) >= n) {
                break
            }
            f2 = f2 + 0.1
        }
        nums = floor(cell1$m*onco$cds*factor*f2 + 0.5)
        cells = NULL
        for (i in 1:n) {
            cells = c(cells, cell_copy(cell1))
        }
        pos = 0
        for (i in 1:length(nums)) {
            if (nums[i] > 0) {
                for (j in 1:nums[i]) {
                    if (pos + j <= n) {
                        cells[[pos + j]]$gene[i] = 1
                    }
                }
                pos = pos + nums[i]
            }
        }
    } else {
        data = read.table(cellfile, sep="\t")
        n <- nrow(data)
        
        for (i in 1:n) {
            cell2 = cell_copy(cell1)
            p <- match(onco$name, str_trim(strsplit(as.character(data[i,2]),",")[[1]]))
            cell2$gene[seq(1,length(onco$name))[!is.na(p)]] = 1
            cells = c(cells, cell2)
        }
    }
    for (i in 1:n) {
        cells[[i]]$id = i
        cells[[i]]$parent = 0
        cells[[i]]$birthday = 0
        cells[[i]]$posdriver = ifelse(cells[[i]]$gene == 1,
                                      paste(ceiling(runif(onco$len)*onco$cds),"0",sep = ":"),
                                      cells[[i]]$posdriver)
        cells[[i]]$calcMutden()
        cells[[i]]$calcApoptosis()
    }
    env$last_id = n
    return(as.list(cells))
}

# Genetic recombination
read_w <- function(file) {
    if (!is.null(file) & !is.na(file)) {
        data <- read.table(file, sep="\t")
        w <- NULL
        r = 1
        start = 1
        while (TRUE) {
            window = ceiling(runif(1)*ncol(data))
            if (start + window > ncol(data)) {
                window = ncol(data) - start
            }
            w = c(w, data[r,start:(start+window)])
            if (start+window == ncol(data)) {
                break
            }
            start = start + window
            r = ifelse(r==1,2,1)                  
        }
        hall$setW(w)
    }
}

model <- function(genefile, cellfile, geneoutfile, celloutfile, logoutfile, E0, F0, m, uo, us, s, k, censore_n, censore_t, rp) {
    write_log(genefile, cellfile, geneoutfile, celloutfile, logoutfile, E0, F0, m, uo, us, s, k, censore_n, censore_t)   # write input parameters
    onco = oncogene$new()        # make the vector onco about the hallmarks
    onco$read(genefile)          # read the input info to the onco from genefile - 'gene_cds2.txt'
    hall = hallmark$new()        # make a vector hall with hallmarks parameters
    hall$read(genefile, onco$name)     # read from the genefile - 'gene_cds2.txt'
    env = environ$new(F0)               # new vector for average values of cells
    assign("env", env, env=.GlobalEnv)
    assign("onco", onco, env=.GlobalEnv)
    assign("hall", hall, env=.GlobalEnv)
    assign("uo", uo, env=.GlobalEnv)
    assign("us", us, env=.GlobalEnv)
    cell1 = cell$new(gene_size=length(onco$cds),
                     m=m, s=s, k=k, E=E0)              # cell1  -  empty object of cell
    cells = init_cells(cellfile, cell1)                # cells - cells with hallmarks from cellfile - cellinit.txt - initial cells  
    write_geneout(geneoutfile, hall)                   # write the geneout.txt file with initial hallmarks 
    if (rp == 1) write_header(celloutfile, env)  
    if (rp == 1) write_header("last_steps.txt", env)        # 
    lapply(cells,update_Hallmarks)                     # to calculate the Hallmarks and probabilities for initial cells
    hall$updateEnviron(env, cells)                     # make averaging for cells and first step (T=T+1)
    isFirst = TRUE
 #   write_cellout("initial_cells.txt", env, cells, isFirst)     #  write initial cells

    while(length(cells) > 0 && censore_n > (env$N + env$M) && env$T < censore_t) {

            ret = unlist(lapply(cells, trial))               # The mark of cells during application of the trial !!! WOTHOUT changing of cells
            cells = c(cells, lapply(cells[ret==2], cell_copy))  # add cells after division (ret=2)            
            ret=c(ret,ret[ret==2])
            lapply(cells[ret==2], trial_mutagenesis)         # apply the mutagenesis trial to the parent and child cells independently
            cells = c(cells[ret>=1])                         # delete cells with ret=0           
            ret=c(ret[ret>=1])
            lapply(cells,update_Hallmarks)                  # to calculate the Hallmarks and probabilities for initial cells
            hall$updateEnviron(env, cells)                   #  make averaging for cells, Hallmarks etc and increase step (T=T+1)  
#            write_cellout(celloutfile, env, cells, isFirst)
            write_average(celloutfile, env, cells, isFirst,rp)
    }
    write_last("last_steps.txt", env, cells, isFirst, rp)
}

# Exchange the 3rd and 4th coloumns in the genefile 
 changeCol <- function(genefile) {
  exchange = read.table(file =genefile, header = FALSE)
  Vec1 = exchange$V1
  Vec2 = exchange$V2
  Vec3 = exchange$V3
  Vec4 = exchange$V4
  Vec5 = exchange$V5
  exch_new = data.frame(Vec1,Vec2,Vec4,Vec3,Vec5)
  genefileNew = substr(genefile,1,4)
  write.table(exch_new, file = genefileNew, append = FALSE, quote = TRUE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = FALSE, qmethod = c("escape", "double"),
              fileEncoding = "")
  
 return(genefileNew) 
}

#args <- commandArgs(trailingOnly=TRUE)
#args <- c('gene_cds2.txt', 'cellinit.txt',
#          'geneout.txt', 'cellout.txt',
#          'log.txt', '1E-30', '1E30', '0.5', '0.5', '0.5',
#          '10', '0', '100000', '1000', '10')

library(stringr)
 
genefile <- 'gene_cds2.txt'     #args[1]                # gene file - "gene_cds2.txt"
cellfile <- 'cellinit.txt'      #args[2]                # initial Cells - "cellinit.txt"
geneoutfile <- 'geneout.txt'    #args[3]             # Gene Out file with Hallmarks - "gene out"
celloutfile <- 'cellout.txt'    #args[4]             # output information of simulation - "cellout.txt" 
logoutfile <-  'log.txt'        #args[5]              # log file to save the input information of simulation - "log.txt"
E0 <- 1E-3                     #as.numeric(args[6])          # (1E-1, 1E-2, 1E-3, 1E-4, 1E-5)
F0 <- 1E1                      #as.numeric(args[7])          # (10^1, 10^2, 10^3, 10^6, 10^9, 10^12)
m <-  1E-8                       #as.numeric(args[8])           # mutation probability (1*10e-11, 5*10e-11, 1*10e-10, 5*10e-10, 1*10e-9, 5*10e-9, 1*10e-8, 5*10e-8)
uo <- 0.5                       #as.numeric(args[9])          # oncogene mutation probability
us <- 0.5                       #as.numeric(args[10])         # suppressor mutation probability
s <-  10                        #as.numeric(args[11])          # (10, 100)
k <-  0.27                         #as.numeric(args[12])          # Environmental death probability
censore_n <- 20000               #as.numeric(args[14])  # Max cell number where the program forcibly stops
censore_t <- 45                 #as.numeric(args[15])  # Max time where the program forcibly stops


# if you have a new format of gene file, please, use change of columns function
# genefile <- changeCol(genefile)

setwd("/home/iurii/REPEAT_SIM/Weights/")

for (rp in 1:100) {
  
model(genefile, cellfile, geneoutfile, celloutfile, logoutfile, E0, F0, m, uo, us, s, k, censore_n, censore_t, rp)
}




