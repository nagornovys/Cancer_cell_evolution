# Diversity of clones

library(stringr)

library(ape)
library(ggplot2)
library(ggtree)

par(xpd=TRUE, cex.lab=2, lwd = 2, mar = c(5, 5, 5, 5), tcl = 0.5, cex.axis = 1.75,  mgp = c(3, 0.6, 0))


source(file = "Code/tugHall_functions.R")
onco <<- oncogene$new()        # make the vector onco about the hallmarks
onco$read(genefile)          # read the input info to the onco from genefile - 'gene_cds2.txt'
hall <<- hallmark$new()        # make a vector hall with hallmarks parameters
hall$read(genefile, onco$name)     # read from the genefile - 'gene_cds2.txt'

source("Code/Functions.R")
# The function is in the Analyze.R file 
analyze_data(celloutfile)

# time_max - number of time steps 

clones <- matrix(nrow = length(data_flow$Time), ncol = 7)
clones[,1] <- as.integer(as.character(data_flow$Time))
clones[,2] <- as.integer(as.character(data_flow$Clone.number))
clones[,3] <- as.integer(as.character(data_flow$Passengers.Clone.number))
clones[,4] <- as.integer(as.character(data_flow$Mix.Clone.number))

x <- as.character(data_flow$ParentID.Birthday) 
x <- unlist(strsplit(x, split = ":") )
x <- as.integer(x)
clones[,5] <- x[c(TRUE,FALSE)]
clones[,6] <- x[c(FALSE,TRUE)]
clones[,7] <- as.integer(as.character(data_flow$ID))

dimnames(clones) <- list(1:length(data_flow$Time),c("Time","Clone_ID","Clone_ID_pass",
                                                    "Clone_ID_mix","Parent_ID","Birthday","Cell_ID") )


#ord_cl <- 10^(length(onco$name))

#for (i in 1:length(clones[,1])) {
#clones[i,4] <- as.numeric(as.character(data_flow$Clone.number[i])) * ord_cl  + as.numeric(as.character(data_flow$Passengers.Clone.number[i])) 
#}


diversity <- matrix(nrow = time_max, ncol = 2)

for (i in 1:time_max) {
diversity[i,1] <- length(table(clones[which(clones[,1]==i),2]))
# diversity[i,2] <- diversity[i,1] / (data_avg$N[i] + data_avg$M[i])
diversity[i,2] <- length(table(clones[which(clones[,1]==i),4]))
}  

# postscript(file =  "Figures/N_clones.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 10, width = 10)

g_range_y <- range(0, max(diversity[,2],diversity[,1]))

plot(1:time_max,diversity[,1],type = "l",col = "blue", lwd=2, ylim=g_range_y, xlab = "Time step", ylab = "Number of clones") # number of clones - diversity
lines(1:time_max,diversity[,2],type = "l",col = "red", lwd=2)

legend(1, 1.15*g_range_y[2], c("Drivers","Drivers+Passengers"), 
       lwd=1,cex=1.2,col=c("blue","red"), lty = 1:1, horiz = TRUE, bty = "n")


save_fig("Figures/N_clones.eps")

rl <-  readline(prompt="This is a plot for Number of clones - Press Enter  ")


#dev.copy(pdf, "Figures/N_clones.pdf")    
# dev.copy2eps(file = "Figures/N_clones.eps", height = 10, width = 10) 

# plot(1:time_max,diversity[,2],type = "l",col = "red", lwd=2,xlab = "Time step", ylab = "Normilized diversity") # normalized diversity

# rl <-  readline(prompt="This is a plot for Numbers of cells in each clone - Press Enter  ")


# Number of cells in each clone

x <- table(clones[,2]) # names of all clones - "numeric" name of clones

total_clones <- length(names(x)) # number of all clones

 evolution_clones <- data.frame(matrix(nrow = time_max, ncol = total_clones))

 names(evolution_clones) <- names(x)
 
for (i in 1:time_max) {
y <-  table(clones[which(clones[,1]==i),2])
evolution_clones[i,match(names(y),names(evolution_clones))] <- y  
}

evolution_clones[is.na(evolution_clones)] <- 0


 tree_cl <- calc_tree(evolution_clones, clones, total_clones)

 tree_cl$label <- names(table(clone_tree$Parent_ID))
 
 # ggtree(tree_cl) + theme_tree2()
 
 # label_nodes = names(table(clone_tree$Parent_ID))
 
 #p <- ggtree(tree_cl, color="blue", size=1.5, linetype=1)     + 
#   geom_nodepoint(mapping = NULL, data = NULL, position = "identity", na.rm = FALSE, show.legend = FALSE, color = "red", fill = "red", size = 3, shape=23)   +   
 #        geom_tippoint(color = "skyblue", size = 3) + 
  #         geom_tiplab(size = 6) + 
   #          geom_rootedge(color="blue", size=1.5, linetype=1) 

 
 
 p <- ggtree(tree_cl, color="blue", size=1.5, linetype=1, ladderize = TRUE)     + 
   geom_nodepoint(mapping = NULL, data = NULL, position = "identity", na.rm = FALSE, show.legend = FALSE, color = "red", fill = "red", size = 3, shape=23)   +   
   geom_tippoint(color = "skyblue", size = 3)  
#    geom_rootedge(color="blue", size=1.5, linetype=1) 
 
 d <- p$data
 d$label[is.na(d$label)] <- tree_cl$label
 
 
 plot( p + geom_text2(data=d, aes(label=label), nudge_x = 0.04, size = 6 )  )
 

 save_fig("Figures/ggtree_clones.eps")
 
 rl <-  readline(prompt="This is a ggtree plot for  clones - Press Enter  ") 
  ######## plot TREE:
  # plot(tree_cl, cex = 1.4)
  
 
  ### APE library:
  plotTreeTime(tree_cl, tip.dates = clone_tree$Time_start, show.tip.label = TRUE, label.offset = 0.01 )
  
  nodelabels(text = tree_cl$label, frame = "circle")
  # tiplabels()
  # dgelabels()
  
  save_fig("Figures/Tree_clones.eps")
 
  rl <-  readline(prompt="This is a Time tree plot for  clones - Press Enter  ") 
  
  #dev.copy(pdf, "Figures/Tree_clones.pdf")    
  #dev.off()
  
  
  # to plot all dependeces 
  
# Make a large number of colors

nm <- length(evolution_clones)

w <- (nm^(1/3)) %/% 1 +1

st <-  w^3 %/% nm

sq <- seq(0,1-1/w,1/w)

cr <- 1:nm

l <- 0
R <- 1:(w^3)
G <- R
B <- R

for (i in 1:w) {
  for (j in 1:w) {
    for (k in 1:w) {
      l <- l+1
      R[l] <- sq[i]
      G[l] <- sq[j]
      B[l] <- sq[k]
      
    } 
  }  
}

# seq(1,w^3,st) # is consequence of each color to make a high diversity of colors
jColor <- data.frame(number = 1:length(seq(1,w^3,st)),color = rgb(R[seq(1,w^3,st)],G[seq(1,w^3,st)],B[seq(1,w^3,st)]))

print("This is a plot for  number of cells in each clone: ") 

for (ll in 1:2) {
  
  if (ll == 1) g_range_y <- range(0, max(evolution_clones)) else g_range_y <- range(0,100)

  plot(1:time_max,evolution_clones[,1],type = "l",col = jColor$color[1],lwd=3,xlab = "Time step",
     ylab = "Number of cells in each clone",ylim=g_range_y)

  for (i in 2:total_clones) {
    lines(1:time_max,evolution_clones[,i],col = jColor$color[i],lwd=3)
    }
  if (ll == 1) rl <-  readline("Press Enter for another scale")
 
  #if (ll == 1) dev.copy(pdf, "Figures/N_cells_in_clones_1.pdf")   else dev.copy(pdf, "Figures/N_cells_in_clones_2.pdf")    
  #dev.off()
  if (ll == 1)   save_fig("Figures/N_cells_in_clones_1.eps")    else save_fig("Figures/N_cells_in_clones_2.eps") 
  }

rl <-  readline("Press Enter")

# the clones at the last time step
cl <- clones[which(clones[,1]==time_max),]

# Most popular clone
barplot(table(cl[,2]),xlab = "The ID of clone", ylab = "Number of cells in the clone", main = "FOR DRIVERS", 
         cex.name = 2, space=0.7, col = "green")

save_fig("Figures/Barplot_N_cells_in_clones.eps") 
#dev.copy(pdf, "Figures/Barplot_N_cells_in_clones.pdf")    
#dev.off()
rl <-  readline(prompt="This is a barplot for Numbers of cells in clones at last timestep - Press Enter ")


# Most popular clone (Drivers and Passengers)
barplot(table(cl[,4]),xlab = "The ID of clone", ylab = "Number of cells in the clone", main = "FOR DRIVERS AND PASSENGERS",
        cex.name = 1.6, space=0.7, col = "green")
save_fig("Figures/Barplot_N_cells_in_clones_DP.eps")
#dev.copy(pdf, "Figures/Barplot_N_cells_in_clones_DP.pdf")    
#dev.off()
rl <-  readline(prompt="This is same barplot for DRIVERS AND PASSENGERS - Press Enter ")


ineq_clones <- matrix(0,nrow = time_max,ncol = 2)

# Inequality measure 
for (k in 1:time_max) {
  cl <- clones[which(clones[,1]==k),2]
ineq_clones[k,1] <- ineq(cl,type = "Gini")

cl2 <- clones[which(clones[,1]==k),4]
ineq_clones[k,2] <- ineq(cl2,type = "Gini")
}

g_range_y <- range(0, 1) # max(ineq_clones[,2], ineq_clones[,1]))

plot(1:time_max,ineq_clones[,1],type = "l", lwd=3, pch = 19, col = "blue", 
      ylim = g_range_y, xlab = "Time step", ylab = "Inequlity coefficient")
lines(1:time_max,ineq_clones[,2],type = "l", lwd=3, pch = 19, col =  "red")


# g_range_x[1]/2+g_range_x[2]/2.5, 1.2*g_range_y[2]
legend(1, 1.15*g_range_y[2], c("Drivers","Drivers+Passengers"), 
       lwd=2,cex=1.1,col=c("blue","red"), lty = 1:1, horiz = TRUE, bty = "n")

save_fig("Figures/Inequality.eps")
#dev.copy(pdf, "Figures/Inequality.pdf")    
#dev.off()
rl <-  readline(prompt="This is a plot for inequality coefficient - Press Enter ")



plot(data_avg$N,ineq_clones[,1],type = "l", lwd=3, pch = 19, col = "blue", ylim = g_range_y, xlab = "Number of primary tumor cells", ylab = "Inequlity coefficient")
points(data_avg$N,ineq_clones[,2],type = "l", lwd=3, pch = 19, col =  "red")
legend(min(data_avg$N)-1, 1.15*g_range_y[2], c("Drivers","Drivers+Passengers"), 
       lwd=2,cex=1.1,col=c("blue","red"), lty = 1:1, horiz = TRUE, bty = "n")

save_fig("Figures/Inequality_primary.eps")
#dev.copy(pdf, "Figures/Inequality_primary.pdf")    
#dev.off()

rl <-  readline(prompt="This is a plot for inequality coefficient for primary tumor cells - Press Enter ")



plot(data_avg$M,ineq_clones[,1],type = "l", lwd=3, pch = 19, col = "blue", ylim = g_range_y, xlab = "Number of metastasis cells", ylab = "Inequlity coefficient")
points(data_avg$M,ineq_clones[,2],type = "l", lwd=3, pch = 19, col =  "red")
legend(min(data_avg$M), 1.15*g_range_y[2], c("Drivers","Drivers+Passengers"), 
       lwd=2,cex=1.1,col=c("blue","red"), lty = 1:1, horiz = TRUE, bty = "n")

save_fig("Figures/Inequality_metastasis.eps")
#dev.copy(pdf, "Figures/Inequality_metastasis.pdf")    
#dev.off()
rl <-  readline(prompt="This is a plot for inequality coefficient for metastasis cells - Press Enter ")



plot(data_avg$M+data_avg$N,ineq_clones[,1],type = "l", lwd=3, pch = 19, col = "blue", ylim = g_range_y, xlab = "Number of all cells", ylab = "Inequlity coefficient")
points(data_avg$M+data_avg$N,ineq_clones[,2],type = "l", lwd=3, pch = 19, col =  "red")
legend(min(data_avg$M+data_avg$N), 1.15*g_range_y[2], c("Drivers","Drivers+Passengers"), 
       lwd=2,cex=1.1,col=c("blue","red"), lty = 1:1, horiz = TRUE, bty = "n")
save_fig("Figures/Inequality_all_cells.eps")
#dev.copy(pdf, "Figures/Inequality_all_cells.pdf")    
#dev.off()

rl <-  readline(prompt="This is a plot for inequality coefficient for all cells - Press Enter ")

##### VAF 

### data_last has information for all mutations in genes of all cells 
### from 22 column         to 21 + length(onco) - for drivers and 
### from 22 + length(onco) to 21 + 2*length(onco) for passangers


VAF <- NULL
N <- data_last$N[1]
M <- data_last$M[1]
N_all <- N + M

for (k in 22:(21 + 2 * onco$len)) {
  
  if (k > 21+onco$len) DriverPasngr <- "P" else DriverPasngr <- "D" 
  if (k > 21+onco$len) Gene <- onco$name[ k - 21 - onco$len] else Gene <- onco$name[ k - 21]
  d <- as.vector( data_last[,k] )
  d <- as.character(d)
  d <- str_replace_all(d,":.....,",",")
  d <- str_replace_all(d,":....,",",")
  d <- str_replace_all(d,":...,",",")
  d <- str_replace_all(d,":..,",",")
  d <- str_replace_all(d,":.,",",")

  d <- str_replace_all(d,":.....","")
  d <- str_replace_all(d,":....","")
  d <- str_replace_all(d,":...","")
  d <- str_replace_all(d,":..","")
  d <- str_replace_all(d,":.","")

  d_Primary    <- d[which(data_last$im  < 1)]
  d_Metastatic <- d[which(data_last$im == 1)]
  
  d <- str_split(d,",")  
  d <- unlist(d)
  d <- as.integer(d)

  d_Primary <- str_split(d_Primary,",")  
  d_Primary <- unlist(d_Primary)
  d_Primary <- as.integer(d_Primary)

  d_Metastatic <- str_split(d_Metastatic,",")  
  d_Metastatic <- unlist(d_Metastatic)
  d_Metastatic <- as.integer(d_Metastatic)  
  
  Z_M <- TRUE
  Z_N <- TRUE
  Z_M <- all(is.na(d_Metastatic))
  Z_N <- all(is.na(d_Primary))
  
  out <- NULL
  out <- as.data.frame( table(d) )
  if (!Z_M | !Z_N) names(out) <- c("pos","Freq")
  
  if (!Z_N && N > 0) {
      out_Prim <- NULL
      out_Prim <- as.data.frame( table(d_Primary) )
      names(out_Prim) <- c("pos","Freq_Prim")
      out <- merge.data.frame(out_Prim, out, by = "pos" , all = TRUE)
      
      } else out["Freq_Prim"] <- 0
  
  
  if (!Z_M && M > 0) {
      out_Met <- NULL
      out_Met <- as.data.frame( table(d_Metastatic) )
      names(out_Met) <- c("pos","Freq_Met")
      out <- merge.data.frame(out_Met, out, by = "pos" , all = TRUE)
      
      } else out["Freq_Met"] <- 0
  
  
      out[is.na(out)] <- 0 
  
      if ( N > 0 ) out$VAF_Prim <- 0.5 * out$Freq_Prim / N   else out$VAF_Prim <- 0
      if ( M > 0 ) out$VAF_Met  <- 0.5 * out$Freq_Met  / M   else out$VAF_Met  <- 0
      if (!Z_M | !Z_N) out$VAF <- 0.5 * out$Freq / N_all     else out$VAF      <- 0
      
            
  # nm <- names(data_last[k])
  VAF1 <- NULL
  VAF1 <- cbind.data.frame(DriverPasngr, Gene, out$pos, 
                out$VAF_Prim,  out$Freq_Prim, N, 
                out$VAF_Met ,  out$Freq_Met,  M,
                out$VAF,  out$Freq, N_all)
  VAF <- rbind.data.frame(VAF,VAF1)
  
  }

header <- c( "DriverPasngr", "Gene", "Position", 
             "VAF_Primary", "Ncells_Primary_wMutation", "Ncells_Primary",
             "VAF_Metastatic", "Ncells_Metastatic_wMutation", "Ncells_Metastatic", 
             "VAF_PriMet", "Ncells_PriMet_wMutation", "Ncells_PriMet" )

names(VAF) <- header
write.table(VAF,file = "Output/VAF.txt", append = FALSE, row.names = FALSE, sep="\t")

print("VAF is saved to the file `Output/VAF.txt` ")
