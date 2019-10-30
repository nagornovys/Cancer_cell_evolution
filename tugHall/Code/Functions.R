# The code to analize the results of simulation

save_fig <- function(file_input) {
  
  file_input <- paste0(substr(file_input,1,nchar(file_input)-4),".jpg")
  #dev.copy2eps(file = file_input, height = 10, width = 10) 
  dev.copy(jpeg,file_input, width = 10, height = 10, units = 'in', res = 300)
  dev.off()
  par(xpd=TRUE, cex.lab=2, lwd = 2, mar = c(5, 5, 5, 5), tcl = 0.5, cex.axis = 1.75,  mgp = c(3, 0.6, 0))
  
}


analyze_data <- function(celloutfile) {
  # Analising of results:
  data_out <- read.csv(celloutfile, sep="\t")
  data_out[is.na(data_out)] <- ""
  # make a readible names
  names(data_out)[5] <-  "c"
  names(data_out)[6] <-  "d"
  names(data_out)[7] <-  "i"
  names(data_out)[8] <-  "im"
  names(data_out)[9] <-  "a"
  names(data_out)[10] <- "k"
  names(data_out)[11] <- "E"
  names(data_out)[13] <- "Nmax"
  names(data_out)

  # average data
  data_avg <<- data_out[which(data_out$AvgOrIndx == "avg"),]

  # data without averaging - flow data
  data_flow <<- data_out[which(!data_out$AvgOrIndx == "avg"),]
  
  # the data of the last time step 
  time_max <<- max(data_flow$Time)
  data_last <<- data_flow[which(data_flow$Time == time_max),]
        
          # let draw graphics 
          # Numbers of Metastasis and primary tumor cells
          g_range_y <- range(0, data_avg$N,data_flow$M)
          g_range_x <- range(min(data_avg$Time),max(data_flow$Time))
          plot(data_avg$Time,data_avg$M,type = "l",xlab = "Time step",
           ylab = "Number of cells",ylim=g_range_y,xlim = g_range_x,col = "red")
          lines(data_avg$Time,data_avg$N,type = "l",col = "blue")
          # g_range_x[1]/2+g_range_x[2]/2.5, 1.2*g_range_y[2]
          legend(g_range_x[1], 1.15*g_range_y[2], c("Primary tumor","Metastasis"), 
             lwd=2,cex=1.3,col=c("blue","red"), lty = 1:1,horiz = TRUE, bty = "n")
  #dev.copy(pdf, "Figures/N_cells.pdf")    
  #dev.off()

          save_fig("Figures/N_cells.eps")
  
  rl <-  readline(prompt="This is a plot for Numbers of primary tumor and Metastasis cells - Press Enter  ")


  # Average values of probabilities
  g_range_y <- range(0, data_avg[6:10])
  g_range_x <- range(min(data_avg$Time),max(data_flow$Time))
  
      
        plot(data_avg$Time,data_avg$d,type = "l", ylim=g_range_y,xlim = g_range_x,
         xlab = "Time step",ylab = "Average probabilities",col = "red")
    
        lines(data_avg$Time,data_avg$i,type = "l",col = "blue")
        lines(data_avg$Time,data_avg$im,type = "l",col = "green")
        lines(data_avg$Time,data_avg$a,type = "l",col = "orange")
        lines(data_avg$Time,data_avg$k,type = "l",col = "black")

        legend(g_range_x[1], 1.15*g_range_y[2], c("d","i","im","a","k"), cex=1.3,
               col=c("red","blue","green","orange","black"), lty = 1:1,lwd = 2,horiz = TRUE, bty = "n")
      

        #dev.copy(pdf, "Figures/Probabilities.pdf")    
        #dev.off()
        
        save_fig("Figures/Probabilities.eps")
        
  rl <-  readline(prompt="This is a plot for Average values of probabilities - Press Enter  ")

    # The averaged values of Hallmarks 
    g_range_y <- range(0, data_avg[15:19])
    g_range_x <- range(min(data_avg$Time),max(data_flow$Time))
  
    plot(data_avg$Time,data_avg$Hd,type = "l", ylim=g_range_y,xlim = g_range_x,
       xlab = "Time step",ylab = "Averaged Hallmarks values",col = "red")
  
    lines(data_avg$Time,data_avg$Hi,type = "l",col = "blue",lwd = 2)
    lines(data_avg$Time,data_avg$Him,type = "l",col = "green",lwd = 2)
    lines(data_avg$Time,data_avg$Ha,type = "l",col = "orange",lwd = 2)
    lines(data_avg$Time,data_avg$Hb,type = "l",col = "black",lwd = 2)

    legend(g_range_x[1], 1.15*g_range_y[2], c("Hd","Hi","Him","Ha","Hb"), cex=1.2,
           col=c("red","blue","green","orange","black"), lty = 1:1,lwd = 2,horiz = TRUE, bty = "n")


    #dev.copy(pdf, "Figures/Hallmarks.pdf")    
    #dev.off()
    save_fig("Figures/Hallmarks.eps")  
  
  rl <-  readline(prompt="This is a plot for Average values of Hallmarks - Press Enter, but next takes a time for calculations")

  onco <<- oncogene$new()        # make the vector onco about the hallmarks
  onco$read(genefile)          # read the input info to the onco from genefile - 'gene_cds2.txt'
  
  order_dysfunction <<- data.frame(data_last[,21:(21+length(onco$name))],stringsAsFactors = FALSE)
  names(order_dysfunction) <<- c("ID",onco$name)
  order_dysfunction[,1] <<- data_last$ID

  for (i in 1:(length(onco$name)+1)) {
    order_dysfunction[,i] <<- as.character(order_dysfunction[,i])
  }
  #str(order_dysfunction)
  order_dysfunction$ID <<- as.character(order_dysfunction$ID)
  #str(order_dysfunction)



  # substr(x, regexpr(":",x)+1, ifelse(isTRUE(grep(",",x)>0),(regexpr(",",x)-1),nchar(x)))


  for (k in 1:length(order_dysfunction[,1])) {

      for (i in 2:(length(onco$name)+1)) {

          x<- order_dysfunction[k,i]
         order_dysfunction[k,i] <<- substr(x, regexpr(":",x)+1, ifelse(isTRUE(grep(",",x)>0),(regexpr(",",x)-1),nchar(x)))
      }
  }


  for (i in 1:length(onco$name)+1) {
      order_dysfunction[,i] <<- as.integer(order_dysfunction[,i])
  }

  outfile <- 'Output/Order_of_dysfunction.txt'
  header <- c('Order of gene dysfunction: from first to last', 'Frequency or number of cells with same order')
  write(header, outfile, append=FALSE, ncolumn=length(header), sep="\t")  
  
  #for (i in 1:length(order_dysfunction[,1])) {
  #  data <- c(order_dysfunction$ID[i],names(sort(order_dysfunction[i,2:(length(onco$name)+1)])))
  #  write(data, outfile, append=TRUE, ncolumn=length(data), sep="\t")
  #    }
  
  #    write_order_of_dysfunction('Order_of_dysfunction.txt', env, cells, isFirst)


  x <- array("",dim = length(order_dysfunction[,1]))
  
  for (i in 1:length(order_dysfunction[,1])) {
    data <- c(names(sort(order_dysfunction[i,2:(length(onco$name)+1)])))
    x[i] <- paste(data,collapse = " -> ")
  }  
  
  print("The order of gene dysfunction for each cell is saved into the file")
  print("This is function to save the order of gene dysfunction to the file `Order_of_dysfunction.txt` ")
  # The order of gene dysfunction 
  
  # find the unique orders of genes dysfunction
  uniq_order <<- table(x)
  uniq_order <<- sort(uniq_order,decreasing = TRUE)
  print("Unique of order of genes dysfunction and it's frequency:")
  for (i in 1:length(uniq_order)) {
    print(c(names(uniq_order)[i],uniq_order[[i]]))
    data <- c(names(uniq_order)[i],uniq_order[[i]])
    write(data, outfile, append=TRUE, ncolumn=length(data), sep="\t")
  }
}  # End of function





# Functions to calculate the inequality coefficients:

ineq <- function(x, parameter = NULL, type=c("Gini", "RS", "Atkinson", "Theil",
                                             "Kolm", "var", "square.var", "entropy"), na.rm = TRUE)
{
  switch(match.arg(type),
         Gini = Gini(x, na.rm = na.rm),
         RS = RS(x, na.rm = na.rm),
         Atkinson = Atkinson(x, parameter = parameter, na.rm = na.rm),
         Theil = Theil(x, parameter = parameter, na.rm = na.rm),
         Kolm = Kolm(x, parameter = parameter, na.rm = na.rm),
         var = var.coeff(x, na.rm = na.rm),
         square.var = var.coeff(x, square=TRUE, na.rm = na.rm),
         entropy = entropy(x, parameter = parameter, na.rm = na.rm))
}

Gini <- function(x, corr = FALSE, na.rm = TRUE)
{
  if(!na.rm && any(is.na(x))) return(NA_real_)
  x <- as.numeric(na.omit(x))
  n <- length(x)
  x <- sort(x)
  G <- sum(x * 1L:n)
  G <- 2 * G/sum(x) - (n + 1L)
  if (corr) G/(n - 1L) else G/n
}

RS <- function(x, na.rm = TRUE)
{
  if(!na.rm && any(is.na(x))) return(NA_real_)
  x <- as.numeric(na.omit(x))
  d <- abs(x - mean(x))
  d <- mean(d)/(2*mean(x))
  d
}

Atkinson <- function(x, parameter = 0.5, na.rm = TRUE)
{
  if(!na.rm && any(is.na(x))) return(NA_real_)
  x <- as.numeric(na.omit(x))
  if(is.null(parameter)) parameter <- 0.5
  if(parameter==1)
    A <- 1 - (exp(mean(log(x)))/mean(x))
  else
  {
    x <- (x/mean(x))^(1-parameter)
    A <- 1 - mean(x)^(1/(1-parameter))
  }
  A
}

var.coeff <- function(x, square=FALSE, na.rm = TRUE)
{
  if(!na.rm && any(is.na(x))) return(NA_real_)
  x <- as.numeric(na.omit(x))
  n <- length(x)
  V <- sqrt((n-1)*var(x)/n)/mean(x)
  if(square) V <- V^2
  V
}

Theil <- function(x, parameter = 0, na.rm = TRUE)
{
  if(!na.rm && any(is.na(x))) return(NA_real_)
  x <- as.numeric(na.omit(x))
  if(is.null(parameter)) parameter <- 0
  if(parameter==0)
  {
    x <- x[!(x==0)]
    Th <- x/mean(x)
    Th <- sum(x*log(Th))
    Th <- Th/sum(x)
  }
  else
  {
    Th <- exp(mean(log(x)))/mean(x)
    Th <- -log(Th)
  }
  Th
}

Kolm <- function(x, parameter = 1, na.rm = TRUE)
{
  if(!na.rm && any(is.na(x))) return(NA_real_)
  x <- as.numeric(na.omit(x))
  x <- as.numeric(x)
  if(is.null(parameter)) parameter <- 1
  KM <- parameter * (mean(x)-x)
  KM <- mean(exp(KM))
  KM <- (1/parameter)*log(KM)
  KM
}

entropy <- function(x, parameter = 0.5, na.rm = TRUE)
{
  if(!na.rm && any(is.na(x))) return(NA_real_)
  x <- as.numeric(na.omit(x))
  x <- as.numeric(x)
  if(is.null(parameter)) parameter <- 0.5
  if(parameter==0)
    e <- Theil(x, parameter = 1)
  else
    if(parameter==1)
      e <- Theil(x, parameter = 0)
  else
  {
    k <- parameter
    e <- (x/mean(x))^k
    e <- mean(e-1)/(k*(k-1))
  }
  e
}

    



######### Function to calulate TREE of CLONES:

#### To plot tree of clones
calc_tree <- function(evolution_clones, clones, total_clones) { 
      clone_tree <<- data.frame(matrix(nrow = total_clones, ncol = 5))
      clone_tree[,] <<- 0
      names(clone_tree) <<- c("Clone_ID","Time_start","Time_end","Parent_ID","Length")
      clone_tree[,1] <<- names(evolution_clones)      # Clone_ID
      
      for ( i in 1:total_clones ) {
        st <- min(which(evolution_clones[,i] > 0) )
        en <- max(which(evolution_clones[,i] > 0) )
        clone_tree[i,2] <<- st
        clone_tree[i,3] <<- en    
        

        # To find parent ID: find first cell of clone -> find parent ID of it -> 
        ##################################        ->   find clone ID of parent cell at timestep  = birthday
        ### find first cell of clone:
        ID <- as.numeric(clone_tree[i,1])
        jk <- min(which(clones[,2] == ID))
        P_ID <- clones[jk,"Parent_ID"]   # cell ID of parent
        BirthDay <- clones[jk,"Birthday"]  # birthday time step
        
        ## find clone ID of parent cell at timestep  = birthday
        k <- min(which( (clones[,"Time"] == BirthDay) & (clones[,"Cell_ID"] == P_ID)))  # cell ID of parent at timestep
        
        clone_tree[i,"Parent_ID"] <<- clones[k,"Clone_ID"]
      }
      
      clone_tree[is.na(clone_tree)] <<- 0 
      
      clone_tree$Clone_ID <<- as.numeric(clone_tree$Clone_ID)
      clone_tree$Length <<- clone_tree$Time_end - clone_tree$Time_start
      
      clone_tree$Length <<- clone_tree$Length / max(clone_tree$Length)
      
      clone_tree$Length[ which(clone_tree$Length == 0) ] <<- 0.01
      
      ### order(clone_tree$Time_start) gives the order of appearence of clones:
      clone_tree[,] <<- clone_tree[order(clone_tree$Time_start),]
      
      library(ape)
      n <- length(table(clone_tree$Parent_ID))
      
      parents <- unlist(dimnames(table(clone_tree$Parent_ID)))
      parents <- as.numeric(parents)
      
      children <- as.numeric(clone_tree$Clone_ID)
      
      m <- length(children)
      
      cl_edge <- data.frame(p=1:m,ch=1:m)
      cl_edge[,] <- -1
      
      k <- 1
      l_edge <- NULL
      
      for (i in 1:n) {
        vec_p <-  which(clone_tree$Parent_ID == parents[i])
        for (j in vec_p) {
          
          if ( children[j] %in% parents ) {
            if ( children[j] != parents[1] )  {
              cl_edge[k,] <- c(m+i, ( m + which(children[j] == parents)) )  
              l_edge[k] <- 1.0
              k <- k + 1
            }
            cl_edge[k,] <- c( (m + which(children[j] == parents)) , j )
            l_edge[k] <- 0.5
          }  else  { 
            cl_edge[k,] <- c(m+i,j) 
            l_edge[k] <- 0.5
          }
          k <- k + 1
        } 
      }
      
      tree_cl <- rtree(n = 2, rooted = TRUE)
      
      tree_cl$Nnode <- n 
      tree_cl$tip.label <- as.character( clone_tree$Clone_ID )
      tree_cl$edge <- as.matrix(cl_edge)
      tree_cl$edge.length <- l_edge
      attr(tree_cl$edge,"dimnames") <- NULL
      
      

      write.tree(tree_cl,file = "null")
      
      tree_cl <- read.tree(file = "null")
      tree_cl$root.edge <- 0.15
      file.remove("null")
      
return(tree_cl)
}
    


