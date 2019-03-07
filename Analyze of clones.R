  # The files to analize the average data for each time step and trial, and hole data of last time step

   celloutfile <- "cancer/cellout.txt"
# celloutfile <- "single/cellout.txt"   # average data for each time step
#   celloutfile <- "multi/cellout.txt"

   cell_last <- "cancer/last_steps.txt"
#  cell_last <- "single/last_steps.txt"  # hole data of last time step for each trial
#   cell_last <- "multi/last_steps.txt"
  ############ To read the onco$name data from file #####################
  
  onco = oncogene$new()        # make the vector onco about the hallmarks
  onco$read("cancer/gene_cds2_third.txt")          # read the input info to the onco from genefile - 'gene_cds2.txt'
#  onco$read("multi/gene_cds2_multi.txt") 
#   onco$read("single/gene_cds2.txt") 
     onco$name
  
  
  # The function is in the Analize.R file 
  analize_data(celloutfile,cell_last)

  # time_max - number of time steps 
  # data_avg - data of average values rom analize_data() function
  # data_last - data of last time steps for each trial from analize_data() function

  # The clones calculation
  # table(data_last$Clone.number)

  clones <- matrix(nrow = length(data_last$Time), ncol = 5)
  clones[,1] <- as.integer(as.character(data_last$Trial))
  clones[,2] <- as.integer(as.character(data_last$Time))
  clones[,3] <- as.integer(as.character(data_last$Clone.number))
  clones[,4] <- as.integer(as.character(data_last$Passengers.Clone.number))
  clones[,5] <- as.integer(as.character(data_last$Mix.Clone.number))


  # Distribution of clone's ID (Drivers only)
  barplot(sort(table(clones[,3]),decreasing = TRUE),
        main = "Data for drivers only", 
        xlab = "The clone's ID", ylab = "Number of cells in a clone", 
        space=0.7, col = "green", log = "y")


  # Distribution of clone's ID (Drivers and Passengers)
  barplot(sort(table(clones[,5]),decreasing = TRUE),
        main = "Data for drivers and passengers", 
        xlab = "The clone's ID", ylab = "Number of cells in a clone", 
        space=0.7, col = "green", xaxt='n', log = "y")


  diversity <- matrix(nrow = 100, ncol = 2)

  for (i in 1:100) {
    diversity[i,1] <- length(table(clones[which(clones[,1]==i),3]))
    # diversity[i,2] <- diversity[i,1] / (data_avg$N[i] + data_avg$M[i])
    diversity[i,2] <- length(table(clones[which(clones[,1]==i),5]))
  }  

  g_range_y <- range(0, max(diversity[,2],diversity[,1]))
  plot(1:100,diversity[,1],type = "l",col = "blue", lwd=2, ylim=g_range_y, xlab = "ID of Trial", 
     ylab = "Number of clones") # number of clones - diversity

  lines(1:100,diversity[,2],type = "l",col = "red", lwd=2)

  legend(1, 1.2*g_range_y[2], c("Drivers","Drivers+Passengers"), 
       lwd=2,cex=1,col=c("blue","red"), lty = 1:1,horiz = TRUE)



  # The scatter plot for Drivers vs Drivers+Passengers 
  g_range_x <- range(min(diversity[,1])-1, max(diversity[,1])+1)
  g_range_y <- range(min(diversity[,2])-1, max(diversity[,2])+1)
  plot(diversity[,1],diversity[,2],type = "p",col = "blue", main = "Number of clones", 
     xlim = g_range_x, ylim = g_range_y, cex = 1.2, pch = 18, 
     xlab = "Drivers only", ylab = "Drivers and Passengers") 

  lines(smooth.spline(x = diversity[,1], y = diversity[,2]), col = "red", lwd = 3, lt = 2)
  

  # Distribution of clone number (Drivers)
  # barplot(table(diversity[,1])/max(table(diversity[,1])),xlab = "Number of clones (Diversity)", ylab = "Freqiency", space=0.7, col = "green")

  for (i in 1:2) {
  
      x<- diversity[,i]
      fit <- fitdistr(x, "normal")
      class(fit)
      # [1] "fitdistr"
      para <- fit$estimate
      #         mean            sd 
      #-0.0002000485  0.9886248515 

      hist(x, breaks = max(x), prob = TRUE, border = TRUE, 
        xlab = "Number of clones", ylab = "Freqiency", col = "green", 
        main = "Histogram and normal distribution of Diversity")

      # barplot(table(x)/max(table(x)),xlab = "Number of clones (Diversity)", ylab = "Freqiency", space=0.7, col = "green")
      curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE, lwd = 3 )


      rl <-  readline(prompt="This is a barplot for Driver (first) and for Drivers (second) only - Press Enter")
  }


  # Distribution of clone number (Drivers and Passengers)
  #barplot(table(diversity[,2])/max(table(diversity[,2])),xlab = "Number of clones (Diversity)", ylab = "Freqiency", space=0.7, col = "green")


  ineq_clones <- matrix(0,nrow = 100,ncol = 2)

  # Inequality measure 
  for (k in 1:100) {
    cl <- clones[which(clones[,1]==k),3]
    ineq_clones[k,1] <- ineq(cl,type = "Gini")

    cl2 <- clones[which(clones[,1]==k),5]
    ineq_clones[k,2] <- ineq(cl2,type = "Gini")
  }

  ineq_clones[is.na(ineq_clones)] <- 0


  layout(matrix(c(1,2), nrow=1, ncol=2),widths=c(1,1))
  layout.show(2)

  hist(ineq_clones[,1],  freq  = TRUE, border = TRUE,
      xlab = "Inequality coefficient", ylab = "Freqiency, %", col = "green", 
      main = "Drivers only")

  hist(ineq_clones[,2],  freq = TRUE, border = TRUE,
      xlab = "Inequality coefficient", ylab = "Freqiency, %", col = "green", 
      main = "Drivers and passengers")
 
  layout(1)  
 

  # The order of gene dysfunction 
  
  order_dysfunction <<- data.frame(data_last[,23:(22+length(onco$name))],stringsAsFactors = FALSE)
  names(order_dysfunction) <- onco$name
  #  order_dysfunction[,1] <- 0
  order_dysfunction[,length(order_dysfunction[1,])+1] <- data_last[,1]
  
  
  # The VAF - Varient allel frequency functions for each gene 
  
  num_genes <- length(onco$name)   # Number ofgenes
  
  
  len_data <- length(order_dysfunction[,1])
  len_data
  
  str(order_dysfunction)
  for (i in 1:length(onco$name)) {
    order_dysfunction[,i] <- as.character(order_dysfunction[,i])
  }
  str(order_dysfunction)
  
  
  
  # Please, choose parameters with accordence with your view!!!!
  # layout(matrix(c(1,2,3,4,5,6), nrow=2, ncol=3)) # ,widths=c(1,1,1))
  # layout.show(6)

  for (k in 1:num_genes) {
    
    VAF_data <- order_dysfunction[,c(k,num_genes+1)]
    # to delete the time step data for VAF analysis 
    VAF_data[1:len_data,1] <- sub(":......,",",",VAF_data[1:len_data,1])
    VAF_data[1:len_data,1] <- sub(":.....,",",",VAF_data[1:len_data,1])
    VAF_data[1:len_data,1] <- sub(":....,",",",VAF_data[1:len_data,1])
    VAF_data[1:len_data,1] <- sub(":...,",",",VAF_data[1:len_data,1])
    VAF_data[1:len_data,1] <- sub(":..,",",",VAF_data[1:len_data,1])
    VAF_data[1:len_data,1] <- sub(":.,",",",VAF_data[1:len_data,1])
  
    # to delete the time step data in the end of value 
    VAF_data[1:len_data,1] <- sub(":......","",VAF_data[1:len_data,1])
    VAF_data[1:len_data,1] <- sub(":.....","",VAF_data[1:len_data,1])
    VAF_data[1:len_data,1] <- sub(":....","",VAF_data[1:len_data,1])
    VAF_data[1:len_data,1] <- sub(":...","",VAF_data[1:len_data,1])
    VAF_data[1:len_data,1] <- sub(":..","",VAF_data[1:len_data,1])
    VAF_data[1:len_data,1] <- sub(":.","",VAF_data[1:len_data,1])
  
     for (p in 1:100) {
      
      VAF <- sort(table(VAF_data[which(VAF_data[,2]==p),1]), decreasing = TRUE)
      VAF[which(names(VAF) == "")]  <- 0 # without mutation = 0 

      s_VAF <- sum(VAF)
      VAF <- VAF/s_VAF/2   # for 
      # barplot(VAF, # xlab = "The allel's ID", ylab = "Variant allele frequency", 
      #             main = paste("VAF function for gene ", onco$name[k]), space=0.7, col = "green",  xaxt='n', log = "y")
      # print(VAF)
      if ((k==1) & (p==1)) {VAF_total <- VAF} else {VAF_total <- c(VAF_total,VAF)}
  
    }
  }
  
  # layout(1)  

  hist(VAF_total, breaks = 20, border = TRUE, 
       xlab = "VAF", ylab = "Freqiency", col = "green", main = "Variant allele frequency")
  
  
  
    
  rm(data_avg)
  rm(data_last)
 
  for (i in 1:length(onco$name)) {
      order_dysfunction[1:len_data,i] <- sub("......:","",order_dysfunction[1:len_data,i])
      order_dysfunction[1:len_data,i] <- sub(".....:","",order_dysfunction[1:len_data,i])
      order_dysfunction[1:len_data,i] <- sub("....:","",order_dysfunction[1:len_data,i])
  }
 
  
  for (i in 1:length(onco$name)) {
    order_dysfunction[1:len_data,i] <- sub("...:","",order_dysfunction[1:len_data,i])
    order_dysfunction[1:len_data,i] <- sub("..:","",order_dysfunction[1:len_data,i])
    order_dysfunction[1:len_data,i] <- sub(".:","",order_dysfunction[1:len_data,i])
  }
  
  
  for (i in 1:length(onco$name)) {
    order_dysfunction[1:len_data,i] <- sub(",.....","",order_dysfunction[1:len_data,i])
    order_dysfunction[1:len_data,i] <- sub(",....","",order_dysfunction[1:len_data,i])
    order_dysfunction[1:len_data,i] <- sub(",...","",order_dysfunction[1:len_data,i])
    order_dysfunction[1:len_data,i] <- sub(",..","",order_dysfunction[1:len_data,i])
    order_dysfunction[1:len_data,i] <- sub(",.","",order_dysfunction[1:len_data,i])
  }
  
   
  re <- length(onco$name) + 1
  re
  
   str(order_dysfunction)

   order_dysfunction[1:len_data,re] <- ""
   
   for (k in 1:length(onco$name)) {
      order_dysfunction[1:len_data,re] <- paste(order_dysfunction[1:len_data,re],order_dysfunction[1:len_data,k], sep = ":")
   }
    
   order_dysfunction[1:len_data,re] <- substring(order_dysfunction[1:len_data,re], 2)
     
   # order_dysfunction[1:len_data,re] <- paste(order_dysfunction[1:len_data,1], order_dysfunction[1:len_data,2], order_dysfunction[1:len_data,3],order_dysfunction[1:len_data,4], sep = ":")
   
   RESULT <- sort(table(order_dysfunction[,re]), decreasing = TRUE)  
   RESULT
   # rm(order_dysfunction)
   
   n_data <- names(RESULT)
   n_data
   # to check - one example:
   strsplit(paste(" ", n_data[2], " "), ":")
   
   
   nam_data <- vector(mode = "character")
   
   for (m in 1:length(n_data)) {
  
        x <- as.integer(unlist(strsplit(paste(" ", n_data[m], " "),":")))
        if ( (sum(is.na(x)) < re-1) & (length(x) == re-1) ) {
             names(x) <- onco$name
             nam_data[m] <- paste(names(sort(x)), collapse = "->" )
        }
   }
   nam_data 
   
   names(RESULT) <- nam_data
   
   nam_uniq <- names(table(nam_data)) # names of all possible orders of gene dysfunctions
   
   TAB_RES <- 1:length(nam_uniq)
   names(TAB_RES) <- nam_uniq 
   
   
   for (k in 1:length(nam_uniq)) {
  
        TAB_RES[k] <- sum(RESULT[which(names(RESULT) == nam_uniq[k])])
   }
   
   for (k in 1:length(nam_uniq)) {
     
        print(paste(names( sort(TAB_RES, decreasing = TRUE)[k]) , 
               sort(TAB_RES, decreasing = TRUE)[k], sep = " - "))
   }
   ##########################################################################
   
   