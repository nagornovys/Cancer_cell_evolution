# The code to analize the results of simulation

  # The function to draw the plot with distribution

  library(MASS)
  library(RColorBrewer)

  # The function to show scale in the image method [from https://gist.github.com/menugget/7689145/raw/dac746aa322ca4160a5fe66c70fec16ebe26faf9/image.scale.2.r]
  
  #This function creates a color scale for use with the image()
  #function. Input parameters should be consistent with those
  #used in the corresponding image plot. The "axis.pos" argument
  #defines the side of the axis. The "add.axis" argument defines
  #whether the axis is added (default: TRUE)or not (FALSE).
  
  image.scale <- function(z, zlim, col = heat.colors(12),
                          breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
    if(!missing(breaks)){
      if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
    }
    if(missing(breaks) & !missing(zlim)){
      breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
    }
    if(missing(breaks) & missing(zlim)){
      zlim <- range(z, na.rm=TRUE)
      zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
      zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
      breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
    }
    poly <- vector(mode="list", length(col))
    for(i in seq(poly)){
      poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
    }
    xaxt <- ifelse(horiz, "s", "n")
    yaxt <- ifelse(horiz, "n", "s")
    if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
    if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
    if(missing(xlim)) xlim=XLIM
    if(missing(ylim)) ylim=YLIM
    plot(1,1,t="n",ann = FALSE, ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
    for(i in seq(poly)){
      if(horiz){
        polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
      }
      if(!horiz){
        polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
      }
    }
  }
  
  
  # The function to plot distributions for each time step in the form of color plot
  distr_plot <- function(X,Y,X_str,Y_str,Ave_str) {
#    data_out <- read.csv(celloutfile, sep="\t")
#    X <-  data_out$Time
#    Y <-  data_out$M  
#    X_str <- "Generation number"
#    Y_str <- "Number of cells"
#    Ave_str <- "Average number of cells"
#     View(X)
#     View(Y)
  
  maxT <- max(X)
  Ymax <- max(Y) * 1.2
  
  data_ave <- matrix(data = 0, ncol = 2,nrow = maxT)
  
  for (k in 2:maxT ) {
    
    data_ave[k,1] = k
    n <- which(X == k)
    data_ave[k,2] = sum(Y[n])/length(Y[n])
  }
  
  #       View(data_ave) 
  df <- data.frame(x = X,y = Y)
  
  # colors
  rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
  r <- rf(40)
  # Adjust binning (interpolate - can be computationally intensive for large datasets)
  k <- kde2d(df$x, df$y, n=c(300,300))
  k$z <- k$z / max(k$z)
  
  breaks <- 1:(length(r)+1)
  breaks <- breaks / length(r)
  
#  layout(matrix(c(1,1,1,1,0,2,0,0), nrow=4, ncol=2),widths=c(11,1),  heights=c(1,1,1,1))   #   ncol=3), widths=c(4,4,1), heights=c(4,1)
   layout(matrix(c(1,2), nrow=1, ncol=2),widths=c(11,1))
  
    layout.show(2)
  
  par(mar=c(5,5,1,1)) 
  image.default(k, col=r, breaks = breaks, xlab = X_str, ylab = Y_str, cex.lab = 1.3, ylim = c(0,Ymax),xlim = c(1,maxT),legend=TRUE)
  
  # For each time step
  for (kk in 2:(maxT-1)) {
    
    numb <- which( xor(df$x==kk,df$x==(kk+1)) )
    if (bandwidth.nrd(df$y[numb]) > 0) {
      k <- kde2d(df$x[numb], df$y[numb], n=c(2,100))
      k$z <- k$z / max(k$z)
    
      breaks <- 1:(length(r)+1)
      breaks <- breaks / length(r)
      image(k, col=r, breaks = breaks, add = TRUE)
    }
  }
  
  lines(data_ave[2:maxT,1],data_ave[2:maxT,2],type = "l",lty = 3, lwd = 3, col = "#000000")
  legend("top", c(Ave_str,""), cex = 1.2, bty = "n", col = c("#000000",NA),lty = c(3,NA),pch = c(NA,NA),pt.cex = 0.15,lwd = c(3,NA), horiz = TRUE)
  
  par(mar=c(20,1.5,7,1))
  #  image.scale(volcano, col=pal.1(length(breaks)-1), breaks=breaks, horiz=FALSE)
  image.scale(volcano, col=r, ylim = c(0,1), breaks=breaks, horiz=FALSE)
  
  box()

  
  layout(1)   # to make a normal screen mode
    
}   # End of distr_plot function

# To check the work of functions:
  
# celloutfile <- "1000_Cells_Single/cellout.txt" 
#  data_out <- read.csv(celloutfile, sep="\t") 
  
#  distr_plot(X = data_out$Time, Y = data_out$Nmax. , X_str = "Generation number", Y_str = "Number of cells", Ave_str = "Average number of Hi" )
  



  analize_data <- function(celloutfile,cell_last) {
# Analising of results:
  data_avg <<- read.csv(celloutfile, sep="\t")
  data_avg[is.na(data_avg)] <<- ""
# make a readible names
  names(data_avg)[6] <<-  "c"
  names(data_avg)[7] <<-  "d"
  names(data_avg)[8] <<-  "i"
  names(data_avg)[9] <<-  "im"
  names(data_avg)[10] <<-  "a"
  names(data_avg)[11] <<- "k"
  names(data_avg)[12] <<- "E"
  names(data_avg)[14] <<- "Nmax"
  names(data_avg)

# the data of the last time step 
  time_max <<- max(data_avg$Time)
  data_last <<- read.csv(cell_last, sep="\t")

  names(data_last) <<- names(data_avg)


# let draw graphics 
# Numbers of Metastasis and normal cells
  g_range_y <- range(0, data_avg$N,data_avg$M)
  g_range_x <- range(min(data_avg$Time),max(data_avg$Time))
  plot(data_avg$Time,data_avg$M,type = "p", cex=0.3,
     xlab = "Generation number", ylab = "Number of cells",
     ylim=g_range_y, xlim = g_range_x,
     col = "red")
  points(data_avg$Time,data_avg$N,type = "p",col = "blue",cex = 0.2, pch = 20)
  par(xpd=TRUE)

# g_range_x[1]/2+g_range_x[2]/2.5, 1.2*g_range_y[2]
  legend(g_range_x[1], 1.2*g_range_y[2], c("Normal","Metastasis"), 
       cex=1,col=c("blue","red"), pch = 20:20, horiz = TRUE)

  rl <-  readline(prompt="This is a plot for Normal and Metastasis Numbers of cells - Press Enter  ")


# Average values of probabilities
  g_range_y <- range(0, data_avg[7:11])

  plot(data_avg$Time,data_avg$d,type = "p", ylim=g_range_y,xlim = g_range_x, cex = 0.2, pch = 20,
     cex.lab=1.2, xlab = "Generation number",ylab = "The average probabilities",col = "red")

  points(data_avg$Time,data_avg$i,type = "p",col = "blue",cex = 0.2, pch = 20)
  points(data_avg$Time,data_avg$im,type = "p",col = "green",cex = 0.2, pch = 20)
  points(data_avg$Time,data_avg$a,type = "p",col = "orange",cex = 0.3, pch = 20)
  points(data_avg$Time,data_avg$k,type = "p",col = "black",cex = 0.2, pch = 20)
  par(xpd=TRUE)
  legend(g_range_x[1], 1.2*g_range_y[2], c("d","i","im","a","k"), cex=1,col=c("red","blue","green","orange","black"), pch = 20:20, horiz = TRUE)

  rl <-  readline(prompt="This is a plot for Average values of probabilities - Press Enter  ")

# The averaged values of Hallmarks 
  g_range_y <- range(0, 1)
  g_range_x <- range(min(data_avg$Time),time_max)

  plot(data_avg$Time,data_avg$Hd,type = "p", ylim=g_range_y,xlim = g_range_x,cex.lab=1.2, cex = 0.2,
     pch = 20, xlab = "Generation number",ylab = "The averaged Hallmarks values",col = "red")

  points(data_avg$Time,data_avg$Hi,type = "p",col = "blue",cex = 0.2, pch = 20)
  points(data_avg$Time,data_avg$Him,type = "p",col = "green",cex = 0.2, pch = 20)
  points(data_avg$Time,data_avg$Ha,type = "p",col = "orange",cex = 0.2, pch = 20)
  points(data_avg$Time,data_avg$Hb,type = "p",col = "black",cex = 0.2, pch = 20)
  par(xpd=TRUE)
  legend(g_range_x[1], 1.2*g_range_y[2], c("Hd","Hi","Him","Ha","Hb"), cex=1,col=c("red","blue","green","orange","black"), pch = 20:20, horiz = TRUE)

  rl <-  readline(prompt="This is a plot for Average values of Hallmarks - Press Enter  ")

}  # \ analize_data <- function {  






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




    
    

