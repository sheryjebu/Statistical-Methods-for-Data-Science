#Set the work directory
setwd("~/MSC/Stat/CW/final_codendoc")
# Load the dataset
library(readxl)
gene_data <- read_excel("gene_data.xlsx",col_names =FALSE)
View(gene_data)

#-----------------------------------#
#Task 1: Preliminary data analysis
#-----------------------------------#
gene_data = as.matrix(gene_data)
X = gene_data[,1]
Y_x1 = gene_data[,2]
Y_x2 = gene_data[,3]
Y_x3 = gene_data[,4]
Y_x4 = gene_data[,5]
Y_x5 = gene_data[,6]

#Time series plots (of each gene against sampling time
plot(X,Y_x1,type= "l",
     xlab="Sampling Time",
     ylab="X1 Gene data", 
     main="Time series plots with X1 gene data"
     )
plot(X,Y_x2,type= "l",
     xlab="Sampling Time",
     ylab="X2 Gene data", 
     main="Time series plots with X2 gene data"
     )
plot(X,Y_x3,type= "l",
     xlab="Sampling Time",
     ylab="X3 Gene data", 
     main="Time series plots with X3 gene data"
     )
plot(X,Y_x4,type= "l",
     xlab="Sampling Time",
     ylab="X4 Gene data", 
     main="Time series plots with X4 gene data"
     )
plot(X,Y_x5,type= "l",
     xlab="Sampling Time",
     ylab="X5 Gene data", 
     main="Time series plots with X5 gene data"
     )


ts_gene_dat = ts(gene_data[,2:6], frequency = 10.0)
plot.ts(ts_gene_dat)

#Distribution for each gene (time-series)
hist(Y_x1, col = "Grey")
abline(v = mean(Y_x1), lwd = 5 , col="Brown")
abline(v = median(Y_x1), lwd = 3 , col = "green")

mode =+ unique(Y_x1)
mode = mode[which.max(tabulate(match(Y_x1, mode)))]
abline(v = mode, lwd = 4 , col="red")

#Plot matrix, consisting of scatter plots corresponding to each data frame
install.packages("corrplot")
install.packages("RColorBrewer")
install.packages("psych")

library(corrplot)
library(psych)
cols_combo = c("#00AFBB","#E7B800","#FC4E07")

pairs(gene_data[,2:6])

pairs.panels(gene_data[,2:6],
             method   = "pearson", #correlation method
             hist.col = "#00AFBB",
             density  = TRUE,     # show density plots
             ellipses = TRUE      #show correlation ellipses
             )

pairs(gene_data[,2:6],                     # Data frame of variables
      labels = colnames(gene_data[,2:6]),  # Variable names
      pch = 21,                            # Pch symbol
      bg = cols_combo,                     # Background color of the symbol (pch 21 to 25)
      col = cols_combo,                    # Border color of the symbol
      main = "Correlation of Scatter Plots",    # Title of the plot
      row1attop = TRUE,                         # If FALSE, changes the direction of the diagonal
      gap = 1,                                  # Distance between subplots
      cex.labels = NULL,                        # Size of the diagonal text
      font.labels = 1)                          # Font style of the diagonal text
