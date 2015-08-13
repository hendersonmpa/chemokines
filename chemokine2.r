#!/usr/bin/env Rscript
library(randomForest)
library(ggplot2)
library(boot)
library(xtable)
library(reshape)
###################################################
## Functions
## Simple function to select the important varables for the RF used
## with sapply on the data frame
var.select <- function (data) {
  threshold <- 3
  if (data > threshold){
    keep = "Y"
  }else{
    keep = "N"
  }
}

pred.class <- function(data, cutoff1, cutoff2, cutoff3) {
  location <- which(data == max(data))
  max.prob <- data[location]
  if (location == 1 && max.prob > cutoff1) {
    class <- "B"
  } else if (location == 2 && max.prob > cutoff2) {
    class <- "H"
  } else if (location == 3 && max.prob > cutoff3) {
    class <- "O"
  } else {
    class <- "U"
  }
  ##  cat("Location:", location,"Prob:", max.prob, "Class:", class ,"\n")
  class
}

## eff:  takes the class in question (ie "B" "H" or "O" ) the
## rf.votes column for that class and true codes for the dataset
truth.table <- function(data, cutoff, class) {
  prob <- data[1:3]
  code <- data[4]
  pred <- pred.class(prob, cutoff1= cutoff, cutoff2 = cutoff, cutoff3 =
                     cutoff)
  if (code == class && pred == class) {
    truth <- "TP"
  } else if (code != class && pred == class) {
    truth <- "FP"
  } else if (code != class && pred != class) {
    truth <- "TN"
  } else if (code == class && pred != class) {
    truth <- "FN"
  }
  ## cat("Code:", code,"Pred:", pred, "Truth:", truth ,"\n")
  truth
}
##truth <- apply(votes,1, truth.table, cutoff = .5, class = "B")

counter <- function(data , cutoff, class) {
  TP <- 0
  FP <- 0
  FN <- 0
  TN <- 0
  truth <- truth.table(data, cutoff = cutoff, class = class)
  truth
  if (truth == "TP"){
    TP=TP+1
  }
  if (truth == "FP"){
    FP=FP+1
  }
  if (truth == "FN"){
    FN=FN+1
  }
  if (truth == "TN"){
    TN=TN+1
  }
  ## cat("TP:", TP,"TN:",TN, "FP:", FP, "FN:",FN,"\n")
  c(TP,FP,FN,TN)
}
##counter(votes[1,], 0.5, "B")
##c <- apply(votes, 1, counter, cutoff=0.5,class= "B")

eff <- function(cutoff,data, class) {
  TP <- 0
  FP <- 0
  FN <- 0
  TN <- 0
  c <- apply(data, 1, counter, cutoff=cutoff, class=class)
  ## c <- counter(data,cutoff=cutoff, class = class)
  s <- colSums(t(c))
  ##cat("s", s ,"\n")
  TP <- s[1]
  FP <- s[2]
  FN <- s[3]
  TN <- s[4]
  E <- (TP + TN )/(TP + TN + FP + FN) * 100
  temp <- c(cutoff,E,TP,FP,TN,FN)
  ## cat("data:", data,"\n")
  temp
}

code.format <- function(data, class) {
  temp <- as.data.frame(t(data))
  colnames(temp) <- c("x","E","TP","FP","TN","FN")
  temp$code <- rep(class, length(temp[,1]))
  temp
}

Error.fun <- function(data,class){
  observed <- as.character(data[,1])
  predicted <- as.character(data[,2])
  freq <- data[,3]
  wrong <- 0
  total <- 0
  for ( i in 1:length(data[,1])){
    if ( observed[i] == class && predicted[i] != class) {
      wrong <- wrong + freq[i]
    }
    if (observed[i] == class){
      total <- total + freq[i]
    }
    ##cat("total:",total,"wrong:", wrong,"\n")
  }
  error <- wrong/total
  error
}

Error.adj <- function(data,class){
  observed <- as.character(data[,1])
  predicted <- as.character(data[,2])
  freq <- data[,3]
  wrong <- 0
  total <- 0
  for ( i in 1:length(data[,1])){
    if ( observed[i] == class && predicted[i] != class && predicted[i] != "U") {
      wrong <- wrong + freq[i]
    }
    if (observed[i] == class){
      total <- total + freq[i]
    }
    ##cat("total:",total,"wrong:", wrong,"\n")
  }
  error <- wrong/total
  error
}


## sens: Determine the sensitivity for use with boot
sens <- function(type, ind) {
  TP <- 0
  TN <- 0
  FP <- 0
  FN <- 0
  x <- type[ind]
  for ( i in 1:length(x)){
    if (x[i] == "TP" ){
      TP <- TP + 1
    }
    if (x[i] == "TN"){
      TN <- TN + 1
    }
    if (x[i] == "FP"){
      FP = FP + 1
    }
    if (x[i] == "FN"){
      FN = FN + 1
    }
  }
  TP/ (TP + FN) * 100
}

## spec: Determine the specificity for use with boot
spec <- function(type,ind) {
  TP <- 0
  TN <- 0
  FP <- 0
  FN <- 0
  x <- type[ind]
  for ( i in 1:length(x)){
    if (x[i] == "TP" ){
      TP <- TP + 1
    }
    if (x[i] == "TN"){
      TN <- TN + 1
    }
    if (x[i] == "FP"){
      FP = FP + 1
    }
    if (x[i] == "FN"){
      FN = FN + 1
    }
  }
  TN/ (TN + FP) * 100
}


########################################################
## Data import and clean-up
setwd("~/Projects/Research/Chemokines/")
chemokines <- read.csv("Chemokine.csv")

## Command used to trim the added modifiers from the TMN code
## cut -f 3 clinical.csv | sed 's/^\(.\{2\}\)./\1/' > temp
##  paste clinical.txt temp > clinical2.txt
clinical <- read.table("clinical2.txt", header = TRUE)

# Change the wd for figure generation
setwd("./figures")

## Merge the clinical and chemokine data
data.all <- merge(chemokines,clinical, by="ID")
data.uniq <- (unique(data.all))
data.uniq$ID <- factor(data.uniq$ID)

## Select the B H and O codes with no metastasis
bho <- with(data.uniq,data.uniq[code.y != "M" & code.y != "C" , ])
bho$code <- factor(bho$code.y)
bho$M <- factor(bho$M)
bho <- na.omit(bho)


#######################################################################
### Main Body of script
## Create a summary of the population
pop <- table(Group = bho$code, Progression = bho$T)
#pop.xtable <- xtable(pop, caption =  "Tumour staging of the study
#population", label = "tab:pop", table.placement = "h",
#                    caption.placement  = "top", rowname = NULL)
#print(pop.xtable, include.rownames= TRUE , file = "pop.tex")
rownames(pop) <- c("Breast","Healthy","Ovarian")
pop.xtable <- xtable(pop)
print(pop.xtable, include.rownames= TRUE, include.colnames=FALSE,
      hline.after=NULL,only.contents= TRUE, file = "pop.tex")

## Create test and training indices
set.seed(808)
ind <- sample(2, nrow(bho), replace = TRUE, prob=c(0.6, 0.4))
data <- bho[, c(-1,-2,-20,-21,-22)]

## Trainining set RF
bho.rf <- randomForest(code ~ .,
                       data=data[ind == 1,],
                       ntree = 1000,
                       proximity = TRUE,
                       importance=TRUE,
                       keep.forest=TRUE
                       )

# Create a summary table
## Table of the confusion matrix
rownames(bho.rf$confusion) <- c("Breast","Healthy","Ovarian")
all_conf.xtable <- xtable(bho.rf$confusion, digits = c(0,0,0,0,2))
print(all_conf.xtable, include.rownames= TRUE, include.colnames=FALSE,
      hline.after=NULL,only.contents= TRUE, file = "all_conf.tex")

## what are the important variables (via permutation)
##varImpPlot(bho.rf, type=1)
## Rank importance for each class
importance <- data.frame(bho.rf$importance)
var_imp <- importance[rev(order(importance$MeanDecreaseGini)),]
var_imp$names <- rownames(var_imp)
var_imp$names <- with( var_imp, factor( var_imp$names,
                                       levels = var_imp[order(MeanDecreaseGini),]$names))


var_imp$include <- sapply(var_imp$MeanDecreaseGini, var.select)

pdf(file='var_imp.pdf')
bp <- ggplot(var_imp, aes( MeanDecreaseGini,names)) +
  geom_point(aes(shape = include), colour = "black", size = 3) +
  xlab("Mean Decrease in Gini") +
  ylab("")
bp
dev.off()

## Predict classes for the test data
bho.pred <- predict(bho.rf, data[ind == 2,], type="response", proximity=TRUE)
pred.table <- table(observed = data[ind==2, "code"], predicted = bho.pred$predicted)
colnames(pred.table) <- c("Breast","Healthy","Ovarian")
rownames(pred.table) <- c("Breast","Healthy","Ovarian")

######################################
## Use only the important analytes for classification
## See if accuracy is affected
## Get a vector of the variable names to keep

data.trim <- data[ c(as.character(var_imp$names[var_imp$include == "Y"]), "code")]

## A second RF model
trim.rf <- randomForest(code ~ .,
                        data=data.trim[ind == 1,],
                        ntree = 1000,
                        proximity = TRUE,
                        importance=TRUE,
                        keep.forest=TRUE
                        )

rownames(trim.rf$confusion) <- c("Breast","Healthy","Ovarian")
## what are the important variables (via permutation)
## MDS plot of training data
trim.mds <-  cmdscale(1-trim.rf$proximity)
trim.data <- data.frame(observed = as.character(data.trim$code[ind == 1]),
                        predicted = as.character(trim.rf$predicted),
                        trim.mds)

trim.xtable <- xtable(trim.rf$confusion,digits = c(0,0,0,0,2))

print(trim.xtable, include.rownames= TRUE, include.colnames=FALSE,
hline.after=NULL, only.contents=TRUE, file = "trim_conf.tex")

colours <- c("royalblue4","springgreen4","firebrick","gray50")
pdf(file="mds.pdf")
hm <- ggplot(trim.data, aes(X1,X2, colour = observed)) +
  scale_colour_manual(values = colours) +
  geom_point(size = 4, aes(shape = predicted))
hm
dev.off()

## Create a votes dataframe
str(trim.rf$votes)
votes <- round(data.frame(trim.rf$votes),3)
votes$code <- as.factor(data.trim$code[ind==1])
votes$.row <- rownames(votes)
## Use the efficiency function for diagnosis code
## Find the point of maximum efficiency
## x is a list of cut-off values
x <- seq(0.0,1,.001)
code.B <- code.format(sapply(x, eff, data = votes, class = "B"), "B")
cut <- code.B[which(code.B$E == max(code.B$E)),]
B.cut <- cut$x[length(cut$x)] ##median(cut$x) #
B.eff <- code.format(eff(B.cut,data = votes, class = "B"), "B")

code.H <- code.format(sapply(x,eff,data = votes, class = "H"), "H")
cut <- code.H[which(code.H$E == max(code.H$E)),]
H.cut <- cut$x[length(cut$x)] ##median(cut$x)
H.eff <- code.format(eff(H.cut,data = votes, class = "H"), "H")

code.O <- code.format(sapply(x,eff,data = votes, class = "O"), "O")
cut <- code.O[which(code.O$E == max(code.O$E)),]
O.cut <-  cut$x[length(cut$x)] ####[median(cut$x)
O.eff <- code.format(eff(O.cut,data = votes, class = "O"), "O")

## Combine into a single data set
eff.data <- rbind(code.B, code.H, code.O)

## Plot the efficiency data
pdf(file="eff.pdf")
eff.plot <- ggplot(eff.data, aes(x,E)) +
  geom_line( aes(colour = code), size = 1) +
  scale_colour_manual(values = colours) +
  xlab("Threshold") +
  ylab("Percent Efficiency")
eff.plot
dev.off()

## Parallel Coordinates plots for the training data
votes.melt <- melt(votes, id = c( ".row", "code"))
votes.cast <- cast(votes.melt, code ~ variable , median)

## Predict the class of the test data
trim.resp <- predict(trim.rf, data.trim[ind == 2,], type="response",
                     proximity=TRUE)
                                        #table(observed = data.trim[ind==2, "code"], predicted = trim.resp$predicted)
## Get the votes
test.votes <- predict(trim.rf, data.trim[ind == 2,], type="vote",
                      proximity=FALSE)

test.votes <- data.frame(test.votes)
test.votes$code <- as.factor(data.trim$code[ind==2])
test.votes$.row <- rownames(test.votes)

test.votes$pred <- apply(test.votes[,1:3], 1, pred.class , cutoff1 =
                         B.cut, cutoff2= H.cut, cutoff3=O.cut)

## Use grouping function to apply the cut-offs predictions

test.votes$group.B <- apply(test.votes,1, truth.table, cutoff=B.cut, class="B")
test.votes$group.H <- apply(test.votes,1, truth.table, cutoff=H.cut, class="H")
test.votes$group.O <- apply(test.votes,1, truth.table, cutoff=O.cut, class="O")

test.melt <- melt(test.votes, id = c( ".row", "code","pred","group.B" , "group.H", "group.O"))

## Parallel co-ordinates plots of the three classes
colours <- c("black", "black", "gray75","royalblue4")
lines <- c( 4, 2, 3, 1 )
pdf(file = "par_corB.pdf")
sp <- ggplot(test.melt, aes(variable, value, group = .row, colour =
                            group.B , linetype = group.B)) +
    geom_line(size = 1) +
    scale_linetype_manual(values = lines)+
    scale_colour_manual(values = colours) +
    geom_hline(yintercept=B.cut,size = 0.5) +
    scale_x_discrete(breaks=c("B", "H", "O"), labels=c("B", "F", "O")) +
    opts(legend.title = theme_blank(), labels ) +
    xlab("") +
    ylab("Probability")
sp
dev.off()

colours <- c("black","gray75","springgreen4")
threelines <- c( 4, 3, 1 )
pdf(file = "par_corH.pdf")
sp <- ggplot(test.melt, aes(variable, value, group = .row, colour =
                            group.H , linetype = group.H)) +
    geom_line( size = 1) +
    scale_linetype_manual(values = threelines )+
    scale_colour_manual(values = colours) +
    geom_hline(yintercept=c(H.cut),size = 0.5) +
    scale_x_discrete(breaks=c("B", "H", "O"), labels=c("B", "F", "O")) +
    opts(legend.title = theme_blank()) +
    xlab("Class") +
    ylab("")
sp
dev.off()

colours <- c("black", "black", "gray75","firebrick")
pdf(file = "par_corO.pdf")
sp <- ggplot(test.melt, aes(variable, value, group = .row, colour =
                            group.O , linetype = group.O)) +
    geom_line(size = 1) +
    scale_colour_manual(values = colours) +
    scale_linetype_manual(values = lines) +
    geom_hline(yintercept=c(O.cut),size = 0.5) +
    scale_x_discrete(breaks=c("B", "H", "O"), labels=c("B", "F", "O")) +
    opts(legend.title = theme_blank()) +
    xlab("") +
    ylab("")
sp
dev.off()



## Create a confusion matrix for the test data
test.table <- table( observed = test.votes$code, predicted = test.votes$pred)
test.table.df <- as.data.frame(test.table)
B.error <- Error.fun(test.table.df, "B")
H.error <- Error.fun(test.table.df, "H")
O.error <- Error.fun(test.table.df, "O")
##U.error <- 1
error <- c(B.error,H.error,O.error) # ,U.error)
# Calculate the "adjusted error"
B.adj <- Error.adj(test.table.df, "B")
H.adj <- Error.adj(test.table.df, "H")
O.adj <- Error.adj(test.table.df, "O")
adjust <- c(B.adj,H.adj,O.adj) # ,U.adj)
test.table <- cbind(test.table,error,adjust)
rownames(test.table) <- c("Breast","Cancer free","Ovarian")
# Create the test data latex table
test.xtable <- xtable(test.table, digits = c(0,0,0,0,0,2,2))
print(test.xtable, include.rownames= TRUE, include.colnames=FALSE,
hline.after=NULL, only.contents=TRUE, file = "test_conf.tex")

## Boot strapped sensitivity and specificity
## B
Bsn.boot <- boot(test.votes$group.B, sens, R=1000)
##plot(Bsn.boot)
Bsn.ci <- boot.ci(Bsn.boot,conf = 0.95,type = "all")
Bsp.boot <- boot(test.votes$group.B, spec, R=1000)
Bsp.ci <- boot.ci(Bsp.boot,conf = 0.95,type = "all")

## H
Hsn.boot <- boot(test.votes$group.H, sens, R=1000)
Hsn.ci <- boot.ci(Hsn.boot,conf = 0.95,type = "all")
Hsp.boot <- boot(test.votes$group.H, spec, R=1000)
Hsp.ci <- boot.ci(Hsp.boot,conf = 0.95,type = "all")

## O
Osn.boot <- boot(test.votes$group.O, sens, R=1000)
Osn.ci <- boot.ci(Osn.boot,conf = 0.95,type = "all")
Osp.boot <- boot(test.votes$group.O, spec, R=1000)
Osp.ci <- boot.ci(Osp.boot,conf = 0.95,type = "all")

## Create a table with n, efficiency, cut-off, sens (ci), spec (ci
Bn.test <- length(test.votes$B[test.votes$code == "B"])
Bn.train <- length(votes$B[votes$code == "B"])
B <- c( Bn.train, Bn.test, B.eff[1:2],  Bsn.boot$t0,
       round(Bsn.ci$bca[4:5],1) ,  Bsp.boot$t0, round(Bsp.ci$bca[4:5],1))

Hn.test <- length(test.votes$H[test.votes$code == "H"])
Hn.train <- length(votes$H[votes$code == "H"])
H <- c( Hn.train, Hn.test, H.eff[1:2], Hsn.boot$t0,
      round(Hsn.ci$bca[4:5],1) ,  Hsp.boot$t0,  c("*","*"))

On.test <- length(test.votes$O[test.votes$code == "O"])
On.train <- length(votes$O[votes$code == "O"])
O <- c(On.train, On.test ,O.eff[1:2],  Osn.boot$t0,
       round(Osn.ci$bca[4:5],1) ,  Osp.boot$t0, round(Osp.ci$bca[4:5],1))

results <- rbind(B, H, O)
senCI <- paste(results[,6],results[,7], sep=" - ")
spCI <- paste(results[,9],results[,10], sep=" - ")
results.table <- cbind(results[,1:5], senCI,unlist(results[,8]), spCI)
colnames(results.table) <- c("Train (n)", "Test (n)", "Threshold",
                              "Efficiency", "Sensitivity", "95% CI", "Specificity", "95% CI" )

rownames(results.table) <- c("Breast cancer", "Cancer free", "Ovarian cancer")

results.xtable <- xtable(results.table ,caption =
                         "Summary of the Random Forest algorithm classification
  accuracy using the optimal classification threshold. Bootstrapped
  confidence intervals are provided. The asterisk indicates that
  intervals could not be calculated as there was no variation between the
  bootstrapped data sets.",label = "tab:summary", align = "lrrrrrrrr", digits = c(0,0,0,3,1,1,1,1,1))

print(results.xtable, include.rownames=TRUE,  table.placement = "h", caption.placement = "bottom",file = "summary_table.tex")

### Examine the characteristics of the misidentified classes
## B = 3 and O = 2

## Create an rda file
save(bho, pop, bho.rf, var_imp, trim.rf, ind , trim.data, eff.data,
     B.cut, H.cut, O.cut, test.melt, test.table, results, results.table, file="chemokine.rda")
