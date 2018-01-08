#Project 3
#Library
library(MASS)
library(mvtnorm)
library(rpart)
library(class)
library(nnet)
#Question 1
## Load data
p3data20<- read.table("http://wwwf.imperial.ac.uk/~eakc07/S7/data20.dat")
#p3data20new <- p3data20[complete.cases(p3data20),]
#Data treatment
#range of different columns
range(p3data20,na.rm=TRUE)
#summary of every feature vector sets
summary(p3data20)
#standard deviation
sd <- apply(p3data20,2,sd,na.rm=TRUE)
#Boxplots of each dimensions of data

boxplot((p3data20), main = "Boxplot of each dimensions",ylab="")
par(mfrow=c(2,2))
for (i in 2:5) {title <- "Boxplot of the dimension "
                boxplot((p3data20[,i]), main = paste(title,i),ylab="")}
for (i in 6:9) {title <- "Boxplot of the dimension "
                boxplot((p3data20[,i]), main = paste(title,i),ylab="Values of observations")}
for (i in 10:13) {title <- "Boxplot of the dimension "
                boxplot((p3data20[,i]), main = paste(title,i),ylab="Values of observations")}
for (i in 14:17) {title <- "Boxplot of the dimension "
                boxplot((p3data20[,i]), main = paste(title,i),ylab="Values of observations")}
for (i in 18:21) {title <- "Boxplot of the dimension "
                  boxplot((p3data20[,i]), main = paste(title,i),ylab="Values of observations")}
for (i in 22:25) {title <- "Boxplot of the dimension "
                  boxplot((p3data20[,i]), main = paste(title,i),ylab="Values of observations")}
for (i in 26:29) {title <- "Boxplot of the dimension "
                  boxplot((p3data20[,i]), main = paste(title,i),ylab="Values of observations")}

#scatter plots of each dimensions of data
par (mfrow=c(2,2))
for (i in 2:5) {title <- "Scatter plot of the dimension"
  plot(p3data20[(which(p3data20[,1]==0)), i], col = "red" , main = paste(title,i),ylab="Values of observations")
  points (p3data20[(which(p3data20[,1]==1)), i], col = "blue")
}
for (i in 6:9) {title <- "Scatter plot of the dimension"
                plot(p3data20[(which(p3data20[,1]==0)), i], col = "red" , main = paste(title,i),ylab="Values of observations")
                points (p3data20[(which(p3data20[,1]==1)), i], col = "blue")
}
for (i in 10:13) {title <- "Scatter plot of the dimension"
                plot(p3data20[(which(p3data20[,1]==0)), i], col = "red" , main = paste(title,i),ylab="Values of observations")
                points (p3data20[(which(p3data20[,1]==1)), i], col = "blue")
}
for (i in 14:17) {title <- "Scatter plot of the dimension"
                plot(p3data20[(which(p3data20[,1]==0)), i], col = "red" , main = paste(title,i),ylab="Values of observations")
                points (p3data20[(which(p3data20[,1]==1)), i], col = "blue")
}
for (i in 18:21) {title <- "Scatter plot of the dimension"
                  plot(p3data20[(which(p3data20[,1]==0)), i], col = "red" , main = paste(title,i),ylab="Values of observations")
                  points (p3data20[(which(p3data20[,1]==1)), i], col = "blue")
}
for (i in 22:25) {title <- "Scatter plot of the dimension"
                  plot(p3data20[(which(p3data20[,1]==0)), i], col = "red" , main = paste(title,i),ylab="Values of observations")
                  points (p3data20[(which(p3data20[,1]==1)), i], col = "blue")
}
for (i in 26:29) {title <- "Scatter plot of the dimension"
                  plot(p3data20[(which(p3data20[,1]==0)), i], col = "red" , main = paste(title,i),ylab="Values of observations")
                  points (p3data20[(which(p3data20[,1]==1)), i], col = "blue")
}

#Histogram for each dimension
par(mfrow=c(2,2))
for (i in 2:5) {title <- "Histogram of the dimension "
                hist((p3data20[which(p3data20[,1]==0),i]), 
                     col="red",main = paste(title,i),ylab="Frequency",
                     xlab="Values of observations")
                hist((p3data20[which(p3data20[,1]==1),i]), 
                     col="blue",main = paste(title,i),ylab="Frequency",
                     xlab="Values of observations",add=TRUE)}
par(mfrow=c(2,2))
for (i in 6:9) {title <- "Histogram of the dimension "
                hist((p3data20[which(p3data20[,1]==0),i]), 
                     col="red",main = paste(title,i),ylab="Frequency",
                     xlab="Values of observations")
                hist((p3data20[which(p3data20[,1]==1),i]), 
                     col="blue",main = paste(title,i),ylab="Frequency",
                     xlab="Values of observations",add=TRUE)}
par(mfrow=c(2,2))
for (i in 10:13) {title <- "Histogram of the dimension "
                hist((p3data20[which(p3data20[,1]==0),i]), 
                     col="red",main = paste(title,i),ylab="Frequency",
                     xlab="Values of observations")
                hist((p3data20[which(p3data20[,1]==1),i]), 
                     col="blue",main = paste(title,i),ylab="Frequency",
                     xlab="Values of observations",add=TRUE)}
par(mfrow=c(2,2))
for (i in 14:17) {title <- "Histogram of the dimension "
                hist((p3data20[which(p3data20[,1]==0),i]), 
                     col="red",main = paste(title,i),ylab="Frequency",
                     xlab="Values of observations")
                hist((p3data20[which(p3data20[,1]==1),i]), 
                     col="blue",main = paste(title,i),ylab="Frequency",
                     xlab="Values of observations",add=TRUE)}
par(mfrow=c(2,2))
for (i in 18:21) {title <- "Histogram of the dimension "
                hist((p3data20[which(p3data20[,1]==0),i]), 
                     col="red",main = paste(title,i),ylab="Frequency",
                     xlab="Values of observations")
                hist((p3data20[which(p3data20[,1]==1),i]), 
                     col="blue",main = paste(title,i),ylab="Frequency",
                     xlab="Values of observations",add=TRUE)}
par(mfrow=c(2,2))
for (i in 22:25) {title <- "Histogram of the dimension "
                hist((p3data20[which(p3data20[,1]==0),i]), 
                     col="red",main = paste(title,i),ylab="Frequency",
                     xlab="Values of observations")
                hist((p3data20[which(p3data20[,1]==1),i]), 
                     col="blue",main = paste(title,i),ylab="Frequency",
                     xlab="Values of observations",add=TRUE)}
par(mfrow=c(2,2))
for (i in 26:29) {title <- "Histogram of the dimension "
                hist((p3data20[which(p3data20[,1]==0),i]), 
                     col="red",main = paste(title,i),ylab="Frequency",
                     xlab="Values of observations")
                hist((p3data20[which(p3data20[,1]==1),i]), 
                     col="blue",main = paste(title,i),ylab="Frequency",
                     xlab="Values of observations",add=TRUE)}
#Question 2
#Treatments for missing values (NA)
#function for replace NA values of a column to its column mean
rpmean <- function(x) 
{x[which(is.na(x))] <- mean(x,na.rm=TRUE)
 x
}
#function for replace NA values of a column to its column mode
rpmode <- function(x) 
{ fre <- as.matrix(summary(factor(x)))
  w <- which(fre==max(fre))
  colmode <- as.numeric(names(fre[w,1]))
  x[which(is.na(x))] <- colmode
  x
}
#define new set of data without any missing values NA 
p3data20new <- matrix(,nrow(p3data20),ncol(p3data20))
#replace NA values in each dimension with its corresponding mean or mode of its dimension's observations
#The classes corresponding to each observations
p3data20new[,1] <- p3data20[,1]
#dimension 5 has discrete datas, which I choose to replace the NA by the mode of this dimension's observations
p3data20new[,5] <- rpmode(p3data20[,5])
#for the rest of the dimensions with continuous datas, replace the NA by the mean of its dimension's observations
for(i in 2:4){
  p3data20new[,i] <- rpmean(p3data20[,i])
}
for(i in 6:29){
  p3data20new[,i] <- rpmean(p3data20[,i])
}

#Split the data into training set and test set
#Test set has size of 100 and is reserved for later sections
#Shuffle the data
set.seed(135)
p3data20newr<-p3data20new[sample(nrow(p3data20new)),]
split  <- sort(sample(nrow(p3data20newr),100,FALSE))
data.tr <- p3data20newr[-split,]
data.te <- p3data20newr[split,]

#Data classification methods
#KNN
#normalise data
p3data20newr.knn <- matrix(NA,nrow(p3data20newr),ncol(p3data20newr))
p3data20newr.knn[,1] <- p3data20newr[,1]
for(j in 2:29){
  p3data20newr.knn[,j] <- (p3data20newr[,j]-mean(p3data20newr[,j]))/sd(p3data20newr[,j])
}
data.tr.knn <- p3data20newr.knn[-split,]
data.te.knn <- p3data20newr.knn[split,]
#Sequential forward feature selection by 10-fold cross validation
set.seed(135)
fold<- cut(1:nrow(data.tr.knn),breaks=10,labels=FALSE)
er.knn.store <-matrix(NA,1,28)
#First, find out the best individual performing feature
ini.feat.knn<-c()
for(j in 2:29){
  temp.er.knn<-c()
  for (i in 1:10){
    data.tr.part <- data.tr.knn[which(fold==i),]
    data.tr.main <- data.tr.knn[-which(fold==i),]
pred <- knn(as.matrix(data.tr.main[,j]),as.matrix(data.tr.part[,j]),as.matrix(data.tr.main[,1]))
temp.er.knn[i] <- (sum(pred != data.tr.part[,1]))/nrow(data.tr.part)
}
er.knn.store[j-1] <- mean(temp.er.knn)
}
min.er.knn <- min(er.knn.store)
#the initial selected feature
ini.feat.knn <- which(er.knn.store==min.er.knn)+1
#initialise algorithm
f.knn <- ini.feat.knn
no.feat <- 1
limit <- 1
FS.knn <- matrix(0,28,29)
er.knn.feat.store <- matrix(NA,2,28)
while(min.er.knn <= limit && no.feat < 29){
  #range.feat <- c(2:29)[-(ini.feat.knn-1)]  
  FS.knn[,no.feat] <- f.knn
  FS.knn[no.feat,29] <- min.er.knn
  limit <- min.er.knn
  #ini.feat.knn.1 <- ini.feat.knn
for(j in 2:29){
  temp.er.knn<-c()
  f <- unique(c(ini.feat.knn,j))
  for(i in 1:10){
data.tr.part <- data.tr.knn[which(fold==i),]
data.tr.main <- data.tr.knn[-which(fold==i),]
pred <- knn(as.matrix(data.tr.main[,f]),as.matrix(data.tr.part[,f]),as.matrix(data.tr.main[,1]))
temp.er.knn[i]<- (sum(pred != data.tr.part[,1]))/nrow(data.tr.part)
}
er.knn.feat.store[1,j-1] <- j
er.knn.feat.store[2,j-1] <- mean(temp.er.knn)
}
for(h in 1:length(ini.feat.knn)){
  er.knn.feat.store <- er.knn.feat.store[,-(which(er.knn.feat.store[1,]==ini.feat.knn[h]))]}
min.er.knn <- min(er.knn.feat.store[2,])
f.knn <- er.knn.feat.store[1,which(er.knn.feat.store[2,]==min.er.knn)]
if(length(f.knn)>1){f.knn <- f.knn[1]}
ini.feat.knn<- c(ini.feat.knn,f.knn)
no.feat <- no.feat +1
er.knn.feat.store <- matrix(NA,2,28)
}
#selected features
ini.feat.knn
table.feat.knn <- cbind(FS.knn[1,1:length(ini.feat.knn)],FS.knn[1:length(ini.feat.knn),29])
#plot(table.feat.knn[,2],xlab="Number of dimensions of the dataset",
 #    ylab="error rate",
  #   main=paste("Performance of the KNN classifier \n"," as number of dimensions of dataset increases"))
#parameter selection by 10-fold cross validation
er.knn.par.store <- matrix(NA,10,10)
k <- c(3,7,11,15,19,23,27,31,35,39)
for(i in 1:10){
  for(d in 1:10){
    data.tr.part <- data.tr.knn[which(fold==i),]
    data.tr.main <- data.tr.knn[-which(fold==i),]
    pred <- knn(as.matrix(data.tr.main[,ini.feat.knn]),as.matrix(data.tr.part[,ini.feat.knn]),as.matrix(data.tr.main[,1]),k=k[d])
    er.knn.par.store[i,d] <- (sum(pred != data.tr.part[,1]))/nrow(data.tr.part)
  }
}
er.knn.par.store <- colMeans(er.knn.par.store)
min.par.er.knn <- min(er.knn.par.store)
ini.par.knn <- which(er.knn.par.store==min.par.er.knn)
k <- k[ini.par.knn]   
k
#Assessment of performance of KNN based on training set 
#with selected features and parameters with 10-fold cross validation
fold<- cut(1:nrow(data.tr.knn),breaks=10,labels=FALSE)
er.perform.tr.knn <- c()
for(i in 1:10){
    data.tr.part <- data.tr.knn[which(fold==i),]
    data.tr.main <- data.tr.knn[-which(fold==i),]
    pred <- knn(as.matrix(data.tr.main[,ini.feat.knn]),as.matrix(data.tr.part[,ini.feat.knn]),as.matrix(data.tr.main[,1]),k=k)
    er.perform.tr.knn[i] <- (sum(pred != data.tr.part[,1]))/nrow(data.tr.part)
  }
er.perform.tr.knn <- mean(er.perform.tr.knn)
#QDA
#Sequential backward feature selection
#Shuffle the data
set.seed(12)
p3data20newr<-p3data20new[sample(nrow(p3data20new)),]
fold <- cut(1:nrow(data.tr),breaks=10,labels=FALSE)
er.qda.store <-matrix(NA,1,28)
#First, take out features from data one by one and do cross validation
for(j in 2:29){
  temp.er.qda <- c()
  for (i in 1:10){
    data.tr.part <- data.tr[which(fold==i),]
    data.tr.main <- data.tr[-which(fold==i),]
    data.tr.part.cl <- data.tr.part[,1]
    data.tr.part.da <- data.tr.part[,2:29]
    data.tr.main.cl <- data.tr.main[,1]
    data.tr.main.da <- data.tr.main[,2:29]
    qda.jj1 <- qda(data.tr.main.da[,-(j-1)],data.tr.main.cl)
    qda.jj1.pred <- predict(qda.jj1,data.tr.part.da[,-(j-1)])$class
    temp.er.qda[i] <-sum(qda.jj1.pred != data.tr.part.cl)/(nrow(data.tr.part))
  }
  er.qda.store[j-1] <- mean(temp.er.qda)
}
min.er.qda <- min(er.qda.store)
#the initial selected feature to be removed from the feature set
ini.feat.qda <- which(er.qda.store==min.er.qda)+1
#compute the corresponding error rate for the feature subset with the initial selected feature
f.qda <- ini.feat.qda
no.feat <- 1
limit <- 1
FS.qda <- matrix(0,28,29)
er.qda.feat.store <- matrix(NA,2,28)
while(min.er.qda <= limit&& no.feat < 29){ 
  FS.qda[,no.feat] <- f.qda
  FS.qda[no.feat,29] <- min.er.qda
  limit <- min.er.qda
  for(j in 2:29){
    temp.er.qda<-c()
    f <- unique(c(ini.feat.qda,j))
    for(i in 1:10){
      data.tr.part <- data.tr[which(fold==i),]
      data.tr.main <- data.tr[-which(fold==i),]
      data.tr.part.cl <- data.tr.part[,1]
      data.tr.part.da <- data.tr.part[,2:29]
      data.tr.main.cl <- data.tr.main[,1]
      data.tr.main.da <- data.tr.main[,2:29]
      qda.jj1 <- qda(data.tr.main.da[,-(f-1)],data.tr.main.cl)
      qda.jj1.pred <- predict(qda.jj1,data.tr.part.da[,-(f-1)])$class
      temp.er.qda[i] <-sum(qda.jj1.pred != data.tr.part.cl)/(nrow(data.tr.part))
    }
    er.qda.feat.store[1,j-1] <- j
    er.qda.feat.store[2,j-1] <- mean(temp.er.qda)
  }
  for(h in 1:length(ini.feat.qda)){
    er.qda.feat.store <- er.qda.feat.store[,-(which(er.qda.feat.store[1,]==ini.feat.qda[h]))]}
  min.er.qda <- min(er.qda.feat.store[2,])
  f.qda <- er.qda.feat.store[1,which(er.qda.feat.store[2,]==min.er.qda)]
  if(length(f.qda)>1){f.qda <- f.qda[1]}
  ini.feat.qda<- c(ini.feat.qda,f.qda)
  no.feat <- no.feat +1
  er.qda.feat.store <- matrix(NA,2,28)
}
#selected features to be get rid of the feature set
ini.feat.qda
#table.feat.qda <- cbind(FS.qda[1,1:length(ini.feat.qda)],FS.qda[1:length(ini.feat.qda),29])
#plot(table.feat.qda[,2],xlab="Number of dimensions ommitted",
     #ylab="error rate",
     #main=paste("Performance of the QDA classifier \n"," as dimensions being ommitted increases"))
#Assessment of performance of QDA analysis based on training set 
#with selected features with 10-fold cross validation
fold <- cut(1:nrow(data.tr),breaks=10,labels=FALSE)
er.perform.tr.qda <- c()
for(i in 1:10){
  data.tr.part <- data.tr[which(fold==i),]
  data.tr.main <- data.tr[-which(fold==i),]
  data.tr.part.cl <- data.tr.part[,1]
  data.tr.part.da <- data.tr.part[,2:29]
  data.tr.main.cl <- data.tr.main[,1]
  data.tr.main.da <- data.tr.main[,2:29]
  qda.jj1 <- qda(data.tr.main.da[,-(ini.feat.qda-1)],data.tr.main[,1])
  qda.jj1.pred <- predict(qda.jj1,data.tr.part.da[,-(ini.feat.qda-1)])$class
  er.perform.tr.qda[i] <- sum(qda.jj1.pred != data.tr.part.cl)/(nrow(data.tr.part.da))
}
er.perform.tr.qda <- mean(er.perform.tr.qda)
dev.off()

#Classification and Regression Tree
set.seed(135)
#Shuffle the data
p3data20newr<-p3data20new[sample(nrow(p3data20new)),]
fold <- cut(1:nrow(data.tr),breaks=10,labels=FALSE)
i <- 1
data.tr.part <- data.tr[which(fold==i),]
data.tr.main <- data.tr[-which(fold==i),]
#er.tree.store <- c()
#Test with different approach in finding cp values by experimenting on partition of training data
#Gini index
jj1 <- rpart(as.factor(data.tr.main[,1])~.,data=as.data.frame(data.tr.main[,2:29]),parms=list(split="gini"),
             control=rpart.control(xval=nrow(data.tr.main[,2:29])))
plotcp(jj1)
plot(jj1)
text(jj1)
#The default cp value is 0.014
jj2 <- prune(jj1,0.014)
plot(jj2)
text(jj2)
pred.tree.cl<-c()
pred.tree <- predict(jj2,as.data.frame(data.tr.part[,2:29]))
pred.tree.cl[which(pred.tree[,2]>0.5)] <-1
pred.tree.cl[which(pred.tree[,1]>0.5)] <-0
er.tree.default <- sum(pred.tree.cl != data.tr.part[,1])/nrow(data.tr.part)
print(er.tree.default)
#error rate is 0.138 for default cp
#Try to find the appropriate cp value by One-SE rule 
min.jj1.tree <-min(jj1$cptable[,4])
std.jj1 <- jj1$cptable[which.min(jj1$cptable[,4]),5]
cp.jj1 <- max(jj1$cptable[which(jj1$cptable[,4] < std.jj1+min.jj1.tree),1])
jj3<- prune(jj1, cp=cp.jj1)
plot(jj3)
text(jj3)
pred.tree.cl<-c()
pred.tree <- predict(jj3,as.data.frame(data.tr.part[,2:29]))
pred.tree.cl[which(pred.tree[,2]>0.5)] <-1
pred.tree.cl[which(pred.tree[,1]>0.5)] <-0
er.tree.cp <- sum(pred.tree.cl != data.tr.part[,1])/nrow(data.tr.part)
print(er.tree.cp)
#choose cp.jj1 which is 0.01
#One-SE rule does not require manuel search of cp-value, which is preferred and more convenient in computation
#Assessment of performance of Tree-based analysis based on training set
#to find the average cp-value by 10-fold cross validation and One-SE rule
jj.tr.tree <- rpart(as.factor(data.tr[,1])~.,data=as.data.frame(data.tr[,2:29]),parms=list(split="gini"),
                         control=rpart.control(xval=nrow(data.tr[,2:29])))
plot(jj.tr.tree)
min.jj.tree <-min(jj.tr.tree$cptable[,4])
std.jj <- jj.tr.tree$cptable[which.min(jj.tr.tree$cptable[,4]),5]
cp.jj.tree <- max(jj.tr.tree$cptable[which(jj.tr.tree$cptable[,4] < std.jj+min.jj.tree),1])
prune.tree<- prune(jj.tr.tree, cp=cp.jj.tree)
plot(prune.tree)
text(prune.tree)
er.perform.tr.tree <- c()
fold <- cut(1:nrow(data.tr),breaks=10,labels=FALSE)
for(i in 1:10){
data.tr.part <- data.tr[which(fold==i),]
data.tr.main <- data.tr[-which(fold==i),]
pred.tree.cl<-c()
pred.tree <- predict(jj.tr.tree,as.data.frame(data.tr.part[,2:29]))
pred.tree.cl[which(pred.tree[,2]>0.5)] <-1
pred.tree.cl[which(pred.tree[,1]>0.5)] <-0
er.perform.tr.tree[i] <- sum(pred.tree.cl != data.tr.part[,1])/nrow(data.tr.part)}
er.perform.tr.tree <- mean(er.perform.tr.tree)
er.perform.tr.tree 
#MLP
#Sequential forward feature selection
set.seed(132)
#scale data
p3data20newr.mlp <- matrix(NA,nrow(p3data20newr),ncol(p3data20newr))
p3data20newr.mlp[,1] <- p3data20newr[,1]
for(j in 2:29){
  p3data20newr.mlp[,j] <- (p3data20newr[,j])/max(p3data20newr[,j])
}
data.tr.mlp <- p3data20newr.mlp[-split,]
data.te.mlp <- p3data20newr.mlp[split,]
#partition for 5 fold cross validation for feature selection since computationally expensive
fold<- cut(1:nrow(data.tr.mlp),breaks=10,labels=FALSE)
er.mlp.store <-matrix(NA,1,28)
ini.feat.mlp<-c()
for(j in 2:29){
  temp.er.mlp<-c()
  for (i in 1:5){
    data.tr.part <- data.tr.mlp[which(fold==i),]
    data.tr.main <- data.tr.mlp[-which(fold==i),]
    data.tr.part.cl <- data.tr.part[,1]
    data.tr.part.da <- as.matrix(data.tr.part[,j])
    data.tr.main.cl <- data.tr.main[,1]
    data.tr.main.da <- as.matrix(data.tr.main[,j])
    score <- numeric(5)
    count <- 0
    while(count < 5) {
      jj1 <- nnet(as.factor(data.tr.main.cl) ~ ., data = data.tr.main.da, size = 6,
                  decay = 0.01, maxit = 1000)
      evs <- eigen(nnetHess(jj1, data.tr.main.da, data.tr.main.cl),
                   T)$values
      if(min(evs) > 1e-005) {
        count <- count + 1
        preds <- predict(jj1, data.tr.part.da)
        score[count] <- sum((preds > 0.5) != data.tr.part.cl)/nrow(data.tr.part.da)
        
      }
      wm <- which.min(score)
      temp.er.mlp[i] <- score[wm]
    }}
  er.mlp.store[j-1] <- mean(temp.er.mlp)}
min.er.mlp <- min(er.mlp.store)
#the initial selected feature
ini.feat.mlp <- which(er.mlp.store==min.er.mlp)+1
#compute the corresponding error rate for the feature subset with the initial selected feature
#CAUTION: long running time 
f.mlp <- ini.feat.mlp
no.feat <- 1
limit <- 1
FS.mlp <- matrix(0,28,29)
er.mlp.feat.store <- matrix(NA,2,28)
#set condition for iteration, limit+0.05 for finding the actual minimum
# to observethe general trend of error rates after each feature is being added
while(min.er.mlp < limit+0.05 && no.feat < 11){
  FS.mlp[,no.feat] <- f.mlp
  FS.mlp[no.feat,29] <- min.er.mlp
  limit <- min.er.mlp
  for(j in 2:29){
    temp.er.mlp<-c()
    f <- unique(c(ini.feat.mlp,j))
    for (i in 1:5){
      data.tr.part <- data.tr.mlp[which(fold==i),]
      data.tr.main <- data.tr.mlp[-which(fold==i),]
      data.tr.part.cl <- data.tr.part[,1]
      data.tr.part.da <- as.matrix(data.tr.part[,f])
      data.tr.main.cl <- data.tr.main[,1]
      data.tr.main.da <- as.matrix(data.tr.main[,f])
      score <- numeric(3)
      count <- 0
      while(count < 3) {
        jj1 <- nnet(as.factor(data.tr.main.cl) ~ ., data = data.tr.main.da, size = 6,
                    decay = 0.01, maxit = 1000)
        evs <- eigen(nnetHess(jj1, data.tr.main.da, data.tr.main.cl),
                     T)$values
        if(min(evs) > 1e-005) {
          count <- count + 1
          preds <- predict(jj1, data.tr.part.da)
          score[count] <- sum((preds > 0.5) != data.tr.part.cl)/nrow(data.tr.part.da)
        }
        wm <- which.min(score)
        temp.er.mlp[i] <- score[wm]
      }}
    er.mlp.feat.store[1,j-1] <- j
    er.mlp.feat.store[2,j-1] <- mean(temp.er.mlp)
  }
  for(h in 1:length(ini.feat.mlp)){
    er.mlp.feat.store <- er.mlp.feat.store[,-(which(er.mlp.feat.store[1,]==ini.feat.mlp[h]))]}
  min.er.mlp <- min(er.mlp.feat.store[2,])
  f.mlp <- er.mlp.feat.store[1,which(er.mlp.feat.store[2,]==min.er.mlp)]
  if(length(f.mlp)>1){f.mlp <- f.mlp[1]}
  ini.feat.mlp <- c(ini.feat.mlp,f.mlp)
  no.feat <- no.feat +1
  er.mlp.feat.store <- matrix(NA,2,28)
}
#features added to the training data one by one with their corresponding error rate
ini.feat.mlp
table.feat.mlp <- cbind(FS.mlp[1,1:length(ini.feat.mlp)],FS.mlp[1:length(ini.feat.mlp),29])
plot(table.feat.mlp[,2],xlab="Number of dimensions of the dataset",
     ylab="error rate",
     main=paste("Performance of the MLP classifier \n"," as number of dimensions of dataset increases"))
#minimum error rate occurs when the 10th feature (feature set 7) is added
#included features: 2 13 10 18 14 17  5 24  4  6
ini.feat.mlp <- c(2,13,10,18,14,17,5,24,4,6)
#parameter selection by 5-fold cross validation
fold <- cut(1:nrow(data.tr),breaks=5,labels=FALSE)
H <- c(3,4,5,6,7,8,9,10,11)
wd <- c(0.1,0.01,0.001)
gr <- as.matrix(expand.grid(wd,H))
"nnet.f" <-
  function(x)
  {
    H <- x[2]
    decay <- x[1]
    temp.er.mlp<-c()
    for (i in 1:5){
      data.tr.part <- data.tr.mlp[which(fold==i),]
      data.tr.main <- data.tr.mlp[-which(fold==i),]
      data.tr.part.cl <- data.tr.part[,1]
      data.tr.part.da <- as.matrix(data.tr.part[,ini.feat.mlp])
      data.tr.main.cl <- data.tr.main[,1]
      data.tr.main.da <- as.matrix(data.tr.main[,ini.feat.mlp])
      score <- numeric(3)
      count <- 0
      while(count < 3) {
        jj1 <- nnet(as.factor(data.tr.main.cl) ~ ., data = data.tr.main.da, size = 6,
                    decay = 0.01, maxit = 1000)
        evs <- eigen(nnetHess(jj1, data.tr.main.da, data.tr.main.cl),
                     T)$values
        if(min(evs) > 1e-005) {
          count <- count + 1
          preds <- predict(jj1, data.tr.part.da)
          score[count] <- sum((preds > 0.5) != data.tr.part.cl)/nrow(data.tr.part.da)
        }
        wm <- which.min(score)
        temp.er.mlp[i] <- score[wm]
      }
    }
    mean(temp.er.mlp)
  }
er.mlp.par.store <- apply(gr,1,nnet.f)
decay <- as.numeric(gr[which.min(er.mlp.par.store),1])
H<- as.numeric(gr[which.min(er.mlp.par.store),2])
#parameter selected: weight decay = 0.01 and number of hidden nodes H = 8
#Assessment of performance of MLP analysis based on training set 
#with selected features with 10-fold cross validation
fold <- cut(1:nrow(data.tr.mlp),breaks=10,labels=FALSE)
er.perform.tr.mlp <- c()
for (i in 1:10){
  data.tr.part <- data.tr.mlp[which(fold==i),]
  data.tr.main <- data.tr.mlp[-which(fold==i),]
  data.tr.part.cl <- data.tr.part[,1]
  data.tr.part.da <- as.matrix(data.tr.part[,ini.feat.mlp])
  data.tr.main.cl <- data.tr.main[,1]
  data.tr.main.da <- as.matrix(data.tr.main[,ini.feat.mlp])
  score <- numeric(5)
  count <- 0
  while(count < 5) {
    jj1 <- nnet(as.factor(data.tr.main.cl) ~ ., data = data.tr.main.da, size = 8,
                decay = 0.01, maxit = 1000)
    evs <- eigen(nnetHess(jj1, data.tr.main.da, data.tr.main.cl),
                 T)$values
    if(min(evs) > 1e-005) {
      count <- count + 1
      preds <- predict(jj1, data.tr.part.da)
      score[count] <- sum((preds > 0.5) != data.tr.part.cl)/nrow(data.tr.part.da)
    }
    wm <- which.min(score)
    er.perform.tr.mlp[i] <- score[wm]
  }}
er.perform.tr.mlp <- mean(er.perform.tr.mlp)
dev.off()
#Question 4
#Two best performing classifiers
#1. QDA
#Assessment of performance of QDA based on test set 
#with deselected features {22  6  9 21 15 28  8}
er.perform.te.qda <- c()
for(i in 1:10){
data.tr.qda <- data.tr[,2:29]
data.te.qda <- data.te[,2:29]
qda.jj1 <- qda(data.tr.qda[,-(ini.feat.qda-1)],data.tr[,1])
pred.qda.te <- predict(qda.jj1,data.te.qda[,-(ini.feat.qda-1)])$class
er.perform.te.qda[i] <- sum(pred.qda.te != data.te[,1])/(nrow(data.te))}
er.perform.te.qda <- mean(er.perform.te.qda)
#2. MLP
#Assessment of performance of MLP analysis based on testing set
#with selected features { 2, 13, 10, 18, 14, 17, 5, 24, 4 and 6} and parameters H=8,decay=0.01
ini.feat.mlp <- c(2,13,10,18,14,17,5,24,4,6)
er.perform.te.mlp <- c()
  score <- numeric(5)
  count <- 0
  while(count < 5) {
    jj1 <- nnet(as.factor(data.tr.mlp[,1]) ~ ., data = data.tr.mlp[,ini.feat.mlp], size = 8,
                decay = 0.01, maxit = 1000)
    evs <- eigen(nnetHess(jj1, data.tr.mlp[,ini.feat.mlp], data.tr.mlp[,1]),
                 T)$values
    if(min(evs) > 1e-005) {
      count <- count + 1
      pred.mlp.te <- predict(jj1, data.te[,ini.feat.mlp])
      score[count] <- sum((pred.mlp.te > 0.5) != data.te[,1])/nrow(data.te)
    }
    wm <- which.min(score)
    er.perform.te.mlp<- score[wm]
  }
er.perform.te.mlp
dev.off()
#McNemar's Test
mcnemar<-cbind(as.numeric(pred.mlp.te > 0.5),as.numeric(pred.qda.te)-1,data.te[,1])
#number of samples misclassified by MLP but not QDA
n01 <- sum((mcnemar[,1] != mcnemar[,3]) & (mcnemar[,2] == mcnemar[,3]))
#number of samples misclassified by QDA but not MLP
n10 <- sum((mcnemar[,1] == mcnemar[,3]) & (mcnemar[,2] != mcnemar[,3]))
#number of samples misclassified by both MLP and QDA
n00 <- sum((mcnemar[,1] != mcnemar[,3]) & (mcnemar[,2] != mcnemar[,3]))
#number of sample misclassified by neither MLP nor QDA
n11 <- sum((mcnemar[,1] == mcnemar[,3]) & (mcnemar[,2] == mcnemar[,3]))
#Under the null hypothesis, that the classifiers have the same error,the test statistic is
z <- (abs(n01 - n10) - 1)/(sqrt(n10 + n01))
#To conduct a test at the 1% significance level, we reject the null hypothesis if |z| > 2.58
