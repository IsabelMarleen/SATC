# Read in data
#filename <- "/home/genis/projects/satc2/leopards/test_leopardwin100kb.txt"
filename <- "examples/test_leopardwin100kb.txt"
raw <- read.table(filename, header=T, as.is=T, row.names = 1)

# Normalisation -----------------------------------------------------------------------------------------------------

# First normalisation step - normalise for scaffold length
norm1 <- apply(raw, 2, function(x) x/raw[, "length"])
stopifnot(all(norm1[, "length"] ==1))

# Second normalisation step - normalise for library size
norm2 <- t(apply(norm1[,-1], 1, function(x) x/median(x)))


# PCA ---------------------------------------------------------------------------------------------------------------

sexDetermine2 <- function(dat,K=2,weight=TRUE,model="gaussian",lengthWeight=FALSE){
  # Run PCA to determine sex for comparison with current SATC
  model <- char.expand(model, c("gaussian","hclust"))
  
  mat_first <- as.matrix(dat)
  noNArow=apply(mat_first,1,function(x) ifelse(sum(is.na(x))==0,TRUE,FALSE))
  #center
  mat <- mat_first[noNArow,]-rowMeans(mat_first[noNArow,])
  
  #    pca <- prcomp(t(mat),scale=F)
  w <- 1
  if(lengthWeight){
    ## weight by squared chromosome length
    w <- sqrt(dat[[1]][noNArow,"Length"])
  }
  maxRank <- min(dim(mat))
  svd <- svd(t(mat*w))
  SIG <- matrix(0,maxRank,maxRank)
  diag(SIG)<-svd$d
  pca <- svd
  pca$x <- svd$u[,1:maxRank]%*%SIG
  
  if(weight){
    d<- pca$x[,1:K]
  } else{
    d <- svd$u[,1:maxRank]
  }
  if(model=="gaussian"){
    group <- Mclust(d,G=2,modelName="EVV")
    if(is.null(group))
      group <- Mclust(d,G=2)
    g <- group$classification
  } else if(model=="hclust"){
    hh<-hclust(dist(d))
    g <- cutree(hh, k=2) # cut
  }
  
  if(sum(g==1) ==1 |  sum(g==2) ==1 ) stop("Dataset can not be used for sex determination (poor clustering). One of the inferred sex groups only have one member.")
  
  beta <- rowMeans(mat[,g==1])-rowMeans(mat[,g==2])
  
  if(!any( abs(beta) > 0.4 & abs(beta) < 0.6)){
    stop("No good candidates for sex scaffold found based on depth of coverage. Try changing the clustering method, or consider SATC might not work for your data.")
  }
  
  if( mean(beta[ abs(beta) > 0.4 & abs(beta) < 0.6]) > 0 ){
    sex <- c("homomorphic","heteromorphic")[g]
  } else{
    sex <- c("heteromorphic","homomorphic")[g]
    beta <- - beta
  }
  homoMedian <- apply( mat_first [,sex=="homomorphic"],1,median)
  heteroMedian <- apply( mat_first [,sex=="heteromorphic"],1,median) 
  pval <- apply(mat,1,function(x) t.test(x~sex)$p.value) #new
  sexAssoScafs <- pval < 0.05/nrow(mat) #new
  
  X_Z_Scaffold <- beta > 0.4 & beta < 0.6  & homoMedian<1.3 & homoMedian>0.7 & heteroMedian<0.7 & heteroMedian>0.3 & sexAssoScafs
  Abnormal <- sexAssoScafs & !X_Z_Scaffold
  outlierScafs<- rowMeans(mat_first) > 1.3
  autoScafs<- as.logical
  
  list(dat=dat,pca=pca,sex=sex,SexScaffolds=data.frame(Name=rownames(dat),X_Z_Scaffolds=X_Z_Scaffold,Abnormal_sex_linked_Scaffolds=Abnormal,Pval=pval,homoMedian=homoMedian,heteroMedian=heteroMedian,stringsAsFactors = FALSE))
}

sex_test <- sexDetermine2(norm2)


## plot the clustering---------------------------------------------------------------------------

plotGroup(sex_test)

## plot the scaffolds depths stratificed by inferred sex including the (abnormal) scaffolds

plotScafs2 <- function(x,ylim,abnormal=FALSE,main=""){
  par(mar=c(4.1,4.1,3.1,2.1))
  #mat <- sapply(x$dat,function(y) y$norm)
  mat <- as.matrix(x$dat)
  rownames(mat) <- as.character(rownames(x$dat))
  sexLinkedScaf <- x$SexScaffolds$Abnormal_sex_linked_Scaffolds|x$SexScaffolds$X_Z_Scaffolds
  XZScaf <- x$SexScaffolds$X_Z_Scaffolds
  
  keep <- as.character(rownames(x$dat[XZScaf]))
  if(abnormal){
    keep<- c(keep, as.character(rownames(x$dat)[sexLinkedScaf & !XZScaf])) # this forces correct order in plot
  }
  mat <- mat[keep,,drop=FALSE]
  #nam <- gsub("NW_0176|NW_0050","",rownames(mat))
  nam <- rownames(mat)
  
  n<- ncol(mat)
  m<-nrow(mat)
  s <- factor(rep(nam,each=n),levels=nam)
  g<-rep(x$sex,m)
  at <- cumsum(rep(c(1.2,1),m))
  if(missing(ylim)){
    ylim <-range(mat)
  }
  b<-boxplot(as.vector(t(mat))~g+s,las=2,col=col12,names=NA,ylab="Normalized Depth",at=at,axes=F,xlab="Sex Associated Scaffolds",lwd=0.5,pch=16,outcol=col12,ylim=ylim,cex.lab=1.2)
  #  abline(v=1:m*2.2+0.6,lty=2)
  axis(2,cex.axis=1.05)
  
  text(1:m*2.2-1.1,y=rep(min(ylim)-diff(ylim)/10,m),nam,xpd=T,srt=45,cex=0.9)
  title(main,adj=0.2,cex.main=1.5)
  gen <- c("Heterogametic ","Homogametic ")
  if(abnormal)
    abline(v=sum(x$SexScaffolds$X_Z_Scaffolds)*2.2+0.6,lty=2)
  #legend("topright",gen,pch=16,cex=1,col=col12,bty="n",horiz=T,)
  legend(2.5*m*0.5,1.1*max(ylim)*1,gen,pch=16,cex=1,col=col12,bty="n",horiz=T,xpd=T)
  abline(h = 0.5, col="#20A387FF", lwd=1, lty=2)
  abline(h = 1, col="#20A387FF", lwd=1, lty=2)
}
plotScafs2(sex_test,abnormal=T)

## plot the sex scaffolds' normalized depth for each individuals ------------------------------------

plotSamples2 <- function(x,ylim=c(0,2),abnormal=TRUE,main=""){
  par(mar=c(8.1,4.1,3.1,2.1))
  # mat <- sapply(x$dat,function(y) y$norm)
  mat <- as.matrix(x$dat)
  scafNames <- as.character(rownames(x$dat))
  indNames<- colnames(x$dat)
  rownames(mat) <- scafNames
  colnames(mat)<- indNames
  sexLinkedScaf <- as.character(scafNames[x$SexScaffolds$Abnormal_sex_linked_Scaffold | x$SexScaffolds$X_Z_Scaffolds])
  XZScaf <- as.character(scafNames[x$SexScaffolds$X_Z_Scaffolds])
  keep <-XZScaf
  if(abnormal){
    keep <-sexLinkedScaf 
  }
  newOrder <- order(x$sex,decreasing=TRUE)
  matOrder<- mat[keep,newOrder]
  namOrder <- gsub("X","",colnames(norm2))
  s <- factor(rep(namOrder,each=nrow(matOrder)),levels=namOrder)
  n<-length(namOrder)
  
  b<-boxplot(as.vector(matOrder)~s,
             col=col12[as.factor(x$sex)[newOrder]]
             ,ylab="Normalized
             Depth",
             axes=F,xlab="",ylim=ylim)
  axis(2,cex.axis=1.05,ylim=ylim)
  text(1:n,y=rep(-0.1,n),namOrder,xpd=T,srt=90,cex=0.7)
  title(main,adj=0.2,cex.main=1.5)
  gen <- c("Heterogametic ","Homogametic ")
  abline(h=0.5,lty=2,lwd=1.5,xlim=c(0,n),col=col12[1])
  abline(h=1,lty=2,lwd=1.5,xlim=c(0,n),col=col12[2])
  abline(v=sum(x$sex=="homomorphic")+0.5,lty=1,lwd=1.5,ylim=c(0,2))
  legend("topright",gen,pch=16,cex=1,col=col12,bty="n",horiz=T,xpd=T)
}

plotSamples2(sex_test)

## visualized the gaussian mixtures
plotUnc(sex_test)

## See the inferred status of each scaffold
head(sex_test$SexScaffolds)

## See the inferred sex of each indiviual
head(cbind(colnames(sex_test$dat),sex_test$sex))
