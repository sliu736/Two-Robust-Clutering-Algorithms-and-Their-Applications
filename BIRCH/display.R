rm(list=ls())
library(factoextra)
library(cluster)
library(tidyverse)
source('BIRCH_v1.R')
#kmeans exp
ggplot(iris,aes(iris$Sepal.Length,iris$Sepal.Width))+ geom_point(size=2.5, aes(colour=factor(iris$Species)))+theme(legend.position='right')+
  guides(colour=guide_legend(title='cluster'))+ labs(title = "Iris")
for(i in 1:100){
km<-kmeans(iris[1:4],3)
ratio<-c()
for(i in 1:length(unique(km$cluster))){
  nam<-iris[5][which(km$cluster==i),]%>%table()%>%which.max()%>%names
  ratio<-c(ratio,sum(iris[5][which(km$cluster==i),]!=nam))
}
kr<-sum(ratio)/150
if(kr>0.3) break
}
ggplot(iris,aes(iris$Sepal.Length,iris$Sepal.Width))+ geom_point(size=2.5, aes(colour=factor(km$cluster)))+theme(legend.position='right')+
  guides(colour=guide_legend(title='cluster'))+ labs(title = "k-means2")

#birch+kmeans exp
for(j in 1:100){
ls<-data.frame(x1=NA,x2=NA,x3=NA,x4=NA)
BIRCH(iris[1:4],3,3,0.01,dist='intercluster')

i=1
ci<-c()
index<-c()
for(cf in CFtree[[1]]){
  ls[i,]<-cf$LS/cf$N
  ci<-c(ci,cf$CI)
  index<-c(index,rep(cf$index,length(cf$CI)))
  i=i+1
}
km1<-kmeans(ls,3)
clust<-as_tibble(data.frame(clu=km1$cluster,index=1:length(km1$cluster)))
dat<-as_tibble(data.frame(ci=ci,index=index))
join<-dat%>%left_join(clust,by='index')%>%arrange(ci)
ratio<-c()
for(i in 1:length(unique(km1$cluster))){
  nam<-iris[5][which(km1$cluster==i),]%>%table()%>%which.max()%>%names
  ratio<-c(ratio,sum(iris[5][which(km1$cluster==i),]!=nam))
}
ratio
br<-sum(ratio)/150
if(br>0.3) break
}

ggplot(iris,aes(iris$Sepal.Length,iris$Sepal.Width))+ geom_point(size=2.5, aes(colour=factor(join$clu)))+theme(legend.position='right')+
  guides(colour=guide_legend(title='cluster'))+ labs(title = "BIRCH")

#============================
library(ggplot2)
library(plyr)
library(gridExtra)
set.seed(665544)
x1 <- seq(0,pi,length.out=100)
y1 <- sin(x1) + 0.1*rnorm(100)
x2 <- 1.5+ seq(0,pi,length.out=100)
y2 <- cos(x2) + 0.1*rnorm(100)
data <- data.frame(c(x1,x2),c(y1,y2))
names(data) <- c('x','y')
qplot(data$x, data$y)

BIRCH(data,8,8,1,dist='intercluster')
ls<-data.frame(x=NA,y=NA)
i=1
ci<-c()
index<-c()
for(cf in CFtree[[1]]){
  ls[i,]<-cf$LS/cf$N
  ci<-c(ci,cf$CI)
  index<-c(index,rep(cf$index,length(cf$CI)))
  i=i+1
}
km1<-kmeans(ls,2)
clust<-as_tibble(data.frame(clu=km1$cluster,index=1:length(km1$cluster)))
dat<-as_tibble(data.frame(ci=ci,index=index))
join<-dat%>%left_join(clust,by='index')%>%arrange(ci)

ggplot(data,aes(x,y))+ geom_point(size=2.5, aes(colour=factor(join$clu)))+theme(legend.position=c(0.85,0.85))+
  guides(colour=guide_legend(title='cluster'))+ labs(title = "B=L=8,T=1")

#=============record===============
record<-data.frame(depth=NA,leaf=NA)
for(i in 1:4){
  BIRCH(data,2*i,2*i,0.01,dist='intercluster')
  depth<-length(CFtree)
  leaf<-length(CFtree[[1]])
  record[i,]$depth=depth
  record[i,]$leaf=leaf
}
write.table(record,'record.txt')
record1<-data.frame(depth=NA,leaf=NA)
se<-c(0.001,0.01,0.1,1)

for(i in 1:4){
  BIRCH(data,8,8,se[i],dist='intercluster')
  depth<-length(CFtree)
  leaf<-length(CFtree[[1]])
  record1[i,]$depth=depth
  record1[i,]$leaf=leaf
}
write.table(record1,'record1.txt')
#===================
set.seed(665544)
n <- 600
data <- data.frame(cbind(runif(10, 0, 10)+rnorm(n, sd=0.2), runif(10, 0, 10)+rnorm(n,sd=0.2)))
names(data) <- c('x','y')
qplot(data$x, data$y)

  BIRCH(data,5,5,0.1,dist='intercluster')
  ls<-data.frame(x=NA,y=NA)
  i=1
  ci<-c()
  index<-c()
  for(cf in CFtree[[1]]){
    ls[i,]<-cf$LS/cf$N
    ci<-c(ci,cf$CI)
    index<-c(index,rep(cf$index,length(cf$CI)))
    i=i+1
  }
  km1<-kmeans(ls,10)
  clust<-as_tibble(data.frame(clu=km1$cluster,index=1:length(km1$cluster)))
  dat<-as_tibble(data.frame(ci=ci,index=index))
  join<-dat%>%left_join(clust,by='index')%>%arrange(ci)

ggplot(data,aes(x,y))+ geom_point(size=2.5, aes(colour=factor(join$clu)))+theme(legend.position='right')+
  guides(colour=guide_legend(title='cluster'))+ labs(title = "B=L=5,T=0.1")

km<-kmeans(data,10)

ggplot(data,aes(x,y))+ geom_point(size=2.5, aes(colour=factor(km$cluster)))+theme(legend.position='right')+
  guides(colour=guide_legend(title='cluster'))+ labs(title = "k-means")
#=====================
data<-read.table("PigNosePoints.txt")
names(data) <- c('x','y')
qplot(data$x, data$y)

BIRCH(data,8,8,0.1,dist='intercluster')
ls<-data.frame(x=NA,y=NA)
i=1
ci<-c()
index<-c()
for(cf in CFtree[[1]]){
  ls[i,]<-cf$LS/cf$N
  ci<-c(ci,cf$CI)
  index<-c(index,rep(cf$index,length(cf$CI)))
  i=i+1
}
km1<-kmeans(ls,3)
clust<-as_tibble(data.frame(clu=km1$cluster,index=1:length(km1$cluster)))
dat<-as_tibble(data.frame(ci=ci,index=index))
join<-dat%>%left_join(clust,by='index')%>%arrange(ci)

ggplot(data,aes(x,y))+ geom_point(size=2.5, aes(colour=factor(join$clu)))+theme(legend.position='right')+
  guides(colour=guide_legend(title='cluster'))+ labs(title = "B=L=8,T=0.1")

km<-kmeans(data,3)

ggplot(data,aes(x,y))+ geom_point(size=2.5, aes(colour=factor(km$cluster)))+theme(legend.position='right')+
  guides(colour=guide_legend(title='cluster'))+ labs(title = "k-means")

