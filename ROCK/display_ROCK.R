rm(list=ls())
source("ROCK.r")
library(ggplot2)
#=========实例数据集1========
set.seed(665544)
x1 <- seq(0,pi,length.out=100)
y1 <- sin(x1) + 0.1*rnorm(100)
x2 <- 1.5+ seq(0,pi,length.out=100)
y2 <- cos(x2) + 0.1*rnorm(100)
data <- data.frame(c(x1,x2),c(y1,y2))
names(data) <- c('x','y')

samp<-sampleROCK(data,0.3,0.1,nrow(data)/2)
ROCK(samp,1,2,funct='L2',type=1)
result<-LableonDisk(data,samp,theta=1,funct='L2',type=1)%>%arrange(index)

ggplot(data,aes(x,y))+ geom_point(size=2.5, aes(colour=factor(result$cluster)))+theme(legend.position=c(0.85,0.85))+
  guides(colour=guide_legend(title='cluster'))+ labs(title = "theta=0.8,N=89")

s<-c()
s1<-c()
s2<-c()
for(i in seq(0.01,1,0.01)){
  s<-c(s,sampleROCK(data,i,0.1,nrow(data)/2)%>%nrow())
  s1<-c(s1,sampleROCK(data,0.5,i,nrow(data)/2)%>%nrow())
  s2<-c(s2,sampleROCK(data,0.5,0.1,nrow(data)/(ceiling(i*10)))%>%nrow())
}
data<-c(s,s1,s2)
group<-rep(c('f','delta','umin'),each=100)
a<-data.frame(s=data,group=group,x=seq(0.01,1,0.01))

ggplot(a, aes(x=x, y=s, group,color=group)) + geom_line(size=1)
#=========模拟数据集2========
set.seed(665544)
n <- 600
data <- data.frame(cbind(runif(10, 0, 10)+rnorm(n, sd=0.2), runif(10, 0, 10)+rnorm(n,sd=0.2)))
names(data) <- c('x','y')

samp<-sampleROCK(data,0.1,0.01,nrow(data)/10)
ROCK(samp,0.2,10,funct='L2',type=1)
result<-LableonDisk(data,samp,theta=0.2,funct='L2',type=1)%>%arrange(index)
ggplot(data,aes(x,y))+ geom_point(size=2.5, aes(colour=factor(result$cluster)))+theme(legend.position='right')+
  guides(colour=guide_legend(title='cluster'))+ labs(title = "theta=0.2,N=194")
#=========模拟数据集3========
data<-read.table("PigNosePoints.txt")
names(data) <- c('x','y')
qplot(data$x, data$y)

samp<-sampleROCK(data,0.2,0.01,nrow(data)/3)
ROCK(samp,0.3,3,funct='L2',type=1)
result<-LableonDisk(data,samp,theta=0.3,funct='L2',type=1)%>%arrange(index)
ggplot(data,aes(x,y))+ geom_point(size=2.5, aes(colour=factor(result$cluster)))+theme(legend.position='right')+
  guides(colour=guide_legend(title='cluster'))+ labs(title = "theta=0.3,N=194")
#=========模拟数据集4========
data1 = matrix(c(1,1,1,1,1,0,0,0,0,0,
                   1,1,1,0,1,0,0,0,0,0,
                   0,1,0,1,1,0,0,0,0,0,
                   0,1,1,0,1,0,0,0,0,0,
                   0,0,1,1,1,0,0,0,0,0,
                   0,0,0,0,0,1,0,1,0,1,
                   0,1,0,1,0,0,1,0,0,0,
                   0,0,0,0,0,1,1,1,1,1,
                   0,0,0,0,0,1,1,0,1,1,
                   0,0,0,0,0,0,1,1,1,1,
                   0,0,0,0,0,1,0,0,0,0,
                   0,0,0,0,0,1,0,0,1,1,
                   1,0,0,0,0,1,0,0,1,0,
                   0,0,0,0,0,0,0,1,0,1),ncol=10,byrow=T)
#
ROCK(data1,0.01,2,funct='Jaccard',type=2)
result<-LableonDisk(data1,data1,theta=0.01,funct='Jaccard',type=2)%>%arrange(index)

#=========实例数据集1========
data<-dataset2
samp<-sampleROCK(data,0.1,0.1,nrow(data)/3)
ROCK(samp,0.01,3,funct='Jaccard',type=2)
result<-LableonDisk(data,samp,theta=0.01,funct='Jaccard',type=2)%>%arrange(index)
miss<-setdiff(1:nrow(data),result$index)
result1<-result
ratio<-c()
for(i in 1:length(unique(result$cluster))){
  nam<-test$Class[which(result$cluster==i)]%>%table()%>%which.max()%>%names
  ratio<-c(ratio,sum(test$Class[which(result$cluster==i)]!=nam))
  result$cluster[which(result$cluster==i)]<-nam
}
sum(ratio)/(150-length(miss))
confusionMatrix(table(result$cluster, test$Class[-miss]%>%as.character()))
#=========实例数据集2========

data2<-iris[1:4]

samp<-sampleROCK(data2,0.5,0.025,nrow(data)/3)
ROCK(data2,0.4,3,funct='L2',type=1)
result<-LableonDisk(data2,data2,theta=0.4,funct='L2',type=1)%>%arrange(index)


ratio<-c()
for(i in 1:length(unique(result$cluster))){
  nam<-iris[5][which(result$cluster==i),]%>%table()%>%which.max()%>%names
  ratio<-c(ratio,sum(iris[5][which(result$cluster==i),]!=nam))
}
ratio
sum(ratio)/150


ggplot(data2,aes(iris$Sepal.Length,iris$Sepal.Width))+ geom_point(size=2.5, aes(colour=factor(result$cluster)))+theme(legend.position='right')+
  guides(colour=guide_legend(title='cluster'))+ labs(title = "ROCK,N=85,theta=0.4")




