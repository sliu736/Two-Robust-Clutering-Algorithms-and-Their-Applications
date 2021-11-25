rm(list=ls())
library(dplyr)
#===========phase1 SAMPLE============
sampleROCK<-function(dataset,f,delta,umin){
  N<-nrow(dataset)
  #minimum sample size
  s<-min(ceiling(f*N+N/umin*log(1/delta)+N/umin*sqrt((log(1/delta))^2+2*f*umin*log(1/delta))),N)
  #sampling the serial number
  sample_id<<-sample(1:N,s,replace = F)
  left_id<<-setdiff(1:N,sample_id)
  dat<-dataset[sample_id,]
  #if no difference in sample size before and after sampling，
  #then return original sample(or not consistent on labeling)
  if(nrow(dat)==nrow(dataset)){
    warning('equal sample size as original data, return with original data')
    return(dataset)
  }
  else return(dat)
}
#=============phase2 ROCK CLUSTERING=============
avesamp<-function(samp){
  n <- nrow(samp)
  comb<-combn(n,2,simplify=T)
  averdis <- 0
  for(i in 1:ncol(comb)){
    averdis <- averdis + sqrt(sum(samp[comb[1,i],]-samp[comb[2,i],])^2)
  }
  averdis <<- averdis/(ncol(comb))
}


Similarity<-function(p1, p2, funct='Jaccard'){
  if(funct=='Jaccard'){
    s1<-sum(p1|p2)
    if(s1==0) s1<-1
    return(sum(p1&p2)/s1)
  }
  if(funct=='cosine'){
    return(sum(p1*p2)/sqrt(sum(p1^2)*sum(p2^2)))
  }
  if(funct=='L2'){
    return(sqrt(sum((p1-p2)^2))/averdis)
  }
}

Goodness<-function(theta, n1, n2, link){
  return(link/((n1+n2)^(1+2*f(theta))-n1^(1+2*f(theta))-n2^(1+2*f(theta))))
}

f<-function(theta){
  return((1-theta)/(1+theta))
}

Compute_link<-function(dataset, theta, funct, type){
  n<-nrow(dataset)
  link<-nbr<-matrix(rep(0,n*n), nrow=n)
  for(i in 1:(n-1))
    for(j in (i+1):n)
      nbr[j,i]<-nbr[i,j]<-Similarity(dataset[i,],dataset[j,],funct)
  nbr[type==1]<-nbr<=theta
  nbr[type!=1]<-nbr>=theta
  #nbr%*%nbr, diagnal=0
  for(i in 1:n){
    if(which(nbr[i,]==T)%>%length>1){
    comb<-combn(which(nbr[i,]==T),2)
      for(j in 1:ncol(comb)){
        link[comb[1,j],comb[2,j]]<-link[comb[1,j],comb[2,j]]+1
        link[comb[2,j],comb[1,j]]<-link[comb[2,j],comb[1,j]]+1
      }
    }
  }
  return(link)
}

#============MAIN FUNCTION=============
Build_localheap<-function(dataset, link, theta){
  localheap<-list()
  for(i in 1:nrow(dataset)){
    goodness<-data.frame(goodness=NA, j=NA)
    count<-1
    for(j in which(link[i,]!=0)){
      goodness[count,1]<-Goodness(theta, 1, 1, link[i,j])
      goodness[count,2]<-j
      count<-count+1
    }
    localheap[[i]]<-list(gninf=goodness, index=i)
  }
  return(localheap)
}

Build_globalheap<-function(localheap){
  globalheap<-data.frame(max_goodness=NA, max_i=NA, max_j=NA)
  count<-1
  for(i in 1:length(localheap)){
    c_i<-localheap[[i]]
    if(is.na(c_i$gninf[1,1])) next
    else{
      globalheap[count,]<-c(c_i$gninf[,1]%>%max,
                            i,
                            c_i$gninf[,2][c_i$gninf[,1]%>%which.max])
      count<-count+1
    }
  }
  return(globalheap)
}


ROCK<-function(dataset, theta, k, funct, type){
  if(funct=='L2') avesamp(dataset)
  link<-Compute_link(dataset, theta,funct,type)
  localheap<-Build_localheap(dataset, link, theta)
  globalheap<-Build_globalheap(localheap)
  while(nrow(globalheap)>k){
    u<-globalheap[which.max(globalheap$max_goodness),3]
    v<-globalheap[which.max(globalheap$max_goodness),2]
    globalheap<-globalheap[-which(globalheap[,2]==u|globalheap[,2]==v),]
    w<-length(localheap)+1
    localheap[[w]]<-list(gninf=data.frame(goodness=NA, j=NA), index=c(localheap[[u]]$index, localheap[[v]]$index)%>%sort)
    inter<-c(localheap[[u]]$gninf[,2],localheap[[v]]$gninf[,2])%>%base::unique()
    inter<-inter[-which(inter==v|inter==u)]
    link<-link%>%cbind(rep(0,nrow(link)))%>%rbind(rep(0,ncol(link)+1))
    if(length(inter)!=0){
      for(x in inter){
        link[w,x]<-link[x,w]<-link[x,u]+link[x,v]#no need to judge x,u(write into symmetric matrix)
        #delete u,v in 'localheap[[x]]'
        localheap[[x]]$gninf<-localheap[[x]]$gninf[-which(localheap[[x]]$gninf[,2]==u|localheap[[x]]$gninf[,2]==v),]
        gn_wx<-Goodness(theta,length(localheap[[w]]$index),length(localheap[[x]]$index),link[w,x])
        #update 'localheap[[x]]' and add 'localheap[[w]]'
        localheap[[x]]$gninf<-rbind(localheap[[x]]$gninf,c(gn_wx,w))
        localheap[[w]]$gninf<-rbind(localheap[[w]]$gninf,c(gn_wx,x))
        #update the max goodness in 'globalheap' (might change)
        c_x<-localheap[[x]]
        globalheap[which(globalheap$max_i==x),]<-c(c_x$gninf[,1]%>%max, x,
                            c_x$gninf[,2][c_x$gninf[,1]%>%which.max])
      }
    localheap[[w]]$gninf<-localheap[[w]]$gninf[-1,]
    }
    globalheap<-rbind(globalheap, c(localheap[[w]]$gninf[,1]%>%max, w,
                                    localheap[[w]]$gninf[,2][localheap[[w]]$gninf[,1]%>%which.max]))
    localheap[[u]]<-w
    localheap[[v]]<-w
    print(nrow(globalheap))
  }
  globalheap<<-globalheap;localheap<<-localheap
  link<<-link
  return(1)
}

#===============phase3 LABEL DATA============
label<-function(dataset,u,newpoint,theta,k,funct,type){
  criteria<-c()
  for(i in 1:length(u)){
    d<-dataset[clu[which(clu[,2]==i),][,1],]
    sim<-c()
    for(j in 1:nrow(d)) sim<-c(sim,Similarity(newpoint,d[j,],funct = funct))
    ni<-ifelse(type==1,sum((sim<=theta)+0),sum((sim>=theta)+0))
    criteria<-c(criteria,ni/(1+nrow(d))^f(theta))
  }
  clu[nrow(clu)+1,]<<-c(index=left_id[k],cluster=which.max(criteria))
}
#main function

LableonDisk<-function(dataset,samp,theta,funct,type){
  u<-globalheap$max_i
  sequ<-seq(1,length(u),1)
  index_sample<-c()
  cl<-c()
  for(i in 1:length(u)){
    l<-localheap[[u[i]]]$index
    index_sample<-c(index_sample,l)
    cl<-c(cl,rep(sequ[i],length(l)))
  }
  #resort the serial no and the belonging cluster after sampling
  clu<<-data.frame(index=index_sample,cluster=cl)  
  if(nrow(dataset)==nrow(samp)) return(clu)
  clu$index<<-sample_id[clu$index]
  if(length(left_id)!=1)  for(k in 1:length(left_id)) {
    label(dataset,u,dataset[left_id,][k,],theta,k,funct,type)
    print(k)
  }
  else  label(dataset,u,dataset[left_id,],theta,1,funct)
  return(clu)
}
