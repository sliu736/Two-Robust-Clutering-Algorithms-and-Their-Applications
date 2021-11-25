library(dplyr)
BIRCH<-function(dataset, B, L, Thre, dist){
  dataset<-as.data.frame(dataset)
  CFtree<<-list()
  for(i in 1:nrow(dataset)){
    Generate_CFtree(newpoint=dataset[i,]%>%as.numeric(), B, L, Thre, i, dist)
    print(i)
  }
}
#=============MAIN FUNCTION==============
Generate_CFtree<-function(newpoint, B=3, L=3, Thre=0.3, i, dist){
  depth<-length(CFtree)
  if(depth==0){ #when tree is empty, build a 2 layer CF tree directly
    Generate_nullnode(newpoint, 0, 1, 1, 1)
  }
  else{
    index<-Find_nearCF(newpoint, depth, CI=CFtree[[depth]][[1]]$CI)
    depth<-2
    cf_leaf<-CFtree[[1]][[index]]
    PN<-cf_leaf$PI
    cf0<-CFtree[[2]][[PN]]
    #whether within the radius of T
    if(distance(newpoint, CFtree[[1]][[index]]$LS/CFtree[[1]][[index]]$N)<=Thre){
      #update the new added leaf node
      CFtree[[1]][[index]]<<-list(N=cf_leaf$N+1, LS=cf_leaf$LS+newpoint, SS=cf_leaf$SS+sum(newpoint^2),
                                PI=cf_leaf$PI, CI=append(cf_leaf$CI,i), index=cf_leaf$index)
    }
    else{ #without the radius, create the new leaf node CF and update root node
      cf_new<-Generate_nullnode(newpoint, 1, i, PN, NULL)
      #update root node, triple group already updated in 'Find_nearCF'
      CFtree[[2]][[PN]]$CI<<-append(cf0$CI, cf_new$index)
      #whether need to split node
      Split_node(B, L, 1, CI=CFtree[[2]][[PN]]$CI, i, dist)
      }
    }
  }
#=======CREATE NEW NODE======
Generate_nullnode<-function(newpoint, depth=1, i, PI=0, CI){
  if(depth==0){ #empty tree
    temp<-list(N=1, LS=newpoint, SS=sum(newpoint^2), PI=1, CI=1, index=1)
    CFtree[[paste0('height_',depth+1)]]<<-list(temp)
    temp<-list(N=1, LS=newpoint, SS=sum(newpoint^2), PI=0, CI=1, index=1)
    CFtree[[paste0('height_',depth+2)]]<<-list(temp)
  }
  else if(depth==1){ #if leaf node, requires PI
    temp<-list(N=1, LS=newpoint, SS=sum(newpoint^2), PI=PI, CI=i, index=length(CFtree[[1]])+1)
    CFtree[[paste0('height_',depth)]][[temp$index]]<<-temp
  }
  else{ #not leaf node, requires PI, CI
    if(is.null(CFtree[[paste0('height_',depth)]])==T){#create root node
      temp<-list(N=i, LS=newpoint, SS=Extract(depth-1, 1)$SS%>%unlist%>%as.numeric(), PI=PI, CI=CI, index=1)
      CFtree[[paste0('height_',depth)]]<<-list(temp)
    }
    else{#create split node, N is the sample size of every CI
      #a<-sapply(Extract(depth-1, CI)$CI,'[',i = 1:max(sapply(Extract(depth-1, CI)$CI, length)))%>%t()%>%as.data.frame()%>%is.na()
      temp<-list(N=Extract(depth-1, CI)$N%>%as.data.frame()%>%sum(), LS=newpoint, SS=sapply(Extract(depth-1, CI)$SS,'[',i = 1)%>%unlist%>%sum(), PI=PI, CI=CI, index=length(CFtree[[depth]])+1)
      CFtree[[paste0('height_',depth)]][[temp$index]]<<-temp
    }
  }
}
#============DISTANCE FUNCTION==========
distance=function(point_1, point_2){
  return(sqrt(sum((point_1-point_2)^2)))
}
distance_CF = function(CF_1, CF_2, dist='intercluster' ){
  if (dist == 'Euclidean'){
    temp = distance(CF_1$LS/CF_1$N, CF_2$LS/CF_2$N)
  }
  else if (dist == 'Manhattan'){
    temp = abs(CF_1$LS/CF_1$N-CF_2$LS/CF_2$N)
  }
  else if (dist == 'intercluster'){
    temp = sqrt(CF_1$SS/CF_1$N-2*sum(CF_1$LS*CF_2$LS)/(CF_1$N*CF_2$N)+CF_2$SS/CF_2$N)
  }
  else{
    stop('The distance method doesn\'t exist.')
  }
  return(temp)
}
#=======AMONG NODE where DEPTH=depth-1, FIND NEAREST NODE to 'newpoint'===
#=======UPDATE ROOT NODE AFTER INSERT=====================================
Find_nearCF<-function(newpoint, depth, CI){
  PN<-CFtree[[depth-1]][[CI[1]]]$PI
  cf0<-CFtree[[depth]][[PN]]
  CFtree[[depth]][[PN]]<<-list(N=cf0$N+1, LS=cf0$LS+newpoint, 
                               SS=cf0$SS+sum(newpoint^2),
                               PI=cf0$PI, CI=cf0$CI, index=cf0$index)
  mindist<-Inf
  for(cf in CFtree[[depth-1]][CI]){
    temp<-distance(newpoint, cf$LS/cf$N)
    if(temp<mindist){
      mindist<-temp
      near<-cf$index
      CN<-cf$CI
    }
  }
  if(depth==2)
    return(near)
  else
    return(Find_nearCF(newpoint, depth-1, CN))
}
#===========SPLIT NODE===========
Split_node<-function(B, L, depth, CI, i, dist){
  BL<-ifelse(depth==1, L, B)
  if(length(CI)>BL){
    PPI<-Find_farCF(depth, CI, i, dist)
    depth<-depth+1
    if(depth<length(CFtree))
      return(Split_node(B, L, depth, CI=CFtree[[depth+1]][[PPI]]$CI, i, dist))
  }
  return(1)
}

#===========FIND FARTHEST 2 CF============
Find_farCF<-function(depth, CI=1:length(CFtree[[1]]), i, dist){
  #compute distance of every 2 CF
  comb<-combn(CI,2)
  maxdist<--Inf
  for(k in 1:ncol(comb)){
    temp<-distance_CF(CFtree[[depth]][[comb[1,k]]], CFtree[[depth]][[comb[2,k]]], dist)
    if(temp>maxdist){
      maxdist<-temp
      far<-comb[,k]
    }
  }
  #find where CF belong except the farthest 2
  which_cf<-c()
  for(cf in CFtree[[depth]][CI[!CI%in%far]]){
    temp<-c(distance_CF(cf, CFtree[[depth]][[far[1]]], dist), distance_CF(cf, CFtree[[depth]][[far[2]]], dist))%>%which.min()
    which_cf<-c(which_cf,temp)
  }
  #judge whether need to increase the layer of the tree, build 2 root nodes and reallocate the child nodes
  tag<-0
  if(length(CFtree[[depth+1]])==1){#if the number of nodes that needs to split is 1
                                  #then need to increase the layer
    Generate_nullnode(CFtree[[depth+1]][[1]]$LS, depth+2, i, PI=0, CI=1:2)
    tag<-1
  }
  PN<-CFtree[[depth]][[CI[1]]]$PI
  CN<-c(far[1],CI[!CI%in%far][which(which_cf==1)])
  PPI<-ifelse(tag==1, 1, CFtree[[depth+1]][[PN]]$PI)
  #build the first node
  cf_new<-Generate_nullnode(newpoint=apply(Extract(depth, CN)[,2]%>%as.data.frame,1,sum), depth+1, i, PI=PPI, CI=CN)
  CN1<-CI[!CI%in%CN]
  ls<-apply(Extract(depth, CN1)[,2]%>%as.data.frame,1,sum)
  #build the second node
  CFtree[[depth+1]][[PN]]<<-list(N=Extract(depth, CN1)$N%>%as.data.frame()%>%sum, LS=ls, SS=Extract(depth, CN1)[,3]%>%as.data.frame%>%sum, PI=PPI, CI=CN1, index=PN)
  #if tree layer doesn't increase, then need to update the CI of the root node
  if(tag==0)
    CFtree[[depth+2]][[PPI]]$CI<<-append(CFtree[[depth+2]][[PPI]]$CI, cf_new$index)
  #update the parent node of the child node
  for(j in 1:length(CI)){
    CFtree[[depth]][[CI[j]]]$PI<<-ifelse(CI[j]%in%CN, cf_new$index, PN)
  }
  return(PPI)
}

Extract<-function(depth, CI=1:length(CFtree[[depth]])){
  CFdf<-as.data.frame(t(sapply(CFtree[[depth]][CI], "[", i = 1:max(sapply(CFtree[[1]], length)))))
  return(CFdf)
}

