coxsplity=function(y, nfolds){
  N=nrow(y)
  tem=data.frame(y, i=seq(N), foldid=0)
  tem=tem[order(y[, "time"], y[, "status"]), ]
  n1=sum(y[, "status"]);n2=N-n1
  
  tem$foldid[tem[, "status"]==1]=sample(rep(seq(nfolds), length=n1))
  tem$foldid[tem[, "status"]==0]=sample(rep(seq(nfolds), length=n2))
  
  foldid=tem$foldid[order(tem$i)]
  return(foldid)
}

coxsplitss=function(y, id, nfolds){
  full = data.frame(y, foldid=0, id=id)
  tem = full %>% dplyr::group_by(.data$id) %>% dplyr::filter(row_number()==n())
  N=nrow(tem)
  tem$i = seq(N)
  tem=tem[order(tem$stop, tem$status), ]
  n1=sum(y[, "status"]);n2=N-n1
  # tem = as.matrix(tem)
  
  tem$foldid[tem$status==1]=sample(rep(seq(nfolds), length=n1))
  tem$foldid[tem$status==0]=sample(rep(seq(nfolds), length=n2))
  
  temif <- tem %>% select(.data$foldid,.data$id)  #data.frame(tem[,c("foldid","id")])
  full <- full %>% select(.data$start,.data$stop,.data$status,.data$id) %>% left_join(temif,by="id")
  
  foldid <- full$foldid
  
  return(foldid)
}

fgfoldid=function(id, foldid){
  idfoldid <- data.frame(id=unique(id),foldid=foldid)
  yfoldid <- data.frame(id=id)
  mergedid <- yfoldid %>% left_join(idfoldid,by="id")
  
  foldid <- mergedid$foldid
  return(foldid)
}