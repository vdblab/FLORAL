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

binsimuar1 <- function(mu,gamma){
  
  # simulate correlated binary variables using method by Qaqish (2003) with AR(1) corstr
  
  len <- length(mu)
  y <- rep(0,len)
  y[1] <- rbinom(1,1,mu[1])
  for (i in 2:len){
    lambda <- mu[i] + gamma*(y[i-1] - mu[i-1])*sqrt((mu[i]*(1-mu[i]))/(mu[i-1]*(1-mu[i-1])))
    if (lambda < 0) lambda = 0.0001
    if (lambda > 1) lambda = 0.9999
    y[i] <- rbinom(1,1,lambda)
  }
  y
}

binsimuexch <- function(mu,gamma){
  
  # simulate correlated binary variables using method by Qaqish (2003) with exchangeable corstr
  
  len <- length(mu)
  y <- rep(0,len)
  y[1] <- rbinom(1,1,mu[1])
  
  for (i in 2:len){
    
    lambda <- mu[i]
    
    for (j in 1:(i-1)){
      
      lambda <- lambda + (gamma/(1+(i-2)*gamma)) * sqrt((mu[i]*(1-mu[i]))/(mu[j]*(1-mu[j]))) * (y[j]-mu[j])
      
    }
    
    if (lambda < 0) lambda = 0.0001
    if (lambda > 1) lambda = 0.9999
    y[i] <- rbinom(1,1,lambda)
  }
  y
}
