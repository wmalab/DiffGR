
#' DiffGR: Differential genomic region detection function
#'
#'
#' @dat1,@dat2 numeric. N*N HiC contact maps which have been preprocessed (2D mean filter smoothing and KR normalization)
#' @tad1 @tad2 numeric. A vector of TAD boundaries of dat1 and dat2 (calling by HiCseg)
#' @res numeric. The resolution of HiC contact maps, eg:100kb will input 100,000
#' @N.perm numeric. The number of iterations in permutation test
#' @speedup.option logical. (True/FALSE) Calculation with or without speed-up algorithm
#' @alpha logical. Significant level of differential region testing 
#'
#'
#' @return a list that contains the tad result and genomic region result
#' 
#' The tad result table contains the following elements:
#' tad.start: the starting locus of TAD
#' tad.end: the end locus of TAD
#' scc: the SCC value of corresponding domain
#' pvalue: the pvalue of differential testing on corresponding domain
#' pvalue.adj: the adjusted pvalue of differential testing on corresponding domain (adjusted by Benjamin-Hochberg)
#' 
#' The genomic result table contains the following elements:
#' genom.start: the starting locus of genomic region
#' genom.end: the end locus of genomic region
#' condition.type: the type if candidate genomic region belonging to. 1:single-TAD, 2: Hierachical-TAD, 3: Alternating-TAD
#' detect.result: the differential testing result for corresponding genomic region. 1:Differential 0:Non-differential 
#'  


DiffGR<- function(dat1,tad1,dat2,tad2,res,N.perm=2000,speedup.option=TRUE,alpha=0.05){
  library(HiCcompare)
  library(hicrep)
  library(R.utils)
  max.distance <- 10000000/res
  
  # Function of calculating SCC 
  scc <- function(dat1.i,dat2.i){
  dat1.i <- as.matrix(dat1.i)
  colnames(dat1.i) <- 1:ncol(dat1.i)
  rownames(dat1.i) <- 1:ncol(dat1.i)
  dat1.vec <- MatToVec(dat1.i)
  dat2.i <- as.matrix(dat2.i)
  colnames(dat2.i) <- 1:ncol(dat2.i)
  rownames(dat2.i) <- 1:ncol(dat2.i)
  dat2.vec <- MatToVec(dat2.i)
  dat.comb <- matrix(0,ncol=4,nrow=nrow(dat1.vec))
  dat.comb[,1:3] <- dat1.vec
  dat.comb[,4] <- dat2.vec[,3]
  posi.0 <- sort(union(which(dat.comb[,3]+dat.comb[,4]==0),which(dat.comb[,1]==dat.comb[,2])))
  sum.1 <- sum(dat.comb[-posi.0,3])
  sum.2 <- sum(dat.comb[-posi.0,4])
  if(length(posi.0)!=0){
    dat.comb <- dat.comb[-posi.0,]
  }
  dd <- dim(dat.comb)
  if(dd[1]==0){
    scc.result <- 1
  }else{
    if(sum.1==0 || sum.2==0){
      scc.result <- -1
    }else{
      scc.output <- get.scc(dat.comb,1,max.distance)
      scc.result <- scc.output$scc[1,1]
    }
  }
  return(scc.result)
}
  
  # Function of permutation test
  perm <- function(dat1,dat2,len,N.perm){
  rho.vec <- c()
  n <- nrow(dat1)
  for(l in 1:N.perm){
    dat1.i <- matrix(0,nrow=len,ncol=len)
    dat2.i <- matrix(0,nrow=len,ncol=len)
    for(i in 1:len){
      posi.i <- sample(1:(n-i+1),(len-i+1),replace=FALSE)
      for(j in 1:(len-i+1)){
        dat1.i[j,(j+i-1)] <- dat1[posi.i[j],(posi.i[j]+i-1)]
        dat1.i[(j+i-1),j] <- dat1.i[j,(j+i-1)]
        dat2.i[j,(j+i-1)] <- dat2[posi.i[j],(posi.i[j]+i-1)]
        dat2.i[(j+i-1),j] <- dat2.i[j,(j+i-1)]
      }
    }
    rho.vec[l] <- scc(dat1.i,dat2.i)
  }
  rho.na <- is.na(rho.vec)
  rho.vec <- rho.vec[rho.na==FALSE]
  rho.vec <- sort(rho.vec)
  rho.vec <- rho.vec[1:(alpha*N.perm)]
  return(list(rho.vec=rho.vec))
}

  # Function of TAD boundary adjustment
  boundary.adj <- function(tad1,tad2){
    tad.boundary.m <- sort(union(tad1,tad2))
    tad.boundary.table <- matrix(0,ncol=3,nrow=(length(tad.boundary.m)))
    tad.boundary.table[,1] <- tad.boundary.m
    for(i in 1:length(tad.boundary.m)){
      if(tad.boundary.m[i] %in% tad1) {tad.boundary.table[i,2] <- 1}
      if(tad.boundary.m[i] %in% tad2) {tad.boundary.table[i,3] <- 1}
    }
    tad.merge.po <- which(tad.boundary.table[,2]+tad.boundary.table[,3]==1)
    tad.boundary.d <-
      tad.boundary.table[which(tad.boundary.table[,2]+tad.boundary.table[,3]==1),1]
    diff.tad <- c()
    same.group <- c()
    for(i in 1:(length(tad.boundary.d)-1)){
      diff.tad[i] <- tad.boundary.d[i+1]-tad.boundary.d[i]
      same.group[i] <-
        tad.boundary.table[tad.merge.po[i+1],2]-tad.boundary.table[tad.merge.po[i],2]
    }
    merge.bin <- intersect(which(diff.tad<3),which(same.group!=0))
    if(length(merge.bin)>0){
      for(j in 1:length(merge.bin)){
        if(tad.boundary.table[tad.merge.po[merge.bin[j]],1]!=0){
          tad.boundary.table[tad.merge.po[merge.bin[j]],1] <-
            floor((tad.boundary.table[tad.merge.po[merge.bin[j]],1]+tad.boundary.table[tad.merge.po[merge.bin[j]]+1,1])/2)
          tad.boundary.table[tad.merge.po[merge.bin[j]],2:3] <- rep(1,2)
          tad.boundary.table[tad.merge.po[merge.bin[j]]+1,1] <- 0
        }
      }
      tad.boundary.table <- tad.boundary.table[-which(tad.boundary.table[,1]==0),]
    }
    tad.bound <- matrix(0,ncol=3,nrow=nrow(tad.boundary.table)+1)
    tad.bound[1,2:3] <- rep(1,2)
    tad.bound[2:nrow(tad.bound),] <- tad.boundary.table
    return(tad.bound)
  }
  
  
  
  tad.bound <- boundary.adj(tad1,tad2)
  if(tad.bound[2,1]==1){
    tad.bound <- tad.bound[-2,1]
  }
  bins <- which(tad.bound[,2]+tad.bound[,3]==2)
  tad.interval <- c()
  genomic.interval <- c()
  condition.type <- c()
  for(i in 1:(length(bins)-1)){
    tad.bound.i <- tad.bound[bins[i]:bins[i+1],1]
    genomic.interval <- rbind(genomic.interval,tad.bound[c(bins[i],bins[i+1]),1])
    tad.interval <- rbind(tad.interval,t(combn(tad.bound.i,2)))
    
    if(bins[i+1]-bins[i]==1){
      condition.type[i] <- 1
    }else{
       if(bins[i+1]-bins[i]==2){
      condition.type[i] <- 2
    }else{
      tad.cutpoint <- tad.bound[(bins[i]+1):(bins[i+1]-1),2:3]
      if(sum(tad.cutpoint[,1])!=0 && sum(tad.cutpoint[,2])!=0){
        condition.type[i] <- 3
      }else{
        condition.type[i] <- 2
      }
    }
    }
  }
  
  len.tad.uniq <- sort(unique(tad.interval[,2]-tad.interval[,1]+1))
  

  
  rho.table <- matrix(0,nrow=length(len.tad.uniq),ncol=((alpha*N.perm+1)))
  
  if(speedup.option==FALSE){
   for(l in 1:length(len.tad.uniq)){
    perm.result  <- perm(dat1,dat2,len.tad.uniq[l],N.perm)
    rho.table[l,] <- c(len.tad.uniq[l],perm.result$rho.vec)
   }
  }else{
    
    rho.table[,1] <- len.tad.uniq
    sample.num <- min(round((length(len.tad.uniq)-15)*0.25),40)
    len.chosen <- c(4:15, sort(sample(16:length(len.tad.uniq),sample.num)))
    for(l in 1:length(len.chosen)){
      perm.result  <- perm(dat1,dat2,len.tad.uniq[len.chosen[l]],N.perm)
      rho.table[len.chosen[l],2:ncol(rho.table)] <- perm.result$rho.vec
    }
    for(k in 1:(alpha*N.perm)){
      spline.l <- smooth.spline(len.tad.uniq[len.chosen],rho.table[len.chosen,k+1])
      pred.result <- predict(spline.l,len.tad.uniq[16:length(len.tad.uniq)])
      rho.table[16:length(len.tad.uniq),k+1] <- pred.result$y
    }
  }
  
  tad.result <- matrix(0,ncol=5,nrow=nrow(tad.interval))
  colnames(tad.result) <- c("tad.start","tad.end","scc","pvalue","pvalue.adj")
  tad.result[,1:2] <- tad.interval
  for(j in 1:nrow(tad.result)){
    tad.s.length <- tad.result[j,2]-tad.result[j,1]
    if(tad.s.length<5) {
    tad.result[j,3] <- 0
    tad.result[j,4] <- 0.5
    }else{
      dat1.s <- dat1[(tad.result[j,1]+1):tad.result[j,2],(tad.result[j,1]+1):tad.result[j,2]]
      dat2.s <- dat2[(tad.result[j,1]+1):tad.result[j,2],(tad.result[j,1]+1):tad.result[j,2]]
      tad.result[j,3] <- scc(dat1.s,dat2.s)
    }
      tad.result[j,4] <- sum(rho.table[which(rho.table[,1]==tad.s.length),-1]<tad.result[j,3])/N.perm
  }
     tad.result[,5] <- p.adjust(tad.result[,4],method="BH")
     
     detection.result <- rep(0,nrow(genomic.interval))
     for(i in 1:nrow(genomic.interval)){
       s <- min(which(tad.result[,1]==genomic.interval[i,1]))
       e <- max(which(tad.result[,2]==genomic.interval[i,2]))
       detect.result[i] <- sum(tad.result[s:e,5]<alpha)
       lenn <- result[tt,2]-result[tt,1]+1
       #tt.l <- which(lenn>4)
       if(sum(tad.result[s:e,5]<alpha)>0){
         diff.tad <- which(tad.result[s:e,5]<alpha)
         diff.len.max  <- max(tad.result[(diff.tad+s-1),2]-tad.result[(diff.tad+s-1),1])
         if(diff.len.max>max(4,(genomic.interval[i,2]-genomic.interval[i,1])/3)){
           detection.result[i] <- 1
         }
       }
       }
       
     genomic.result <- as.matrix(cbind(genomic.interval,condition.type,detection.result))
     colnames(genomic.result) <- c("genom.start","genom.end","condition.type","detect.result")
  
     return(list(tad.result=tad.result,genomic.result=genomic.result))
}
  
  
  
