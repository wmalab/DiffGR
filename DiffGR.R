
#' DiffGR: Differential genomic region detection function
#'
#'
#' @dat1,@dat2 numeric. raw N*N HiC contact maps
#' @res numeric. The resolution of HiC contact maps, eg:100kb will input 100,000
#' @smooth.size numeric. The size varied with different resolutions
#' @N.perm numeric. The number of iterations in permutation test
#' @cutoff.default logical. logical. Whether set the SCC cutoff (meaningful SCC between the two Hi-C datasets that mustbe reached in order to call a differential TAD truly significant) with self-defined value(True) or with automatic computed value (False)
#' @speedup.option logical. (True/FALSE) Calculation with or without speed-up algorithm
#' @alpha numeric. Significant level of differential region testing 
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

DiffGR<- function(dat1,dat2,res,smooth.size,N.perm=2000,cutoff.default=TRUE,speedup.option=TRUE,alpha=0.05){
  library(HiCcompare)
  library(HiCseg)
  library(hicrep)
  library(R.utils)
  
  #Function of KR normalization
  normalize <- function(dat){
    n <- nrow(dat)
    M.A <- as.matrix(dat)
    dat1.sum <-  colSums(M.A)
    posi1 <- which(dat1.sum!=0)
    M.A.s <- KRnorm(M.A)
    M.A.n <- matrix(0,ncol=n,nrow=n)
    M.A.n[posi1,posi1] <- M.A.s
    return(M.A.n)
  }
  
  #Function of generating pseudo replicates
  generate_replicate <- function(dat11){
    colnames(dat11) <- 1:nrow(dat11)
    dat <- as.matrix(full2sparse(dat11))
    dat.g <- dat
    dat.g[,3] <- 0
    n <- nrow(dat)
    prob <- dat[,3]/sum(dat[,3])
    sample.result <- sample(n,sum(dat[,3]),replace=TRUE,prob=prob)
    posi.chose <- sort(unique(sample.result))
    result.sum <- as.vector(table(sample.result))
    dat.g[posi.chose,3] <- result.sum
    pseudo.mat <- sparse2full(dat.g)
    s <- min(as.numeric(rownames(pseudo.mat)))
    e <- max(as.numeric(rownames(pseudo.mat)))
    pseudo.rep <- matrix(0,ncol=nrow(dat11),nrow=nrow(dat11))
    pseudo.rep[s:e,s:e] <- pseudo.mat
    return(pseudo.rep=pseudo.rep)
  }
  
  # Function of calculating SCC 
  scc <- function(dat1,dat2){
    dat1 <- as.matrix(dat1)
    colnames(dat1) <- 1:ncol(dat1)
    rownames(dat1) <- 1:ncol(dat1)
    dat1.vec <- MatToVec(dat1)
    dat2 <- as.matrix(dat2)
    colnames(dat2) <- 1:ncol(dat2)
    rownames(dat2) <- 1:ncol(dat2)
    dat2.vec <- MatToVec(dat2)
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
  
  max.distance <- 10000000/res
  tad <- HiCseg_linkC_R(nrow(dat1),round(nrow(dat1)/3),"P",dat1,"D")
  tad1 <- tad$t_hat[tad$t_hat!=0]
  tad <- HiCseg_linkC_R(nrow(dat2),round(nrow(dat2)/3),"P",dat2,"D")  
  tad2 <- tad$t_hat[tad$t_hat!=0]
  if(min(length(tad1),length(tad1))==1){
    stop("The TADs of HiC maps aren't effectively detected by HiCSeg.\n", call. = FALSE)
  }
  
  
  dat1 <- smoothMat(dat1,smooth.size)
  dat2 <- smoothMat(dat2,smooth.size)
  
  res <- NULL
  tryCatch({
    res <- withTimeout({
      dat1 <- normalize(dat1) 
    }, timeout = 300, onTimeout = "silent")
  })
  if(is.null(res)==TRUE) {
    stop("KR normalization for dat1 did not converge.\n", call. = FALSE)
  }
  
  res <- NULL
  tryCatch({
    res <- withTimeout({
      dat2 <- normalize(dat2) 
    }, timeout = 300, onTimeout = "silent")
  })
  if(is.null(res)==TRUE) {
    stop("KR normalization for dat2 did not converge.\n", call. = FALSE)
  }
  
  
  tad.bound <- boundary.adj(tad1,tad2)
  if(tad.bound[2,1]==1){
    tad.bound <- tad.bound[-2,]
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
  
  len.tad.uniq <- sort(unique(tad.interval[,2]-tad.interval[,1]))
  
  
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
  
  #cutoff
  rho.cutoff <- matrix(0,ncol=2,nrow=length(len.tad.uniq))
  rho.cutoff[,1] <- len.tad.uniq
  if(cutoff.default=TRUE){
    rho.cutoff[,2] <- rep(0.85,length(len.tad.uniq))
  }else{
    
    res <- NULL
    for(kk in 1:5){
    tryCatch({
      res <- withTimeout({
        dat1.rep <- generate_replicate(dat1)
        dat1.rep <- smoothMat(dat1.rep,smooth.size)
        dat1.rep <- normalize(dat1.rep)
      }, timeout = 1800,onTimeout = "silent")
    })
      if(sum(res=="Tac")==1) {break}
    } 
    
    res <- NULL
    for(kk in 1:5){
      tryCatch({
        res <- withTimeout({
          dat2.rep <- generate_replicate(dat2) 
          dat2.rep <- smoothMat(dat2.rep,smooth.size)
          dat2.rep <- normalize(dat2.rep)
        }, timeout = 1800,onTimeout = "silent")
      })
      if(sum(res=="Tac")==1) {break}
    } 
    
    
   
    rho.cutoff1 <- rho.cutoff
    rho.cutoff2 <- rho.cutoff
    
    if(speedup.option==FALSE){
      for(l in 1:length(len.tad.uniq)){
        perm.result  <- perm(dat1,dat1.rep,len.tad.uniq[l],N.perm)
        rho.cutoff1[l,2] <- perm.result$rho.vec[length(perm.result$rho.vec)]
        perm.result  <- perm(dat2,dat2.rep,len.tad.uniq[l],N.perm)
        rho.cutoff2[l,2] <- perm.result$rho.vec[length(perm.result$rho.vec)]
      }
    }else{
      sample.num <- min(round((length(len.tad.uniq)-15)*0.25),40)
      len.chosen <- c(4:15, sort(sample(16:length(len.tad.uniq),sample.num)))
      for(l in 1:length(len.chosen)){
        perm.result  <- perm(dat1,dat1.rep,len.tad.uniq[len.chosen[l]],N.perm)
        rho.cutoff1[len.chosen[l],2] <- perm.result$rho.vec[length(perm.result$rho.vec)]
        perm.result  <- perm(dat12,dat2.rep,len.tad.uniq[len.chosen[l]],N.perm)
        rho.cutoff2[len.chosen[l],2] <- perm.result$rho.vec[length(perm.result$rho.vec)]
      }
        
      spline.l <- smooth.spline(len.tad.uniq[len.chosen],rho.cutoff1[len.chosen,2])
      pred.result <- predict(spline.l,len.tad.uniq[16:length(len.tad.uniq)])
      rho.cutoff1[16:length(len.tad.uniq),2] <- pred.result$y
      spline.l <- smooth.spline(len.tad.uniq[len.chosen],rho.cutoff2[len.chosen,2])
      pred.result <- predict(spline.l,len.tad.uniq[16:length(len.tad.uniq)])
      rho.cutoff2[16:length(len.tad.uniq),2] <- pred.result$y
    }
    rho.cutoff[,2] <- 0.25*(rho.cutoff1[,2]+rho.cutoff2[,2])+0.5*rho.table[,ncol(rho.table)]
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
      if(tad.result[j,3]>rho.cutoff[which(rho.cutoff[,1]==tad.s.length),2]){
        tad.result[j,4] <- 0.5
      }else{
        tad.result[j,4] <- sum(rho.table[which(rho.table[,1]==tad.s.length),-1]<tad.result[j,3])/N.perm
      }
    }
  }
  tad.result[,5] <- p.adjust(tad.result[,4],method="BH")
  
  detection.result <- rep(0,nrow(genomic.interval))
  for(i in 1:nrow(genomic.interval)){
    s <- min(which(tad.result[,1]==genomic.interval[i,1]))
    e <- max(which(tad.result[,2]==genomic.interval[i,2]))
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
  
  return(list(tad.result=tad.result,genomic.result=genomic.result,rho.table=rho.table))
} 
