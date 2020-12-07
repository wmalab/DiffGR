# simulation for noise matrix

generate_noise <- function(dat11){
  library(HiCcompare)
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





# simulation for single-TADs condition (simulation1)

#' @dat numeric. N*N HiC contact map 
#' @len.threshold numeric. The minimum TAD size to be considered (default value=5)
#' @tad.change.prop numeric. The proportion of altered TADs (default value=50%)
#' @random.prop numeric. The proportion of TAD alternation in each selected altered TADs (default value=100%)
#' @noise.level numeric. The proportion of noisy counts containing in simulated contacts (defalut value=10%)
#' @return a list that contains the altered TAD regions, TAD boundaries and simulated contact map


simulation1 <-  function(dat,len.threshold=5,tad.change.prop=0.5,random.prop=1,noise.level=0.1){
  tad <- HiCseg_linkC_R(nrow(dat),round(nrow(dat)/3),"P",dat,"D")
  tad <- tad$t_hat[tad$t_hat!=0]
  n <- round(length(tad)*tad.change.prop)
  tad.start  <- c(0,tad[-length(tad)])
  tad.length <- tad-tad.start
  tad.remain <- intersect(which(tad.length>(len.thre-1)),which(tad.length<max(tad.length)))
  tad.change.num <- sample(tad.remain,n)
  tad.change <- matrix(0,ncol=2,nrow=n)
  dat2 <- as.matrix(dat)
  for(i in 1:n){
    start <- tad.start[tad.change.num[i]]+1
    end <- tad[tad.change.num[i]]
    tad.change[i,] <- c(start,end)
    dat.s <- dat[start:end,start:end]
    len <- end-start+1
    m <- len*(len+1)/2
    if(round(random.prop*m)>0){
      posi <- sample(1:m,round(random.prop*m))
      for(j in 1:length(posi)){
        loci <- find.position(posi[j])
        diag.num <- abs(loci[1]-loci[2])
        diag.count <-  dat[row(dat)==(col(dat)-diag.num)]
        dat.s[loci[1],loci[2]] <- sample(diag.count,1) 
        dat.s[loci[2],loci[1]] <- dat.s[loci[1],loci[2]]
      }
    }
    dat2[start:end,start:end] <- dat.s
  }
  
  dat.noise <- generate_noise(dat)
  dat.sim <- (1-noise.level)*dat2+noise.level*dat.noise
  return(list(tad.change=tad.change,tad.boundary=tad,dat.sim=dat.sim))
}


# simulation for hierarchical-TADs condition (simulation2)

#' @dat numeric. N*N HiC contact map 
#' @tad numeric. TAD boundaries of dat(calling by HiCseg)
#' @len.threshold numeric. The minimum TAD size to be considered (default value=5)
#' @tad.change.prop numeric. The proportion of altered TADs (default value=50%)
#' @random.prop numeric. The proportion of TAD alternation in each selected altered TADs (default value=100%)
#' @noise.level numeric. The proportion of noisy counts containing in simulated contacts (defalut value=10%)
#' @return a list that contains the altered TAD regions, updated TAD boundaries and simulated contact map


simulation2 <- function(dat,len.threshold=5,random.prop=1,tad.change.prop=0.5,noise.level=0.1){
  tad <- HiCseg_linkC_R(nrow(dat),round(nrow(dat)/3),"P",dat,"D")
  tad <- tad$t_hat[tad$t_hat!=0]
  n <- round(length(tad)*tad.change.prop)
  tad.start  <- c(0,tad[-length(tad)])
  tad.length <- tad-tad.start
  tad.remain <- intersect(which(tad.length>(2*len.threshold-1)),which(tad.length<max(tad.length)))
  tad.change.num <- sample(tad.remain,n)
  tad.change <- matrix(0,ncol=2,nrow=n)
  add.boundaries <- c()
  dat2 <- as.matrix(dat)
  for(i in 1:n){
    start <- tad.start[tad.change.num[i]]+1
    end <- tad[tad.change.num[i]]
    tad.change[i,] <- c(start,end)
    dat.s <- dat[start:end,start:end]
    len <- end-start+1
    add.boundaries[i] <- sample(max(len.threshold,round(len/3)):min(len-len.threshold,round(2*len/3)),1)+tad1[tad.change.num[i]]
    m <- len*(len+1)/2
    if(round(random.prop*m)>0){
      posi <- sample(1:m,round(random.prop*m))
      for(j in 1:length(posi)){
        loci <- find.position(posi[j])
        diag.num <- abs(loci[1]-loci[2])
        diag.count <-  dat[row(dat)==(col(dat)-diag.num)]
        dat.s[loci[1],loci[2]] <- sample(diag.count,1) 
        dat.s[loci[2],loci[1]] <- dat.s[loci[1],loci[2]]
      }
    }
    dat2[start:end,start:end] <- dat.s
    dat2[start:add.boundaries[i],start:add.boundaries[i]] <- dat[start:add.boundaries[i],start:add.boundaries[i]]
    dat2[(add.boundaries[i]+1):end,(add.boundaries[i]+1):end] <- dat[(add.boundaries[i]+1):end,(add.boundaries[i]+1):end]
  }
  
  dat.noise <- generate_noise(dat)
  dat.sim <- (1-noise.level)*dat2+noise.level*dat.noise
  tad.boundary <- sort(c(tad,add.boundaries))
  return(list(tad.change=tad.change,tad.boundary=tad.boundary,dat.sim=dat.sim))
}




# simulation with different coverages

#' @dat numeric. N*N HiC contact map 
#' @sparse.level numeric. The relative coverage level of down-sampled contact matrix
#' @return a list that contains the simulated contact map


simulation3 <- function(dat,coverage.level){
  n <- nrow(dat)
  datt <- as.matrix(dat)
  datt[upper.tri(datt)] <- 0
  avail.posi <- which(datt!=0)
  for(j in 1:length(avail.posi)){
    count <- datt[avail.posi[j]]
    loci <- c()
    loci[1] <- avail.posi[j]%%n
    if(loci[1]==0){
      loci[1] <- n
      loci[2] <- avail.posi[j]%/%n
    }else{
      loci[2] <- 1+(avail.posi[j]%/%n)}
    datt[loci[1],loci[2]] <- rbinom(1,count,p=coverage.level)
    datt[loci[2],loci[1]] <- datt[loci[1],loci[2]]
  }
  dat.noise <- generate_noise(dat)
  dat.sim <- (1-0.1)*datt+0.1*dat.noise*coverage.level
  return(list(dat.sim=dat.sim))
}


