
rejection=function(boxes,lines,n){
  # initial setting
  sample=matrix(NA,nrow=0,ncol=3)
  # draw box indices
  area = (boxes[,2]-boxes[,1])*(boxes[,4]-boxes[,3])
  prob = exp(boxes[,6]-min(boxes[,6]))
  
  N = 1.5*n
  
  while (N>0) {
    # for remaining samples
    boxindex = sample(seq(1:dim(boxes)[1]),N,replace=T,prob = area * prob)
    
    # calculate sample frequency for each index
    indexfreq=data.frame(table(boxindex))
    indexfreq$boxindex=as.numeric(levels(indexfreq$boxindex))[indexfreq$boxindex]
    for (j in 1:dim(indexfreq)[1]){
      # sample uniformly from the box
      sigmas=with(boxes,
                  runif(n=indexfreq[j,2],min=sigsqs.lo[indexfreq[j,1]],max=sigsqs.hi[indexfreq[j,1]]))
      sigmae=with(boxes,
                  runif(n=indexfreq[j,2],min=sigsqe.lo[indexfreq[j,1]],max=sigsqe.hi[indexfreq[j,1]]))
      eval.sigma=NULL
      for (i in 1:length(sigmas)){
        sigma = with(lines, a*sigmas[i]+b*sigmae[i])
        eval.sigma[i] = with(lines, sum(- 0.5*(multiplier.log * log(sigma) + multiplier.inv/sigma)))
      }
      acceptprob=exp(eval.sigma-boxes[indexfreq[j,1],6])
      
      # accept if alpha > u
      u=runif(length(sigmas))
      sample=rbind(sample,matrix(c(sigmas[acceptprob>u],sigmae[acceptprob>u],acceptprob[acceptprob>u]),ncol=3))
      
    }  
    # calculate remaining sample size
    N=n-dim(sample)[1]
    
  }
  output=as.data.frame(sample)
  colnames(output)= c("sigmas","sigmae","acceptance")
  return(output[sample(1:dim(output)[1],n),])
}