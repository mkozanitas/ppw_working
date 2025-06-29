# Adjust UTM to mazimize trees inside corner coordinates
# assuming the following files are available
# tmp - Plot,Num,UTM.E.17,UTM.N.17
# bounds - xlow,xhigh,ylow,yhigh for corner posts

howManyInside <- function(xy,b=bounds,xadj=0,yadj=0) {
  Ntot <- nrow(xy)
  inTrees <- which(xy[,1]+xadj>=b[1]
          & xy[,1]+xadj <= b[2]
          & xy[,2]+yadj >= b[3]
          & xy[,2]+yadj <= b[4])
  inNum <- length(inTrees)
  return(list(Ntot,inNum,inTrees))
}

# xadj <- 5
# yadj <- 5
# res <- howManyInside(tmp[,c('UTM.E.17','UTM.N.17')],bounds,xadj,yadj)
# print(c(res[[1]],res[[2]],res[[2]]/res[[1]]))
# 
# drawBox(b,exlim=50)
# points(tmp$UTM.E.17+xadj,tmp$UTM.N.17+yadj)
# points(tmp$UTM.E.17[res[[3]]]+xadj,tmp$UTM.N.17[res[[3]]]+yadj,pch=19)


adjRes <- data.frame(xadj=rep(0:20,21),yadj=rep(0:20,each=21),nTot=NA,inNum=NA)

i=399
for (i in 1:nrow(adjRes)) {
  res <- howManyInside(tmp[,c('UTM.E.17','UTM.N.17')],b,adjRes[i,1],adjRes[i,2])
  adjRes[i,3] <- res[[1]]
  adjRes[i,4] <- res[[2]]
}
head(adjRes)
(optRes <- adjRes[which.max(adjRes[,4]),])

# drawBox(b,exlim=50)
# points(tmp$UTM.E.17+as.numeric(optRes[1]),tmp$UTM.N.17+as.numeric(optRes[2]))
