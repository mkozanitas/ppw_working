allequal <- function(x) {
  x <- x[-which(is.na(x))]
  x <- x[-which(x=='<NA>')]
  x <- x[-which(x=='')]
  sp <- x[1]
  all(x==sp)
}
x <- c('x','x',NA)
allequal(x)

tAllm[25,c('Species.13','Species.17','Species.18','Species.20')]

speq <- apply(tAllm[,c('Species.13','Species.17','Species.18','Species.20')],1,allequal)
which(!speq)