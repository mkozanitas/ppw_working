# row numbers of column x which are NAs
which(is.na(x))

# how many NAs are in row z
howmanyNAs <- function(z) length(which(is.na(z)))

xx <- c(1,1,NA,2,NA,NA,NA,3)
howmanyNAs(xx)

xx <- data.frame(z1=1:5,z2=2:6)
xx

# second position: 1 for rows, 2 for cols
apply(xx,2,sum)
apply(xx,1,howmanyNAs)
