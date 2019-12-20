#clean
rm(list=ls())


get1Dfrom2D <- function( DD, SEP="~" ){

    GN <- colnames(DD)
    N  <- length(GN)

    S  = (N*(N-1))/2
    sa = 1
    M  = vector(length=S)
    
    for( i in 1:N ){
        for( j in i:N ){
            if( i != j ){
                names(M)[sa] = sprintf("%s%s%s",rownames(DD)[i],SEP,colnames(DD)[j])
                M[sa] = DD[i,j]
                sa=sa+1
                                                      
            }
        }
    }

    return(M)
    
}

firstOrder <- function(GG=NULL, CC=NULL, Ns=NULL, NCPUs=0, print=FALSE){
    library(Rcpp)
    dyn.load("CIParCor.so")
    res = .Call("firstOrder",GG=GG,CC=CC,N=Ns,nCPU=NCPUs,verbose=print)
    res
}


secondOrder <- function(GG=NULL, CC=NULL, Ns=NULL, NCPUs=0, print=FALSE){
    library(Rcpp)
    dyn.load("CIParCor.so")
    res = .Call("secondOrder",GG=GG,CC=CC,N=Ns,nCPU=NCPUs,verbose=print)
    res
}

#---Script required inputs
args  <- commandArgs(TRUE);
N     <- as.numeric(args[1]) #dummy network size 
cpus  <- as.numeric(args[2]) #no: of cpus

if( N < 0 || is.na(N) || is.null(N) ){ N = 100 }

if( cpus < 0 || is.na(cpus) || is.null(cpus) ){ cpus = 0 }

#print args
cat("args...\n")
cat(sprintf("N    = %d \n", N))
cat(sprintf("CPUs = %d \n\n", cpus))


TT = matrix(nrow=N,ncol=N)
GG = matrix(1, nrow=N,ncol=N)

colnames(TT) = sprintf("t%d",seq(1,N,1))
rownames(TT) = sprintf("t%d",seq(1,N,1))


for(i in 1:N){
    for(j in i:N){
        val = rnorm(1,0,0.25)
        if( abs(val) > 1 ){ val = sign(val) * 1 }
        TT[i,j] = val
        TT[j,i] = val
    }
}
#---

cat(sprintf("run firstOrder... \n"))

#calculate Conditional Independent first Order Partial Correlations
print(system.time(res.fO <- firstOrder(GG=GG, CC=TT, Ns=NULL, NCPUs=cpus, print=TRUE)))

cat(sprintf("...done.\n\n"))
cat(sprintf("run secondOrder... \n"))

#calculate Conditional Independent second Order Partial Correlations
print(system.time(res.sO <- secondOrder(GG=GG, CC=TT, Ns=NULL, NCPUs=cpus, print=TRUE)))

cat(sprintf("...done.\n\n"))

#store Conditional Independent Partial Correlations results
colnames(res.fO$PC) = colnames(TT)
rownames(res.fO$PC) = colnames(TT)

colnames(res.sO$PC) = colnames(TT)
rownames(res.sO$PC) = colnames(TT)

#
oo = cbind(get1Dfrom2D(TT),get1Dfrom2D(res.fO$PC),get1Dfrom2D(res.sO$PC))
colnames(oo) = c("TT","fO.PC","sO.PC")

#mean
cat(sprintf("print mean.\n\n"))
print(apply(oo, 2, mean))

#sd
cat(sprintf("print sd.\n\n"))
print(apply(oo, 2, sd))
