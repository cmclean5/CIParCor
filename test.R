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

vineBeta <- function(d, betaparam){
    # Reference:
    # https://stats.stackexchange.com/questions/124538/how-to-generate-a-large-full-rank-random-correlation-matrix-with-some-strong-cor
    
    P = matrix(0,d,d) #storing partial correlations
    S = diag(d);

    for(k in 1:d-1){
        for(i in (k+1):d){
            P[k,i] = rbeta(1,betaparam,betaparam) #sampling from beta
            P[k,i] = (P[k,i]-0.5)*2;     #linearly shifting to [-1, 1]
            p = P[k,i];
            for(l in (k):1){ #converting partial correlation to raw correlation
                p = p * sqrt((1-P[l,i]^2)*(1-P[l,k]^2)) + P[l,i]*P[l,k];
            }
            S[k,i] = p;
            S[i,k] = p;
        }
    }

    #permuting the variables to make the distribution permutation-invariant
    nr <- nrow(S)
    nc <- ncol(S)
    S  <- S[sample(nr), sample(nc)] 

    return(S)
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
Bs    <- as.numeric(args[3]) #beta dist. shape param.

if( N < 0 || is.na(N) || is.null(N) ){ N = 100 }

if( cpus < 0 || is.na(cpus) || is.null(cpus) ){ cpus = 0 }

if( Bs < 0 || is.na(Bs) || is.null(Bs) ){ Bs = 2 }

#print args
cat("args...\n")
cat(sprintf("N    = %d \n", N))
cat(sprintf("CPUs = %d \n", cpus))
cat(sprintf("Bs   = %d \n\n", Bs))


TT = vineBeta(N,2)
colnames(TT) = sprintf("t%d",seq(1,N,1))
rownames(TT) = sprintf("t%d",seq(1,N,1))


GG = matrix(1, nrow=N,ncol=N)


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
