#################################################
#################################################
#################################################
##
## In the following document, we will
## define a method for generating positive
## definite matrices with hubs, which is 
## based on first generating the eigen-
## structure.
##
#################################################
#################################################

examples = FALSE

#################################################
#################################################
## Function that generates the adjacency matrix of
## a network with hubs.
.adjmat <- function(p, t, H1, H2, ph1, ph2, pnh){
  .A = matrix(rep(0, p*p), ncol = p)
  
  for(.i in 2:t){
    for(.j in 1:(.i-1)) {
      if(.i %in% H1 || .j %in% H1) 
        .A[.i, .j] = rbinom(1,1,ph1)
      else if(.i %in% H2 || .j %in% H2) 
        .A[.i, .j] = rbinom(1,1,ph2)
      else {
        .A[.i, .j] = rbinom(1,1,pnh)
      }
    }
  }
  .A = .A + t(.A)
  return(.A)
}


.rsign <- function(n = 1) {
  return(2*rbinom(n,1,0.5)-1)
}

## Generates a random variable with distribution:
## Unif( [-max,-min]U[min,max] )
.rsignunif <- function(n = 1 , min = 4, max = 5) {
  .cont = runif(n, min, max)
  .sign = .rsign(n)
  return(.cont*.sign )
}

## Generates a symmetric matrix with entries given
## by rsignunif.
.rsymmmatrix <- function(p, t, H1, H2, 
                         type = c("unif", "gaussian"),
                         hmin1 = 0.5, hmax1 = 0.8,
                         hmin2 = 0.5, hmax2 = 0.8,
                         nhmin = 0.5, nhmax = 0.8,
                         hsd1 = 1, hsd2 = 1, nhsd = 1){
  
  .theta = matrix(rep(0, p*p), ncol = p)
  
  if(type == "unif"){
    for(.i in 2:t){
      for(.j in 1:(.i-1)) {
        if(.i %in% H1 || .j %in% H1) {
          .theta[.i, .j] = .rsignunif(n = 1, min = hmin1, max = hmax1)
        } else if(.i %in% H2 || .j %in% H2) {
          .theta[.i, .j] = .rsignunif(n = 1, min = hmin2, max = hmax2)
        } else{
          .theta[.i, .j] = .rsignunif(n = 1, min = nhmin, max = nhmax)
        } 
      }
    }
  }
  if(type == "gaussian"){
    for(.i in 2:t){
      for(.j in 1:(.i-1)) {
        if(.i %in% H1 || .j %in% H1) {
          .theta[.i, .j] = rnorm(n = 1, sd = hsd1)
        } else if(.i %in% H2 || .j %in% H2) {
          .theta[.i, .j] = rnorm(n = 1, sd = hsd2)
        } else{
          .theta[.i, .j] = rnorm(n = 1, sd = nhsd)
        } 
      }
    }
  }
  .theta = .theta + t(.theta)
  return(.theta)
}

## Shuffles the position of the hubs from the
## first couple entries to any random entry.
.shufflemat <- function(A,p){
  .neworder = sample(1:p, p, replace = FALSE)
  .M = A[.neworder, .neworder]
  return(.M)
}

## Generates a symmetric matrix that is not
## positive definite. It has zero diagonals.
.rhubmat <- function(p, t, H1, H2, ph1, ph2, pnh, shuffle = FALSE,
                     type = c("unif", "gaussian"),
                     hmin1 = 0.5, hmax1 = 0.8,
                     hmin2 = 0.5, hmax2 = 0.8,
                     nhmin = 0.5, nhmax = 0.8,
                     hsd1 = 5, hsd2 = 5, nhsd = 1){
  .A = .adjmat(p, t, H1, H2, ph1, ph2, pnh)
  .theta = .A * .rsymmmatrix(p, t, H1, H2,
                             type = type,
                             hmin1 = hmin1, hmax1 = hmax1,
                             hmin2 = hmin2, hmax2 = hmax2,
                             nhmin = nhmin, nhmax = nhmax,
                             hsd1 = hsd1, hsd2 = hsd2, nhsd = nhsd)
  if( shuffle ){
    .theta = .shufflemat(.theta, p)
  }
  return( .theta )
}

## With this functions, we are able to select the 
## network structure of the precision matrix. Now,
## we have to design methods that select the value
## of the diagonal.

######################
######################
### Example:
if(examples){

  adjmat = .adjmat(p = 20, t = 10, 
                   H1 = c(10), H2 = c(4,5), 
                   ph1 = 1, ph2 = 0.5, pnh = 0)
  par(oma = c(0,0,0,2))
  plot(adjmat) 
  
  ## We generate two conventional matrices with
  ## two types of hubs. 
  rhubmat = .rhubmat(p = 20, t = 10, 
                     H1 = c(5), H2 = c(1,2), 
                     ph1 = 1, ph2 = 0.5, pnh = 1, 
                     shuffle = FALSE,
                     type = "unif", 
                     hmin1 = 100, hmax1 = 101, 
                     hmin2 = 50, hmax2 = 51, 
                     nhmin = 1, nhmax = 2)
  plot(rhubmat)
  rhubmat = .rhubmat(p = 10, t = 8, 
                     H1 = c(8), H2 = c(1,2), 
                     ph1 = 1, ph2 = 0.5, pnh = 1, 
                     shuffle = FALSE,
                     type = "gaussian", 
                     hsd1 = 100, hsd2 = 100, nhsd = 1) 
  plot(abs(rhubmat))
  
  
  
  ## Lets see what happens if r1 = 0
  rhubmat = .rhubmat(p = 20, t = 10, 
                     H1 = NULL, H2 = c(1,2), 
                     ph1 = 1, ph2 = 0.5, pnh = 1, 
                     shuffle = FALSE,
                     type = "unif", 
                     hmin1 = 100, hmax1 = 101, 
                     hmin2 = 50, hmax2 = 51, 
                     nhmin = 1, nhmax = 2)
  par(mar = c(4,4,4,4))
  plot(rhubmat)

  ## Lets see what happens if r2 = 0
  rhubmat = .rhubmat(p = 20, t = 10, 
                     H1 = c(4,5), H2 = NULL, 
                     ph1 = 1, ph2 = 0.5, pnh = 1, 
                     shuffle = FALSE,
                     type = "unif", 
                     hmin1 = 100, hmax1 = 101, 
                     hmin2 = 50, hmax2 = 51, 
                     nhmin = 1, nhmax = 2)
  par(mar = c(4,4,4,4))
  plot(rhubmat)
  
  rm(adjmat, rhubmat, examples)
  
}



#################################################
#################################################

## This function generates precision matrices by adding a 
## single uniform diagonal to all variables.
r.sparse.pdhubmat <- function(p, t, H1, H2, ph1, ph2, pnh, 
                              diagonal_shift = 1, 
                              shuffle = FALSE,
                              type = c("unif", "gaussian"),
                              hmin1 = 0.5, hmax1 = 0.8,
                              hmin2 = 0.5, hmax2 = 0.8,
                              nhmin = 0.5, nhmax = 0.8,
                              hsd1 = 1, hsd2 = 1, nhsd = 1,
                              verbose = FALSE){
  .pm = .rhubmat(p, t, H1, H2, ph1, ph2, pnh, shuffle = shuffle,
                 type = type,
                 hmin1 = hmin1, hmax1 = hmax1,
                 hmin2 = hmin2, hmax2 = hmax2,
                 nhmin = nhmin, nhmax = nhmax,
                 hsd1 = hsd1, hsd2 = hsd2, nhsd = nhsd)

  .lambda = eigen(.pm)$value[p]
  .I =  (-.lambda + diagonal_shift) * diag(p)
  
  if(verbose){
    print(paste("The dimension of PM is:",dim(.pm)[1],"x", dim(.pm)[2]))
    print(paste("The dimension of .I is:",dim(.I)[1],"x", dim(.I)[2]))
  }
  
  .pm = .pm + .I
  return(.pm)
}

r.laplacian.pdhubmat <- function(p, t, H1, H2, ph1, ph2, pnh, 
                           diagonal_shift = 1, 
                           shuffle = FALSE,
                           type = c("unif", "gaussian"),
                           hmin1 = 0.5, hmax1 = 0.8,
                           hmin2 = 0.5, hmax2 = 0.8,
                           nhmin = 0.5, nhmax = 0.8,
                           hsd1 = 1, hsd2 = 1, nhsd = 1,
                           verbose = FALSE){
  .pm = .rhubmat(p, t, H1, H2, ph1, ph2, pnh, shuffle = shuffle,
                 type = type,
                 hmin1 = hmin1, hmax1 = hmax1,
                 hmin2 = hmin2, hmax2 = hmax2,
                 nhmin = nhmin, nhmax = nhmax,
                 hsd1 = hsd1, hsd2 = hsd2, nhsd = nhsd)
  .D = apply(abs(.pm), MARGIN = 1, sum)
  
  .pm = diag(.D) - .pm + diagonal_shift* diag(p)
  return(.pm)
}


######################
######################
### Example:
if(examples){
  Theta = r.sparse.pdhubmat(p = 20, t = 10, r1 = 1, r2 = 2, ph1 = 1, ph2 = 0.5, pnh = 1, 
                            diagonal_shift = 1,  shuffle = FALSE,
                            type = "unif", 
                            hmin1 = 100, hmax1 = 101, 
                            hmin2 = 50, hmax2 = 51, 
                            nhmin = 1, nhmax = 2)

  plot(abs(Theta - diag(diag(Theta))))
  
  rm(Theta)
  
}

rm(examples)
