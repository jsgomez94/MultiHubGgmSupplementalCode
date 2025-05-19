##########################################################
##########################################################
##########################################################
##
## Relevant function: BICemp
##  
##########################################################
##########################################################
##########################################################

examples = FALSE

##########################################################
##########################################################

## Trace function
Trace = function(M){
  return(sum(diag(M)))
}

## BIC function:
BICemp <- function(s, pm.est = NULL, w.est = NULL, p, n,
                   threshold = 10^(-7)){
  .bic <- 0
  .pm.est <- NULL
  
  #### FIRST PORTION: LOGDET PART.
  ## Case 1: W is provided:
  if(!is.null(w.est)){
    .pm.est = solve(w.est)
    
    .eigw <- eigen(w.est, symmetric = TRUE)
    .bic  <- .bic + n * sum( log(.eigw$values) )
  }
  ## Case 2: Only Theta is provided:
  if(is.null(w.est) & !is.null(pm.est)){
    .pm.est <- pm.est
    
    .eigt <- eigen(.pm.est)
    .bic <- .bic - n * sum( log(.eigt$values) )
  }
  
  #### SECOND PORTION: TRACE PART.
  .bic <- .bic + n * Trace(.pm.est %*% s)
  
  #### THIRD PORTION: COMPLEXITY QUANTIFICATION.
  .edgenum = sum( abs(.pm.est) > threshold )/2
  .bic <- .bic + log(n)*.edgenum
  
  #### Conclude.
  return(.bic)
}



##########################################################
##########################################################

if(examples){
  theta = r.sparse.pdhubmat(p = 40, t = 20, 
                            r1 = 2, r2 = 2, 
                            ph1 = 0.6, ph2 = 0.6, pnh = 0.05, 
                            diagonal_shift = 10, shuffle = TRUE, 
                            type = "unif", 
                            hmin1 = 4, hmax1 = 5, 
                            hmin2 = 4, hmax2 = 5, 
                            nhmin = 4, nhmax = 5)
  sigma = solve(theta)
  
  X = rmvnorm(n = 200, sigma = sigma)
  
  s = cov(X)
  maxlambda = max( abs(s - diag(diag(s))) )
  rholist = maxlambda * 10^seq(from = -10, to =0, length.out = 100)
  bic1 = rep(0,100)
  bic2 = rep(0,100)
  bic3 = rep(0,100)
  bic4 = rep(0,100)
  for(.i in 1:100){
    
    .model = glasso(s = cov(X), rho = rholist[.i], nobs = 200, 
                    penalize.diagonal = TRUE)
    
    bic1[.i] = BICemp(s = cov(X), w.est = .model$w, 
                      p = 40, n = 200, 
                      threshold = 10^(- 5))
    bic2[.i] = BICemp(s = cov(X), w.est = .model$w, 
                      p = 40, n = 200, 
                      threshold = 10^(- 7))
    bic3[.i] = BICemp(s = cov(X), w.est = .model$w, 
                      p = 40, n = 200, 
                      threshold = 10^(- 10))
    bic4[.i] = BICemp(s = cov(X), pm.est = .model$wi, 
                      p = 40, n = 200, 
                      threshold = 10^(- 7))
  }
  plot(x = rholist, y = bic1, type = "l", col = "red")
  abline(v = rholist[which.min(bic1)], col = "red")
  lines(x = rholist, y = bic2, col = "blue")
  abline(v = rholist[which.min(bic2)], col = "blue")
  lines(x = rholist, y = bic3, col = "green")
  abline(v = rholist[which.min(bic3)], col = "green")
  lines(x = rholist, y = bic4, col = "black")
  abline(v = rholist[which.min(bic4)], col = "black")
  
  ## As we can see from this example, the BIC does become
  ## the same if we choose an appropriate threshold. The
  ## question is: what is the most appropriate threshold?
  
  rm(bic1, bic2, bic3, bic4,
     theta, sigma, X, s, .i, maxlambda, rholist, .model)
}

rm(examples)

## Rank BIC (not useful.)
# BICrank = function(s, pm.est, p, n, threshold){
#   .pm.est = (pm.est + t(pm.est)) / 2
#   .eig = eigen(.pm.est, symmetric = TRUE)
# 
#   .bic <- 0
#   .bic <- .bic - n * sum( log(.eig$values) )
#   .bic <- .bic + n * Trace(.pm.est %*% s)
# 
#   print(.bic)
# 
#   .rank = 0
#   .pseudoL = PseudoLaplacian(Theta = pm.est, p = p)
#   .pseudoL = (.pseudoL + t(.pseudoL))/(2)
#   .eigL = eigen(.pseudoL, symmetric = TRUE)
# 
#   .K = NULL
#   if(.eigL$values[1] == 0)
#     .K = 0
#   if(.eigL$values[1] != 0)
#     .K = sum(.eigL$values)/.eigL$values[1]
# 
#   print(.K)
# 
#   .bic <- .bic + log(n)* (2 * p - .K)* .K
#   return(.bic)
# 
# }


