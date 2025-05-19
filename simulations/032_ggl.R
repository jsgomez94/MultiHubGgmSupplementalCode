##########################################################
##########################################################
##########################################################
##
## Relevant function: cv.hwgl
##  
##########################################################
##########################################################
##########################################################

examples = FALSE

##########################################################
##########################################################
## Function that chooses optimal Tuning Parameters for the
## GGL.
cv.ggl = function(X, K, p, n, penalize.diagonal = TRUE,
                  range = -5, rholength1 = 20, rholength2 = 20,
                  maxit = 10000, threshold = 10^(-7),
                  normalize = TRUE, AIC = TRUE){
  
  ##########################################################
  ##########################################################
  ## Step 1: prepare inputs:
  
  ######################
  ## Calculate normalized data:
  .X = NULL
  if(normalize){
    .normX = .multiple.normalize(X, K, p, n) 
    .X = .normX
  } else{
    .X = X
  }
  
  ######################
  ## Find upper limits 
  ## for tuning parameters:
  .rhos = .max.tuningpars(.X, K, p, n)
  
  ######################
  ## Generate estimates:
  
  .rholist1 = .rhos[[1]] * 10^seq(from = range, to = 0, length.out = rholength1)
  .rholist2 = .rhos[[2]] * 10^seq(from = range, to = 0, length.out = rholength2)
  
  .InfCrit1 = rep(0, rholength1)
  .InfCrit2 = rep(0, rholength2)

  ##########################################################
  ##########################################################
  ## Step 2: select the right value for rho1
  for(.ind1 in 1:rholength1){
    .ggl.temp = JGL(Y = .X, penalty = "group", 
                    lambda1 = .rholist1[.ind1], lambda2 = .rholist2[1], 
                    rho = 1, penalize.diagonal = penalize.diagonal, 
                    maxiter = maxit, tol = 10^(-7), return.whole.theta = TRUE)
    
    if(AIC){
      .InfCrit1[.ind1] = .AICemp.multi(s = .rhos[[3]], pm.est = .ggl.temp$theta, 
                                       K = K, p = p, n = n, 
                                       threshold = threshold)
    } else{
      .InfCrit1[.ind1] = .BICemp.multi(s = .rhos[[3]], pm.est = .ggl.temp$theta, 
                                       K = K, p = p, n = n, 
                                       threshold = threshold)
    }
  }
  .opt.index1 = which.min(.InfCrit1)
  .opt.rho1 = .rholist1[.opt.index1]
  
  ##########################################################
  ##########################################################
  ## Step 2: Select the right value of rho2
  for(.ind2 in 1:rholength2){
    .ggl.temp = JGL(Y = .X, penalty = "group", 
                    lambda1 = .rholist1[.opt.index1], lambda2 = .rholist2[.ind2], 
                    rho = 1, penalize.diagonal = penalize.diagonal, 
                    maxiter = maxit, tol = 10^(-7), 
                    return.whole.theta = TRUE)
    
    if(AIC){
      .InfCrit2[.ind2] = .AICemp.multi(s = .rhos[[3]], pm.est = .ggl.temp$theta, 
                                    K = K, p = p, n = n, 
                                    threshold = threshold)
    } else{
      .InfCrit2[.ind2] = .BICemp.multi(s = .rhos[[3]], pm.est = .ggl.temp$theta, 
                                    K = K, p = p, n = n, 
                                    threshold = threshold)
    }
  }
  .opt.index2 = which.min(.InfCrit2)
  .opt.rho2 = .rholist1[.opt.index2]

  ##########################################################
  ##########################################################
  ## Step 4: Save outputs:
  .opt.rho = c(.opt.rho1, .opt.rho2)
  .opt.index = c(.opt.index1,  .opt.index2)
  .opt.model = JGL(Y = .X, penalty = "group", 
                   lambda1 = .opt.rho[1], lambda2 = .opt.rho[2], 
                   rho = 1, penalize.diagonal = penalize.diagonal, 
                   maxiter = maxit, tol = 10^(-7), return.whole.theta = TRUE)
  
  .OUTPUT = list(opt.model = .opt.model,
                 opt.rho = .opt.rho,
                 opt.index = .opt.index,
                 rholist1 = .rholist1,
                 rholist2 = .rholist2,
                 InfCrit1 = .InfCrit1,
                 InfCrit2 = .InfCrit2)
  
  return(.OUTPUT)
}
##########################################################
##########################################################
## Function that chooses optimal Tuning Parameters for the
## GGL. It performs full-grid search, and returns the 
## pathway solution.
pathway.ggl = function(X, K, p, n, penalize.diagonal = TRUE,
                  range = -5, rholength1 = 20, rholength2 = 20,
                  maxit = 10000, threshold = 10^(-7),
                  normalize = TRUE, AIC = TRUE){
  
  ##########################################################
  ##########################################################
  ## Step 1: prepare inputs:
  
  ######################
  ## Calculate normalized data:
  .X = NULL
  if(normalize){
    .normX = .multiple.normalize(X, K, p, n) 
    .X = .normX
  } else{
    .X = X
  }
  
  ######################
  ## Find upper limits 
  ## for tuning parameters:
  .rhos = .max.tuningpars(.X, K, p, n)
  
  ######################
  ## Generate estimates:
  
  .rholist1 = .rhos[[1]] * 10^seq(from = range, to = 0, length.out = rholength1)
  .rholist2 = .rhos[[2]] * 10^seq(from = range, to = 0, length.out = rholength2)
  
  .evalmeasure = matrix(0, nrow = rholength1, ncol = rholength2)

  ##########################################################
  ##########################################################
  ## Step 2: Calculate the optimal path:
  .opt.models <- list()
  for (.ind1 in 1:rholength1) {

    for (.ind2 in 1:rholength2) {
      print(c(.ind1, .ind2))
      ## Calculate JGL solution:
      .ggl.temp = JGL(
        Y = .X, penalty = "group", 
        lambda1 = .rholist1[.ind1], lambda2 = .rholist2[.ind2], 
        rho = 1, penalize.diagonal = penalize.diagonal, 
        maxiter = maxit, tol = 10^(-7), return.whole.theta = TRUE)
      ## Calculate AIC/BIC of model:
      if (AIC) {
      .evalmeasure[.ind1, .ind2] <- .AICemp.multi(s = .rhos[[3]], pm.est = .ggl.temp$theta, 
                                                  K = K, p = p, n = n, 
                                                  threshold = threshold)
      } else{
      .evalmeasure[.ind1, .ind2] <- .BICemp.multi(s = .rhos[[3]], pm.est = .ggl.temp$theta, 
                                                  K = K, p = p, n = n, 
                                                  threshold = threshold)
      }
      ## Save/update best model:
      if (.evalmeasure[.ind1, .ind2] == min(.evalmeasure[.ind1, 1:.ind2])) {
        .opt.models[[.ind1]] <- list(
          model = .ggl.temp, 
          rho1 = .rholist1[.ind1], 
          rh2 = .rholist2[.ind2])
      }
    }
  }
  print(.evalmeasure)
  .opt.index = which(.evalmeasure == min(.evalmeasure), arr.ind = TRUE)
  .opt.rho = c(.rholist1[.opt.index[1]], .rholist2[.opt.index[2]])
  
  ##########################################################
  ##########################################################
  ## Step 4: Save outputs:
  
  .OUTPUT = list(opt.models = .opt.models,
                 opt.rhos = .opt.rho,
                 opt.index = .opt.index,
                 rholist1 = .rholist1,
                 rholist2 = .rholist2,
                 evalmeasure = .evalmeasure)
  
  return(.OUTPUT)
}












##########################################################
## Normalize data:
.multiple.normalize = function(X, K, p, n){
  .normX = list()
  for(.k in 1:K){
    .normX[[.k]] = scale(x = X[[.k]], center = TRUE, scale = TRUE)
    .normX[[.k]] = .normX[[.k]][1:n[.k], 1:p]
  }
  return(.normX)
}

# library(mvtnorm)
# 
# X = list()
# X[[1]] = matrix(rnorm(200)+ 10, ncol = 10, nrow = 20)
# X[[2]] = matrix(rnorm(200)+ 10, ncol = 10, nrow = 20)
# X[[3]] = matrix(rnorm(200)+ 10, ncol = 10, nrow = 20)
# 
# Xnorm = .multiple.normalize(X = X, K = 3, p = 10, n = 20)
# Xnorm[[1]][1:20, 1:10]
# 
# b = scale(X) 
# b

##########################################################
## Find maximal choices of tuning parameters:
.max.tuningpars = function(X, K, p, n){
  
  ## Generate scaled correlation:
  .cov = list()
  .covN = array(data = 0, dim = c(p,p,K))
  for(.k in 1:K){
    .cov[[.k]] = cov(X[[.k]])
    .covN[, , .k] = n[.k] * .cov[[.k]]
    .covN[, , .k] = abs( .covN[, , .k] - diag(diag(.covN[, , .k])) )
  }
  
  ## Find max val of rho1
  .rho1 = max( .covN ) 
  
  ## Find max val for rho2
  .covnormN = matrix(0, ncol = p, nrow = p)
  for(.i in 2:p){
    for(.j in 1:(.i-1)){
      .covnormN[.i, .j] = sqrt(sum( .covN[.i, .j, ]^2 ))
    }
  }
  .rho2 = max(.covnormN)
  
  ## Save outputs:
  .output = list(maxrho1 = .rho1, maxrho2 = .rho2, cov = .cov)
  
  return(.output)

}

##########################################################
## This function computes the BIC of a given estimation:
.BICemp.multi = function(s, pm.est = NULL, w.est = NULL, p, n,
                     threshold = 10^(-7)){
  .bic <- 0
  .pm.est <- list()
  
  #### FIRST PORTION: LOGDET PART.
  ## Case 1: W is provided:
  if(!is.null(w.est)){
    
    for(.k in 1:K){
      .pm.estk[[.k]] = solve(w.est[[.k]])
      .eigwk = eigen(w.est[[.k]], symmetric = TRUE)
      
      .bic  <- .bic + n[.k] * sum( log(.eigw$values) )
    }
  }
  ## Case 2: Only Theta is provided:
  if(is.null(w.est) & !is.null(pm.est)){
    .pm.est <- pm.est
    
    for(.k in 1:K){
      .eigt <- eigen(.pm.est[[.k]])
      .bic <- .bic - n[.k] * sum( log(.eigt$values) )
    }
  }
  
  for(.k in 1:K){
    #### SECOND PORTION: TRACE PART.
    .bic <- .bic + n[.k] * Trace(.pm.est[[.k]] %*% s[[.k]])
    
    #### THIRD PORTION: COMPLEXITY QUANTIFICATION.
    .edgenumk = sum( abs(.pm.est[[.k]]) > threshold )/2
    .bic <- .bic + log(n[.k]) * .edgenumk 
  }
  
  #### Conclude.
  return(.bic)

}

##########################################################
## This function computes the AIC of a given estimation:
## BIC function:
.AICemp.multi <- function(s, pm.est = NULL, w.est = NULL, K, p, n,
                   threshold = 10^(-7)){
  .aic <- 0
  .pm.est <- list()
  
  #### FIRST PORTION: LOGDET PART.
  ## Case 1: W is provided:
  if(!is.null(w.est)){
    
    for(.k in 1:K){
      .pm.estk[[.k]] = solve(w.est[[.k]])
      .eigwk = eigen(w.est[[.k]], symmetric = TRUE)
      
      .aic  <- .aic + n[.k] * sum( log(.eigw$values) )
    }
  }
  ## Case 2: Only Theta is provided:
  if(is.null(w.est) & !is.null(pm.est)){
    .pm.est <- pm.est
    
    for(.k in 1:K){
      .eigt <- eigen(.pm.est[[.k]])
      .aic <- .aic - n[.k] * sum( log(.eigt$values) )
    }
  }
  
  for(.k in 1:K){
    #### SECOND PORTION: TRACE PART.
    .aic <- .aic + n[.k] * Trace(.pm.est[[.k]] %*% s[[.k]])

    #### THIRD PORTION: COMPLEXITY QUANTIFICATION.
    .edgenumk = sum( abs(.pm.est[[.k]]) > threshold )/2
    .aic <- .aic + 2*.edgenumk
    
  }
  
  #### Conclude.
  return(.aic)
}




##########################################################
##########################################################
if(examples){
  source("002_GeneratingMultipleMatrixSparse.R")
  source("003_UsefulMatrixTransforms.R")
  source("013_BIC_OLD.R")
  
  library(mvtnorm)
  library(JGL)
  
  #################################################
  #################################################
  ## Establishing parameters.
  H1 = c(1,2,3,4,5)
  H2.1 = NULL
  H2.2 = NULL
  H2.3 = NULL
  p = 100
  T0 = 100
  
  ph1 = 0.7
  ph2 = 0
  pnh = 0.05
  pneff = 0.01
  
  diagonal_shift = 3
  
  n1 = 150
  n2 = 150
  n3 = 150
  n = c(n1, n2, n3)
  
  
  #################################################
  #################################################
  ## Generate the three precision matrices with hubs.
  Theta1 = r.sparse.pdhubmat(p = p, T0 = T0, 
                             H1 = H1, H2 = H2.1, 
                             ph1 = ph1, ph2 = ph2, pnh = pnh, pneff = pneff,
                             diagonal_shift = diagonal_shift, shuffle = FALSE, type = "unif",
                             hmin1 = 4, hmax1 = 5, 
                             hmin2 = 4, hmax2 = 5, 
                             nhmin = 4, nhmax = 5)
  Theta2 = r.sparse.pdhubmat(p = p, T0 = T0, 
                             H1 = H1, H2 = H2.2, 
                             ph1 = ph1, ph2 = ph2, pnh = pnh,  pneff = pneff,
                             diagonal_shift = diagonal_shift, shuffle = FALSE, type = "unif",
                             hmin1 = 4, hmax1 = 5, 
                             hmin2 = 4, hmax2 = 5, 
                             nhmin = 4, nhmax = 5)
  Theta3 = r.sparse.pdhubmat(p = p, T0 = T0,  
                             H1 = H1, H2 = H2.3, 
                             ph1 = ph1, ph2 = ph2, pnh = pnh,  pneff = pneff,
                             diagonal_shift = diagonal_shift, shuffle = FALSE, type = "unif",
                             hmin1 = 4, hmax1 = 5, 
                             hmin2 = 4, hmax2 = 5, 
                             nhmin = 4, nhmax = 5)
  
  ## Obtain covariance matrix
  cov1 = solve(Theta1)
  cov2 = solve(Theta2)
  cov3 = solve(Theta3)
  
  #################################################
  #################################################
  ## Generate data and estimate.
  X1 = rmvnorm(n = n1, mean = rep(0,p), sigma = cov1)
  X2 = rmvnorm(n = n2, mean = rep(0,p), sigma = cov2)
  X3 = rmvnorm(n = n3, mean = rep(0,p), sigma = cov3)
  
  X = list(X1, X2, X3)
  
  ########################
  ## Step 1: 
  ## Verify normalization:
  exp = .multiple.normalize(X = X, K = 3, p = p, n = n)
  
  ########################
  ## Step 1: 
  ## Verify normalization:
  maxtp = .max.tuningpars(X = exp, K = 3, p = p, n = n)
    
  maxtp[[1]]
  maxtp[[2]]
  diag(maxtp[[3]][[1]])
  
  

  
  ########################
  ## Step 1: 
  ## Verify estimation method:
  BICmodel = cv.ggl(X = X, K = 3, p = p, n = n, 
                    penalize.diagonal = TRUE, range = -7.5, 
                    rholength1 = 20, rholength2 = 20, 
                    maxit = 1000, threshold = 10^(-7), 
                    normalize = TRUE, AIC = TRUE)
  AICmodel = cv.ggl(X = X, K = 3, p = p, n = n, 
                    penalize.diagonal = TRUE, range = -7.5, 
                    rholength1 = 20, rholength2 = 20, 
                    maxit = 1000, threshold = 10^(-7), 
                    normalize = TRUE, AIC = TRUE)
  
  plot(BICmodel$opt.model$theta[[1]]!= 0)
  plot(BICmodel$opt.model$theta[[2]]!= 0)
  plot(BICmodel$opt.model$theta[[3]]!= 0)
  
  plot(AICmodel$opt.model$theta[[1]]!= 0)
  plot(AICmodel$opt.model$theta[[2]]!= 0)
  plot(AICmodel$opt.model$theta[[3]]!= 0)
  
  
  BICpathway = pathway.ggl(X = X, K = 3, p = p, n = n, 
                           penalize.diagonal = TRUE, range = -5, 
                           rholength1 = 5, rholength2 = 5, 
                           maxit = 10, threshold = 10^(-7), 
                           normalize = TRUE, AIC = TRUE)
  AICpathway <- pathway.ggl(X = X, K = 3, p = p, n = n, 
                            penalize.diagonal = TRUE, range = -5, 
                            rholength1 = 5, rholength2 = 5, 
                            maxit = 10, threshold = 10^(-7), 
                            normalize = TRUE, AIC = TRUE)

  ## BICpathway$opt.models:
    ## List of rholength1 GGL models.
    ## BICpathway$opt.models[[t]]           is the optimal model corresponding to (lambda1[t], lambda2[opt]).
    ## (BICpathway$opt.models[[t]])$theta   are the K matrices resulted from GGL with (lambda1[t], lambda2[opt]).

  names(BICpathway)                           ## Full output object
  print((BICpathway$opt.models[[1]]))         ## Object corresponding to lambda1[1] 
  ((BICpathway$opt.models[[1]])$model)$theta  ## Theta estimators for lambda1[1]
  
  print((BICpathway$opt.models[[1]]))         ## Object corresponding to lambda1[2] 
  ((BICpathway$opt.models[[1]])$model)$theta  ## Theta estimators for lambda1[2]


  
  rm(H1, H2.1, H2.2, H2.3, p, t, ph1, ph2, pnh, 
     diagonal_shift, n1, n2, n3, n,
     Theta1, Theta2, Theta3, cov1, cov2, cov3, X1, X2, X3, X,
     exp, maxtp, BICmodel, AICmodel)
}



rm(examples)

