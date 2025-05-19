##########################################################
##########################################################
##########################################################
##
## Relevant function: cv.glasso.
##  
##########################################################
##########################################################
##########################################################

examples = FALSE

##########################################################
## Method for estimating hubs
##########################################################

cv.glasso = function(s, n, p, penalize.diagonal = TRUE,
                     range = -5, rholength = 100,
                     maxit = 10000, bic.threshold = 10^(-5),
                     is.corr = TRUE){
  ## Step 1: find the minimum rho for which wi is diagonal:
  .maxrho = maxtuning.glasso(s = s, p = p, n = n, 
                             penalize.diagonal = penalize.diagonal,
                             is.corr = is.corr)
  
  ## Step 2: Create the glasso fit for the range provided:
  .rholist = .maxrho * 10^(seq(from = range, to = 0, length.out = rholength))
  .BIC = rep(0, rholength)
  .glasso.temp = NULL
  for(.i in 1:rholength){
    
    if(.i == 1){
      .glasso.temp = glasso(s = s, rho = .rholist[.i], 
                            nobs = n, zero = NULL, 
                            thr = 1.0e-4, maxit = maxit,  approx = FALSE, 
                            penalize.diagonal = penalize.diagonal, 
                            start = "cold", w.init = NULL, wi.init = NULL, 
                            trace = FALSE)
      
    } else {
      .glasso.temp = glasso(s = s, rho = .rholist[.i], 
                            nobs = n, zero = NULL, 
                            thr = 1.0e-4, maxit = maxit,  approx = FALSE, 
                            penalize.diagonal = penalize.diagonal, 
                            start = "warm",
                            w.init  = .glasso.temp$w, 
                            wi.init = .glasso.temp$wi, 
                            trace = FALSE)
      
    }
    .BIC[.i] = BICemp(s = s, w.est = .glasso.temp$w,
                      n = n, p = p, 
                      threshold = bic.threshold)
  }
  
  ## Step 3: Refit to obtain the initializing estimator.
  .opt.ind = which.min(.BIC)
  .opt.rho = .rholist[.opt.ind]
  
  .opt.model = glasso(s = s, rho = .opt.rho, nobs = n, 
                      penalize.diagonal = penalize.diagonal)
  
  .OUTPUT = list(opt.model = .opt.model,
                 opt.rho = .opt.rho,
                 opt.index = .opt.ind,
                 rholist = .rholist,
                 bic = .BIC)
  return(.OUTPUT)
}



##########################################################
##########################################################
## This function finds the smallest tuning parameter for
## which the GLASSO estimate is diagonal.

maxtuning.glasso = function(s = NULL, true.pm = NULL, use.cov = TRUE,
                            p, n, 
                            penalize.diagonal = TRUE,
                            is.corr = TRUE) {
  .s = NULL
  if(is.null(s) & !is.null(true.pm)){
    .sigma = solve(pm)
    .X = rmvnorm(n = n, 
                 sigma = .sigma,
                 method = "svd")
    if(use.cov) cov(.X) 
    if(!use.cov) cor(.X)
  } else {
    .s = s
  }
  .rho = 10^30
  if(is.corr){
    .rho = 1
  } 
  
  ##############################################
  ##############################################
  ## Estimate the value for which the model is 
  ## first diagonal.
  .output.glasso = glasso(s = .s, rho = .rho, 
                        nobs = n, maxit = 200, 
                        penalize.diagonal = penalize.diagonal)
  .is.diagonal = is.diagonal.matrix(.output.glasso$wi)
  
  .count = 1
  while(!.is.diagonal & .count < 5){
    .count = .count + 1
    .rho = 10^5 * .rho
    .output.hwgl = glasso(s = .s, rho = .rho*.W, 
                          nobs = n, maxit = 200, 
                          penalize.diagonal = penalize.diagonal)
    .is.diagonal = is.diagonal.matrix(.output.glasso$wi)
  }
  if(!.is.diagonal){
    print("Warning: No diagonal/non-diagonal threshold was found: lower bound only.")
    return(.rho)
  }
  
  .count = 1
  while(.is.diagonal & .count < 201){
    ## Reduce the value of the tuning parameter
    .count = .count + 1 
    .rho = .rho/2
    
    ## find the new model:
    .output.glasso = glasso(s = .s, rho =  .rho, 
                          nobs = n, maxit = 200, 
                          penalize.diagonal = penalize.diagonal)
    .is.diagonal = is.diagonal.matrix(.output.glasso$wi)
    
  }   
  if(.is.diagonal){
    print("Warning: No diagonal/non-diagonal threshold was found: upper bound only.")
    return(.rho)
  }  
  
  
  print(paste("The diagonal/non-diagonal treshold is:", .rho))
  return(.rho)
}


##########################################################
##########################################################
if(examples){
  
  theta = r.sparse.pdhubmat(p = 40, t = 20, 
                            r1 = 2, r2 = 2, 
                            ph1 = 0.6, ph2 = 0.6, pnh = 0.05, 
                            diagonal_shift = 10, shuffle = FALSE, 
                            type = "unif", 
                            hmin1 = 4, hmax1 = 5, 
                            hmin2 = 4, hmax2 = 5, 
                            nhmin = 4, nhmax = 5)
  sigma = solve(theta)
  
  X = rmvnorm(n = 200, sigma = sigma)
  
  s = cov(X)
  
  model = cv.glasso(s = s, n = 200, p = 40, 
                    penalize.diagonal = TRUE,
                    range = -10, rholength = 100, 
                    maxit = 10000)
  
  plot(x = log(model$rholist, base = 10), y = model$bic)
  
  par(mfrow = c(2,1), oma = c(0,0,0,4))
  plot(theta!= 0)
  plot(model$.opt.model$wi != 0)

  rm(theta, sigma, X, s, model)    
}



rm(examples)
