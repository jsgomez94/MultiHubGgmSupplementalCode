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


cv.hwgl = function(s, p, n, penalize.diagonal = TRUE,
                   range = -5, rholength = 100, 
                   maxit = 10000, bic.threshold = 10^(-5)){
  
  ##########################################################
  ## Find penalty weight matrix
  .inv_cov = solve(s + 0.1 * diag(p) )
  .W  <- matrix(rep(0, p * p), ncol = p)
  for (.i in (2:p)) {
    for(.j in (1:(.i-1))) {
      .a  <- abs(.inv_cov[.i, .j])
      .ai <- sum(abs(.inv_cov[.i, -.i])) 
      .aj <- sum(abs(.inv_cov[.j, -.j])) 
      
      .W[.i, .j] <- (.i!= .j)/(.a*.ai*.aj)
    }
  }
  .W = .W + t(.W)
  
  ##########################################################
  ## Find rho vector:
  .maxrho = maxtuning.hwgl(s = s, p = p, n = n, 
                           penalize.diagonal = penalize.diagonal)
  .rholist = .maxrho * 10^seq(from = range, to = 0, length.out = rholength)
  .BIC = rep(NA, rholength)
  
  .hwgl.temp = NULL
  .start_time = Sys.time()
  for(.i in 1:rholength){
    .hwgl.temp = glasso(s = s, rho = .rholist[.i]*.W, 
                        nobs = n, zero = NULL, 
                        thr = 1.0e-4, maxit = maxit,  approx = FALSE,
                        penalize.diagonal = penalize.diagonal, 
                        start = "cold", w.init = NULL, wi.init = NULL, 
                        trace = FALSE)
    .BIC[.i] = BICemp(s = s, 
                      pm.est = .hwgl.temp$wi, 
                      n = n, 
                      p = p,
                      threshold = bic.threshold)
  }
  .end_time = Sys.time()
  .total_time_HWGL = .end_time-.start_time
  
  ## Obtain optimal model:
  .opt.rho = .rholist[which.min(.BIC)]
  .opt.index = which.min(.BIC)
  .opt.model = glasso(s = s, rho = .opt.rho * .W, 
                      nobs = n, zero = NULL, 
                      thr = 1.0e-4, maxit = maxit,  approx = FALSE,
                      penalize.diagonal = penalize.diagonal, 
                      start = "cold", w.init = NULL, wi.init = NULL, 
                      trace = FALSE)
  
  .OUTPUT = list(opt.model = .opt.model,
                 opt.rho = .opt.rho,
                 opt.index = .opt.index,
                 rholist = .rholist,
                 bic = .BIC)
  
  return(.OUTPUT)
}


##########################################################
##########################################################
## This function finds the smallest tuning parameter for
## which the HWGL estimate is diagonal.

maxtuning.hwgl = function(s = NULL, true.pm = NULL, use.cov = TRUE,
                            p, n, 
                            penalize.diagonal = TRUE) {
  .s = NULL
  if(is.null(s) & !is.null(true.pm)){
    .sigma = solve(true.pm)
    .X = rmvnorm(n = n, 
                 sigma = .sigma,
                 method = "svd")
    if(use.cov) cov(.X) 
    if(!use.cov) cor(.X)
  } else {
    .s = s
  }
  .rho = 1000

  ##############################################
  ##############################################
  ## Find weight matrix: 
  .inv_cov = solve(.s + 0.1 * diag(p) )
  .W  <- matrix(rep(0, p * p), ncol = p)
  for (.i in (2:p)) {
    for(.j in (1:(.i-1))) {
      .a  <- abs(.inv_cov[.i, .j])
      .ai <- sum(abs(.inv_cov[.i, -.i])) 
      .aj <- sum(abs(.inv_cov[.j, -.j])) 
      
      .W[.i, .j] <- (.i!= .j)/(.a*.ai*.aj)
    }
  }
  .W = .W + t(.W)
  
  ##############################################
  ##############################################
  ## Estimate the value for which the model is 
  ## first diagonal.
  .output.hwgl = glasso(s = .s, rho = .rho*.W, 
                        nobs = n, maxit = 200, 
                        penalize.diagonal = penalize.diagonal)
  .is.diagonal = is.diagonal.matrix(.output.hwgl$wi)
  
  .count = 1
  while(!.is.diagonal & .count < 5){
    .count = .count + 1
    .rho = 10^5 * .rho
    .output.hwgl = glasso(s = .s, rho = .rho*.W, 
                          nobs = n, maxit = 200, 
                          penalize.diagonal = penalize.diagonal)
    .is.diagonal = is.diagonal.matrix(.output.hwgl$wi)
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
    .output.hwgl = glasso(s = .s, rho =  .rho*.W, 
                          nobs = n, maxit = 200, 
                          penalize.diagonal = penalize.diagonal)
    .is.diagonal = is.diagonal.matrix(.output.hwgl$wi)
    
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
  
  model = cv.hwgl(s = s, p = 40, n = 200, penalize.diagonal = TRUE, 
                  range = -13, rholength = 100, 
                  maxit = 500, bic.threshold = 10^(-5))
  
  plot(x = log(model$ rholist, base = 10), y = model$bic)
  
  par(mfrow = c(2,1), oma = c(0,0,0,4))
  plot(theta!= 0)
  plot(model$opt.model$wi != 0)
  
  rm(theta, sigma, X, s, model)    
}



rm(examples)

