#################################################
#################################################
##
## In the following document, we introduce the
## procedure for variable screening for multiple
## covariance matrices.
##
#################################################
#################################################

######################
## Reduce dimension by choosing larger max-corr.
.reduce.dim = function (empcovs, p, t, K, n, method = "maxcor") {
  
  ## Measure variable strength:
  .maxcorK = matrix(0, ncol = K, nrow = p)
  for(.k in 1:K){
    if (method == "maxcor") {
      .maxcorK[, .k] = apply(empcov - diag(diag(empcov)), 
                            MARGIN = 1,
                            function(x) max(abs(x)))
    } else if (method == "l2") {
      .maxcorK[, .k] = apply(empcov - diag(diag(empcov)), 
                            MARGIN = 1,
                            function(x) sum(x^2))
    } else if (method == "l1") {
      .maxcorK[, .k] = apply(empcov - diag(diag(empcov)), 
                            MARGIN = 1,
                            function(x) sum(abs(x)))
    }    
  }
  .maxcor = apply(X = .maxcorK,
                  MARGIN = 1,
                  max)
  
  ## Choose the variables with top values
  .order = order(.maxcor, decreasing = TRUE)
  
  ## TRUE on selected vars/ 
  ## FALSE on dropped vars.
  .vars = (1:p) %in% .order[1:t]
  
  ## New-index-num on selected vars/ 
  ## zero on dropped vars.
  .position = rep(0, p)
  for(.i in 1:p){
    .position[.i] = sum(.vars[1:(.i)])* .vars[.i]
  }
  
  ## Reduced matrices.
  .red.empcov = list()
  for(.k in 1:K){
    .red.empcov[[.k]] = empcov[.vars, .vars]
  }
  
  ## Produce outcomes:
  .OUTPUT = list(red.empcov = .red.empcov,
                 vars = .vars,
                 strength = .maxcor, 
                 position = .position)
  return(.OUTPUT)
}
