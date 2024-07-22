
#' Variable importance for Distributional Random Forests
#'
#' @param X Matrix with input training data.
#' @param Y Matrix with output training data.
#' @param X_test Matrix with input testing data. If NULL, out-of-bag estimates are used.
#' @param num.trees Number of trees to fit DRF. Default value is 500 trees.
#' @param silent If FALSE, print variable iteration number, otherwise nothing is print. Default is FALSE.
#'
#' @return The list of importance values for all input variables.
#' @export
#'
#' @examples
compute_causaldrf_vimp <- function(X, Y, W, X_test = NULL, num.trees = 500, silent = FALSE){


  X<-as.matrix(X)
  Y<-as.matrix(Y)
  W<-as.matrix(W)

  if (is.null(X_test)){

    X_test<-X

  }

  # Step 1: fit initial DRF to estimate mu(x)=\E[k(Y,.) \mid x]
  bandwidth_Y <- drf:::medianHeuristic(Y)
  k_Y <- rbfdot(sigma = bandwidth_Y)
  K <- matrix(kernelMatrix(k_Y, Y, Y), ncol = nrow(Y))
  DRF <- drf(X, Y, num.trees = num.trees)
  wall <- predict(DRF, newdata=X_test, newtreatment=NULL)$weights
  #wall <- matrix(wall, ncol = ncol(wall))
  
  

  

  # compute normalization constant
  wbar <- colMeans(wall)
  wall_wbar <- sweep(wall, 2, wbar, "-")
  I0 <- as.numeric(sum(diag(wall_wbar %*% K %*% t(wall_wbar))))

  # compute drf importance dropping variables one by one
  I <- sapply(1:ncol(X), function(j) {
    
    # Step 2: Estimate X_j | X_{-j}
    DRFj <- drf(X=X[,-j,drop=F], Y=X[,j,drop=F], num.trees = num.trees)
    wallj <- predict(DRFj, newdata=X_test[,-j,drop=F], newtreatment=NULL)$weights
    
    
   Xjsample<- sapply(1:ncol(Xtest), function(l){
      
      X[sample(1:nrow(Xtest), size = n, replace = T, prob = wallj[l,]),j,drop=F]
      
    })
    
   sapply(1:ncol(Xtest), function(l){
     ## Step 3: Obtain new weights!!
    
     
   })
    
    if (!silent){print(paste0('Running importance for variable X', j, '...'))}
    DRFj <- drf(X = X[, -j, drop=F], Y = Y, W=W, num.trees = num.trees)


    wj <- predict(DRFj, X_test[, -j])$weights


    wj <- matrix(wj, ncol = ncol(wj))
    Ij <- sum(diag((wj - wall) %*% K %*% t(wj - wall)))/I0
    return(Ij)
  })

  # compute retraining bias
  DRF0 <- drf(X = X, Y = Y, W=W, num.trees = num.trees)
  w0 = predict(DRF0, X_test)$weights
  w0 <- matrix(w0, ncol = ncol(w0))
  vimp0 <- sum(diag((w0 - wall) %*% K %*% t(w0 - wall)))/I0

  # compute final importance (remove bias & truncate negative values)
  vimp <- sapply(I - vimp0, function(x){max(0,x)})

  names(vimp)<-colnames(X)

  return(vimp)

}
