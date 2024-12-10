## Moran's Basis
## Copied from https://github.com/JonathanBradley28/CM/blob/master/R/MoransIBasis.R

#' Moran's I Basis Function
#'
#' This code computes the Moran's I basis functions
#' @param X An nxp matrix of covariates. Each column represents a covariate. Each row corresponds to one of the n replicates.
#' @param r The number of basis functions.
#' @param A An nxn adjacency matrix.
#' @return Psi nxr matrix of basis functions
#' @export
MoransI.Basis<-function(X,r,A){


  n = dim(X)[1]

  PsiOper = (diag(n) - X%*%solve(t(X)%*%X)%*%t(X))%*%A%*%(diag(n) - X%*%solve(t(X)%*%X)%*%t(X))
  output2<-eigen(PsiOper)
  Psi = output2$vectors[,1:r]

  return(Psi)
}


## BSpline Basis
## Based on R package "splines2": https://www.rdocumentation.org/packages/splines2/versions/0.2.5/topics/bSpline
#' B-spline Basis Function
#'
#' This code computes the B-spline basis functions
#' @param x A length n vector of predictor values.
#' @param knots A vector of internal knots. E.g., seq(0,1,length.out=r-degree+1)[-c(1,r-degree+1)]
#' @param boundary A vector of boundary knots. E.g., c(0,1)
#' @param J The number of basis functions.
#' @param degree The degress of basis functions. The default value is 3 for cubic splines. 
#' @return Psi nxr matrix of basis functions
#' @export
Bspline.Basis<-function(x, knots = seq(0,1,length.out=8)[-c(1,8)], boundary = c(0,1),r = 10, degree = 3){

  Phi = bSpline(x = x ,knots = knots, df = r, degree = degree, Boundary.knots = boundary, intercept = TRUE)

  return(Phi)
}
