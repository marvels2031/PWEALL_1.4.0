#' Define weight function
#'
#' This function defines and returns weight function..
#'
#' @param degree The degree of splines.
#' @param inner.knots The inner knots.
#' @param boundary.knots The boundary knots.
#' @param bt A vector of coefficients for the splines.
#' @param base A character determining which splines are used, T: T-splines, B: B-splines.
#' @param tau A constant denoting the upper bound of the observed data.
#' @param type Type of the transformation, 1: Gamma type, 2: Weibull type, and other values: uniform.
#' @param a Parameter a in the transformation.
#' @param b Parameter b in the transformation.
#' @param c Parameter c in the transformation.
#'
#' @return Returns a weight function. 
#'
#' @examples
#' n.weights=4 #will define 4 weight functions
#' degree=3    #with degree of 3
#' inner.knots=c(1,2,2.5) #inner.knots
#' boundary.knots=c(0,6)  #boundary knots
#' np=degree+1+length(inner.knots) #total number of parameters
#' btmatrix=matrix(0,nrow=n.weights,ncol=np) #coefficient matrix
#' btmatrix[1,]=1
#' btmatrix[2,]=rnorm(np)
#' btmatrix[3,]=rnorm(np)
#' btmatrix[4,]=rnorm(np)
#' ## Define the 4 weight functions
#' weightfuns <- vector("list", n.weights)
#' for (i in 1:n.weights) {
#'  weightfuns[[i]] <- fun1(degree=3,inner.knots=c(1,2,2.5),boundary.knots=c(0,6),bt=btmatrix[i,],base='B')
#' }
#'
#' @export
#'
fun1 <- function(degree=3,inner.knots=c(1,2,2.5),boundary.knots=c(0,6),bt=rep(0,1+degree+length(inner.knots)),base='T',tau=1,type=-1,a=tau,b=1,c=0) {
  force(degree)
  force(inner.knots)
  force(boundary.knots)
  force(bt)
  force(base)
  force(tau)
  force(type)
  force(a)
  force(b)
  force(c)
  fun2 <- function(x) {
    tx=trany(x,tau=tau,type=type,a=a,b=b,c=c)$ty
    txinner.knots=trany(inner.knots,tau=tau,type=type,a=a,b=b,c=c)$ty
    txboundary.knots=trany(boundary.knots,tau=tau,type=type,a=a,b=b,c=c)$ty
    if (base=='T'){
      XB=tbase(tx,iknots=txinner.knots,lk=txboundary.knots[1],hk=txboundary.knots[2],deg=degree)
    }
    else if (base=='B'){
      XB=splines2::bSpline(tx,knots=txinner.knots, degree=degree,intercept=TRUE,Boundary.knots=txboundary.knots)
    }
    XB%*%bt
  }
  return(fun2)
}