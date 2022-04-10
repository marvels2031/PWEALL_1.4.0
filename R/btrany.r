#' Back transform the y values to the original scale
#'
#' This function back transforms the y values to the original scale.
#'
#' @param y A vector of transformed y values, y must be in interval [0,1].
#' @param tau A constant denoting the upper bound of the observed data.
#' @param type Type of the transformation, 1: Gamma type, 2: Weibull type, and other values: uniform.
#' @param a Parameter a in the transformation.
#' @param b Parameter b in the transformation.
#' @param c Parameter c in the transformation.
#'
#' @return Returns the back-transformed values in the original scale.
#'
#' @examples
#' y=runif(100)
#' btrany(y=y,tau=10,type=1,a=10,b=2,c=0.5)
#'
#' @export
#'
btrany=function(y,tau=100,type=1,a=tau,b=1,c=0){
  #y must be in interval [0,1]
  if (type==1){temp=y*pgamma(tau+c,shape=b,scale=a);ty=qgamma(temp,shape=b,scale=a)-c}##Gamma
  else if (type==2){ ##Weibull
    tt=(tau+c)/a
    temp=(1-exp(-tt^b))
    ty=a*(-log(1-y*temp))^(1/b)-c
  }
  else {ty=y*(tau+c)-c}
  list(ty=ty)
}
