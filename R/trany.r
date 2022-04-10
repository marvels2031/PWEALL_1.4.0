#' Transform the observed survival times into [0,1]
#'
#' This function transforms the observed survival times into [0,1].
#'
#' @param y A vector of observed times to be transformed.
#' @param tau A constant denoting the upper bound of the observed data.
#' @param type Type of the transformation, 1: Gamma type, 2: Weibull type, and other values: uniform.
#' @param a Parameter a in the transformation.
#' @param b Parameter b in the transformation.
#' @param c Parameter c in the transformation.
#'
#' @return Returns the transformed values in [0,1].
#'
#' @examples
#' y=rexp(100)
#' trany(y=y,tau=10,type=1,a=10,b=2,c=0.5)
#'
#' @export
#'
trany=function(y,tau=max(y)+0.1,type=1,a=tau,b=1,c=0){
  if (type==1){ty=pgamma(y+c,shape=b,scale=a)/pgamma(tau+c,shape=b,scale=a)}##Gamma
  else if (type==2){ytemp=(y+c)/a;temp=(tau+c)/a;ty=(1-exp(-ytemp^b))/((1-exp(-temp^b)))}##Weibull
  else {ty=(y+c)/(tau+c)}
  list(ty=ty)
}
