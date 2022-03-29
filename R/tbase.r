#' Generate Truncated Spline Base
#' 
#' @param x A vector.
#' @param iknots A vector of inner knots.
#' @param lk lower boundary knot.
#' @param hk higher boundary knot.
#' @param deg degree of the splines, default is 3.
#' @return The truncated spline base.
#' @examples
#' tbase(x=c(0.25,0.35,0.5,0.7),iknots=c(0.2,0.4,0.6,0.8),lk=0,hk=1,deg=3)
#' 
#' @export
#' 
tbase=function(x,iknots,lk,hk,deg=3){
  knots=c(lk,iknots)
  B=cbind(outer(x-lk,0:(deg-1),"^"),
          outer(x,knots,function(x,y)pmax(x-y,0)^deg))
  B
}


