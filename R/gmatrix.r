#' Generate Truncated Spline Base
#' 
#' @param r Degree of the splines, must be greater than zero.
#' @param iknots A vector of inner knots.
#' @param lk lower boundary knot.
#' @param hk higher boundary knot.
#' @return The $G$ matrix that links the truncated base (T-base) and the B-base.
#' @details $TG^'=B$, note that the matrix $G$ only depends on the degree \code{r}, 
#' inner knots \code{iknots} and outer knots \code{lk} and \code{hk}.
#' @examples
#' gmatrix(r=3,iknots=c(0.2,0.3),lk=0,hk=1)
#' x=sort(runif(1000)) # random numbers to ensure no cherry-picking
#' knots=c(0.04,0.14,0.26,0.38,0.5,0.62,0.74,0.86,0.96)
#' lk=0;hk=1 #boundary knots
#' r=3 #order of splines, i.e. cubic splines
#' ###generate B-base using the splines2 package
#' B2=splines2::bSpline(x,knots=knots, degree=r, intercept=TRUE,Boundary.knots=c(lk,hk))
#' ###generate B-base using splines package
#' B1=splines::bs(x,knots=knots, degree=r, intercept = TRUE,Boundary.knots=c(lk,hk))
#' ###t-base
#' TB0=tbase(x,iknots=knots,lk=lk,hk=hk,deg=r)
#' ###compute the linkage matrix GM
#' abc=gmatrix(r=r,iknots=knots,lk=lk,hk=hk)
#' ###generate B-base from T-base using the linkage function
#' TB1=TB0%*%t(abc$GM)
#' ### Check the differences
#' c(max(abs(B2-B1)),max(abs(B2-TB1)))
#' 
#' @export
#' 
gmatrix=function(r=3,iknots=c(0.2,0.3),lk=0,hk=1){
  ###r must be greater than zero
  ###lk; ###lower boundary knot
  ###hk; ###upper boundary knot, PLEASE make sure no point is exceeding it
  r1=r+1
  
  lm2=length(iknots) ###number of inner knots 
  
  lx=lm2+2 ### number of knots (including two boundary knots)
  p=lx+r-1 ### number of parameters for both B-splines and T-splines
  tall=c(rep(lk,r),c(lk,iknots,hk),rep(hk,r)) ####repeat r times of both the boundary knots
  l2r=length(tall) ###total number of knots in action is l+2*r
  
  t1=tall[1+r];t2=tall[2+r];t3=tall[3+r]
  H2=matrix(c(1,-2/(t2-t1),0,2/(t2-t1)),nrow=2,byrow=T)
  H3=matrix(c(1,-3/(t2-t1),3/(t2-t1)^2,0,3/(t2-t1),-3/(t2-t1)/(t3-t1)-3/(t2-t1)^2,0,0,3/(t2-t1)/(t3-t1)),nrow=3,byrow=T)
  
  G0=matrix(0,nrow=l2r,ncol=lx)
  G0[1:r,1]=1
  G0[(lx+r+1):l2r,lx]=1
  G0[(1+r):(lx+r),]=diag(lx)
  j=1
  while (j<=r){
    la=l2r-j+1
    gd=tall[-(1:j)]-tall[-(la:l2r)]
    ii=(gd>0)
    ngd=length(gd)
    ww=rep(0,ngd)
    ww[ii==1]=1/gd[ii==1];#length(ww)
    G1=diag(ww)%*%diff(diag(la),diff=1)%*%G0;dim(G1)
    G0=G1
    j=j+1
  }
  G0=(-1)^(r+1)*diff(diag(la-1),diff=1)%*%G0
  if (r==0){GM=G0[,-lx]}
  else if (r==1){GM=cbind(rep(0,p),G0[,-lx]);GM[1:r,1:r]=1}
  else if (r==2){GM=cbind(rep(0,p),rep(0,p),G0[,-lx]);GM[1:r,1:r]=H2}
  else if (r==3){GM=cbind(rep(0,p),rep(0,p),rep(0,p),G0[,-lx]);GM[1:r,1:r]=H3}
  list(GM=GM,G0=G0,H2=H2,H3=H3)
}
