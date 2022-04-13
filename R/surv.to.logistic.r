#' Transform the two-sample survival data into logistic regression data
#' 
#' The following function transforms the two-sample survival data into a logistic regression model with an offset. 
#' 
#' @param y A vector of observed times.
#' @param d A vector of event indicators with 1 denoting event and 0 denoting censoring.
#' @param z A vector of group indicators with 1 being the treatment and 0 the control.
#' @param tcut A constant as the cut-point.
#' @return A dataframe for logistic regression with z denoting 1/0 outcome, y the (univariate) covariate and yset the offset.
#' @examples
#' n1=100;n0=100;n=n1+n0
#' zz=rep(0,n);zz[1:n1]=1
#' tt=c(rexp(n1)/0.1,rexp(n0)/0.2)
#' cc=c(rexp(n1)/0.05,rexp(n0)/0.06)
#' yy=pmin(tt,cc)
#' dd=(tt<=yy)
#' df=surv.to.logistic(y=yy,d=dd,z=zz,tcut=5)$adata
#'  
#' @export
#' 
surv.to.logistic=function(y,d,z,tcut=max(y[d==1])){
  abc=PWEALL::wlrutil(y=y,d=d,z=z,te=y)
  aindex=as.numeric(abc$mfunc[,4]>0&abc$mfunc[,5]>0)
  atemp=rep(0,length(y))
  atemp[aindex==1]=log(abc$mfunc[aindex==1,4]/abc$mfunc[aindex==1,5])
  bindex=as.numeric(aindex==1&d==1&y<=tcut)
  ax=cbind(y,d,z,atemp,abc$mfunc[,4],abc$mfunc[,5],abc$mfunc[,6])[bindex==1,]
  adata=data.frame(z=ax[,3],y=ax[,1],yset=ax[,4])
  list(adata=adata)
}
