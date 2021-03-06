% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gmatrix.r
\name{gmatrix}
\alias{gmatrix}
\title{Generate Truncated Spline Base}
\usage{
gmatrix(r = 3, iknots = c(0.2, 0.3), lk = 0, hk = 1)
}
\arguments{
\item{r}{Degree of the splines, must be greater than zero.}

\item{iknots}{A vector of inner knots.}

\item{lk}{lower boundary knot.}

\item{hk}{higher boundary knot.}
}
\value{
The $G$ matrix that links the truncated base (T-base) and the B-base.
}
\description{
Generate Truncated Spline Base
}
\details{
$TG^'=B$, note that the matrix $G$ only depends on the degree \code{r}, 
inner knots \code{iknots} and outer knots \code{lk} and \code{hk}.
}
\examples{
gmatrix(r=3,iknots=c(0.2,0.3),lk=0,hk=1)
x=sort(runif(1000)) # random numbers to ensure no cherry-picking
knots=c(0.04,0.14,0.26,0.38,0.5,0.62,0.74,0.86,0.96)
lk=0;hk=1 #boundary knots
r=3 #order of splines, i.e. cubic splines
###generate B-base using the splines2 package
B2=splines2::bSpline(x,knots=knots, degree=r, intercept=TRUE,Boundary.knots=c(lk,hk))
###generate B-base using splines package
B1=splines::bs(x,knots=knots, degree=r, intercept = TRUE,Boundary.knots=c(lk,hk))
###t-base
TB0=tbase(x,iknots=knots,lk=lk,hk=hk,deg=r)
###compute the linkage matrix GM
abc=gmatrix(r=r,iknots=knots,lk=lk,hk=hk)
###generate B-base from T-base using the linkage function
TB1=TB0\%*\%t(abc$GM)
### Check the differences
c(max(abs(B2-B1)),max(abs(B2-TB1)))

}
