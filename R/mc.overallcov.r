#' Calculate the covariance for NPH MCP-Mod
#'
#' This function combines functions "mc.overallcovp1" and "mc.overallcovp1" to calculate the covariance of max-combo tests.
#' 
#' @param wfunctions A vector of the indexes of selected weight functions for calculation made at time \code{tfix}.
#' @param wfunctions0 A vector of the indexes of selected weight functions for calculation made at time \code{tfix0}.
#' @param tfix A constant for the time for weighted log-rank using weights \code{wfunctions}.
#' @param tfix0 A constant for the time for weighted log-rank using weights \code{wfunctions0}.
#' @param taur A constant denoting the recruitment period.
#' @param u A vector of recruitment rates.
#' @param ut A vector of time-points when recruitment rate changes. 
#' @param pi1 A constant denoting the proportion of patients randomized to treatment arm. 
#' @param rate11 Hazard before crossover for the treatment group.
#' @param rate21 Hazard after crossover for the treatment group.
#' @param rate31 Hazard for time to crossover for the treatment group.
#' @param rate41 Hazard after crossover for the treatment group for complex case.
#' @param rate51 Hazard after crossover for the treatment group for complex case.
#' @param ratec1 Hazard for time to censoring for the treatment group.
#' @param rate10 Hazard before crossover for the control group.
#' @param rate20 Hazard after crossover for the control group.
#' @param rate30 Hazard for time to crossover for the control group.
#' @param rate40 Hazard after crossover for the control group for complex case.
#' @param rate50 Hazard after crossover for the control group for complex case.
#' @param ratec0 Hazard for time to censoring for the control group.
#' @param tchange A strictly increasing sequence of time points at which the event rates changes. 
#' The first element of tchange must be zero. It must have the same length as \code{rate11}, \code{rate21}, \code{rate31}, etc.
#' @param type1 Type of crossover in the treatment group.
#' @param type0 Type of crossover in the control group.
#' @param rp21 re-randomization prob for the treatment group.
#' @param rp20 re-randomization prob for the control group.
#' @param eps A small number representing the error tolerance when calculating the utility function 
#'  \deqn{\Phi_l(x)=\frac{\int_0^x s^l e^{-s}ds}{x^{l+1}}} with \eqn{l=0,1,2}.
#' @param veps A small number representing the error tolerance when calculating the Fisher information.
#' @param beta The value at which the covaraince is computed.
#' @param beta0 The value at which the covaraince is computed.
#'
#' @return Returns the covariance at the designated time-points \code{tfix} and \code{tfix0}. 
#' 
#' Extremely important! Please make sure "mc.weightfuns" have been defined in the global environment. The following provides a simple example. 
#' Please refer to the vignette file for more detail.  
#' 
#' @examples
#' #Define 'mc.weightfuns'
#' n.weights=4
#' degree=3
#' inner.knots=c(0.1,0.2,0.3,0.5,0.7)
#' boundary.knots=c(0,1)
#' np=degree+1+length(inner.knots)
#' btmatrix=matrix(0,nrow=n.weights,ncol=np)
#' btmatrix[1,1]=1
#' btmatrix[2,1:3]=c(3,-2,0)/6
#' btmatrix[3,1:3]=c(3,2,-2)/4
#' btmatrix[4,1:3]=c(1,-1,1)
#' mc.weightfuns <- vector("list", n.weights)
#' for (i in 1:n.weights) {
#'   mc.weightfuns [[i]] <- mc.fun1(degree=3,inner.knots=inner.knots,boundary.knots=boundary.knots,
#'   bt=btmatrix[i,],base='T',type=-1,tau=6)
#' }
#' #Calculate the covariances 
#' mc.overallcov(wfunctions=c(1,2),wfunctions0=c(1,2),tfix = 2, tfix0 = 1)
#' 
#' @export
#'
mc.overallcov=function(wfunctions=c(1,2),wfunctions0=c(1,2),tfix = 2, tfix0 = 1, 
                       taur = 5, u = c(1/taur, 1/taur), ut = c(taur/2, taur), pi1 = 0.5, 
                       rate11 = c(1, 0.5), rate21 = rate11, rate31 = c(0.7, 0.4), 
                       rate41 = rate21, rate51 = rate21, ratec1 = c(0.5, 0.6), 
                       rate10 = rate11, rate20 = rate10, rate30 = rate31, 
                       rate40 = rate20, rate50 = rate20, ratec0 = ratec1, 
                       tchange = c(0, 1), type1 = 1, type0 = 1, rp21 = 0.5, rp20 = 0.5, eps = 0.01, 
                       veps = 0.01, beta = 0, beta0 = 0) 
{
  part1 <- mc.overallcovp1(wfunctions=wfunctions,wfunctions0=wfunctions0,tfix = tfix, tfix0 = tfix0, taur = taur, 
                           u = u, ut = ut, pi1 = pi1, rate11 = rate11, rate21 = rate21, 
                           rate31 = rate31, rate41 = rate41, rate51 = rate51, ratec1 = ratec1, 
                           rate10 = rate10, rate20 = rate20, rate30 = rate30, rate40 = rate40, 
                           rate50 = rate50, ratec0 = ratec0, tchange = tchange, 
                           type1 = type1, type0 = type0, rp21 = rp21, rp20 = rp20, 
                           eps = eps, veps = veps, beta = beta, beta0 = beta0)
  part2 <- mc.overallcovp2(wfunctions=wfunctions,wfunctions0=wfunctions0,tfix = tfix, tfix0 = tfix0, taur = taur, 
                           u = u, ut = ut, pi1 = pi1, rate11 = rate11, rate21 = rate21, 
                           rate31 = rate31, rate41 = rate41, rate51 = rate51, ratec1 = ratec1, 
                           rate10 = rate10, rate20 = rate20, rate30 = rate30, rate40 = rate40, 
                           rate50 = rate50, ratec0 = ratec0, tchange = tchange, 
                           type1 = type1, type0 = type0, rp21 = rp21, rp20 = rp20, 
                           eps = eps, veps = veps, beta = beta, beta0 = beta0)
  EA1 <- part1$EA1
  EA2 <- part2$EA2
  covbeta <- part1$covbeta1 + part2$cov234-EA1%*%t(EA2)
  covbeta1 <- part1$covbeta1
  covbeta2 <- part2$covbeta2
  covbeta3 <- part2$covbeta3
  covbeta4 <- part2$covbeta4
  list(covbeta = covbeta, covbeta1 = covbeta1, covbeta2 = covbeta2, 
       covbeta3 = covbeta3, covbeta4 = covbeta4, EA1 = EA1, 
       EA2 = EA2)
}
