#' Calculate the information function for NPH MCP-Mod
#'
#' This function is to calculate the information function for the max-combo test statistics.
#' 
#' @param wfunctions A vector of the indexes of selected weight functions for calculation made at time \code{tfix}.
#' @param tfix A constant for the time for weighted log-rank using weights \code{wfunctions}.
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
#' @param beta The value at which the covariance is computed.
#'
#' @return Returns two integrations at the designated time-points \code{tfix}. 
#' 
#' @details
#' This function is to calculate as qf1[j]
#' \deqn{\int_0^{t} w^2_j(s)G_E(t-s)q_1(s)q_0(s)\zeta(s)ds;}
#' calculate as EA1[j]
#' \deqn{\int_0^{t} w_j(s)G_E(t-s)a(s)ds.}
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
#'   mc.weightfuns[[i]] <- mc.fun1(degree=3,inner.knots=inner.knots,
#'                               boundary.knots=boundary.knots,
#'                               bt=btmatrix[i,],base='T',type=-1,tau=6)
#' }
#' #Calculate the covariances 
#' mc.inf(wfunctions=c(1,2),tfix = 2)
#' 
#' @export
#'
mc.inf=function(wfunctions=c(1,2),tfix = 2, taur = 5, u = c(1/taur, 1/taur), ut = c(taur/2,taur), pi1 = 0.5, 
                rate11 = c(1, 0.5), rate21 = rate11, rate31 = c(0.7, 0.4), 
                rate41 = rate21, rate51 = rate21, ratec1 = c(0.5, 0.6), 
                rate10 = rate11, rate20 = rate10, rate30 = rate31, rate40 = rate20, 
                rate50 = rate20, ratec0 = c(0.6, 0.5), tchange = c(0, 1), 
                type1 = 1, type0 = 1, rp21 = 0.5, rp20 = 0.5, eps = 0.01, 
                veps = 0.01, beta = 0) 
{
  ratemax <- max(abs(rate11 - rate10)) + max(abs(rate21 - rate20)) + 
    max(abs(rate31 - rate30)) + max(abs(rate41 - rate40)) + 
    max(abs(rate51 - rate50)) + max(abs(ratec1 - ratec0))
  rateb <- max(0.01, min(ratemax, 1))
  err <- veps/rateb
  tmax <- max(c(tfix, tchange, taur)) + err
  nr <- length(rate11)
  tplus <- rep(0, nr)
  tplus[nr] <- tmax
  if (nr > 1) 
    tplus[-nr] <- tchange[-1]
  nn <- rep(1, nr)
  nn[1] <- ceiling((tplus[1] - tchange[1])/err)
  atchange <- rep(0, nn[1])
  atchange <- seq(tchange[1], tplus[1], by = (tplus[1] - tchange[1])/nn[1])[1:nn[1]]
  if (nr >= 2) {
    for (i in 2:nr) {
      nn[i] <- ceiling((tplus[i] - tchange[i])/err)
      atchange <- c(atchange, seq(tchange[i], tplus[i], 
                                  by = (tplus[i] - tchange[i])/nn[i])[1:nn[i]])
    }
  }
  atchange1 <- sort(unique(c(atchange, tfix), fromLast = T))
  anr <- length(atchange1) + 1
  atplus <- rep(0, anr)
  atplus[anr] <- tmax
  atplus[-anr] <- atchange1
  ats <- atplus[atplus < (tfix - 0.1 * err)]
  atsupp <- c(ats, tfix)
  nsupp <- length(atsupp)
  BigK1 <- pwecxpwuforvar(tfix = tfix, t = atsupp, taur = taur, 
                          u = u, ut = ut, rate1 = rate11, rate2 = rate21, rate3 = rate31, 
                          rate4 = rate41, rate5 = rate51, ratec = ratec1, tchange = tchange, 
                          type = type1, rp2 = rp21, eps = eps)
  BigK0 <- pwecxpwuforvar(tfix = tfix, t = atsupp, taur = taur, 
                          u = u, ut = ut, rate1 = rate10, rate2 = rate20, rate3 = rate30, 
                          rate4 = rate40, rate5 = rate50, ratec = ratec0, tchange = tchange, 
                          type = type0, rp2 = rp20, eps = eps)
  dk1 <- BigK1$f0[-1] - BigK1$f0[-nsupp]
  dk0 <- BigK0$f0[-1] - BigK0$f0[-nsupp]
  adk1 <- (dk1 > 1e-08)
  adk0 <- (dk0 > 1e-08)
  tk1 <- tk0 <- atsupp[-nsupp]
  tk1[adk1 == 1] <- (BigK1$f1[-1] - BigK1$f1[-nsupp])[adk1 == 
                                                        1]/dk1[adk1 == 1]
  tk0[adk0 == 1] <- (BigK0$f1[-1] - BigK0$f1[-nsupp])[adk0 == 
                                                        1]/dk0[adk0 == 1]
  ST11 <- pwecx(t = tk1, rate1 = rate11, rate2 = rate21, rate3 = rate31, 
                rate4 = rate41, rate5 = rate51, tchange = tchange, type = type1, 
                rp2 = rp21, eps = eps)$surv
  ST10 <- pwecx(t = tk1, rate1 = rate10, rate2 = rate20, rate3 = rate30, 
                rate4 = rate40, rate5 = rate50, tchange = tchange, type = type0, 
                rp2 = rp20, eps = eps)$surv
  SC11 <- pwe(t = tk1, rate = ratec1, tchange = tchange)$surv
  SC10 <- pwe(t = tk1, rate = ratec0, tchange = tchange)$surv
  ST01 <- pwecx(t = tk0, rate1 = rate11, rate2 = rate21, rate3 = rate31, 
                rate4 = rate41, rate5 = rate51, tchange = tchange, type = type1, 
                rp2 = rp21, eps = eps)$surv
  ST00 <- pwecx(t = tk0, rate1 = rate10, rate2 = rate20, rate3 = rate30, 
                rate4 = rate40, rate5 = rate50, tchange = tchange, type = type0, 
                rp2 = rp20, eps = eps)$surv
  SC01 <- pwe(t = tk0, rate = ratec1, tchange = tchange)$surv
  SC00 <- pwe(t = tk0, rate = ratec0, tchange = tchange)$surv
  bb1 <- (1 - pi1) * ST10 * SC10
  bb0 <- (1 - pi1) * ST00 * SC00
  aa1 <- pi1 * exp(beta) * ST11 * SC11
  aa0 <- pi1 * exp(beta) * ST01 * SC01
  q1bs <- aa1/(aa1 + bb1)
  q0bs <- aa0/(aa0 + bb0)
  
  nw=length(wfunctions)
  qf1=EA1=rep(0,nw)
  for (i in 1:nw){
    ai=wfunctions[i]
    wf0i=mc.weightfuns[[ai]](tk0)
    wf1i=mc.weightfuns[[ai]](tk1)
    qf1[i] <- pi1 * sum(q1bs * (1 - q1bs) * dk1*wf1i^2) + (1 - pi1) * 
      sum(q0bs * (1 - q0bs) * dk0*wf0i^2)
    EA1[i] <- pi1 * sum((1 - q1bs) * dk1*wf1i) - (1 - pi1) * 
      sum(q0bs * dk0*wf0i)
  }
  list(qf1 = qf1,EA1 = EA1)
}