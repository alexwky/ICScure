Gpx<-function(x,rho){
  ifelse(rho==0,(1+x)^(-1),exp(-rho^(-1)*((1+x)^rho-1)))
}

#' Evaluate a user-provided Bernstein polynomial
#'
#' @description  Evaluate the Bernstein polynomial at given time points using user-provided coefficients.
#'
#' @param t Values at which the Bernstein polynomial is evaluated
#' @param psi Vector of coefficients for the Bernstein polynomial (not including the intercept). The dimension of the vector is equal to the degree of the polynomial. The intercept is assumed to be 0
#' @param mint Lower boundary of the support of the polynomial
#' @param maxt Upper boundary of the support of the polynomial
#'
#' @details The vector \code{psi} contains the coefficients of the Bernstein basis polynomials. The polynomial is supported on (\code{mint} , \code{maxt}) and is equal to zero at \code{mint}. The intercept should not be included in \code{psi}. The function returns the values of the polynomial at the time points in \code{t}.
#' @examples
#' t <- seq(1,10,1)
#' psi <- seq(0.25,1,0.25)
#' BernsteinPolynomial(t = t, psi = psi, mint = 0, maxt = 5)
#' #0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0
#'
#' @return A vector of values of the Bernstein polynomial at \code{t}.
#' @export

BernsteinPolynomial<-function(t,psi,mint,maxt){
  degree <- length(psi)
  t.scaled<-(t-mint)/(maxt-mint)
  res <- rep(NA,length(t))
  for (j in 1:length(t)){
    res[j] <- sum(psi*sapply(1:(degree),function(i){choose((degree),i)*t.scaled[j]^i*(1-t.scaled[j])^(degree-i)}))
  }
  return(res)
  #t.scaled<-(t-mint)/(maxt-mint)
  #return(
  #  sum(psi*sapply(1:(degree),function(i){choose((degree),i)*t.scaled^i*(1-t.scaled)^(degree-i)}))
  #)
}


LOG.Lik.csw<-function(par,rho,degree,data,x0num,t.min,t.max){
  nn<-nrow(data)
  id<-c(1:nn)
  wi<-data[,2]
  Ci<-data[,3]
  cen<-data[,4]
  if(ncol(data)<5){
    xij<-rep(1,nrow(data))
  }else{
    xij<-cbind(1,data[,5:ncol(data)])
  }

  lambda <- par[1:degree]
  phi <- exp(lambda)
  psi <-sapply(c(1:degree),function(a){sum(phi[1:a])})/sum(phi)  #x.scaled<-(x-min(x))/(max(x)-min(x))
  beta <- par[(degree+1):(degree+x0num)]

  mu1<-exp(apply(xij*matrix(rep(beta,nn),nrow=nn,byrow=T),1,sum))

  SiC<-ifelse(Ci<t.max,sapply(1:nn, function(o){Gpx(BernsteinPolynomial(t=Ci[o], psi = psi, mint=t.min, maxt=t.max)*mu1[o],rho)}),
                        sapply(1:nn, function(o){Gpx(sum(exp(lambda))*mu1[o],rho)}))

  Sum.Log.Lik<-0
  for(i in c(1:nn)){
    Log.Lik.attribute <- wi[i]^(-1)*(cen[i]*log(1-SiC[i])+(1-cen[i])*log(SiC[i]))
    Sum.Log.Lik<-Sum.Log.Lik+Log.Lik.attribute
  }
  -(Sum.Log.Lik)
}




LOG.Lik.IJ.csw<-function(par,rho,degree,datline,x0num,t.min,t.max){ #data is in matrix format

  wi<-datline[2]
  Ci<-datline[3]
  cen<-datline[4]
  if(length(datline)<5){
    xij<- 1
  }else{
    xij<-c(1,datline[5:length(datline)])
  }


  lambda <- par[1:degree]
  phi <- exp(lambda)
  psi <-sapply(c(1:degree),function(a){sum(phi[1:a])})/sum(phi) #x.scaled<-(x-min(x))/(max(x)-min(x))
  beta <- par[(degree+1):(degree+x0num)]

  mu1<-exp(sum(xij*beta))

  SiC<-ifelse(Ci<t.max,Gpx(BernsteinPolynomial(t=Ci, psi = psi, mint=t.min, maxt=t.max)*mu1,rho),
                       Gpx(sum(exp(lambda))*mu1,rho))

  cen*log(1-SiC)+(1-cen)*log(SiC)
}

#' Data simulator
#'
#' @description Generate a data set based on the simulation setting in Lam et al. (2021) <DOI: \href{https://doi.org/10.1002/sim.8910}{10.1002/sim.8910}>
#'
#' @param seed Seed of the random generator (optional)
#' @param n Cluster index
#' @param beta00 True parameter value of the first component of \eqn{\beta}
#' @param beta10 True parameter value of the second component of \eqn{\beta}
#' @param beta20 True parameter value of the third component of \eqn{\beta}
#' @param rho Transformation parameter \eqn{\rho} in the Box-Cox transformation
#' @param gamma Parameter \eqn{\gamma} in the Gumbel Copula
#' @param cs Maximum cluster size
#'
#' @details See Section 4 of Lam et al. (2021) <DOI: \href{https://doi.org/10.1002/sim.8910}{10.1002/sim.8910}> for the data generation mechanism.
#'
#' @return A data frame consisting of the following columns: \itemize{
#'   \item \strong{family} : Cluster index
#'   \item \strong{Ci} : Observation time
#'   \item \strong{delta} : Event indicator
#'   \item \strong{x1} : 1st covariate
#'   \item \strong{x2} : 2nd covariate
#' }
#'
#' @examples data <- ICDASim(seed = 1942, n = 100, beta00 = 0.5, beta10 = -1, beta20 = 1,
#' cs = 40, rho = 0.5,gamma = 0.5)
#' ## Generate 100 clusters with maximum size 40
#' @export


ICDASim<-function(seed=NA,n,beta00,beta10,beta20,rho,gamma,cs){

  if(!is.na(seed)){
    set.seed(seed)
  }
  iu <- complex(real=0, imaginary=1)
  xi<-rstable(n = n*3, alpha = gamma, beta = 1, gamma = abs(1-iu*tan(pi*gamma/2))^(-1/gamma), delta = 0, pm = 1)

  clustersize<-ifelse(xi<=qstable(p = 0.5,alpha=gamma,beta = 1,gamma = abs(1-iu*tan(pi*gamma/2))^(-1/gamma),delta = 0,pm = 1),
                      rbinom(n*3,cs,0.75),rbinom(n*3,cs,(1-0.75)))

  cs.reduced<-clustersize[which((clustersize>1)&(clustersize<cs))][1:n]
  xi.reduced<-xi[which((clustersize>1)&(clustersize<cs))][1:n]

  data<-NULL
  for(i in c(1:n)){

    Sij<- exp(-(-log(runif(n = cs.reduced[i],0,1))/xi.reduced[i])^gamma)

    x1<-rbinom(n = cs.reduced[i],size = 1,prob = 0.5)
    x2<-rnorm(n = cs.reduced[i],mean = 0,sd = 1)

    etaij <-exp(beta00+beta10*x1+beta20*x2)

    if(rho==0){
      cureij<-(1+etaij)^(-1)
    }else if(rho>0){
      cureij<-exp(-rho^(-1)*((1+etaij)^rho-1))
    }

    if(rho==0){
      t1<-ifelse(Sij<cureij, 999, -log(1-(Sij^(-1)-1)*(1-exp(-4))/etaij))
    }else if(rho>0){
      t1<-ifelse(Sij<cureij, 999, -log(1-((-log(Sij)*rho+1)^(1/rho)-1)*(1-exp(-4))/etaij))
    }

    Ci<-pmin(runif(n = cs.reduced[i],min = 0,max = 8), 4)

    delta<-ifelse(t1<=Ci,1,0)

    data<-rbind(data,cbind(i,Ci,delta,x1,x2))
  }
  data<-data.frame(data)
  names(data)<-c("family","Ci","delta","x1","x2")
  return(data)
}

#' Estimate model parameters and standard errors
#'
#' @description Perform the cluster-weighted GEE or GEE estimation in Lam et al. (2021) <DOI: \href{https://doi.org/10.1002/sim.8910}{10.1002/sim.8910}>
#'
#' @param data An \code{n x (p+3)} matrix, where \code{n} is the sample size and \code{p} is the number of covariates. The first column of data consists of cluster indices, the second column consists of the observation time, the third column consists of the event indicator, and the fourth to the last columns consist of the covariates (not including the intercept). The set of covariates can be empty.
#' @param rho Transformation parameter \eqn{\rho} in the Box-Cox transformation
#' @param degree Degree of the Bernstein polynomial
#' @param weighting Logical; If TRUE, then CWGEE is adopted; otherwise, GEE is adopted. Default is TRUE
#' @param t.min Left end point of the support of the distribution function \emph{F}
#' @param t.max Cure threshold. If \code{NA}, then it is set to be the largest observation with \eqn{\Delta = 1}
#' @param reltol Relative tolerance for the function \code{optim} (under the package \code{stat}) used for the maximization of the pseudo likelihood function
#' @param maxit Maximum number of iterations for \code{optim}
#' @param calSE Logical; If TRUE, then the standard error of the regression parameter estimators would be estimated
#'
#' @return A list of the following components: \itemize{
#'   \item \strong{degree} : Degree of the Bernstein polynomial for the estimated distribution function \emph{F}
#'   \item \strong{psi} : Vector of the coefficients of the Bernstein basis polynomial for the estimated distribution function \emph{F}
#'   \item \strong{beta} : Vector of the regression parameters \eqn{\beta}
#'   \item \strong{betaSE} : Estimated standard error of the estimator for \eqn{\beta}; NA if calSE = FALSE
#'   \item \strong{iteration} : Numbers of iteration used in \code{optim}
#'   \item \strong{covergence} : Logical; If TRUE, then the optimization algorithm has converged
#' }
#'
#' @author Kin Yau (Alex) Wong <kin-yau.wong@polyu.edu.hk>
#' @references Lam, K. F., Lee, C. Y., Wong, K. Y., & Bandyopadhyay, D. (2021). Marginal analysis of current status data with informative cluster size using a class of semiparametric transformation cure models. Statistics in Medicine [online], DOI:10.1002/sim.8910
#' @seealso \code{BernsteinPolynomial}
#' @examples Dataset <- ICDASim(seed = 1942, n = 100, beta00 = 0.5, beta10 = -1, beta20 = 1,
#' cs = 40, rho=0.5,gamma = 0.5)
#' Result <- Est.ICScure(data = Dataset, rho = 0.5, degree = 3, weighting = FALSE)
#' @export

Est.ICScure <- function(data, rho, degree, weighting=TRUE, t.min=0, t.max=NA, reltol=10^-6, maxit=1000, calSE=TRUE){

  if(is.na(t.max)){
  t.max<-max(data$Ci*data$delta)
  }

  data <- cbind(data[,1],rep(0,nrow(data)),data[,2:ncol(data)])
  names(data) <- c("family","cs",names(data)[3:ncol(data)])
  data <- data[order(data[,1],decreasing = F),]
  for ( i in 1:nrow(data)){
    data$cs[i] <- length(data$family[data$family==data$family[i]])
  }

  if(!weighting){data$cs<-1}

  x0num <- ncol(data) - 4 + 1

  Maximization<-optim(c(rep(0,degree),rep(0,x0num)),
                      f = LOG.Lik.csw, NULL, hessian = FALSE, method = "BFGS",
                      rho=rho, degree=degree,data=as.matrix(data),x0num=x0num, t.min=t.min,t.max=t.max,
                      control=list(maxit=maxit,reltol=reltol,ndeps=rep(10^-8,times=degree+x0num)))

  #AIC<-2*Maximization$value+2*length(c(rep(1,degree),rep(0,x0num)))

  psi <- sapply(c(1:degree),function(a){sum(exp(Maximization$par[1:a]))})/sum(exp(Maximization$par[1:degree]))

  if(calSE){
    Hess <- hessian(func = LOG.Lik.csw, x = Maximization$par, rho=rho, degree=degree,data=as.matrix(data), x0num=x0num, t.min=t.min, t.max=t.max) #Hessian
  }else{
    Hess <- NA
  }

  if((sum(is.na(Hess))>=1)){
    results<-list(degree,psi,
             Maximization$par[(degree+1):(degree+x0num)],
             rep(NA,x0num),
             as.numeric(Maximization$counts[1]),as.numeric(Maximization$convergence),AIC,t.max)
  }
  else{
    matrix.store<-matrix(0,nrow=length(c(rep(0,degree),rep(0,x0num))),
                         ncol=length(c(rep(0,degree),rep(0,x0num))))

    for(k in c(1:length(unique(data$family))) ){
      clusteri<-subset(data,data$family==k)
      score_i<-numeric(length(c(rep(0,degree),rep(0,x0num))))
      for(p in c(1:nrow(clusteri))){
        score_pj<-grad(LOG.Lik.IJ.csw, rho=rho, datline=as.numeric(clusteri[p,]), Maximization$par,
                       degree=degree,t.min=t.min,t.max=t.max, x0num=x0num)
        score_i<-score_i+(as.numeric(clusteri$cs[1]))^(-1)*score_pj
      }
      matrix.store<-matrix.store+as.matrix(score_i)%*%t(as.matrix(score_i))
    }

    sandwich<-ginv(Hess,tol = 10^-30)%*%matrix.store%*%ginv(Hess,tol = 10^-30) #Sandwich estimator

    convergence <- ifelse(as.numeric(Maximization$convergence)==0,TRUE,FALSE)
    results<-list(degree,psi,
                  Maximization$par[(degree+1):(degree+x0num)],
                  sqrt(diag(sandwich))[(degree+1):(degree+x0num)],
                  as.numeric(Maximization$counts[1]),as.character(convergence))
  }

  names(results)<-c("degree","psi","beta","betaSE","iteration","covergence")

return(results)
}


