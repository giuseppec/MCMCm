\name{MCMCm}
\alias{MCMCm}
\title{STAR LVMs}
\description{
  blabla
}
\usage{
MCMCm(ind.form, 
      fixed.form = NA, 
      covar.form = NA, 			###	3 Formeln fuer Response, Fixed Eff., Covariates
      offset.form=NA,
 			data=parent.environment(), 
      factors=1, 
      lambda.constraints=list(), 
      burnin = 500, 
      mcmc = 1000,
      thin=1, 
      verbose = TRUE, 
      seed = NA,
      beta_stern = 0, 
      Beta_stern = 0, 						## Priori fuer Param.vektor beta
      gamma_stern = 0, 
      Gamma_stern = 0, 					## Priori fuer Param.vektor gamma
      v_stern = 0, 
      s_stern = 0,	## Priori fuer Param.vektor phi* (metrische Var.)
      beta.start = 0, 
      gamma.start = 0,						## Startwerte fuer Param.vektoren beta und gamma
      store.scores = FALSE, 
      store.z = FALSE,
      mh=FALSE, 
      tune=0.10, 
      gm = TRUE,						## Angabe, ob Grouped Move Step 1 durchgefuehrt werden soll
      DIC=FALSE, 
      sim = 0, 
      ... )
}
\arguments{
  \item{ind.form}{...} 
  \item{fixed.form}{...} 
  \item{covar.form}{...} 
	\item{offset.form}{...} 
 	\item{data}{...} 
  \item{factors}{...} 
 	\item{lambda.constraints}{...} 
  \item{burnin}{...}  
  \item{mcmc}{...} 
  \item{thin}{...} 
  \item{verbose}{...} 
  \item{seed}{...} 
  \item{beta_stern}{...} 
  \item{Beta_stern}{...}
  \item{gamma_stern}{...} 
  \item{Gamma_stern}{...} 
  \item{v_stern}{...} 
  \item{s_stern}{...} 
  \item{beta.start}{...} 
  \item{gamma.start}{...} 	
  \item{store.scores}{...} 
  \item{store.z}{...} 
  \item{mh}{...} 
  \item{tune}{...} 
  \item{gm}{...} 
  \item{DIC}{...} 
  \item{sim}{...} 
  \item{\dots}{Additional arguments.}
}
\details{
blabla.
}
\value{
   An mcmc object.
}
\keyword{models}
