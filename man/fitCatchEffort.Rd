\name{fitCatchEffort}
\alias{fitCatchEffort}
\title{
  Fit a Catch-Effort Model on Wild Boar Data
}
\description{
  This function fits a catch-effort model as described in Vajas et
  al. (in prep).
}
\usage{
fitCatchEffort(listData, n.chains = 4, thin = 10, n.iter = 1e+05, n.burnin = 5000)
}
\arguments{
  \item{listData}{
    a list formatted like \code{bagArdeche} or \code{bagHerault}
  }
  \item{n.chains}{
    Integer. the number of parallel chains for the MCMC (see the same argument in
    the function \code{jags.model} in the package \code{rjags}).
  }
  \item{thin}{
    Integer. Thinning interval for monitors (see the same
    argument in the function \code{coda.samples} in the package
    \code{coda})
  }
  \item{n.iter}{
    Integer. Number of iterations monitor (see the same argument in the
    function \code{coda.samples} in the package \code{coda}).
  }
  \item{n.burnin}{
    Integer. Number of burnin iterations (thrown away before monitoring iterations).
  }
}
\details{
This function fits the catch-effort model to the list of hunting logs
data, through "rjags" and "coda" packages. 
}
\value{
  This function returns an object of class \code{mcmc.list} (this is the
  object returned by the function \code{coda.sample} used by
  \code{fitCatchEffort} to fit the model. See this function for more details on the
  format of such objects).
}
\references{
  Vajas, P., Calenge, C., Gamelon, M., Girard, F., Melac, O., Chandosne,
  C., Richard, E., Said, S., Baubet, Eric. in prep : Catch-effort model
  used as a management tool in exploited populations: wild boar as a
  case study. 
}
\author{
  Pablo Vajas \email{pablo.vajas@ofb.gouv.fr}
}
\seealso{
  \code{\link[coda]{coda.sample}}, \code{\link[rjags]{jags.model}},
  \code{\link{bagArdeche}}, \code{\link{bagHerault}}
}
\examples{
\dontrun{

## Loads the dataset containing the hunting bags in Ardeche
data(bagArdeche)

## fit the model (WARNING : can take several hours !)
mod <- fitCatchEffort(bagArdeche, n.chains = 4, thin = 1, n.iter = 1000,
                      n.burnin = 50)

## In the article, the values used were
## n.chains = 4
## thin = 10
## n.iter = 100000
## n.burnin = 1000


}
}
\keyword{ model }
