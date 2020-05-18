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
    In \code{jags.model} function in "rags" package : " the number of parallel chains for the model "
  }
  \item{thin}{
    In \code{coda.samples} function in "coda" package : " thinning interval for monitors "
}
  \item{n.iter}{
    In \code{coda.samples} function in "coda" package : " number of iterations to monitor "
}
  \item{n.burnin}{
    In \code{coda.samples} function in "coda" package : number of iteration throwing away at the beginning of an MCMC run
  }
}
\details{
This function fits the list of hunting logs data with the catch-effort model through "rjags" package and "coda" package arguments.
}
\value{
  This function returns an list of Markov Chain Monte Carlo (MCMC) output returned by the function
  \code{coda.sample} used internally to fit the model in the function. See the help page of this function for more details on the
  format of such objects.
}
\references{
  Vajas, P., Calenge, C., Gamelon, M., Girard, F., Melac, O., Chandosne, C., Richard, E., Said, S., Baubet, Eric. in prep : Catch-effort model used as a management tool in exploited populations: wild boar as a case study.
}
\author{
  Pablo Vajas \email{pablo.vajas@ofb.gouv.fr}
}
\seealso{
  \code{\link[coda]{coda.sample}}, \link[rjags]{jags.model}}, \code{\link{bagArdeche}}, \code{\link{bagHerault}}
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
