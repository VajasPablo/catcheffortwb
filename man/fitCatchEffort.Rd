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
    XX%%     ~~Describe \code{n.chains} here~~
  }
  \item{thin}{
    XX%%     ~~Describe \code{thin} here~~
}
  \item{n.iter}{
    XX%%     ~~Describe \code{n.iter} here~~
}
  \item{n.burnin}{
    XX%%     ~~Describe \code{n.burnin} here~~
  }
}
\details{
This function fits XX%%  ~~ If necessary, more details than the description above ~~
}
\value{
  This function returns an object of class xxx returned by the function
  \code{coda.sample} used internally to fit the model in the
  function. See the help page of this function for more details on the
  format of such objects.
}
\references{
  Vajas, P. etc.
}
\author{
  Pablo Vajas \email{pablo.vajas@ofb.gouv.fr}
}
\seealso{
  \code{\link[coda]{coda.sample}}, \code{\link{bagArdeche}}, \code{\link{bagHerault}}
}
\examples{
\dontrun{

## Loads the dataset containing the hunting bags in Ardeche
data(bagArdeche)

## fit the model (WARNING : can take several hours)
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
