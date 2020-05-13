\name{bagArdeche}
\alias{bagArdeche}
\docType{data}
\title{
  Hunting Bag for the Ardeche Department
}
\description{
    This list contains the data required to fit the model in Vajas et
    al. (in prep.) for the Ardeche department.  For anonymity reasons,
    the Management Units number have been randomized.
}
\usage{data("bagArdeche")}
\format{
  A list with the following elements:

  \describe{
    \item{Nua}{Integer value. Number of xxxx}
    
    \item{cull}{a numeric vector containing the hunting bag resulting
      from each day of hunting. XXX}

    \item{nhunters}{Number of hunters}

    \item{mu}{Number of XXXX}

    \item{idstart}{Number of XXXX}

    \item{month}{Number of XXXX}

    \item{mug}{Number of management units}

    \item{type}{Number of XXXX}

    \item{nmus}{Number of XXXX}

    \item{area}{Number of XXXX}

    \item{nmonths}{Number of XXXX}

    \item{ntypes}{Number of XXXX}
  }
}
\details{
  See the help page of \code{fitCatchEffort} for an example of
  catch-effort model fit on this dataset
}
\source{
  Fédération Départementale des Chasseurs de l'Ardèche. Adresse.XXX
}
\references{
  Vajas, P., Calenge, C. etc. Titre. in prep.XXXX
}
\examples{
data(bagArdeche)
str((bagArdeche)
}
\keyword{datasets}
