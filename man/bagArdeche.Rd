 \name{bagArdeche}
\alias{bagArdeche}
\docType{data}
\title{
  Hunting Bag for the Ardeche Department
}
\description{
    This list contains the data required to fit the model in Vajas et
    al. (in prep.) for the Ardeche department.  For anonymity reasons,
    the Management Units numbers have been randomized.
}
\usage{data("bagArdeche")}
\format{
  A list with the following elements:

  \describe{ \item{Nua}{Integer vector value corresponding to the length
      of the vector \code{nmus}: here, equal to 308}
    
    \item{cull}{Integer vector containing the hunting bag resulting from
      each day of hunting. For a given combination of management unit
      and year, all the hunting days correspond to a batch of successive
      values in this vector. \code{integer[1:28573]}}
    
    \item{nhunters}{Integer vector containing the number of hunters from
      each day of hunting. \code{integer[1:28573]}} 
    
    \item{mu}{Integer vector giving the Management Unit ID for each
      hunting day in \code{cull}.  Note that these IDs have been
      randomized. \code{integer[1:28573]}}

    \item{idstart}{Integer vector giving the index of the first hunting
      day in the batch for a given combination of management unit / year
      in \code{cull}. \code{integer[1:308]}}

    \item{month}{Integer vector giving the month number for each day of
      hunting in \code{cull}. \code{numeric[1:28573]}}

    \item{nmus}{Integer vector giving the number of hunting days in the
      batch of values in \code{cull} for each combination of management
      units and year. \code{integer[1:308]}}

    \item{mug}{Integer vector giving the ID of the combination of
      management units and year in each element of
      \code{nmus}. \code{numeric[1:308]}}

    \item{type}{Integer value vector of the type of territory associated
      with each management unit. \code{numeric [1:28]}} 

    \item{area}{Numeric vector giving the area of the management units
      (hectares). \code{numeric[1:28]}}

    \item{nmonths}{Integer value. Number of months during a hunting
    season: 6.}

    \item{ntypes}{Integer value. Number of types of management units:
      3.}
  }
}
\details{
  See the help page of \code{fitCatchEffort} for an example of
  catch-effort model fit on this dataset
}
\source{
  Fédération Départementale des Chasseurs de l'Ardèche. L’Escrinet,
    07200 Saint-Etienne-de-Boulogne. France.
}
\references{
  Vajas, P., Calenge, C., Gamelon, M., Girard, F., Melac, O., Chandosne,
    C., Richard, E., Said, S., Baubet, Eric. in prep : Catch-effort
    model used as a management tool in exploited populations: wild boar
    as a case study.
}
\examples{
data(bagArdeche)
str((bagArdeche)
}
\keyword{datasets}
