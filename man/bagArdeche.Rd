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
    \item{Nua}{An unique integer vector value corresponding to the length of "nmus": here 308}
    
    \item{cull}{Integer value vector containing the hunting bag resulting
      from each day of hunting. int [1:28573]}

    \item{nhunters}{Integer value vector containing the number of hunters from each day of hunting. int [1:28573]}

    \item{mu}{Integer value vector of le MU number from each day of hunting. int [1:28573]}

    \item{idstart}{Integer value vector of the line index in this table corresponding to the first hunting day in a MU / year. in [1:308]}

    \item{month}{Numerical values vector of the month number from each day of hunting. num [1:28573]}

    \item{nmus}{Integer value vector of the number of hunting days for each MU/year. int [1:308]}

    \item{mug}{Numerical values vector corresponding to the MU-year in each component of "nmus". [1:308]}

    \item{type}{Numeric value vector of the type of territory associated with each MU. num [1:28]}

    \item{area}{Numerical value vector of the MU size (hectare). num [1:28]}

    \item{nmonths}{Unique numerical value vector of the month number. 6}

    \item{ntypes}{Unique numerical value vector of the type number. 3}
  }
}
\details{
  See the help page of \code{fitCatchEffort} for an example of
  catch-effort model fit on this dataset
}
\source{
  Fédération Départementale des Chasseurs de l'Ardèche. Adresse: Fédération Départementale des Chasseurs de l’Ardèche, l’Escrinet, 07200 Saint-Etienne-de-Boulogne
}
\references{
  Vajas, P., Calenge, C., Gamelon, M., Girard, F., Melac, O., Chandosne, C., Richard, E., Said, S., Baubet, Eric. in prep : Catch-effort model used as a management tool in exploited populations: wild boar as a case study.
}
\examples{
data(bagArdeche)
str((bagArdeche)
}
\keyword{datasets}
