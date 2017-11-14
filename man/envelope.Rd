\name{envelope}
\alias{envelope}
\alias{envelope.gstatVariogram}
\alias{envelope.variogram}
\title{Compute Pointwise Variogram Envelope}
\description{
Apply Monte Carlo permutation on the data values over different locations, and compute the variogram estimates for each permutation. The pointwise variogram envelope is taken based on quantiles of variogram estimates at each distance.
}
\usage{
\method{envelope}{gstatVariogram}(vario, data, locations = coordinates(data), formula = NULL,
  nsim = 999, conf.level = 0.95, save.sim = FALSE, ...)

\method{envelope}{variogram}(vario, data, locations = data$coords, trend = NULL,
  nsim = 999, conf.level = 0.95, save.sim = FALSE, ...)
}

\arguments{
\item{vario}{Variogram of class \code{gstatVariogram} from the \code{gstat} package or \code{variogram} from the \code{geoR} package}
\item{data}{data frame where the names in formula are to be found}
\item{locations}{ spatial data locations.  For variogram.formula: a
formula with only the coordinate variables in the right hand (explanatory
variable) side e.g. \code{~x+y}; see examples.}
\item{formula}{only work with \code{gstat} package, a formula defining the response vector and possible covariates.}
\item{trend}{only work with \code{geoR} package, definng the response vector and possible covariates. See documentation of \code{\link[geoR]{trend.spatial}} for further details.}
\item{nsim}{integer. The number of permutations to be conducted to compute the envelope. Defaults to 999.}
\item{conf.level}{confidence level of the envelope. Defaults to 0.95}
\item{save.sim}{logical. Whether to save all the permutated data.}
\item{...}{any other arguments that will be passed to variogram computation using either \code{\link[gstat]{variogram}} from \code{gstat} package or \code{\link[geoR]{variog}} from \code{geoR} package, based on the class of the object argument.}
}

\value{
A list with the following elements:
\item{data}{contains all the permutated data, only shown if \code{save.sim=TRUE}}
\item{variogram}{the variogram estimates from each of the permutations}
\item{upper}{upper quantile of the variogram estimates according to the confidence level}
\item{lower}{lower quantile of the variogram estimates according to the confidence level}
\item{variogram0}{the original variogram supplied via \code{vario} argument}
}

\references{
Diggle, P.J., Ribeiro, P.J. (2007), Model-Based Geostatistics, Springer Series in Statistics. Springer. doi:10.1007/978-0-387-48536-2}

\author{ Craig Wang }
\seealso{
\code{\link{envplot}} to plot an empirical variogram based on the data and its envelope
}
\examples{

library(sp)
data(meuse)

## for gstat variogram

library(gstat)
data(meuse)
coordinates(meuse) = ~x+y
vario0 <- gstat::variogram(log(zinc)~1, meuse)
varioEnv <- envelope(vario0, data = meuse, formula = log(zinc)~1, nsim = 99)

## for geoR variog

library(geoR)
data(meuse)
meuse <- as.geodata(meuse, coords.col = 1:2, data.col = 6, covar.col = 3:5)
meuse$data <- log(meuse$data)
vario0 <- variog(meuse,max.dist=1500)
varioEnv <- envelope(vario0, data = meuse, trend = trend.spatial(~lead, meuse), nsim = 99)
}

\keyword{models}
