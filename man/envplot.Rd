\name{envplot}
\alias{envplot}
\title{ Plot Empirical Variogram and its Pointwise Envelope}
\description{
Plot the empirical variogram along with pointwise variogram envelope based on Monte Carlo permutations.
}
\usage{
envplot(envlist, shade = TRUE, shade.color = "lightgrey",
        xlim = NULL, ylim = NULL, main = NULL,
        xlab = "Distance", ylab = "Semivariance")
}

\arguments{
\item{envlist}{ Output from \code{envelope} function.}
\item{shade}{ logical. If the envelope should be shaded or not.}
\item{shade.color}{ Color of the envelope if it is shaded.}
\item{xlim, ylim}{ x-axis and y-axis range. If \code{NULL}, default values are chosen.}
\item{main}{ Title.}
\item{xlab, ylab}{ x-axis and y-axis labels.}
}

\value{ Nothing is returned, generates a plot.}

\author{ Craig Wang }
\seealso{
\code{\link{envelope}} to generate Monte Carlo permutations for computing the envelope.
}
\examples{
library(sp)
data(meuse)
coordinates(meuse) = ~x+y
vario0 <- gstat::variogram(log(zinc)~1, meuse)
varioEnv <- envelope(vario0, data = meuse, formula = log(zinc)~1, nsim = 99)
envplot(varioEnv)
}
