\name{envsig}
\alias{envsig}
\title{ Determine the Significance of Pointwise Variogram Envelope}
\description{
Determine the significance of pointwise variogram envelope based on permutation test.
}
\usage{
envsig(envlist, adjustment = c("holm", "hochberg", "hommel",
       "bonferroni", "BH", "BY", "fdr", "none"))
}

\arguments{
\item{envlist}{ Output from \code{envelope} function.}
\item{adjustment}{ One of p value adjustment methods.}
}

\value{ A list contains:
\item{p.pointwise}{ Adjusted pointwise p-values.}
\item{p.overall}{ Overall variogram p-value using Fisher's method.}
}

\references{Diblasi, A. and Bowman, A. W. (2001), On the Use of the Variogram in Checking for Independence in Spatial Data. Biometrics, 57: 211-218. doi:10.1111/j.0006-341X.2001.00211.x

Diggle, P.J., Ribeiro, P.J. (2007), Model-Based Geostatistics, Springer Series in Statistics. Springer. doi:10.1007/978-0-387-48536-2}

\author{ Craig Wang }
\seealso{
\code{\link{envelope}} to generate Monte Carlo permutations for computing the envelope.
}
\examples{
\dontrun{
library(sp)
data(meuse)
coordinates(meuse) = ~x+y
vario0 <- gstat::variogram(log(zinc)~1, meuse)
varioEnv <- envelope(vario0, data = meuse, formula = log(zinc)~1, nsim = 999)
envplot(varioEnv)
envsig(varioEnv, adjustment = "holm")
}
}
