envelope <- function(vario, ...){
  if (!class(vario)[1] %in% c("gstatVariogram","variogram")) stop(
"The method 'envelope' must be applied to an object of class either 'gstatVariogram' or'variogram'")
  UseMethod("envelope",vario)
}

envelope.gstatVariogram <- function(vario, data, locations = coordinates(data), formula = NULL,
                                    nsim = 999, conf.level = 0.95, save.sim = FALSE, ...) {
  dots <- list(...)
  variogramDefault <- vario
  if (!is.null(formula)){
    if (inherits(formula, "formula")){
      dataValues <- lm(formula, data = data)$residuals
      dots$object <- NULL
    }}
  simulation <- list()
  simulation$data <- sapply(1:nsim, function(i) dataValues[sample(1:length(dataValues))])

  simulation$variogram <- sapply(1:nsim, function(i){
      data.temp <- data.frame(locations, simulation$data[,i])
      names(data.temp) <- c("x","y","residuals")
      sp::coordinates(data.temp) = ~ x + y
      gstat::variogram(residuals~1, locations, data.temp, ...)$gamma
    })

  simulation$upper <- apply(simulation$variogram, 1, quantile, probs = 1-(1-conf.level)/2)
  simulation$lower <- apply(simulation$variogram, 1, quantile, probs = (1-conf.level)/2)
  simulation$variogram0 <- variogramDefault
  simulation$conf.level <- conf.level
  if (!save.sim){simulation$data <- NULL}
  return(simulation)
}

envelope.variogram <- function(vario, data, locations = data$coords, trend = NULL,
                               nsim = 999, conf.level = 0.95, save.sim = FALSE, ...) {
  dots <- list(...)
  dataValue <- data$data
  variogramDefault <- vario
  if (!is.null(trend)){
    dataValues <- lm(dataValue ~ trend-1)$residuals
  } else {
    dataValues <- dataValue
  }
  simulation <- list()
  simulation$data <- sapply(1:nsim, function(i) dataValues[sample(1:length(dataValues))])

  simulation$variogram <- sapply(1:nsim, function(i){
      data.temp <- data.frame(locations, simulation$data[,i])
      data.temp <- as.geodata(data.temp, coords.col = 1:2, data.col = 3)
      geoR::variog(data.temp, ...)$v
  })

  simulation$upper <- apply(simulation$variogram, 1, quantile, probs = 1-(1-conf.level)/2)
  simulation$lower <- apply(simulation$variogram, 1, quantile, probs = (1-conf.level)/2)
  simulation$variogram0 <- variogramDefault
  simulation$conf.level <- conf.level
  if (!save.sim){simulation$data <- NULL}
  return(simulation)
}

