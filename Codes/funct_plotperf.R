get.error.rate <- function (x,
                            dist = c("all","max.dist","centroids.dist","mahalanobis.dist"),
                            measure = c("all","overall","BER"),
                            weighted = TRUE,
                            xlab = NULL,
                            ylab = NULL,
                            legend.position=c("vertical","horizontal"),
                            sd = TRUE,
                            ...) {
  measure.input = measure
  measure = NULL
  if(any(measure.input == "all"))
    measure.input = c("BER", "overall")
  
  if(any(measure.input == "BER"))
    measure = c(measure, "Overall.BER")
  
  if (any(measure.input == "overall"))
    measure = c(measure, "Overall.ER")
  
  if(!all(measure.input %in% c("all", "overall", "BER")))
    stop("'measure' must be 'all', 'overall' or 'BER'")
  
  if (any(dist == "all"))
    dist = colnames(x$error.rate[[1]])
  
  if (is.null(dist) || !any(dist %in% colnames(x$error.rate[[1]])))
    stop("'dist' should be among the ones used in your call to 'perf': ", paste(colnames(x$error.rate[[1]]),collapse = ", "),".")
  
  
  
  if(weighted == TRUE)
  {
    perfo = "WeightedVote.error.rate"
    perfo.sd = "WeightedVote.error.rate.sd"
  } else {
    perfo = "MajorityVote.error.rate"
    perfo.sd = "MajorityVote.error.rate.sd"
  }
  
  if(sd == TRUE)
  {
    if(is.null(x[[perfo.sd]]))
      sd = FALSE
  }
  
  error.rate = error.rate.sd = list()
  for(mea in measure)
  {
    error.temp = error.temp.sd = NULL
    for(di in dist)
    {
      temp = t(x[[perfo]][[di]][mea, , drop=FALSE])
      colnames(temp) = di
      error.temp = cbind(error.temp, temp)
      if(sd)
      {
        temp.sd = t(x[[perfo.sd]][[di]][mea, , drop=FALSE])
        colnames(temp.sd) = di
        error.temp.sd = cbind(error.temp.sd, temp.sd)
      }
      
    }
    error.rate[[mea]] = error.temp
    if(sd)
    {
      error.rate.sd[[mea]] = error.temp.sd
    } else {
      error.rate.sd = NULL
    }
  }
  
  if (is.null(ylab))
    ylab = 'Classification error rate'
  
  if (is.null(xlab))
    xlab = 'Component'
  
  error.rate <- list(error.rate, error.rate.sd)
  names(error.rate) <- c("ER", "ER_sd")
  
  return(error.rate)
}

internal_graphic.perf <- function(error.rate.list, ylab = NULL, xlab = NULL) {
  require(tidyr)
  
  if (is.null(ylab))
    ylab = 'Classification error rate'
  
  if (is.null(xlab))
    xlab = 'Component'
  
  # Error rate
  dt.error <- as.data.frame(do.call("rbind", error.rate.list[[1]]))
  dt.error$type <- rep(c("BER", "ER"), each = 5)
  dt.error$comp <- rep(1:(nrow(dt.error)/2), 2)
  rownames(dt.error) <- NULL
  dt.error.long <- gather(dt.error, distance, error.rate, max.dist:mahalanobis.dist, factor_key = TRUE)
  
  # SD error rate
  dt.error.sd <- as.data.frame(do.call("rbind", error.rate.list[[2]]))
  dt.error.sd$type <- rep(c("BER", "ER"), each = 5)
  dt.error.sd$comp <- rep(1:(nrow(dt.error.sd)/2), 2)
  rownames(dt.error.sd) <- NULL
  dt.error.sd.long <- gather(dt.error.sd, distance, sd, max.dist:mahalanobis.dist, factor_key = TRUE)
  
  dt <- merge(dt.error.long, dt.error.sd.long, by = c("row.names", "type", "distance", "comp"))
  dt <- subset(dt, select = -Row.names)
  
  ylim = range(c(dt$error.rate + dt$sd), c(dt$error.rate - dt$sd))
  
  p <- ggplot(dt, aes(x=comp, y=error.rate, color=distance, linetype = type)) + 
    geom_line() +
    geom_point() +
    labs(x = xlab,
         y = ylab,
         color = "Distance",
         linetype = "Error Type") +
    scale_color_manual(values = c("grey47", "steelblue4", "orangered2")) +
    geom_errorbar(aes(ymin = error.rate-sd, ymax = error.rate + sd), width=.2,
                  position=position_dodge(0.05)) +
    theme_bw()
  return(p)
}