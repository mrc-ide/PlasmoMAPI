
#### Member functions for objects of class pm_sim

#------------------------------------------------
# Overload print()
#' @method print pm_sim
#' @export
print.pm_sim <- function(x, ...) {
  
  # print summary
  summary(x)
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
# Overload summary()
#' @method summary pm_sim
#' @export
summary.pm_sim <- function(object, ...) {
  
  # print data details
  message("Simulation details:")
  message(sprintf("  demes = %s", length(object$daily_values)))
  message(sprintf("  max time = %s", nrow(object$daily_values[[1]])))
  
}

#------------------------------------------------
# Overload plot()
#' @method plot rmapi_project
#' @export
plot.pm_sim <- function(x, y, ...) {
  plot_daily_states(x)
}

#------------------------------------------------
#' @title Determine if object is of class pm_sim
#'
#' @description Determine if object is of class \code{pm_sim}.
#'
#' @param x object to query.
#'
#' @export

is.pm_sim <- function(x) {
  inherits(x, "pm_sim")
}

