
#### Member functions for objects of class pm_project

#------------------------------------------------
# Overload print()
#' @method print pm_project
#' @export
print.pm_project <- function(x, ...) {
  
  # print summary
  summary(x)
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
# Overload summary()
#' @method summary pm_project
#' @export
summary.pm_project <- function(object, ...) {
  
  # print data details
  message("# DATA:")
  if (is.null(object$data$coords)) {
    message("no data loaded")
  } else {
    message(sprintf("n = %s samples", nrow(object$data$coords)))
  }
  message("")
  
  # print map details
  message("# HEX MAP:")
  if (is.null(object$map)) {
    message("no map generated")
  } else {
    message(sprintf("h = %s hexagons", length(object$map$hex)))
  }
  message("")
  
  # print output details
  message("# OUTPUT:")
  if (length(object$output) == 0) {
    message("no output")
  } else {
    message("TODO - some details of output")
  }
}

#------------------------------------------------
# Overload plot()
#' @method plot pm_project
#' @export
plot.pm_project <- function(x, y, ...) {
  plot_map(x)
}

#------------------------------------------------
#' @title Determine if object is of class pm_project
#'
#' @description Determine if object is of class \code{pm_project}.
#'
#' @param x object to query.
#'
#' @export

is.pm_project <- function(x) {
  inherits(x, "pm_project")
}

