
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
#' @importFrom stats median
#' @export
summary.pm_project <- function(object, ...) {
  
  # print data details
  message("# DATA:")
  if (is.null(object$data$coords)) {
    message("no data loaded")
  } else {
    message(sprintf("samples = %s", nrow(object$data$coords)))
  }
  message("")
  
  # print map details
  message("# HEX MAP:")
  if (is.null(object$map)) {
    message("no map generated")
  } else {
    message(sprintf("eccentricity = %s", object$map$eccentricity))
    message(sprintf("hexagons = %s", length(object$map$hex)))
  }
  message("")
  
  # print output details
  message("# OUTPUT:")
  if (length(object$output) == 0) {
    message("no output")
  } else {
    hex_coverage <- object$output$hex_coverage
    message(sprintf("coverage (median) = %s", median(hex_coverage, na.rm = TRUE) ))
    message(sprintf("coverage (range) = %s - %s", min(hex_coverage, na.rm = TRUE), max(hex_coverage, na.rm = TRUE) ))
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

