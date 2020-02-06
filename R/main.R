#------------------------------------------------
#' @title Square a vector of values
#'
#' @description Simple test function that demonstrates some of the features of
#'   this package by squaring an input vector of values.
#'
#' @param x vector of values.
#'
#' @export
#' @examples
#' # Find square of first 100 values
#' square(1:100)

square <- function(x = 1:5) {

  # print message to console
  message("running R square function")

  # get arguments in list form
  args <- list(x = x)

  # run C++ function with these arguments
  output_raw <- square_cpp(args)

  # some optional processing of output
  message("processing output")
  ret <- output_raw$x_squared

  # return
  return(ret)
}
