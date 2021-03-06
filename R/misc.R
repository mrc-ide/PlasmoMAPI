
#------------------------------------------------
# if NULL then replace with chosen value, otherwise keep original value
#' @noRd
define_default <- function(x, default) {
  if (is.null(x)) {
    return(default)
  } else {
    return(x)
  }
}

#------------------------------------------------
# if a single value is provided then expand to a vector of length n
#' @noRd
force_vector <- function(x, n) {
  if (length(x) == 1) {
    return(rep(x,n))
  } else {
    return(x)
  }
}

#------------------------------------------------
# calculate midpoints of a vector
#' @noRd
midpoints <- function(x) {
  return((x[-1] + x[-length(x)])/2)
}

# -----------------------------------
# takes matrix as input, converts to list format for use within Rcpp code
#' @noRd
matrix_to_rcpp <- function(x) {
  return(split(x, f = 1:nrow(x)))
}

# -----------------------------------
# takes list format returned from Rcpp and converts to matrix
#' @noRd
rcpp_to_matrix <- function(x) {
  ret <- matrix(unlist(x), nrow = length(x), byrow = TRUE)
  return(ret)
}

# -----------------------------------
# takes list format returned from Rcpp and converts to three-dimensional array.
# Array indexing is in the same order as the underlying list, for example
# x[i,j,k] is equivalent to l[[i]][[j]][[k]]
#' @noRd
rcpp_to_array <- function(x) {
  ret <- array(unlist(x), dim = c(length(x[[1]][[1]]), length(x[[1]]), length(x)))
  ret <- aperm(ret, perm = c(3,2,1))
  return(ret)
}

#------------------------------------------------
# return 95% quantile
#' @importFrom stats quantile
#' @noRd
quantile_95 <- function(x) {
  ret <- quantile(x, probs = c(0.025, 0.5, 0.975))
  names(ret) <- c("Q2.5", "Q50", "Q97.5")
  return(ret)
}

#------------------------------------------------
# sum logged values without underflow, i.e. do log(sum(exp(x)))
#' @noRd
log_sum <- function(x) {
  if (all(is.na(x))) {
    return(rep(NA, length(x)))
  }
  x_max <- max(x, na.rm = TRUE)
  ret <- x_max + log(sum(exp(x - x_max)))
  return(ret)
}

#------------------------------------------------
# update progress bar
# pb_list = list of progress bar objects
# name = name of this progress bar
# i = new value of bar
# max_i = max value of bar (close when reach this value)
# close = whether to close when reach end
#' @importFrom utils setTxtProgressBar
#' @noRd
update_progress <- function(pb_list, name, i, max_i, close = TRUE) {
  setTxtProgressBar(pb_list[[name]], i)
  if (i == max_i & close) {
    close(pb_list[[name]])
  }
}

# -----------------------------------
# ask user a yes/no question. Return TRUE/FALSE.
#' @noRd
user_yes_no <- function(x = "continue? (Y/N): ") {
  user_choice <- NA
  while (!user_choice %in% c("Y", "y" ,"N", "n")) {
    user_choice <- readline(x)
  }
  return(user_choice %in% c("Y", "y"))
}
