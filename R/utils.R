#' Function to peform min max scale
#'
#' This function perform min max scale
#' @param x input vector
#' @return a vector in the interval of 0 between 1
#' @export
#' @examples
#' MinMaxScale(rnorm(100))
MinMaxScale <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#' Function to remove NA
#'
#' This function replace NA by 0
#' @param df data.frame
#' @return a data.frame without NA
#' @export
#' @examples
#' myRemoveNA(data.frame(x = c(0,1,NA),y = c(1,2,3)))

myRemoveNA <- function(df){
  df[is.na(df)] <- 0
  return(df)
}


