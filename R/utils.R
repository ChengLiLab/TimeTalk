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
