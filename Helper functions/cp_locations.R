# Determine change point locations.
#' Function to return segment_start and segment_end locations given delta.
#'
#' @param delta delta beta vector (current)
#' @param Y observed data
#'
#' @return segment_start and segment_end.
#' @export
#'
cp_locations <- function(delta, Y) {

  segment_start <- c(0)
  segment_start <- append(segment_start, which(delta == 1))
  segment_end <- c(Y[["t_end"]])
  segment_end <- append(which(delta == 1), segment_end)

  return(list("segment_start" = segment_start, "segment_end"=segment_end))
}
