
#' Draws a traceplot
#'

#'
#' @param df data frame
#' @param var variable under consideration
#' @param thin
#'
#' @return a traceplot
#' @export
#'
draw_tp <- function(df, var, thin){
  tp <- df  %>%
    filter(row_number() %% thin == 1) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$Iteration, y = .data[[var]])) +
    ggplot2::geom_line() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 15),
      axis.text.x = ggplot2::element_text(size = 15)
    )
  return(tp)
}


#' Generates a histogram
#'
#' @inheritParams draw_tp
#'
#' @param bins number of bins in the histogram
#'
#' @return a histogram
#' @export
#'
draw_histogram <- function(df, var, bins = 10) {

  g <- df %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[[var]])) +
    ggplot2::geom_histogram(ggplot2::aes(y = .data$..count../sum(.data$..count..)), bins = bins) +
    ggplot2::theme(text = ggplot2::element_text(size = 15)) +
    ggplot2::labs(y = "Proportion")

  return(g)
}


#' Generates a histogram with marking lines for quantile locations
#'
#' @inheritParams draw_tp
#'
#' @param bins number of bins in the histogram
#' @param lower_quantile
#' @param upper_quantile
#'
#' @return a histogram
#' @export
#'
draw_histogram2 <- function(df, var, bins = 10, lower_quantile, upper_quantile) {

  g <- df %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[[var]])) +
    ggplot2::geom_histogram(ggplot2::aes(y = .data$..count../sum(.data$..count..)), bins = bins) +
    ggplot2::theme(text = ggplot2::element_text(size = 15)) +
    ggplot2::labs(y = "Proportion")+
    geom_vline(xintercept = lower_quantile, color = "red") +
    geom_vline(xintercept = upper_quantile, color = "red")

  return(g)
}


#' Generate a plot of the auto-correlation function
#'
#' @inheritParams draw_traceplot
#'
#' @param lag_max largest lag shown
#'
#' @return Barplot of the ACF
#' @export
#'
draw_acf <- function(df, var, lag_max = 200) {

  x     <- dplyr::pull(df, .data[[var]])
  x_acf <- stats::acf(x, lag.max = lag_max, plot = FALSE)
  df    <- tibble::tibble(ACF = as.numeric(x_acf[["acf"]]), Lag = as.numeric(x_acf[["lag"]]))

  g <- df %>%
    ggplot2::ggplot(ggplot2::aes(.data$Lag, .data$ACF)) +
    ggplot2::geom_col(width = .25) +
    ggplot2::theme(text = ggplot2::element_text(size = 15)) +
    ggplot2::labs(title =  var)

  return(g)
}
