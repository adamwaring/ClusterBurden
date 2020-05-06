#' Find contiguous ranges of values in a vector of integers
#'
#' @author Adam Waring - adam.waring@@msdtc.ox.ac.uk
#'
#' @param x a vector of integers
#' @param output whether to output a data.table with columns start and end ("dt") or string ("str")
#'
#' @keywords RVAT, distribution, cluster, genetics, case-control, gene
#'
#' @details if output = "dt" then a data.table with the start and end positions of contiguous ranges of integers in the vector.
#' if output = "str" then each range is presented in a single string separated by commas e.g. "1-3,6-10,17-19"
#'
#' @export
#' @examples
#' x = sample(1:10, 5)
#' findRanges(x, "dt")
#' findRanges(x, "str")

findRanges = function(x, output="dt"){

  if(length(x) == 0){

    starts = numeric(0)
    stops = numeric(0)

  } else if(length(x) == 1){

    starts = x
    stops = x

  } else {

    x = sort(unique(x))

    starts = x[1]
    stops = NULL
    for(i in 2:length(x)){
      if(x[i] != x[i-1] + 1){
        starts = c(starts, x[i])
        stops = c(stops, x[i-1])
      }
    }
    stops = c(stops, x[length(x)])

  }

  x = data.table(start=starts, end=stops)

  if(output=="dt"){
    return(x)
  } else if(output=="str"){
    x = paste(apply(x, 1, paste, collapse="-"), collapse=",")
    return(x)
  }

}
