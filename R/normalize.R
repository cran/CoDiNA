#' @title normalize
#' @description Normalize a given variable.
#' @param m variable to be normalized in the interval [0,1]
#'
#' @export
#' @author Deisy Morselli Gysi <deisy at bioinf.uni-leipzig.de>
#' @examples
#' Z = runif(10,-10,10)
#' normalize(Z)
normalize<-function(m){
  if(length(m)>1){
    m =  (m - min(m))/(max(m)-min(m))
  }
  else{ m = 1}
}
