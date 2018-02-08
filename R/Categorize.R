#' @title Categorize
#' @description Categorize the links into -1, 0 and 1 given a cutoff.
#' @param M A data.frame to be categorized.
#' @param cutoff By default, the cutoff is 0.33. If the user wants to use another value, it has to be cited on the description of the used methodology that the cutoff was changed.
#' @author Deisy Morselli Gysi <deisy at bioinf.uni-leipzig.de>
#' @keywords internal

Categorize = function(M, cutoff = 0.33){
  message('Coding correlations.')
  M[M > cutoff]<- 1
  M[M < -cutoff]<- -1
  M[M >= -cutoff & M <= cutoff]<-0
  return(M)
}

