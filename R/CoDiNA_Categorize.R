#' @title PhiCategory
#' @description Categorize the links into Phi and Phi tilda categories.
#' @param Base data.frame to be categorized
#' @param n number of networks to be compared
#' @author Deisy Morselli Gysi <deisy at bioinf.uni-leipzig.de>
#' @keywords internal
#'
PhiCategory = function(Base, n){
  message('Starting Phi categorization.')
  Base$Phi_tilda = ifelse(test = (abs(Base$sum) == Base$sum_abs & Base$sum_abs == 1) , yes = paste('g',Base$pos,Base$neg, sep = '.'),
                           no = ifelse(test = (abs(Base$sum) == Base$sum_abs & Base$sum_abs == n), yes = paste('a',Base$pos,Base$neg, sep = '.'),
                                       no = ifelse(test = (abs(Base$sum_abs)) == n, yes = paste('b',Base$pos,Base$neg, sep = '.'),
                                                   no = paste('g',Base$pos,Base$neg, sep = '.'))))
  Base$Phi = ifelse(test = (abs(Base$sum) == Base$sum_abs & Base$sum_abs == 1) , yes = 'g',
                   no = ifelse(test = (abs(Base$sum) == Base$sum_abs & Base$sum_abs == n), yes = 'a',
                               no = ifelse(test = (abs(Base$sum_abs)) == n, yes = 'b',
                                           no = 'g')))

  return(Base)
}
