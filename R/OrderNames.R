
#' @title OrderNames
#' @description Sorts each link's Nodes by the smallest value. Removes links that both nodes are the same.
#' @param M data.frame to have the names ordered. Node.1, Node.2 and correlation value.
#'
#' @return a data.table whith Node.1 and Node.2, sorted by the smallest value between both.
#' @author Deisy Morselli Gysi <deisy at bioinf.uni-leipzig.de>
#' @export
#' @importFrom data.table data.table
#' @examples
#' Nodes = LETTERS[1:10]
#' Z = data.frame(Node.1 = sample(Nodes) ,
#' Node.2 = sample(Nodes), cor = runif(10,-1,1))
#' OrderNames(Z)

OrderNames=function (M){
  # message('Ordering names.')
  M[,1]= as.character(M[,1])
  M[,2]= as.character(M[,2])
  n1 = apply(M[,1:2], 1, min)
  n2 = apply(M[,1:2], 1, max)
  df = data.table::data.table(Node.1 = n1,
                              Node.2 = n2,
                              cor = M[,3])
  names(df)[3] = 'cor'
  df = subset(df, df$Node.1 != df$Node.2)
  df = subset(df, !is.na(df$Node.1))
  df = subset(df, !is.na(df$Node.2))
  df = subset(df, !is.na(df$cor))
  return(df)
}
