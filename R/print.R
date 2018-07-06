#' @title  print.CoDiNA
#' @aliases print.CoDiNA
#' @description Print on the screen the number of nodes and edges. To see the data.frame, call: data.frame().
#' @param x Output from MakeDiffNet
#' @param \dots Additional plotting parameters.

#' @author Deisy Morselli Gysi <deisy at bioinf.uni-leipzig.de>
#' @return Print on the screen the number of nodes and edges.
#' @method print CoDiNA
#' @export
#' @examples

#' Nodes = LETTERS[1:10]
#' Net1 = data.frame(Node.1 = sample(Nodes) , Node.2 = sample(Nodes), wTO = runif(10,-1,1))
#' Net2 = data.frame(Node.1 = sample(Nodes) , Node.2 = sample(Nodes), wTO = runif(10,-1,1))
#' Net3 = data.frame(Node.1 = sample(Nodes) , Node.2 = sample(Nodes), wTO = runif(10,-1,1))

#' DiffNet = MakeDiffNet (Data = list(Net1,Net2,Net3), Code = c('Net1', 'Net2', 'Net3') )
#' print(DiffNet)

# print <- function(x){
#   UseMethod('CoDiNA')
# }

print.CoDiNA <- function(x, ...) {
  Nodes_n = length(unique(c(as.character(x$Node.1),
                            as.character(x$Node.2))))
  cat('Nodes', Nodes_n, "\n")
  cat('Links', nrow(x), "\n")
  as.data.frame(x)
}
