#' @title  as.igraph
#' @aliases as.igraph
#' @description Converts the CoDiNA.plot into an igraph object.
#' @param x the output from the function plot.
#' @author Deisy Morselli Gysi <deisy at bioinf.uni-leipzig.de>
#' @importFrom igraph graph_from_data_frame
#' @return the CoDiNA plot as an igraph object.
#' @export
#'
#' @examples

#' Nodes = LETTERS[1:10]
#' Net1 = data.frame(Node.1 = sample(Nodes) , Node.2 = sample(Nodes), wTO = runif(10,-1,1))
#' Net2 = data.frame(Node.1 = sample(Nodes) , Node.2 = sample(Nodes), wTO = runif(10,-1,1))
#' Net3 = data.frame(Node.1 = sample(Nodes) , Node.2 = sample(Nodes), wTO = runif(10,-1,1))
#' DiffNet = MakeDiffNet (Data = list(Net1,Net2,Net3), Code = c('Net1', 'Net2', 'Net3') )

#' Graph = plot(x = DiffNet,
#'  layout = NULL, smooth.edges = TRUE,
#'  path = NULL, MakeGroups = FALSE, Cluster = FALSE,
#'  legend = TRUE, manipulation = FALSE, sort.by.Phi = FALSE)
#' x = as.igraph(Graph)
#'
#' plot(x)

#'
as.igraph <- function(x){
  stopifnot(any(class(x) %in% 'CoDiNA.plot'))
  x =  igraph::graph_from_data_frame(x$network$x$edges, directed = FALSE, vertices = x$network$x$nodes)
  igraph::V(x)$shape = 'circle'
  igraph::V(x)$size = igraph::V(x)$size
  return(x)
}

