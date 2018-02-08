#' @title CreateFullBase
#' @description Joins a set of data.frames, order the nodes names by it's smaller value.
#' @param x List of data.frames containig Node.1, Node.2 and the correlation value
#' @param Code Name of each one of the networks.
#' @return Returns a list contating: The nodes names and it's correlation values in all networks, 0 if this node is absent.
#' @importFrom plyr join_all
#' @keywords internal
#' @author Deisy Morselli Gysi <deisy at bioinf.uni-leipzig.de>
CreateFullBase = function(x, Code){
  message('Starting now.')
  FULL_Bases = list()
  NODES_Bases = list()
  for ( i in 1:length(x)){
    message(i)
    X = x[[i]]

      names(X)[1:2] = c('Node.1', 'Node.2')

    NODES_Bases[[i]] = data.frame(Nodes = unique(c(as.character(X$Node.1),
                                                   as.character(X$Node.2))))
    message(paste('Edges:', nrow(X)))
    if(nrow(X)>0){
      FULL_Bases[[i]] = OrderNames(X)
      names(FULL_Bases[[i]])[3] = as.character(Code[i])
    }
  }
  Bases = suppressMessages(plyr::join_all(FULL_Bases, type = "full"))
  Bases[is.na(Bases)] <- 0
  Nodes = suppressMessages(plyr::join_all(NODES_Bases, type = "inner"))
  if(any(duplicated(data.frame(Bases$Node.1, Bases$Node.2)))){
    warning('You have duplicated links. Results will be presented, but you should double check your input.')
  }
  message(paste('Total of nodes:', nrow(Nodes)))
  message(paste('Total of edges:', dim(Bases)[1]))
  return(list(Bases, Nodes))
}
