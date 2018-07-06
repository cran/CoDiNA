#' @title ClusterNodes
#' @aliases ClusterNodes
#' @description Categorize the Nodes into Phi and Phi tilde.
#'
#' @param DiffNet The Differential network from MakeDiffNet
#' @param cutoff.external The cut-off between the clusters (delta from the center to the edge coordinates), the closer to 1, the better.
#' @param cutoff.internal The cut-off inside the clusters (delta from the theoretical cluster to the edge coordinates), the closer to zero, the better.
#' @importFrom stats p.adjust
#' @importFrom magrittr "%>%" "%<>%"
#' @importFrom data.table dcast
#' @export
#'
#' @examples

#' DiffNet = MakeDiffNet (Data = list(CTR,  AST), Code = c('CTR', 'AST') )
#' Genes_Phi = ClusterNodes(DiffNet, cutoff.external = 0.5, cutoff.internal = 0.25)
#' table(Genes_Phi$Phi_tilde)
#'
ClusterNodes <- function(DiffNet, cutoff.external = 0.8, cutoff.internal = 0.5){
  clean = subset(DiffNet, DiffNet$Score_Phi_tilde > cutoff.external & DiffNet$Score_internal < cutoff.internal )

  if(nrow(clean)==0){
    stop('Please, choose different cut-offs. You have zero nodes to cluster.')
  }
  Nodes = data.frame(Node = clean$Node.1, Phi_tilde = clean$Phi_tilde, Phi = clean$Phi)
  Nodes2 = data.frame(Node = clean$Node.2, Phi_tilde = clean$Phi_tilde, Phi = clean$Phi)

  Nodes = rbind(Nodes, Nodes2)

  NodesPhi_tilde =  data.table::dcast(Nodes, Node~Phi_tilde, fun.aggregate = length,
                                 value.var = 'Phi_tilde')
  # message('NodesPhi_tilde')
  NodesPhi =  data.table::dcast(Nodes, Node~Phi, fun.aggregate = length,
                                value.var = 'Phi')
  Chi2 = 1
  Chi = 1


  if (ncol(NodesPhi) > 2){
    #p = apply(NodesPhi[,-1], 2, sum)#/sum(NodesPhi[,-1])
    p = rep((1/(ncol(NodesPhi)-1)), ncol(NodesPhi)-1)
    GP1 = apply(NodesPhi[,-1], 1, function(x){
      y = names(x)[x==max(x)]
      y = ifelse(length(y)>1, 'U', y)
      return(y)
    })
    Chi = apply(NodesPhi[,-1], 1, function(x){
      u = suppressMessages( suppressWarnings( chisq.test(x, p = p, rescale.p = TRUE)$p.value))
      u = ifelse(is.na(u), 0.99, u)
      return(u)
    })

  }

  if (ncol(NodesPhi_tilde) > 2){
    GP2 = apply(NodesPhi_tilde[,-1], 1, function(x){
      y = names(x)[x==max(x)]
      y = ifelse(length(y)>1, 'U', y)
      return(y)
    })
    #p = apply(NodesPhi_tilde[,-1], 2, sum)#/sum(NodesPhi_tilde[,-1])
    p = rep((1/(ncol(NodesPhi_tilde)-1)), ncol(NodesPhi_tilde)-1)
    Chi2 = apply(NodesPhi_tilde[,-1], 1, function(x){
      u = suppressWarnings(suppressMessages( chisq.test(x, p=p, rescale.p = TRUE)$p.value))
      u = ifelse(is.na(u), 0.99, u)
      return(u)
    })

  }


  NodesPhi_tilde$Phi_tilde = ifelse(p.adjust(Chi2, method = 'BH') < 0.05, GP2, 'U')
  NodesPhi$Phi  = ifelse(p.adjust(Chi, method = 'BH') < 0.05,  GP1, 'U')
  NodesPhi = data.frame(Node = NodesPhi$Node,
                        Phi = NodesPhi$Phi)

  NodesPhi_tilde = data.frame(Node = NodesPhi_tilde$Node,
                              Phi_tilde = NodesPhi_tilde$Phi_tilde
                             )

  Out = suppressMessages( plyr::join(NodesPhi, NodesPhi_tilde))

  gp = character(length = length(Out))
  for( i in 1:nrow(Out)){
    gp[i] = strsplit(as.character(Out$Phi_tilde), "[.]")[[i]][1]
  }

  Out$Phi = ifelse(Out$Phi == 'U' & Out$Phi_tilde != 'U', gp, as.character(Out$Phi))

  return(Out)
}

