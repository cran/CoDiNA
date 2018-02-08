#' @title  summary.CoDiNA
#' @aliases summary.CoDiNA
#' @description summary of the CoDiNA network.
#' @param object Output from MakeDiffNet
#' @param \dots Additional plotting parameters.
#' @author Deisy Morselli Gysi <deisy at bioinf.uni-leipzig.de>
#' @return Returns a summary describing the network.
#' @method summary CoDiNA
#' @export
#' @examples
#' set.seed(123)
#' Nodes = LETTERS[1:10]
#' Net1 = data.frame(Node.1 = sample(Nodes) , Node.2 = sample(Nodes), wTO = runif(10,-1,1))
#' Net2 = data.frame(Node.1 = sample(Nodes) , Node.2 = sample(Nodes), wTO = runif(10,-1,1))
#' Net3 = data.frame(Node.1 = sample(Nodes) , Node.2 = sample(Nodes), wTO = runif(10,-1,1))

#' DiffNet = makeDiffNet (x = list(Net1,Net2,Net3), Code = c('Net1', 'Net2', 'Net3') )
#' summary(DiffNet)

# summary <- function(object){
#   UseMethod('CoDiNA')
# }

summary.CoDiNA <- function(object,...) {
  x = object
  Nodes_n = length(unique(c(as.character(x$Node.1),
                            as.character(x$Node.2))))
  cat('Nodes', Nodes_n, "\n")
  cat('Links', nrow(x), "\n")
cat('\n')
  PHI = as.data.frame(table(x$Phi))
  names(PHI)[1] = 'Phi'

  PHI_tilda = as.data.frame(table(x$Phi_tilda))
  names(PHI_tilda)[1] = 'Phi_tilda'

  Group = as.data.frame(table(x$Group))
  names(Group)[1] = 'Group'

  # print(PHI)
  # print(PHI_tilda)
  # print(Group)
  return(list(Phi = PHI, Phi_tilda = PHI_tilda, Group = Group))
}
