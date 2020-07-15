#' @title  plot.CoDiNA
#' @aliases plot.CoDiNA
#' @description Categorize the Nodes into Phi and groups categories. Also, creates an interactive view of the CoDiNA network.
#' @param x Output from MakeDiffNet
#' @param cutoff.external The cut-off between the clusters (delta from the center to the edge coordinates), the closer to 1, the better.
#' @param cutoff.internal The cut-off inside the clusters (delta from the theoretical cluster to the edge coordinates), the closer to zero, the better.
#' @param cutoff.ratio The cut-off for the ratio of both scores. Default is set to 1. The greater, the better.

#' @param layout a layout from the igraph package.
#' @param smooth.edges If the edges should be smoothed or not.
#' @param sort.by.Phi if the graph should be plotted in the Phi order
#' @param path If the graph should be saved specify the name of the file.
#' @param Cluster TRUE or FALSE if the nodes should be clustered (double click to uncluster).
#' @param MakeGroups algorithm to find clusters. One of the followings: walktrap, optimal, spinglass, edge.betweenness, fast_greedy, infomap, louvain, label_prop, leading_eigen. Default to FALSE.
#' @param legend TRUE or FALSE if the legend should appear.
#' @param manipulation TRUE or FALSE if the graph should be editable.
#' @param \dots Additional plotting parameters.

#' @author Deisy Morselli Gysi <deisy at bioinf.uni-leipzig.de>
#' @return Returns a list contatining: The nodes description, the Edges description and the network graph.
#' @method plot CoDiNA
#' @importFrom grDevices colorRampPalette x11
#' @importFrom graphics plot
#' @importFrom stats aggregate kmeans chisq.test
#' @importFrom visNetwork visNetwork visInteraction visEdges visOptions visClusteringByGroup visLegend visPhysics visIgraphLayout visOptions visSave
#' @importFrom plyr arrange join join_all
#' @importFrom igraph graph_from_data_frame degree E plot.igraph
#' @importFrom data.table as.data.table
#' @importFrom magrittr "%>%" "%<>%"
#' @importFrom reshape2 dcast
#' @export
#' @export plot.CoDiNA
#' @examples
#' suppressWarnings(RNGversion("3.5.0"))
#'
#' Nodes = LETTERS[1:10]
#' Net1 = data.frame(Node.1 = sample(Nodes) , Node.2 = sample(Nodes), wTO = runif(10,-1,1))
#' Net2 = data.frame(Node.1 = sample(Nodes) , Node.2 = sample(Nodes), wTO = runif(10,-1,1))
#' Net3 = data.frame(Node.1 = sample(Nodes) , Node.2 = sample(Nodes), wTO = runif(10,-1,1))
#' DiffNet = MakeDiffNet (Data = list(Net1,Net2,Net3), Code = c('Net1', 'Net2', 'Net3') )

#' Graph = plot(x = DiffNet,
#'  layout = NULL, smooth.edges = TRUE,
#'  path = NULL, MakeGroups = FALSE, Cluster = FALSE,
#'  legend = TRUE, manipulation = FALSE, sort.by.Phi = FALSE)
#' Graph
#'

plot.CoDiNA = function(x, cutoff.external = 0, cutoff.internal = 1, cutoff.ratio = 1,
                        layout = NULL, smooth.edges = TRUE,
                        path = NULL, MakeGroups = FALSE,
                        Cluster = FALSE, legend = TRUE,
                        manipulation = FALSE,
                        sort.by.Phi = FALSE , ...)
{
  `%ni%` <- Negate(`%in%`)
  `%>%` <- magrittr::`%>%`
  `%<>%` <- magrittr::`%<>%`
  clean = subset(x, x$Score_Phi_tilde > cutoff.external &
                   x$Score_internal < cutoff.internal &
                   x$Score_ratio > cutoff.ratio)

  Vars = c('Node.1', 'Node.2', 'Score_Phi_tilde',
           'Score_internal', 'Phi', 'Phi_tilde', 'Score_ratio')
  # message(Vars)
  if (any(Vars %ni% names(x))) {
    stop("x input is not complete.")
  }
  Score = by =  'Score_Phi_tilde'

  input_vis = data.frame(subset(clean,!is.na(clean$Node.1) ,select = Vars))


  if (is.numeric(cutoff.external) == FALSE | is.numeric(cutoff.internal) == FALSE) {
    stop("cutoff value must be numeric.")
  }
  if (Cluster %ni% c(TRUE, FALSE)) {
    stop("Cluster must be TRUE or FALSE.")
  }
  if (smooth.edges %ni% c(TRUE, FALSE)) {
    stop("smooth.edges must be TRUE or FALSE.")
  }
  if(cutoff.external<0.01){
      input_vis = droplevels(subset(input_vis, abs(input_vis$Score_Phi_tilde) > 0.01))
  }
  if (nrow(input_vis) <= 2) {
    stop("Not enough nodes on your network. Choose a another cutoff.")
  }
  if (smooth.edges == TRUE) {
    smooth.edges = "enabled"
  }

  Nodes_Phi = ClusterNodes(DiffNet = x, cutoff.external = cutoff.external, cutoff.internal = cutoff.internal)

  input_vis = input_vis[!is.na(input_vis$Score_Phi_tilde), ]
  input_vis = input_vis[!is.na(input_vis$Score_internal), ]

  ####
  #### getting colors for Phi groups
  ####
  input_vis$GROUP_FULL = as.factor(input_vis$Phi_tilde)


  colsC = data.frame(GROUP_FULL = levels(droplevels(subset(input_vis$GROUP_FULL, input_vis$Phi == 'a'))))
  colsS = data.frame(GROUP_FULL = levels(droplevels(subset(input_vis$GROUP_FULL, input_vis$Phi == 'g'))))
  colsD = data.frame(GROUP_FULL = levels(droplevels(subset(input_vis$GROUP_FULL, input_vis$Phi == 'b'))))



  if(nrow(colsD) == 0){
    colsD = data.frame(GROUP_FULL = NULL,
                       Phi = NULL, color = NULL, shape = NULL)
  }
  if(nrow(colsS) == 0){
    colsD = data.frame(GROUP_FULL = NULL, Phi = NULL,color = NULL, shape = NULL)
  }
  if(nrow(colsC) == 0){
    colsC = data.frame(GROUP_FULL = NULL, Phi = NULL,color = NULL, shape = NULL)
  }
  if(nrow(colsC) > 0){
    colsC$color = colorRampPalette(c('limegreen', 'olivedrab1'))(nrow(colsC))
    colsC$shape = 'triangle'
    colsC$Phi = 'a'
  }
  if(nrow(colsS) > 0){
    colsS$color = colorRampPalette(c('lightskyblue', 'royalblue'))(nrow(colsS))
    colsS$shape = 'star'
    colsS$Phi = 'g'
  }
  if(nrow(colsD) > 0){
    colsD$color = colorRampPalette(c('salmon','violetred3'))(nrow(colsD))
    colsD$shape = 'square'
    colsD$Phi = 'b'
  }

  colsI = data.frame(GROUP_FULL = 'U', color = '#bdbdbd',
                     shape = 'diamond', Phi = 'U')

  colormap = rbind(colsC, colsD, colsS, colsI)
  input_vis = suppressMessages(plyr::join(input_vis, colormap))

  input_vis = suppressMessages(plyr::arrange(input_vis, input_vis$Node.1, input_vis$Node.2))
  input_vis = droplevels(input_vis)
  nodes <- data.frame(id = sort(unique(c(as.character(input_vis$Node.1),
                                         as.character(input_vis$Node.2)))))
  gg = igraph::graph_from_data_frame(data.frame(input_vis$Node.1, input_vis$Node.2, weights = input_vis$Score_Phi_tilde), directed = FALSE)

  DEGREE = as.data.frame(igraph::degree(gg))
  if(MakeGroups == FALSE){
    group = 1
  }
  if (MakeGroups == 'infomap'){

    group = igraph::cluster_infomap(gg)$membership
  }
  else if (MakeGroups == 'walktrap'){
    group = igraph::cluster_walktrap(gg)$membership
  }
  else if (MakeGroups == 'leading_eigen'){
    group = igraph::cluster_leading_eigen(gg)$membership
  }
  else if (MakeGroups == 'louvain'){
    group = igraph::cluster_louvain(gg)$membership
  }
  else if (MakeGroups == 'label_prop'){
    group = igraph::cluster_label_prop(gg)$membership
  }
  else if (MakeGroups == 'fast_greedy'){
    group = igraph::cluster_fast_greedy(gg)$membership
  }
  else if (MakeGroups == 'optimal'){
    group = igraph::cluster_optimal(gg)$membership
  }
  else if (MakeGroups == 'spinglass'){
    group = igraph::cluster_spinglass(gg)$membership
  }
  else if (MakeGroups == 'edge.betweenness'){
    group = igraph::edge.betweenness.community(gg)$membership
  }
  nodes = suppressMessages(plyr::join(nodes, data.frame(id = igraph::V(gg)$name, cluster = group)))

  igraph::E(gg)$weight = abs(input_vis$Score_Phi_tilde)
  names(DEGREE) = "degree"
  DEGREE$id = row.names(DEGREE)
  nodes = suppressMessages(plyr::join(nodes, DEGREE))
  nodes$value = (nodes$degree - min(nodes$degree))/(max(nodes$degree) -
                                                      min(nodes$degree))
  nodes$value = nodes$value^2 + 5
  igraph::V(gg)$size = nodes$value
  igraph::V(gg)$color = nodes$cluster

  nodes$size = nodes$value
  names(Nodes_Phi) [1]= 'id'
  nodes = suppressWarnings(suppressMessages( plyr::join(nodes, Nodes_Phi)))

  Nodes1 = data.frame(id = clean$Node.1, Phi_tilde = clean$Phi_tilde, Phi = clean$Phi)
  Nodes2 = data.frame(id = clean$Node.2, Phi_tilde = clean$Phi_tilde, Phi = clean$Phi)
  Nodes = rbind(Nodes1, Nodes2)


  Map2 =  reshape2::dcast(Nodes, id~Phi, fun.aggregate = length,
                                value.var = 'Phi')
  rm(Nodes1)
  rm(Nodes2)
  rm(Nodes)



  names(Map2)[1] = 'id'
  if(is.null(Map2$a)){
    Map2$a = 0
  }
  if(is.null(Map2$g)){
    Map2$g = 0
  }
  if(is.null(Map2$b)){
    Map2$b = 0
  }

  nodes = suppressMessages(plyr::join(nodes, Map2))


  nodes$label = nodes$id
  nodes$title = paste0("<p> Node ID: ", nodes$label,
                       "<br>Degree: ", nodes$degree,
                       "<br>Degree a: ", nodes$a,
                       "<br>Degree b: ", nodes$b,
                       "<br>Degree g: ", nodes$g,
                       '<br>Phi:', nodes$Phi,
                       '<br>Phi_tilde:', nodes$Phi_tilde,
                       "</p>")
  if(sort.by.Phi == TRUE){
    nodes$id = paste(nodes$groupPhi, nodes$group, nodes$id, sep ='_')

  }

  nodes$GROUP_FULL <- nodes$Phi <- as.character(nodes$Phi)

  nodes$color = ifelse(nodes$Phi == 'a',  'green', 'gray30')
  nodes$color = ifelse(nodes$Phi == 'b',  'red', nodes$color)
  nodes$color = ifelse(nodes$Phi == 'g',  'blue', nodes$color)
  nodes$frame.color = nodes$color

  nodes = nodes[order(nodes$id),]
  node1_ID =  data.frame(Node.1 = nodes$label, from = nodes$id)
  node2_ID = data.frame(Node.2 = nodes$label, to = nodes$id)
  # input2 = data.frame(Node.1 = input_vis$Node.1, Node.2= input_vis$Node.2)
  edges_ID = suppressMessages(plyr::join_all(list(input_vis, node1_ID, node2_ID) ))
  edges_ID$label = apply(edges_ID[,1:2], 1, paste, collapse = '<->')
  edges_ID$L1 = edges_ID[,1]
  edges_ID$L2 = edges_ID[,2]
  edges <- data.frame(from = edges_ID$from, to = edges_ID$to,
                      Label = edges_ID$label, L1 = edges_ID$L1, L2 = edges_ID$L2,
                      group = edges_ID$Phi_tilde,
                      Score = edges_ID$Score_Phi_tilde,
                      Phi = edges_ID$Phi, color = edges_ID$color)
  wto = abs(edges_ID$Score_ratio)
  edges$width = 3*abs((wto - min(wto))/(max(wto) -
                                          min(wto)))^4
  nodes = droplevels(nodes)
  edges = droplevels(edges)
  ledges <- data.frame(color = colormap$color,
                       label = colormap$GROUP_FULL,
                       arrows = '')
  ledges2 <- data.frame(color = c('green', 'red', 'blue', 'grey'),
                        label = c('a', 'b', 'g', 'U'),
                        arrows = c("", "", '', ''))
  ledges2 = rbind(ledges, ledges2)
  edges$title = paste0("<p>Edge: ", edges$Label,
                       "<br>Score: ", round(edges$Score, 2),
                       '<br>Group:', edges$group,
                       '<br>Phi:', edges$Phi,"</p>")

  NODESIZE = length(nodes$id)
  EDGESIZE = length(edges$Label)
  main = paste('Network contains:', NODESIZE, 'nodes and', EDGESIZE, 'edges.')

  nodes$color.border= nodes$color
  network <- visNetwork::visNetwork(nodes, edges, main = main) %>%

    visNetwork::visInteraction(navigationButtons = TRUE, hover = TRUE, multiselect = FALSE) %>%
    # visNetwork::visEdges(smooth = smooth.edges )%>%
    visNetwork::visOptions(highlightNearest = list(enabled = TRUE,
                                                   degree = 1, hover = FALSE),
                           nodesIdSelection = list(enabled = TRUE, useLabels =TRUE,
                                                   style = "width: 200px; height: 26px;\\n   background: #f8f8f8;\\n   color: darkblue;\\n   border:none;\\n   outline:none;"),
                           manipulation = F,
                           # selectedBy = list(variable = 'cluster', multiple = FALSE),
                           selectedBy = list(variable = 'Phi_tilde', multiple = FALSE),
                           # collapse = list(enabled = TRUE, clusterOptions =list(Phi = nodes$groupPhi),resetHighlight = TRUE)) %>%
    ) %>%
    visNetwork::visPhysics(enabled = F) %>%
    visNetwork::visExport(type = "png",
                          name = "networkpng",
                          float = "right",
                          label = "Save png",
                          background = "transparent",
                          style= "")

  if (Cluster == T) {
    network <- network %>%
      visNetwork::visClusteringByGroup(groups = unique((nodes$group)))
  }
  if (legend == T) {

    network <- network %>%
      visNetwork::visLegend(width = 0.3, useGroups = FALSE,
                            position = "right", main = "Group", addNodes = ledges2,
                            ncol = 1)

  }
  if (!is.null(layout)) {
    network <- network %>% visNetwork::visIgraphLayout(layout = layout)
  }
  if (manipulation == T) {
    network <- network %>% visNetwork::visOptions(manipulation = TRUE)
  }
  if (is.null(path)) {
    network
  }
  else if (!is.null(path)) {
    visNetwork::visSave(network, file = path)
    message(path)
  }



  nodesout = data.frame(id = nodes$label,
                        Phi_tilde = nodes$Phi_tilde,
                        Phi = nodes$Phi,
                        Degree_Total = nodes$degree,
                        Degree_a= nodes$a,
                        Degree_b= nodes$b,
                        Degree_g= nodes$g)

  edgesout = data.frame(id1 = edges$L1, id2 = edges$L2, Group = edges$group,
                        Phi = edges$Phi, Score = edges$Score)
  U = (list(Nodes = nodesout, Edges = edgesout, network = network))
  class(U)<- append('CoDiNA.plot', class(U))
  return(U)
}

