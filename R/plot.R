#' @title  plot.CoDiNA
#' @aliases plot.CoDiNA
#' @description Categorize the Nodes into Phi and groups categories. Also, creates an interactive view of the CoDiNA network.
#' @param x Output from MakeDiffNet
#' @param distance distance to be used can be: Group, all, Phi or Phi_tilda
#' @param cutoff cutoff to be applied on the distance. By default, this number is 0.33.
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
#' @importFrom magrittr "%>%"
#' @export
#' @export plot.CoDiNA
#' @examples
#' set.seed(123)
#' Nodes = LETTERS[1:10]
#' Net1 = data.frame(Node.1 = sample(Nodes) , Node.2 = sample(Nodes), wTO = runif(10,-1,1))
#' Net2 = data.frame(Node.1 = sample(Nodes) , Node.2 = sample(Nodes), wTO = runif(10,-1,1))
#' Net3 = data.frame(Node.1 = sample(Nodes) , Node.2 = sample(Nodes), wTO = runif(10,-1,1))
#' DiffNet = makeDiffNet (x = list(Net1,Net2,Net3), Code = c('Net1', 'Net2', 'Net3') )

#' Graph = plot(x = DiffNet,
#'  cutoff = 0.3, layout = NULL, smooth.edges = TRUE,
#'  path = NULL, MakeGroups = FALSE, Cluster = FALSE,
#'  legend = TRUE, manipulation = FALSE, sort.by.Phi = FALSE)
#' Graph
#'

plot.CoDiNA = function(x, distance = 'Phi_tilda', cutoff = 0.33,
                        layout = NULL, smooth.edges = TRUE,
                        path = NULL, MakeGroups = FALSE,
                        Cluster = FALSE, legend = TRUE,
                        manipulation = FALSE,
                        sort.by.Phi = FALSE , ...)
{
  `%ni%` <- Negate(`%in%`)
  `%>%` <- magrittr::`%>%`

  Vars = c('Node.1', 'Node.2', paste('distance',distance, sep = '_'), 'Group', 'Phi', 'Phi_tilda')
  # message(Vars)
  if (any(Vars %ni% names(x))) {
    stop("x input is not complete.")
  }
  distance = by =  paste('distance',distance, sep = '_')

  input_vis = data.frame(subset(x,!is.na(x$Node.1) ,select = Vars))

  if (is.numeric(cutoff) == FALSE) {
    stop("cutoff value must be numeric.")
  }
  if (Cluster %ni% c(TRUE, FALSE)) {
    stop("Cluster must be T / F.")
  }
  if (smooth.edges %ni% c(TRUE, FALSE)) {
    stop("smooth.edges must be T / F.")
  }
  if(cutoff<0.01){
    input_vis = droplevels(subset(input_vis, abs(input_vis$distance) > 0.01))
  }
  if(cutoff>0.01){
    input_vis = droplevels(subset(input_vis,
                                  abs(input_vis$distance) > cutoff))
  }
  if (nrow(input_vis) <= 2) {
    stop("Not enough nodes on your network. Choose a lower cutoff.")
  }
  if (smooth.edges == TRUE) {
    smooth.edges = "enabled"
  }

  input_vis = input_vis[!is.na(input_vis$distance), ]
  ####
  #### getting colors for Phi groups
  ####
  input_vis$GROUP_FULL = as.factor(apply(cbind(as.character(input_vis$Phi),
                                               as.character(input_vis$Group)),
                                         1, paste, collapse = '.'))

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
  gg = igraph::graph_from_data_frame(data.frame(input_vis$Node.1, input_vis$Node.2, weights = input_vis$distance), directed = FALSE)
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

  igraph::E(gg)$weight = abs(input_vis$distance)
  names(DEGREE) = "degree"
  DEGREE$id = row.names(DEGREE)
  nodes = suppressMessages(plyr::join(nodes, DEGREE))
  nodes$value = (nodes$degree - min(nodes$degree))/(max(nodes$degree) -
                                                      min(nodes$degree))
  nodes$value = nodes$value^2 + 5
  igraph::V(gg)$size = nodes$value
  igraph::V(gg)$color = nodes$cluster

  nodes$size = nodes$value

  Df = rbind(data.frame(ID = input_vis$Node.1, Phi = input_vis$Phi),
             data.frame(ID = input_vis$Node.2, Phi = input_vis$Phi))

  Map2 = suppressMessages(data.table::dcast(data = Df, ID ~ Phi))
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


  Map = data.frame(table(Df$ID, Df$Phi))
  Z = suppressMessages(aggregate(Map$Freq, by=list(Category=Map$Var1), FUN=max))
  names(Z) = c('Var1', 'Freq')
  Group_Node = suppressMessages(plyr::join(Z, Map, match = 'first'))
  Group_Node = data.frame(id = Group_Node$Var1, groupPhi = Group_Node$Var2)
  nodes = suppressMessages(plyr::join(nodes, Group_Node))


  QUI2 = data.frame(Degree_a= nodes$a,   Degree_b= nodes$b,  Degree_g= nodes$g)
  V = suppressWarnings(apply(QUI2, 1, chisq.test, p = c(1/3,1/3,1/3)))

  for(no in 1:length(V)){
    # nodes$p.value_phi[no] = (U[[no]]$p.value)
    nodes$p.value_phi[no] = (V[[no]]$p.value)
  }
  nodes$groupPhi = as.character(nodes$groupPhi)
  # nodes$p.value.adj = p.adjust(nodes$p.value_phi, method = 'BH')
  nodes$groupPhi = ifelse(nodes$p.value_phi > 0.05, 'U', nodes$groupPhi)

  Df = rbind(data.frame(ID = input_vis$Node.1, Phi = input_vis$GROUP_FULL),
             data.frame(ID = input_vis$Node.2, Phi = input_vis$GROUP_FULL))

  Map = data.frame(table(Df$ID, Df$Phi))
  Map2 = suppressMessages(data.table::dcast(data = Df, ID ~ Phi))
  names(Map2)[1] = 'id'

  Z = suppressMessages(aggregate(Map$Freq, by=list(Category=Map$Var1), FUN=max))
  names(Z) = c('Var1', 'Freq')
  Group_Node = suppressMessages(plyr::join(Z, Map, match = 'first'))
  Group_Node = data.frame(id = Group_Node$Var1, group = Group_Node$Var2)
  nodes = suppressMessages(plyr::join(nodes, Group_Node))

  QUI2 = data.frame(Map2[,-1])
  V = suppressWarnings(apply(QUI2, 1, chisq.test))

  for(no in 1:length(V)){
    nodes$p.value_group[no] = (V[[no]]$p.value)
  }

  nodes$group = as.character(nodes$group)
  nodes$group = ifelse(nodes$p.value_group > 0.05, 'U', nodes$group)

  nodes$GROUP_FULL = nodes$group
  nodes = suppressMessages(plyr::join(nodes, colormap))
  nodes$color = ifelse(nodes$groupPhi == 'U',  'grey50', nodes$color)
  nodes$shape = ifelse(nodes$groupPhi == 'U', 'diamond', nodes$shape)

  # nodes$ID = nodes$id
  nodes$label = nodes$id
  nodes$title = paste0("<p> Node ID: ", nodes$label,
                       "<br>Degree: ", nodes$degree,
                       "<br>Degree a: ", nodes$a,
                       "<br>Degree b: ", nodes$b,
                       "<br>Degree g: ", nodes$g,
                       '<br>Phi:', nodes$groupPhi,
                       '<br>Group:', nodes$group,
                       "</p>")
  if(sort.by.Phi == TRUE){
    nodes$id = paste(nodes$groupPhi, nodes$group, nodes$id, sep ='_')

  }

  nodes$groupPhi <- nodes$Phi <- as.character(nodes$groupPhi)

  nodes$color = ifelse(nodes$GROUP_FULL == 'U' & nodes$Phi == 'a',  'green', nodes$color)
  nodes$color = ifelse(nodes$GROUP_FULL == 'U' & nodes$Phi == 'b',  'red', nodes$color)
  nodes$color = ifelse(nodes$GROUP_FULL == 'U' & nodes$Phi == 'g',  'blue', nodes$color)

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
                      group = edges_ID$Group,
                      distance = edges_ID$distance,
                      Phi = edges_ID$Phi, color = edges_ID$color)
  wto = abs(edges_ID$distance)
  edges$width = 3*abs((wto - min(wto))/(max(wto) -
                                          min(wto)))^4
  nodes = droplevels(nodes)
  edges = droplevels(edges)
  ledges <- data.frame(color = colormap$color,
                       label = colormap$GROUP_FULL,
                       arrows = '', shape = 'elipse')
  ledges2 <- data.frame(color = c('green', 'red', 'blue', 'grey'),
                        label = c('a', 'b', 'g', 'U'),
                        arrows = c("", "", '', ''),
                        shape = c('triangle', 'square', 'star', 'diamond'))
  ledges2 = rbind(ledges, ledges2)
  edges$title = paste0("<p>Edge: ", edges$Label,
                       "<br>Distance: ", round(edges$distance, 2),
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
                           selectedBy = list(variable = 'group', multiple = FALSE),
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
                        group = nodes$group,
                        Phi = nodes$groupPhi,
                        Degree_Total = nodes$degree,
                        Degree_a= nodes$a,
                        Degree_b= nodes$b,
                        Degree_g= nodes$g,
                        p.value_phi = nodes$p.value_phi,
                        p.value_group = nodes$p.value_group)
  edgesout = data.frame(id1 = edges$L1, id2 = edges$L2, Group = edges$group,
                        Phi = edges$Phi, distance = edges$distance)
  U = (list(Nodes = nodesout, Edges = edgesout, network = network))
  class(U)<- append('CoDiNA.plot', class(U))
  return(U)
}

