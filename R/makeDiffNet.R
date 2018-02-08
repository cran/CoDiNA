#' @title makeDiffNet
#' @description Categorize links into Phi categories, calculate the distance to the center and also normlize the distance into some categories: Phi and Phi tilda, group and all.
#' @param x List of data.frames containig Node.1, Node.2 and the correlation value
#' @param Code Name of each one of the networks.
#' @param cutoff By default, the cutoff is 0.33. If the user wants to use another value, it has to be cited on the description of the used methodology that the cutoff was changed.
#' @author Deisy Morselli Gysi <deisy at bioinf.uni-leipzig.de>
#' @return Returns a data.table contating: Nodes names, correlation value for each network (the input values), the k means cluster that link belongs, the Phi groups (Phi and Phi tilda), the signed group that link belongs to, the unsigned group. The distance to the center, and the distance normalized by: Phi_tilda, Phi, signed group or all data.
#' @export
#' @import data.table
#' @importFrom stats aggregate kmeans
#' @examples
#' set.seed(123)
#' Nodes = LETTERS[1:10]
#' Net1 = data.frame(Node.1 = sample(Nodes) , Node.2 = sample(Nodes), wTO = runif(10,-1,1))
#' Net2 = data.frame(Node.1 = sample(Nodes) , Node.2 = sample(Nodes), wTO = runif(10,-1,1))
#' Net3 = data.frame(Node.1 = sample(Nodes) , Node.2 = sample(Nodes), wTO = runif(10,-1,1))

#' DiffNet = makeDiffNet (x = list(Net1,Net2,Net3), Code = c('Net1', 'Net2', 'Net3') )
#' print(DiffNet)

makeDiffNet = function(x, Code, cutoff = 0.33){
  `:=` <- data.table::`:=`

  BASE = CreateFullBase(x = x, Code = Code)
  Nodes = BASE[[2]]
  if(length(Nodes$Nodes)< 2){
    stop('Not enough Node. Please use a higher cutoff.')
  }
  BASE = BASE[[1]]
  BASE = subset(BASE, (BASE$Node.1 %in% Nodes$Nodes & BASE$Node.2 %in% Nodes$Nodes ))
  CAT_BASE = Categorize(M = BASE[,3:ncol(BASE)])
  Stay =rowSums(abs(CAT_BASE))
  CAT_BASE = subset(CAT_BASE, Stay>0)
  BASE = subset(BASE, Stay>0)
  Z = apply(CAT_BASE, 1, paste, collapse = '.')
  # length(table(Z))
  # message('K means!')
  # K = kmeans(CAT_BASE, centers = length(table(Z)))
  BASE_COMPLETE = data.frame(cbind(BASE, Z))

  BASE_COMPLETE$sum = apply(CAT_BASE, 1, sum)
  BASE_COMPLETE$sum_abs = apply(abs(CAT_BASE), 1, sum)
  neg = pos = CAT_BASE
  neg[neg == 1] <- 0
  pos[pos == -1] <- 0

  BASE_COMPLETE$pos = abs(rowSums(pos))
  BASE_COMPLETE$neg = abs(rowSums(neg))

  BASE_COMPLETE$min = abs(apply(cbind(BASE_COMPLETE$neg, BASE_COMPLETE$pos), 1, min, na.rm = T))
  y = PhiCategory(Base = BASE_COMPLETE, n = ncol(CAT_BASE))

  #### Now I want to put the meaning in the code:
  #### if -1.0.0 means: -Con
  #### if -1.-1.-1 means: -con-aut-alz
  message('Coding the groups.')
  CAT_NEW = CAT_BASE
  CAT_NEW <- matrix(names(CAT_BASE), nrow = nrow(CAT_BASE),
                    ncol = ncol(CAT_BASE), byrow = T)
  Y = CAT_NEW
  Y[CAT_BASE == 0]<- 'No'
  Y[CAT_BASE == 1]<- '+'
  Y[CAT_BASE == -1]<- '-'
  CAT_NEW =matrix(paste(t(Y), t(CAT_NEW), sep = ''),ncol(CAT_NEW))

  Names = apply(CAT_NEW, 2, paste, collapse = ',')
  # message('Recode is done!')
  y$Group = Names

  CAT_NEW = CAT_BASE
  CAT_NEW <- matrix(names(CAT_BASE), nrow = nrow(CAT_BASE),
                    ncol = ncol(CAT_BASE), byrow = T)
  Y = CAT_NEW
  Y[CAT_BASE == 0]<- 'No'
  Y[CAT_BASE == 1]<- ''
  Y[CAT_BASE == -1]<- ''
  CAT_NEW =matrix(paste(t(Y), t(CAT_NEW), sep = ''),ncol(CAT_NEW))

  Names = apply(CAT_NEW, 2, paste, collapse = ',')
  message('Recode is done!')
  y$Group_Simple = Names

  y$dist_center = sqrt(apply(BASE[,3:ncol(BASE)]^2, 1, sum))/sqrt(rowSums(abs(CAT_BASE)))
  y$dist_center[rowSums(abs(CAT_BASE)) == 0]<- 0
  y$dist_center[y$dist_center > 1]<- 1
  y$sum = y$sum_abs = y$Z = y$pos=  y$neg = y$min = NULL

  # Phi_gp = y$Phi
  # GP = y$Group

  y = data.table::data.table(y)

  y[,distance_Phi_tilda:= normalize(dist_center), by = y$Phi_tilda]
  y[,distance_Phi:= normalize(dist_center), by = y$Phi]
  y[,distance_Group:= normalize(dist_center), by = y$Group]
  y$distance_all = normalize(y$dist_center)
  distance_Phi_tilda = dist_center = Phi_tilda = distance_Phi = distance_Group = Group = NULL
  # class(y)<- append(class(y), 'data.frame')
  class(y)<- append('CoDiNA', class(y))
# class(y)<- 'CoDiNA'
  return(y)
}
