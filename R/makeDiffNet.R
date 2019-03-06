#' @title MakeDiffNet
#' @description Categorize links into Phi categories, calculate the distance to the center and also normlize the distance into some categories: Phi and Phi tilda, group and all.
#' @param Data List of data.frames containig Node.1, Node.2 and the correlation value
#' @param Code Name of each one of the networks.
#' @param stretch Should the input data be normalized? Default to TRUE.
#' @param cutoff By default, the cutoff is 0.33. If the user wants to use another value, it has to be cited on the description of the used methodology that the cutoff was changed.
#' @author Deisy Morselli Gysi <deisy at bioinf.uni-leipzig.de>
#' @return Returns a data.table contating: Nodes names, correlation value for each network (the input values), the k means cluster that link belongs, the Phi groups (Phi and Phi tilda), the signed group that link belongs to, the unsigned group. The distance to the center, and the distance normalized by: Phi_tilda, Phi, signed group or all data.
#' @export
#' @import data.table
#' @importFrom stats aggregate kmeans
#'
#' @examples
#' suppressWarnings(RNGversion("3.5.0"))
#' Nodes = LETTERS[1:20]
#' Net1 = data.frame(Node.1 = sample(Nodes) , Node.2 = sample(Nodes), wTO = runif(10,-1,1))
#' Net2 = data.frame(Node.1 = sample(Nodes) , Node.2 = sample(Nodes), wTO = runif(10,-1,1))
#' Net3 = data.frame(Node.1 = sample(Nodes) , Node.2 = sample(Nodes), wTO = runif(10,-1,1))

#' DiffNet = MakeDiffNet (Data = list(Net1,Net2,Net3), Code = c('Net1', 'Net2', 'Net3') )
#' print(DiffNet)

MakeDiffNet = function(Data, Code, cutoff = 0.33,
                       stretch = TRUE){
  `:=` <- data.table::`:=`

  BASE = CreateFullBase(x = Data, Code = Code)
  Nodes = BASE[[2]]
  if(length(Nodes$Nodes)< 2){
    stop('Not enough Nodes. Please use a higher cutoff.')
  }
  BASE = BASE[[1]]
  BASE = subset(BASE, (BASE$Node.1 %in% Nodes$Nodes & BASE$Node.2 %in% Nodes$Nodes ))
  if(stretch == TRUE){
    mini_base =  BASE[,3:ncol(BASE)]
    mini_base = apply(mini_base, 2, function(x) normalize(abs(x))*ifelse(x <0,-1,1))
    BASE = cbind.data.frame(BASE[,1:2], mini_base)
    rm(mini_base)
  }
  CAT_BASE = Categorize(M = BASE[,3:ncol(BASE)])
  Stay =rowSums(abs(CAT_BASE))
  CAT_BASE = subset(CAT_BASE, Stay>0)
  BASE = subset(BASE, Stay>0)
  message(paste('Total of edges (inside the cutoff):',nrow(BASE)) )
  Z = apply(CAT_BASE, 1, paste, collapse = '.')


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

  n = ncol(CAT_BASE)
  y$Phi2<-NA
  #### Assumes the first code as ref level

    y$Phi2 = Recode(Phi = y$Phi, values = CAT_BASE)

  #### Now I want to put the meaning in the code:
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
  message('Recode is done!')
  y$Group = Names


  y$dist_center = sqrt(apply(BASE[,3:ncol(BASE)]^2, 1, sum))/ sqrt(rowSums(abs(CAT_BASE)))# sqrt(rowSums(abs(round(BASE[,-c(1:2)],1))))
  Bla = (BASE[,3:ncol(BASE)] - CAT_BASE)^2

  y$distance_to_group = sqrt(apply(Bla, 1, sum))/sqrt(rowSums(abs(CAT_BASE)))


  y$distance_to_group[y$distance_to_group > 1]<- 1


  y$dist_center[rowSums(abs(CAT_BASE)) == 0]<- 0
  y$dist_center[y$dist_center > 1]<- 1
  y$sum = y$sum_abs = y$Z = y$pos=  y$neg = y$min = NULL

  # Phi_gp = y$Phi
  # GP = y$Group

  y = data.table::data.table(y)


  y[,distance_Phi:= normalize(dist_center), by = y$Phi]
  y[,distance_Phi2:= normalize(dist_center), by = y$Phi2]
  y$distance_all = normalize(y$dist_center)
  y$distance_center = y$dist_center
  distance_Phi_tilda =y$dist_center= dist_center = y$Phi_tilda = distance_Phi= distance_Phi2= distance_Group = Group = Phi_tilda = NULL

  y$Phi_tilde= y$Phi2
  y$distance_Phi_tilde = y$distance_Phi2

  y$Phi2 = y$distance_Phi2 = NULL

  y$distance_internal=y$distance_to_group


  y$Score_center = y$distance_center
  y$Score_Phi = y$distance_Phi
  y$Score_Phi_tilde = y$distance_Phi_tilde
  y$Score_internal = y$distance_internal
  y$Score_ratio = y$Score_Phi_tilde/y$Score_internal

  y = subset(y, select = c('Node.1', 'Node.2', Code,
                           'Phi', 'Phi_tilde',
                           'Group',
                           'Score_center',
                           'Score_Phi',
                           'Score_Phi_tilde',
                           'Score_internal', 'Score_ratio'))


  # class(y)<- append(class(y), 'data.frame')
  class(y)<- append('CoDiNA', class(y))
  # class(y)<- 'CoDiNA'
  return(y)
}
