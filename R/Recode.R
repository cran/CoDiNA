#' Title Recode
#'
#' @param Phi Phi categories
#' @param values Categorical base
#' @importFrom magrittr "%>%" "%<>%"
#'
#' @keywords internal
#' 
Recode = function(Phi, values){
  Phi = Phi
  values = values
  
  Beta_in = matrix(NA, nrow = nrow(values), ncol = length(values)-1)
  Gama_in = matrix(NA, nrow = nrow(values), ncol = length(values))
  
  for ( i in 1:ncol(values)){
    Gama_in[,i] = ifelse(values[,i]!=0, names(values)[i], '')
  }
  
  for ( i in 2:ncol(values)){
    Beta_in[,i-1] = ifelse(values[,1] != values[,i], names(values)[i], '')
  }
  g1 = apply(Gama_in, 1, function(x) paste(x, collapse = '.'))
  b1 = apply(Beta_in, 1, function(x) paste(x, collapse = '.'))
  
  Phi2 = NULL
  Phi2 = ifelse(Phi == 'a', 'a', ifelse(Phi == 'g', paste('g', g1, sep = '.'),
                                        paste('b', b1, sep = '.')))
  
  Phi2 = gsub(pattern = '\\.{2,}', replacement = '.', x = Phi2) 
  Phi2 = gsub(pattern = '\\.$', replacement = '', x = Phi2)
  return(Phi2)
}

