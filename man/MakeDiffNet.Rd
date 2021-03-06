% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeDiffNet.R
\name{MakeDiffNet}
\alias{MakeDiffNet}
\title{MakeDiffNet}
\usage{
MakeDiffNet(Data, Code, cutoff = 0.33, stretch = TRUE)
}
\arguments{
\item{Data}{List of data.frames containig Node.1, Node.2 and the correlation value}

\item{Code}{Name of each one of the networks.}

\item{cutoff}{By default, the cutoff is 0.33. If the user wants to use another value, it has to be cited on the description of the used methodology that the cutoff was changed.}

\item{stretch}{Should the input data be normalized? Default to TRUE.}
}
\value{
Returns a data.table contating: Nodes names, correlation value for each network (the input values), the k means cluster that link belongs, the Phi groups (Phi and Phi tilda), the signed group that link belongs to, the unsigned group. The distance to the center, and the distance normalized by: Phi_tilda, Phi, signed group or all data.
}
\description{
Categorize links into Phi categories, calculate the distance to the center and also normlize the distance into some categories: Phi and Phi tilda, group and all.
}
\examples{
suppressWarnings(RNGversion("3.5.0"))
Nodes = LETTERS[1:20]
Net1 = data.frame(Node.1 = sample(Nodes) , Node.2 = sample(Nodes), wTO = runif(10,-1,1))
Net2 = data.frame(Node.1 = sample(Nodes) , Node.2 = sample(Nodes), wTO = runif(10,-1,1))
Net3 = data.frame(Node.1 = sample(Nodes) , Node.2 = sample(Nodes), wTO = runif(10,-1,1))
DiffNet = MakeDiffNet (Data = list(Net1,Net2,Net3), Code = c('Net1', 'Net2', 'Net3') )
print(DiffNet)
}
\author{
Deisy Morselli Gysi <deisy at bioinf.uni-leipzig.de>
}
