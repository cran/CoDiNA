% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summary.CoDiNA}
\alias{summary.CoDiNA}
\title{summary.CoDiNA}
\usage{
\method{summary}{CoDiNA}(object, ...)
}
\arguments{
\item{object}{Output from MakeDiffNet}

\item{\dots}{Additional plotting parameters.}
}
\value{
Returns a summary describing the network.
}
\description{
summary of the CoDiNA network.
}
\examples{
Nodes = LETTERS[1:10]
Net1 = data.frame(Node.1 = sample(Nodes) , Node.2 = sample(Nodes), wTO = runif(10,-1,1))
Net2 = data.frame(Node.1 = sample(Nodes) , Node.2 = sample(Nodes), wTO = runif(10,-1,1))
Net3 = data.frame(Node.1 = sample(Nodes) , Node.2 = sample(Nodes), wTO = runif(10,-1,1))
DiffNet = MakeDiffNet (Data = list(Net1,Net2,Net3), Code = c('Net1', 'Net2', 'Net3') )
summary(DiffNet)
}
\author{
Deisy Morselli Gysi <deisy at bioinf.uni-leipzig.de>
}
