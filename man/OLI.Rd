\name{OLI}
\alias{OLI}
\docType{data}
\title{
OLI
}
\description{
This data.table contains node and the weighted topological overlap (wTO) of Transcription Factors (TFs), from GSE4290 (Sun, 2006) for 50 brain samples with oligodendrogliomas. The wTO was calculated using the package wTO.
}
\usage{data("OLI")}
\format{
  A data frame with 3488761 observations on the following 3 variables.
  \describe{
    \item{\code{Node.1}}{a factor with levels. TF names }
    \item{\code{Node.2}}{a factor with levels. TF names }
    \item{\code{cor}}{a numeric vector. wTO values calculated using only the TFs}
  }
}

\source{
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse4290
}
\references{
Sun L, Hui AM, Su Q, Vortmeyer A et al. Neuronal and glioma-derived stem cell factor induces angiogenesis within the brain. Cancer Cell 2006 Apr;9(4):287-300. PMID: 16616334
}
\examples{
data(OLI)
str(OLI)
}
\keyword{datasets}
