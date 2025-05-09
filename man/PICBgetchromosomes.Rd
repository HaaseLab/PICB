% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PICBgetchromosomes.R
\name{PICBgetchromosomes}
\alias{PICBgetchromosomes}
\title{Get SeqInfo object from standard non-circular chromosome names from your genome}
\usage{
PICBgetchromosomes(REFERENCE.GENOME, SEQ.LEVELS.STYLE = "UCSC")
}
\arguments{
\item{REFERENCE.GENOME}{name of genome. For example "BSgenome.Dmelanogaster.UCSC.dm6", or directly a SeqInfo object.}

\item{SEQ.LEVELS.STYLE}{naming of chromosomes style. "UCSC" by default.}
}
\value{
SeqInfo object with standard non-circular chromosome names
}
\description{
Get SeqInfo object from standard non-circular chromosome names from your genome
}
\examples{
library(BSgenome.Dmelanogaster.UCSC.dm6)
mySI <- PICBgetchromosomes("BSgenome.Dmelanogaster.UCSC.dm6", "UCSC")
}
\author{
Aleksandr Friman

Franziska Ahrend
}
