% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClussCluster.R
\name{plot_ClussCluster}
\alias{plot_ClussCluster}
\alias{top.m.hm}
\title{Plots the results of \code{ClussCluster}}
\usage{
plot_ClussCluster(object, m = 10, snames = NULL, gnames = NULL, ...)

top.m.hm(object, m, snames = NULL, gnames = NULL, ...)
}
\arguments{
\item{object}{An object that is obtained by applying the ClussCluster function to the data set.}

\item{m}{The number of top signature genes selected to produce the heatmap.}

\item{snames}{The names of the cells.}

\item{gnames}{The names of the genes}

\item{...}{Addtional parameters, sent to the method}
}
\value{
a ggplot2 object of the heatmap with top signature genes selected by \code{ClussCluster}
}
\description{
Plots the number of signature genes against the tuning parameters if multiple tuning parameters are evaluated in the object. If only one is included, then \code{plot_ClussCluster} returns a venn diagram and a heatmap at this particular tuning parameter.
}
\details{
Takes the normalized and log transformed number of reads mapped to genes (e.g., log(RPKM+1) or log(TPM+1) where RPKM stands for Reads Per Kilobase of transcript per Million mapped reads and TPM stands for transcripts per million) but NOT centered.

If multiple tuning parameters are evaluated in the object, the number of signature genes is computed for each cluster and is plotted against the tuning parameters. Each color and line type corresponds to a cell type.

If only one tuning parameter is evaluated, two plots will be produced. One is the venn diagram of the cell-type-specific genes, the other is the heatmap of the data with the cells and top m signature genes. See more details in the paper.
}
\examples{
data(Hou_sim)
run.cc <- ClussCluster(Hou_sim$x, nclust = 3, ws = c(2.4, 5, 8.8))
plot_ClussCluster(run.cc, m = 5, snames=Hou$snames, gnames=Hou$gnames)

}
