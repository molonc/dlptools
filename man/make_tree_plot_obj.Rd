% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_tree_heatmap.R
\name{make_tree_plot_obj}
\alias{make_tree_plot_obj}
\title{create a ggplot object of a phylogenetic tree}
\usage{
make_tree_plot_obj(phylo_tree, clones = NULL)
}
\arguments{
\item{phylo_tree}{a phylo object}

\item{clones}{(optional) adds clones as OTUs (TODO is to allow this to then
make colors)}
}
\value{
ggplot object
}
\description{
Will make a plot of a tree that can be added to a heatmap object. The plot
itself is just a simple tree with no tip labels.
}
