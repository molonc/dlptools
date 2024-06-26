% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_tree_heatmap.R
\name{get_clone_id_label_positions}
\alias{get_clone_id_label_positions}
\title{get a row number of where to place each clone label.}
\usage{
get_clone_id_label_positions(clones, only_largest_clone_group = TRUE)
}
\arguments{
\item{clones}{table of cell_id, clone_id}

\item{only_largest_clone_group}{boolean. TRUE means only largest clone clade
will receive the clone ID label. FALSE means all clades wil get the clone
label.}
}
\value{
tibble of clone_id, and a row_number to place it, which
ComplexHeatmap understands as the position to place the label on the left
side annotation.
}
\description{
A single clone can appear on multiple clades of the tree or in one
monophyletic group. On the heatmap, we need to determine were to put a clone
ID label (e.g., A, B, C, etc), as we don't want a label for each cell id.
Typically, we label the largest clade with the clone ID, and then provide a
legend mapping clone colors to IDs. This function will find the "middle" cell
ID of the clade being labelled as the position to place the clone ID.
}
