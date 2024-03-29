% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CELESTA_functions.R
\name{GetNeighborInfo}
\alias{GetNeighborInfo}
\title{GetNeighborInfo}
\usage{
GetNeighborInfo(coords, number_of_neighbors = 5, bandwidth = 100)
}
\arguments{
\item{coords}{the x, y coordinates of each cell}

\item{number_of_neighbors}{the number of cells in a single neighborhood}

\item{bandwidth}{the upper distance bound used when calculating neighborhoods
by distance}
}
\value{
a list of the following information
\describe{
\item{\code{nb_list}}{the list of N-nearest neighbors}
\item{\code{all_cell_nb_in_bandwidth}}{the cells located within a bandwidth to
cell \emph{c}}
\item{\code{cell_nb_dist}}{the distance of each cell to cell \emph{c} within a
bandwidth}
}
}
\description{
Gets the neighborhood information, including neighborhoods by
number and distance.
}
