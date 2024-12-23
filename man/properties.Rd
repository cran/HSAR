% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{properties}
\alias{properties}
\title{Dataset of properties in the municipality of Athens}
\format{
An object of class \code{sf} (inherits from \code{data.frame}) with 1000 rows and 7 columns.
}
\usage{
properties
}
\description{
A dataset of apartments in the municipality of Athens for 2017. Point location of the properties is given
together with their main characteristics and the distance to the closest metro/train station.
}
\details{
An sf object of 1000 points with the following 6 variables.
\describe{
\item{id}{An unique identifier for each property.}
\item{size}{The size of the property (unit: square meters)}
\item{price}{The asking price (unit: euros) }
\item{prpsqm}{The asking price per squre meter (unit: euroes/square meter).}
\item{age}{Age of property in 2017 (unit: years).}
\item{dist_metro}{The distance to closest train/metro station (unit: meters).}
}
}
\keyword{datasets}
