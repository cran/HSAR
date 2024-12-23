#' @title Boundaries of districts in Beijing
#' @name Beijingdistricts
"Beijingdistricts"

#' @title The spatial locations of the Beijing land price data
#' @name land
"land"

#' @title Leased residential land parcels, from 2003 to 2009 in Beijing, China
#' @name landprice
#' @details
#' A `data.frame` with 1117 observations on the following 11 variables.
#'  \describe{
#'    \item{obs}{An unique identifier for each land parcel.}
#'    \item{lnprice}{The log of the leasing price per square metre of each residential land parcel (unit: RMB, Chinese yuan) }
#'    \item{dsubway}{The log of the distance of each land parcel to the nearest railway station (unit:meters)}
#'    \item{dele}{The log of the distance of each land parcel to the nearest elementary school (unit:meters) }
#'    \item{dpark}{The log of the distance of each land parcel to the nearest green park (unit:meters) }
#'    \item{lnarea}{The log of the size of each land parcel (unit: square meters).}
#'    \item{lndcbd}{The log of the distance of each land parcel to the CBD (centre business district) in Beijing (unit:meters) }
#'    \item{year}{The year when each land parcel was leased with values of 0,1,2,3,4,5,6 representing year 2003,2004,2005,2006,2007,2008,2009}
#'    \item{popden}{The population density of each district (unit: 1000 persons per square kilometers) }
#'    \item{crimerate}{The number of reported serious crimes committed in each district per 1000 persons.}
#'    \item{district.id}{The identifier of the district where each land parcel is located.}
#' }
"landprice"

#' @title Municipality departments of Athens
#' @name depmunic
#' @details
#' An sf object of 7 polygons with the following 7 variables:
#'  \describe{
#'     \item{num_dep}{An unique identifier for each municipality department.}
#'     \item{airbnb}{The number of airbnb properties in 2017}
#'     \item{museums}{The number of museums}
#'     \item{population}{The population recorded in census at 2011.}
#'     \item{pop_rest}{The number of citizens that the origin is a non european country.}
#'     \item{greensp}{The area of green spaces (unit: square meters).}
#'     \item{area}{The area of the polygon (unit: square kilometers). }
#'  }
"depmunic"

#' @title Dataset of properties in the municipality of Athens
#' @description
#' A dataset of apartments in the municipality of Athens for 2017. Point location of the properties is given
#' together with their main characteristics and the distance to the closest metro/train station.
#' @details
#' An sf object of 1000 points with the following 6 variables.
#' \describe{
#'   \item{id}{An unique identifier for each property.}
#'   \item{size}{The size of the property (unit: square meters)}
#'   \item{price}{The asking price (unit: euros) }
#'   \item{prpsqm}{The asking price per squre meter (unit: euroes/square meter).}
#'   \item{age}{Age of property in 2017 (unit: years).}
#'   \item{dist_metro}{The distance to closest train/metro station (unit: meters).}
#' }
#' @name properties
"properties"
