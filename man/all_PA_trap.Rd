% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MLMSfunctions.R
\name{all_PA_trap}
\alias{all_PA_trap}
\title{all_PA_trap: function to calculate all peak areas using the integration trapezoidal rule (via trapz)}
\usage{
all_PA_trap(start.vec, end.vec, time.vec, int.vec, pk.Nrs)
}
\arguments{
\item{start.vec}{numeric vector of peak start times}

\item{end.vec}{numeric vector of peak end times}

\item{time.vec}{numeric vector of times from the raw data}

\item{int.vec}{numeric vector of intensities from the raw data}
}
\value{
a numeric vector containing all the peak areas in V*s
}
\description{
all_PA_trap: function to calculate all peak areas using the integration trapezoidal rule (via trapz)
}
\examples{
Usage Example
all_PA_trap(start.v1,end.v1,time.s,v44,pk.Nrs)
}
