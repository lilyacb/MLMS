% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MLMSfunctions.R
\name{iso_ratio_similarity}
\alias{iso_ratio_similarity}
\title{iso_ratio_similarity: function that checks whether isotope ratios are within similarity criteria for specific peaks}
\usage{
iso_ratio_similarity(vend.df, peakNr.vec, sdC.thresh = 0.1, sdO.thresh = 0.1)
}
\arguments{
\item{vend.df}{dataframe of vendor table data}

\item{peakNr.vec}{vector of peak numbers to be analyzed}

\item{sdC.thresh}{threshold for the standard deviation of d13C/12C (default=0.1)}

\item{sdO.thresh}{threshold for the standard deviation of d18O/16O (default=0.1)}
}
\value{
list with boolean values if ratios are within thresholds and the standard deviation values
}
\description{
iso_ratio_similarity: function that checks whether isotope ratios are within similarity criteria for specific peaks
}
\examples{
Usage Example
refIsoRatios<-iso_ratio_similarity(vi.df,c(1,2,4,5,16),sdC.thresh=0.1,sdO.thresh=0.1)
}
