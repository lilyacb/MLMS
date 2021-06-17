% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MLMSfunctions.R
\name{plot_ms}
\alias{plot_ms}
\title{plot_ms:}
\usage{
plot_ms(vendor_info.df, x_name = "Rt", y_name = "rIntensity_All")
}
\arguments{
\item{vendor_info.df}{dataframe of vendor info for only one experiment (may need to parse output from select_vendor_info())}

\item{x_name}{name for desired x units for ms plot from vendor_info.df (default Rt)}

\item{y_name}{name for desired y units for ms plot from vendor_info.df (default rIntensity_All)}
}
\value{
PlotSpec plot of the specified columns from vendor_info.df
}
\description{
plot_ms:
}
\examples{
Usgae example
}