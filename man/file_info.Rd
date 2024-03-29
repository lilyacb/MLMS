% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MLMSfunctions.R
\name{file_info}
\alias{file_info}
\title{file_info: select specific information from a collection of .dxf files (Identifier 1, Analysis, Preparation, file_datetime)}
\usage{
file_info(files)
}
\arguments{
\item{files}{vector containing character strings of .dxf file names}
}
\value{
dataframe of file information - file_id, Identifier_1, Analysis, Preparation, Date_and_Time
}
\description{
file_info: select specific information from a collection of .dxf files (Identifier 1, Analysis, Preparation, file_datetime)
}
\examples{
Usage example
file_info(data_files)
}
