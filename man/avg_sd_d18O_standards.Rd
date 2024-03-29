% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MLMSfunctions.R
\name{avg_sd_d18O_standards}
\alias{avg_sd_d18O_standards}
\title{avg_sd_d18O_standards: function to find d18O/16O averages and standard deviations for each set of files represented in input list}
\usage{
avg_sd_d18O_standards(
  allStandards_d18O.list,
  standNames = c("L1", "H1", "LW"),
  standAcceptedVals.vec = c(-8.55, 4.85, -3.85),
  accStandRatioSD = c(0.2, 0.2, 0.2)
)
}
\arguments{
\item{allStandards_d18O.list}{list built from output of d18O_samples_list(), which retrieves sample peaks and d18O/16O data from specified files}

\item{standNames}{character vector containing the names of the standards}

\item{standAcceptedVals.vec}{numeric vector containing the accepted values of the standards}

\item{accStandRatioSD}{numeric vector containing the accepted standard deviation values for each of the standards (in the same order as names in standNames)}
}
\value{
list of length 3 with the following elements:
[\link{1}]: dataframes of d18O/16O averages and standard deviations for each set of files from input
[\link{2}]: accepted and measured d18O/16O dataframe
[\link{3}]: acceptable and calculated standard deviations for d18O/16O
}
\description{
avg_sd_d18O_standards: function to find d18O/16O averages and standard deviations for each set of files represented in input list
}
\examples{
Usage Example
L1_d18O.list<-d18O_samples_list(L1fileNames.vec)
H1_d18O.list<-d18O_samples_list(H1fileNames.vec)
LW_d18O.list<-d18O_samples_list(LWfileNames.vec)
allStand_d18O.list<-list()
allStand_d18O.list[[1]]<-L1_d18O.list
allStand_d18O.list[[2]]<-H1_d18O.list
allStand_d18O.list[[3]]<-LW_d18O.list
AvgSD_d18O<-avg_sd_d18O_standards(allStand_d18O.list)
}
