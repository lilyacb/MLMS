#### Functions for processing .dxf MS data files

#' read_summ: read and print a summary of mass spec data from a .dxf file
#' @param filename character string of the name of the .dxf file of data
#' @return summary table of file contents
#' @example
#' Usage example
#' @export
read_summ<-function(filename){
  msdat<-iso_read_continuous_flow(filename)
  summ<-summary(msdat)
  print(summ)
}
