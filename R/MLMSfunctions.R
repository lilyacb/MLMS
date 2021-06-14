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


#' select_file_info:
#' @param filename character string of the name of the .dxf file of data
#' @return file information - Name, Analysis, Peak Center, H3 Factor, Date & Time
#' @example
#' Usage example
#' @export
select_file_info<-function(filename){
  msdat<-iso_read_continuous_flow(filename)
  file_info<-msdat %>%
    iso_get_file_info(
      select = c(
      # rename sample id columns from the different file types to a new ID column
        ID = `Identifier 1`, ID = `Name`,
      # select columns without renaming
        `Analysis`, `Peak Center`, `H3 Factor`,
      # select the time stamp and rename it to `Date & Time`
        `Date & Time` = file_datetime
      ),
      # explicitly allow for file specific rename (for the new ID column)
      #file_specific = TRUE #mostly useful with data from different instruments
    )
  }

