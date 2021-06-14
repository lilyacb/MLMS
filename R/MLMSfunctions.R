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


#' select_file_info: select specific info from collection of .dxf files (Identifier 1, Analysis, Peak Center, file_datetime)
#' @param files vector containing character strings of .dxf filenames
#' @return dataframe of file information - file_idName, Analysis, Peak Center, Date_and_Time
#' @example
#' Usage example
#' @export
select_file_info<-function(files){
  num_files<-length(files)
  msdat<-iso_read_continuous_flow(files[1:num_files])
  file_info<-msdat %>%
    iso_get_file_info(
      select = c(
        #rename
        Name = `Identifier 1`,
      # select columns without renaming
        `Analysis`, `Peak Center`,
      # select the time stamp and rename it to `Date & Time`
        Date_and_Time = file_datetime
      ),
      # explicitly allow for file specific rename (for the new ID column)
      file_specific = TRUE #mostly useful with data from different instruments
    )
  # convert from tibble to df
  file_info.df<-as.data.frame(file_info)
  }

#' get_raw_df: get the raw data from .dxf files as a dataframe
#' @param files vector containing character strings of .dxf filenames
#' @return dataframe containing all raw data in the .dxf files
#' @examples
#' Usage example
#' @export
get_raw_df<-function(files){
  num_files<-length(files)
  msdat<-iso_read_continuous_flow(files[1:num_files])
  raw_dat<-iso_get_raw_data(msdat)
  raw_dat.df<-as.data.frame(raw_dat)
}
