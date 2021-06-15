#### Functions for processing .dxf MS data files

#' read_summ: read and print a summary of mass spec data from a .dxf file
#' @param filename character string of the name of the .dxf file of data
#' @return summary table of file contents
#' @examples
#' Usage example
#' read_summ("170525_NaHCO3 L + NaCl L_.dxf")
#' @export
read_summ<-function(filename){
  msdat<-iso_read_continuous_flow(filename)
  summ<-summary(msdat)
  print(summ)
}


#' select_file_info: select specific info from collection of .dxf files (Identifier 1, Analysis, Peak Center, file_datetime)
#' @param files vector containing character strings of .dxf file names
#' @return dataframe of file information - file_idName, Analysis, Peak Center, Date_and_Time
#' @examples
#' Usage example
#' data_files<-c("170525_NaHCO3 L + NaCl L_.dxf","170525_NaHCO3 L + NaCl U_.dxf","170525_NaHCO3 L_.dxf","170525_NaHCO3 U + NaCl L_.dxf",
#' "170525_NaHCO3 U + NaCl U_.dxf","170525_NaHCO3 U_.dxf")
#' select_file_info(data_files)
#' @export
select_file_info<-function(files){
  num_files<-length(files)
  msdat<-iso_read_continuous_flow(files[1:num_files])
  file_info<-msdat %>%
    iso_get_file_info(
      select = c(
        #rename
        Identifier_1 = `Identifier 1`,
      # select columns without renaming
        `Analysis`, `Preparation`,
      # select the time stamp and rename it to `Date & Time`
        Date_and_Time = file_datetime
      ),
      # explicitly allow for file specific rename (for the new ID column)
      file_specific = TRUE #mostly useful with data from different instruments
    )
  # convert from tibble to df
  file_info.df<-as.data.frame(file_info)
  }


#' select_vendor_info: get specific vendor info with labeled experiment names for a collection of .dxf files
#' (Identifier 1, Nr., Start, Rt, End, Area All, d13C/12C, d18O/16O)
#' @param files vector containing character strings of .dxf file names
#' @return dataframe of vendor information with rows labeled with experiment name (Identifier 1)
#' @examples
#' Usage Example
#' data_files<-c("170525_NaHCO3 L + NaCl L_.dxf","170525_NaHCO3 L + NaCl U_.dxf","170525_NaHCO3 L_.dxf","170525_NaHCO3 U + NaCl L_.dxf",
#' "170525_NaHCO3 U + NaCl U_.dxf","170525_NaHCO3 U_.dxf")
#' select_vendor_info(data_files)
#' @export
select_vendor_info<-function(files){
  num_files<-length(files)
  msdat<-iso_read_continuous_flow(files[1:num_files])
  file.info<-msdat %>% iso_get_file_info()
  #file.id<-file.info$`file_id`
  #print(file.id)
  ident1<-file.info$`Identifier 1`
  print(ident1)
  vendor_info<-msdat %>% iso_get_vendor_data_table()
  peak_num<-vendor_info$Nr.
  # create name col for data (not part of vendor data table)
  m=1
  i=1
  name.vec<-c()
  for(m in seq(1:length(data_files))){
    while(peak_num[i+1]>peak_num[i]){
      # add name for this experiment
      name<-ident1[m]
      #print(name)
      name.vec<-c(name.vec,name)
      i=i+1
      #print(i)
      if(i==length(peak_num)){
        break
      }
    }
    i=i+1
    name<-ident1[m]
    name.vec<-c(name.vec,name)
  }
  vendor_info_select<-cbind(name.vec,peak_num,
                 vendor_info$Start,
                 vendor_info$Rt,
                 vendor_info$End,
                 vendor_info$`Area All`,
                 vendor_info$`d 13C/12C`,
                 vendor_info$`d 18O/16O`)
  vendor_info_select.df<-as.data.frame(vendor_info_select)
  colnames(vendor_info_select.df)<-c("Name","Peak_Nr","Start","Rt","End","Area_All","d13C/12C","d18O/16O")
  return(vendor_info_select.df)
}


#' get_raw_df: get the raw data from .dxf files as a dataframe
#' @param files vector containing character strings of .dxf file names
#' @return dataframe containing all raw data in the .dxf files
#' @examples
#' Usage example
#' data_files<-c("170525_NaHCO3 L + NaCl L_.dxf","170525_NaHCO3 L + NaCl U_.dxf","170525_NaHCO3 L_.dxf","170525_NaHCO3 U + NaCl L_.dxf",
#' "170525_NaHCO3 U + NaCl U_.dxf","170525_NaHCO3 U_.dxf")
#' get_raw_df(data_files)
#' @export
get_raw_df<-function(files){
  num_files<-length(files)
  msdat<-iso_read_continuous_flow(files[1:num_files])
  raw_dat<- msdat %>% iso_get_raw_data()
  raw_dat.df<-as.data.frame(raw_dat)
}


#' get_resistor_df: get resistor info for collection of .dxf files as a dataframe
#' @param files vector containing character strings of .dxf file names
#' @return dataframe of resistor info in the .dxf files
#' @examples
#' Usage example
#' data_files<-c("170525_NaHCO3 L + NaCl L_.dxf","170525_NaHCO3 L + NaCl U_.dxf","170525_NaHCO3 L_.dxf","170525_NaHCO3 U + NaCl L_.dxf",
#' "170525_NaHCO3 U + NaCl U_.dxf","170525_NaHCO3 U_.dxf")
#' resist<-get_resistor_df(data_files)
#' @export
get_resistor_df<-function(files){
  num_files<-length(files)
  msdat<-iso_read_continuous_flow(files[1:num_files])
  resistors<-iso_get_resistors(msdat)
  resistors.df<-as.data.frame(resistors)
}


#' get_reference_values_df: get isotopic reference values as a dataframe
