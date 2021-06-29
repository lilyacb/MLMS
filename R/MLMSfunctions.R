#### Functions for processing .dxf MS data files

#' cha: function to find the area of a convex hull
#' @param x dataframe of x values
#' @param y dataframe of y values
#' @return area of the polygon
#' @examples
#' Usage example
#' cha(peak1t,peak1v)
#' @export
cha<-function(x,y){
  i<-chull(x,y)
  return(areapl(cbind(x[i],y[i])))
}


#' extract_rintensity_all_tsfeatures: extract time series features from rIntensity_All using tsfeatures
#' @param rintensity_all.num numeric vector containing the rIntensity_All data
#' @return dataframe containing the extracted tsfeatures of the rIntensity_All data
#' @examples
#' Usage Example
#' extract_rintensity_all_tsfeatures(int_all_num)
#' @export
extract_rintensity_all_tsfeatures<-function(rintensity_all.num){
  features.tib<-tsfeatures(rintensity_all.num,
                           features=c("acf_features","arch_stat","crossing_points",
                                      "entropy","flat_spots","heterogeneity",
                                      "holt_parameters","hurst",
                                      "lumpiness","max_kl_shift","max_level_shift",
                                      "max_var_shift","nonlinearity",#"pacf_features",
                                      #"stability",
                                      #"stl_features",
                                      #"unitroot_kpss",
                                      #"unitroot_pp",
                                      #"ac_9",
                                      #"firstmin_ac",
                                      #"firstzero_ac",
                                      #"fluctanal_prop_r1",
                                      "histogram_mode","localsimple_taures","motiftwo_entro3",
                                      "outlierinclude_mdrmd","sampenc","sampen_first",
                                      "std1st_der","trev_num", #,"spreadrandomlocal_meantaul"
                                      "walker_propcross"))
  features.df<-as.data.frame(features.tib)
}


#' file_move: function to move a file from one directory to another, and to create that directory if it does not exist
#' @param from the original file directory
#' @param to the new file directory
#' @examples
#' Usage example
#' file_move(orig_dir,new_dir)
#' @export
file_move<-function(from, to){
  todir <- dirname(to)
  if (!isTRUE(file.info(todir)$isdir)){
    dir.create(todir)
  }
  file.rename(from,to)
}


#' get_all_filenames: function to get all file names from a directory of .dxf files
#' @param path path to the directory of .dxf files
#' @return vector of filenames for all .dxf files in the specified directory
#' @examples
#' Usage example
#' get_all_filenames("~/path_to_files")
#' @export
get_all_filenames<-function(path){ #path to the directory of .dxf files
  all_iso<-iso_read_continuous_flow(path)
  # get just file names
  all_file_info<-iso_get_file_info(all_iso)
  all_file_names<-all_file_info$file_id
}


#' get_all_peak_areas: function to get all peak areas in an experiment
#' @param start.vec numeric vector containing all peak start times
#' @param end.vec numeric vector containing all peak end times
#' @param time.vec numeric vector containing all the raw time.s data
#' @param int.vec numeric vector continaing all the raw intentisty data (i.e. v44.mV, etc,)
#' @return a vector containing the areas of each peak
#' @examples
#' Usage example
#' get_all_peak_areas(start.v1,end.v1,allt.s,allv44)
#' @export
get_all_peak_areas<-function(start.vec,end.vec,time.vec,int.vec){
  all_areas<-c()
  for(i in seq(1:length(start.vec))){
    all_areas<-c(all_areas,peak_area(start.vec[i],end.vec[i],time.vec,int.vec))
  }
  return(all_areas)
}


#' get_identifier_1_files: function to get filenames whose Identifier_1 data matches the one specified
#' @param files vector of file names
#' @param identifier_1 the name of the desired Identifier_1
#' @param cores number of cores to use for the grepl function to search through the files vector
#' @return dataframe of file names with the specified Identifier_1
#' @examples
#' Usage example
#' get_identifier_1_files(my_filenames,my_identifier_1)
#' @export
get_identifier_1_files<-function(files,identifier_1,cores=2){
  identifier_1_files_ind<-pvec(seq_along(files),function(i)
    grepl(identifier_1,files[i],fixed=T),mc.cores=cores)
  # get the indices
  identifier_1_files_ind<-which(identifier_1_files_ind) #which are TRUE
  # get the identifier_1 file names
  identifier_1_files<-files[identifier_1_files_ind]
  i1.df<-as.data.frame(identifier_1_files)
  colnames(i1.df)<-c(identifier_1)
  return(i1.df)
}


#' get_raw_df: get the raw data from .dxf files as a dataframe
#' @param files vector containing character strings of .dxf file names
#' @return dataframe containing all raw data in the .dxf files
#' @examples
#' Usage example
#' get_raw_df(data_files)
#' @export
get_raw_df<-function(files){
  num_files<-length(files)
  msdat<-iso_read_continuous_flow(files[1:num_files])
  raw_dat<- msdat %>% iso_get_raw_data()
  raw_dat.df<-as.data.frame(raw_dat)
}


#' get_reference_values_no_ratio: get isotopic reference values as a dataframe
#' @param files vector containing character strings of .dxf file names
#' @return dataframe of reference values in the .dxf files
#' @examples
#' Usage example
#' get_reference_values_no_ratio(data_files)
#' @export
get_reference_values_no_ratio <- function(files){
  num_files<-length(files)
  msdat<-iso_read_continuous_flow(files[1:num_files])
  # reference delta values without ratio values
  delta_no_ratio<-msdat %>% iso_get_standards(file_id:reference)
  delta_no_ratio.df<-as.data.frame(delta_no_ratio)
}


#' get_reference_values_ratio: get isotopic reference values as a dataframe
#' @param files vector containing character strings of .dxf file names
#' @return dataframe of reference values in the .dxf files
#' @examples
#' Usage example
#' get_reference_values_ratio(data_files)
#' @export
get_reference_values_ratio <- function(files){
  num_files<-length(files)
  msdat<-iso_read_continuous_flow(files[1:num_files])
  # reference values with ratios
  ref_with_ratios<-msdat %>% iso_get_standards()
  ref_with_ratios.df<-as.data.frame(ref_with_ratios)
}


#' get_resistor_df: get resistor info for collection of .dxf files as a dataframe
#' @param files vector containing character strings of .dxf file names
#' @return dataframe of resistor info in the .dxf files
#' @examples
#' Usage example
#' get_resistor_df(data_files)
#' @export
get_resistor_df<-function(files){
  num_files<-length(files)
  msdat<-iso_read_continuous_flow(files[1:num_files])
  resistors<-iso_get_resistors(msdat)
  resistors.df<-as.data.frame(resistors)
}


#' move_Identifier_1_files: function that copies all files with a specified identifier into a folder labeled with the Identifier_1 value
#' @param identifier_1_files files that contain the desired identifier
#' @param path location of the file to be copied
#' @param identifier_1 the identifier_1 of the files to be copied
#' @examples
#' Usage example
#' move_Identifier_1_files(identifier_1_d_files,pathc,identifier_1_d)
#' @export
move_Identifier_1_files<-function(identifier_1_files,path,identifier_1){
  num_files<-dim(identifier_1_files)[1]
  for(i in seq(1:num_files)){
    orig_dir<-paste(path,"/",identifier_1_files[i,],sep="")
    #orig.dir
    new_dir<-paste(path,identifier_1,"/",identifier_1_files[i,],sep="")
    #new.dir
    file_move(orig_dir,new_dir)
  }
}


#' peak_area: Function to calculate the area under a peak
#' @param start.t start time of the peak
#' @param end.t end time of the peak
#' @param time.vec the raw time.s data
#' @param int.vec the raw intensity data for a specified mass
#' @return the area under the specified peak
#' @examples
#' Usage example
#' peak_area(start1,end1,allt.s,allv44)
#' @export
peak_area<-function(start.t,end.t,time.vec,int.vec){
  # get times for peak
  peak.t<-c()
  time.ind<-c()
  # get peak times and indices
  for(i in seq(1:length(time.vec))){
    if((time.vec[i]>=start.t) && (time.vec[i]<=end.t)){
      peak.t<-c(peak.t,time.vec[i])
      time.ind<-c(time.ind,i)
    }
  }
  # get peak intensities
  peak.mv<-c()
  for(i in seq(1:length(time.ind))){
    peak.mv<-c(peak.mv,int.vec[time.ind[i]])
  }
  #plot(peak.t,peak.mv,type="l")
  peak.df<-as.data.frame(cbind(peak.t,peak.mv))
  colnames(peak.df)<-c("time","intensity")
  peak.area<-cha(peak.t,peak.mv)
  return(peak.area)
}


#' plot_individual_peak: function to plot an individual peak in an experiment
#' @param start.t start time of the peak
#' @param end.t end time of the peak
#' @param time.vec vector containing all time.s raw data
#' @param int.vec vectoring containing the raw intensity data for a specified mass (i.e. v44.mV, etc.)
#' @param peak_num the peak number
#' @param v.mv vector containing the raw mV intensity data
#' @examples
#' Usage example
#' plot_individual_peak(start1,end1,allt.s,allv44,"1","v44.mV")
#' @export
plot_individual_peak<-function(start.t,end.t,time.vec,int.vec,peak_num,v.mv){
  peak.t<-c()
  time.ind<-c()
  # get peak times and indices
  for(i in seq(1:length(time.vec))){
    if((time.vec[i]>=start.t) && (time.vec[i]<=end.t)){
      peak.t<-c(peak.t,time.vec[i])
      time.ind<-c(time.ind,i)
    }
  }
  # get peak intensities
  peak.mv<-c()
  for(i in seq(1:length(time.ind))){
    peak.mv<-c(peak.mv,int.vec[time.ind[i]])
  }
  plot(peak.t,peak.mv,main=paste("peak_",peak_num,"_",v.mv,sep=""),type="l")
  peak.df<-as.data.frame(cbind(peak.t,peak.mv))
  colnames(peak.df)<-c("time","intensity")
  return(peak.df)
}


#' plot_ms: Function to plot mass spec data
#' @param vendor_info.df dataframe of vendor info for only one experiment
#' @param x_name name for desired x units for ms plot from vendor_info.df (default Rt)
#' @param y_name name for desired y units for ms plot from vendor_info.df (default rIntensity_All)
#' @return PlotSpec plot of the specified columns from vendor_info.df
#' @examples
#' Usage example
#' plot_ms("Rt","rIntensity_All")
#' @export
plot_ms<-function(vendor_info.df,x_name="Rt",y_name="rIntensity_All"){
  x_ind<-which(colnames(vendor_info.df)==x_name)
  x<-as.numeric(vendor_info.df[,x_ind])
  y_ind<-which(colnames(vendor_info.df)==y_name)
  y<-as.numeric(vendor_info.df[,y_ind])
  plot_dat.df<-as.data.frame(cbind(x,y))
  colnames(plot_dat.df)<-c(x_name,y_name)
  PlotSpec(plot_dat.df)
}


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


#' select_file_info: select specific information from a collection of .dxf files (Identifier 1, Analysis, Preparation, file_datetime)
#' @param files vector containing character strings of .dxf file names
#' @return dataframe of file information - file_id, Identifier_1, Analysis, Preparation, Date_and_Time
#' @examples
#' Usage example
#' select_file_info(data_files)
#' @export
select_file_info<-function(files){
  num_files<-length(files)
  msdat<-iso_read_continuous_flow(files[1:num_files])
  file_info<-msdat %>%
    iso_get_file_info(
      select = c(
        #rename?
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


#' select_vendor_info: get selected vendor info with labeled experiment names for a collection of .dxf files
#' (Identifier 1, Nr., Start, Rt, End, Intensity All, rIntensity All, d13C/12C, d18O/16O)
#' @param files vector containing character strings of .dxf file names
#' @return dataframe of vendor information with rows labeled with experiment name (Identifier 1)
#' @examples
#' Usage Example
#' select_vendor_info(data_files)
#' @export
select_vendor_info<-function(files){
  num_files<-length(files)
  msdat<-iso_read_continuous_flow(files[1:num_files])
  file.info<-msdat %>% iso_get_file_info()
  ident1<-file.info$`Identifier 1`
  #print(ident1)
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
                            vendor_info$`Intensity All`,
                            vendor_info$`rIntensity All`,
                            vendor_info$`d 13C/12C`,
                            vendor_info$`d 18O/16O`)
  vendor_info_select.df<-as.data.frame(vendor_info_select)
  colnames(vendor_info_select.df)<-c("Identifier_1","Peak_Nr","Start","Rt","End","Intensity_All","rIntensity_All","d13C/12C","d18O/16O")#"Area_All"
  return(vendor_info_select.df)
}
