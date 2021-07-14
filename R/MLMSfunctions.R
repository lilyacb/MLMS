#### Functions for processing .dxf MS data files
# Function in this file:
# (1) all_filenames(path)
# (2) all_peak_areas_poly(start.vec,end.vec,time.vec,int.vec)
# (3) all_peak_areas_trap()
# (4) cha(x,y)
# (5) extract_tsfeatures(intensity.data)
# (6) file_info(files)
# (7) file_move(from, to)
# (8) groupPickedTimes(timesPicked)
# (9) identifier_1_files(files,identifier_1,cores=2)
# (10) move_identifier_1_files(identifier_1_files,path,identifier_1)
# (11) peak_area_poly(start.t,end.t,time.vec,int.vec)
# (12) peak_area_trap(start.t,end.t,time.vec,int.vec)
# (13) peakTimes(Int.mat,time.interval)
# (14) plot_all_peaks(start.v,end.v,t.s,v.mV,v.name)
# (15) plot_individual_peaks(start.t,end.t,time.vec,int.vec,peak_num,v.mv)
# (16) plot_ms(vendor_info.df,x_name="Rt",y_name="Intensity_All")
# (17) qc_num_vendPeaks(vend.df)
# (18) qc_peaks_present(start.times)
# (19) qc_samplePeaks(intensity.mat,time.interval,expectedNum.samplePeaks,expectedStartSamples,samplePeakTimes)
# (20) unique_identifiers(path,filenames)
# (21) raw_data(file)
# (22) read_summary(filename)
# (23) reference_values_no_ratio(files)
# (24) reference_values_ratio(files)
# (25) remove_276(vend.df)
# (26) resistor_data(files)
# (27) sample_peaks_vend(vend.df)
# (28) sort_by_identifier_1(path)
# (29) ThresholdingAlgo(y,lag,threshold,influence)
# (30) vendor_info(files)


# (1)
#' all_filenames: function to get all file names from a directory of .dxf files
#' @param path path to the directory of .dxf files
#' @return vector of filenames for all .dxf files in the specified directory
#' @examples
#' Usage example
#' all_filenames("~/path_to_files")
#' @export
all_filenames<-function(path){ #path to the directory of .dxf files
  all_iso<-iso_read_continuous_flow(path)
  # get just file names
  all_file_info<-iso_get_file_info(all_iso)
  all_file_names<-all_file_info$file_id
}


# (2)
#' all_peak_areas_poly: function to get all peak areas in an experiment
#' @param start.vec numeric vector containing all peak start times
#' @param end.vec numeric vector containing all peak end times
#' @param time.vec numeric vector containing all the raw time.s data
#' @param int.vec numeric vector continaing all the raw intentisty data (i.e. v44.mV, etc,)
#' @return a vector containing the areas of each peak in V*s
#' @examples
#' Usage example
#' all_peak_areas_poly(start.v1,end.v1,allt.s,allv44)
#' @export
all_peak_areas_poly<-function(start.vec,end.vec,time.vec,int.vec){
  all_areas<-c()
  for(i in seq(1:length(start.vec))){
    all_areas<-c(all_areas,peak_area_poly(start.vec[i],end.vec[i],time.vec,int.vec))
  }
  all_areas_Vs<-all_areas/1000
  return(all_areas_Vs)
}


# (3)
#' all_peak_areas_trap: function to calculate all peak areas using the integration trapezoidal rule via trapz
#' @param start.vec numeric vector of peak start times
#' @param end.vec numeric vector of peak end times
#' @param time.vec numeric vector of times from the raw data
#' @param int.vec numeric vector of intensities from the raw data
#' @return a numeric vector containing all the peak areas in V*s
#' @examples
#' Usage Example
#' all_peak_areas_trap<-all_peak_areas_trap(start.v1,end.v1,time.s,v44)
#' @export
all_peak_areas_trap<-function(start.vec,end.vec,time.vec,int.vec){
  all_areas<-c()
  for(i in seq(1:length(start.vec))){
    all_areas<-c(all_areas,peak_area_trap(start.vec[i],end.vec[i],time.vec,int.vec))
  }
  all_areas_Vs<-all_areas
  return(all_areas_Vs)
}


# (4)
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


# (5)
#' extract_tsfeatures: extract time series features from intensity data using tsfeatures
#' @param intensity.data numeric vector containing the rIntensity_All data
#' @return dataframe containing the extracted tsfeatures of the rIntensity_All data
#' @examples
#' Usage Example
#' extract_tsfeatures(int_all_num)
#' @export
extract_tsfeatures<-function(intensity.data){
  features.tib<-tsfeatures(intensity.data,
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


# (6)
#' file_info: select specific information from a collection of .dxf files (Identifier 1, Analysis, Preparation, file_datetime)
#' @param files vector containing character strings of .dxf file names
#' @return dataframe of file information - file_id, Identifier_1, Analysis, Preparation, Date_and_Time
#' @examples
#' Usage example
#' file_info(data_files)
#' @export
file_info<-function(files){
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


# (7)
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


# (8)
#' groupPickedTimes: function that organizes the times picked from the PeakTimes function
#' into single average times for peaks
#' @param timesPicked the times returned in a list from ThresholdingAlgo
#' @return the average times from each "group" (groups returned from ThresholdingAlgo$signals that needs to be separated)
#' @examples
#' Usage Example
#' groupPickedTimes(timesPicked.list)
#' @export
groupPickedTimes<-function(timesPicked){
  group.list<-list()
  listIndex<-1
  timeGroup<-c()
  for(i in seq(1,length(timesPicked))){
    if(i==length(timesPicked)){
      #put the last element in the last group
      inGroup<-abs(timesPicked[i]-timesPicked[(i-1)])<3
      timeGroup<-c(timeGroup,timesPicked[i])
      #print(paste("end of Group",listIndex,sep=""))
      group.list[[listIndex]]<-timeGroup
      break
    } else if(abs(timesPicked[(i+1)]-timesPicked[i])<3){
      timeGroup<-c(timeGroup,timesPicked[i])
      #timeGroup
    } else{
      #print(paste("end of Group",listIndex,sep=""))
      group.list[[listIndex]]<-timeGroup
      #print(group.list)
      listIndex<-listIndex+1
      timeGroup<-c()
    }
  }
  avgTime.vec<-c()
  for(m in seq(1,length(group.list))){
    avgTime<-mean(group.list[[m]])
    avgTime.vec<-c(avgTime.vec,avgTime)
  }
  return(avgTime.vec)
}


# (9)
#' identifier_1_files: function to get filenames whose Identifier_1 data matches the one specified
#' @param files vector of file names
#' @param identifier_1 the name of the desired Identifier_1
#' @param cores number of cores to use for the grepl function to search through the files vector (default=2)
#' @return dataframe of file names with the specified Identifier_1
#' @examples
#' Usage example
#' identifier_1_files(my_filenames,my_identifier_1)
#' @export
identifier_1_files<-function(files,identifier_1,cores=2){
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


# (10)
#' move_identifier_1_files: function that moves all files with a specified identifier into a folder labeled with the Identifier_1 value
#' @param identifier_1_files files that contain the desired identifier
#' @param path location of the file to be copied
#' @param identifier_1 the identifier_1 of the files to be copied
#' @examples
#' Usage example
#' move_Identifier_1_files(identifier_1_d_files,pathc,identifier_1_d)
#' @export
move_identifier_1_files<-function(identifier_1_files,path,identifier_1){
  num_files<-dim(identifier_1_files)[1]
  for(i in seq(1:num_files)){
    orig_dir<-paste(path,"/",identifier_1_files[i,],sep="")
    #orig.dir
    new_dir<-paste(path,identifier_1,"/",identifier_1_files[i,],sep="")
    #new.dir
    file_move(orig_dir,new_dir)
  }
}


# (11)
#' peak_area_poly: Function to calculate the area under a peak
#' @param start.t start time of the peak
#' @param end.t end time of the peak
#' @param time.vec the raw time.s data
#' @param int.vec the raw intensity data for a specified mass
#' @return the area under the specified peak
#' @examples
#' Usage example
#' peak_area_poly(start1,end1,allt.s,allv44)
#' @export
peak_area_poly<-function(start.t,end.t,time.vec,int.vec){
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


# (12)
#' peak_area_trap: function that returns the area of a peak using the integration trapezoidal rule via trapz
#' @param start.t start time of the peak
#' @param end.t end time of the peak
#' @param time.vec time vector of times in the interval for the peak
#' @param int.vec intensity vector of intensity values for the peak
#' @return peak area in V*s
#' @examples
#' Usage example
#' peak_area_trap(start.v1[1],end.v1[1],time.s,v45)
#' @export
peak_area_trap<-function(start.t,end.t,time.vec,int.vec){
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
  #peak.area<-cha(peak.t,peak.mv)
  peak.area<-trapz(peak.t,peak.mv)/1000
  return(peak.area)
}


# (13)
#' peakTimes: function
#' @param Int.mat matrix of intensity values within the desired interval
#' @param z.thresh z-value threshold for ThresholdAlgo
#' @param time.interval time values in the desired interval
#' @return vector of average peak start times as found by ThresholdingAlgo and groupPickedTimes
#' @examples
#' Usage Example
#' PeakTimes(v44Int.mat,interval.times)
#' @export
# Int.mat: intensity values within the desired time interval
# time.interval: time values in the desired interval
peakTimes<-function(Int.mat,time.interval,z.thresh,reactionLabel){
  testAlgo<-ThresholdingAlgo(Int.mat,lag=10,threshold=z.thresh,influence=0.5)
  #(testAlgo$signals)
  plot(Int.mat,type="l")
  title(main=paste("Peak times detected:",reactionLabel,sep=""))
  x<-which(testAlgo$signals==1) # what about -1?
  y<-rep(0,length(x))
  points(x,y,col="red",pch=20)
  timesPicked<-time.interval[x]
  timeGroups<-groupPickedTimes(timesPicked)
  return(timeGroups)
}


# (14)
#' plot_all_peaks: Function that plots every peak in an experiment in a separate graph
#' @param start.v numeric vector of start times for all the peaks in the experiment
#' @param end.v numeric vector of end times for all the peaks in the experiment
#' @param t.s numeric vector time.s from raw data
#' @param v.mV numeric vector of raw voltages (i.e. v44.mV, v45.mV, v46.mV, etc.)
#' @param v.name the name of raw voltage column
#' @examples
#' Usage example
#' plot_all_peaks(start.v1,end.v1,allt.s,allv44,"v44.mV")
#' @export
plot_all_peaks<-function(start.v,end.v,t.s,v.mV,v.name){
  for(i in seq(1:length(start.v))){
    plot_individual_peak(start.v[i],end.v[i],t.s,v.mV,i,v.name)
  }
}


# (15)
#' plot_individual_peaks: function to plot an individual peak in an experiment
#' @param start.t start time of the peak
#' @param end.t end time of the peak
#' @param time.vec vector containing all time.s raw data
#' @param int.vec vectoring containing the raw intensity data for a specified mass (i.e. v44.mV, etc.)
#' @param peak_num the peak number
#' @param v.mv vector containing the raw mV intensity data
#' @examples
#' Usage example
#' plot_individual_peaks(start1,end1,allt.s,allv44,"1","v44.mV")
#' @export
plot_individual_peaks<-function(start.t,end.t,time.vec,int.vec,peak_num,v.mv){
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


# (16)
#' plot_ms: Function to plot mass spec data from vendor table
#' @param vendor_info.df dataframe of vendor info for only one experiment
#' @param x_name name for desired x units for ms plot from vendor_info.df (default Rt)
#' @param y_name name for desired y units for ms plot from vendor_info.df (default rIntensity_All)
#' @return PlotSpec plot of the specified columns from vendor_info.df
#' @examples
#' Usage example
#' plot_ms("Rt","Intensity_All")
#' @export
plot_ms<-function(vendor_info.df,x_name="Rt",y_name="Intensity_All"){
  x_ind<-which(colnames(vendor_info.df)==x_name)
  x<-as.numeric(vendor_info.df[,x_ind])
  y_ind<-which(colnames(vendor_info.df)==y_name)
  y<-as.numeric(vendor_info.df[,y_ind])
  plot_dat.df<-as.data.frame(cbind(x,y))
  colnames(plot_dat.df)<-c(x_name,y_name)
  PlotSpec(plot_dat.df)
}


# (17)
#' qc_num_vendPeaks: Function that checks for the correct number of peaks
#' @param vend.df dataframe of vendor data that contains Peak_Nr
#' @param expectedPeakNum1 first value of expected peak numbers (i.e. if a peak can be missing - set both to the same number if no variation in peak numbers is expected)
#' @param expectedPeakNum2 second value of expected peak numbers
#' @return Boolean that indicates whether the expected number of peaks is present
#' @examples
#' Usage example
#' qc_num_vendPeaks(vend1.df,15,16)
#' @export
qc_num_vendPeaks<-function(vend.df,expectedPeakNum1=15,expectedPeakNum2=16){
  num_peaks<-as.numeric(vend.df$Peak_Nr)
  print(paste("Number of peaks in vendor table: ",length(num_peaks),sep=""))
  if(length(num_peaks==expectedPeakNum1) | length(num_peaks==expectedPeakNum2)){
    peaks_qc<-TRUE
    #print(paste("number of peaks:",length(num_peaks)))
  }
  else{
    peaks_qc<-FALSE
  }
  return(peaks_qc)
}


# (18)
#' qc_peaks_present: Function to check for peak presence
#' @param start.times vector of start times for each peak in the data
#' @return Boolean that indicates whether peaks are present
#' @examples
#' Usage example
#' qc_peaks_present(start_times.vec)
#' @export
qc_peaks_present<-function(start.times){
  start.times<-as.numeric(start.times)
  if(length(start.times)>=0){ # What do start and end look like if there are no peaks? blank or NA?
    peaks=TRUE
  }
  else{
    peaks=FALSE
  }
  return(peaks)
}


# (19)
#' qc_samplePeaks: function that tests whether the expected number of sample peaks are present and whether they occur at the expected times
#' @param intensity.mat matrix of intensity values for the desired interval
#' @param time.interval numeric vector of time values for the desired interval
#' @param expectedNum.samplePeaks number of expected sample peaks
#' @param expectedStartSamples numeric vector of expected start times
#' @param samplePeakTimes numeric vector of times that sample peaks occur, such as that returned by PeakTimes()
#' @return Boolean dataframe with columns "Expected_Number_Sample_Peaks" and "Expected_Start_Times" which indicate whether the expected number of sample peaks and expected start times are found
#' @examples
#' Usage Example
#' qc_samplePeaks(intensity.mat=v44Int.mat,time.interval=above_below.t,expectedNum.samplePeaks=9,
#'                 expectedStartSamples=start.v1[7:15],samplePeakTimes=peakTimes)
#' @export
qc_samplePeaks<-function(intensity.mat,time.interval,expectedNum.samplePeaks,expectedStartSamples,samplePeakTimes){
  peakTimes<-peakTimes(intensity.mat,time.interval)
  samplePeaks<-c()
  if(length(peakTimes)==expectedNum.samplePeaks){
    peaks_qc<-TRUE
    #print(peaks_qc)
    samplePeaks<-c(peaks_qc)
  } else{
    peaks_qc<-FALSE
    #print(peaks_qc)
    samplePeaks<-c(peaks_qc)
  }
  print(paste("Number of sample peaks detected:",length(peakTimes)))
  # check that values are close to reported start times[7:15]?
  absT<-abs(expectedStartSamples-samplePeakTimes)
  #absT
  withinExpected<-(absT/expectedStartSamples)<0.01
  #withinExpected
  if(sum(withinExpected)==expectedNum.samplePeaks){
    allWithinExpected<-TRUE
    print("All sample peak start times within expected time intervals.")
  } else{
    allWithinExpected<-FALSE
    print("Not all sample peak start times within expected time intervals")
    notIn<-which(withinExpected==FALSE)
    print(paste("Sample peak number not within expected start times: ",notIn,sep=""))
  }
  if(allWithinExpected==TRUE){
    samplePeaks<-c(samplePeaks,allWithinExpected)
  } else{
    samplePeaks<-c(samplePeaks,allWithinExpected)
  }
  samplePeaks.mat<-matrix(samplePeaks,ncol=2)
  samplePeaks.mat
  samplePeaks.df<-as.data.frame(samplePeaks.mat)
  samplePeaks.df
  colnames(samplePeaks.df)<-c("Expected_Number_Sample_Peaks","Expected_Start_Times")
  samplePeaks.df
  return(samplePeaks.df)
}


# (20)
#' unique_identifiers: Function to find all of the different Identifier_1 labels in a directory of .dxf files
#' @param filenames vector of filenames from which to search for unique Identifier_1 labels
#' @return vector of each different Identifier_1 in the given files
#' Usage example
#' unique_identifiers(double_component_files)
#' @export
unique_identifiers<-function(filenames){
  #filenames<-get_all_filenames(path)
  vendor.info<-select_file_info(filenames)
  all_identifiers<-vendor.info$Identifier_1
  num_all_identifiers<-length(all_identifiers)
  #num_all_identifiers
  unique_identifiers<-c()
  for(i in seq(1:num_all_identifiers)){
    curr_identifier<-all_identifiers[i]
    curr_identifier
    if(sum(which(unique_identifiers==curr_identifier))==0){
      unique_identifiers<-c(unique_identifiers,curr_identifier)
    }
  }
  return(unique_identifiers)
}


# (21)
#'raw_data: get the raw data from the specified .dxf file as a dataframe
#' @param file character string defining the file to extract raw data from
#' @return dataframe containing the raw data in the .dxf file
#' @examples
#' Usage example
#' raw_data(data_files)
#' @export
raw_data<-function(file){
  #num_files<-length(files)
  msdat<-iso_read_continuous_flow(file)#files[1:num_files]
  raw_dat<- msdat %>% iso_get_raw_data()
  raw_dat.df<-as.data.frame(raw_dat)
}


# (22)
#' read_summary: read and print a summary of mass spec data from a .dxf file
#' @param filename character string of the name of the .dxf file of data
#' @return summary table of file contents
#' @examples
#' Usage example
#' read_summ("170525_NaHCO3 L + NaCl L_.dxf")
#' @export
read_summary<-function(filename){
  msdat<-iso_read_continuous_flow(filename)
  summ<-summary(msdat)
  print(summ)
}


# (23)
#' reference_values_no_ratio: get isotopic reference values as a dataframe
#' @param files vector containing character strings of .dxf file names
#' @return dataframe of reference values in the .dxf files
#' @examples
#' Usage example
#' reference_values_no_ratio(data_files)
#' @export
reference_values_no_ratio <- function(files){
  num_files<-length(files)
  msdat<-iso_read_continuous_flow(files[1:num_files])
  # reference delta values without ratio values
  delta_no_ratio<-msdat %>% iso_get_standards(file_id:reference)
  delta_no_ratio.df<-as.data.frame(delta_no_ratio)
}


# (24)
#' reference_values_ratio: get isotopic reference values as a dataframe
#' @param files vector containing character strings of .dxf file names
#' @return dataframe of reference values in the .dxf files
#' @examples
#' Usage example
#' reference_values_ratio(data_files)
#' @export
reference_values_ratio <- function(files){
  num_files<-length(files)
  msdat<-iso_read_continuous_flow(files[1:num_files])
  # reference values with ratios
  ref_with_ratios<-msdat %>% iso_get_standards()
  ref_with_ratios.df<-as.data.frame(ref_with_ratios)
}


# (25)
#' remove_276: function to remove the (1st sample) peak at 276 which often has some contamination for CO2
#' @param vend.df dataframe of full vendor table data
#' @return new vend.df without the 276 peak data
#' @examples
#' Usage Example
#' remove_276(vend1.df)
#' @export
remove_276<-function(vend.df){
  num_peaks<-as.numeric(length(vend.df$Peak_Nr))
  start.times<-as.numeric(vend.df$Start)
  Rts<-as.numeric(vend.df$Rt)
  peak_nums<-as.numeric(vend.df$Nr)
  end.times<-as.numeric(vend.df$End)
  for(i in seq(1:length(Rts))){
    if(abs(Rts[i]-276)<5){
      peak_276_Nr<-i
      print("276 peak detected")
      newVend.df<-vend.df[-(peak_276_Nr),]
      start.276<-start.times[i]
      end.276<-end.times[i]
      print(paste("Peak number", peak_276_Nr,"removed from",round(start.276,2),
                  "s to",round(end.276,2),"s with Rt:",round(Rts[i],2)))
    }
  }
  return(newVend.df)
}


# (26)
#' resistor_data: get resistor info for collection of .dxf files as a dataframe
#' @param files vector containing character strings of .dxf file names
#' @return dataframe of resistor info in the .dxf files
#' @examples
#' Usage example
#' resistor_df(data_files)
#' @export
resistor_data<-function(files){
  num_files<-length(files)
  msdat<-iso_read_continuous_flow(files[1:num_files])
  resistors<-iso_get_resistors(msdat)
  resistors.df<-as.data.frame(resistors)
}


# (27)
#' sample_peaks_vend: function to grab the vendor info for the peaks at expected sample peak times
#' @param vend.df vendor info dataframe of all peaks
#' @param start.sample approximate expected start time for sample peaks  (default 326s)
#' @param stop.sample approximate expected stop time for sample peaks (default 825s)
#' @return vendor table for just the (sample) peaks from the given start and stop times
#' @examples
#' Usage Example
#' sample_peaks_vend(no276Vend1.df,326,825)
#' @export
sample_peaks_vend<-function(vend.df,start.sample=326,stop.sample=825){
  #vend.df<-newVend1.df
  num_peaks<-as.numeric(length(vend.df$Peak_Nr))
  #num_peaks
  start.times<-as.numeric(vend.df$Start)
  #start.times
  #peak_nums<-as.numeric(vend.df$Nr)
  sample.indices<-c()
  for(i in seq(1:num_peaks)){
    if((start.times[i]>start.sample) && (start.times[i]<stop.sample)){
      sample.indices<-c(sample.indices,i)
    }
  }
  #sample.indices
  sampleVend.df<-vend.df[sample.indices,]
  return(sampleVend.df)
}


# (28)
#' sort_by_identifier_1: Function that sorts all .dxf files in a path into folders labeled with the Identifier_1 values
#' @param path path to files to be sorted
#' @examples
#' Usage example
#' sort_by_identifier_1("~lily/Desktop/Europa MLMS/Europa_Data/Europa double comp")
#' @export
sort_by_identifier_1<-function(path){
  all_filenames<-get_all_filenames(path)
  un_identifiers<-get_unique_identifiers(path,all_filenames)
  num_identifiers<-length(un_identifiers)
  path_ID<-paste(path,"/",sep="")
  for(i in seq(1:num_identifiers)){
    # get all files with un_identifier[i]
    I1_files<-get_identifier_1_files(all_filenames,un_identifiers[i])
    move_identifier_1_files(I1_files,path_ID,un_identifiers[i])
  }
}


# (29)
#' ThresholdingAlgo: function based on dispersion
#' Code from: Robust peak detection algorithm (using z-scores)
#' https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data?page=1&tab=votes#tab-top
#' @param y intensity matrix (1 col)
#' @param lag lag of the moving window
#' @param threshold z-score at which the algorithm signals
#' @param influence the influence (between 0 and 1) of new signals on the mean and standard deviation
#' @return list of "signals", "avgFilter", "stdFilter" - signal times (signals) are of the most interest for MS purposes
#' @examples
#' Usage Example
#' ThresholdingAlgo(v44Int.mat,lag=10,threshold=5,influence=0.5)
#' @export
ThresholdingAlgo <- function(y,lag,threshold,influence) {
  signals <- rep(0,length(y))
  filteredY <- y[0:lag]
  avgFilter <- NULL
  stdFilter <- NULL
  avgFilter[lag] <- mean(y[0:lag], na.rm=TRUE)
  stdFilter[lag] <- sd(y[0:lag], na.rm=TRUE)
  for (i in (lag+1):length(y)){
    if (abs(y[i]-avgFilter[i-1]) > threshold*stdFilter[i-1]) {
      if (y[i] > avgFilter[i-1]) {
        signals[i] <- 1;
      } else {
        signals[i] <- -1;
      }
      filteredY[i] <- influence*y[i]+(1-influence)*filteredY[i-1]
    } else {
      signals[i] <- 0
      filteredY[i] <- y[i]
    }
    avgFilter[i] <- mean(filteredY[(i-lag):i], na.rm=TRUE)
    stdFilter[i] <- sd(filteredY[(i-lag):i], na.rm=TRUE)
  }
  return(list("signals"=signals,"avgFilter"=avgFilter,"stdFilter"=stdFilter))
}


# (30)
#' vendor_info: get selected vendor info with labeled experiment names for a collection of .dxf files
#' (Identifier 1, Nr., Start, Rt, End, Intensity All, rIntensity All, d13C/12C, d18O/16O)
#' @param files vector containing character strings of .dxf file names
#' @return dataframe of vendor information with rows labeled with experiment name (Identifier 1)
#' @examples
#' Usage Example
#' select_vendor_info(data_files)
#' @export
vendor_info<-function(files){
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
                            vendor_info$`Ampl 44`,
                            vendor_info$`Ampl 45`,
                            vendor_info$`Ampl 46`,
                            #vendor_info$`Intensity 44`,
                            #vendor_info$`rIntensity 44`,
                            #vendor_info$`Intensity 45`,
                            #vendor_info$`rIntensity 45`,
                            #vendor_info$`Intensity 46`,
                            #vendor_info$`rIntensity 46`,
                            vendor_info$`d 13C/12C`,
                            vendor_info$`d 18O/16O`)
  vendor_info_select.df<-as.data.frame(vendor_info_select)
  colnames(vendor_info_select.df)<-c("Identifier_1","Peak_Nr","Start","Rt","End","Intensity_All",
                                     "rIntensity_All",#"Intensity_44","rIntensity_44",
                                     #"Intensity_45","rIntensity_45","Intensity_46","rIntensity_46",
                                     "Ampl_44","Ampl_45","Ampl_46",
                                     "d13C/12C","d18O/16O")#"Area_All"
  return(vendor_info_select.df)
}

