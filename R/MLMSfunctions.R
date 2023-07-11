################################## Functions for MLMS ############################
# Functions in this file
#
# (1) all_filenames(path) -
#         get all .dxf filenames in directory
#
# (2) all_PA_trap(start.vec,end.vec,time.vec,int.vec)
#         calculate all peak areas using trapz
#
# (3) avg_sd_d18O_standards(allStandards_d18O.list,standNames=c("L1","H1","LW"),standAcceptedVals.vec=c(-8.55,4.85,-3.85),accStandRatioSD=c(0.2,0.2,0.2)) -
#         find average and sd values for d18O/16O standards data in a list built from output of d18O_samples_list()
#
# (4) average_times(peakTimes.list,nPeaks) -
#         find the averages of each peak time for vendor data of different experiments
#
# (5) build_expected_times(ref.ind,sample.ind,avgs.df) -
#         build expected times data for reference and sample peaks
#
# (6) check_picked_and_vendor_peaks(vend.df,int.mat,time.vec,z.thresh=5,graphLab,timeThresh=2) -
#         check peak times and numbers from vendor table against picked peaks
#
# (7) d18O_samples_list(standardFileNames) -
#         get list of sample peaks and d18O data for the specified files
#
# (8) extract_tsfeatures(intensity.data) -
#         extract time series features from intensity data using tsfeatures
#
# (9) file_info(files) -
#         select info from .dxf files (Identifier 1, Analysis, Preparation, file_datetime)
#
# (10) file_move(from, to) -
#         move a file from one directory to another (creates directory if it does not exist)
#
# (11) filterSpecSignals(threshAlgoList) -
#         filters out possible false positives in peak detection
#
# (12) find_files_by_analysis_num(standards.Iso,analysisNum.vec) -
#         get filenames of files with specified analysis numbers
#
# (13) find_expected_times(peaks_times.list,expNumPks=c(15,16)) -
#         find average values for Start, Rt and End times
#
# (14) generic_plot_all_raw(raw.list) -
#         plot all the raw data in a list containing raw data from different experiments
#
# (15) generic_raw_plot(raw.df,title) -
#         plot raw data using generic plot
#
# (16) group_picked_times(timesPicked) -
#         organizes the times picked from the PeakTimes function
#
# (17) identifier_1_files(files,identifier_1,cores=2) -
#         get filenames for files with specified Identifier_1 labels
#
# (18) intensity_similarity_check(vendAmpl,amplName,peakNr.vec,relDiff.thresh=0.1) -
#         check specified peaks for similarity in intensity values
#
# (19) iso_ratio_similarity(vend.df,peakNr.vec,sdC.thresh=0.1,sdO.thresh=0.1) -
#         check whether isotope ratios are within similarity criteria for specific peaks
#
# (20) move_identifier_1_files(identifier_1_files,path,identifier_1)
#         move all files with a specified identifier_1 into a folder labeled with the identifier_1
#
# (21) peak_area_trap(start.t,end.t,time.vec,int.vec) -
#         find the area of a peak using the integration trapezoidal rule (via trapz)
#
# (22) peak_times_check<-function(vend.df,expectedPeak.num=5, diff.t=10,expectedStart,expectedRt,expectedEnd) -
#         check whether peaks occur at the expected times
#
# (23) peak_times_pick(Int.mat,time.interval,z.thresh,reactionLabel) -
#         find the approximate peak start times
#
# (24) peak_vend_times(vend.df) -
#         get the peak numbers and Start, Rt, and End times for vendor data
#
# (25) peak_vend_times_all(vend.list) -
#         get peak times from a list of vendor data from several experiments
#
# (26) plot_all_peaks(start.v,end.v,t.s,v.mV,v.name) -
#          plot every peak in an experiment in a separate graph (in a shared panel)
#
# (27) plot_individual_peaks(start.t,end.t,time.vec,int.vec,peak_num,v.mv) -
#          plot an individual peak in an experiment
#
# (28) plot_ms(vendor_info.df,x_name="Rt",y_name="Intensity_All") -
#         plot mass spec data from vendor table
#
# (29) qc_num_vendPeaks(vend.df,expectedPeakNum1=15,expectedPeakNum2=16) -
#         check for the correct number of peaks
#
# (30) qc_peaks_vend(start.times) -
#         check for peak presence
#
# (31) qc_samplePeaks(intensity.mat,time.interval,expectedNum.samplePeaks,expectedStartSamples,samplePeakTimes,label) -
#         test whether expected number of sample peaks are present and if they occur at expected times
#
# (32) raw_data(file) -
#         get raw data from the specified .dxf file as a dataframe
#
# (33) raw_data_all(files) -
#         get all raw data from multiple files in a list
#
# (34) read_summary(filename) -
#         read and print a summary of data from a .dxf file
#
# (35) reference_times_check(vend.df,expectedPeak.num=5, diff.t=10, expectedStart,expectedRt,expectedEnd) -
#         determine if the expected number and peak times of reference peaks are present; return reference vendor data
#
# (36) reference_values_no_ratio(files) -
#         get isotopic reference values without ratios
#
# (37) reference_values_ratio(files) -
#         get isotopic reference values including ratios
#
# (38) remove_276(vend.df) -
#         remove the (1st sample) peak at 276 which often has some contamination for CO2
#
# (39) resistor_data(files) -
#          get resistor info for .dxf files
#
# (40) sample_peaks_pre_process(refTimesOutput,vend.df)
#         separate out sample peaks data (without any processing like flush peak removal)
#
# (41) sample_peaks_process<-function(refTimesOutput,vend.df,flushExpT=135,flushTint=10,firstSampExpT=275,firstSampTint=10) -
#         remove flush peak if present and 1st sample peak
#
# (42) sample_peaks_vend(vend.df,start.sample=326,stop.sample=825) -
#         grab the vendor info for the peaks at expected sample peak times
#
# (43) separate_by_analysis_num<-function(vend.df) -
#         separate out data by analysis number
#
# (44) sort_by_identifier_1(path) -
#         sort all .dxf files in a path into folders labeled with Identifier_1 labels
#
# (45)  standards_input_list(standardPaths.vec,standardAnalysisNums.df)
#         create input list for internal standards quality check
#
# (46) stand_lm(acceptedMeas.df) -
#         perform the linear regression for the internal standards quality check (and plot)
#
# (47) standardInputProc(path,standardAnalysisNums.vec) -
#         gets sample peak data for internal standards quality check
#
# (48) ThresholdingAlgo(y,lag,threshold,influence) -
#         peak detection algorithm (uses z-scores) to get peak (signal) start times
#
# (49) trap_area_allPks(raw.df,vend.df,mV.rawName) -
#         calculate area of all peaks in an experiment via trapz
#
# (50)  unique_identifiers(filenames) -
#         get all of the different identifier_1 labels in a set of files
#
# (51) vendor_info(file) -
#         get selected vendor info for a .dxf file: Identifier_1, Peak_Nr, Start, Rt, End, Intensity_All,rIntensity_All,
#         Ampl_44,Ampl_45,Ampl_46,d13C/12C,d18O/16O
#
# (52) vendor_info_all(files) -
#         get vendor table data from several .dxf files
#
#########################################################################

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
#' all_PA_trap: function to calculate all peak areas using the integration trapezoidal rule (via trapz)
#' @param start.vec numeric vector of peak start times
#' @param end.vec numeric vector of peak end times
#' @param time.vec numeric vector of times from the raw data
#' @param int.vec numeric vector of intensities from the raw data
#' @return a numeric vector containing all the peak areas in V*s
#' @examples
#' Usage Example
#' all_PA_trap(start.v1,end.v1,time.s,v44,pk.Nrs)
#' @export
all_PA_trap<-function(start.vec,end.vec,time.vec,int.vec,pk.Nrs){
  all_areas<-c()
  for(i in seq(1:length(start.vec))){
    all_areas<-c(all_areas,peak_area_trap(start.vec[i],end.vec[i],time.vec,int.vec))
  }
  all_areas_Vs<-all_areas
  areas.mat<-matrix(c(pk.Nrs,all_areas_Vs),ncol=2)
  areas.df<-as.data.frame(areas.mat)
  colnames(areas.df)<-c("Pk_Nr","trap_area")
  return(areas.df)
}


# (3)
#' avg_sd_d18O_standards: function to find d18O/16O averages and standard deviations for each set of files represented in input list
#' @param allStandards_d18O.list list built from output of d18O_samples_list(), which retrieves sample peaks and d18O/16O data from specified files
#' @param standNames character vector containing the names of the standards
#' @param standAcceptedVals.vec numeric vector containing the accepted values of the standards
#' @param accStandRatioSD numeric vector containing the accepted standard deviation values for each of the standards (in the same order as names in standNames)
#' @return list of length 3 with the following elements:
#'          [[1]]: dataframes of d18O/16O averages and standard deviations for each set of files from input
#'          [[2]]: accepted and measured d18O/16O dataframe
#'          [[3]]: acceptable and calculated standard deviations for d18O/16O
#' @examples
#' Usage Example
#' L1_d18O.list<-d18O_samples_list(L1fileNames.vec)
#' H1_d18O.list<-d18O_samples_list(H1fileNames.vec)
#' LW_d18O.list<-d18O_samples_list(LWfileNames.vec)
#' allStand_d18O.list<-list()
#' allStand_d18O.list[[1]]<-L1_d18O.list
#' allStand_d18O.list[[2]]<-H1_d18O.list
#' allStand_d18O.list[[3]]<-LW_d18O.list
#' AvgSD_d18O<-avg_sd_d18O_standards(allStand_d18O.list)
#' @export
avg_sd_d18O_standards<-function(allStandards_d18O.list,standNames=c("L1","H1","LW"),standAcceptedVals.vec=c(-8.55,4.85,-3.85),accStandRatioSD=c(0.2,0.2,0.2)){
  # make df for avgs (avgs of all d18O in all files for particular standard)
  standAccepted.mat<-matrix(standAcceptedVals.vec,ncol=1)
  standAccepted.df<-as.data.frame(standAccepted.mat)
  # start df -- will add average of all measured values
  d18OstandVals_Avgs.df<-data.frame(matrix(rep(NA,6),ncol=2))
  colnames(d18OstandVals_Avgs.df)<-c("accepted_d18O16O","measured_d18O16O")
  rownames(d18OstandVals_Avgs.df)<-standNames
  d18OstandVals_Avgs.df[,1]<-standAccepted.df

  # make df for SDs
  SDAccepted.mat<-matrix(accStandRatioSD,ncol=1)
  SDAccepted.df<-as.data.frame(SDAccepted.mat)
  # start df -- will add SD of all measured values
  d18OstandVals_SDs.df<-data.frame(matrix(rep(NA,6),ncol=2))
  colnames(d18OstandVals_SDs.df)<-c("acceptable_SD_d18O16O","calculated_SD_d18O16O")
  rownames(d18OstandVals_SDs.df)<-standNames
  d18OstandVals_SDs.df[,1]<-accStandRatioSD
  # initialize lists for return vals
  ret_files_d18.list<-list()
  ret_d18.list<-list()
  d18O.list<-allStandards_d18O.list
  list.len<-length(d18O.list)
  alld18OFile.list<-list()
  alld18OFile.index<-1
  for(i in seq(1,list.len)){
    # for the different standards
    d18OFile.mat<-matrix(rep(NA,list.len*3),ncol=3)
    d18OFile.df<-as.data.frame(d18OFile.mat)
    colnames(d18OFile.df)<-c("file_id","avg_d18O16O","sd_d18O16O")
    numSampDFs<-length(d18O.list[[i]])
    alld18O.vec<-c()
    for(j in seq(1,numSampDFs)){
      # for each file in standard set
      dlist<-(d18O.list[[i]][j])
      # get sample peaks d18O/16O for one file
      dlist.df<-as.data.frame(dlist)
      # get the filename
      fileN<-dlist.df$file_id[1]
      # get the d18O data for the file
      d18O.vec<-as.numeric(dlist.df$d18O.d16O)
      # put in vector of all d18O for the given standard (i iteration, over the 3 files for a standard)
      alld18O.vec<-c(alld18O.vec,d18O.vec)
      # get avg d18O for each file
      avg<-mean(d18O.vec)
      # get sd for each file
      sd<-sd(d18O.vec)
      # put file_id, avg, and sd for each file in df
      d18OFile.df[j,]<-data.frame(fileN,avg,sd)
    }
    # add to list for dfs of file sets
    alld18OFile.list[[alld18OFile.index]]<-d18OFile.df
    alld18OFile.index<-alld18OFile.index+1

    # avg of all d18O in files for given standard
    d18OstandAvg<-mean(alld18O.vec)
    # put in df to return
    d18OstandVals_Avgs.df$`measured_d18O16O`[i]<-d18OstandAvg

    # SD
    sdAllStandFiles<-sd(alld18O.vec)
    d18OstandVals_SDs.df$`calculated_SD_d18O16O`[i]<-sdAllStandFiles

    # return d18O.df in list
    ret_files_d18.list[[i]]<-d18OFile.df
  }
  # add values to return list
  ret_d18.list[[1]]<-ret_files_d18.list
  ret_d18.list[[2]]<-d18OstandVals_Avgs.df # all avgs
  ret_d18.list[[3]]<-d18OstandVals_SDs.df # all SDs

  return(ret_d18.list)
}


# (4)
#' average_times: function that finds the averages of each peak time for vendor data (called in find_expected_times())
#' @param peakTimes.list list of peak numbers and times from vendor table (output from peak_vend_times_all())
#' @param nPeaks number of peaks
#' @return dataframe of mean values for each peak number in the input vendor data
#' @examples
#' Usage Example
#' peak15avg.df<-average_times(peaks15.list,15)
#' @export
average_times<-function(peakTimes.list,nPeaks){
  nTimes<-length(peakTimes.list)
  nCols<-4
  mean.df<-as.data.frame(matrix(rep(NA,nCols*nPeaks),ncol=nCols))
  colnames(mean.df)<-c("Peak_Nr","Start_Avg","Rt_Avg","End_Avg")
  # iterate over peaks
  for(m in seq(1,nPeaks)){
    # peak number
    peak.num<-round(m,0)
    # initialize time vecs
    start.v<-c()
    rt.v<-c()
    end.v<-c()
    # get all the times for that peak
    for(k in seq(1,nTimes)){
      # for each peak, grab all starts, rt, and ends
      # get start, add to vec
      currSt<-as.numeric(peakTimes.list[[k]]$Start[m])
      start.v<-c(start.v,currSt)
      # get Rt, add to vec
      currRt<-as.numeric(peakTimes.list[[k]]$Rt[m])
      rt.v<-c(rt.v,currRt)
      # get End, add to vec
      currEnd<-as.numeric(peakTimes.list[[k]]$End[m])
      end.v<-c(end.v,currEnd)
    }
    # take averages for each peak, store
    start.avg<-mean(start.v)
    rt.avg<-mean(rt.v)
    end.avg<-mean(end.v)

    peak.v<-as.character(c(peak.num,round(start.avg,4),round(rt.avg,4),round(end.avg,4)))
    peak.m<-matrix(peak.v,ncol=4)
    peak.df<-as.data.frame(peak.m)
    mean.df[peak.num,]<-peak.df
  }
  return(mean.df)
}


# (5)
#' build_expected_times: function to build expected times data for reference and sample peaks
#' @param ref.ind numeric vector containing the indices of reference peaks in avgs.df
#' @param sample.ind numeric vector containing the indices of sample peaks in avgs.df
#' @param avgs.df dataframe of average values for peaks times
#' @return list of 2 elements: first for reference peaks expected times and the second for sample peaks
#' @examples
#' Usage Example
#' refPk.vec<-c(1,2,4,5,16)
#' samplePk.vec<-c(3,6,7:15)
#' pkAvg.df<-average_times(peaks15.list,15)
#' expTimes.list<-build_expected_times(ref.ind,sample.ind,avgs.df)
#' @export
build_expected_times<-function(ref.ind,sample.ind,avgs.df){
  # get expected start times for reference and sample peaks
  expectedRefStart.vec<-as.numeric(avgs.df$Start_Avg[ref.ind])
  expectedSampleStart.vec<-as.numeric(avgs.df$Start_Avg[sample.ind])
  # get expected Rts for reference and sample peaks
  expectedRefRt.vec<-as.numeric(avgs.df$Rt_Avg[ref.ind])
  expectedSampleRt.vec<-as.numeric(avgs.df$Rt_Avg[sample.ind])
  # get expected end times for reference and sample peaks
  expectedRefEnd.vec<-as.numeric(avgs.df$End_Avg[ref.ind])
  expectedSampleEnd.vec<-as.numeric(avgs.df$End_Avg[sample.ind])
  # make reference df
  expectedRefs.mat<-matrix(c(ref.ind,expectedRefStart.vec,expectedRefRt.vec,expectedRefEnd.vec),ncol=4)
  expectedRefs.df<-as.data.frame(expectedRefs.mat)
  colnames(expectedRefs.df)<-c("Ref_Peak_Nr","Expected_Start","Expected_Rt","Expected_End")
  # make sample df
  expectedSamps.mat<-matrix(c(sample.ind,expectedSampleStart.vec,expectedSampleRt.vec,expectedSampleEnd.vec),ncol=4)
  expectedSamps.df<-as.data.frame(expectedSamps.mat)
  colnames(expectedSamps.df)<-c("Sample_Peak_Nr","Expected_Start","Expected_Rt","Expected_End")
  # put ref and sample df into list
  expected.list<-list()
  expected.list[[1]]<-expectedRefs.df
  expected.list[[2]]<-expectedSamps.df
  # return list of dfs
  return(expected.list)
}


# (6)
#' check_picked_and_vendor_peaks: function that checks peak times and numbers from vendor table against picked peaks
#' @param vend.df dataframe of vendor data
#' @param int.mat 1-col matrix of intensity values
#' @param time.vec numeric vector of time values
#' @param z.thresh z-value threshold for ThresholdingAlgo
#' @param graphLab label for graph of picked times
#' @param timeThresh threshold for time differences between picked and vendor peaks
#' @param expPeaks vector containg two values for expected numbers of peaks
#' @return list of 3 elements: [[1]] dataframe indicating whether the number of peaks in
#'         vendor table match number picked; [[2]] dataframe containing number of vendor peaks
#'         and picked peaks; [[3]] dataframe of peak times for vendor and picked peaks
#' @examples
#' Usage Example
#' int.mat<-matrix(v44,ncol=1) # for ThresholdingAlgo
#' time.vec<-time.s
#' graphLab<-"170525_NaHCO3 L + NaCl L"
#' numPeakCheck<-check_picked_and_vendor_peaks(vend1.df,int.mat,time.s,5,graphLab,2)
#' @export
check_picked_and_vendor_peaks<-function(vend.df,int.mat,time.vec,z.thresh=5,graphLab,timeThresh=2,expPeaks=c(15,16)){
  peaks.list<-list()
  numPeaks<-FALSE
  times<-FALSE
  pickedPeakTimes<-peak_times_pick(int.mat,time.vec,z.thresh,graphLab)
  numPickedPeaks<-length(pickedPeakTimes)
  vendTimes<-peak_vend_times(vend.df)
  numVendTimes<-length(vendTimes$Start)
  # check lengths the same
  if(numPickedPeaks==numVendTimes){
    numPeaks<-TRUE
    print("Number of peaks in vendor table and number of peaks detected match")
  }
  # check times similar (within 2 secs?)
  timeDiff.vec<-abs(pickedPeakTimes-vendTimes$Start)
  numDiffs<-length(timeDiff.vec)
  indexOff<-c()
  for(i in seq(1,numDiffs)){
    if(timeDiff.vec[i]>2){
      print(paste("significant time difference detected at index ", i,sep=""))
      indexOff<-c(indexOff,i)
    }
  }
  if(length(indexOff)==0){ # no significant time differences between peaks
    times<-TRUE
    print("No siginificant time differences between vendor peaks and picked peaks")
  }
  # boolean vector for whether number of peaks detected matches vendor table and if times are similar
  peaks_times.mat<-matrix(c(numPeaks,times),ncol=2)
  peaks_times.df<-as.data.frame(peaks_times.mat)
  colnames(peaks_times.df)<-c("Num_peaks_match","Times_match")
  peaks_times.df
  # vector for peak numbers in vendor table and picked out
  peak_nums.mat<-matrix(c(numVendTimes,numPickedPeaks),ncol=2)
  peak_nums.df<-as.data.frame(peak_nums.mat)
  colnames(peak_nums.df)<-c("Num_Peaks_Vendor","Num_Peaks_Picked")
  peak_nums.df
  # df for peak times
  peak_times.mat<-matrix(c(vendTimes$Start,pickedPeakTimes),ncol=2)
  peak_times.df<-as.data.frame(peak_times.mat)
  colnames(peak_times.df)<-c("Vendor_peak_start","Picked_peak_start")
  #########
  # check if num peaks expected
  expPeaks<-c(15,16)
  NPV<-peak_nums.df$Num_Peaks_Vendor
  NPP<-peak_nums.df$Num_Peaks_Picked
  npv.exp<-NPV %in% expPeaks
  npp.exp<-NPP %in% expPeaks
  ret.exp<-matrix(c(npv.exp,npp.exp),ncol=2)
  retExp.df<-as.data.frame(ret.exp)
  colnames(retExp.df)<-c("Vendor_Exp_Pk_Nums","Picked_Exp_Pk_Nums")
  ########
  peaks.list[[1]]<-peaks_times.df
  peaks.list[[2]]<-peak_nums.df
  peaks.list[[3]]<-peak_times.df
  peaks.list[[4]]<-retExp.df

  return(peaks.list)
}


# (7)
#' d18O_samples_list: function to get sample peak and d18O data for the specified files
#' @param standardFilesNames character vector of filenames containing standards (of same type)
#' @return list of dataframes containing sample peaks and d18O data
#' @examples
#' Usage example
#' LWstandFiles<-find_files_by_analysis_num(LWiso,LWstandardAnalysisNum)
#' LW_d18O.list<-d18O_samples_list(LWstandardFileNames)
#' @export
d18O_samples_list<-function(standardFileNames,expPk.num=5,exp.df){
  # get expected times
  expSt.vec<-exp.df$Expected_Start
  expEnd.vec<-exp.df$Expected_End
  expRt.vec<-exp.df$Expected_Rt

  # get vendor data and d1O values
  isoVend<-vendor_info_all(standardFileNames)
  d18Osample.list<-list()
  numSampleVend<-length(isoVend)
  for(i in seq(1,numSampleVend)){
    # get sample peaks
    refTC<-reference_times_check(isoVend[[i]],expectedPeak.num=expPk.num,expectedStart=expSt.vec,
                                 expectedRt=expRt.vec,expectedEnd=expEnd.vec)
    sampleVend<-sample_peaks_process(refTC,isoVend[[i]])
    # get d18O/d16O data
    sample.d18O<-as.numeric(sampleVend$`d18O/16O`)
    # get peak numbers
    sample.pkNr<-as.numeric(sampleVend$Peak_Nr)
    sampleNumPeaks<-length(sample.pkNr)
    # put data together
    sample.mat<-matrix(c(rep(standardFileNames[i],sampleNumPeaks),sample.pkNr,sample.d18O),ncol=3)
    d18O.df<-as.data.frame(sample.mat)
    colnames(d18O.df)<-c("file_id","Pk_Nr","d18O/d16O")
    # add to list
    d18Osample.list[[i]]<-d18O.df
  }
  return(d18Osample.list)
}


# (8)
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


# (9)
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


# (10)
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


# (11)
#' filterSpecSignals: function that filters out possible false positives in peak detection
#' @param threshAlgoList list returned by ThresholdingAlgo()
#' @return list that has been filtered for possible false positive signals (3 consecutive 1s or less in signal matrix)
#' @examples
#' Usage Example
#' filterSpecSignals(ThreshAlgoList)
#' @export
filterSpecSignals<-function(threshAlgoList){
  specSignals<-threshAlgoList$signals
  newSpecSignals<-specSignals
  ones<-which(specSignals==1)
  ones
  # loop through to find consecutive ones, then count and keep only peaks with 4 or more consecutive ones
  num_ones<-length(ones)
  consec.vec<-c()
  # find consecutive indices, make groups representing signals
  list.index<-1
  consec.list<-c()
  for(i in seq(1,num_ones-1)){
    if(ones[i+1]-ones[i]==1){
      #print("consecutive 1s found")
      consec.vec<-c(consec.vec,ones[i])
    }else{
      # end of list
      consec.vec<-c(consec.vec,ones[i])
      #print(consec.vec)
      # add to list and reset vector
      consec.list[[list.index]]<-consec.vec
      consec.vec<-c()
      list.index<-list.index+1
    }
  }
  num_consecs<-length(consec.list)
  num_consecs
  for(i in seq(1,num_consecs)){
    consec.mat<-matrix(consec.list[[i]],nrow=1)
    consec.mat
    if(length(consec.list[[i]])<3){
      # replace 1s with 0 in signals
      newSpecSignals[consec.mat[1,]]<-0
    }
  }
  newSpecSignals
  threshAlgoList$signals<-newSpecSignals
  return(threshAlgoList)
}


# (12)
#' find_files_by_analysis_num: function that gets filenames of files with specified analysis numbers
#' @param standards.Iso continuous flow data for standards files (of the same standard)
#' @param analysisNum.vec: numeric vector of analysis numbers identifying which standard files to pull
#' @return character vector of file names containing the specified analyses
#' @examples
#' Usage Example
#' LWstandardFileNames<-find_files_by_analysis_num(LWiso,LWstandardAnalysisNum)
#' @export
find_files_by_analysis_num<-function(standards.Iso,analysisNum.vec){
  fileInfo<-iso_get_file_info(standards.Iso)
  files<-fileInfo$file_id
  AnalysisNums<-fileInfo$Analysis
  numFiles<-length(files)
  numAnalyses<-length(analysisNum.vec)
  fileNames<-c()
  for(i in seq(1,numAnalyses)){
    if(analysisNum.vec[i] %in% AnalysisNums){
      matchInd<-match(analysisNum.vec[i],AnalysisNums)
      # get filename
      standMatch<-files[matchInd]
      fileNames<-c(fileNames,standMatch)
    }
  }
  return(fileNames)
}


# (13)
#' find_expected_times: function to find average values for start, Rt and end times
#' @param peaks_times.list a list of dataframes of peak numbers and times data
#' @return if data for both peak numbers present, a list of dataframes of expected times for data
#' with both peak numbers, if data of just one of the expected peak numbers is found, it returns a dataframe of expected times
#' @examples
#' Usage Example
#' pT.df<-find_expected_times(peaksTimes.list)
#' @export
find_expected_times<-function(peaks_times.list,expNumPks=c(15,16)){
  nTimes<-length(peaks_times.list)
  nPeaks.vec<-c()
  # find any elements that have any number but 15 or 16 peaks
  not15or16<-c()
  corrNumPeaks<-c()
  peaks15<-c()
  peaks16<-c()
  for(i in seq(1,nTimes)){
    nPeaks<-length(peaks_times.list[[i]]$Peak_Nr)
    if((nPeaks!=expNumPks[1]) && (nPeaks!=expNumPks[2])){ #expected number of peaks---add as arg**
      not15or16<-c(not15or16,i)
      nPeaks.vec<-c(nPeaks.vec,nPeaks)
      print(paste("element with unexpected number of peaks found. List index ",i," will be removed from the average calculation",sep=""))
    } else{
      corrNumPeaks<-c(corrNumPeaks,nPeaks)
      # list indices with 15 peaks
      if(nPeaks==expNumPks[1]){
        peaks15<-c(peaks15,i)
      }
      # list indices with 16 peaks
      if(nPeaks==expNumPks[2]){
        peaks16<-c(peaks16,i)
      }
    }
  }
  # separate peaks_times.list by number of peaks to find the averages
  p15present<-length(peaks15)>0
  p16present<-length(peaks16)>0

  # for 15 peaks
  if(length(peaks15)>0){
    peaks15.list<-list()
    for(i in seq(1,length(peaks15))){
      ind15<-peaks15[i]
      peaks15.list[[i]]<-peaks_times.list[[ind15]]
    }
    # call average_times() function to average the times
    peaks15.list
    peak15avg.df<-average_times(peaks15.list,expNumPks[1])
  }

  # for 16 peaks
  if(length(peaks16)>0){
    peaks16.list<-list()
    for(i in seq(1,length(peaks16))){
      ind16<-peaks16[i]
      peaks16.list[[i]]<-peaks_times.list[[ind16]]
    }
    # call average_times() function to average the times
    peak16avg.df<-average_times(peaks16.list,expNumPks[2])
  }

  # put them together?
  if(p16present && p15present){
    # combine them
    list15and16<-list()
    list15and16[[1]]<-peak15avg.df
    list15and16[[2]]<-peak16avg.df
    return(list15and16)
  } else{
    if(p16present){
      # return only that avg
      return(peak16avg.df)
    } else{ #p15present
      # return only that avg
      return(peak15avg.df)
    }
  }
}


# (14)
#' generic_plot_all_raw: function that plots all the raw data in a list containing raw data from different experiments
#' @param raw.list list whose elements are the raw.df to be plotted
#' @examples
#' Usage Example
#' generic_plot_all_raw(rawList)
#' @export
generic_plot_all_raw<-function(raw.list){
  raw.length<-length(raw.list)
  par(mfrow=c(2,3))
  for(i in seq(1,raw.length)){
    # get title
    raw.dat<-raw.list[[i]]
    raw.title<-raw.dat$file_id[1]
    raw.title<-paste("Index: ",i," Raw Data ",raw.title,sep="")
    raw.title
    # plot
    generic_raw_plot(raw.dat,raw.title)
  }
  par(mfrow=c(1,1))
}


# (15)
#' generic_raw_plot: function to plot raw data using generic plot
#' @param raw.df full dataframe of raw data
#' @param title title of the plot
#' @examples
#' Usage Example
#' generic_raw_plot(rawDat.df,"170525_NaHCO3 L + NaCl L_.dxf")
#' @export
generic_raw_plot<-function(raw.df,title){
  #par(mfrow=c(1,1))
  # get intensity data for each mass
  v44<-raw.df$v44.mV
  v45<-raw.df$v45.mV
  v46<-raw.df$v46.mV
  # get time data
  time.s<-raw.df$time.s
  # this way, have to plot the tallest one first to get all in same window
  plot(time.s,v46,type="l",col="magenta",ylab="Intensity (mV)",xlab="Time(s)")#v44_v45_v46.mV
  lines(time.s,v45,type="l",col="green")
  lines(time.s,v44,type="l",col="blue")
  title(main=title)
}

##### FIXME
# get files by analysis numbers from iso file directory
get_analysis_nums<-function(path,analysis.vec){
  allFiles<-all_filenames(path)
  fileInfo.df<-file_info(allFiles)
  retFilesInd<-which(fileInfo.df$Analysis %in% analysis.vec)
  retFiles<-fileInfo.df$file_id[retFilesInd]
  return(retFiles)
}

# (16)
#' group_picked_times: function that organizes the times picked from the PeakTimes function
#' into single average times for peaks
#' @param timesPicked the times returned in a list from ThresholdingAlgo
#' @return the average times from each "group" (groups returned from ThresholdingAlgo$signals that needs to be separated)
#' @examples
#' Usage Example
#' group_picked_times(timesPicked.list)
#' @export
group_picked_times<-function(timesPicked){
  group.list<-list()
  listIndex<-1
  timeGroup<-c()
  for(i in seq(1,length(timesPicked))){
    if(i==length(timesPicked)){
      #put the last element in the last group
      #inGroup<-abs(timesPicked[i]-timesPicked[(i-1)])<3
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


# (17)
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
  #identifier_1_files_ind<-which(identifier_1_files_ind) #which are TRUE
  # get the identifier_1 file names
  identifier_1_files<-files[identifier_1_files_ind]
  i1.df<-as.data.frame(identifier_1_files)
  colnames(i1.df)<-c(identifier_1)
  return(i1.df)
}


# (18)
#' intensity_similarity_check: function that checks specified peaks for similarity in intensity values
#' @param vendAmpl vector of amplitudes (as from vendor table)
#' @param amplName label of the amplitude (as from vendor table)
#' @param peakNr.vec vector of peak numbers of interest
#' @param relDiff.thresh threshold for relative difference between intensities
#' @return list that contains TRUE/FALSE if intensities are within specified threshold and a dataframe of the peak numbers and intensities analyzed
#' @examples
#' Usage Example
#' intensity_similarity_check(vend.df$Ampl_44,"Ampl_44",c(1,2,4,5,16),0.1)
#' @export
intensity_similarity_check<-function(vendAmpl,amplName,peakNr.vec,relDiff.thresh=0.1,ID1){
  intCheck.list<-list()
  # analyze intensity similarity by relative difference
  int.check<-FALSE
  relDiff<-c()
  for(i in seq(2,length(vendAmpl))){
    curr.relDiff<-abs(vendAmpl[i]-vendAmpl[i-1])/mean(vendAmpl)
    if(curr.relDiff<-relDiff.thresh){
      relDiff<-c(relDiff,TRUE)
    }
  }
  sumRelDiff<-sum(relDiff)
  if(sumRelDiff==length(relDiff)){
    int.check<-TRUE
    #print("intensities within similarity criteria")
    intCheck.list[[1]]<-int.check
  } else{
    int.check<-FALSE
    print("intensities not within similarity criteria")
    print(paste("ID1: ",ID1,sep=""))
    #return(int.check)
    intCheck.list[[1]]<-int.check
  }
  # dataframe of intensity values at peaks
  intCheck.mat<-matrix(c(peakNr.vec,vendAmpl[peakNr.vec]),ncol=2)
  intCheck.df<-as.data.frame(intCheck.mat)
  colnames(intCheck.df)<-c("Peak_Nr",amplName)
  # add to list
  intCheck.list[[2]]<-intCheck.df
  return(intCheck.list)
}


# (19)
#' iso_ratio_similarity: function that checks whether isotope ratios are within similarity criteria for specific peaks
#' @param vend.df dataframe of vendor table data
#' @param peakNr.vec vector of peak numbers to be analyzed
#' @param sdC.thresh threshold for the standard deviation of d13C/12C (default=0.1)
#' @param sdO.thresh threshold for the standard deviation of d18O/16O (default=0.1)
#' @return list with boolean values if ratios are within thresholds and the standard deviation values
#' @examples
#' Usage Example
#' refIsoRatios<-iso_ratio_similarity(vi.df,c(1,2,4,5,16),sdC.thresh=0.1,sdO.thresh=0.1)
#' @export
iso_ratio_similarity<-function(vend.df,peakNr.vec,sdC.thresh=0.1,sdO.thresh=0.1){
  isoR.list<-list()
  d13C<-vend.df$`d13C/12C`
  d13C.ref<-as.numeric(d13C[peakNr.vec])
  d13C.ref
  d18O<-vend.df$`d18O/16O`
  d18O.ref<-as.numeric(d18O[peakNr.vec])

  # sd of d13C/12C
  sd.13C<-sd(d13C.ref)
  sd.C<-FALSE
  if(sd.13C<sdC.thresh){
    sd.C<-TRUE
    #print("SD 13C/12C within accepted threshold")
  } else{
    print(paste("SD 13C/12C not within accepted threshold for Analysis ",
                analysis,sep=""))

  }

  # sd of 18O/16O
  sd.18O<-sd(d18O.ref)
  sd.O<-FALSE
  if(sd.18O<sdO.thresh){
    sd.O<-TRUE
    #print("SD 18O/16O within accepted threshold")
  } else{
    print("SD 18O/16O not within accepted threshold for Analysis ",
          analysis,sep=",")
  }

  # dataframe for boolean vals for within threshold
  sdBool.mat<-matrix(c(sd.C,sd.O),ncol=2)
  sdBool.df<-as.data.frame(sdBool.mat)
  colnames(sdBool.df)<-c("d13C/12C","d18O/16O")
  isoR.list[[1]]<-sdBool.df

  # dataframe for numeric value of sds
  sdnum.mat<-matrix(c(sd.13C,sd.18O),ncol=2)
  sdnum.df<-as.data.frame(sdnum.mat)
  colnames(sdnum.df)<-c("SD d13/12C","SD d18O/16O")
  isoR.list[[2]]<-sdnum.df

  return(isoR.list)
}


# (20)
#' move_identifier_1_files: function that moves all files with a specified identifier into a folder labeled with the Identifier_1 value
#' @param identifier_1_files files that contain the desired identifier
#' @param path location of the file to be copied
#' @param identifier_1 the identifier_1 of the files to be copied
#' @examples
#' Usage example
#' move_Identifier_1_files(identifier_1_d_files,pathc,identifier_1_d)
#' @export
move_identifier_1_files<-function(identifier_1files,path,identifier_1){
  num_files<-length(identifier_1files)
  for(i in seq(1:num_files)){
    orig_dir<-paste(path,"/",identifier_1files[i],sep="")
    #orig.dir
    new_dir<-paste(path,"/",identifier_1,"/",identifier_1files[i],sep="")
    new_dir
    file_move(orig_dir,new_dir)
  }
}

# (21)
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


# (22)
#' peak_times_check: function that checks whether peaks occur at the expected times
#' @param vend.df dataframe of vendor table info
#' @param expectedPeak.num expected number of reference peaks (default=5)
#' @param diff.t acceptable time difference between expected peak times and actual peak times (default=5s)
#' @param expectedStart vector of expected peak Start times (default for CO2)
#' @param expectedRt vector of expected peak Rt (default for CO2)
#' @param expectedEnd vector of expected peak End times (default for CO2)
#' @return list containing boolean values indicating whether the expected number of reference peaks was found at the expected times and a dataframe
#' containing the peak numbers and times
#' @examples
#' Usage Example
#' expTimes.list<-build_expected_times(refInd.vec,sampInd.vec,avgs.df)
#' peak_times_check(vendor.list[[10]],expectedStart=expTimes.list[[1]]$Expected_Start,expectedRt=expTimes.list[[1]]$Expected_Rt,expectedEnd=expTimes.list[[1]]$Expected_End) # using default values
#' @export
peak_times_check<-function(vend.df,expectedPeak.num, diff.t=10,expectedStart,expectedRt,expectedEnd){
  refs.list<-list()
  peak.nums<-as.numeric(vend.df$Peak_Nr)
  num.peaks<-length(peak.nums)
  # times from vendor table
  start.times<-as.numeric(vend.df$Start)
  start.times
  Rts<-as.numeric(vend.df$Rt)
  end.times<-as.numeric(vend.df$End)
  # intialize loop data
  start.ind<-c()
  rt.ind<-c()
  end.ind<-c()
  ref.ind<-c()
  expectedIndex<-1
  for(i in seq(1,num.peaks)){ # check times and grab reference peaks if detected
    # initialize ref bools
    start.ref<-FALSE
    rt.ref<-FALSE
    end.ref<-FALSE
    # check start
    if(abs(start.times[i]-expectedStart[expectedIndex])<diff.t){ # ??
      # matching start time
      start.ref<-TRUE
      #print(paste("Matching Start time detected at index ",i,", time ",round(start.times[i],3)," s",sep=""))
      start.ind<-c(start.ind,i)
    }
    # check Rt
    if(abs(Rts[i]-expectedRt[expectedIndex])<diff.t){
      # matching Rt
      rt.ref<-TRUE
      #print(paste("Matching Rt detected at index ",i,", time ",round(Rts[i],3)," s",sep=""))
      rt.ind<-c(rt.ind,i)
    }
    # check End
    if(abs(end.times[i]-expectedEnd[expectedIndex])<diff.t){
      # matching end time
      end.ref<-TRUE
      #print(paste("Matching End time detected at index ",i,", time ",round(end.times[i],3)," s",sep=""))
      end.ind<-c(end.ind,i)
    }
    if(sum(start.ref,rt.ref,end.ref==3)){#all times match for peak as ref
      print(paste("Start, Rt, and End times match an expected peak for Peak_Nr ",peak.nums[i],sep=""))
      # add to ref peak ind
      ref.ind<-c(ref.ind,i)
      expectedIndex<-expectedIndex+1
    }
  } # end for
  # check if number of reference peaks found matches the expected number
  numRefs<-FALSE
  foundAll<-(length(ref.ind)==expectedPeak.num)
  if(foundAll){   # return reference peak numbers and times
    numRefs<-TRUE
    print("Expected number of peaks detected at expected times.")
    refs.list[[1]]<-numRefs
  } else{
    print("Did not detected the expected number of peaks at expected times")
    refs.list[[1]]<-numRefs
  }
  ret<-cbind(vend.df$Peak_Nr[ref.ind],start.times[ref.ind],Rts[ref.ind],end.times[ref.ind])
  ret.df<-as.data.frame(ret)
  colnames(ret.df)<-c("Pk_Nr","Start","Rt","End")
  refs.list[[2]]<-ret.df
  return(refs.list)
}


# (23)
#' peak_times_pick: function that gives the approximate peak start times
#' @param Int.mat matrix of intensity values within the desired interval
#' @param time.interval time values in the desired interval
#' @param z.thresh z-value threshold for ThresholdAlgo
#' @param reactionLabel label for the reaction (for plot generated showing peak start times)
#' @return vector of average peak start times as found by ThresholdingAlgo and groupPickedTimes
#' @examples
#' Usage Example
#' peak_times_pick(v44Int.mat,interval.times)
#' @export
# Int.mat: intensity values within the desired time interval
# time.interval: time values in the desired interval
peak_times_pick<-function(Int.mat,time.interval,z.thresh,reactionLabel){
  firstAlgo<-ThresholdingAlgo(Int.mat,lag=10,threshold=z.thresh,influence=0.5)
  #(testAlgo$signals)
  filteredAlgo<-filterSpecSignals(firstAlgo)
  #print(filteredAlgo$signals)
  plot(Int.mat,type="l")
  title(main=paste("Peak times detected:",reactionLabel,sep=""))
  x<-which(filteredAlgo$signals==1) # what about -1?
  y<-rep(0,length(x))
  points(x,y,col="red",pch=20)
  timesPicked<-time.interval[x]
  timesPicked
  timeGroups<-groupPickedTimes(timesPicked)
  return(timeGroups)
}


#(24)
#' peak_vend_times: function that returns the peak numbers and times for a given vendor table
#' @param vend.df dataframe of vendor data
#' @return dataframe containing peak numbers, start, retention times, and end times
#' @examples
#' Usage example
#' pkTimes.df<-peak_vend_times(v1.df)
#' @export
peak_vend_times<-function(vend.df){
  # get peak nums and start, end, rt
    peakNr<-as.numeric(vend.df$Peak_Nr)
    start<-as.numeric(vend.df$Start)
    Rt<-as.numeric(vend.df$Rt)
    end<-as.numeric(vend.df$End)
    peak.mat<-matrix(c(peakNr,start,Rt,end),ncol=4)
    peak.df<-as.data.frame(peak.mat)
    colnames(peak.df)<-c("Peak_Nr","Start","Rt","End")
  return(peak.df)
}


# (25)
#' peak_vend_times_all: function that returns all peak times from the vendor list of several experiments
#' @param vend.list list of vendor tables from each experiment
#' @return list of peaks and times for each experiment
#' @examples
#' Usage Example
#' peaks_and_times_all(vi.list)
#' @export
peak_vend_times_all<-function(vend.list){
  # get all peak nums and start, end, rt
  peaks_and_times.list<-list()
  numData<-length(vend.list)
  for(i in seq(1,numData)){
    peakNr<-as.numeric(vend.list[[i]]$Peak_Nr)
    start<-as.numeric(vend.list[[i]]$Start)
    Rt<-as.numeric(vend.list[[i]]$Rt)
    end<-as.numeric(vend.list[[i]]$End)
    intRef.mat<-matrix(c(peakNr,start,Rt,end),ncol=4)
    intRef.mat
    intRef.df<-as.data.frame(intRef.mat)
    colnames(intRef.df)<-c("Peak_Nr","Start","Rt","End")
    intRef.df
    peaks_and_times.list[[i]]<-intRef.df
  }
  return(peaks_and_times.list)
}


# (26)
#' plot_all_peaks: Function that plots every peak in an experiment in a separate graph (in shared panel)
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
  par(mfrow=c(4,4))
  for(i in seq(1:length(start.v))){
    plot_individual_peaks(start.v[i],end.v[i],t.s,v.mV,i,v.name)
  }
  par(mfrow=c(1,1))
}


# (27)
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
  #return(peak.df)
}


# (28)
#' plot_ms: Function to plot mass spec data from vendor table
#' @param vendor_info.df dataframe of vendor info for only one experiment
#' @param x_name name for desired x units for ms plot from vendor_info.df (default Rt)
#' @param y_name name for desired y units for ms plot from vendor_info.df (default rIntensity_All)
#' @return PlotSpec plot of the specified columns from vendor_info.df
#' @examples
#' Usage example
#' plot_ms(vend.df,"Rt","Intensity_All")
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


# (29)
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


# (30)
#' qc_peaks_vend: Function to check for peak presence
#' @param start.times vector of start times for each peak in the data
#' @return Boolean that indicates whether peaks are present
#' @examples
#' Usage example
#' qc_peaks_vend(start_times.vec)
#' @export
qc_peaks_vend<-function(start.times){
  start.times<-as.numeric(start.times)
  if(length(start.times)>=0){ # What do start and end look like if there are no peaks? blank or NA?
    peaks=TRUE
  }
  else{
    peaks=FALSE
  }
  return(peaks)
}


# (31)
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
qc_samplePeaks<-function(intensity.mat,time.interval,expectedNum.samplePeaks,expectedStartSamples,samplePeakTimes,label){
  peakTimes<-peakTimes(intensity.mat,time.interval,z.thresh=5,label)
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


# (32)
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


# (33)
#' raw_data_all: function to get all raw data from multiple files
#' @param files vector of filenames to get raw data from
#' @return a list containing raw data for each file
#' @examples
#' Usage example
#' rawList<-raw_data_all(files)
#' @export
raw_data_all<-function(files){
  raw.list<-list()
  for(i in seq(1,length(files))){
    raw.dat<-raw_data(files[i])
    raw.list[[i]]<-raw.dat
  }
  return(raw.list)
}


# (34)
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


# (35)
#' reference_times_check: function to determine if the expected number and start,end rt times of reference peaks are present
#' @param vend.df dataframe of vendor data
#' @param expectedPeak.num number of expected peaks
#' @param diff.t acceptable difference between peak times and expected times
#' @param expectedStart vector of expected start times
#' @param expectedRt vector of expected retention times
#' @param expectedEnd vector of expected end times
#' @return list of 3 elements: [[1]] boolean-if peaks are at expected times; [[2]] peak numer and times vendor data
#'        [[3]] full vendor data for reference peaks
#' @examples
#' Usage example
#' refTC<-reference_times_check(vi.list[[1]])
#' @export
reference_times_check<-function(vend.df,expectedPeak.num, diff.t=10,expectedStart=c(27.1700,67.0193,166.5653,206.4920,843.3227),
                                expectedRt=c(47.3888,87.1840,186.7608,226.5096,863.5029),expectedEnd=c(50.4929,90.3422,189.8572,229.6987,866.6069)){
  refs.list<-list()
  peak.nums<-as.numeric(vend.df$Peak_Nr)
  num.peaks<-length(peak.nums)
  # times from vendor table
  start.times<-as.numeric(vend.df$Start)
  Rts<-as.numeric(vend.df$Rt)
  end.times<-as.numeric(vend.df$End)
  # intialize loop data
  start.ind<-c()
  rt.ind<-c()
  end.ind<-c()
  ref.ind<-c()
  expectedIndex<-1
  for(i in seq(1,num.peaks)){ # check times and grab reference peaks if detected
    # initialize ref bools
    start.ref<-FALSE
    rt.ref<-FALSE
    end.ref<-FALSE
    # check start
    if(abs(start.times[i]-expectedStart[expectedIndex])<diff.t){ # ??
      # matching start time
      start.ref<-TRUE
      start.ind<-c(start.ind,i)
    }
    # check Rt
    if(abs(Rts[i]-expectedRt[expectedIndex])<diff.t){
      # matching Rt
      rt.ref<-TRUE
      rt.ind<-c(rt.ind,i)
    }
    # check End
    if(abs(end.times[i]-expectedEnd[expectedIndex])<diff.t){
      # matching end time
      end.ref<-TRUE
      end.ind<-c(end.ind,i)
    }
    if(sum(start.ref,rt.ref,end.ref==3)){#all times match for peak as ref
      print(paste("Start, Rt, and End times match an expected peak for Peak_Nr ",peak.nums[i],sep=""))
      # add to ref peak ind
      ref.ind<-c(ref.ind,i)
      expectedIndex<-expectedIndex+1
    }
  } # end for
  # check if number of reference peaks found matches the expected number
  numRefs<-FALSE
  foundAll<-(length(ref.ind)==expectedPeak.num)
  if(foundAll){   # return reference peak numbers and times
    numRefs<-TRUE
    print("Expected number of peaks detected at expected times.")
    refs.list[[1]]<-numRefs
  } else{
    print("Did not detected the expected number of peaks at expected times")
    refs.list[[1]]<-numRefs
  }
  ret<-cbind(vend.df$Peak_Nr[ref.ind],start.times[ref.ind],Rts[ref.ind],end.times[ref.ind])
  ret.df<-as.data.frame(ret)
  colnames(ret.df)<-c("Pk_Nr","Start","Rt","End")
  refs.list[[2]]<-ret.df
  # vendor data
  refs.list[[3]]<-vend.df[ref.ind,]
  return(refs.list)
}


# (36)
#' reference_values_no_ratio: get isotopic reference values without ratios
#' @param files vector containing character strings of .dxf file names
#' @return dataframe of reference values in the .dxf files
#' @examples
#' Usage example
#' reference_values_no_ratio(data_files)
#' @export
reference_values_no_ratio<-function(files){
  num_files<-length(files)
  msdat<-iso_read_continuous_flow(files[1:num_files])
  # reference delta values without ratio values
  delta_no_ratio<-msdat %>% iso_get_standards(file_id:reference)
  delta_no_ratio.df<-as.data.frame(delta_no_ratio)
}


# (37)
#' reference_values_ratio: get isotopic reference values including ratios
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


# (38)
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


# (39)
#' resistor_data: get resistor info for .dxf files as a dataframe
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


# (40)
#' sample_peaks_pre_process: function to separate out sample peaks (all peaks that are not reference peaks)
#' @param refTimesOutput output from reference_times_check (a list)
#' @param vend.df dataframe of vendor data
#' @return dataframe of vendor data for sample peaks
#' @examples
#' Usage Example
#' samp.df<-sample_peaks_pre_process(RTO.list,v1.df)
#' @export
sample_peaks_pre_process<-function(refTimesOutput,vend.df){
  allPeaks<-seq(1,length(vend.df$Peak_Nr))
  # start with all as samples then remove reference, flush and 1st sample peak
  samplePeaks<-allPeaks
  # get sample peaks
  refVend<-refTimesOutput[[3]]
  refPeaks<-as.numeric(refTimesOutput[[2]]$Pk_Nr)
  numR<-length(refPeaks)
  # check for no samples
  anySamples<-length(refPeaks)<length(allPeaks)
  if(anySamples){
    for(i in seq(1,numR)){
      if(refPeaks[i] %in% samplePeaks){
        matchInd<-match(refPeaks[i],samplePeaks)
        # remove from samplePeaks
        samplePeaks[matchInd]
        samplePeaks<-samplePeaks[-(matchInd)]
      }
    }
    # get sample peak vendor info
    sample.vend<-vend.df[samplePeaks,]
  } else{
    print("no sample peaks detected")
    return()#NA
  }
  return(sample.vend)
}


# (41)
#' sample_peaks_process: function to remove flush peak if present and 1st sample peak
#' @param refTimesOutput output from reference_times_check (a list)
#' @param vend.df dataframe of vendor data
#' @param flushExpT expected Rt of the flush peak
#' @param flushTint acceptable difference in Rt between flush peak Rt and expected time
#' @param firstSampleT expected Rt of the first sample peak
#' @param firstSampTint acceptable difference in Rt between first sample peak Rt and expected time
#' @return dataframe of processed sample peak vendor data
#' @examples
#' Usage Example
#' samplePeaksVend.df<-sample_peaks_process(refTC,vend1.df)
#' @export
sample_peaks_process<-function(refTimesOutput,vend.df,flushExpT=135,flushTint=10,
                               firstSampExpT=275,firstSampTint=10){
  allPeaks<-seq(1,length(vend.df$Peak_Nr))
  # start with all as samples then remove reference, flush and 1st sample peak
  samplePeaks<-allPeaks
  # after 4th reference peak
  pkNrFirstSampleAfterRefs<-as.numeric(refTimesOutput[[2]]$Pk_Nr[4])+1
  # grab sample peak vendor info using reference output
  refPeaks<-as.numeric(refTimesOutput[[2]]$Pk_Nr)
  numR<-length(refPeaks)
  # check for samples
  anySamples<-length(refPeaks)<length(allPeaks)
  if(anySamples){ # samples may have flush peak but lack other samples
    for(i in seq(1,numR)){
      if(refPeaks[i] %in% samplePeaks){
        matchInd<-match(refPeaks[i],samplePeaks)
        # remove from samplePeaks
        samplePeaks[matchInd]
        samplePeaks<-samplePeaks[-(matchInd)]
      }
    }
    # get sample peak vendor info
    sample.vend<-vend.df[samplePeaks,]
  }else{
    if(length(refTimesOutput[[3]]$Analysis)>0){
      analysisNum<-refTimesOutput[[3]]$Analysis[1]
      print(paste("no sample peaks detected for Analysis ",analysisNum,sep=""))
    }else{
      print("no sample peaks detected")
    }
    return()
  }
  if(refTimesOutput[[1]]==FALSE){
    # check if analysis number present in vend.df, return if so
    if(length(refTimesOutput[[3]]$Analysis)>0){
      analysisNum<-refTimesOutput[[3]]$Analysis[1]
      print(paste("FALSE at reference_times_check[[1]]: check reference peak times for quality for Analysis ",analysisNum,sep=""))
    }else if(length(refTimesOutput[[3]]$Date)<0){
      date<-refTimesOutput[[3]]$Date
      print(paste("FALSE at reference_times_check[[1]]: check reference peak times for quality for date ",date,sep=""))
    } else{
      print("FALSE at reference_times_check[[1]]: check reference peak times for quality")
    }
  }

  # remove flush peak and 1st sample peak
  procSample<-sample.vend
  procSampleStart<-as.numeric(procSample$Start)
  # find flush peak by time**
  for(i in seq(1,length(procSampleStart))){
    if(abs(flushExpT-procSampleStart[i])<flushTint){ #within 10 secs? or use a different way to ID flush peak?
      flushPeak<-procSample[i,]
      flushInd<-i
    }
  }
  flushNotPresent<-length(flushInd)==0
  if(flushNotPresent){
    print("no flush peak detected within expected time interval")
  } else{ # remove flush peak
    #print("flush peak detected within expected time interval")
    procSample<-procSample[-flushInd,]
  }
  # first sample
  firstSampleInd<-which(procSample$Peak_Nr==pkNrFirstSampleAfterRefs)
  firstSampleVend<-procSample[firstSampleInd,]
  firstSampleTime<-as.numeric(firstSampleVend$Rt)
  # check if other samples present (other than flush)
  if(length(firstSampleInd)==0){
    print("no sample peaks after flush detected, returning()")
    return()#NA
  }
  # check that the time for the first sample peak is in the expected interval
    sample1withinTime<-abs(firstSampExpT-firstSampleTime)<firstSampTint
    if(sample1withinTime){
      #print("First sample peak detected within expected time interval")
      procSample<-procSample[-c(firstSampleInd),]
    } else{ # dont remove first sample peak
      print("First sample peak not within expected time interval")
  }
  # remove flush peak and 1st sample peak
  #print(paste("Removed peaks at Rts ",round(as.numeric(flushPeak$Rt),3), "s (flush peak) and ",round(firstSampleTime,3),"s (1st sample peak)",sep=""))
  return(procSample)
}


# (42)
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


# (43)
#' separate_by_analysis_num: function to separate out data by analysis number
#' @param vend.df dataframe of vendor data for multiple experiments
#' @return list of vendor data separated out by analysis number
#' @examples
#' Usage Example
#' sep_analyses<-separate_by_analysis_num(CO2_1p)
#' @export
separate_by_analysis_num<-function(vend.df){
  # make lists of data separated by analysis numbers
  analysis_nums<-vend.df$Analysis
  # change "Peak.Nr" to "Peak_Nr" to use functions ***remember to change back***
  origNames<-colnames(vend.df)
  newNames<-origNames
  peakNrIndex<-which(origNames=="Peak.Nr")
  newNames<-origNames
  newNames[peakNrIndex]<-"Peak_Nr"
  colnames(vend.df)<-newNames
  # separate CO2_1p by analysis number
  totNumRows<-length(analysis_nums)
  differentAnalyses<-unique(analysis_nums)
  analysisData.list<-list()
  analysisSetInd.vec<-c()
  uniqueAnalysisInd<-1
  # check for end of analysis numbers
  for(i in seq(2,(totNumRows+1))){
    prevAnalysis<-analysis_nums[i-1]
    currAnalysis<-analysis_nums[i]
    # check for end of analysis nums
    if(i==(totNumRows+1)){
      #add last element to vec and then to list
      print("End of dataset")
      analysisSetInd.vec<-c(analysisSetInd.vec,i-1)
      analysisData.list[[uniqueAnalysisInd]]<-vend.df[analysisSetInd.vec,]
      break
    }
    # check if the analysis numbers are the same
    if(currAnalysis==prevAnalysis){
      analysisSetInd.vec<-c(analysisSetInd.vec,i-1)
    } else{ # different analysis numbers
      print("Reached end of analysis set")
      # add i-1 to vec
      analysisSetInd.vec<-c(analysisSetInd.vec,i-1)
      # add analysis set to list
      analysisData.list[[uniqueAnalysisInd]]<-vend.df[analysisSetInd.vec,]
      # reset analysis index vector and set next list index
      analysisSetInd.vec<-c()
      uniqueAnalysisInd<-uniqueAnalysisInd+1
    }
  }
  return(analysisData.list)
}


# (44)
#' sort_by_identifier_1: Function that sorts all .dxf files in a path into folders labeled with the Identifier_1 labels
#' @param path path to files to be sorted
#' @examples
#' Usage example
#' sort_by_identifier_1("~lily/Desktop/Europa MLMS/Europa_Data/Europa double comp")
#' @export
sort_by_identifier_1<-function(path){
  all_filenames<-all_filenames(path)
  un_identifiers<-unique_identifiers(path,all_filenames)
  num_identifiers<-length(un_identifiers)
  path_ID<-paste(path,"/",sep="")
  for(i in seq(1:num_identifiers)){
    # get all files with un_identifier[i]
    I1_files<-identifier_1_files(all_filenames,un_identifiers[i])
    move_identifier_1_files(I1_files,path_ID,un_identifiers[i])
  }
}


# (45)
#' standards_input_list: function to create input list for internal standards quality check
#' @param standardPaths.vec character vector containing the paths to the .dxf files for internal standards experiments
#' @param standardAnalysisNums.df dataframe containing the analysis numbers for internal standards experiments and colnames are the standard names
#' @return list of sample peak data for each analysis
#' @examples
#' Usage Example
#' standAnalysisNums.df<-read.table("standards_Analysis_Nums.txt")
#'standPaths.vec<-c(L1path,H1path,LWpath)
#' allStand_d18O.list<-standards_input_list(standPaths.vec,standAnalysisNums.df)
#' @export
standards_input_list<-function(standardPaths.vec,standardAnalysisNums.df,exp.df){ # orders of elems need to corresp in args
  allStand_d18O.list<-list()
  numStands<-length(standardAnalysisNums.df[,1])
  for(i in seq(1:numStands)){
    stand_d18O.list<-standardInputProc(standardPaths.vec[i],standardAnalysisNums.df[,i],exp.df=exp.df)
    # add to return list
    allStand_d18O.list[[i]]<-stand_d18O.list
  }
  return(allStand_d18O.list)
}


# (46)
#' stand_lm: function to perform the linear regression for the internal standards quality check and plot the results
#' @param acceptedMeas.df dataframe with accepted and measured values for internal standards (rownames are thestandard names, colnames=c(accepted,measured))
#' @return list whose first element is the lm model and the second is the model summary
#' @examples
#' Usage Example
#' avgSD_d18O<-avg_sd_d18O_standards(allStand_d18O.list)
#' standDat<-testAvg_d18O[[2]]
#' stand.lm<-stand_lm(standDat)
#' @export
stand_lm<-function(acceptedMeas.df){
  ret.list<-list()
  # linear regression
  standard.fit<-lm(accepted_d18O16O~measured_d18O16O,data=acceptedMeas.df)
  standard.fit
  fit.summ<-summary(standard.fit)
  # plot
  plot(acceptedMeas.df$measured,acceptedMeas.df$accepted,main="lm for accepted vs measured d18O/16O standards",
       xlab="measured",ylab="accepted",pch=20,col="blue")
  abline(standard.fit,col="light blue")
  ret.list[[1]]<-standard.fit
  ret.list[[2]]<-fit.summ
  return(ret.list)
}


# (47)
#' standardInputProc: function that gets sample peak and d18O/16O data for internal standards quality check
#' @param path path to the .dxf internal standards experiment files
#' @param standardAnalysisNums.vec numeric vector of analysis numbers
#' @return list of dataframes containing d18O/16O vendor data for each file
#' @examples
#' Usage Example
#' stand_d18O.list<-standardInputProc(standPath,standAnalysisNum)
#' @export
standardInputProc<-function(path,standardAnalysisNums.vec,exp.df){
  standardFiles<-all_filenames(path)
  setwd(path)
  iso<-iso_read_continuous_flow(standardFiles)
  standFileNames<-find_files_by_analysis_num(iso,standardAnalysisNums.vec)
  # get sample peak d18O data from standard files
  d18O.list<-d18O_samples_list(standFileNames,exp.df=exp.df)
  return(d18O.list)
}


# (48)
#' ThresholdingAlgo: function based on dispersion to peak peak start times in time series data
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
ThresholdingAlgo<-function(y,lag,threshold,influence){
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


# (49)
#' trap_area_allPks: function to calculate area of all peaks in an experiment via trapz (wrapper function for all_PA_trap())
#' @param raw.df dataframe of raw data
#' @param vend.df dataframe of vendor data
#' @param mV.rawName character string of the name of intensity column to use (in raw.df)
#' @return dataframe containing peak numbers and areas
#' @examples
#' Usage Example
#' pkAreas.df<-trap_area_allPks(raw.df,vi.df,"v44.mV")
#' @export
trap_area_allPks<-function(raw.df,vend.df,mV.rawName){
  # get mass voltage index
  massInd<-which(colnames(raw.df)==mV.rawName)
  massV<-raw.df[,massInd]
  massT<-raw.df$time.s
  massSt<-as.numeric(vend.df$Start)
  massEnd<-as.numeric(vend.df$End)
  massPkNr<-as.numeric(vend.df$Peak_Nr)

  rawMass.mat<-matrix(c(massT,massV),ncol=2)
  rawMass.df<-as.data.frame(rawMass.mat)
  colnames(rawMass.df)<-c("time.s",mV.rawName)

  vendMass.mat<-matrix(c(massSt,massEnd,massPkNr),ncol=3)
  vendMass.df<-as.data.frame(vendMass.mat)
  colnames(vendMass.df)<-c("Start","End","Pk_Nr")

  ret.list<-list()
  ret.list[[1]]<-rawMass.df
  ret.list[[2]]<-vendMass.df

  rawIntTime<-ret.list[[1]]
  vendTimePk<-ret.list[[2]]

  start<-vendTimePk$Start
  end<-vendTimePk$End
  PkNr<-vendTimePk$Pk_Nr
  time<-rawIntTime$time.s
  mV<-rawIntTime[,2]

  # call area func
  allA<-all_PA_trap(start,end,time,mV,PkNr)
  return(allA)
}


# (50)
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


# (51)
#' vendor_info: get selected vendor info for a .dxf file:
#' (Identifier 1, Nr., Start, Rt, End, Intensity All, rIntensity All, Ampl_44, Ampl_45, Ampl_46, d13C/12C, d18O/16O)
#' @param file file name
#' @return dataframe of vendor information
#' @examples
#' Usage Example
#' vendor_info(file)
#' @export
vendor_info<-function(file){
    msdat<-iso_read_continuous_flow(file)
    file.info<-msdat %>% iso_get_file_info()
    ident1<-file.info$`Identifier 1`
    vendor_info<-msdat %>% iso_get_vendor_data_table()
    peak_num<-vendor_info$Nr.

    vendor_info_select<-cbind(ident1,peak_num,
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
    colnames(vendor_info_select.df)<-c("Identifier_1","Peak_Nr","Start","Rt","End",
                                     "Intensity_All","rIntensity_All",
                                     "Ampl_44","Ampl_45","Ampl_46",
                                     "d13C/12C","d18O/16O")
  return(vendor_info_select.df)
}


# (52)
#' vendor_info_all: function to get vendor table data from several .dxf files
#' @param files vector of filenames
#' @return list that contains vendor table info for each file
#' @examples
#' Usage example
#' vendor_info_all(files.vec)
#' @export
vendor_info_all<-function(files){
  vi.list<-list()
  for(i in seq(1,length(files))){
    vi.df<-vendor_info(files[i])
    vi.list[[i]]<-vi.df
  }
  return(vi.list)
}

