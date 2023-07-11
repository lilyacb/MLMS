library(isoreader)
source("/home/lily/Documents/git/MLMS/R/MLMSfunctions.R")
source("/home/lily/Documents/git/MLMS/R/QCfunctions.R")



######### IRMS automated QC based on Isodat output

# qc.config.list: default params for CO2 IRMS data but can be changed by the user
qc.config.list<-list()
# checks: checks and tasks to perform - dataframe of booleans
qc.checks.df<-as.data.frame(matrix(c(T,T,T,T,T,T,T,T,T),nrow=1))
colnames(qc.checks.df)<-c("checkMaxPkNum","checkRefTimes",
                         "checkRefPkNum","checkSDrefIso",
                         "checkSDsampIso","checkRefPkIntensity",
                         "checkSampPkIntensity","doInternalStandardsStats",
                         "sortOutBadData")
qc.config.list[[1]]<-qc.checks.df

# testValues: maxPkNum, expected number of ref peaks, relDiffs in intensity thresholds
qc.pkNums.df<-as.data.frame(matrix(c(5, 18, 0.1, 0.1),nrow=1))
colnames(qc.pkNums.df)<-c("expRefPkNr","maxPkNum","relDiffIntRefs",
                          "relDiffIntSamps") # TODO: make sep df for inensity check
# add to qc.config.list
qc.config.list[[2]]<-qc.pkNums.df

# isoSDvalues
# reference peaks
refIsoCheckSD.df<-as.data.frame(matrix(c(0.1, 0.1),nrow=1))
colnames(refIsoCheckSD.df)<-c("d18O16O", "d13C12C")
refIsoCheckSD.df

# sample peaks
sampIsoCheckSD.df<-as.data.frame(matrix(c(0.2,0.2),nrow=1))
colnames(sampIsoCheckSD.df)<-c("d18O16O", "d13C12C")
sampIsoCheckSD.df

# add both to qc.config
qc.config.list[[3]]<-refIsoCheckSD.df
qc.config.list[[4]]<-sampIsoCheckSD.df

# refPkExpTimes.df
# get expected times for reference peaks
expRef.df<-read.table("~/Desktop/EuropaMLMS/Europa_Data/exampleDXFdata/SupplementalFiles/referencePeaks_expectedTimes.txt")

# add to qc.config.list
qc.config.list[[5]]<-expRef.df
qc.config.list

# standardsInfo.df
# user-specified:
standardsLab.vec<-c("H1","L1","LW")
# create mapping for user labels
standInfo.df<-as.data.frame(matrix(rep(NA,2*length(standardsLab.vec)),ncol=2))
colnames(standInfo.df)<-c("standSerial","standLab")
# build df
standInfo.df$standLab<-standardsLab.vec
standSerial.vec<-rep(NA,length(standInfo.df$standLab))
for(i in seq(1,length(standSerial.vec))){
  standSerial.vec[i]<-paste("S",i,sep="")
}
standInfo.df$standSerial<-standSerial.vec
# add to qc.config.list
qc.config.list[[6]]<-standInfo.df

# add names to qc.config.list
names(qc.config.list)<-c("checks", "testValues","refIsoSD",
                         "sampIsoSD", "refPkExpTimes","standInfo")
qc.config.list
# $checks
#  checkMaxPkNum checkRefTimes checkRefPkNum checkSDrefIso checkSDsampIso
#1          TRUE          TRUE          TRUE          TRUE           TRUE
#  checkRefPkIntensity checkSampPkIntensity doInternalStandardsStats sortOutBadData
#1                TRUE                 TRUE                     TRUE           TRUE
#
# $testValues
#   expRefPkNr maxPkNum relDiffIntRefs relDiffIntSamps
#1          5       18            0.1             0.1
#
# $refIsoSD
#  d18O16O d13C12C
#1     0.1     0.1
#
# $sampIsoSD
#  d18O16O d13C12C
#1     0.6     0.6
#
# $refPkExpTimes
#. Ref_Peak_Nr Expected_Start Expected_Rt Expected_End
#1           1        27.1700     47.3888      50.4929
#2           2        67.0193     87.1840      90.3422
#3           4       166.5653    186.7608     189.8572
#4           5       206.4920    226.5096     229.6987
#5          16       843.3227    863.5029     866.6069
#
# $standInfo
#  standSerial standLab
#1          S1       H1
#2          S2       L1
#3          S3       LW



### wrapQC: wrapper function for QC and internal standard analysis

# NOTE: copy orig directories to demo sorting process

# initialize QC
# europa bacteria
#data.path<-"/media/lily/My Book/Lilys_NASA_microbe_data/Europa_bacteria1/Europa_bacteria/"
#setwd(data.path)
# microbial mud
#data.path<-"/media/lily/My Book/Lilys_NASA_microbe_data/Europa_microbial_mud1/Europa_microbial_mud/"
# europa_bacteria
data.path<-"/media/lily/My Book/Lilys_NASA_microbe_data/Europa_bacteria1/Europa_bacteria/"
setwd(data.path)
# abiotic

# specify dataset label for processing
#datasetLabel<-"microbial_mud"
datasetLabel<-"europa_bacteria"
# datasetLabel<-"abiotic"

# read in data and prep for QC
dxfFiles<-all_dxf_files()


# get file info
allFilesInfo.list<-file_info_all(dxfFiles)
head(allFilesInfo.list)[1:3]
# microbial mud
# [[1]]
#file_id Identifier1 Analysis Preparation       Date_and_Time
#1 170418_H1_.dxf          H1     2959   CO2 in He 2017-04-18 13:45:36
#
#[[2]]
#file_id Identifier1 Analysis Preparation       Date_and_Time
#1 170418_H1_(1).dxf          H1     2969   CO2 in He 2017-04-18 16:35:06
#
#[[3]]
#file_id Identifier1 Analysis Preparation       Date_and_Time
#1 170418_L1_.dxf          L1     2958   CO2 in He 2017-04-18 13:28:39

# change default values in qc.config.list for microbes
qc.config.list$sampIsoSD$d18O16O<-0.6
qc.config.list$sampIsoSD$d13C12C<-0.6

# for europa_bacteria, change standInfo.df in qc.config.list
qc.config.list$standInfo
# TODO: standLabEuropa<-c() *** check with Bethany about "acceptedVals"


# give output directory for QC'd data and QC report
#qcOut.path<-"/home/lily/Desktop/EuropaMLMS/microbe/QC_SU22/bio_abio/microbial_mud/"
qcOut.path<-"/home/lily/Desktop/EuropaMLMS/microbe/QC_SU22/bio_abio/europa_bacteria/"
#qcOut.path<-"/home/lily/Desktop/EuropaMLMS/microbe/QC_SU22/bio_abio/abiotic/"

# wrapQC<-function(qc.config=qc.config.list,data.path,datasetLabel="IRMSdataset",qcOut.path){
#c("microbial_mud","europa_data","abiotic","TitanAnalogs") - loop through different batches
  # TODO: input for user-specified data columns
  dataCols.vec<-c( "fileId","Identifier1","Analysis","Preparation","DateTime",
                    "PeakNr","Start","Rt","End","Ampl44","Ampl45",
                    "Ampl46","BGD44","BGD45","BGD46","rIntensity44","rIntensity45",
                    "rIntensity46","rIntensityAll","Intensity44","Intensity45",
                    "Intensity46","IntensityAll","ListFirstPeak","rR45CO244CO2",
                    "rR46CO244CO2","IsRef","R45CO244CO2","RefName","rd45CO244CO2",
                    "d45CO244CO2", "R46CO244CO2", "rd46CO244CO2","d46CO244CO2",
                    "R13C12C","d13C12C","AT13C12C","R18O16O","d18O16O", "AT18O16O",
                    "R17O16O","d17O16O","Rps45CO244CO2","Rps46CO244CO2")

  # make combined vendor and file info df for processing
  combList<-combineVendFileInfo(data.path,dataCols.vec,datasetLabel)
  # original num experiments
  # microbial_mud: 86
  # europa_bacteria: 216

  # get params from qc.config.list
  # checks boolean
  checks.df<-qc.config.list$checks
  # testValues
  testVals.df<-qc.config.list$testValues
  # ref pk iso ratios
  refPkIso.df<-qc.config.list$refIsoSD
  # samp pk iso ratios
  sampPkIso.df<-qc.config.list$sampIsoSD
  # reference peak expected times
  expRefT.df<-qc.config.list$refPkExpTimes
  # standards names
  standInfo.df<-qc.config.list$standInfo

  # run QC
  # TODO: *generalize args further so dfs are the args instead
  # so that the names of the iso ratios dont matter -- pass the dfs in qc.config.list
  currProc<-removeFailedAnalysesDXF(combList, expRef.df = expRefT.df,
                                    maxPkNum = testVals.df$maxPkNum,
                                    expRefPkNr = testVals.df$expRefPkNr,
                                    sdCrefIso.thresh = refPkIso.df$d13C12C, # *pass df
                                    sdOrefIso.thresh = refPkIso.df$d18O16O, # * pass df
                                    relDiffInt.thresh = testVals.df$relDiffIntRefs, #*sep check for samps, optional
                                    sdCsampIso.thresh = sampPkIso.df$d13C12C, #*pass df
                                    sdOsampIso.thresh = sampPkIso.df$d18O16O) #*pass df

  # adjust reference data for any samples removed
  # list of 2 lists: one for ref peak dfs, one for samps
  currFiltered<-rmRefDatDXF(currProc)
  # want same number of analyses in both ref and samp sets
  lengthEqual<-length(currFiltered[[1]])==length(currFiltered[[2]])

  # reference and sample peak intensity similarity - dont run for microbes
  # TODO: make this process into QC function
  # checkPkIntRelDiff()
  # function:
  # # reference and sample peak intensity similarity - dont run for microbes
  #if(lengthEqual){
  #  # separate dfs into lists
  #  refFiltered.df<-currFiltered[[1]]
  #  sampFiltered.df<-currFiltered[[2]]
  # separate by experiment
  # refFiltered.list<-separate_by_analysis_numDXF(refFiltered.df)
  #  sampFiltered.list<-separate_by_analysis_numDXF(sampFiltered.df)
  # check peak similarities
  # rmAnalysis.vec<-c()
  #  rmInd.vec<-c()
  #  for(k in seq(1,length(refFiltered.list))){
  #    # make vecs of intensities
  #   analysis<-refFiltered.list[[k]]$Analysis[1]
  #    ID1<-refFiltered.list[[k]]$Identifier1[1]
  #    refInt.vec<-refFiltered.list[[k]]$Ampl44
  #    refPkNr.vec<-refFiltered.list[[k]]$PeakNr
  #    sampInt.vec<-sampFiltered.list[[k]]$Ampl44
  #    sampPkNr.vec<-sampFiltered.list[[k]]$PeakNr
  #    # put data together
  #    refSampPkNr.vec<-sort(c(refPkNr.vec,sampPkNr.vec))
  #    refSampVendAmpl.vec<-c()
  #    refInd=0
  #    sampInd=0
  #    for(j in seq(1,length(refSampPkNr.vec))){
  #      # match peak numbers with intensity amplitudes
  #     if(refSampPkNr.vec[j] %in% refPkNr.vec){
  #    refInd<-refInd+1
  #     refSampVendAmpl.vec<-c(refSampVendAmpl.vec,refInt.vec[refInd])
  # }else if(refSampPkNr.vec[j] %in% sampPkNr.vec){
  #    sampInd<-sampInd+1
  #    refSampVendAmpl.vec<-c(refSampVendAmpl.vec,sampInt.vec[sampInd])
  # }
  #  }
  # combine in df
  # refSampAmpl.mat<-matrix(c(refSampPkNr.vec,refSampVendAmpl.vec),ncol=2)
  #  refSampAmpl.df<-as.data.frame(refSampAmpl.mat)
  #  colnames(refSampAmpl.df)<-c("PeakNr","Ampl44")
  # peak intensity similarity check using intensity_similarityDXF
  # refSampIntCheck<-intensity_similarityDXF(vendAmpl=refSampAmpl.df$Ampl44,amplName="Ampl44",peakNr.vec=refSampAmpl.df$PeakNr,relDiffInt.thresh=0.75)
  #  maxRefSampRelDiff<-refSampIntCheck[[3]]
  #  if(!refSampIntCheck[[1]]){
  #   rmInd=k
  #  print(paste("Significant difference in intensities of reference and sample peaks in Analysis ",analysis , " detected.",sep=""))
  #  # put that analysis number in vec to be removed
  # rmAnalysis.vec<-c(rmAnalysis.vec,analysis)
  #  rmInd.vec<-c(rmInd.vec,rmInd)
  #  failedRefSampIntAnalysis<-c(failedRefSampIntAnalysis,analysis)
  #  failedRefSampIntID1<-c(failedRefSampIntID1,ID1)
  #failedRefSampRelDiff.vec<-c(failedRefSampRelDiff.vec,maxRefSampRelDiff)
  # }
  #}
  #}
  # update currFiltered to write new results to file
  #if(length(rmAnalysis.vec)>0){
  # remove those analyses from the lists
  #print("removing analyses...")
  #refFiltered.list[rmInd.vec]<-NULL
  #sampFiltered.list[rmInd.vec]<-NULL
  # change lists back to dfs
  # if(length(refFiltered.list)!=0){
  #ret.list<-threadRefSampListsToDF(refFiltered.list,sampFiltered.list)
  ## currFiltered[[1]]<-ret.list[[1]]
  # currFiltered[[2]]<-ret.list[[2]]
  #}else{ # all data processed out!
  # currFiltered[[1]]<-as.data.frame(matrix(NA))
  # currFiltered[[2]]<-as.data.frame(matrix(NA))
  #}
  #}

  # write to file in qc output dir
  setwd(qcOut.path)
  # reference data
  #head(currFiltered[[1]])
  refFileName<-paste(datasetLabel,"_refs_QC_",Sys.Date(),".csv",sep="")
  # sample peak data
  sampFileName<-paste(datasetLabel,"_samps_QC_",Sys.Date(),".csv",sep="")
  # get ref and sample peak data
  refsFiltered.df<-currFiltered[[1]]
  sampsFiltered.df<-currFiltered[[2]]
  # write both QC'd datasets to file
  write.table(refsFiltered.df,file=refFileName,quote=F,row.names=F,sep=",")
  # sample data
  write.table(sampsFiltered.df,file=sampFileName,quote=F,row.names=F,sep=",")
  # write QC summary
  writeFailStats(currProc,datasetLabel)

  # internal standards statistics
  #unique(sampsFiltered.df$Identifier1)
  # vec of internal standard IDs
  standLabs<-standInfo.df$standLab

  # pick out experiements that are internal standards and get d18O/16O and d13C/12C data for sample peaks
  standIsoR.list<-DXFvendListFromID1(sampsFiltered.df,standName.vec=standLabs)
  #standIsoR.list
  standAvgSD<-avgSD_d18O_standards(standIsoR.list,standNames=standLabs,
                                   #TODO: add standAcceptedVals and accStandRatioSD as args to QC wrapper
                                   standAcceptedVals.vec=c(4.85,-8.55,-3.85),accStandRatioSD=c(0.2,0.2,0.2))
  # write data to file
  # summary of avg d18O/16O and SD d18O/16O for all internal standards
  # TODO: generalize with iso ratios
  intStand1.fileName<-paste(datasetLabel,"_intStand_d180160Summ_",".txt",sep="")
  standList<-standAvgSD[[1]]
  restAvgSD<-c(2:length(standList))
  write.table(data.frame(standList[[1]]),intStand1.fileName,col.names=T,quote=F,row.names=F)
  lapply(standList[restAvgSD],function(x) write.table(data.frame(x),intStand1.fileName,append=T,
                                                      col.names=F,quote=F,row.names=F))
  # df of accepted and measured standard isotope ratios
  stand18O16O<-standAvgSD[[2]]
  intStand2.fileName<-paste(datasetLabel,"_intStand_AccMeas_d18O16O.txt",sep="")
  write.table(stand18O16O,intStand2.fileName,row.names=F,quote=F)

  # df of sds of iso ratios
  standSD<-standAvgSD[[3]]
  intStand3.fileName<-paste(datasetLabel,"_intStand_SDd18O16O.txt")
  write.table(standSD,intStand3.fileName,row.names=T,quote=F)

  # df of values and sds for iso ratios
  standFullSumm<-standAvgSD[[4]]
  intStand4.fileName<-paste(datasetLabel,"_intStand_fullSumm.txt")
  write.table(standFullSumm,intStand4.fileName,row.names=T,quote=F)

  ## save graph to file
  pdf(file=paste(datasetLabel,"_int_stand_d18O_lm.pdf",sep=""),width=6,height=4)
  standlm<-stand_lm(stand18O16O)
  dev.off()

  ## write stats
  # slope, intercept, r^2
  lmr2<-standlm[[2]]$r.squared
  lmresid<-standlm[[2]]$residuals
  lmCoeff.df<-as.data.frame(matrix(c(standlm[[1]]$coefficients,lmr2,lmresid),ncol=6))
  colnames(lmCoeff.df)<-c("intercept","slope","r^2","L1_resid","H1_resid","LW_resid")
  # write lm stats to file
  write.table(lmCoeff.df,file=paste(datasetLabel,"_intStand_lm_d18O16O.txt",sep=""),row.names=F,quote=F)
  # std error, t val and p val
  lmCoeffStat<-standlm[[2]]$coefficients
  write.table(lmCoeffStat,file=paste(datasetLabel,"instStand_coefficients.txt",sep=""),row.names=T,quote=F)

  # graph passed chromatograms
  passedFiles<-unique(sampsFiltered.df$fileId)
  setwd(data.path)
  rawList<-raw_data_all(passedFiles)
  #length(rawList)
  # microbial mud: 58
  # europa_bacteria: 118

  # analysis.vec,ID1.vec,date.vec - for chromatograms of passed experiments
  passedFileInfo.df<-file_info(passedFiles)
  passedAnalyses.vec<-passedFileInfo.df$Analysis
  passedID1.vec<-passedFileInfo.df$Identifier1
  passedDate.vec<-as.character(passedFileInfo.df$Date_and_Time)
  newDate.vec<-c()
  for(i in seq(1, length(passedDate.vec))){
    currDate<-substr(passedDate.vec[i],1,10)
    # update with yr-mo-day
    newDate.vec[i]<-currDate
  }
  passedDate.vec<-newDate.vec
  # plot passed chromatograms
  setwd(qcOut.path)
  generic_plot_all_raw(rawList,fileName=paste(datasetLabel,"_passedQC_chromatograms.pdf",sep=""),
                        analysis.vec=passedAnalyses.vec, ID1.vec=passedID1.vec,date.vec=passedDate.vec)

  # graph failed chromatograms
  # read in failed analysis nums then get file names
  # dxfFiles = all dxf files in directory before QC
  passedInd<-which(dxfFiles %in% passedFiles)
  failedFiles<-dxfFiles[-passedInd]
  #analysis.vec,ID1.vec,date.vec - for chromatogram plots
  setwd(data.path)
  failedFileInfo.df<-file_info(failedFiles)
  failedAnalyses.vec<-failedFileInfo.df$Analysis
  failedID1.vec<-failedFileInfo.df$Identifier1
  failedDate.vec<-as.character(failedFileInfo.df$Date_and_Time)
  newDate.vec<-c()
  for(i in seq(1, length(failedDate.vec))){
    currDate<-substr(failedDate.vec[i],1,10)
    # update with yr-mo-day
    newDate.vec[i]<-currDate
  }
  failedDate.vec<-newDate.vec
  # get raw data
  rawList<-raw_data_all(failedFiles)

  # plot all raw data
  setwd(qcOut.path)
  generic_plot_all_raw(rawList,fileName=paste(datasetLabel,"_failedQC_chromatograms.pdf",sep=""),
                       analysis.vec=failedAnalyses.vec, ID1.vec=failedID1.vec,date.vec=failedDate.vec)

  # TODO: add code to generate plotly graphs for passed/failed chromatograms
  # TODO: optional sorting of passed/failed experiments into different folders


#} # end wrapper func


