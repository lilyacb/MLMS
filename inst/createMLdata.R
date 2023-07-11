#### create ML data

## check out previously-processed abiotic data
setwd("~/Desktop/EuropaMLMS/microbe/QC_SU22/bio_abio/abiotic")
abio.dat<-read.csv("dxfDataLarge_samps_allChecks.csv")
head(abio.dat)
dim(abio.dat)
# [1] 2529   44
unique(abio.dat$Identifier1)
# process Preparation labels?
# replace spaces in Preparation with underscores
abio.dat<-replace_space_prep(abio.dat)
unique(abio.dat$Preparation)
# [1] "He_only"                 "CO2_in_He"               "2%_CO2_in_He_48hrs"
# [4] "2%_CO2_in_He_24hrs"      "0.3%_CO2_in_He_72_hr"    "1%_CO2_in_He_72hrs"
# [7] "0.115ml_CO2_in_He_24hrs" "0.115ml_CO2_in_He_48hrs"
# 1% = 0.115 mL
# process Identifier1 labels
abio.dat<-replace_space_Identifier1(abio.dat)
unique(abio.dat$Identifier1)
#[1] "H1"                  "L1"                  "LW"
#[4] "MgSO4_L_+_NaCl_L"    "MgSO4_L_+_NaCl_U"    "MgSO4_L_+_NaHCO3_L"
#[7] "MgSO4_L"             "Na2SO4_L_+_NaCl_L"   "Na2SO4_L_+_NaCl_U"
#[10] "Na2SO4_L_+_NaHCO3_L" "Na2SO4_L_+_NaHCO3_U" "Na2SO4_L"
#[13] "Na2SO4_U"            "NaCl_L"              "NaCl_U"
#[16] "NaHCO3_L_+_NaCl_L"   "NaHCO3_L_+_NaCl_U"   "NaHCO3_L"
#[19] "NaHCO3_U_+_NaCl_L"   "NaHCO3_U_+_NaCl_U"   "NaHCO3_U"


# **** add areas of peaks
# TODO: make into function for data processing
# dir to abiotic dxf data
setwd("~/Desktop/EuropaMLMS/Europa_Data/exampleDXFdata/CopyOfUnsortedData")
# get raw data for all files
#abioDxfFiles<-all_dxf_files()

# get only the ones in abio.dat
length(unique(abio.dat$Analysis))
# 281
# get these analyses that passed QC
abioPassedFiles<-unique(abio.dat$fileId)

# get raw intensity vs time data
abioRawList<-raw_data_all(abioPassedFiles)
head(abioRawList[[1]])
#         file_id tp time.s   v44.mV    v45.mV   v46.mV
#1 170404_H1_.dxf  1  0.209 1.088951 0.3888736 1.906578
#2 170404_H1_.dxf  2  0.418 1.083217 0.4060438 1.893181
#3 170404_H1_.dxf  3  0.627 1.090862 0.3907814 1.939115
#4 170404_H1_.dxf  4  0.836 1.087039 0.4060438 1.852996
#5 170404_H1_.dxf  5  1.045 1.096596 0.4556530 1.898922
#6 170404_H1_.dxf  6  1.254 1.079395 0.4041359 1.898922

# get file info for passed files
abioFileInfo.df<-file_info(abioPassedFiles)
head(abioFileInfo.df)
# file_id Identifier1 Analysis        Preparation       Date_and_Time
#1           170404_H1_.dxf          H1     2729            He-Only 2017-04-04 13:51:15
#2        170404_H1_(1).dxf          H1     2745            He-Only 2017-04-04 18:22:28
#3        170404_H1_(2).dxf          H1     2763            He-Only 2017-04-04 23:27:36
#4 170503_H1__CO2_in_He.dxf          H1     3248          CO2 in He 2017-05-03 20:34:24
#5 170504_H1__CO2_in_He.dxf          H1     3275          CO2 in He 2017-05-04 04:11:56
#6        170518_H1_(1).dxf          H1     3832 2% CO2 in He 48hrs 2017-05-18 23:51:23

# get vendor data
abioVendInfo.list<-vendor_info_all(abioPassedFiles)
head(abioVendInfo.list[[1]])

# input for pk area - label for intensity to use
rawVname<-"v46.mV"
library(pracma)
# loop and calc peak areas in list then add to df
abioPkAreas.df<-trap_area_allPks(raw.df=abioRawList[[1]],vend.df=abioVendInfo.list[[1]],
                              mV.rawName = rawVname)
abioPkAreas.df
# PkNr  trapArea
#1     1 60.693675
#2     2 60.780505
#3     3  6.572158
#4     4 60.750938
#5     5 60.530529
#6     6 28.689129
#7     7 27.308097
#8     8 25.956695
#9     9 24.657210
#10   10 23.377436
#11   11 22.222060
#12   12 21.111828
#13   13 20.051748
#14   14 19.006667
#15   15 18.069748
#16   16 60.632867

# will need to match areas to pkNr in data for sample peaks

# file name is the name of the data frame
abioPkArea.list<-list()
VIfiles.vec<-c()
for(i in seq(1,length(abioRawList))){
  abioPkArea<-trap_area_allPks(raw.df=abioRawList[[i]],vend.df=abioVendInfo.list[[i]],
                                   mV.rawName = rawVname)
  abioPkArea.list[[i]]<-abioPkArea
  fileId<-abioVendInfo.list[[i]]$file_id[1]
  #names(abioPkArea.list[[i]])<-fileId
  VIfiles.vec<-c(VIfiles.vec,fileId)
}
length(abioPkArea.list)
head(abioPkArea.list[[1]])
names(abioPkArea.list)<-VIfiles.vec
head(names(abioPkArea.list))

#abioAnalyses.vec<-abioFileInfo.df$Analysis
FIfiles<-abioFileInfo.df$file_id
head(FIfiles)
VIfiles<-names(abioPkArea.list)
head(VIfiles)
# same order!
head(unique(abio.dat$fileId))

# add additoinal values to data
# separate by analysis then add data for each experiment

abioSep.list<-separate_by_analysis_numDXF(abio.dat)
for(i in seq(1,length(abioSep.list))){
      currDf<-abioSep.list[[i]]
      # add peak areas -- add avg peak areas during tsfeature extraction
      # get sample peak numbers
      sampPkNrs<-abioSep.list[[i]]$PeakNr
      currDf$pkArea<-abioPkArea.list[[i]]$trapArea[sampPkNrs]
      # update abioSep.list[[i]]
      abioSep.list[[i]]$pkArea<-currDf$pkArea
      # add iso ration d18O/13C
      d18O13C.lab<-d18O13Clab(currDf)$d18O13C
      abioSep.list[[i]]$d18O13C<-d18O13C.lab
      # add biotic label
      abioLab<-rep("abiotic",dim(currDf)[1])
      abioSep.list[[i]]$biotic<-abioLab
      # add avgs of deltas
      deltaAvg.df<-avg3Deltas(currDf)
      d1<-rep(deltaAvg.df$avg_d13C12C,dim(currDf)[1])
      d2<-rep(deltaAvg.df$avg_d18O13C,dim(currDf)[1])
      d3<-rep(deltaAvg.df$avg_d18O16O,dim(currDf)[1])
      abioSep.list[[i]]$avg_d13C12C<-d1
      abioSep.list[[i]]$avg_d18O13C<-d2
      abioSep.list[[i]]$avg_d18O16O<-d3
      # add intensity avg values
      avg_Ampl46<-rep(mean(abioSep.list[[i]]$Ampl46),dim(currDf)[1])
      abioSep.list[[i]]$avg_Ampl46<-avg_Ampl46
      # add intensityAll avg values
      avg_IntensityAll<-rep(mean(abioSep.list[[i]]$IntensityAll),dim(currDf)[1])
      abioSep.list[[i]]$avg_IntensityAll<-avg_IntensityAll
      # add rIntensityAll avg values
      avg_rIntensityAll<-rep(mean(abioSep.list[[i]]$rIntensityAll),dim(currDf)[1])
      abioSep.list[[i]]$avg_rIntensityAll<-avg_rIntensityAll
}
head(abioSep.list[[1]])








## microbial mud
setwd("~/Desktop/EuropaMLMS/microbe/QC_SU22/bio_abio/microbial_mud")
mmbio.dat<-read.csv("microbial_mud_samps_QC_2022-06-23.csv",header=T)
head(mmbio.dat)
dim(mmbio.dat)
# [1] 522  44

# process Preparation labels
# replace spaces in Preparation with underscores
mmbio.dat<-replace_space_prep(mmbio.dat)
unique(mmbio.dat$Preparation)
# [1] "CO2_in_He" "He_only"

# process Identifier1 labels
unique(mmbio.dat$Identifier1)
#[1] "H1"                 "L1"                 "LW"                 "MgSO4 L"
#[5] "MgSO4 L+NaCl L"     "MgSO4 L+NaCl U"     "MgSO4 L+NaHCO3 L"   "MgSO4 L+NaHCO3 U"
#[9] "NACL L"             "NACL U"             "NaHCO3 L"           "NaHCO3 L+NaCl L"
#[13] "NaHCO3 L+NaCl U"    "NaHCO3 U"           "NaHCO3 U+NaCl L"    "NaHCO3 U+NaCl U"
#[17] "MgSO4 L B"          "MgSO4 L+NaCl L B"   "MgSO4 L+NaHCO3 L B" "MgSO4 L+NaHCO3 U B"
#[21] "NACL U B"           "NaHCO3 L+NaCl L B " "NaHCO3 L+NaCl U B"

# fix Identifer1s
mmID1.map<-read.csv("microbialMudIdentifier1map.csv",header=T)
head(mmID1.map)
# oldIdentifier1   newIdentifier1
#1             H1               H1
#2             L1               L1
#3             LW               LW
#4        MgSO4 L          MgSO4_L
#5 MgSO4 L+NaCl L MgSO4_L_+_NaCl_L
#6 MgSO4 L+NaCl U MgSO4_L_+_NaCl_U


newMM.df<-mmbio.dat
newID1.vec<-c()
for(i in seq(1,length(mmbio.dat$Identifier1))){
  currID1<-mmbio.dat$Identifier1[i]
  # match to mapping
  id1Ind<-which(mmID1.map$oldIdentifier1==currID1)
  newID1<-mmID1.map$newIdentifier1[id1Ind]
  newID1.vec<-c(newID1.vec,newID1)
}
length(newID1.vec)
newMM.df$Identifier1<-newID1.vec
unique(newMM.df$Identifier1)
#[1] "H1"                 "L1"                 "LW"                 "MgSO4_L"
#[5] "MgSO4_L_+_NaCl_L"   "MgSO4_L_+_NaCl_U"   "MgSO4_L_+_NaHCO3_L" "MgSO4_L_+_NaHCO3_U"
#[9] "NaCl_L"             "NaCl_U"             "NaHCO3_L"           "NaHCO3_L_+_NaCl_L"
#[13] "NaHCO3_L_+_NaCl_U"  "NaHCO3_U"           "NaHCO3_U_+_NaCl_L"  "NaHCO3_U_+_NaCl_U"

mmbio.dat<-newMM.df



### peak areas and data building
##### get raw intensity vs time data
abioRawList<-raw_data_all(abioPassedFiles)
head(abioRawList[[1]])
#         file_id tp time.s   v44.mV    v45.mV   v46.mV
#1 170404_H1_.dxf  1  0.209 1.088951 0.3888736 1.906578
#2 170404_H1_.dxf  2  0.418 1.083217 0.4060438 1.893181
#3 170404_H1_.dxf  3  0.627 1.090862 0.3907814 1.939115
#4 170404_H1_.dxf  4  0.836 1.087039 0.4060438 1.852996
#5 170404_H1_.dxf  5  1.045 1.096596 0.4556530 1.898922
#6 170404_H1_.dxf  6  1.254 1.079395 0.4041359 1.898922

# get file info for passed files
abioFileInfo.df<-file_info(abioPassedFiles)
head(abioFileInfo.df)
# file_id Identifier1 Analysis        Preparation       Date_and_Time
#1           170404_H1_.dxf          H1     2729            He-Only 2017-04-04 13:51:15
#2        170404_H1_(1).dxf          H1     2745            He-Only 2017-04-04 18:22:28
#3        170404_H1_(2).dxf          H1     2763            He-Only 2017-04-04 23:27:36
#4 170503_H1__CO2_in_He.dxf          H1     3248          CO2 in He 2017-05-03 20:34:24
#5 170504_H1__CO2_in_He.dxf          H1     3275          CO2 in He 2017-05-04 04:11:56
#6        170518_H1_(1).dxf          H1     3832 2% CO2 in He 48hrs 2017-05-18 23:51:23

# get vendor data
abioVendInfo.list<-vendor_info_all(abioPassedFiles)
head(abioVendInfo.list[[1]])

# input for pk area - label for intensity to use
rawVname<-"v46.mV"
library(pracma)
# loop and calc peak areas in list then add to df
abioPkAreas.df<-trap_area_allPks(raw.df=abioRawList[[1]],vend.df=abioVendInfo.list[[1]],
                                 mV.rawName = rawVname)
abioPkAreas.df
# PkNr  trapArea
#1     1 60.693675
#2     2 60.780505
#3     3  6.572158
#4     4 60.750938
#5     5 60.530529
#6     6 28.689129
#7     7 27.308097
#8     8 25.956695
#9     9 24.657210
#10   10 23.377436
#11   11 22.222060
#12   12 21.111828
#13   13 20.051748
#14   14 19.006667
#15   15 18.069748
#16   16 60.632867

# will need to match areas to pkNr in data for sample peaks

# file name is the name of the data frame
abioPkArea.list<-list()
VIfiles.vec<-c()
for(i in seq(1,length(abioRawList))){
  abioPkArea<-trap_area_allPks(raw.df=abioRawList[[i]],vend.df=abioVendInfo.list[[i]],
                               mV.rawName = rawVname)
  abioPkArea.list[[i]]<-abioPkArea
  fileId<-abioVendInfo.list[[i]]$file_id[1]
  #names(abioPkArea.list[[i]])<-fileId
  VIfiles.vec<-c(VIfiles.vec,fileId)
}
length(abioPkArea.list)
head(abioPkArea.list[[1]])
names(abioPkArea.list)<-VIfiles.vec
head(names(abioPkArea.list))

#abioAnalyses.vec<-abioFileInfo.df$Analysis
FIfiles<-abioFileInfo.df$file_id
head(FIfiles)
VIfiles<-names(abioPkArea.list)
head(VIfiles)
# same order!
head(unique(abio.dat$fileId))

# add additoinal values to data
# separate by analysis then add data for each experiment

abioSep.list<-separate_by_analysis_numDXF(abio.dat)
for(i in seq(1,length(abioSep.list))){
  currDf<-abioSep.list[[i]]
  # add peak areas -- add avg peak areas during tsfeature extraction
  # get sample peak numbers
  sampPkNrs<-abioSep.list[[i]]$PeakNr
  currDf$pkArea<-abioPkArea.list[[i]]$trapArea[sampPkNrs]
  # update abioSep.list[[i]]
  abioSep.list[[i]]$pkArea<-currDf$pkArea
  # add iso ration d18O/13C
  d18O13C.lab<-d18O13Clab(currDf)$d18O13C
  abioSep.list[[i]]$d18O13C<-d18O13C.lab
  # add biotic label
  abioLab<-rep("abiotic",dim(currDf)[1])
  abioSep.list[[i]]$biotic<-abioLab
  # add avgs of deltas
  deltaAvg.df<-avg3Deltas(currDf)
  d1<-rep(deltaAvg.df$avg_d13C12C,dim(currDf)[1])
  d2<-rep(deltaAvg.df$avg_d18O13C,dim(currDf)[1])
  d3<-rep(deltaAvg.df$avg_d18O16O,dim(currDf)[1])
  abioSep.list[[i]]$avg_d13C12C<-d1
  abioSep.list[[i]]$avg_d18O13C<-d2
  abioSep.list[[i]]$avg_d18O16O<-d3
  # add intensity avg values
  avg_Ampl46<-rep(mean(abioSep.list[[i]]$Ampl46),dim(currDf)[1])
  abioSep.list[[i]]$avg_Ampl46<-avg_Ampl46 # abs intensities can change dramatically based on volatge
  # add intensityAll avg values
  avg_IntensityAll<-rep(mean(abioSep.list[[i]]$IntensityAll),dim(currDf)[1])
  abioSep.list[[i]]$avg_IntensityAll<-avg_IntensityAll
  # add rIntensityAll avg values
  avg_rIntensityAll<-rep(mean(abioSep.list[[i]]$rIntensityAll),dim(currDf)[1])
  abioSep.list[[i]]$avg_rIntensityAll<-avg_rIntensityAll # ratio of intensities?
  # avg peak areas
  avg_pkArea<-rep(mean(abioSep.list[[i]]$pkArea),dim(currDf)[1])
  abioSep.list[[i]]$avg_pkArea<-avg_pkArea
}
head(abioSep.list[[1]])
dim(abioSep.list[[1]])
# [1]  9 54


## europa bacteria
setwd("~/Desktop/EuropaMLMS/microbe/QC_SU22/bio_abio/europa_bacteria")
eurbio.dat<-read.csv("europa_bacteria_samps_QC_2022-06-23.csv",header=T)
head(eurbio.dat)
dim(eurbio.dat)
# [1] 1061   44

# process Preparation labels
unique(eurbio.dat$Preparation)
# [1] "Bacteria Inoculation 0.3% CO2 in He 5 day"
# [2] "0.5% CO2 in He 5 day"

# replace spaces in Preparation with underscores
eurbio.dat<-replace_space_prep(eurbio.dat)
unique(eurbio.dat$Preparation)
# [1] "Bacteria_Inoculation_0.3%_CO2_in_He_5_day"
# [2] "0.5%_CO2_in_He_5_day"

# process Identifier1 labels --- ** ask Bethany
unique(eurbio.dat$Identifier1)
#[1] "CA"                              "DT-A"
#[3] "H1"                              "KCl L"
#[5] "KCl L + MgCl2 U"                 "KCl U + MgCl2 U"
#[7] "LW"                              "MgCl2 U"
#[9] "MgSo4 L"                         "MgSO4 L + Na2SO4 L + MgCl2 L"
#[11] "MgSo4 L + NaCl L"                "MgSO4 L + NaCl U"
#[13] "MgSO4 L + NaHCO3 L"              "Na2SO4 L"
#[15] "Na2SO4 L + NaCl L"               "Na2SO4 L + NaCl L + NaHCO3 L"
#[17] "Na2SO4 L + NaCl U"               "Na2SO4 L + NaHCO3 L"
#[19] "Na2SO4 U"                        "NaCl L"
#[21] "NaCl L + KCl U + MgCl2 L"        "NaCl U"
#[23] "NaCl U + KCl L + MgCl2 L"        "NaHCO3 L + NaCl L"
#[25] "C-A 50ul"                        "DT-A 300ul"
#[27] "DT-A 50ul"                       "KCL U C 50ul"
#[29] "KCL U DT 300ul"                  "MgCl2 L C 100ul"
#[31] "MgCl2 L C 300ul"                 "MgCl2 L C 50ul"
#[33] "MgCl2 U DT 300ul"                "MgSO4 L DT 100ul"
#[35] "MgSO4 L DT 50ul"                 "Na2SO4 L C 50ul"
#[37] "Na2SO4 L DT 300ul"               "Na2SO4 U C 300ul"
#[39] "Na2SO4 U DT 100ul"               "NaCl L C 100ul"
#[41] "NaCl L DT 300ul"                 "NaHCO3 L DT 100ul"
#[43] "C-A 100ul"                       "C-A 300ul"
#[45] "DT-A 100ul"                      "KCL L C 100ul"
#[47] "KCL L C 300ul"                   "KCL L C 50ul"
#[49] "KCL L DT 100ul"                  "KCL L DT 50ul"
#[51] "KCL U C 100ul"                   "KCL U C 300ul"
#[53] "KCL U DT 100ul"                  "KCL U DT 50ul"
#[55] "L1"                              "MgCl2 L DT 100ul"
#[57] "MgCl2 L DT 300ul"                "MgCl2 U C 100ul"
#[59] "MgCl2 U C 300ul"                 "MgCl2 U C 50ul"
#[61] "MgCl2 U DT 50ul"                 "MgSO4 L C 100ul"
#[63] "MgSO4 L C 300ul"                 "MgSO4 L C 50ul"
#[65] "MgSO4 L DT 300ul"                "Na2SO4 L C 100ul"
#[67] "Na2SO4 L C 300ul"                "Na2SO4 L DT 100ul"
#[69] "Na2SO4 L DT 50ul"                "Na2SO4 U C 100ul"
#[71] "Na2SO4 U C 50ul"                 "Na2SO4 U DT 300ul"
#[73] "Na2SO4 U DT 50ul"                "NaCl L C 300ul"
#[75] "NaCl L C 50ul"                   "NaCl L DT 100ul"
#[77] "NaCl L DT 50ul"                  "NaCl U C 100ul"
#[79] "NaCl U C 300ul"                  "NaCl U C 50ul"
#[81] "NaCl U DT 100ul"                 "NaCl U DT 300ul"
#[83] "NaCl U DT 50ul"                  "NaHCO3 L C 300ul"
#[85] "NaHCO3 L DT 300ul"               "NaHCO3 L DT 50ul"
#[87] "H1 C"                            "KCl U + MgCl2 L DT"
#[89] "KCl U C"                         "L1 C"
#[91] "MgCl2 U DT"                      "MgSO4 L + MgCl2 L DT"
#[93] "MgSO4 L + Na2SO4 L + MgCl2 L DT" "MgSO4 L + NaCl L DT"
#[95] "MgSO4 L + NaCl U C"              "MgSO4 L + NaHCO3 L DT"
#[97] "NaHCO3 L C"

# TODO: fix Identifier1 - not ready
#eurID1.map<-read.csv("europa_bacteria_ID1map.csv",header=T)

# just take out spaces for now
eurbioCopy.dat<-replace_space_Identifier1(eurbio.dat)
unique(eurbioCopy.dat$Identifier1)


# get raw intensity vs time data
setwd("/media/lily/My Book/Lilys_NASA_microbe_data/Europa_microbial_mud1/Europa_microbial_mud/")
mmRawList<-raw_data_all(mmPassedFiles)
head(mmRawList[[1]])
#    file_id tp time.s   v44.mV    v45.mV   v46.mV
#1 170418_H1_.dxf  1  0.209 1.017767 0.3216353 1.651687
#2 170418_H1_.dxf  2  0.418 1.019678 0.2892148 1.561420
#3 170418_H1_.dxf  3  0.627 1.031142 0.3368935 1.614990
#4 170418_H1_.dxf  4  0.836 1.025410 0.3502452 1.712875
#5 170418_H1_.dxf  5  1.045 1.023499 0.4093825 1.886962
#6 170418_H1_.dxf  6  1.254 1.027321 0.3388008 1.775992

# get file info for passed files
mmFileInfo.df<-file_info(mmPassedFiles)
head(mmFileInfo.df)
#             file_id Identifier1 Analysis Preparation       Date_and_Time
#1      170418_H1_.dxf          H1     2959   CO2 in He 2017-04-18 13:45:36
#2   170418_H1_(1).dxf          H1     2969   CO2 in He 2017-04-18 16:35:06
#3      170418_L1_.dxf          L1     2958   CO2 in He 2017-04-18 13:28:39
#4   170418_L1_(1).dxf          L1     2980   CO2 in He 2017-04-18 19:41:32
#5   170418_LW_(1).dxf          LW     2991   CO2 in He 2017-04-18 22:47:58
#6 170418_MgSO4_L_.dxf     MgSO4 L     2972   CO2 in He 2017-04-18 17:25:56

# get vendor data
mmVendInfo.list<-vendor_info_all(mmPassedFiles)
head(mmVendInfo.list[[1]])

# find peak areas

# file name is the name of the data frame
mmPkArea.list<-list()
VIfiles.vec<-c()
for(i in seq(1,length(mmRawList))){
  mmPkArea<-trap_area_allPks(raw.df=mmRawList[[i]],vend.df=mmVendInfo.list[[i]],
                             mV.rawName = rawVname)
  mmPkArea.list[[i]]<-mmPkArea
  fileId<-mmVendInfo.list[[i]]$file_id[1]
  VIfiles.vec<-c(VIfiles.vec,fileId)
}
names(mmPkArea.list)<-VIfiles.vec
mmPkArea.list[[1]]
#  PkNr  trapArea
#1    1 60.631344
#2    2 60.731316
#3    3  6.642827
#4    4 60.669124
#5    5 60.387140
#6    6 27.986229
#7     7 26.524341
#8     8 25.132668
#9     9 23.803885
#10   10 22.520561
#11   11 21.320711
#12   12 20.196341
#13   13 19.143684
#14   14 18.077922
#15   15 17.155240
#16   16 60.561875


# add additional values to data
# separate by analysis then add data for each experiment
mmSep.list<-separate_by_analysis_numDXF(mmbio.dat)
for(i in seq(1,length(mmSep.list))){
  currDf<-mmSep.list[[i]]
  # add peak areas -- add avg peak areas during tsfeature extraction
  # get sample peak numbers
  sampPkNrs<-mmSep.list[[i]]$PeakNr
  currDf$pkArea<-mmPkArea.list[[i]]$trapArea[sampPkNrs]
  # update abioSep.list[[i]]
  mmSep.list[[i]]$pkArea<-currDf$pkArea
  # add iso ration d18O/13C
  d18O13C.lab<-d18O13Clab(currDf)$d18O13C
  mmSep.list[[i]]$d18O13C<-d18O13C.lab
  # add biotic label
  mmLab<-rep("biotic",dim(currDf)[1])
  mmSep.list[[i]]$biotic<-mmLab
  # add avgs of deltas
  deltaAvg.df<-avg3Deltas(currDf)
  d1<-rep(deltaAvg.df$avg_d13C12C,dim(currDf)[1])
  d2<-rep(deltaAvg.df$avg_d18O13C,dim(currDf)[1])
  d3<-rep(deltaAvg.df$avg_d18O16O,dim(currDf)[1])
  mmSep.list[[i]]$avg_d13C12C<-d1
  mmSep.list[[i]]$avg_d18O13C<-d2
  mmSep.list[[i]]$avg_d18O16O<-d3
  # add intensity avg values
  avg_Ampl46<-rep(mean(mmSep.list[[i]]$Ampl46),dim(currDf)[1])
  mmSep.list[[i]]$avg_Ampl46<-avg_Ampl46
  # add intensityAll avg values
  avg_IntensityAll<-rep(mean(mmSep.list[[i]]$IntensityAll),dim(currDf)[1])
  mmSep.list[[i]]$avg_IntensityAll<-avg_IntensityAll
  # add rIntensityAll avg values
  avg_rIntensityAll<-rep(mean(mmSep.list[[i]]$rIntensityAll),dim(currDf)[1])
  mmSep.list[[i]]$avg_rIntensityAll<-avg_rIntensityAll
  # avg peak areas
  avg_pkArea<-rep(mean(mmSep.list[[i]]$pkArea),dim(currDf)[1])
  mmSep.list[[i]]$avg_pkArea<-avg_pkArea
}
head(mmSep.list[[1]])
dim(mmSep.list[[1]])
# [1]  9 54



##### europa bacteria
setwd("~/Desktop/EuropaMLMS/microbe/QC_SU22/bio_abio/europa_bacteria")

# get analyses that passed QC
eurPassedFiles<-unique(eurbio.dat$fileId)
head(eurPassedFiles)
#[1] "170725_CA_(1).dxf"   "170725_CA_(2).dxf"   "170725_DT-A_.dxf"    "170725_DT-A_(1).dxf"
#[5] "170725_DT-A_(2).dxf" "170725_H1_.dxf"

# get raw intensity vs time data
setwd("/media/lily/My Book/Lilys_NASA_microbe_data/Europa_bacteria1/Europa_bacteria/")
eurRawList<-raw_data_all(eurPassedFiles)
head(eurRawList[[1]])
# file_id tp time.s    v44.mV    v45.mV   v46.mV
#1 170725_CA_(1).dxf  1  0.209 0.9189074 0.1468496 2.122827
##2 170725_CA_(1).dxf  2  0.418 0.9284573 0.2668090 2.145792
#3 170725_CA_(1).dxf  3  0.627 0.9093578 0.2649022 2.143878
#4 170725_CA_(1).dxf  4  0.836 0.9093578 0.2706226 2.084629
#5 170725_CA_(1).dxf  5  1.045 0.9169974 0.2534618 2.105697
#6 170725_CA_(1).dxf  6  1.254 0.9112677 0.2420219 2.050157

# get file info for passed files
eurFileInfo.df<-file_info(eurPassedFiles)
head(eurFileInfo.df)
#             file_id Identifier1 Analysis  Preparation
#1   170725_CA_(1).dxf          CA     5759 Bacteria Inoculation 0.3% CO2 in He 5 day
#2   170725_CA_(2).dxf          CA     5770 Bacteria Inoculation 0.3% CO2 in He 5 day
#3    170725_DT-A_.dxf        DT-A     5765 Bacteria Inoculation 0.3% CO2 in He 5 day
#4 170725_DT-A_(1).dxf        DT-A     5778 Bacteria Inoculation 0.3% CO2 in He 5 day
#5 170725_DT-A_(2).dxf        DT-A     5783 Bacteria Inoculation 0.3% CO2 in He 5 day
#6      170725_H1_.dxf          H1     5754 Bacteria Inoculation 0.3% CO2 in He 5 day
#.       Date_and_Time
#1 2017-07-25 10:51:48
#2 2017-07-25 13:58:17
#3 2017-07-25 12:33:31
#4 2017-07-25 16:13:54
#5 2017-07-25 17:38:39
#6 2017-07-25 09:27:03

# get vendor data
eurVendInfo.list<-vendor_info_all(eurPassedFiles)
head(eurVendInfo.list[[1]])


# find peak areas
# file name is the name of the data frame
eurPkArea.list<-list()
VIfiles.vec<-c()
for(i in seq(1,length(eurRawList))){
  eurPkArea<-trap_area_allPks(raw.df=eurRawList[[i]],vend.df=eurVendInfo.list[[i]],
                               mV.rawName = rawVname)
  eurPkArea.list[[i]]<-eurPkArea
  fileId<-eurVendInfo.list[[i]]$file_id[1]
  VIfiles.vec<-c(VIfiles.vec,fileId)
}
names(eurPkArea.list)<-VIfiles.vec
head(eurPkArea.list[[1]])
#  PkNr trapArea
#1     1 51.18976
#2     2 51.28452
#3     3 15.54117
#4     4 51.23127
#5     5 50.98622
#6     6 66.56839


# TODO: make into a function
# add additional values to data
# separate by analysis then add data for each experiment
eurSep.list<-separate_by_analysis_numDXF(eurbio.dat)
for(i in seq(1,length(eurSep.list))){
  currDf<-eurSep.list[[i]]
  # add peak areas -- add avg peak areas during tsfeature extraction
  # get sample peak numbers
  sampPkNrs<-eurSep.list[[i]]$PeakNr
  currDf$pkArea<-eurPkArea.list[[i]]$trapArea[sampPkNrs]
  # update abioSep.list[[i]]
  eurSep.list[[i]]$pkArea<-currDf$pkArea
  # add iso ration d18O/13C
  d18O13C.lab<-d18O13Clab(currDf)$d18O13C
  eurSep.list[[i]]$d18O13C<-d18O13C.lab
  # add biotic label
  eurLab<-rep("biotic",dim(currDf)[1])
  eurSep.list[[i]]$biotic<-eurLab
  # add avgs of deltas
  deltaAvg.df<-avg3Deltas(currDf)
  d1<-rep(deltaAvg.df$avg_d13C12C,dim(currDf)[1])
  d2<-rep(deltaAvg.df$avg_d18O13C,dim(currDf)[1])
  d3<-rep(deltaAvg.df$avg_d18O16O,dim(currDf)[1])
  eurSep.list[[i]]$avg_d13C12C<-d1
  eurSep.list[[i]]$avg_d18O13C<-d2
  eurSep.list[[i]]$avg_d18O16O<-d3
  # add intensity avg values
  avg_Ampl46<-rep(mean(eurSep.list[[i]]$Ampl46),dim(currDf)[1])
  eurSep.list[[i]]$avg_Ampl46<-avg_Ampl46
  # add intensityAll avg values
  avg_IntensityAll<-rep(mean(eurSep.list[[i]]$IntensityAll),dim(currDf)[1])
  eurSep.list[[i]]$avg_IntensityAll<-avg_IntensityAll
  # add rIntensityAll avg values
  avg_rIntensityAll<-rep(mean(eurSep.list[[i]]$rIntensityAll),dim(currDf)[1])
  eurSep.list[[i]]$avg_rIntensityAll<-avg_rIntensityAll
  # avg peak areas
  avg_pkArea<-rep(mean(eurSep.list[[i]]$pkArea),dim(currDf)[1])
  eurSep.list[[i]]$avg_pkArea<-avg_pkArea
}
head(eurSep.list[[1]])
dim(eurSep.list[[1]])
# [1]  9 54



# all 3 datsets have the same dims now
# add salt and co2 labels after feature extraction?

# flatten lists to dfs, combine and write to file
abioProc.dat<-listToDF(abioSep.list)

mmbioProc.dat<-listToDF(mmSep.list)

eurbioProc.dat<-listToDF(eurSep.list)


colnames(eurbioProc.dat)
# [1] "fileId"            "Identifier1"       "Analysis"          "Preparation"
#[5] "DateTime"          "PeakNr"            "Start"             "Rt"
#[9] "End"               "Ampl44"            "Ampl45"            "Ampl46"
#[13] "BGD44"             "BGD45"             "BGD46"             "rIntensity44"
#[17] "rIntensity45"      "rIntensity46"      "rIntensityAll"     "Intensity44"
#[21] "Intensity45"       "Intensity46"       "IntensityAll"      "ListFirstPeak"
#[25] "rR45CO244CO2"      "rR46CO244CO2"      "IsRef"             "R45CO244CO2"
#[29] "RefName"           "rd45CO244CO2"      "d45CO244CO2"       "R46CO244CO2"
#[33] "rd46CO244CO2"      "d46CO244CO2"       "R13C12C"           "d13C12C"
#[37] "AT13C12C"          "R18O16O"           "d18O16O"           "AT18O16O"
#[41] "R17O16O"           "d17O16O"           "Rps45CO244CO2"     "Rps46CO244CO2"
#[45] "pkArea"            "d18O13C"           "biotic"            "avg_d13C12C"
#[49] "avg_d18O13C"       "avg_d18O16O"       "avg_Ampl46"        "avg_IntensityAll"
#[53] "avg_rIntensityAll" "avg_pkArea"


abioProc.dat$dataset<-rep("abiotic",dim(abioProc.dat)[1])
head(abioProc.dat)
unique(abioProc.dat$Preparation)
unique(abioProc.dat$Identifier1)

mmbioProc.dat$dataset<-rep("microbial_mud",dim(mmbioProc.dat)[1])
head(mmbioProc.dat)
unique(mmbioProc.dat$Preparation)
unique(mmbioProc.dat$Identifier1)

eurbioProc.dat$dataset<-rep("europa_bacteria",dim(eurbioProc.dat)[1])
head(eurbioProc.dat)
unique(eurbioProc.dat$Preparation)
eurBCopy<-replace_space_Identifier1(eurbioProc.dat)
eurbioProc.dat<-eurBCopy
unique(eurbioProc.dat$Identifier1)

###### ML labels
# process Identifier1 labels and make salt content and ion labels
setwd("~/Desktop/EuropaMLMS/microbe/QC_SU22/bio_abio/bioAbio/MLdataMappingFiles")

# read in mapping files
meta.dat<-read.csv("meta_data.csv",header=T)
head(meta.dat)

ts.dat<-read.csv("ms_ts_data.csv",header=T)
head(ts.dat)

head(bioAbio.dat)

abioID1.map<-read.csv("abiotic_fixIdentifier1_rxn_mapping.csv",header=T)
abioID1.map

mmID1.map<-read.csv("microbial_mud_fixIdentifier1_rxn_mapping.csv",header=T)
mmID1.map

eurID1.map<-read.csv("europa_bacteria_fixIdentifier1_rxn_mapping.csv",header=T)
head(eurID1.map)


# replace old Identifier1s with newIdentifier1
# Identifier1s in abiotic data already good
# fix Identifier1s in microbial mud data
newMMsamps.df<-mmbioProc.dat
unique(mmbioProc.dat$Identifier1)
#[1] "H1"                 "L1"                 "LW"                 "MgSO4_L"
#[5] "MgSO4_L_+_NaCl_L"   "MgSO4_L_+_NaCl_U"   "MgSO4_L_+_NaHCO3_L" "MgSO4_L_+_NaHCO3_U"
#[9] "NaCl_L"             "NaCl_U"             "NaHCO3_L"           "NaHCO3_L_+_NaCl_L"
#[13] "NaHCO3_L_+_NaCl_U"  "NaHCO3_U"           "NaHCO3_U_+_NaCl_L"  "NaHCO3_U_+_NaCl_U"

mmID1.map$newIdentifier1
# [1] "H1"                 "L1"                 "LW"                 "MgSO4_L"
#[5] "MgSO4_L_+_NaCl_L"   "MgSO4_L_+_NaCl_U"   "MgSO4_L_+_NaHCO3_L" "MgSO4_L_+_NaHCO3_U"
#[9] "NaCl_L"             "NaCl_U"             "NaHCO3_L"           "NaHCO3_L_+_NaCl_L"
#[13] "NaHCO3_L_+_NaCl_U"  "NaHCO3_U"           "NaHCO3_U_+_NaCl_L"  "NaHCO3_U_+_NaCl_U"

#[17] "MgSO4_L"            "MgSO4_L_+_NaCl_L"   "MgSO4_L_+_NaHCO3_L" "MgSO4_L_+_NaHCO3_U"
#[21] "NaCl_U"             "NaHCO3_L_+_NaCl_L"  "NaHCO3_L_+_NaCl_U"

inData<-which(unique(mmbioProc.dat$Identifier1)%in%mmID1.map$newIdentifier1)
unique(mmbioProc.dat$Identifier1)[-inData]

newMMsamps.df<-mmbioProc.dat
newID1.vec<-c()
for(i in seq(1,length(mmbioProc.dat$Identifier1))){
  #print(i)
  currID1<-mmbioProc.dat$Identifier1[i]
  currID1
  # match to mapping
  id1Ind<-which(mmID1.map$oldIdentifier1==currID1)
  if(length(id1Ind)==0){
    newID1<-currID1
    newID1
  }else{
    newID1<-mmID1.map$newIdentifier1[id1Ind]
    newID1
  }

  newID1.vec<-c(newID1.vec,newID1)
  newID1.vec
}
head(newID1.vec)
length(newID1.vec)
newMMsamps.df$Identifier1<-newID1.vec
head(newMMsamps.df)
unique(newMMsamps.df$Identifier1)
mmbioProc.df<-newMMsamps.df

# add ML labels
# add salt content label
# microbial mud
newMMsamps.df<-saltContentLab(mmbioProc.df)
unique(newMMsamps.df$saltContent)
# [1] "no_salt"      "MgSO4"        "MgSO4_NaCl"   "MgSO4_NaHCO3" "NaCl"
# [6] "NaHCO3"       "NaHCO3_NaCl"

head(newMMsamps.df)
tail(newMMsamps.df)

# ion labels
# TODO: make into func add to newFuncs file
cations.map<-c("Mg[2+]","Na[+]")
anions.map<-c("SO4[2-]","Cl[-]","CO3H[-]")
catVals.map<-c("Mg","Na") # for string matching
anVals.map<-c("SO4","Cl","HCO3")
catCont.vec<-c()
anCont.vec<-c()
for(i in seq(1,length(newMMsamps.df$saltContent))){
  currSaltCont<-newMMsamps.df$saltContent[i]
  if(currSaltCont=="no_salt"){
    catCont.vec<-c(catCont.vec,NA)
    anCont.vec<-c(anCont.vec,NA)
  }else{
    # split by underscore
    saltSplit<-strsplit(currSaltCont,"_",fixed=T)
    saltSplit<-saltSplit[[1]]
    # get salts, cations, and anions
    numSalts<-length(saltSplit)
    salts.mat<-matrix(rep(0,numSalts),ncol=numSalts)
    saltnames.vec<-c()
    # add to salts.mat
    for(j in seq(1,numSalts)){
      salts.mat[1,j]<-saltSplit[j]
      saltnames.vec<-c(saltnames.vec,paste("salt",j,sep=""))
    }
    salts.df<-as.data.frame(salts.mat)
    colnames(salts.df)<-saltnames.vec
    # ion matrix
    ions.mat<-matrix(rep(0,2*dim(salts.df)[2]),ncol=2)
    an.vec<-c()
    cat.vec<-c()
    # split salts into cations and anions
    for(k in seq(1,dim(salts.df)[2])){
      anion<-substr(salts.df[1,k],3,nchar(salts.df[1,k]))
      if(substr(anion,1,1)=="2"){ # if == "2SO4"
        anion<-substr(anion,2,nchar(anion))
      }
      anInd<-which(grepl(anion,anVals.map))
      catInd<-which(grepl(substr(salts.df[1,k],1,2),catVals.map))
      currAn<-anions.map[anInd]
      currCat<-cations.map[catInd]
      an.vec<-c(an.vec,anions.map[anInd])
      cat.vec<-c(cat.vec,cations.map[catInd])
    }
    ions.mat[,1]<-cat.vec
    ions.mat[,2]<-an.vec
    ions.df<-as.data.frame(ions.mat)
    colnames(ions.df)<-c("cations","anions")
    # paste cations and anions together for labels
    anString<-""
    catString<-""
    for(l in seq(1,length(an.vec))){
      #ind.vec<-seq(1,l)
      anString<-paste(anString,"_",an.vec[l],sep="")
      catString<-paste(catString,"_",cat.vec[l],sep="")
      #rm first underscore
      #anString<-substr(anString,2,nchar(anString))
      #catString<-substr(catString,2,nchar(catString))
    }
    # add to catCont and anCont vecs
    catCont.vec<-c(catCont.vec,substr(catString,2,nchar(catString)))
    anCont.vec<-c(anCont.vec,substr(anString,2,nchar(anString)))
  }

}
head(anCont.vec)
tail(anCont.vec)
length(anCont.vec)
unique(anCont.vec)
#[1] NA                "SO4[2-]"         "SO4[2-]_Cl[-]"   "SO4[2-]_CO3H[-]"
#[5] "Cl[-]"           "CO3H[-]"         "CO3H[-]_Cl[-]"
head(catCont.vec)
tail(catCont.vec)
length(catCont.vec)
unique(catCont.vec)
#[1] NA             "Mg[2+]"       "Mg[2+]_Na[+]" "Na[+]"        "Na[+]_Na[+]"
# add cation and anion labels to df
newMMsamps.df$anion<-anCont.vec
newMMsamps.df$cation<-catCont.vec
head(newMMsamps.df)
tail(newMMsamps.df)

mmbioProc.dat<-newMMsamps.df
head(mmbioProc.dat)
tail(mmbioProc.dat)


## abiotic data
# add salt content label
newAbioSamps.df<-saltContentLab(abioProc.dat)
head(newAbioSamps.df)
unique(newAbioSamps.df$saltContent)
# [1] "no_salt"       "MgSO4_NaCl"    "MgSO4_NaHCO3"  "MgSO4"         "Na2SO4_NaCl"
# [6] "Na2SO4_NaHCO3" "Na2SO4"        "NaCl"          "NaHCO3_NaCl"   "NaHCO3"

# ion labels - Identifier1 is fine
cations.map<-c("Mg[2+]","Na[+]")
anions.map<-c("SO4[2-]","Cl[-]","CO3H[-]") #TODO: switch back to HCO3[-]
catVals.map<-c("Mg","Na") # for string matching
anVals.map<-c("SO4","Cl","HCO3")# TODO: CO3H
catCont.vec<-c()
anCont.vec<-c()
for(i in seq(1,length(newAbioSamps.df$saltContent))){
  currSaltCont<-newAbioSamps.df$saltContent[i]
  if(currSaltCont=="no_salt"){
    catCont.vec<-c(catCont.vec,NA)
    anCont.vec<-c(anCont.vec,NA)
  }else{
    # split by underscore
    saltSplit<-strsplit(currSaltCont,"_",fixed=T)
    saltSplit<-saltSplit[[1]]
    # get salts, cations, and anions
    numSalts<-length(saltSplit)
    salts.mat<-matrix(rep(0,numSalts),ncol=numSalts)
    saltnames.vec<-c()
    # add to salts.mat
    for(j in seq(1,numSalts)){
      salts.mat[1,j]<-saltSplit[j]
      saltnames.vec<-c(saltnames.vec,paste("salt",j,sep=""))
    }
    salts.df<-as.data.frame(salts.mat)
    colnames(salts.df)<-saltnames.vec
    # ion matrix
    ions.mat<-matrix(rep(0,2*dim(salts.df)[2]),ncol=2)
    an.vec<-c()
    cat.vec<-c()
    # split salts into cations and anions
    for(k in seq(1,dim(salts.df)[2])){
      anion<-substr(salts.df[1,k],3,nchar(salts.df[1,k]))
      if(substr(anion,1,1)=="2"){ # if == "2SO4"
        anion<-substr(anion,2,nchar(anion))
      }
      anInd<-which(grepl(anion,anVals.map))
      catInd<-which(grepl(substr(salts.df[1,k],1,2),catVals.map))
      currAn<-anions.map[anInd]
      currCat<-cations.map[catInd]
      an.vec<-c(an.vec,anions.map[anInd])
      cat.vec<-c(cat.vec,cations.map[catInd])
    }
    ions.mat[,1]<-cat.vec
    ions.mat[,2]<-an.vec
    ions.df<-as.data.frame(ions.mat)
    colnames(ions.df)<-c("cations","anions")
    # paste cations and anions together for labels
    anString<-""
    catString<-""
    for(l in seq(1,length(an.vec))){
      #ind.vec<-seq(1,l)
      anString<-paste(anString,"_",an.vec[l],sep="")
      catString<-paste(catString,"_",cat.vec[l],sep="")
      #rm first underscore
      #anString<-substr(anString,2,nchar(anString))
      #catString<-substr(catString,2,nchar(catString))
    }
    # add to catCont and anCont vecs
    catCont.vec<-c(catCont.vec,substr(catString,2,nchar(catString)))
    anCont.vec<-c(anCont.vec,substr(anString,2,nchar(anString)))
  }

}
head(anCont.vec)
tail(anCont.vec)
length(anCont.vec)
unique(anCont.vec)
#[1] NA                "SO4[2-]_Cl[-]"   "SO4[2-]_CO3H[-]" "SO4[2-]"
#[5] "Cl[-]"           "CO3H[-]_Cl[-]"   "CO3H[-]"
head(catCont.vec)
tail(catCont.vec)
length(catCont.vec)
unique(catCont.vec)
#[1] NA             "Mg[2+]_Na[+]" "Mg[2+]"       "Na[+]_Na[+]"  "Na[+]"
# add cation and anion labels to df
newAbioSamps.df$anion<-anCont.vec
newAbioSamps.df$cation<-catCont.vec
head(newAbioSamps.df)
tail(newAbioSamps.df)

abioProc.dat<-newAbioSamps.df
head(abioProc.dat)
tail(abioProc.dat)


### europa bacteria
setwd("~/Desktop/EuropaMLMS/microbe/QC_SU22/bio_abio/bioAbio/MLdataMappingFiles")
# read in ID1 mapping file
eurID1.map<-read.csv("europa_bacteria_fixIdentifier1_rxn_mapping.csv",header=T)
head(eurID1.map)

# fix Identifier1
newEurSamps.df<-eurbioProc.dat
dim(newEurSamps.df)
newID1.vec<-c()
eurRxn.vec<-c()
for(i in seq(1,length(newEurSamps.df$Identifier1))){#length(newEurSamps.df$Identifier1)
  currID1<-newEurSamps.df$Identifier1[i]
  # match to mapping
  id1Ind<-which(eurID1.map$oldIdentifier1==currID1)
  if(length(id1Ind)==0){
    newID1<-currID1
    newID1.vec<-c(newID1.vec,newID1)
    newID1.vec
  }else{
    newID1<-eurID1.map$newIdentifier1[id1Ind]
    newID1
    newID1.vec<-c(newID1.vec,newID1)
    newID1.vec
  }

}
length(newID1.vec)
unique(newID1.vec)
write.table(unique(newID1.vec),"eurNewID1.tab",quote=F,row.names=F)

newEurSamps.df$Identifier1<-newID1.vec
head(newEurSamps.df)
unique(newEurSamps.df$Identifier1)
eurbioProc.dat<-newEurSamps.df
head(eurbioProc.dat)

# add salt content and cation/anion labels -- doing manually for this dataset
setwd("~/Desktop/EuropaMLMS/microbe/QC_SU22")
saltIonMapEur<-read.csv("europaBactID1saltMap.csv",header=T)
head(saltIonMapEur)
unique(saltIonMapEur$Identifier1)
unique(saltIonMapEur$rxn)
# [1] "C"                            "DT"
#[3] "H1"                           "KCl_L"
#[5] "KCl_L_+_MgCl2_U"              "KCl_U_+_MgCl2_U"
#[7] "LW"                           "MgCl2_U"
#[9] "MgSO4_L"                      "MgSO4_L_+_Na2SO4_L_+_MgCl2_L"
#[11] "MgSO4_L_+_NaCl_L"             "MgSO4_L_+_NaCl_U"
#[13] "MgSO4_L_+_NaHCO3_L"           "Na2SO4_L"
#[15] "Na2SO4_L_+_NaCl_L"            "Na2SO4_L_+_NaCl_L_+_NaHCO3_L"
#[17] "Na2SO4_L_+_NaCl_U"            "Na2SO4_L_+_NaHCO3_L"
#[19] "Na2SO4_U"                     "NaCl_L"
#[21] "NaCl_L_+_KCl_U_+_MgCl2_L"     "NaCl_U"
#[23] "NaCl_U_+_KCl_L_+_MgCl2_L"     "NaHCO3_L_+_NaCl_L"
#[25] "KCl_U"                        "MgCl2_L"
#[27] "NaHCO3_L"                     "L1"
#[29] "KCl_U_+_MgCl2_L"              "MgSO4_L_+_MgCl2_L"
#[31] "MgSO4_+_Na2SO4_L_+_MgCl2_L"

# add data to eurbioProc.dat
rxn.vec<-c()
saltContent.vec<-c()
cation.vec<-c()
anion.vec<-c()
for(i in seq(1,dim(eurbioProc.dat)[1])){
  # match Identifier1 with rxn and add to df
  currID1<-eurbioProc.copy$Identifier1[i]
  currID1
  # use map
  # get rxn
  mapRxnInd<-which(saltIonMapEur$Identifier1==currID1)
  mapRxnInd
  rxnMap<-saltIonMapEur$rxn[mapRxnInd]
  rxnMap
  # add to vec
  rxn.vec<-c(rxn.vec,rxnMap)
  # get saltContent, cation and anion
  salt<-saltIonMapEur$saltContent[mapRxnInd]
  cation<-saltIonMapEur$cation[mapRxnInd]
  anion<-saltIonMapEur$anion[mapRxnInd]
  # add to vecs
  saltContent.vec<-c(saltContent.vec,salt)
  cation.vec<-c(cation.vec,cation)
  anion.vec<-c(anion.vec,anion)
}
length(rxn.vec) # all there!
length(saltContent.vec)
length(cation.vec)
length(anion.vec)

unique(anion.vec)
#[1] "no_anion"                "Cl[-]"
#[3] "Cl[-]_Cl2[2-]"           "Cl2[2-]"
#[5] "SO4[2-]"                 "SO4[2-]_SO4[2-]_Cl2[2-]"
#[7] "SO4[2-]_Cl[-]"           "SO4[2-]_HCO3[-]"
#[9] "SO4[2-]_Cl[-]_HCO3[-]"   "Cl[-]_Cl[-]_Cl2[2-]"
#[11] "HCO3[-]_Cl[-]"          "HCO3[-]"
#[13] "SO4[2-]_Cl2[2-]"

unique(cation.vec)
# [1] "no_cation"             "K[+]"                  "K[+]_Mg[2+]"
#[4] "Mg[2+]"                "Mg[2+]_Na2[2+]_Mg[2+]" "Mg[2+]_Na[+]"
#[7] "Na2[2+]"               "Na2[2+]_Na[+]"         "Na2[2+]_Na[+]_Na[+]"
#[10] "Na[+]"                 "Na[+]_K[+]_Mg[2+]"     "Na[+]_Na[+]"
#[13] "Mg[2+]_Mg[2+]"         "Mg[2+]_Na[+]_Mg[2+]"

# add to Proc df
eurbioProc.dat$rxn<-rxn.vec
eurbioProc.dat$saltContent<-saltContent.vec
eurbioProc.dat$cation<-cation.vec
eurbioProc.dat$anion<-anion.vec
head(eurbioProc.dat)


# check other datasets salt cols
# abiotic
head(abioProc.dat) # replace NA with "no_cation" and "no_anion"
abioProc.dat$anion[which(is.na(abioProc.dat$anion))]<-"no_anion"
head(abioProc.dat$anion)
abioProc.dat$cation[which(is.na(abioProc.dat$cation))]<-"no_cation"

# microbial mud
head(mmbioProc.dat)
mmbioProc.dat$anion[which(is.na(mmbioProc.dat$anion))]<-"no_anion"
mmbioProc.dat$cation[which(is.na(mmbioProc.dat$cation))]<-"no_cation"
head(mmbioProc.dat)

# mm and abio need their rxn cols added
# for these datasets rxn is the same as Identifier1
mmRxn.vec<-c()
for(i in seq(1,dim(mmbioProc.dat)[1])){#dim(mmbioProc.dat)[1]
  currID1<-mmbioProc.copy$Identifier1[i]
  #rxnInd<-which(mmID1.map$rxn==currID1)
  #rxn<-mmID1.map$rxn[rxnInd]
  rxn<-currID1
  print(rxn)
  mmRxn.vec<-c(mmRxn.vec,rxn)
  mmRxn.vec
}
head(mmRxn.vec)
dim(mmbioProc.dat)[1]
length(mmRxn.vec)
tail(mmRxn.vec)
# add to df
mmbioProc.dat$rxn<-mmRxn.vec
head(mmbioProc.dat)

unique(mmbioProc.dat$rxn) # same as Identifier1
# [1] "H1"                 "L1"                 "LW"
#[4] "MgSO4_L"            "MgSO4_L_+_NaCl_L"   "MgSO4_L_+_NaCl_U"
#[7] "MgSO4_L_+_NaHCO3_L" "MgSO4_L_+_NaHCO3_U" "NaCl_L"
#[10] "NaCl_U"             "NaHCO3_L"           "NaHCO3_L_+_NaCl_L"
#[13] "NaHCO3_L_+_NaCl_U"  "NaHCO3_U"           "NaHCO3_U_+_NaCl_L"
#[16] "NaHCO3_U_+_NaCl_U"


abioRxn.vec<-c()
for(i in seq(1,dim(abioProc.dat)[1])){
  currID1<-abioProc.dat$Identifier1[i]
  #rxnInd<-which(abioID1.map$rxn==currID1)
  #rxn<-abioID1.map$rxn[rxnInd]
  rxn<-currID1
  abioRxn.vec<-c(abioRxn.vec,rxn)
}
length(abioRxn.vec)
dim(abioProc.dat)
abioProc.dat$rxn<-abioRxn.vec
unique(abioProc.dat$rxn) # same as Identifier1
#[1] "H1"                  "L1"                  "LW"
#[4] "MgSO4_L_+_NaCl_L"    "MgSO4_L_+_NaCl_U"    "MgSO4_L_+_NaHCO3_L"
#[7] "MgSO4_L"             "Na2SO4_L_+_NaCl_L"   "Na2SO4_L_+_NaCl_U"
#[10] "Na2SO4_L_+_NaHCO3_L" "Na2SO4_L_+_NaHCO3_U" "Na2SO4_L"
#[13] "Na2SO4_U"            "NaCl_L"              "NaCl_U"
#[16] "NaHCO3_L_+_NaCl_L"   "NaHCO3_L_+_NaCl_U"   "NaHCO3_L"
#[19] "NaHCO3_U_+_NaCl_L"   "NaHCO3_U_+_NaCl_U"   "NaHCO3_U"


# need to reorder columns so they are all in the same order for the 3 datasets
# these datasets are the same
colnames(mmbioProc.dat)==colnames(abioProc.dat)
# match these then determine a final order after all processing
colnames(eurbioProc.dat)==colnames(mmProc.dat)
# not equal
# eurbiop
# [55] "dataset"           "rxn"               "saltContent"
#[58] "cation"            "anion"
# mmbio
# [55] "dataset"           "saltContent"       "anion"
#[58] "cation"            "rxn"


eurbio1<-eurbioProc.copy[,seq(1:55)] # these are fine for now
eurbioRxn<-eurbioProc.dat$rxn
eurbioSalt<-eurbioProc.dat$saltContent
eurbioCat<-eurbioProc.dat$cation
eurbioAn<-eurbioProc.dat$anion

eurbio2<-cbind.data.frame(eurbioSalt,eurbioAn,eurbioCat,eurbioRxn)
colnames(eurbio2)<-c("saltContent","anion","cation","rxn")
head(eurbio2)

eurbioProc.copy<-cbind.data.frame(eurbio1,eurbio2)
head(eurbioProc.copy)
colnames(eurbioProc.copy)==colnames(mmbioProc.dat)
eurbioProc.dat<-eurbioProc.copy
# now they all match


# merge these datasets, then process for individual ion labels
bioAbio.dat<-rbind.data.frame(abioProc.dat,mmbioProc.dat,eurbioProc.dat)
dim(bioAbio.dat)
#[1] 4112   59
# now have about 175 biotic experiments and 281 abiotic

bioAbio.dat<-read.csv("bioAbio_samp_ML_labels_newPrepBiCarb.csv",header=T)
head(bioAbio.dat)
#### need to actually remove the cation and anion labels!
unique(bioAbio.dat$allAnion)
unique(bioAbio.dat$allCation)
which(colnames(bioAbio.dat)=="cation")
which(colnames(bioAbio.dat)=="anion")
# make cols for each unique cation/anion in the datasets
# need to make vecs for cations and anions -- these will be input for processing
# cations
unique(mmbioProc.dat$cation)
unique(abioProc.dat$cation)
unique(eurbioProc.dat$cation)

allCats.vec<-c("no_cation", "Mg[2+]_Na[+]", "Mg[2+]","Na[+]_Na[+]",
              "Na[+]","K[+]", "K[+]_Mg[2+]", "Mg[2+]_Na2[2+]_Mg[2+]",
              "Na2[2+]", "Na2[2+]_Na[+]","Na2[2+]_Na[+]_Na[+]",
              "Na[+]_K[+]_Mg[2+]","Mg[2+]_Mg[2+]","Mg[2+]_Na[+]_Mg[2+]")

# anions
unique(mmbioProc.dat$anion)
# change CO3H's to HCO3's
# microbial mud
mmbioProc.dat$anion[which(mmbioProc.dat$anion=="CO3H[-]")]<-"HCO3[-]"
mmbioProc.dat$anion[which(mmbioProc.dat$anion=="SO4[2-]_CO3H[-]")]<-"SO4[2-]_HCO3[-]"
mmbioProc.dat$anion[which(mmbioProc.dat$anion=="CO3H[-]_Cl[-]")]<-"HCO3[-]_Cl[-]"
unique(mmbioProc.dat$anion)
# [1] "no_anion"        "SO4[2-]"         "SO4[2-]_Cl[-]"   "SO4[2-]_HCO3[-]"
#[5] "Cl[-]"           "HCO3[-]"         "HCO3[-]_Cl[-]"

# abiotic
unique(abioProc.dat$anion)
abioProc.dat$anion[which(abioProc.dat$anion=="CO3H[-]")]<-"HCO3[-]"
abioProc.dat$anion[which(abioProc.dat$anion=="CO3H[-]_Cl[-]")]<-"HCO3[-]_Cl[-]"
abioProc.dat$anion[which(abioProc.dat$anion=="SO4[2-]_CO3H[-]")]<-"SO4[2-]_HCO3[-]"
unique(abioProc.dat$anion)
#[1] "no_anion"        "SO4[2-]_Cl[-]"   "SO4[2-]_HCO3[-]" "SO4[2-]"
#[5] "Cl[-]"           "HCO3[-]_Cl[-]"   "HCO3[-]"

unique(eurbioProc.dat$anion)
#[1] "no_anion"                "Cl[-]"
#[3] "Cl[-]_Cl2[2-]"           "Cl2[2-]"
#[5] "SO4[2-]"                 "SO4[2-]_SO4[2-]_Cl2[2-]"
#[7] "SO4[2-]_Cl[-]"           "SO4[2-]_HCO3[-]"
#[9] "SO4[2-]_Cl[-]_HCO3[-]"   "Cl[-]_Cl[-]_Cl2[2-]"
#[11] "HCO3[-]_Cl[-]"           "HCO3[-]"
#[13] "SO4[2-]_Cl2[2-]"

allAns.vec<-c("no_anion", "SO4[2-]", "SO4[2-]_Cl[-]", "SO4[2-]_HCO3[-]",
              "Cl[-]", "HCO3[-]", "HCO3[-]_Cl[-]", "Cl[-]_Cl2[2-]",
              "Cl2[2-]", "SO4[2-]_SO4[2-]_Cl2[2-]",
              "SO4[2-]_Cl[-]_HCO3[-]","Cl[-]_Cl[-]_Cl2[2-]",
              "SO4[2-]_Cl2[2-]")

# all individual cations -- need cols to match for all datasets
indivCats.vec<-c("Mg[2+]", "Na[+]", "K[+]", "Na2[2+]")

# all individual anions
indivAns.vec<-c("SO4[2-]", "Cl[-]", "HCO3[-]","Cl2[2-]")

# change colnames


# create df for ion vals to be added to bioAbio data
ionColNames<-c("sulfate","chloride","carbonate","magnesium","sodium",
               "potassium","anion","cation")
ion.df<-as.data.frame(matrix(rep(NA,length(ionColNames)*dim(bioAbio.dat)[1]),
                      ncol=length(ionColNames)))
colnames(ion.df)<-ionColNames

# create ion mappings
# cations
indivCatsLab.vec<-c("magnesium","sodium","potassium","sodium")
indivAns.df<-cbind.data.frame(indivAns.vec,indivAnsLab.vec)
colnames(indivAns.df)<-c("ion","name")
indivAns.df #mapping of ion to colname in ion.df
#     ion      name
#1 SO4[2-]   sulfate
#2   Cl[-]  chloride
#3 HCO3[-] carbonate
#4 Cl2[2-]  chloride

# anions
indivAnsLab.vec<-c("sulfate","chloride","carbonate","chloride")
indivCats.df<-cbind.data.frame(indivCats.vec,indivCatsLab.vec)
colnames(indivCats.df)<-c("ion","name")
indivCats.df # mappin gof ion to colname in ion.df
# ion      name
#1  Mg[2+] magnesium
#2   Na[+]    sodium
#3    K[+] potassium
#4 Na2[2+]    sodium

catColNames<-c("magnesium","sodium","potassium")
anColNames<-c("sulfate","chloride","carbonate")
long.dat<-bioAbio.dat
# add vals to df
for(i in seq(1,dim(bioAbio.dat)[1])){
  currFullCat<-bioAbio.dat$anion[i]
  currFullCat
  currFullAn<-bioAbio.dat$cation[i]
  currFullAn
  # check if no_cation
  if(currFullCat=="no_cation"){ # add neg values to each cation column
    # add no_cation to cation col
    ion.df[i,which(colnames(ion.df)=="cation")]<-"no_cation"
    # and NA to every other cation col -- already NA!
    #ion.df[i,which(colnames(ion.df)%in%catColNames)]<-NA
    ion.df[i,]
  }else{ # search for matching cations and put them in the proper column
    splitCats<-strsplit(currFullCat,"_")
    splitCats[[1]]
    for(j in seq(1,length(splitCats[[1]]))){# loop through all cations
      currCat<-splitCats[[1]][j]
      currCat
      # find which cation column it belongs in
      catMatchColInd<-which(indivCats.df$ion==currCat)
      catMatchColName<-indivCats.df$name[catMatchColInd]
      catMatchColName

      ionDFind<-which(colnames(ion.df)==catMatchColName)
      # NA for all other ions - don't need to change
      ion.df[i,ionDFind]<-catMatchColName
      ion.df[i,]

      # add to cation col
      ion.df[i,8]<-currCat
      }
  }
  if(currFullAn=="no_anion"){ # add neg values to each anion column
    ion.df[i,which(colnames(ion.df)=="anion")]<-"no_anion"
    ion.df[i,]
  }else{ # search for matching anions and put them in the proper col
    splitAns<-strsplit(currFullAn,"_")
    splitAns[[1]]
    for(k in seq(1,length(splitAns[[1]]))){#loop through all anions
      currAn<-splitAns[[1]][k]
      currAn
      # find which anion col it belongs in
      anMatchColInd<-which(indivAns.df$ion==currAn)
      anMatchColName<-indivAns.df$name[anMatchColInd]
      anMatchColName

      ionDFind<-which(colnames(ion.df)==anMatchColName)
      # NA for all other ions - don't need to change
      ion.df[i,ionDFind]<-anMatchColName
      ion.df[i,]

      # add to anion col
      ion.df[i,7]<-currAn
    }
  }
}
head(ion.df)
tail(ion.df)

# use mappings to replace NAs with "no_x"
sulfateNA.ind<-which(is.na(ion.df$sulfate))
ion.df$sulfate[sulfateNA.ind]<-"no_sulfate"
head(ion.df)

chlorideNA.ind<-which(is.na(ion.df$chloride))
ion.df$chloride[chlorideNA.ind]<-"no_chloride"
head(ion.df)

carbNA.ind<-which(is.na(ion.df$carbonate))
ion.df$carbonate[carbNA.ind]<-"no_carbonate"
head(ion.df)

magNA.ind<-which(is.na(ion.df$magnesium))
ion.df$magnesium[magNA.ind]<-"no_magnesium"
head(ion.df)

sodNA.ind<-which(is.na(ion.df$sodium))
ion.df$sodium[sodNA.ind]<-"no_sodium"
head(ion.df)

potNA.ind<-which(is.na(ion.df$potassium))
ion.df$potassium[potNA.ind]<-"no_potassium"
head(ion.df)
tail(ion.df)

unique(ion.df$sulfate)
unique(ion.df$chloride)
unique(ion.df$carbonate)
unique(ion.df$magnesium)
unique(ion.df$sodium)
unique(ion.df$potassium)
unique(ion.df$anion)
unique(ion.df$cation)


# add ion.df to bioAbiodata
bioAbio.copy<-bioAbio.dat
bioAbio.copy<-cbind.data.frame(bioAbio.dat,ion.df)
head(bioAbio.copy)
tail(bioAbio.copy)
bioAbio.dat<-bioAbio.copy

# and add ph/ionic strength--Brett's code -- need to send him these files!
# need to fix DateTime again -- Brett's code, sent him the file already

# write this data to file
setwd("~/Desktop/EuropaMLMS/microbe/QC_SU22/bio_abio/bioAbio")
write.table(bioAbio.dat,"bioAbio_samps_ML_labels_2022-06-27.csv",quote=F,row.names=F,sep=",")


# process out internal standards and put into separate file...
# use Brett's code? and rxn column
inStandID<-c("H1","L1","LW")




#### ***
# TODO: pH and ionic strength for abioitc data -- process for Brett's code!!
setwd("~/Desktop/EuropaMLMS/microbe/QC_SU22/bio_abio/bioAbio/MLdataMappingFiles")
pHionicAbio.dat<-read.csv("pHionicAbio.csv",header=T)
colnames(pHionicAbio.dat)<-c("pH","Identifier1","ionicStrength[M]")
head(pHionicAbio.dat)
# fix Identifier1 spaces
pHionicAbio.dat<-replace_space_Identifier1(pHionicAbio.dat)
pHionicAbio.dat$dataset<-rep("abiotic",dim(pHionicAbio.dat)[1])
pHionicAbio.dat<-pHionicAbio.dat[,c(1,2,3,5)]
head(pHionicAbio.dat)
#   pH Identifier1 ionicStrength[M] dataset
#1 4.5       KCl_L             0.07 abiotic
#2 4.5       KCl_U             0.94 abiotic
#3 4.5     MgCl2_L             0.16 abiotic
#4 5.0     MgCl2_U             2.21 abiotic
#5 3.5     MgSO4_L             7.98 abiotic
#6 3.5    Na2SO4_L             5.28 abiotic


# write to file
write.table(pHionicAbio.dat,"pHionic_abiotic.csv",quote=F,row.names=F,sep=",")

# europa_bacteria
pHionicBio.dat<-read.csv("pHionicBio.csv",header=T)
colnames(pHionicBio.dat)<-c("Identifier1","pH")
pHionicBio.dat<-replace_space_Identifier1(pHionicBio.dat)
head(pHionicBio.dat)
dim(pHionicBio.dat)

#colnames(pHionicBio.dat)<-c("pH","rxn","ionicStrength[M]")
# reorder this data to match the data above; add ionicStrength col
pHionicBioCopy<-pHionicBio.dat
pHionicBioCopy[,1]<-pHionicBio.dat$pH
pHionicBioCopy[,2]<-pHionicBio.dat$Identifier1
pHionicBioCopy[,3]<-rep(0,dim(pHionicBioCopy)[1])
colnames(pHionicBioCopy)<-c("pH","Identifier1","ionicStrength[M]")
pHionicBioCopy$dataset<-rep("europa_bacteria",dim(pHionicBio.dat)[1])
pHionicBio.dat<-pHionicBioCopy

head(pHionicBio.dat)
#   pH    Identifier1 ionicStrength[M] #TODO:0
#1 8.5             CA         unkIonic
#2 9.5             DT         unkIonic
#3 5.5  KCL_L_C_100ul         unkIonic
#4 6.5  KCL_L_C_300ul         unkIonic
#5 6.0   KCL_L_C_50ul         unkIonic
#6 6.5 KCL_L_DT_100ul         unkIonic

# pH    Identifier1 ionicStrength[M]         dataset
#1 8.5             CA                0 europa_bacteria
#2 9.5             DT                0 europa_bacteria
#3 5.5  KCL_L_C_100ul                0 europa_bacteria
#4 6.5  KCL_L_C_300ul                0 europa_bacteria
#5 6.0   KCL_L_C_50ul                0 europa_bacteria
#6 6.5 KCL_L_DT_100ul                0 europa_bacteria

# write to file
write.table(pHionicBio.dat,"pHionic_europa_bacteria.csv",row.names=F,quote=F,sep=",")


mm_pHionic.df<-as.data.frame(matrix(rep(NA,4*length(unique(mmbioProc.dat$Identifier1))),ncol=4))
colnames(mm_pHionic.df)<-c("pH","Identifier1","ionicStrength[M]","dataset")
mm_pHionic.df$pH<-"unk"
mm_pHionic.df$Identifier1<-unique(mmbioProc.dat$Identifier1)
mm_pHionic.df$`ionicStrength[M]`<-0
mm_pHionic.df$dataset<-"microbial_mud"
head(mm_pHionic.df)
# write to file
write.table(mm_pHionic.df,"pHionic_microbial_mud.csv",row.names=F,quote=F,sep=",")





# TODO: npdr - plot chromatograms with high scores for important features
#abioID1.vec<-abioFileInfo.df$Identifier1
#abioDate.vec<-as.character(abioFileInfo.df$Date_and_Time)
#newDate.vec<-c()
#for(i in seq(1, length(abioDate.vec))){
#  currDate<-substr(abioDate.vec[i],1,10)
# update with yr-mo-day
#  newDate.vec[i]<-currDate
#}
#abioDate.vec<-newDate.vec
# plot all raw data
#setwd(qcOut.path)
#generic_plot_all_raw(rawList,fileName=paste(datasetLabel,"_failedQC_chromatograms.pdf",sep=""),
#                     analysis.vec=failedAnalyses.vec, ID1.vec=failedID1.vec,date.vec=failedDate.vec)

# dplyr to select sample IDs
# subject.attrs

# TODO: serialize labels for general data per Steve's suggestion
# makes it easier to generalize code
# makes it easier to plot, cluster, etc...


### fix dateTimes

# read in dxf files and make df of Analyses and DateTime
setwd("~/Desktop/EuropaMLMS/microbe/QC_SU22/bio_abio/bioAbio/bioAbioDXF")
dxfFiles<-all_dxf_files()
dxfFI<-file_info_all(dxfFiles)
head(dxfFI[[1]])
class(dxfFI[[1]]$Date_and_Time)

passAnalyses<-unique(bioAbio.dat$Analysis)
analysis.vec<-c()
dateTime.vec<-c()
for(i in seq(1,length(dxfFI))){
  currAn<-dxfFI[[i]]$Analysis
  if(currAn %in% passAnalyses){
    analysis.vec<-c(analysis.vec,currAn)
    currDate<-as.character(dxfFI[[i]]$Date_and_Time)
    newDate<-gsub(" ","_",currDate)
    dateTime.vec<-c(dateTime.vec,newDate)
  }

}
head(dateTime.vec)
length(dateTime.vec)

anDate.df<-as.data.frame(cbind(analysis.vec,dateTime.vec))
head(anDate.df)
setwd("~/Desktop/EuropaMLMS/microbe/QC_SU22/bio_abio/bioAbio/MLdataMappingFiles")
write.table(anDate.df,"analysisDateTimeBioAbio.csv",row.names=F,quote=F,sep=",")


## update pH and ionic strength file


############# TODO: update Preparations and have Brett fix the prepTime col for all datastes
# prepTimes = 5 days need to be 7 days
# CO2 conc for microbial mud is 0.3pCO2 for all
# CA and DT stay at 0 ionic strength, all others match rxn essentially....
# internal standards L1,H1,LW stay at 0 also
# in future maybe have a column for prepGas
setwd("~/Desktop/EuropaMLMS/microbe/QC_SU22")
pHionic.dat<-read.csv("pHionic_bioAbio.csv",header=T)
head(pHionic.dat)

length(unique(pHionic.dat$Identifier1))
# 116
length(unique(bioAbio.dat$Identifier1))
# 104

# duplicates in microbial mud, and there are 12 Identifier1s in pH data that
# are not in bioAbio?
inBiodat.ind<-which(unique(pHionic.dat$Identifier1)%in%unique(bioAbio.dat$Identifier1))
unique(pHionic.dat$Identifier1)[-inBiodat.ind]
# none of these passed QC? ... probably include in QC report and do from now on
#[1] "KCl_L_DT_300uL"               "MgCl2_L_DT_50uL"
#[3] "MgCl2_U_DT_100uL"             "NaHCO3_L_C_100uL"
#[5] "NaHCO3_L_C_50uL"              "KCl_U"
#[7] "MgCl2_L"                      "KCl_L_+_MgCl2_L"
#[9] "KCl_U_+_MgCl2_L"              "MgSO4_L_+_MgCl2_L"
#[11] "MgSO4_L_+_Na2SO4_L_+_MgCl2_U" "Na2SO4_L+_NaCl_L_+_NaHCO3_L"
#[13] "NaCl_L_+_KCl_L_+_MgCl2_L"     "NaCl_L_+_KCl_L_+_MgCl2_U"

which(bioAbio.dat$Identifier1=="NaCl_L_+_KCl_L_+_MgCl2_U")
# integer(0)

# there are some duplicates in the pH data also
pHionic.dat$Identifier1[which(duplicated(pHionic.dat$Identifier1))]
#[1] "MgSO4_L"            "MgSO4_L_+_NaCl_L"   "MgSO4_L_+_NaCl_U"   "MgSO4_L_+_NaHCO3_L"
#[5] "MgSO4_L_+_NaHCO3_U" "NaCl_L"             "NaCl_U"             "NaHCO3_L"
#[9] "NaHCO3_L_+_NaCl_L"  "NaHCO3_L_+_NaCl_U"  "NaHCO3_U"           "NaHCO3_U_+_NaCl_L"
#[13] "NaHCO3_U_+_NaCl_U"

# dataset of duplicates -- in microbial mud, and don't have vals for these for pH or ionic strength
#[1] "microbial_mud" "microbial_mud" "microbial_mud" "microbial_mud"
#[6] "microbial_mud" "microbial_mud" "microbial_mud" "microbial_mud" "microbial_mud"
#[11] "microbial_mud" "microbial_mud" "microbial_mud" "microbial_mud"




###############
bioAbio.dat<-read.csv("bioAbio_samps_ML_labels_2022-06-27.csv")
head(bioAbio.dat)

## fix Preparation labels
unique(bioAbio.dat$Preparation)
#[1] "He_only"
#[2] "CO2_in_He"
#[3] "2%_CO2_in_He_48hrs"
#[4] "2%_CO2_in_He_24hrs"
#[5] "0.3%_CO2_in_He_72_hr"
#[6] "1%_CO2_in_He_72hrs"
#[7] "0.115ml_CO2_in_He_24hrs"
#[8] "0.115ml_CO2_in_He_48hrs"
#[9] "Bacteria_Inoculation_0.3%_CO2_in_He_5_day"
#[10] "0.5%_CO2_in_He_5_day"


# all these are in the europa_bacteria dataset
bioAbio.dat$Preparation[which(bioAbio.dat$Preparation=="0.5%_CO2_in_He_5_day")]<-"0.5%_CO2_in_He_7_days"
bioAbio.dat$Preparation[which(bioAbio.dat$Preparation=="Bacteria_Inoculation_0.3%_CO2_in_He_5_day")]<-"Bacteria_Inoculation_0.3%_CO2_in_He_7_days"
# change the values in 0.115mL to be 1%
bioAbio.dat$Preparation[which(bioAbio.dat$Preparation=="0.115ml_CO2_in_He_48hrs")]<-"1%_CO2_in_He_2_days"
bioAbio.dat$Preparation[which(bioAbio.dat$Preparation=="0.115ml_CO2_in_He_24hrs")]<-"1%_CO2_in_He_1_day"

# change the values in Preparation to days if not done already -- use the same units!
# change back if Bethany doesn't like it, but Brett's code uses Preparation to build prepTime
bioAbio.dat$Preparation[which(bioAbio.dat$Preparation=="2%_CO2_in_He_48hrs")]<-"2%_CO2_in_He_2_days"
bioAbio.dat$Preparation[which(bioAbio.dat$Preparation=="2%_CO2_in_He_24hrs")]<-"2%_CO2_in_He_1_day"
bioAbio.dat$Preparation[which(bioAbio.dat$Preparation=="0.3%_CO2_in_He_72_hr")]<-"0.3%_CO2_in_He_3_days"
bioAbio.dat$Preparation[which(bioAbio.dat$Preparation=="1%_CO2_in_He_72hrs")]<-"1%_CO2_in_He_3_days"

unique(bioAbio.dat$Preparation)
#[1] "He_only"
#[2] "CO2_in_He"
#[3] "2%_CO2_in_He_2_days"
#[4] "2%_CO2_in_He_1_day"
#[5] "0.3%_CO2_in_He_3_days"
#[6] "1%_CO2_in_He_3_days"
#[7] "1%_CO2_in_He_1_day"
#[8] "1%_CO2_in_He_2_days"
#[9] "Bacteria_Inoculation_0.3%_CO2_in_He_7_days"
#[10] "0.5%_CO2_in_He_7_days"

# change "CO2_in_He" for microbial mud to 0.3pCO2
unique(bioAbio.dat$dataset[which(bioAbio.dat$Preparation=="CO2_in_He")])
# [1] "abiotic"       "microbial_mud"
# there are some abiotic Preparation labels without CO2 concentration information
#      -- will need to check with Bethany
missingCO2.ind<-which(bioAbio.dat$Preparation=="CO2_in_He")
unique(bioAbio.dat$Preparation[missingCO2.ind])
mm_missingCO2.ind<-which(bioAbio.dat$dataset[missingCO2.ind]=="microbial_mud")
unique(bioAbio.dat$dataset[missingCO2.ind][mm_missingCO2.ind])
bioAbio.dat$Preparation[missingCO2.ind][mm_missingCO2.ind]<-"0.3%_CO2_in_He_7_days"

unique(bioAbio.dat$Preparation)
#[1] "He_only"
#[2] "CO2_in_He"
#[3] "2%_CO2_in_He_2_days"
#[4] "2%_CO2_in_He_1_day"
#[5] "0.3%_CO2_in_He_3_days"
#[6] "1%_CO2_in_He_3_days"
#[7] "1%_CO2_in_He_1_day"
#[8] "1%_CO2_in_He_2_days"
#[9] "0.3%_CO2_in_He_7_days"
#[10] "Bacteria_Inoculation_0.3%_CO2_in_He_7_days"
#[11] "0.5%_CO2_in_He_7_days"

datasetMissingCO2concInd<-which(bioAbio.dat$Preparation=="CO2_in_He")
unique(bioAbio.dat$dataset[datasetMissingCO2concInd])
# [1] "abiotic"
unique(bioAbio.dat$fileId[datasetMissingCO2concInd])
#[1] "170503_H1__CO2_in_He.dxf"
#[2] "170504_H1__CO2_in_He.dxf"
#[3] "170503_L1__CO2_in_He.dxf"
#[4] "170503_L1_CO2_in_He.dxf"
#[5] "170504_L1__CO2_in_He.dxf"
#[6] "170503_LW__CO2_in_He.dxf"
#[7] "170503_LW_CO2_in_He.dxf"
#[8] "170504_LW__CO2_in_He.dxf"
#[9] "170503_MgSO4_L_+_NaCl_L__CO2_in_He.dxf"
#[10] "170503_MgSO4_L_+_NaCl_L__CO2_in_He(1).dxf"
#[11] "170504_MgSO4_L_+_NaCl_L__CO2_in_He.dxf"
#[12] "170503_MgSO4_L_+_NaCl_U__CO2_in_He.dxf"
#[13] "170503_MgSO4_L_+_NaCl_U__CO2_in_He(1).dxf"
#[14] "170504_MgSO4_L_+_NaCl_U__CO2_in_He.dxf"
#[15] "170504_MgSO4_L_+_NaHCO3_L__CO2_in_He.dxf"
#[16] "170503_MgSO4_L__CO2_in_He.dxf"
#[17] "170503_MgSO4_L__CO2_in_He(1).dxf"
#[18] "170503_NaCl_L__CO2_in_He.dxf"
#[19] "170503_NaCl_L_CO2_in_He.dxf"
#[20] "170504_NaCl_L__CO2_in_He.dxf"
#[21] "170503_NaCl_U__CO2_in_He.dxf"
#[22] "170503_NaCl_U_CO2_in_He.dxf"
#[23] "170504_NaCl_U__CO2_in_He.dxf"
#[24] "170503_NaHCO3_L_+_NaCl_L__CO2_in_He.dxf"
#[25] "170503_NaHCO3_L_+_NaCl_L_CO2_in_He.dxf"
#[26] "170504_NaHCO3_L_+_NaCl_L__CO2_in_He.dxf"
#[27] "170503_NaHCO3_L_+_NaCl_U_CO2_in_He.dxf"
#[28] "170504_NaHCO3_L_+_NaCl_U__CO2_in_He.dxf"
#[29] "170504_NaHCO3_L_+_NaCl_U__CO2_in_He(1).dxf"
#[30] "170503_NaHCO3_L__CO2_in_He.dxf"



# add gasPrep column and put He in all
bioAbio.dat$gasPrep<-rep("He",dim(bioAbio.dat)[1])
head(bioAbio.dat)

# give Brett data with "correct" column order
colnames(bioAbio.dat)
sampleMScolOrder<-c("PeakNr", "Start", "Rt", "End", "Ampl44", "Ampl45", "Ampl46",
                    "BGD44", "BGD45", "BGD46", "rIntensity44",
                    "rIntensity45", "rIntensity46", "rIntensityAll", "Intensity44",
                    "Intensity45", "Intensity46", "IntensityAll", "ListFirstPeak",
                    "rR45CO244CO2", "rR46CO244CO2", "IsRef", "R45CO244CO2",
                    "RefName", "rd45CO244CO2", "d45CO244CO2", "R46CO244CO2",
                    "rd46CO244CO2", "d46CO244CO2", "R13C12C", "d13C12C",
                    "AT13C12C", "R18O16O", "d18O16O", "AT18O16O",
                    "R17O16O", "d17O16O", "Rps45CO244CO2", "Rps46CO244CO2",
                    "pkArea", "d18O13C", "biotic", "avg_d13C12C",
                    "avg_d18O13C", "avg_d18O16O", "avg_Ampl46", "avg_IntensityAll",
                    "avg_rIntensityAll", "avg_pkArea", "ionicStrength", "pH",
                    "fileId","DateTime","dataset", "Identifier1", "Analysis","Preparation",
                    "prepTime", "CO2conc", "gasPrep", "rxn","saltContent","allAnion",
                    "allCation", "sulfate","chloride", "carbonate", "magnesium",
                    "sodium", "potassium", "anion", "cation")
length(sampleMScolOrder)
length(colnames(bioAbio.dat))

# write to file!
write.table(bioAbio.dat,"bioAbio_samps_ML_labels_2022-06-28.csv",row.names=F,quote=F,sep=",")


## extra Prep processing: keep 7 days and change the rest to hours
# all the CO2_in_He in abiotic dataset are 0.3%_CO2_in_He_7_days
unique(bioAbio.dat$Preparation)
#[1] "He_only"
#[2] "CO2_in_He"
#[3] "2%_CO2_in_He_2_days"
#[4] "2%_CO2_in_He_1_day"
#[5] "0.3%_CO2_in_He_3_days"
#[6] "1%_CO2_in_He_3_days"
#[7] "1%_CO2_in_He_1_day"
#[8] "1%_CO2_in_He_2_days"
#[9] "0.3%_CO2_in_He_7_days"
#[10] "Bacteria_Inoculation_0.3%_CO2_in_He_7_days"
#[11] "0.5%_CO2_in_He_7_days"

setwd("~/Desktop/EuropaMLMS/microbe/QC_SU22/bio_abio/bioAbio")

long.dat<-read.csv("long_data_withStandards.csv")
ms.dat<-read.csv("bioAbio_samps_ML_labels_2022-06-28.csv",header=T)

long.dat<-ms.dat
head(long.dat)

longDat.prep<-long.dat$Preparation
unique(longDat.prep)

CO2change<-which(longDat.prep=="CO2_in_He")
longDat.prep[CO2change]<-"0.3%_CO2_in_He_7_days"

p2dayChange<-which(longDat.prep=="2%_CO2_in_He_2_days")
longDat.prep[p2dayChange]<-"2%_CO2_in_He_48_hours"

p2day1<-which(longDat.prep=="2%_CO2_in_He_1_day")
longDat.prep[p2day1]<-"2%_CO2_in_He_24_hours"

change<-which(longDat.prep=="0.3%_CO2_in_He_3_days")
longDat.prep[change]<-"0.3%_CO2_in_He_72_hours"

change<-which(longDat.prep=="1%_CO2_in_He_3_days")
longDat.prep[change]<-"1%_CO2_in_He_72_hours"

change<-which(longDat.prep=="1%_CO2_in_He_1_day")
longDat.prep[change]<-"1%_CO2_in_He_24_hours"

change<-which(longDat.prep=="1%_CO2_in_He_2_days")
longDat.prep[change]<-"1%_CO2_in_He_48_hours"

change<-which(longDat.prep=="Bacteria_Inoculation_0.3%_CO2_in_He_7_days")
longDat.prep[change]<-"0.3%_CO2_in_He_7_days"

orig<-bioAbio.dat$Preparation
changed<-longDat.prep

prep.df<-as.data.frame(matrix(c(orig,changed),ncol=2))
colnames(prep.df)<-c("Prep_days","Prep_final_change")
head(prep.df)

long.dat$Preparation<-longDat.prep
unique(long.dat$Preparation)

#write.table(prep.df,"prepMapping.csv",row.names=F,quote=F,sep=",")

## change carbonate to bicarbonate
carbs<-long.dat$carbonate
which(colnames(long.dat)=="carbonate")
colnames(long.dat)[62]<-"bicarbonate"
carbInd<-which(carbs=="carbonate")
carbs[carbInd]<-"bicarbonate"
carbs[-carbInd]<-"no_bicarbonate"

unique(carbs)
long.dat$bicarbonate<-carbs

head(long.dat)
tail(long.dat)
unique(long.dat$Preparation)
long.dat$prepDays<-orig
unique(long.dat$prepDays)
head(long.dat)


BI.ind<-which(orig=="Bacteria_Inoculation_0.3%_CO2_in_He_7_days")
long.dat$prepDays[BI.ind]<-"0.3%_CO2_in_He_7_days"

CO2ind<-which(orig=="CO2_in_He")
long.dat$prepDays[CO2ind]<-"0.3%_CO2_in_He_7_days"

unique(long.dat$prepDays)
unique(long.dat$Preparation)

write.table(long.dat,"bioAbio_samp_ML_labels_newPrepBiCarb.csv",row.names=F,quote=F,sep=",")

unique(long.dat$anion)
unique(long.dat$cation)

##

