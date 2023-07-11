### functions for IRMS QC
#'@export
all_dxf_files<-function(){
  all_files<-list.files()
  dxfFiles<-c()
  for(i in seq(1:length(all_files))){
    currFile<-all_files[i] # get file name
    if(grepl(".dxf", currFile)){
      dxfFiles<-c(dxfFiles,currFile)
    }
  }
  return(dxfFiles)
}


analysisNumsEqualDXF<-function(refsList,sampsList){
  lenRefs<-length(refsList)
  refsAnNums<-c()
  for(i in seq(1,lenRefs)){
    refsAnNums<-c(refsAnNums,as.numeric(refsList[[i]]$Analysis[1]))
  }

  lenSamps<-length(sampsList)
  sampsAnNums<-c()
  for(i in seq(1,lenSamps)){
    sampsAnNums<-c(sampsAnNums,as.numeric(sampsList[[i]]$Analysis[1]))
  }
  # if they are equal
  lenEq<-(length(refsAnNums)==length(sampsAnNums))
  if(lenEq){
    analysisNumsEq<-sum(refsAnNums-sampsAnNums)==0
    if(analysisNumsEq){
      return(analysisNumsEq)
    }
  }else{ # if they are not equal, return analysis nums that are different?
    anNum.list<-list()
    anRef.df<-as.data.frame(matrix(refsAnNums,ncol=1))
    colnames(anRef.df)<-c("ref_an")
    anS.df<-as.data.frame(matrix(sampsAnNums,ncol=1))
    colnames(anS.df)<-c("samps_an")
    anNum.list[[1]]<-FALSE
    anNum.list[[2]]<-anRef.df
    anNum.list[[3]]<-anS.df
    return(anNum.list)
  }

}

# debugging
#allStandards_d18O.list<-standIsoR.list
# standNames<-internalStandID
#standAcceptedVals.vec<-c(-8.55,4.85,-3.85)
#accStandRatioSD<-c(0.2,0.2,0.2)
avgSD_d18O_standards<-function(allStandards_d18O.list,standNames,standAcceptedVals.vec,accStandRatioSD){
  # TODO: add isoRatioNames as arg to make FileDf and other labels
  # isoRatioNames<-c("d18O16O","d13C12C")
  # debugging
  #allStandards_d18O.list<-subStandIsoR.list
  #standNames=c("H1","LW")
  #standAcceptedVals.vec=c(-8.55,-3.85)
  #accStandRatioSD=c(0.2,0.2)

  # make df for avgs (avgs of all d18O in all files for particular standard)
  standAccepted.mat<-matrix(standAcceptedVals.vec,ncol=1)
  standAccepted.df<-as.data.frame(standAccepted.mat)
  # start df -- will add average of all measured values
  # TODO: add d13C
  d18OstandVals_Avgs.df<-data.frame(matrix(rep(NA,dim(standAccepted.df)[1]*2),ncol=2))
  colnames(d18OstandVals_Avgs.df)<-c("acceptAvg_d18O16O","measAvg_d18O16O")
  rownames(d18OstandVals_Avgs.df)<-standNames
  d18OstandVals_Avgs.df[,1]<-standAccepted.df
  d18OstandVals_Avgs.df
  # d13C
  #d13CstandVals_Avgs.df<-data.frame(matrix(rep(NA,dim(standAccepted.df)[1]*2),ncol=2))
  #colnames(d13CstandVals_Avgs.df)<-c("acceptAvg_d13C12C","measAvg_d13C12C")
  #rownames(d13CstandVals_Avgs.df)<-standNames
  #d13CstandVals_Avgs.df[,1]<-standAccepted.df

  # make df for SDs
  SDAccepted.mat<-matrix(accStandRatioSD,ncol=1)
  SDAccepted.df<-as.data.frame(SDAccepted.mat)
  # start df -- will add SD of all measured values
  d18OstandVals_SDs.df<-data.frame(matrix(rep(NA,dim(SDAccepted.df)[1]*2),ncol=2))
  colnames(d18OstandVals_SDs.df)<-c("accept_sd_d18O16O","calc_sd_d18O16O")
  rownames(d18OstandVals_SDs.df)<-standNames
  d18OstandVals_SDs.df[,1]<-accStandRatioSD
  d18OstandVals_SDs.df

  # initialize lists for return vals
  ret_d18.list<-list()
  d18O.list<-allStandards_d18O.list
  list.len<-length(d18O.list)
  alld18OFile.list<-list()
  alld18OFile.index<-1

  d18OFile.mat<-matrix(rep(NA,list.len*4),ncol=4)
  d18OFile.df<-as.data.frame(d18OFile.mat)
  colnames(d18OFile.df)<-c("Analysis","rxn","avg_d18O16O","sd_d18O16O")
  d18OFile.df
  #**
  #d13CFile.mat<-matrix(rep(NA,list.len*4),ncol=4)
  #d13CFile.df<-as.data.frame(d13CFile.mat)
  #colnames(d13CFile.df)<-c("Analysis","Identifier1","avg_d13C12C","sd_d13C")
  #d13CFile.df

  stand.list<-list()
  stand.vec<-c(rep(NA,2*length(standNames)))

  for(i in seq(1,length(stand.vec))){
    stand.list[[i]]<-stand.vec[i]
    #stand.list[[i+1]]<-stand.vec[i]
  }
  stand.list

  # j<-1
  # for(i in seq(1,length(stand.list),by=2)){
  #   currStandName<-standNames[j]
  #   currStandName
  #   names(stand.list)[i]<-c(paste(standNames[j],"_tot_d18O16O",sep=""))
  #   #names(stand.list)[i+1]<-c(paste(standNames[j],"tot_d13C12C",sep=""))
  #   stand.list
  #   j<-j+1
  # }

  names(stand.list)<-c("L1_tot_d18O16O","L1_tot_d13C12C",
                       "H1_tot_d18O16O","H1_tot_d13C12C",
                       "LW_tot_d18O16O","LW_tot_d13C12C")
  stand.list

  d18Otot.vec<-c()
  #d13Ctot.vec<-c()
  for(i in seq(1,length(allStandards_d18O.list))){
    # for each file in standard set
    dlist<-(d18O.list[[i]])
    dlist
    # get sample peaks d18O/16O for one file
    dlist.df<-as.data.frame(dlist)
    dlist.df
    # get the analysis number
    anNum<-dlist.df$Analysis[1]
    anNum
    # get Identifier1
    Id1<-dlist.df$Identifier1[1]
    #Id1<-dlist.df$rxn[1]
    Id1
    # get the d18O data for the file
    d18O.vec<-as.numeric(dlist.df$d18O16O)
    d18O.vec
    #d13C.vec<-as.numeric(dlist.df$d13C12)
    #d13C.vec
    # add to total vecs
    # find which standard, add to that vec

    ind<-which(grepl(Id1,names(stand.list),fixed=T))
    ind
    # there should be two for each Identifier1
    firstInd<-ind[1]
    firstInd
    #secInd<-ind[2]
    #secInd

    if(length(stand.list[[firstInd]])<=1){
      d18Otot.vec<-c(d18O.vec)
      d18Otot.vec
      #d13Ctot.vec<-c(d13C.vec)
      #d13Ctot.vec
    }else{
      d18Otot.vec<-c(stand.list[[firstInd]],d18O.vec)
      d18Otot.vec
      #d13Ctot.vec<-c(stand.list[[secInd]],d13C.vec)
      #d13Ctot.vec
    }
    stand.list[[firstInd]]<-d18Otot.vec
    #stand.list[[secInd]]<-d13Ctot.vec
    stand.list

    # get avg d18O for each file
    avg_d18O<-mean(d18O.vec)
    avg_d18O
    # get sd for each file
    sd_d18O<-sd(d18O.vec)
    sd_d18O

    # d13C
    #avg_d13C<-mean(d13C.vec)
    #avg_d13C
    #sd
    #sd_d13C<-sd(d13C.vec)
    #sd_d13C

    # TODO: put file_id, avg, and sd for each file in df
    # TODO: makeDeltaFile.df and add cols for d13C avg and sd
    d18OFile.df[i,]<-data.frame(anNum,Id1,avg_d18O,sd_d18O)
    d18OFile.df
    #d13CFile.df[i,]<-data.frame(anNum,Id1,avg_d13C,sd_d13C)
    #d13CFile.df

  }

  # find the avgs and sd of the totals
  # want to do d13C and d18O separately!
  # d18O/16O
  tot_d18O.df<-as.data.frame(matrix(rep(NA,2*length(standNames)),ncol=2))
  rownames(tot_d18O.df)<-standNames
  colnames(tot_d18O.df)<-c("avg_d18O16O","sd_d18O16O") #,"avg_d13C12C","sd_d13C12C"
  tot_d18O.df
  # d13C/12C
  #tot_d13C.df<-as.data.frame(matrix(rep(NA,2*length(standNames)),ncol=2))
  #rownames(tot_d13C.df)<-standNames
  #colnames(tot_d13C.df)<-c("avg_d13C12C","sd_d13C12C")
  #tot_d13C.df
  # find avgs and sds and add to df
  for(n in seq(1,dim(tot_d18O.df)[1])){
    stand.ind<-which(grepl(rownames(tot_d18O.df)[n],names(stand.list)))
    stand.ind
    stand_d18O.vals<-stand.list[[stand.ind[1]]]
    #stand_d13C.vals<-stand.list[[stand.ind[2]]]

    d18O.avg<-mean(stand_d18O.vals)
    d18O.sd<-sd(stand_d18O.vals)

    #d13C.avg<-mean(stand_d13C.vals)
    #d13C.sd<-sd(stand_d13C.vals)

    tot_d18O.df[n,]<-c(d18O.avg,d18O.sd)#,d13C.avg,d13C.sd
    #tot_d13C.df[n,]<-c(d13C.avg,d13C.sd)
  }

  tot_d18O.df
  #tot_d13C.df


  #*****
  # build d18O/16O return data
  d18OstandVals_Avgs.df$measAvg_d18O16O<-tot_d18O.df$avg_d18O16O
  d18OstandVals_SDs.df$calc_sd_d18O16O<-tot_d18O.df$sd_d18O16O

  # d13C/12C
  #d13CstandVals_Ags.df
  #d13CstandVals_SDs.df

  # sep by identifier1
  numGroups<-length(standNames)
  IDGroup.list<-list()

  for(i in seq(1,numGroups)){
    IDind<-which(d18OFile.df$rxn==standNames[i]) #****
    # get that data
    IDdat<-d18OFile.df[IDind,]
    IDGroup.list[[i]]<-IDdat
  }

  # add values to return list
  ret_d18.list[[1]]<-IDGroup.list
  ret_d18.list[[2]]<-d18OstandVals_Avgs.df # all avgs
  ret_d18.list[[3]]<-d18OstandVals_SDs.df # all SDs
  ret_d18.list[[4]]<-tot_d18O.df # total avg and sd for d18O/16O and d13C/12C

  return(ret_d18.list)
}

# debugging
#path<-unsortedPath
#combColNames
#outputID<-dataName
combineVendFileInfo<-function(path,combColNames,outputID){
  # make a combined dataframe for vendor and file info
  files<-all_dxf_files()
  vend<-vendor_info_all(files)
  colnames(vend[[1]])


  fileInfo<-file_info(files)
  #head(fileInfo)
  dim(fileInfo)
  # make a list of combined vendor data and file info dfs for each experiment
  combList<-list()
  failedInd<-c()
  failedFiles<-c()
  for(i in seq(1,length(files))){
    # make file info df to add to vend df
    mat<-matrix(rep(NA,nrow(vend[[i]])*ncol(fileInfo)),#ncol(fileInfo[[i]])),
                nrow=nrow(vend[[i]]))
    df<-as.data.frame(mat)
    df[1:nrow(df),]<-fileInfo[i,]
    # need to check if there is data
    isVendData<-!(dim(vend[[i]])[1]==0)
    if(isVendData){
      # make a combined df
      combDF<-data.frame(cbind(df,vend[[i]][,2:ncol(vend[[i]])]))
      colnames(combDF)
      # fileId, Identifier1, Analysis, DateTime, PeakNr,
      # Start, Rt, End, Ampl44,
      combColNames[1:10]
      colnames(combDF)<-combColNames
      # add to list
      combList[[i]]<-combDF
    }else{
      fileName<-files[i]
      print(paste("No data found in file: ",fileName,". File name inserted at index ",i,sep=""))
      failedFiles<-c(failedFiles,fileName)
      failedInd<-c(failedInd,i)
      combList[[i]]<-fileName
    }
  }
  # if there is missing data remove that file
  if(length(failedInd)>0){
    print("Removing files with no data from the analysis...")
    onlyDat<-combList[c(-failedInd)]
    failedFiles.df<-as.data.frame(matrix(failedFiles),ncol=1)
    colnames(failedFiles.df)<-c("No_Data")
    write.table(failedFiles.df,file=paste(outputID,"_NoData.txt",sep=""),row.names=F,quote=F)
    return(onlyDat)
  }else{
    return(combList)
  }
}


DXFvendListFromID1<-function(all.df,standName.vec){#=c("L1","H1","LW")
  vendRet.list<-list()
  listIndex<-1
  # separate into different experiments using analysis numbers
  sepList<-separate_by_analysis_numDXF(all.df)
  # look just for identifier1 that are in standName.vec
  lenList<-length(sepList)
  for(i in seq(1,lenList)){
    currID<-sepList[[i]]$Identifier1[1]
    currID
    foundStand<-currID %in% standName.vec
    foundStand
    if(foundStand){ # get vend data and add to output list: Analysis,Identifier1,PkNr,d18O16O,d13C/12C
      standVend<-sepList[[i]]
      currAn.vec<-standVend$Analysis
      currID.vec<-standVend$Identifier1
      currPkNr.vec<-standVend$PeakNr
      currd18O<-standVend$d18O16O
      currd13C<-standVend$d13C12C
      curr.df<-as.data.frame(matrix(c(currAn.vec,currID.vec,currPkNr.vec,currd18O,currd13C),ncol=5))
      colnames(curr.df)<-c("Analysis","Identifier1","PkNr","d18O16O","d13C12C")
      # add to list
      vendRet.list[[listIndex]]<-curr.df
      listIndex<-listIndex+1
    }
  }
  return(vendRet.list)
}


file_info<-function(files){
  num_files<-length(files)
  msdat<-iso_read_continuous_flow(files[1:num_files])
  file_info<-msdat %>%
    iso_get_file_info(
      select = c(
        #rename?
        Identifier1 = `Identifier 1`,
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
  file_info.df
  return(file_info.df)
}


file_info_all<-function(files){
  file.list<-list()
  for(i in seq(1,length(files))){
    file.df<-file_info(files[i])
    file.list[[i]]<-file.df
    #print(i)
  }
  file.list[[i]]
  return(file.list)
}


fix_spaces_fileInfo<-function(path){
  ret.list<-list() # new file names, new dataframe with correct labels
  filenames<-list.files(path)
  #head(filenames)

  # need to get only dxf files -- TODO: make into separate function!
  dxfFiles.logical<-grepl(".dxf",filenames,fixed = TRUE)
  dxfFile.ind<-which(dxfFiles.logical==TRUE)
  dxfFiles<-filenames[dxfFile.ind]
  # store files in return list
  ret.list[[1]]<-dxfFiles

  # read in dxf data
  iso.dat<-iso_read_continuous_flow(path)
  #head(iso.dat)
  #class(iso.dat) # continuous flow

  # get file info
  iso.info<-iso_get_file_info(iso.dat)
  #class(iso.info)
  isoInfo.df<-as.data.frame(iso.info)
  #head(isoInfo.df)

  # get original labels and colnames
  identifier1.orig<-isoInfo.df$`Identifier 1`
  colnames.orig<-colnames(isoInfo.df)
  prep.orig<-isoInfo.df$Preparation

  # set up new labels vecs
  newID1.vec<-identifier1.orig
  newColnames.vec<-colnames.orig
  ### testing code
  #testCol.vec<-c("test  1", "test 2", "test3")
  #colnames.orig<-newColnames.vec
  newPrep.vec<-prep.orig

  # use grepl for string matching - check for two spaces :(
  # colnames
  # check for two and fix first, then check for one because if two spaces
  # are found in a row, then the two spaces will be replaced with two underscores
  twoSpacesCols<-grepl("  ",colnames.orig,fixed=T)
  twoSpacesCols
  if(sum(twoSpacesCols)>0){
    # find, fix, add to vec
    containsTwoSpacesCols.ind<-which(twoSpacesCols==T)
    newNames2.vec<-gsub("  ","_",colnames.orig[containsTwoSpacesCols.ind])
    newColnames.vec[containsTwoSpacesCols.ind]<-newNames2.vec
  }
  # now check for one space
  oneSpaceCols<-grepl(" ",newColnames.vec,fixed=T)
  oneSpaceCols
  if(sum(oneSpaceCols) > 0){
    # find, fix, add to vec
    containsSpacesCols.ind<-which(oneSpaceCols==T)
    newNames1.vec<-gsub(" ","_",newColnames.vec[containsSpacesCols.ind])
    newColnames.vec[containsSpacesCols.ind]<-newNames1.vec
    newColnames.vec
  }
  # add new colnames to data
  colnames(isoInfo.df)<-newColnames.vec

  # fix identifier1 labels
  # again, check for two spaces first and fix, then one
  twoSpacesID1<-grepl("  ",identifier1.orig)
  twoSpacesID1
  if(sum(twoSpacesID1)>0){
    # find, fix, add to vec
    contains2SpacesID1.ind<-which(twoSpacesID1==T)
    newID1_2.vec<-gsub("  ","_",identifier1.orig[contains2SpacesID1.ind])
    newID1.vec[contains2SpacesID1.ind]<-newID1_2.vec
  }
  # now check for one space
  oneSpaceID1<-grepl(" ",newID1.vec,fixed=T)
  oneSpaceID1
  if(sum(oneSpaceID1) > 0){
    # find, fix, add to vec
    containsSpacesID1.ind<-which(oneSpaceID1==T)
    newID1_1.vec<-gsub(" ","_",newID1.vec[containsSpacesID1.ind])
    newID1.vec[containsSpacesID1.ind]<-newID1_1.vec
    newID1.vec
  }
  # replace identifier1 one in df with new ved
  isoInfo.df$Identifier_1<-newID1.vec

  # preparation
  unique(isoInfo.df$Preparation)
  twoSpacesPrep<-grepl("  ",prep.orig)
  twoSpacesPrep
  if(sum(twoSpacesPrep)>0){
    # find, fix, add to vec
    contains2SpacesPrep.ind<-which(twoSpacesPrep==T)
    newPrep2.vec<-gsub("  ","_",prep.orig[contains2SpacesPrep.ind])
    newPrep.vec[contains2SpacesPrep.ind]<-newPrep2.vec
  }
  # now check for one space
  oneSpacePrep<-grepl(" ",newPrep.vec,fixed=T)
  oneSpacePrep
  if(sum(oneSpacePrep) > 0){
    # find, fix, add to vec
    containsSpacesPrep.ind<-which(oneSpacePrep==T)
    newPrep1.vec<-gsub(" ","_",newPrep.vec[containsSpacesPrep.ind])
    newPrep.vec[containsSpacesPrep.ind]<-newPrep1.vec
    newPrep.vec
  }
  # add to df
  isoInfo.df$Preparation<-newPrep.vec
  head(isoInfo.df)

  # make df of file_id, Identifier1, Analysis, Preparation
  ret.df<-cbind.data.frame(isoInfo.df$file_id,isoInfo.df$Identifier_1,
                           isoInfo.df$Analysis, isoInfo.df$Preparation
  )
  colnames(ret.df)<-c("file_id","Identifier1","Analysis","Preparation")
  head(ret.df)
  ret.list[[2]]<-ret.df # return smaller df
  ret.list[[3]]<-isoInfo.df # return full fileInfo df
  names(ret.list)<-c("file_id","noSpaceInfo","fullInfo")
  return(ret.list)
}
intensity_similarityDXF<-function(vendAmpl,amplName,peakNr.vec,
                                  relDiffInt.thresh=0.1,ID1){
  intCheck.list<-list()
  # analyze intensity similarity by relative difference
  int.check<-FALSE
  relDiff<-c()
  rel.vec<-c()
  for(i in seq(2,length(vendAmpl))){
    curr.relDiff<-abs(vendAmpl[i]-vendAmpl[i-1])/max(vendAmpl[i],vendAmpl[i-1])
    rel.vec<-c(rel.vec,curr.relDiff)
    #print(curr.relDiff)
    if(curr.relDiff<relDiffInt.thresh){
      relDiff<-c(relDiff,TRUE)
    }else{
      relDiff<-c(relDiff,FALSE)
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
  intCheck.mat<-matrix(c(peakNr.vec,vendAmpl),ncol=2)
  intCheck.df<-as.data.frame(intCheck.mat)
  colnames(intCheck.df)<-c("PeakNr",amplName)
  # add to list
  intCheck.list[[2]]<-intCheck.df
  intCheck.list[[3]]<-max(rel.vec)
  return(intCheck.list)
}

# refSampIntCheck<-intensity_similarityDXF(vendAmpl=refSampAmpl.df$Ampl44,
#                                          amplName="Ampl44",
#                                          peakNr.vec=refSampAmpl.df$PeakNr,
#                                          relDiffInt.thresh=0.75)


isoR_similarityDXF<-function(vend.df,peakNr.vec,sdC.thresh=0.1,sdO.thresh=0.1,ID1){
  isoR.list<-list()
  d13C<-vend.df$d13C12C
  d13C.ref<-as.numeric(d13C[peakNr.vec])
  d18O<-vend.df$d18O16O
  d18O.ref<-as.numeric(d18O[peakNr.vec])
  analysis<-vend.df$Analysis
  #ID1<-vend.df$Identifier1
  # sd of d13C/12C
  sd.13C<-sd(d13C.ref)
  sd.C<-FALSE
  sdnac<-is.na(sd.13C)
  if(!sdnac){
    if(sd.13C<sdC.thresh){
      sd.C<-TRUE
      #print("SD 13C/12C within accepted threshold")
    } else{
      print(paste("SD 13C/12C not within accepted threshold: ",
                  round(sd.13C,2),sep=""))
      sd.18O<-sd(d18O.ref)
      print(paste("SD 18O/16O: ",
                  round(sd.18O,2),sep=""))
      print(paste("Analysis: ",analysis[1],sep=""))
      print(paste("ID1=",ID1[1],sep=""))
      print(paste("num peaks = ",length(peakNr.vec),sep=""))
    }
  } else{
    print("SD 13C/12C = NA")
  }

  # sd of 18O/16O
  sd.18O<-sd(d18O.ref)
  sd.O<-FALSE
  isnasdo<-is.na(sd.18O)
  if(!isnasdo){
    if(sd.18O<sdO.thresh){
      sd.O<-TRUE
      #print("SD 18O/16O within accepted threshold")
    } else{
      print(paste("SD 18O/16O not within accepted threshold: ",round(sd.18O,2),sep=""))
      print(paste("SD 13C/12C: ",
                  round(sd.13C,2),sep=""))
      print(paste("Analysis: ",analysis[1]))
      print(paste("ID1: ",ID1[1],sep=""))
      print(paste("num peaks = ",length(peakNr.vec),sep=""))
    }
  }else{
    print("SD 18O/16O = NA")
  }

  # dataframe for boolean vals for within threshold
  sdBool.mat<-matrix(c(sd.C,sd.O),ncol=2)
  sdBool.df<-as.data.frame(sdBool.mat)
  colnames(sdBool.df)<-c("d13C12C","d18O16O")
  isoR.list[[1]]<-sdBool.df

  # dataframe for numeric value of sds
  sdnum.mat<-matrix(c(sd.13C,sd.18O),ncol=2)
  sdnum.df<-as.data.frame(sdnum.mat)
  colnames(sdnum.df)<-c("SD_d13C12C","SD_d18O16O")
  isoR.list[[2]]<-sdnum.df

  return(isoR.list)
}


listToDF<-function(dat.list){
  # thread lists into dataframes
  # need total number of elements in the lists
  # sample list
  numSamp<-length(dat.list)
  sampSum<-0
  for(i in seq(1,numSamp)){
    sampIndexLength<-dim(dat.list[[i]])[1]
    sampSum<-sampSum+sampIndexLength
  }
  # initalize df
  numSampCols<-length(colnames(dat.list[[1]]))
  samp.df<-as.data.frame(matrix(rep(NA,numSampCols*sampSum),ncol=numSampCols))
  colnames(samp.df)<-colnames(dat.list[[1]])
  # loop through list and build dfs
  sampIndex<-1
  for(i in seq(1,numSamp)){
    # grab ref vend
    sampDat<-dat.list[[i]]
    numSampPts<-length(sampDat$PeakNr)
    # set indices
    sampInd.vec<-seq(sampIndex,(sampIndex+numSampPts-1))
    # add to df
    samp.df[sampInd.vec,]<-sampDat
    # update ref index
    sampIndex<-sampInd.vec[length(sampInd.vec)]+1
  }
  return(samp.df)
}


lm_graph_eqn <- function(lm){
  r2<-summary(lm)$r.squared
  intercept<-lm$coefficients[1]
  slope<-lm$coefficients[2]
  eqTitle<-paste("\n y = ",round(slope,3)," x + ",round(intercept,3),"\n r^2 = ",round(r2,4),sep="")
  return(eqTitle)
}


maxNumPeaksDXF<-function(vend.df,maxExpectedPks=18){

  maxRet.list<-list()
  maxRet.df<-as.data.frame(matrix(rep(NA,2),ncol=2))
  colnames(maxRet.df)<-c("NumPks","Analysis")
  num_peaks<-length(vend.df$PeakNr)

  # check if data has more than the max expected number of peaks
  weirdNum<-( (num_peaks-7 != 9) && (num_peaks-7 != 10) && (num_peaks-7 != 11))
  if(weirdNum){
    print(paste("weird num peaks Analysis ",vend.df$Analysis[1],sep=""))
  }
  if((num_peaks<=maxExpectedPks) && (num_peaks>0) && (!weirdNum)){
    withinMaxPkNum<-TRUE
    maxRet.list[[1]]<-withinMaxPkNum
    anNum<-vend.df$Analysis[1]
    maxRet.df[1,]<-c(num_peaks,anNum)
    maxRet.list[[2]]<-maxRet.df
    #return(maxRet.list)
  }else{
    if(num_peaks==0){
      print(paste("detected 0 peaks for Analysis ",vend.df$Analysis[1],sep=""))
    }
    withinMaxPkNum<-FALSE
    maxRet.list[[1]]<-withinMaxPkNum
    failedAnNum<-vend.df$Analysis[1]
    print(paste("Anomalous number of peaks detected:",num_peaks,"peaks in Analysis",failedAnNum))
    maxRet.df[1,]<-c(num_peaks,failedAnNum)
    maxRet.list[[2]]<-maxRet.df
  }
  return(maxRet.list)
}


reference_times_checkDXF<-function(vend.df,expectedPeak.num=5, diff.t=10,expectedStart,expectedRt,expectedEnd){
  refs.list<-list()
  peak.nums<-as.numeric(vend.df$PeakNr)
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
      #print(paste("Start, Rt, and End times match an expected peak for Peak_Nr ",peak.nums[i],sep=""))
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
    #print("Expected number of peaks detected at expected times.")
    refs.list[[1]]<-numRefs
  } else{
    print("Did not detected the expected number of peaks at expected times")
    refs.list[[1]]<-numRefs
  }
  ret<-cbind(vend.df$PeakNr[ref.ind],start.times[ref.ind],Rts[ref.ind],end.times[ref.ind])
  ret.df<-as.data.frame(ret)
  colnames(ret.df)<-c("PkNr","Start","Rt","End")
  refs.list[[2]]<-ret.df
  # vendor data
  refs.list[[3]]<-vend.df[ref.ind,]
  return(refs.list)
}


removeFailedAnalysesDXF<-function(sepList,expRef.df,maxPkNum=18,expRefPkNr=5,sdCrefIso.thresh=0.1,
                                  sdOrefIso.thresh=0.1,relDiffInt.thresh=0.1,
                                  sdCsampIso.thresh=0.3,sdOsampIso.thresh=0.2){
  sep.list<-list()
  samps.list<-list()
  refs.list<-list()
  # get reference expected times
  expRefStart<-expRef.df$Expected_Start
  expRefRt<-expRef.df$Expected_Rt
  expRefEnd<-expRef.df$Expected_End
  # set up loop
  numCheck<-length(sepList)
  # initialize list indices
  refIndex<-1
  sampIndex<-1
  peakIndex<-1
  # intialize vecs for each check-analyis nums and identifier1
  # peak number check-0 peaks and over max
  failedPkNum.vec<-c()
  failedPkID.vec<-c()
  failedNumPks.vec<-c()
  # reference times
  failedRefAnNum.vec<-c()
  failedRefID.vec<-c()
  failedRefDat.list<-list()
  # reference isotope ratios
  failedRefIso.vec<-c()
  failedRefIsoID.vec<-c()
  failedRefIsoSD.vec<-c() #d18O/16O
  failedRefIsoSDC.vec<-c()
  failedRefIsoReason.vec<-c()
  # reference intensity similiarity
  failedRefInt.vec<-c()
  failedRefIntID.vec<-c()
  failedRefIntRD.vec<-c()
  # sample processing
  failedSampAnNum.vec<-c()
  failedSampID.vec<-c()
  failedSampReas.vec<-c()
  # sample isotope ratio check
  failedSampIso.vec<-c()
  failedSampIsoID.vec<-c()
  failedSampIsoSD.vec<-c()
  failedSampIsoSDC.vec<-c()
  failedSampIsoReason.vec<-c()

  for(i in seq(1,numCheck)){
    #if(sepList[[i]]$Analysis[1]=="5772"){
    #  print(i)
    #}
    # check number of peaks
    numPeaksCheck<-maxNumPeaksDXF(sepList[[i]],maxExpectedPks=maxPkNum)
    #numPeaksCheck<-maxNumPeaksDXF(vend.df)
    numPks<-numPeaksCheck[[2]]$NumPks
    numPksPass<-numPeaksCheck[[1]]
    if(!numPksPass){ # failed num peak check
      failedPkNum.vec<-c(failedPkNum.vec,sepList[[i]]$Analysis[1])
      #failedPkNum.vec<-c(eurSep[[8]]$Analysis[1])
      failedPkID.vec<-c(failedPkID.vec,sepList[[i]]$Identifier1[1])
      failedNumPks.vec<-c(failedNumPks.vec,numPks)
    } else{ # add to data list for further checks
      sep.list[[peakIndex]]<-sepList[[i]]
      peakIndex<-peakIndex+1
    }
  }
  # update numCheck
  numCheck<-length(sep.list)
  if(is.null(failedPkNum.vec)){
    print("no analyses failed max peak number check")
    failedPkNum.vec<-c(NA)
    failedPkID.vec<-c(NA)
    failedNumPks.vec<-c(NA)
  }
  #*
  if(sum(grepl(" ", failedPkID.vec))>0){
    failedPkID.vec<-gsub(" ","_",failedPkID.vec)
  }
  #*
  failedPkNum.df<-as.data.frame(matrix(c(failedPkNum.vec,failedPkID.vec,failedNumPks.vec),ncol=3))
  colnames(failedPkNum.df)<-c("failed_PkNr_Analysis","failed_PkNr_Identifier1","failed_PkNr_NumPeaks")
  ## reference and sample peak checks
  failedRefDatIndex<-1
  for(i in seq(1,numCheck)){
    sampRefCheck<-reference_times_checkDXF(sep.list[[i]],expectedPeak.num=expRefPkNr,expectedStart=expRefStart,expectedRt=expRefRt,expectedEnd=expRefEnd)
    if(!sampRefCheck[[1]]){ # if failed ref time check
      failedRefAnNum.vec<-c(failedRefAnNum.vec,sep.list[[i]]$Analysis[1])
      failedRefID.vec<-c(failedRefID.vec,sep.list[[i]]$Identifier1[1])
      #*
      if(sum(grepl(" ", failedRefID.vec))>0){
        failedRefID.vec<-gsub(" ","_",failedRefID.vec)
      }
      #*
      failedRefDat.list[[failedRefDatIndex]]<-sampRefCheck[[2]]
      failedRefDatIndex<-failedRefDatIndex+1
    } else{
      # add ref iso ratio check and intensity similarity check
      # need pknr.vec for iso_ratio_similarity
      refPkNr.vec<-sampRefCheck[[2]]$PkNr
      refIsoCheck<-isoR_similarityDXF(vend.df=sep.list[[i]],peakNr.vec=refPkNr.vec,
                                      sdC.thresh=sdCrefIso.thresh,sdO.thresh=sdOrefIso.thresh,
                                      ID1=sep.list[[i]]$Identifier1[1])
      if((!refIsoCheck[[1]]$d18O16O)|(!refIsoCheck[[1]]$d13C12C)){
        refReasInd.vec<-which(refIsoCheck[[1]][1,]==FALSE)
        refReason<-colnames(refIsoCheck[[1]])[refReasInd.vec]
        if(length(refReason)==2){
          # both isotope SDs failed
          refReason<-c("d13C12C_d18O16O")
        }
        failedRefIso.vec<-c(failedRefIso.vec,sep.list[[i]]$Analysis[1])
        failedRefIsoID.vec<-c(failedRefIsoID.vec,sep.list[[i]]$Identifier1[1])
        #*
        if(sum(grepl(" ", failedRefIsoID.vec)>0)){
          failedRefIsoID.vec<-gsub(" ","_",failedRefIsoID.vec)
        }
        #*
        failedRefIsoSD.vec<-c(failedRefIsoSD.vec,refIsoCheck[[2]]$SD_d18O16O)
        failedRefIsoSDC.vec<-c(failedRefIsoSDC.vec,refIsoCheck[[2]]$SD_d13C12C)
        failedRefIsoReason.vec<-c(failedRefIsoReason.vec,refReason)
      } else{ # if it passes check intensity similarities
        refIntCheck<-intensity_similarityDXF(sep.list[[i]]$Ampl44[refPkNr.vec],
                                             amplName="Ampl44",refPkNr.vec,
                                             relDiffInt.thresh=relDiffInt.thresh,
                                             ID1=sep.list[[i]]$Identifier1[1])
        #maxRefIntRelDiff<-refIntCheck[[3]] #**
        if(!refIntCheck[[1]]){
          failedRefInt.vec<-c(failedRefInt.vec,sep.list[[i]]$Analysis[1])
          failedRefIntID.vec<-c(failedRefIntID.vec,sep.list[[i]]$Identifier1[1])
          #*
          if(sum(grepl(" ", failedRefIntID.vec)>0)){
            failedRefIntID.vec<-gsub(" ","_",failedRefIntID.vec)
          }
          #*
          failedRefIntRD.vec<-c(failedRefIntRD.vec,refIntCheck[[3]])
        }else{# only if passes all above 3 ref checks add to ref data
          refs.list[[refIndex]]<-sampRefCheck[[3]]
          refIndex<-refIndex+1
          # ** only check samples if ref data was added
          # process the sample peaks: remove flush peak and 1st sample peak (contamination)
          sampProcess<-sample_peaks_processDXF(sampRefCheck,sep.list[[i]])
          if(is.null(sampProcess)){
            failedSampAnNum.vec<-c(failedSampAnNum.vec,sep.list[[i]]$Analysis[1])
            failedSampID.vec<-c(failedSampID.vec,sep.list[[i]]$Identifier1[1])
            #*
            if(sum(grepl(" ", failedSampID.vec)>0)){
              failedSampID.vec<-gsub(" ","_",failedSampID.vec)
            }
            #*
            failedSampReas.vec<-c(failedSampReas.vec,"NoSamps")
          }else if(length(sampProcess$PeakNr)==0){
            print(paste("No sample peaks detected after processing for Analysis ",sep.list[[i]]$Analysis[1] ))
            # add to remove vector
            failedSampAnNum.vec<-c(failedSampAnNum.vec,sep.list[[i]]$Analysis[1])
            failedSampID.vec<-c(failedSampID.vec,sep.list[[i]]$Identifier1[1])
            #*
            if(sum(grepl(" ", failedSampID.vec)>0)){
              failedSampID.vec<-gsub(" ","_",failedSampID.vec)
            }
            #*
            failedSampReas.vec<-c(failedSampReas.vec,"NoSamps_afterProc")
          }else{
            # check sample iso ratios for d18O/16O
            sampIsoCheck<-isoR_similarityDXF(vend.df=sep.list[[i]],peakNr.vec=sampProcess$PeakNr,
                                             sdC.thresh=sdCsampIso.thresh,sdO.thresh=sdOsampIso.thresh,
                                             ID1=sep.list[[i]]$Identifier1[1])
            if((!sampIsoCheck[[1]]$d13C12C)|(!sampIsoCheck[[1]]$d18O16O)){
              sampReasInd.vec<-which(sampIsoCheck[[1]][1,]==FALSE)
              sampReason<-colnames(sampIsoCheck[[1]])[sampReasInd.vec]
              if(length(sampReason)==2){
                # both isotope SDs failed
                sampReason<-c("d13C12C_d18O16O")
              }
              failedSampIso.vec<-c(failedSampIso.vec,sep.list[[i]]$Analysis[1])
              failedSampIsoID.vec<-c(failedSampIsoID.vec,sep.list[[i]]$Identifier1[1])
              #*
              if(sum(grepl(" ", failedSampIsoID.vec)>0)){
                failedSampIsoID.vec<-gsub(" ","_",failedSampIsoID.vec)
              }
              #*
              failedSampIsoSD.vec<-c(failedSampIsoSD.vec,sampIsoCheck[[2]]$`SD_d18O16O`)
              failedSampIsoSDC.vec<-c(failedSampIsoSDC.vec,sampIsoCheck[[2]]$`SD_d13C12C`)
              failedSampIsoReason.vec<-c(failedSampIsoReason.vec,sampReason)
            }else{ # add to sample data if it passes
              samps.list[[sampIndex]]<-sampProcess
              sampIndex<-sampIndex+1
            }
          }
        }
      }
    }
  }
  # return data
  ret.list<-list()
  # create refT check output
  if(is.null(failedRefAnNum.vec)){
    print("No analyses failed reference times check")
    failedRefAnNum.vec<-c(NA)
    failedRefID.vec<-c(NA)
    failedRefDat.list<-NA
  }
  failedRefDatRetList<-list()
  failedRefAnNum.df<-as.data.frame(matrix(c(failedRefAnNum.vec,failedRefID.vec),ncol=2))
  colnames(failedRefAnNum.df)<-c("failed_refT_Analysis","failed_refT_Identifier1")
  failedRefDatRetList[[1]]<-failedRefAnNum.df
  failedRefDatRetList[[2]]<-failedRefDat.list
  # create refIso check output
  if(is.null(failedRefIso.vec)){
    print("No analyses failed reference isotope ratio check")
    failedRefIso.vec<-c(NA)
    failedRefIsoID.vec<-c(NA)
    failedRefIsoSD.vec<-c(NA)
    failedRefIsoSDC.vec<-c(NA)
    failedRefIsoReason.vec<-c(NA)
  }
  failedRefIso.df<-as.data.frame(matrix(c(failedRefIso.vec,failedRefIsoID.vec,failedRefIsoSD.vec,failedRefIsoSDC.vec,failedRefIsoReason.vec),ncol=5))
  colnames(failedRefIso.df)<-c("failed_refIso_Analysis","failed_refIso_Identifier1","failed_refIsoSD_d18O16O","failed_refIsoSD_d13C12C","failed_refIsoSD_reason")
  # create refIntensity output
  if(is.null(failedRefInt.vec)){
    print("No analyses failed reference intensity check")
    failedRefInt.vec<-c(NA)
    failedRefIntID.vec<-c(NA)
    failedRefIntRD.vec<-c(NA)
  }
  failedRefInt.df<-as.data.frame(matrix(c(failedRefInt.vec,failedRefIntID.vec,failedRefIntRD.vec),ncol=3))
  colnames(failedRefInt.df)<-c("failed_refInt_Analysis","failed_refInt_Identifier1","failed_refInt_relDiff")
  # create sample peak processing output
  if(is.null(failedSampAnNum.vec)){
    print("No analyses failed sample processing")
    failedSampAnNum.vec<-c(NA)
    failedSampID.vec<-c(NA)
    failedSampReas.vec<-c(NA)
  }
  failedSampAnNum.df<-as.data.frame(matrix(c(failedSampAnNum.vec,failedSampID.vec,failedSampReas.vec),ncol=3))
  colnames(failedSampAnNum.df)<-c("failed_sampProc_Analysis","failed_sampProc_Identifier1","failed_SampProc_Reason")
  # create sampIso output
  if(is.null(failedSampIso.vec)){
    print("No analyses failed sample isotope ratio check")
    failedSampIso.vec<-c(NA)
    failedSampIsoID.vec<-c(NA)
    failedSampIsoSD.vec<-c(NA)
    failedSampIsoSDC.vec<-c(NA)
    failedSampIsoReason.vec<-c(NA)
  }
  failedSampIso.df<-as.data.frame(matrix(c(failedSampIso.vec,failedSampIsoID.vec,failedSampIsoSD.vec,failedSampIsoSDC.vec,failedSampIsoReason.vec),ncol=5))
  colnames(failedSampIso.df)<-c("failed_sampIso_Analysis","failed_sampIso_Identifier1","failed_sampIsoSD_d18O16O","failed_sampIsoSD_d13C12C","failed_sampIsoSD_reason")

  # build return data
  ret.list[[1]]<-failedPkNum.df
  ret.list[[2]]<-failedRefDatRetList#failedRefAnNum.df
  ret.list[[3]]<-failedRefIso.df
  ret.list[[4]]<-failedRefInt.df
  ret.list[[5]]<-failedSampAnNum.df
  ret.list[[6]]<-failedSampIso.df
  ret.list[[7]]<-refs.list
  ret.list[[8]]<-samps.list
  return(ret.list)
}


removeRefAnalysisDXF<-function(filtered.list,refInd=7,sampInd=8,sampProcFailedInd=5,sampIsoFailInd=6){
  retSR.list<-list()
  refs<-filtered.list[[refInd]]
  samps<-filtered.list[[sampInd]]

  failedSampIsoAnNum.vec<-filtered.list[[sampIsoFailInd]][,1]
  failedSampProcAnNum.vec<-filtered.list[[sampProcFailedInd]][,1]
  failedSampAn.vec<-c(failedSampIsoAnNum.vec,failedSampProcAnNum.vec)

  removeInd<-c()
  lenRefsList<-length(refs)
  for(i in seq(1,lenRefsList)){
    anNum<-refs[[i]]$Analysis[1]
    if(anNum %in% failedSampAn.vec){
      removeInd<-c(removeInd,i)
    }
  }
  ## remove all at once
  refs[removeInd]<-NULL #closes the hole

  retSR.list[[1]]<-refs
  retSR.list[[2]]<-samps

  return(retSR.list)
}


replace_space_filenames<-function(files,path){
  #files<-list.files(path)
  #path<-getwd()
  for(i in seq(1:length(files))){
    currFile<-files[i] # get file name
    # check for two underscores in original name!
    if(grepl("  ",currFile,fixed=T)){
      newName<-gsub("  ","_",currFile)
      fullNewPath<-paste(path,"/",newName,sep="")
      fullOldPath<-paste(path,"/",currFile,sep="")
    } else if(grepl(" ", currFile)){
      newName<-gsub(" ","_",currFile) # substitute spaces in filenames with underscores
      fullNewPath<-paste(path,"/",newName,sep="")
      fullOldPath<-paste(path,"/",currFile,sep="")
      file_move(fullOldPath,fullNewPath)
    }

  }
}


rmRefDatDXF<-function(procList){
  refsFiltered<-procList[[7]]
  sampsFiltered<-procList[[8]]
  analysesEq<-analysisNumsEqualDXF(refsFiltered,sampsFiltered)
  #if(!analysesEq[[1]]){#remove ref data
  rm<-removeRefAnalysisDXF(procList)
  #}
  # check if analysis numbers are all the same
  finalProcRefList<-rm[[1]]
  finalProcSampList<-rm[[2]]
  filter<-threadRefSampListsToDF(finalProcRefList,finalProcSampList)
  return(filter)
}


sample_peaks_pre_processDXF<-function(refTimesOutput,vend.df){
  allPeaks<-seq(1,length(vend.df$PeakNr))
  # start with all as samples then remove reference, flush and 1st sample peak
  samplePeaks<-allPeaks
  # get sample peaks
  refVend<-refTimesOutput[[3]]
  refPeaks<-as.numeric(refTimesOutput[[2]]$PkNr)
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


sample_peaks_processDXF<-function(refTimesOutput,vend.df,flushExpT=135,flushTint=15,
                                  firstSampExpT=275,firstSampTint=15){
  allPeaks<-seq(1,length(vend.df$PeakNr))
  # start with all as samples then remove reference, flush and 1st sample peak
  samplePeaks<-allPeaks
  # after 4th reference peak
  pkNrFirstSampleAfterRefs<-as.numeric(refTimesOutput[[2]]$PkNr[4])+1
  # grab sample peak vendor info using reference output
  refPeaks<-as.numeric(refTimesOutput[[2]]$PkNr)
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
      print(paste("No sample peaks detected for Analysis ",analysisNum,sep=""))
    }else{
      print("No sample peaks detected")
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
  flushInd=0

  # find flush peak by time**
  for(i in seq(1,length(procSampleStart))){
    if(abs(flushExpT-procSampleStart[i])<flushTint){ #within 10 secs? or use a different way to ID flush peak?
      flushPeak<-procSample[i,]
      flushInd<-i
    }
  }
  flushNotPresent<-flushInd==0
  if(flushNotPresent){
    #print("No flush peak detected within expected time interval")
  } else{ # remove flush peak
    #print("Flush peak detected within expected time interval")
    procSample<-procSample[-flushInd,]
  }
  # first sample
  firstSampleInd<-which(procSample$PeakNr==pkNrFirstSampleAfterRefs)
  firstSampleVend<-procSample[firstSampleInd,]
  firstSampleTime<-as.numeric(firstSampleVend$Rt)
  # check if other samples present (other than flush)
  if(length(firstSampleInd)==0){
    if(length(refTimesOutput[[3]]$Analysis)>0){ # if there's an analysis number
      analysisNum<-refTimesOutput[[3]]$Analysis[1]
      print(paste("No sample peaks after flush detected in Analysis ",analysisNum,sep=""))
    }else{
      print("No sample peaks after flush detected")
    }
    return()
  }
  # check that the time for the first sample peak is in the expected interval
  sample1withinTime<-abs(firstSampExpT-firstSampleTime)<firstSampTint
  if(sample1withinTime){
    #print("First sample peak detected within expected time interval")
    procSample<-procSample[-firstSampleInd,]
  } else{ # dont remove first sample peak
    print("First sample peak not within expected time interval")
  }
  # check if there are any sample peaks left after processing out the first 2
  pkNrProc.len<-length(procSample$PeakNr)
  if(pkNrProc.len==0){
    if(length(refTimesOutput[[3]]$Analysis)>0){ # if there's an analysis number
      analysisNum<-refTimesOutput[[3]]$Analysis[1]
      print(paste("No sample peaks left after processing in Analysis ",analysisNum,sep=""))
    }else{
      print("No sample peaks left after processing ")
    }
    return()
  }
  #print(paste("Removed peaks at Rts ",round(as.numeric(flushPeak$Rt),3), "s (flush peak) and ",round(firstSampleTime,3),"s (1st sample peak)",sep=""))
  return(procSample)
}


separate_by_analysis_numDXF<-function(vend.df){
  # make lists of data separated by analysis numbers
  analysis_nums<-vend.df$Analysis
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
      #print("End of dataset")
      analysisSetInd.vec<-c(analysisSetInd.vec,i-1)
      analysisData.list[[uniqueAnalysisInd]]<-vend.df[analysisSetInd.vec,]
      break
    }
    # check if the analysis numbers are the same
    if(currAnalysis==prevAnalysis){
      analysisSetInd.vec<-c(analysisSetInd.vec,i-1)
    } else{ # different analysis numbers
      #print("Reached end of analysis set")
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


stand_lm<-function(acceptedMeas.df,isoRname="d18O",dataName){
  ret.list<-list()
  # linear regression
  standard.fit<-lm(acceptAvg_d18O16O~measAvg_d18O16O,data=acceptedMeas.df)
  fit.summ<-summary(standard.fit)
  ptitle<-lm_graph_eqn(standard.fit)
  # plot
  #?plot
  plot(acceptedMeas.df$measAvg_d18O16O,acceptedMeas.df$acceptAvg_d18O16O,
       xlab=paste("measured ",isoRname,sep=""),ylab=paste("accepted ",isoRname,sep=""),pch=20,
       col=topo.colors(dim(acceptedMeas.df)[1]),cex=2)
  title(main=paste(dataName,ptitle,sep=""))
  abline(standard.fit,col="light blue")
  legend("bottomright", inset=.02, title="internal standard",
         rownames(acceptedMeas.df), fill=topo.colors(dim(acceptedMeas.df)[1]),
         horiz=TRUE, cex=0.8)
  ret.list[[1]]<-standard.fit
  ret.list[[2]]<-fit.summ
  return(ret.list)
}


threadRefSampListsToDF<-function(ref.list,sample.list){
  # thread lists into dataframes
  # need total number of elements in the lists
  # sample list
  numSamp<-length(sample.list)
  sampSum<-0
  for(i in seq(1,numSamp)){
    sampIndexLength<-dim(sample.list[[i]])[1]
    sampSum<-sampSum+sampIndexLength
  }
  # reference peaks
  numRef<-length(ref.list)
  refSum<-0
  for(i in seq(1,numRef)){
    refIndexLength<-dim(ref.list[[i]])[1]
    refSum<-refSum+refIndexLength
  }
  # df dims: num cols = num colnames; num row=sums
  # initalize references df
  numRefCols<-length(colnames(ref.list[[1]]))
  ref.df<-as.data.frame(matrix(rep(NA,numRefCols*refSum),ncol=numRefCols))
  colnames(ref.df)<-colnames(ref.list[[1]])
  # initalize samples df
  numSampCols<-length(colnames(sample.list[[1]]))
  samp.df<-as.data.frame(matrix(rep(NA,numSampCols*sampSum),ncol=numSampCols))
  colnames(samp.df)<-colnames(sample.list[[1]])
  # loop through lists and build dfs
  # reference df
  numRefs<-length(ref.list)
  refIndex<-1
  for(i in seq(1,numRefs)){
    # grab ref vend
    refDat<-ref.list[[i]]
    numRefPts<-length(refDat$PeakNr)
    # set indices
    refInd.vec<-seq(refIndex,(refIndex+numRefPts-1))
    # add to df
    ref.df[refInd.vec,]<-refDat
    # update ref index
    refIndex<-refInd.vec[length(refInd.vec)]+1
  }
  # sample df
  numSamps<-length(sample.list)
  sampIndex<-1
  for(i in seq(1,numSamps)){
    # grab ref vend
    sampDat<-sample.list[[i]]
    numSampPts<-length(sampDat$PeakNr)
    # set indices
    sampInd.vec<-seq(sampIndex,(sampIndex+numSampPts-1))
    # add to df
    samp.df[sampInd.vec,]<-sampDat
    # update ref index
    sampIndex<-sampInd.vec[length(sampInd.vec)]+1
  }
  refSamp.list<-list()
  refSamp.list[[1]]<-ref.df
  refSamp.list[[2]]<-samp.df
  return(refSamp.list)
}


vendor_info<-function(file){
  msdat<-iso_read_continuous_flow(file)
  file.info<-msdat %>% iso_get_file_info()
  ident1<-file.info$`Identifier 1`
  vendor_info<-msdat %>% iso_get_vendor_data_table()
  vendor_info.df<-as.data.frame(vendor_info)
  return(vendor_info.df)
}


vendor_info_all<-function(files){
  vi.list<-list()
  for(i in seq(1,length(files))){
    vi.df<-vendor_info(files[i])
    vi.list[[i]]<-vi.df
  }
  return(vi.list)
}


writeFailStats<-function(processed.list,outputID){
  # analyses that failed the peak number check
  pkNrFileName<-paste(outputID,"_failedPkNr.txt",sep="")
  write.table(processed.list[[1]],pkNrFileName,row.names=F,quote=F)
  # analyses that failed the reference peak times check
  refTfileName<-paste(outputID,"_failedRefT.txt",sep="")
  write.table(processed.list[[2]][[1]],refTfileName,row.names=F,quote=F)
  # analyses that failed the reference peak isotope ratio test
  refIsoFileName<-paste(outputID,"_failedRefIso.txt",sep="")
  write.table(processed.list[[3]],refIsoFileName,row.names=F,quote=F)
  # analyses that failed the reference peak intensity similarity test
  refIntFileName<-paste(outputID,"_failedRefInt.txt",sep="")
  write.table(processed.list[[4]],refIntFileName,row.names=F,quote=F)
  # analyses that failed sample peak processing (removal of flush peak and first sample)
  sampProcFileName<-paste(outputID,"_failedSampProc.txt",sep="")
  write.table(processed.list[[5]],sampProcFileName,row.names=F,quote=F)
  # analyses that failed sample peak isotope ratio sd test
  sampIsoFileName<-paste(outputID,"_failedSampIso.csv",sep="")
  write.table(processed.list[[6]],sampIsoFileName,row.names=F,quote=F)
  if(length(processed.list)==9){
    refSampIntFileName<-paste(outputID,"_failedRefSampInt.csv",sep="")
    write.table(processed.list[[9]],refSampIntFileName,row.names=F,quote=F)
  }
}


stand_lm<-function(acceptedMeas.df){
  ret.list<-list()
  # linear regression
  standard.fit<-lm(acceptAvg_d18O16O~measAvg_d18O16O,data=acceptedMeas.df)
  standard.fit
  fit.summ<-summary(standard.fit)
  # plot
  plot(acceptedMeas.df$measAvg_d18O16O,acceptedMeas.df$acceptAvg_d18O16O,
       main=paste("lm for d18O/16O standards:","\nr^2 = ",round(fit.summ$r.squared,6),
                  "\ny = ",round(standard.fit$coefficients[2],2),"x + ",
                  round(standard.fit$coefficients[1],2),
                  sep=""),
       xlab="measured",ylab="accepted",pch=20,col=c("blue","orange","magenta"))
  abline(standard.fit,col="light blue")
  # add legend
  legend("topleft", legend=rownames(acceptedMeas.df),
         col=c("blue","orange","magenta"),pch=20)
  ret.list[[1]]<-standard.fit
  ret.list[[2]]<-fit.summ
  return(ret.list)
}


generic_plot_all_raw<-function(raw.list,fileName,analysis.vec,ID1.vec,date.vec){
  raw.length<-length(raw.list)
  pdf(file=fileName,width=6,height=4)
  par(mfrow=c(2,3))

  for(i in seq(1,raw.length)){
    # get title
    raw.dat<-raw.list[[i]]
    raw.title<-raw.dat$file_id[1]
    raw.title<-paste(ID1.vec[i]," ","\nAnalysis: ",analysis.vec[i],
                     "\nDate: ",date.vec[i],sep="")
    #raw.title
    generic_raw_plot(raw.dat,raw.title)
  }
  par(mfrow=c(1,1))
  dev.off()
}
