# Utility functions to quantify genotype-specific treatment responses
# @author: Chuan Li (chuanli@stanford.edu)

#############################
# Section 1. bootstrap mouse and tumor to control for variations across mice and across tumors
#############################
### Usng bootstrap to calculate sd, CI and calculate p-values
# indexBoot: a subfunction that samples sgID by the corresponding fractions
indexBoot <- function(Group){ # Group <- paste0(DataAll$sgID, DataAll$Lib)
  GroupNames <- unique(Group)
  IndexAll <- c()
  for (igroup in GroupNames){
    IndexOri <- which(Group==igroup)
    IndexBoot <- sample(IndexOri, length(IndexOri), replace = T)
    IndexAll <- c(IndexAll, IndexBoot)
  }
  return(sort(IndexAll))
}

# this function allows bootstraping by mouse
indexBootMouse <- function(DataAll, Group){ # Group <- DataAll$Lib
  GroupNames <- unique(Group)
  Mouse <- table(sample(GroupNames, size = length(GroupNames), replace = T))
  DataPseudo <- c()
  for (iMouse in names(Mouse)){ # iMouse <- names(Mouse)[1]
    for (iSamp in 1:Mouse[iMouse]){
      Data1 <- DataAll[DataAll$Lib == iMouse,]
      Data1$Lib <- paste0(Data1$Lib, "-", iSamp)
      DataPseudo <- rbind(DataPseudo, Data1)
    }
  }
  return(DataPseudo)
}

#######################
# Section 2. ScoreRTN to quantify genotype-specific treatment responses
#######################

# Calculate the tumor number relative to Inert
getRelTumorNum <- function(DataAll){
  # calculate the fraction of sgID changes relative to the total
  TumorNum <- c()
  for (sgID in sgIDAll){
    TumorNum <- c(TumorNum, sum(DataAll$sgID == sgID))
  }
  Ind <- which(sgIDAll %in% Inert)
  return(TumorNum/sum(TumorNum[Ind]))
}

# Calculate the tumor burden relative to Inert
getRelTumorBurden <- function(DataAll){
  # calculate the fraction of sgID changes relative to the total
  CellNum <- c()
  for (sgID in sgIDAll){
    CellNum <- c(CellNum, sum(DataAll$CellNum[DataAll$sgID == sgID]))
  }
  Ind <- which(sgIDAll %in% Inert)
  return(CellNum/sum(CellNum[Ind]))
}

#### Calculate the relative tumor number and tumor burden
getRelTNTB <- function(DataTr, DataUn, CellCutoff, Adjust="ON"){
  Shift = binarySearchShift(DataTr, DataUn)
  CellCutoff2 = CellCutoff*Shift
  
  # calculate the relative tumor number and tumor burden above the cutoff.
  DataTr1 <- DataTr[DataTr$CellNum >= CellCutoff2,]
  DataUn1 <- DataUn[DataUn$CellNum >= CellCutoff,]
  TN <- getRelTumorNum(DataTr1)/getRelTumorNum(DataUn1)
  TB <- getRelTumorBurden(DataTr1)/getRelTumorBurden(DataUn1)
  
  return(list(TN,TB))
}

# Estimate the drug effect using the binary search algorithm
binarySearchShift <- function(DataTr, DataUn){
  Upper = 1 # when estimating the drug effect, Upper = 4
  Lower = 0.01
  Shift = 1
  DataTrInert <- DataTr[DataTr$sgID %in% Inert, ]
  DataUnInert <- DataUn[DataUn$sgID %in% Inert, ]
  NumTr <- median(table(DataTrInert$Lib[DataTrInert$CellNum >= CellCutoff]))
  NumUn <- median(table(DataUnInert$Lib[DataUnInert$CellNum >= CellCutoff]))

  if (NumTr >= NumUn){
    return (1)
  }
  iter = 0
  while (Upper-Lower >= 0.001){
    iter = iter + 1
    Shift = (Upper + Lower)/2
    
    # shift untreated to match
    DataUnSft <- DataUnInert
    DataUnSft$CellNum <- DataUnInert$CellNum*Shift
    NumUn <- median(table(DataUnSft$Lib[DataUnSft$CellNum >= CellCutoff]))
    if (NumTr >= NumUn){ # shifted too much
      Lower = Shift
    }else{ # Not enough shift
      Upper = Shift
    }

    if (iter > 50){
      break
    }
  }
  return (Shift)
}

getBootStrapTNTB <- function(DataTr,DataUn, CellCutoff, Mode="Mouse", nSimu=100){
  ResultTN <- c()
  ResultTB <- c()
  for(iSimu in 1:nSimu){
    if (Mode == "Mouse"){ # we no longer use Tumor mode, which only bootstrap tumors
      # we always first bootstrap mouse and then bootstrap tumor
      DataTrPseudo <- indexBootMouse(DataTr, Group = DataTr$Lib)
      DataUnPseudo <- indexBootMouse(DataUn, Group = DataUn$Lib)
      DataTrPseudo <- DataTrPseudo[indexBoot(paste0(DataTrPseudo$sgID, DataTrPseudo$Lib)), ]
      DataUnPseudo <- DataUnPseudo[indexBoot(paste0(DataUnPseudo$sgID, DataUnPseudo$Lib)), ]
    }
    OutPseudo <- getRelTNTB(DataTrPseudo, DataUnPseudo, CellCutoff = CellCutoff)
    ResultTN <- cbind(ResultTN, OutPseudo[[1]])
    ResultTB <- cbind(ResultTB, OutPseudo[[2]])
  }
  return(list(ResultTN, ResultTB))
}

#######################
# Section 3. ScoreRGM or ScoreRLM to quantify genotype-specific treatment responses
#######################
## For choice of relative geometric mean or relative lognormal mean
## it depends on the tumor size distribution
## can run power analysis to figure out which metric works better

getAdaptive1 <- function(DataTr, DataUn, CellCutoff = 500, Adjust="ON"){

  Shift = binarySearchShift(DataTr, DataUn) # shift to 0.303
  CellCutoff2 = CellCutoff*Shift
  
  # So find the new cutoff for DataSet2
  # Figure out the shift in DataUn to achieve the same tumor number for Inert
  DataTrInert <- DataTr[DataTr$sgID %in% Inert, ]
  DataUnInert <- DataUn[DataUn$sgID %in% Inert, ]
  
  # calculate the ratio relative to Inert in each Untreated mouse
  ResultUn <- aggregate(DataUn["sgID"], DataUn["Lib"], table, simplify = T)
  ResultUn <- ResultUn[2][[1]]
  ResultUnFin <- ResultUn/rowSums(ResultUn[,colnames(ResultUn) %in% Inert])
  Ratio <- apply(ResultUnFin, 2, median)
  Ratio <- Ratio[names(Ratio) %in% sgIDAll]
  
  # take out the number of Inert Palbociclib
  NumTrInert <- table(DataTr$Lib[DataTr$CellNum >= CellCutoff2 & DataTr$sgID %in% Inert])
  NumUnInert <- table(DataUn$Lib[DataUn$CellNum >= CellCutoff & DataUn$sgID %in% Inert])
  
  DataTr2 <- c()
  for (sgID in sgIDAll){
    NumTrsgID <- round(NumTrInert*Ratio[sgID])
    
    # calculate the current number for treated
    for (iLib in names(NumTrsgID)){
      DataTrsgID <- DataTr[DataTr$sgID == sgID & DataTr$Lib == iLib,]
      ExpNum = NumTrsgID[iLib]
      ObsNum = dim(DataTrsgID)[1]
      if (ExpNum <= ObsNum){
        DataTr2 <- rbind(DataTr2, DataTrsgID[1:ExpNum,])
      }else{
        # when we don't have enough tumors, fill in tumors at the detection limit
        # but given the current cutoff, this never happens.
        Fake <- DataTrsgID[rep(1,ExpNum - ObsNum),]
        Fake$Count <- 0
        Fake$CellNum <- DataTrsgID$CellNum[ObsNum]
        DataTr2 <- rbind(DataTr2, DataTrsgID)
        DataTr2 <- rbind(DataTr2, Fake)
      }
    }
  }
  
  DataUn2 <- c()
  for (sgID in sgIDAll){
    NumUnsgID <- round(NumUnInert*Ratio[sgID])
    
    # calculate the current number for treated
    for (iLib in names(NumUnsgID)){
      DataUnsgID <- DataUn[DataUn$sgID == sgID & DataUn$Lib == iLib,]
      ExpNum = NumUnsgID[iLib]
      ObsNum = dim(DataUnsgID)[1]
      if (ExpNum <= ObsNum){
        DataUn2 <- rbind(DataUn2, DataUnsgID[1:ExpNum,])
      }else{
        Fake <- DataUnsgID[rep(1,ExpNum - ObsNum),]
        Fake$Count <- 0
        Fake$CellNum <- DataUnsgID$CellNum[ObsNum]
        DataUn2 <- rbind(DataUn2, DataUnsgID)
        DataUn2 <- rbind(DataUn2, Fake)
      }
    }
  }
  return(list(DataTr2,DataUn2))
}

getBootStrapResult <- function(DataTr, DataUn, Mode="Mouse", nSimu = 100, Shift="Tr"){
  # return LN mean of each bootstrapped result
  ResultLNTr <- c()
  ResultLNUn <- c()
  
  for (iSimu in 1:nSimu){ # iSimu=1
    # generate a pseudo sample
    
    if (Mode == "Mouse"){ # we no longer use Tumor mode
      # first bootstrap mouse and then bootstrap tumor
      DataTrPseudo <- indexBootMouse(DataTr, Group = DataTr$Lib)
      DataUnPseudo <- indexBootMouse(DataUn, Group = DataUn$Lib)
      DataTrPseudo <- DataTrPseudo[indexBoot(paste0(DataTrPseudo$sgID, DataTrPseudo$Lib)), ]
      DataUnPseudo <- DataUnPseudo[indexBoot(paste0(DataUnPseudo$sgID, DataUnPseudo$Lib)), ]
    }
    
    Output <- getAdaptive1(DataTrPseudo, DataUnPseudo, CellCutoff = CellCutoff, Adjust="ON")
    
    DataTr2 <- Output[[1]]
    DataUn2 <- Output[[2]]
    
    # calculate the percentile
    LNPseudo <- getLNMean(DataTr2)
    ResultLNTr <- cbind(ResultLNTr, LNPseudo)
    LNPseudo <- getLNMean(DataUn2)
    ResultLNUn <- cbind(ResultLNUn, LNPseudo)
  }
  return(list(ResultLNTr, ResultLNUn))
}


#######################
# Section 4. Calculate p-values and confidence intervals
#######################
getLower <- function(Data){
  return (quantile(Data, 0.025, na.rm = T))
}

getUpper <- function(Data){
  return (quantile(Data, 0.975, na.rm = T))
}

getEmpiricalP <- function(Data, Baseline=1, NoAdjust = FALSE, Side=2){
  # Side can be adjusted to return one sided or two sided p-values
  # NoAdjust can either be false or the value you want to put in
  # which corresponds to the observed value of that metric before bootstrapping
  Data <- Data[is.finite(Data)]
  Counts = (sum(Data > Baseline) + sum(Data == Baseline)/2)/length(Data)
  
  if (Side == 2){ # Two sided test
    return(min(Counts,1-Counts)*2)
  }else{
    
    if (!NoAdjust){ # One sided test
      return(min(Counts,1-Counts))
    }else{
      if (NoAdjust == Baseline){
        return (0.5)
      }else if(NoAdjust > Baseline){
        return (1-Counts)
      }else{
        return (Counts)
      }
    }
  }
}

getEmpiricalPTwoSets <- function(Array1, Array2){
  # The following two subfunctions calculates the P-value by comparing two bootstrap samples
  # Return how often the first argument is larger than the second argument
  Counts = 0
  for (Baseline in Array2){ # User Array2 untreated as the baseline
    Counts = Counts + sum(Array1 > Baseline) + sum(Array1 == Baseline)/2
  }
  P = Counts/length(Array1)/length(Array2)
  return (P)
}

getEmpiricalCIDivide <- function(Array1, Array2){
  ArrayDivide = c()
  for (Baseline in Array2){
    ArrayDivide = c(ArrayDivide, Array1/Baseline)
  }
  ArrayDivide <- ArrayDivide[is.finite(ArrayDivide)]
  return(c(getLower(ArrayDivide), getUpper(ArrayDivide)))
}

getPByCompare <- function(BootTr, BootUn){
  PAll <- c()
  for (isgID in 1:dim(BootTr)[1]){
    PAll <- c(PAll, getEmpiricalPTwoSets(BootTr[isgID,], BootUn[isgID,]))
  }
  return(PAll)
}

getCIByCompare <- function(BootTr, BootUn){
  CI <- c()
  for (isgID in 1:dim(BootTr)[1]){
    CI <- rbind(CI, getEmpiricalCIDivide(BootTr[isgID,], BootUn[isgID,]))
  }
  return(CI)
}

