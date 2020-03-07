# 2Cluster tumor and remove spurious tumors.R
# @author: Chuan Li (chuanli@stanford.edu)

### Cluster reads and remove spurious tumors
# Any "tumours" that are within a Hamming distance of two from a larger tumour is assigned as "spurious tumours", 
# which are likely to be resulting from sequencing or PCR error, and are removed from subsequent analysis. 
# we only focus on the ones with the correct length 

FileName <- paste0("InputFile", "_sgIDCounts.txt") # replace InputFile with the input
FileOut <- sub(pattern = "(*)sgIDCounts.txt", "\\1BarcodeClean.txt", FileName)

SamplesFull <- read.csv(file = FileName, header = F, sep=",", stringsAsFactors = F)
names(SamplesFull) <- c("sgID","BC","Count")
SamplesFull <- SamplesFull[order(-SamplesFull$Count),]
SamplesFull <- SamplesFull[SamplesFull$Count>1,]

SamplesClean <- c()
sgIDAll <- unique(SamplesFull$sgID)

# find whether the sequence contain any N
# remove these ones
findN <- function(String){
  return(grepl('N',String))
}

BoolfindN <- sapply(SamplesFull$BC, FUN = findN)
SamplesFull <- SamplesFull[!BoolfindN, ]

## keep only the one that are more than 2 steps away from the focal sequence
myDist <- function(Str1){
  seq1<-unlist(strsplit(Str1,split=""))
  seq2<-unlist(strsplit(myBC,split=""))
  Dist=0
  for (i in 1:length(seq1)){
    if (seq1[i] != seq2[i]){
      Dist = Dist + 1
    }
    if (Dist > 2){
      return(TRUE) # keep the sequence
    }
  }
  return(FALSE) # filter the sequence
}

for (sgID in sgIDAll){ 
  Samples <- SamplesFull[SamplesFull$sgID == sgID,]
  N <- dim(Samples)[1]
  BoolAll <- rep(TRUE, N)

  for (iBC in 1:N){
    if (Samples$Count[iBC] < 5){
      break
    }
    
    # For current barcode
    myBC <- Samples$BC[iBC]
    
    # Filter through the rest of barcode
    if (iBC==N){
      break # if it reaches the end, break
    }
    
    if (iBC == N-1){
      DistAll <- myDist(Samples$BC[(iBC+1):N])
    }else{
      DistAll <- sapply(Samples$BC[(iBC+1):N], FUN = myDist)
    }
    
    DistAll <- c(rep(TRUE, iBC), DistAll)
    BoolAll[!DistAll] <- FALSE
  }
  
  Samples <- Samples[BoolAll, ]
  SamplesClean <- rbind(SamplesClean, Samples)
}

write.table(SamplesClean, file = FileOut, append = FALSE, quote = FALSE, sep = "\t", row.names = F,
            col.names = TRUE)
