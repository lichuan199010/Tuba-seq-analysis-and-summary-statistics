# 3 Convert read number to cell number based on the 3 spike-in cell lines
# @author: Chuan Li (chuanli@stanford.edu)

### Make a summary the data ####
rm(list=ls())

# read in sgID information and sample information
sgIDList <- read.csv(file = "SgIDList.txt", sep = "\t", stringsAsFactors = F)

SampleFinal <- read.csv(file = "SampleInfo_All.txt", sep = "\t", stringsAsFactors = F)

# read in number of spike-ins for each sample
SpiInfo <- read.csv(file = "BarcodeCountSummary.txt", sep = "\t", stringsAsFactors = F)
SpiBCCore <- SpiInfo$Barcode

FakeBC <- c()
Read2Eql <- c() # sequencing depth
ReadNum <- c()
CellNum <- c()
sgIDNum <- c()
FileID <- c()

# Parse through each sgID-Count files
# convert read number to cell number and add mouse information
for (i in 1:dim(SampleFinal)[1]){
  FileNames <- paste0(SampleFinal$UniqueCode[i], "_BarcodeClean.txt")
  Data <- read.csv(file = FileNames, header = T, sep="\t", stringsAsFactors = F)
  Spi <- Data$Count[Data$sgID == "Spi" & Data$BC %in% SpiBCCore]
  
  FakeBC <- c(FakeBC, sum(Data$sgID == "Spi" & !Data$BC %in% SpiBCSet))
  Correct <- sum(Spi[1:3])/3/500000
  Read2Eql <- c(Read2Eql, 2/Correct)
  
  MouseInfo <- SampleFinal[i, c(5,13:18)]
  Data$CellNum <- Data$Count/Correct
  
  Data <- cbind(Data, MouseInfo)
  ReadNum <- c(ReadNum, sum(Data$Count))
  CellNum <- c(CellNum, sum(Data$CellNum))
  write.table(Data, paste0(SampleFinal$UniqueCode[i], "_BarcodeFinal.txt"), append = FALSE, quote = FALSE, sep = "\t", row.names = F,
              col.names = TRUE)
  sgIDNum <- c(sgIDNum, dim(Data)[1])
  FileID <- c(FileID, SampleFinal$UniqueCode[i])
}

# output for the conversion and sequencing related metrics etc.
Output <- data.frame(FileID, FakeBC, Read2Eql, ReadNum, CellNum, BCNum=sgIDNum)
Output <- cbind(Output, SampleFinal)
write.table(Output, "SampleInfo.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = F,
            col.names = TRUE)