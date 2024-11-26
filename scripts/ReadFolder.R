setwd("C:/Users/aaron/Desktop/Github Projects/MasterProject")

library(plyr)
library(tidyverse)
library(purrr)
library(magrittr)

TestFunction1<-function(SingleMarkerTable){
  EmptyTable<-data.frame(Allele=NA,Count=0,Total=0,RelativeFrequency=0)
  if( empty(SingleMarkerTable))  return(EmptyTable)
  CleanMarkerTable<-  SingleMarkerTable[!is.na(rownames(SingleMarkerTable)),,drop=F]
  if( empty(CleanMarkerTable))  return(EmptyTable)
  SumsVector<-rowSums(CleanMarkerTable)
  OutputDataFrame<-data.frame(Allele=names(SumsVector),Count=SumsVector)
  Total<-sum(OutputDataFrame$Count)      
  OutputDataFrame$Total<-Total   
  OutputDataFrame$RelativeFrequency<- OutputDataFrame$Count/Total
  OutputDataFrame         
}

TestFunction2<-function(MarkerList){
  filter.list <- MarkerList %>% keep(grepl("_SNP1", names(MarkerList)))
  filter.list <- MarkerList %>% keep(!(names(MarkerList) %in% c("KRAS_SNP1", "BRAF_SNP1"))) %>% discard(~ is_empty(.))
  filter.list.processed <- lapply(MarkerList, TestFunction1)
  filter.list.processed <- map2(filter.list.processed, names(MarkerList), ~ mutate(.x, Marker = .y))
  bind_rows(filter.list.processed)
}

ReadGenotypeFolder <- function(FolderPath){
  
  FileList <- lapply(list.files(FolderPath, full.names=TRUE), function(file) {
    source(file, local = TRUE)
    return(Result)
  })
  
  FullTableList <- lapply(FileList, TestFunction2)
  SampleNames <- sapply(strsplit(list.files(FolderPath), "_"), function(x){x[1]})
  FullTableList <- map2(FullTableList, SampleNames, ~ mutate(.x, Sample = .y))
  FullTable <- bind_rows(FullTableList)
  rownames(FullTable) <- NULL
  FullTable$Marker <- gsub("_SNP1", "", FullTable$Marker)
  
  FullTable
}

# CMMRD.byUmiGenotypesRearranged

  samples <- readRDS("RawData/CMMRD/SampleData.20220204.v2.rds")
  samples <- samples %>% select("Sample", "Gene")
  
  Attach <- function(data){
    
    data <- left_join(data, samples, by = "Sample") %>% rename("SampleName" = "Sample")
    data <- data %>% mutate("L1" = if_else(Gene == "Control", "Control", "CMMRD")) %>% rename("L2" = "Gene")
    data <- data %>% mutate(L2 = if_else(L2 == "Control", NA, L2))
    data$L3 <- NA
    
    return(data)
  }
  
  ## CMMRD.byUmiGenotypesRearranged.D0
  
  D0 <- ReadGenotypeFolder("RawData/CMMRD/byUmiGenotypesRearranged/D0")
  D0final <- Attach(D0)
  
  ## CMMRD.byUmiGenotypesRearranged.D1
  
  D1 <- ReadGenotypeFolder("RawData/CMMRD/byUmiGenotypesRearranged/D1")
  D1final <- Attach(D1)
  
  ## CMMRD.byUmiGenotypesRearranged.D2
  
  D2 <- ReadGenotypeFolder("RawData/CMMRD/byUmiGenotypesRearranged/D2")
  D2final <- Attach(D2)
  
  ## full CMMRD data
  
  D0final$UMI <- "D0"
  D1final$UMI <- "D1"
  D2final$UMI <- "D2"
  
  CMMRD <- bind_rows(D0final, D1final)
  CMMRD <- bind_rows(CMMRD, D2final)
  
  write.csv(CMMRD, file = "Data/withoutSNP/CMMRD.csv", row.names = FALSE)

# NEWCRC

  ## Loading phenotype data
  
  samples <- read.csv("RawData/NEWCRC/Tumours.20190923.csv")
  samples <- samples %>% select(Sample,Type)
  
  replacement <- c("MSS" = "MMRp", "MSI" = "MMRd")
  
  ## NEWCRC/20170501
  CRC1_D0 <- ReadGenotypeFolder("RawData/NEWCRC/20170501/D0")
  CRC1_D1 <- ReadGenotypeFolder("RawData/NEWCRC/20170501/D1")
  CRC1_D2 <- ReadGenotypeFolder("RawData/NEWCRC/20170501/D2")
  
  CRC1_D0$UMI <- "D0"
  CRC1_D1$UMI <- "D1"
  CRC1_D2$UMI <- "D2"
  
  CRC1 <- bind_rows(CRC1_D0, CRC1_D1)
  CRC1 <- bind_rows(CRC1, CRC1_D2)
  
  CRC1final <- left_join(CRC1, samples, by = "Sample") %>% rename("SampleName" = "Sample") %>% rename("L1" = "Type") 
  CRC1final$L2 <- NA
  CRC1final$L3 <- NA
  CRC1final$L1 <- replacement[CRC1final$L1]
  
  ## Removing samples which for which Type is unknown
  CRC1final <- CRC1final %>% filter(!SampleName %in% c("E44", "L0303"))
  CRC1final %<>% unique()
  CRC1final %<>% filter(!is.na(RelativeFrequency))
  
  ## Write CSV
  write.csv(CRC1final, file = "Data/withoutSNP/CRC1.csv", row.names = FALSE)
  
  
  ## NEWCRC/20170809
  CRC2_D0 <- ReadGenotypeFolder("RawData/NEWCRC/20170809/D0")
  CRC2_D1 <- ReadGenotypeFolder("RawData/NEWCRC/20170809/D1")
  CRC2_D2 <- ReadGenotypeFolder("RawData/NEWCRC/20170809/D2")
  
  CRC2_D0$UMI <- "D0"
  CRC2_D1$UMI <- "D1"
  CRC2_D2$UMI <- "D2"
  
  CRC2 <- bind_rows(CRC2_D0, CRC2_D1)
  CRC2 <- bind_rows(CRC2, CRC2_D2)
  
  CRC2final <- left_join(CRC2, samples, by = "Sample") %>% rename("SampleName" = "Sample") %>% rename("L1" = "Type")
  CRC2final$L2 <- NA
  CRC2final$L3 <- NA
  CRC2final$L1 <- replacement[CRC2final$L1]
  
  ## Write CSV
  write.csv(CRC2final, file = "Data/withoutSNP/CRC2.csv", row.names = FALSE)
  
# OGCRC
  
  CRC0 <- read.csv("Data/withoutSNP/CRC0.csv")
  CRC0 %<>% rename("L1" = "MMR", "L2" = "Type", "L3" = "Gene")
  write.csv(CRC0, file = "Data/withoutSNP/CRC0.csv", row.names = FALSE)  
  