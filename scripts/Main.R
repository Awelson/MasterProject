setwd("C:/Users/aaron/Desktop/Github Projects/MasterProject")

library(tidyverse)
library(purrr)
library(magrittr)
library(reshape2)
library(plyr)
library(dplyr)

source("Scripts/ReadFolder.R")

# CRC2, CRC3 Phenotype

x <- read.csv("RawData/NEWCRC/Tumours.20190923.csv")
x %<>% select(Sample, Type)
x %<>% rename(c("Sample" = "SampleName", "Type" = "L1"))

replacement <- c("MSS" = "MMRp", "MSI" = "MMRd")

x$L1 <- replacement[x$L1]

x %<>% distinct # there are duplicates in the phenotype csv for some reason

x <- x[!is.na(x$L1), ]

# CRC2

CRC2 <- ReadFolder2("RawData/NEWCRC/20170501")
CRC2$Dataset <- "CRC2"
CRC2$Type <- "CRC"

  ## Missing Sample Check

  test1 <- CRC2$SampleName %>% unique()
  test2 <- x$SampleName %>% unique()
  
  intersect(test1, test2) %>% length() # 96
  test1 %>% length() # 98, 2 samples w/o phenotype

  setdiff(test1, test2) # samples w/o phenotype

  ## Remove Missing Samples
  
  CRC2 <- CRC2[!(CRC2$SampleName %in% setdiff(test1, test2)), ] # remove samples w/o phenotype

  ## Attach Phenotype
  
  CRC2 <- left_join(CRC2, x, by = "SampleName")

# CRC3

CRC3 <- ReadFolder2("RawData/NEWCRC/20170809")
CRC3$Dataset <- "CRC3"
CRC3$Type <- "CRC"

CRC3 <- CRC3[CRC3$SampleName != "D227036", ] # D227036 is present in both CRC2 and CRC3

  ## Missing Sample Check
  
  test1 <- CRC2$SampleName %>% unique()
  test2 <- x$SampleName %>% unique()
  
  intersect(test1, test2) %>% length() # 98
  test1 %>% length() # 98, MATCH!

  ## Attach Phenotype
  
  CRC3 <- left_join(CRC3, x, by = "SampleName")

# CMMRD Phenotype

x <- readRDS("RawData/CMMRD/SampleData.20220204.v2.rds")
x %<>% select(Sample, Gene)
x %<>% rename(c("Sample" = "SampleName", "Gene" = "L1"))
x %<>% mutate(L1 = if_else(L1 == "Control", "MMRp", "MMRd"))

# CMMRD

CMMRD <- ReadFolder2("RawData/CMMRD/byUmiGenotypesRearranged")
CMMRD$Dataset <- "CMMRD"
CMMRD$Type <- "CMMRD"

  ## Missing Sample Check
  
  test1 <- CMMRD$SampleName %>% unique()
  test2 <- x$SampleName %>% unique()
  
  intersect(test1, test2) %>% length() # 139
  test1 %>% length() # 139, MATCH!
  
  ## Attach Phenotype

  CMMRD <- left_join(CMMRD, x, by = "SampleName")

# EC1 (MAN) Phenotype

x <- read.csv("RawData/EC/ManUMI/Manchester.EC.Samples.20240111.csv")
x %<>% select(SequencingID, Type)
x %<>% rename(c("SequencingID" = "SampleName", "Type" = "L1"))

# EC1 (MAN)

EC1 <- ReadFolder2("RawData/EC/MANUMI/20220622")
EC1$Dataset <- "EC1"
EC1$Type <- "EC"

  ## Missing Sample Check

  test1 <- EC1$SampleName %>% unique()
  test2 <- x$SampleName %>% unique()
  
  intersect(test1, test2) %>% length() # 151
  test1 %>% length() # 167, 16 missing

  setdiff(test1, test2) # samples w/o phenotype

  EC1 <- EC1[!(EC1$SampleName %in% setdiff(test1, test2)), ] # remove samples w/o phenotype

  ## Attach Phenotype
  
  EC1 <- left_join(EC1, x, by = "SampleName")

# EC (Ohio) Phenotype

x <- read.csv("RawData/EC/OhioUMI/Ohio.SampleData.20240111.csv")
x %<>% select(SequencingID, Type)
x %<>% rename(c("SequencingID" = "SampleName", "Type" = "L1"))
x$L1 <- replacement[x$L1]

# EC (Ohio)

ECOhio1 <- ReadFolder2("RawData/EC/OhioUMI/20220318")

ECOhio2 <- ReadFolder2("RawData/EC/OhioUMI/20220401")

EC2 <- bind_rows(ECOhio1, ECOhio2)

EC2$Dataset <- "EC2"
EC2$Type <- "EC"

  ## Missing Sample Check

  test1 <- EC2$SampleName %>% unique()
  test2 <- x$SampleName %>% unique()

  intersect(test1, test2) %>% length() # 184
  test1 %>% length() # 184, MATCH!

  EC2 <- left_join(EC2, x, by = "SampleName")

# CRC1 (Jenny) Phenotype

x <- read.csv("RawData/OGCRC/LS.CRCs.BRAF.csv")
x %<>% select(Sample, Type)
x %<>% rename(c("Sample" = "SampleName", "Type" = "L1"))

x$L1 <- gsub("LS MSI-H", "MMRd", x$L1)
x$L1 <- gsub("MSI-H", "MMRd", x$L1)
x$L1 <- gsub("MSS", "MMRp", x$L1)

# CRC1 (JENNY)

CRC1 <- ReadFolder2("RawData/OGCRC/UMImsiMNRsJenny210421")
CRC1$Dataset <- "CRC1"
CRC1$Type <- "CRC"

  ## Renaming Samples
  
  CRC1$SampleName <- gsub("258136A4", "258136", CRC1$SampleName)
  CRC1$SampleName <- gsub("269965B7", "269965", CRC1$SampleName)

  ## Missing Sample Check
  
  test1 <- CRC1$SampleName %>% unique()
  test2 <- x$SampleName %>% unique()
  
  intersect(test1, test2) %>% length() # 128
  test1 %>% length() # 128, MATCH!
  
  ## Attach Phenotype

  CRC1 <- left_join(CRC1, x, by = "SampleName")

# All

All <- bind_rows(CRC1, CRC2, CRC3, CMMRD, EC1, EC2)
All <- All[, c("SampleName", "Marker", "Allele (length)", "Allele (SNP)", "Counts", "UMI", "L1", "Dataset", "Type")]


  ## Check for samples with no classification
  
  All[is.na(All$L1), ] # it is empty, GOOD

  ## Check for samples with > 12 missing tables
  
  MissingTables2("RawData/NewCRC/20170501") # Empty
  MissingTables2("RawData/NewCRC/20170809") # Empty
  MissingTables2("RawData/CMMRD/byUmiGenotypesRearranged") # Empty
  MissingTables2("RawData/EC/ManUMI/20220622") # Empty
  MissingTables2("RawData/EC/OhioUMI/20220318") # OhioECT181
  MissingTables2("RawData/EC/OhioUMI/20220401") # OhioECV81
  MissingTables2("RawData/OGCRC/UMImsiMNRsJenny210421") # Empty
  
  ## Removal of identified samples

  All <- All[!(All$SampleName %in% c("OhioECV81", "OhioECT181")), ]

# Write CSV

write.csv(All, "Data/All.csv", row.names = FALSE)
