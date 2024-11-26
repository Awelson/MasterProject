setwd("C:/Users/aaron/Desktop/Github Projects/MasterProject")

library(tidyverse)
library(purrr)
library(magrittr)
library(reshape2)
library(plyr)
library(dplyr)

source("Scripts/ReadFolder2.R")

# CRC1, CRC2 Phenotype

x <- read.csv("RawData/NEWCRC/Tumours.20190923.csv")
x %<>% select(Sample, Type)
x %<>% rename(c("Sample" = "SampleName", "Type" = "L1"))

replacement <- c("MSS" = "MMRp", "MSI" = "MMRd")

x$L1 <- replacement[x$L1]

x %<>% distinct # there are duplicates in the phenotype csv for some reason

x <- x[!is.na(x$L1), ]

# CRC1

CRC1 <- ReadFolder2("RawData/NEWCRC/20170501")
CRC1$L2 <- NA
CRC1$L3 <- NA
CRC1$Dataset <- "CRC1"
CRC1$Type <- "CRC"

## Missing Sample Check

test1 <- CRC1$SampleName %>% unique()
test2 <- x$SampleName %>% unique()

intersect(test1, test2) %>% length() # 96
test1 %>% length() # 98, 2 samples w/o phenotype

setdiff(test1, test2) # samples w/o phenotype

CRC1 <- CRC1[!(CRC1$SampleName %in% setdiff(test1, test2)), ] # remove samples w/o phenotype

CRC1 <- left_join(CRC1, x, by = "SampleName")

# CRC2

CRC2 <- ReadFolder2("RawData/NEWCRC/20170809")

CRC2 <- CRC2[CRC2$SampleName != "D227036", ] # removing D227036 in CRC2 (it is present in both CRC1 and CRC2)

## Missing Sample Check

test1 <- CRC2$SampleName %>% unique()
test2 <- x$SampleName %>% unique()

intersect(test1, test2) %>% length() # 98
test1 %>% length() # 98, MATCH!

CRC2 <- left_join(CRC2, x, by = "SampleName")
CRC2$L2 <- NA
CRC2$L3 <- NA

CRC2$Dataset <- "CRC2"
CRC2$Type <- "CRC"

# CMMRD Phenotype

x <- readRDS("RawData/CMMRD/SampleData.20220204.v2.rds")
x %<>% select(Sample, Genotype, Gene)
x %<>% rename(c("Sample" = "SampleName", "Gene" = "L1", "Genotype" = "L2"))
x %<>% mutate(L1 = if_else(L1 == "Control", "MMRp", "MMRd"), L2 = if_else(L2 == "Control", NA, L2))

# CMMRD

CMMRD <- ReadFolder2("RawData/CMMRD/byUmiGenotypesRearranged")
CMMRD <- left_join(CMMRD, x, by = "SampleName")
CMMRD$L3 <- NA

## Missing Sample Check

test1 <- CMMRD$SampleName %>% unique()
test2 <- x$SampleName %>% unique()

intersect(test1, test2) %>% length() # 139
test1 %>% length() # 139, MATCH!


CMMRD$Dataset <- "CMMRD"
CMMRD$Type <- "CMMRD"

# EC (MAN) Phenotype

x <- read.csv("RawData/EC/ManUMI/Manchester.EC.Samples.20240111.csv")
x %<>% select(SequencingID, Type)
x %<>% rename(c("SequencingID" = "SampleName", "Type" = "L1"))

# EC (MAN)

ECMAN <- ReadFolder2("RawData/EC/MANUMI/20220622")
ECMAN$L2 <- NA
ECMAN$L3 <- NA
ECMAN$Dataset <- "Man"
ECMAN$Type <- "EC"

test1 <- ECMAN$SampleName %>% unique()
test2 <- x$SampleName %>% unique()

intersect(test1, test2) %>% length() # 151
test1 %>% length() # 167, 16 missing

setdiff(test1, test2) # samples w/o phenotype

ECMAN <- ECMAN[!(ECMAN$SampleName %in% setdiff(test1, test2)), ] # remove samples w/o phenotype

ECMAN <- left_join(ECMAN, x, by = "SampleName")

# EC (Ohio) Phenotype

x <- read.csv("RawData/EC/OhioUMI/Ohio.SampleData.20240111.csv")
x %<>% select(SequencingID, Type)
x %<>% rename(c("SequencingID" = "SampleName", "Type" = "L1"))
x$L1 <- replacement[x$L1]

# EC (Ohio)

ECOhio1 <- ReadFolder2("RawData/EC/OhioUMI/20220318")

ECOhio2 <- ReadFolder2("RawData/EC/OhioUMI/20220401")

ECOhio <- bind_rows(ECOhio1, ECOhio2)

ECOhio$L2 <- NA
ECOhio$L3 <- NA
ECOhio$Dataset <- "Ohio"
ECOhio$Type <- "EC"

test1 <- ECOhio$SampleName %>% unique()
test2 <- x$SampleName %>% unique()

intersect(test1, test2) %>% length() # 184
test1 %>% length() # 184, MATCH!

ECOhio <- left_join(ECOhio, x, by = "SampleName")

# Jenny Phenotype

x <- read.csv("RawData/OGCRC/LS.CRCs.BRAF.csv")
x %<>% select(Sample, Type)
x %<>% rename(c("Sample" = "SampleName", "Type" = "L1"))

x$L1 <- gsub("LS MSI-H", "MMRd", x$L1)
x$L1 <- gsub("MSI-H", "MMRd", x$L1)
x$L1 <- gsub("MSS", "MMRp", x$L1)

# JENNY

Jenny <- ReadFolder2("RawData/OGCRC/UMImsiMNRsJenny210421")
Jenny$L2 <- NA
Jenny$L3 <- NA
Jenny$Dataset <- "Jenny"
Jenny$Type <- "CRC"

## Renaming Samples

Jenny$SampleName <- gsub("258136A4", "258136", Jenny$SampleName)
Jenny$SampleName <- gsub("269965B7", "269965", Jenny$SampleName)

## Missing Sample Check

test1 <- Jenny$SampleName %>% unique()
test2 <- x$SampleName %>% unique()

intersect(test1, test2) %>% length() # 128
test1 %>% length() # 128, MATCH!

Jenny <- left_join(Jenny, x, by = "SampleName")

# All

All <- bind_rows(Jenny, CRC1, CRC2, CMMRD, ECMAN, ECOhio)
All <- All[, c("SampleName", "Marker", "Allele (length)", "Allele (SNP)", "Counts", "UMI", "L1", "L2", "L3", "Dataset", "Type")]

NoL1 <- All[is.na(All$L1), ] # Check Samples with no MMRp / MMRd classification, NoL1 is empty, GOOD

## Identifying "Sample, Marker" combinations with zero (total) counts in the "zero length" allele 

test <- All %>% group_by(SampleName, `Allele (length)`, Marker, UMI, L1, L2, L3, Dataset, Type) %>%
  reframe(
    Counts = sum(Counts))

test <- test %>% filter(test$`Allele (length)` == 0 & test$Counts == 0)

test$SampleName %>% unique # weird samples

OhioECV81 <- All %>% filter(All$SampleName == "OhioECV81") # basically empty, BAD
OhioECT159 <- All %>% filter(All$SampleName == "OhioECT159") # Good
OhioECV48 <- All %>% filter(All$SampleName == "OhioECV48") # Good
OhioECT181 <- All %>% filter(All$SampleName == "OhioECT181") # basically empty, BAD
L0408 <- All %>% filter(All$SampleName == "L0408") # Good
PET103 <- All %>% filter(All$SampleName == "PET103") # Good

# Check for other samples with low number of observations :

SampleList <- All$SampleName %>% unique()

checkobservation <- function(Sample){
  table <- All %>% filter(All$SampleName == Sample)
  return(nrow(table))
} 

df <- data.frame(
  SampleName = SampleList,
  Observations = unlist(lapply(SampleList, checkobservation))
) # Just OhioECV81 and OhioECT181 have low observations


## filtering the above out

All <- All[!(All$SampleName %in% c("OhioECV81", "OhioECT181")), ]

# Write CSV

write.csv(All, "Data/All.csv", row.names = FALSE)
