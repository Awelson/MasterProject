setwd("C:/Users/aaron/Desktop/Github Projects/MasterProject")

library(tidyverse)
library(purrr)
library(magrittr)
library(reshape2)
library(plyr)
library(dplyr)

# Functions

check_string <- function(string) {
  # Split the string at the underscore
  parts <- strsplit(string, "_")[[1]]
  
  # Check if the part after the underscore is exactly "SNP1"
  if (parts[length(parts)] != "SNP1") {
    return(FALSE)
  }
  
  # If we haven't returned FALSE, return TRUE
  return(TRUE)
}

convert_and_replace <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "NA"  # Replace NA with "NA" or any other value you prefer
  factor(x)
}

tempfunction4 <- function(SingleMarkerTable){
  
  SingleMarkerTable <- as.data.frame(SingleMarkerTable)
  
  if( empty(SingleMarkerTable)) return(NA)
  
  SingleMarkerTable[, 1:3] <- lapply(SingleMarkerTable[, 1:3], convert_and_replace)
  
  colnames(SingleMarkerTable) <- c("Allele (length)", "Allele (SNP)", "Counts")
  
  return(SingleMarkerTable)
} 


tempfunction5 <- function(list){
  
  result <- mapply(function(df, name) {
    df$Marker <- name
    return(df)
  }, list, names(list), SIMPLIFY = FALSE)
  
  return(result)
}

tempfunction6 <- function(list){
  
  test <- lapply(list, tempfunction4)
  
  test <- test[!is.na(test)]
  
  test <- tempfunction5(test)
  
  tofilterout <- unlist(lapply(names(test), check_string))
  
  test <- test[tofilterout]
  
  test <- test[!(names(test) %in% c("KRAS_SNP1", "BRAF_SNP1", "NR21_SNP1"))]
  
  df <- bind_rows(test)
  
  return(df)
  
}

ReadFolder3 <- function(path){
  
  FileList <- lapply(list.files(path, full.names=TRUE), function(file) {
    source(file, local = TRUE)
    return(Result)
  })
  
  dflist <- lapply(FileList, tempfunction6)
  
  SampleNames <- sapply(strsplit(list.files(path), "_"), function(x){x[1]})
  
  dflist2 <- mapply(function(df, name) {
    df$SampleName <- name
    return(df)
  }, dflist, SampleNames, SIMPLIFY = FALSE)
  
  df <- bind_rows(dflist2)
  
  return(df)
}