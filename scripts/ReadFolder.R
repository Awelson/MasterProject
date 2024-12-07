setwd("C:/Users/aaron/Desktop/Github Projects/MasterProject")

library(tidyverse)
library(purrr)
library(magrittr)
library(reshape2)
library(plyr)
library(dplyr)

tempfunction1 <- function(SingleMarkerTable){
  
  SingleMarkerTable <- as.data.frame(SingleMarkerTable)
  
  if( empty(SingleMarkerTable)) return(NA)
  
  #returns column number of the NA column
  nacol <- which(is.na(colnames(SingleMarkerTable)))
  
  #renames the NA column to "NA" (string)
  colnames(SingleMarkerTable)[nacol] <- "NA"
  
  SingleMarkerTable$"Allele (length)" <- rownames(SingleMarkerTable)
  rownames(SingleMarkerTable) <- NULL
  
  melted <- melt(SingleMarkerTable,
                 id.vars = "Allele (length)",
                 variable.name = "Allele (SNP)",
                 value.name = "Counts"
  )
  
  return(melted)
}

tempfunction2 <- function(list){
  
  result <- mapply(function(df, name) {
    df$Marker <- name
    return(df)
  }, list, names(list), SIMPLIFY = FALSE)
  
  return(result)
}

tempfunction3 <- function(list){
  list2 <- lapply(list, tempfunction1)
  list3 <- list2[!is.na(list2)]
  list4 <- tempfunction2(list3)
  df <- bind_rows(list4)
  return(df)
}

ReadFolder <- function(path){
  
  FileList <- lapply(list.files(path, full.names=TRUE), function(file) {
    source(file, local = TRUE)
    return(Result)
  })
  
  dflist <- lapply(FileList, tempfunction3)
  
  SampleNames <- sapply(strsplit(list.files(path), "_"), function(x){x[1]})
  
  dflist2 <- mapply(function(df, name) {
    df$SampleName <- name
    return(df)
  }, dflist, SampleNames, SIMPLIFY = FALSE)
  
  df <- bind_rows(dflist2)
  
  return(df)
}

ReadFolder2 <- function(path){
  
  x <- lapply(list.files(path, full.names=TRUE), ReadFolder)
  
  D0 <- x[[1]]
  D1 <- x[[2]]
  D2 <- x[[3]]
  
  D0$UMI <- "D0"
  D1$UMI <- "D1"
  D2$UMI <- "D2"
  
  df <- bind_rows(D0, D1, D2)
  
  return(df)
}

tempfunction4 <- function(sample){
  x <- source(sample)$value
  return(24 - length(x))
}

MissingTables1 <- function(path){
  
  SampleNames <- sapply(strsplit(list.files(path), "_"), function(x){x[1]})
  
  MissingTables <- unlist(lapply(list.files(path, full.names = TRUE), tempfunction4))
  
  df <- data.frame(SampleName = SampleNames, MissingTable = MissingTables)

  return(df)
}

MissingTables2 <- function(path){
  
  x <- lapply(list.files(path, full.names=TRUE), MissingTables1)
  
  D0 <- x[[1]]
  D1 <- x[[2]]
  D2 <- x[[3]]
  
  D0$UMI <- "D0"
  D1$UMI <- "D1"
  D2$UMI <- "D2"
  
  df <- bind_rows(D0, D1, D2)
  
  df <- df %>% filter(MissingTable >= 12)
  
  return(df)
}