FilePathD0 <- list.files("RawData/EC/ManUMI/20220622/D0", full.names=TRUE)
FilePathD1 <- list.files("RawData/EC/ManUMI/20220622/D1", full.names=TRUE)
FilePathD2 <- list.files("RawData/EC/ManUMI/20220622/D2", full.names=TRUE)

NewNames <- list.files("RawData/EC/Man/20220425.24", full.names=FALSE)
NewPathD0 <- file.path(dirname(FilePathD0), NewNames)
NewPathD1 <- file.path(dirname(FilePathD1), NewNames)
NewPathD2 <- file.path(dirname(FilePathD2), NewNames)

file.rename(FilePathD0, NewPathD0)
file.rename(FilePathD1, NewPathD1)
file.rename(FilePathD2, NewPathD2)