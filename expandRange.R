# Purpose:to write out values between start and end positions and apply chromosome number and signal to each
# input: dataframe from bedgraph (works for both roman numerals and numbers)
# output: expanded, dataframe with chromososome number, all values between "start" and "end", and "signal"
# install.packages("dplyr")

expandRange <- function(bedgrph) {
  library(dplyr)
  y <- list()
  expanded <- data.frame()
  colnames(spcl) <- c("chr","start","end","signal")
  spcl <- spcl[order(spcl$chr),]
  for(i in 1:16) {
    tmp <- spcl[which(spcl$chr==unique(spcl$chr)[i]),]
    u <- list()
    for(l in 1:nrow(tmp)) {
      t <- data.frame()
      tmp1 <- seq(tmp$start[l], tmp$end[l], by=1)
      t[1:length(tmp1),1] <- tmp[l,1]  #add chr name to column 1
      t[1:length(tmp1),2] <- tmp1   #add vector containing sequence of numbers to column 2
      t[1:length(tmp1),3] <- tmp[l,4]   #add signal value to column 3
      u[[l]] <- t
    }
    y[[i]] <- do.call(bind_rows, u)
  }
  expanded <- do.call(bind_rows, y)
  colnames(expanded) <- c("chr","position", "signal")
  return(expanded)
}