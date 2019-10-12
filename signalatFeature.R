# Purpose:to average sequencing or microarray data around a feature (ex nucleosome occupancy around hotspots)
# input: 1) output from expandRange function, dataframe with chromososome number, all values between "start" and "end", and "signal"; 2) feature, dataframe from bed file (columns: chr, start, end); 3) range around centered position
# output: avgCntr, dataframe with averaged signal values in specified range around central value
# saving and plotting not included
# install.packages("dplyr")

signalatFeature <- function(expanded, feature, range) {
  if(missing(range)) {
    range=2000}
  colnames(feature) <- c("chr", "start", "end")
  colnames(expnd) <- c("chr","position", "signal")
  feature <- feature[order(feature$chr),]
  expnd <- expnd[order(expnd$chr),]
  ## to find means of start and stop positions in centering file (spo11 file)
  mid <- (feature$start + feature$end) /2
  midR <- round(mid)
  feature[,"position"] <- c(midR)
  
  z <- data.frame()
  u <- list()
  for(i in 1:length(unique(expnd[,1]))) {
    temp=feature[which(feature$chr==unique(feature$chr)[i]),]
    temp4=expnd[which(expnd$chr==unique(expnd$chr)[i]),]
    w <- list()
    for(j in 1:nrow(temp)) {
      temp5 <- data.frame()
      up <- temp$position[j]-range
      down <- temp$position[j]+range
      temp5 <- temp4[temp4$position>=up & temp4$position<=down,]
      temp5$position <- temp5$position-temp$position[j]
      w[[j]] <- temp5
    }
    u[[i]] <- do.call(dplyr::bind_rows,w)
  }
  z <- do.call(dplyr::bind_rows,u)
  avgCntr=aggregate(z$signal, by=list(z$position), FUN=mean)
  colnames(avgCntr) <- c("position", "signal")
  return(avgCntr)
}
