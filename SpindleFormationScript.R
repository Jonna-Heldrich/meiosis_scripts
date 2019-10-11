# Spindle formation plots
# works for 2 or 3 strains
# first run line 22 to import the function

# insert your strain names, time points, and spindle quantification data

# Example:
#
# strain1 <- 'wt'
# strain2 <- 'top1'
# strain3 <- 'top3'
# timepoints<-c(0,3,4,5,6,7,8)
# nospindle_1<-c(200, 194,174,45,11,1,0) #strain1 no spindles
# spindles_1 <- c(0,6,26,155,189,199,200) #strain1 1+2 spindles
# nospindle_2<-c(200, 200, 176,45,3,2,3) #strain2 no spindles
# spindles_2 <- c(0,0,24,155,197,198,197) #strain2 1+2 spindles
# nospindle_3<-c(199,200,200,200,198,191,200) #strain3 no spindles
# spindles_3 <- c(1,0,0,0,2,9,15) #strain3 1+2 spindles

# to plot the data, run the function with your variable names
# example:
# plotspindles(timepoints,strain1,nospindle_1,spindles_1,strain2,nospindle_2,spindles_2,strain3,nospindle_3,spindles_3)
###################

plotspindles <- function(timepoints,strain1,MTOC1,MTOC2,strain2,MTOC3,MTOC4,strain3,MTOC5,MTOC6) {
  if(missing(strain3)){
    sampnum=2 
    MTOC5 <- NULL
    MTOC6 <- NULL
    } else sampnum=3

  spindles1<-data.frame(timepoints,MTOC1,MTOC2,stringsAsFactors=F)
  spindles2<-data.frame(timepoints,MTOC3,MTOC4,stringsAsFactors=F)
  if(sampnum==2){
    spindles3<-NULL
  } else {
    spindles3 <- data.frame(timepoints,MTOC5,MTOC6,stringsAsFactors=F)}

  if (sampnum==2){
    allspindles <- list(spindles1,spindles2)
  } else {
    allspindles <- list(spindles1,spindles2,spindles3)
  } 
  
  spindlesPercent <- list()
  for(j in allspindles) {
    spind <- j
    for(i in 1:nrow(spind)){
      spind[i,2] <- spind[i,2]/sum(spind[i,2:3])*100
      spind[i,3] <- spind[i,3]/sum(spind[i,2:3])*100
    }
    spindlesPercent<- append(spindlesPercent,list(spind))
  }

  par(las=1)
  if(sampnum==2){
  #for 2 strains
    p3 <- plot(x=spindlesPercent[[1]]$timepoints,y=spindlesPercent[[1]]$MTOC2,pch=19, col='blue',type='b', lwd=5, ylim=c(0,100), xlab='Time Points', ylab='Cells with spindles (Percentage)', frame.plot=F);
        legend(0,y=40, legend=c(strain1,strain2), col=c('blue','red'), pch=19);
        points(x=spindlesPercent[[2]]$timepoints,y=spindlesPercent[[2]]$MTOC4, type='b',pch=19, lwd=5, col='red')
  } else if(sampnum==3){
    #for 3 strains
    p3 <- plot(x=spindlesPercent[[1]]$timepoints,y=spindlesPercent[[1]]$MTOC2,pch=19, col='blue',type='b', lwd=5, ylim=c(0,100), xlab='Time Points', ylab='Cell with spindles (Percentage)', frame.plot=F);
        legend(0,y=40, legend=c(strain1,strain2,strain3), col=c('blue','red','green'), pch=19);
        points(x=spindlesPercent[[2]]$timepoints,y=spindlesPercent[[2]]$MTOC4, type='b',pch=19, lwd=5, col='red');
        points(x=spindlesPercent[[3]]$timepoints,y=spindlesPercent[[3]]$MTOC6, type='b',pch=19, lwd=5, col='green')
  }
}





