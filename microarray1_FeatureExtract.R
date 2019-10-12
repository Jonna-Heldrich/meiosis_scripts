#ssDNA analysis function #1

## Generate log2 ratios
## Input: Microarray .txt data file
## Example: FE2_R64("FileName", cy0, "strain") #run Cy# based on sample of interest (ie 5h)
## Optional: average 2 log2 files (dye-swap)
## ssDNA_average("file1","file2","strain")

# function FE2 to generate log2 ratios

#Example: FE2_R64("US93503718_251523910360_S01_ChIP_107_Sep09_1_3.txt",5,"118_Cy5_161005.txt") #run Cy# based on control (ie 0h)

FE2_R64<-function(FileName, cy0, strain)
{
  library(marray) 
  plat<-read.table("SK1rosetta_R64_5.txt", header=T, sep="\t")
  maData<-read.Agilent(fnames=FileName, name.Rf="rMeanSignal", name.Gf="gMeanSignal",
                       name.Rb="rBGMeanSignal", name.Gb="gBGMeanSignal", sep="\t")
  maNorm<-maNorm(maData,norm="printTipLoess")
  if (cy0==5)
  {
    lr<--maNorm@maM
  }
  if (cy0==3)
  {
    lr<-maNorm@maM
  }
  data<-as.data.frame(matrix(0,nrow=nrow(plat),ncol=4))
  data[,1]<-plat[,'chr']
  data[,2]<-plat[,'start']
  data[,3]<-plat[,'stop']
  data[,4]<-lr[,1]
  data = data[order(data[,2]),]
  data = data[order(data[,1]),]
  data = data[which(data[,1]!=0),]
  data = data[which(data[,2]!=0),]
  data = data[which(data[,3]!=0),]
  colnames(data)<-c('chr','start', 'end', "Log2Ratio")
  write.table(data,paste0(strain,"_cy",cy0,"_FE2_R64.txt"),col.names=T, row.names=F, quote=F, sep="\t")
  # save the signal data into a txt file in your working directory.
}


##########################################################
##########################################################
# Average 2 log2 files
#Example: ssDNA_average("file1.txt","file2.txt","WT")

ssDNA_average <- function(file1,file2,strain){
  a=read.table(file1,header=TRUE,sep="\t")
  a2=read.table(file2,header=TRUE,sep="\t")

  #change to dataframe
  a=data.frame(a)
  a2=data.frame(a2)
  #name columns
  colnames(a)=c("chr", "start","end","LogRatio")
  colnames(a2)=c("chr", "start","end","LogRatio")

  #average between the two files
  avg=data.frame(matrix(nrow=nrow(a),ncol=2))
  colnames(avg)=c("a","a2")
  avg[,1]=a$LogRatio
  avg[,2]=a2$LogRatio
  avgs=rowMeans(subset(avg, select=c(a,a2)), na.rm = TRUE)
  a$LogRatio=avgs

  #create and output averaged files
  write.table(a,paste0(strain,"_avg.txt"),col.names=T, row.names=F, quote=F, sep="\t")
}
