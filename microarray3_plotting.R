#ssDNA analysis function #3

## Draw plots for to analyze signal across all chromosomes
## ssDNA_plots("wt","wtpeak","mutant","mutantpeak","mutantname")

#Example: ssDNA_plots("WT_avg.txt","WT_avg_peaks.bed","top1.txt","top1_avg_peaks.bed","top1D")

ssDNA_plots <- function(wt,wtpeak,mutant,mutantpeak,mutantname) {
  #open microarray txt files
  a=read.table(wt,header=TRUE,sep="\t")
  b=read.table(mutant,header=TRUE,sep="\t")
  a=data.frame(a)
  b=data.frame(b)
  colnames(a)=c("chr", "start","end","LogRatio")
  colnames(b)=c("chr", "start","end","LogRatio")
  
  #open peak files
  wtpeaks=read.table(wtpeak, header=F)
  colnames(wtpeaks)=c("chr", "start","end","peak","hits")
  mtpeaks=read.table(mutantpeak, header=F)
  colnames(mtpeaks)=c("chr", "start","end","peak","hits")
  spo=read.table("/Volumes/LabShare/HTGenomics/Microarray_database/array_database/spo11_map.bed", header=F)
  colnames(spo)=c("chr", "start","end","peak","hits")
  max_spo=data.frame()
  max_spo2=data.frame()
  
  for(i in 1:16){
    chrom=c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI")
    row=which(wtpeaks$chr==chrom[i])
    #wt microarray data for this chrom
    temp=wtpeaks[row,]
    row2=which(mtpeaks$chr==chrom[i])
    #mutant microarray data for this chrom
    temp2=mtpeaks[row2,]
    
    #spo11 data for mutant for this chrom
    row=which(spo$chr==chrom[i])
    temps=spo[row,]
    for (j in 1:nrow(temp)){
      #find max spo11 oligo hits in each peak, extract this position for use
      start=temp[j,2]
      end=temp[j,3]
      spo_peaks=temps[which(temps$start>=start&temps$start<=end),]
      row.max=which.max(spo_peaks[,5])
      max_spo=rbind(max_spo,spo_peaks[row.max,])
      
    }
    
    for (j in 1:nrow(temp2)){
      #find max spo11 oligo hits in each peak, extract this position for use
      start=temp2[j,2]
      end=temp2[j,3]
      spo_peaks=temps[which(temps$start>=start&temps$start<=end),]
      row.max=which.max(spo_peaks[,5])
      max_spo2=rbind(max_spo2,spo_peaks[row.max,])
    }
  }
  
  for( i in 1:16){
    rows=which(a$chr==i)
    aa=a[rows,]
    #order by start position
    aa=aa[order(aa[,2]),]
    x1=aa[,2]/1000
    y1=aa[,4]
    
    rows=which(b$chr==i)
    bb=b[rows,]
    #order by start position
    bb=bb[order(bb[,2]),]
    x2=bb[,2]/1000
    y2=bb[,4]
    
    #sliding window of 3
    window=3
    step=1
    slideFunct <- function(data, window, step){
      total <- length(data)
      spots <- seq(from = 1, to = (total - window), by = step)
      result <- vector(length = length(spots))
      for(i in 1:length(spots)){
        result[i] <- mean(data[spots[i]:(spots[i] + window - 1)])
      }
      return(result)
    }
    
    data=y1
    y1=slideFunct(data, window, step)
    
    data=x1
    x1=slideFunct(data, window, step)
    
    data=x2
    x2=slideFunct(data, window, step)
    
    data=y2
    y2=slideFunct(data, window, step)
    
    peaks=wtpeaks[which(wtpeaks[,1]==chrom[i]),]
    peak_wt=(peaks[,2]+((peaks[,3]-peaks[,2])/2))/1000
    peaks=mtpeaks[which(mtpeaks[,1]==chrom[i]),]
    peak_mt=(peaks[,2]+((peaks[,3]-peaks[,2])/2))/1000
    
    
    x=paste0("position on chr ",i," (kb)")
    lab=paste0("WTvs",mutantname,"_1chrom",i,".pdf")
    pdf(lab,width=24,height=8)
    plot(x1,2^y1,col="blue",xlab=x,ylab="ssDNA enrichment",cex.lab=1.5,cex.axis=1.5,pch=20,ylim=c(0,5),type="l")
    points(x2,2^y2,col='red',pch=20, type="l")
    points(x=peak_wt,y=rep(5,length(peak_wt)),col="blue",pch=25,bg="blue",cex=1.5)
    points(x=peak_mt,y=rep(4.5,length(peak_mt)),col="red",pch=25,bg="red",cex=1.5)
    legend("bottomright",legend=c("WT",mutantname),cex=1.5,lty=c(1,1),lwd=2,col=c("blue","red"))
    dev.off()
    
    if(i==(1||3||6)){move=1}else{move=7}
    #     max_spo[which(max_spo[,1]==chrom[i]),2]/1000
    lab2=paste0("WTvs",mutantname,"_2chrom",i,".pdf")
    pdf(lab2)
    # par(mfrow = c(4, 1))
    par(mfrow = c(4, 1))
    for(i in 1:4){
      z=i*(length(x1)/4)
      w=((i-1)*(length(x1)/4))+1
      lab=c("","","",x)
      ax=c("n","n","n","s")
      bottom=c(0.1,1.7,3.4,5.1)
      top=c(5.1,3.5,1.8,0.1)
      par(mar=c(bottom[i],4.1,top[i],2.1))
      
      plot(x1[w:z],2^y1[w:z],col="cyan3",xlab=lab[i],cex=1.5,ylab="ssDNA enrichment",pch=20,ylim=c(0,6.5),xlim=c(x1[w]+move,x1[z]-move),type="l")
      points(x2[w:z],2^y2[w:z],col='deeppink',pch=20, type="l")
      points(x=peak_wt,y=rep(6,length(peak_wt)),col="cyan3",pch=25,bg="cyan3")
      points(x=peak_mt,y=rep(6.5,length(peak_mt)),col="deeppink",pch=25,bg="deeppink")
      legend("topright",legend=c("WT",mutantname),lty=c(1,1),lwd=2,col=c("cyan3","deeppink"))
    }
    dev.off()
    
  }
}
