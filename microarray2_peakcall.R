#ssDNA analysis function #2

## Function to call peaks (highly enriched regions)
## ssDNA_peakcall("log2file","sampleName")

# Example: ssDNA_peakcall("top1_avg.txt","top1D")

ssDNA_peakcall <- function(log2file,sampleName) {
  source('/Volumes/LabShare/HTGenomics/Microarray_database/array_database/Peak_call_microarray.R')
  fe=read.table(log2file,header=T)
  
  #in order to use the peak_call function the signal must be changed from log2
  sig=2^fe[,4]
  fe[,4]=sig
  #save file for peak_call("input","output file")
  peak_call(fe, paste0(sampleName,"_peaks.bed"))
  
  # chromosome names must be changed beause the microarray data is mapped to S288C so the peak files must use roman numerals in order for IGV to read them correctly
  
  fe=read.table(paste0(sampleName,"_peaks.bed"), header=F, sep="\t")
  for(i in 1:16){
    out=c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI")
    fe[which(fe[,1]==i),1] = out[i]
  }
  write.table(fe,file=paste0(sampleName,"_peaks.bed"),quote=F,sep="\t",row.names=F,col.names=F)
  
}
