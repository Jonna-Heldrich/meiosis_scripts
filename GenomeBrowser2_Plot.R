# Genome Browser View plotting
# Purpose: Plot whole chromosome or region of chromosome as genome browser view 
# input: 2 outputs from GenomeBrowser1_Chr.R

plot_genomeview <- function(df_samp1,df_samp2,position1,position2,name1,name2,chrnum,color1,color2) {
  par(las=1)
  par(mfrow=c(2,1))
  ax=c(min(df_samp1[df_samp1$position>=position1 & df_samp1$position<=position2,1]),df_samp1[df_samp1$position>=position1 & df_samp1$position<=position2,1],max(df_samp1[df_samp1$position>=position1 & df_samp1$position<=position2,1]))
  ay=c(0,df_samp1[df_samp1$position>=position1 & df_samp1$position<=position2,2],0)
  bx=c(min(df_samp2[df_samp2$position>=position1 & df_samp2$position<=position2,1]),df_samp2[df_samp2$position>=position1 & df_samp2$position<=position2,1],max(df_samp2[df_samp2$position>=position1 & df_samp2$position<=position2,1]))
  by=c(0,df_samp2[df_samp2$position>=position1 & df_samp2$position<=position2,2],0)
  plot(df_samp1[df_samp1$position>=position1 & df_samp1$position<=position2,], xlab=paste0('Position on chromosome ',chrnum,' (kb)'), ylab=name1, type='h',col=color1,frame.plot=F)
  polygon(ax,ay,col=color1,border = NA)
  plot(df_samp2[df_samp2$position>=position1 & df_samp2$position<=position2,], xlab=paste0('Position on chromosome ',chrnum,' (kb)'), ylab=name2, type='h',col=color2,frame.plot=F)
  polygon(bx,by,col=color2,border = NA)
  par(mfrow=c(1,1))
  
}