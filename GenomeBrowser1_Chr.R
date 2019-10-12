# Genome Browser View: Choose chromosome number and binning
# Purpose: To specify which chromosome and what bin size to use to prepare data for plotting.
# input: dataframe from bedgraph (works for both roman numerals and numbers)

genomeView <- function(samp1,chrnum,tile) {
  genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
  samp1 <- sort(GenomeInfoDb::sortSeqlevels(samp1))
  GenomeInfoDb::seqlengths(samp1) <- GenomeInfoDb::seqlengths(genome_info)
  bins_samp1 <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(samp1),
                                          tilewidth=tile,
                                          cut.last.tile.in.chrom=TRUE)
  score_samp1 <- GenomicRanges::coverage(samp1, weight="score")
  bins_samp1 <- GenomicRanges::binnedAverage(bins_samp1, score_samp1, "binned_score")
  bins_samp1 <- GenomeInfoDb::keepSeqlevels(bins_samp1, paste0("chr",chrnum))
  positions_samp1 <- bins_samp1@ranges@start + floor(bins_samp1@ranges@width / 2)
  df_samp1 <- data.frame(position=positions_samp1/1000, signal=bins_samp1$binned_score)
  return(df_samp1)
}