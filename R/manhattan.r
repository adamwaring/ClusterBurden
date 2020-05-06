#' A Manhattan plot for p-values derived from a ClusterBurden WES analysis
#'
#' @author Adam Waring - adam.waring@@msdtc.ox.ac.uk
#' @return Returns a ggplot object
#'
#' @param pvals pvals data.table from ClusterBurden WES analysis
#' @param test either bin-test, burden or clusterburden (not case sensitive)
#'
#' @keywords RVAT, distribution, cluster, genetics, case-control, gene
#'
#' @details A Manhattan plot for p-values resulting from a ClusterBurden analysis
#'
#' @export

manhattan = function(pvals, test, n_genes, SCALE=1){

  test = c("bin-test"="BIN.test_pvalue", "burden"="burden_pvalue", "clusterburden"="ClusterBurden")[tolower(test)]

  ytest = c("BIN.test_pvalue"="BIN-test", "burden_pvalue"="Fishers-exact", "ClusterBurden"="ClusterBurden")
  ylabel = paste0(ytest[test], " p-value (-log10)")


  point_size = 2 * SCALE
  segment_size = 1.1 * SCALE
  text_size = 6 * SCALE
  text_offset = 50 * SCALE
  theme_size = 25 * SCALE


  data("gene_locations")
  pvals = merge(pvals, gene_locations, all.x=T)

  colnames(pvals)[which(colnames(pvals)==test)] = "pval"
  pvals[,pval:=-log10(pval)]

  pvals[,chrom_color := ifelse(chr %in% c(as.character(seq(1, 22, 2)), "X"), "a", "b")]
  pvals[,chr := factor(chr, levels = c(as.character(1:22), "X"))]
  pvals = pvals[!is.na(pval)][order(chr, pos)]
  pvals[,index:=1:.N]

  maxp = pvals[pval!=Inf, max(pval)]
  pvals[, pval := ifelse(pval==Inf, ifelse(maxp>20, maxp, 20), pval)]

  nps = sum(!is.na(pvals$pval))
  bsignif_threshold = -log10(0.05/nps)

  ymax = ceiling(max(pvals$pval, na.rm=T)/5)*5
  if(ymax < bsignif_threshold) ymax = ceiling(bsignif_threshold/5)*5

  ymin = -(ymax/100)

  midpoints = pvals[,.(mid=mean(index)), by=chr]
  most_signif = pvals[order(-pval)][1:n_genes]

  ggplot(pvals, aes(x = index, y = pval, color = chrom_color)) + geom_point(size=point_size) +
    geom_text(data = midpoints, aes(x = mid, y = ymin*2, label=chr), color="black", size=text_size) +
    geom_segment(x = 1, xend = nrow(pvals), y = -log10(0.05), yend = -log10(0.05), lty="dashed", size=segment_size, color="grey") +
    geom_segment(x = 1, xend = nrow(pvals), y = bsignif_threshold, yend = -log10(0.05/nps), lty="dashed", size=segment_size, color="grey") +
    geom_text(data = most_signif, aes(y = pval + ymax/text_offset, label = symbol), size=text_size) +
    scale_y_continuous(limits = c(ymin*4, ymax), expand = c(0, 0)) +
    scale_x_continuous(limits = c(-10, nps + 10), expand = c(0, 0)) +
    scale_color_manual(values=c("a"="#E69F00", "b"="#009E73")) +
    guides(color=F) + theme_bw(theme_size) +
    labs(x = "Chromosome", y = ylabel) +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.line.y = element_line(),
          panel.border = element_blank())

}



