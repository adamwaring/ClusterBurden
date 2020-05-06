#' Plot the distributions of case and control variant residues across multiple genes
#'
#' @author Adam Waring - adam.waring@@msdtc.ox.ac.uk
#' @return Returns a ggplot object
#'
#' @param pvals pvals data.table from ClusterBurden WES analysis with at least symbol and p-value column for test of interest
#' @param cases case data in format: data.table(aff, symbol, protein_position, ac)
#' @param controls control data in format: data.table(aff, symbol, protein_position, ac)
#' @param case_coverage optional coverage data for cases in format: data.table(symbol, protein_position, over_10)
#' @param control_coverage optional coverage data for controls in format: data.table(symbol, protein_position, over_10)
#' @param test either bin-test, burden or clusterburden (not case sensitive)
#' @param cov_threshold threshold at which to exclude a residue position from the analysis (choose 0 to keep all residues)
#' @param n_genes number of genes to plot e.g. n_genes = 5 plots the 5 most significant genes
#' @param SCALE scale of the plot
#'
#' @keywords RVAT, distribution, cluster, genetics, case-control, gene
#'
#' @details For each gene (n_genes=?): produces a stripchart using geom_jitter for both cases and controls as well as a gaussian density line for each cohort. If coverage files are provided then these are included as a rug underneath and above the case and control variants and are
#' colored by their coverage level (< cov_threshold, < 80% 10X, < 90% 10X, > 90% 10X)
#'
#' @export



plot_signif_distribs = function(pvals, cases, controls, case_coverage, control_coverage, test="BIN-test", cov_threshold=0.5, n_genes=20, SCALE=1){

  test = c("bin-test"="BIN.test_pvalue", "burden"="burden_pvalue", "clusterburden"="ClusterBurden")[tolower(test)]

  point_size = 2 * SCALE
  line_size = 1.5 * SCALE
  cov_size = 3 * SCALE
  theme_size = 25 * SCALE


  top_pvals = pvals[order(pvals[[test]])][1:n_genes]
  top_genes = top_pvals$symbol

  dataset = rbind(cases, controls)[symbol%in%top_genes]
  dataset[,y:=aff+1]

  case_prov = attributes(cases)$provenance
  control_prov = attributes(controls)$provenance

  case_auto = !is.null(case_prov) && case_prov == "auto"
  control_auto = !is.null(control_prov) && control_prov == "auto"

  if(case_auto) case_coverage = attributes(cases)$coverage
  if(control_auto) control_coverage = attributes(controls)$coverage


  null1 = !is.null(case_coverage)
  null2 = !is.null(control_coverage)
  coverage = NULL
  if(null1){
    coverage = rbind(coverage, cbind(case_coverage[symbol%in%top_genes], aff=1, y=3))
  }
  if(null2){
    coverage = rbind(coverage, cbind(control_coverage[symbol%in%top_genes], aff=0, y=0))
  }


  dataset[,color := ifelse(aff==2, "Cases X10 > 80%", "Controls X10 > 80%")]
  COLORS = setNames(c("#D8B70A", "#81A88D"),
                    c("Cases X10 > 80%", "Controls X10 > 80%"))

  if(null1|null2){
    blackname = paste0("Blacklisted (X10 < ", cov_threshold, ")")
    coverage[,color := ifelse(over_10 < cov_threshold, blackname, ifelse(over_10 < 0.8, "X10 < 80%", ifelse(over_10 < 0.9, "X10 < 90%", ifelse(aff==1, "Cases X10 > 80%", "Controls X10 > 80%"))))]
    COLORS = c(COLORS, setNames(c("black", "darkred", "red"),
                                c(blackname, "X10 < 80%",  "X10 < 90%")))

    dataset = merge(dataset, coverage[,.(symbol, protein_position, aff, over_10)], by=c("symbol", "protein_position", "aff"), all.x=T)
    dataset = merge(dataset, unique(coverage[color == blackname,.(symbol, protein_position, black=over_10)]), by=c("symbol", "protein_position"), all.x=T)
    dataset[,color := ifelse(!is.na(black), blackname, ifelse(over_10 < 0.8, "X10 < 80%", ifelse(over_10 < 0.9, "X10 < 90%", ifelse(aff==1, "Cases X10 > 80%", "Controls X10 > 80%"))))]
  }


  g1 = ggplot(dataset, aes(y = y, x = protein_position, color=color)) +

    geom_jitter(alpha=0.3, size=point_size) +

    facet_wrap(~symbol, scales = "free_x") +

    stat_density(geom="line", data=dataset[aff==1], aes(y=1.5 + ..scaled..), size=line_size, color="#D8B70A") +
    stat_density(geom="line", data=dataset[aff==0], aes(y=0.5 + ..scaled..), size=line_size, color="#81A88D") +

    scale_color_manual(name="", values=COLORS) +

    labs(x = "Index in linear protein sequence", y = "") +

    theme_bw(theme_size) + theme(legend.position="top")

  if(null1|null2){
    g1 = g1 +  geom_point(data=coverage, aes(color=color), size=cov_size, shape=15) +
      geom_point(data=coverage[color=="X10 < 80%"], aes(color=color), size=cov_size, shape=15) +
      geom_point(data=coverage[color=="X10 < 90%"], aes(color=color), size=cov_size, shape=15) +
      geom_point(data=coverage[color==blackname], aes(color=color), size=cov_size, shape=15)
  }

  if(null1 & null2){

    g1 = g1 + scale_y_continuous(breaks = 0:3, limits = c(-0.25, 3.25), labels = c("Control coverage", "Control variants", "Case variants", "Case coverage"))

  } else if(null2 & !null1){

    g1 = g1 + scale_y_continuous(breaks = 0:2, limits = c(-0.25, 2.5), labels = c("Control coverage", "Control variants", "Case variants"))

  } else if(null1 & !null2){

    g1 = g1 + scale_y_continuous(breaks = 1:3, limits = c(0.5, 3.25), labels = c("Control variants", "Case variants", "Case coverage"))

  } else {

    g1 = g1 + scale_y_continuous(breaks = 1:2, limits = c(0.5, 2.5), labels = c("Control variants", "Case variants"))

  }

  g1
  return(g1)


}
