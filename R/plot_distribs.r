#' Plot the distributions of case and control variant residues
#'
#' @author Adam Waring - adam.waring@@msdtc.ox.ac.uk
#' @return Returns a ggplot object
#'
#' @param case_residues vector of case variant residue positions
#' @param control_residues vector of control variant residue positions
#' @param case_coverage optional coverage data for cases in format: data.table(symbol, protein_position, over_10)
#' @param control_coverage optional coverage data for controls in format: data.table(symbol, protein_position, over_10)
#' @param cov_threshold threshold at which a residue position is excluded from the analysis
#'
#' @keywords RVAT, distribution, cluster, genetics, case-control, gene
#'
#' @details Produces a stripchart using geom_jitter for both cases and controls as well as a gaussian density line for each cohort. If coverage files are provided then these are included as a rug underneath and above the case and control variants and are
#' colored by their coverage level (< cov_threshold, < 80% 10X, < 90% 10X, > 90% 10X)
#'
#' @export
#' @examples
#' nresidues = 1000 # length of the protein
#' probs = rexp(nresidues)^2
#' case_probs = probs * rep(c(1, 3, 1), c(200, 200, 600))
#' control_probs = probs * rep(c(2, 1, 2), c(200, 200, 600))
#'
#' case_residues = sample(1:nresidues, 100, rep=T, case_probs)
#' control_residues = sample(1:nresidues, 100, rep=T, control_probs)
#'
#' plot_distribs(case_residues, control_residues)


plot_distribs = function(case_residues, control_residues, case_coverage=NULL, control_coverage=NULL, cov_threshold=0.5){

  dataset = data.table(aff = c(rep(1, length(case_residues)), rep(0, length(control_residues))),
                       protein_position = c(case_residues, control_residues))
  dataset[,y:=aff+1]

  null1 = !is.null(case_coverage)
  null2 = !is.null(control_coverage)
  coverage = NULL
  if(null1){
    coverage = rbind(coverage, cbind(case_coverage, aff=1, y=3))
  }
  if(null2){
    coverage = rbind(coverage, cbind(control_coverage, aff=0, y=0))
  }


  dataset[,color := ifelse(aff==1, "Cases X10 > 80%", "Controls X10 > 80%")]
  COLORS = setNames(c("#D8B70A", "#81A88D"),
                    c("Cases X10 > 80%", "Controls X10 > 80%"))

  if(null1|null2){
    blackname = paste0("Blacklisted (X10 < ", cov_threshold, ")")
    coverage[,color := ifelse(over_10 < cov_threshold, blackname, ifelse(over_10 < 0.8, "X10 < 80%", ifelse(over_10 < 0.9, "X10 < 90%", ifelse(aff==1, "Cases X10 > 80%", "Controls X10 > 80%"))))]
    COLORS = c(COLORS, setNames(c("black", "darkred", "red"),
                                c(blackname, "X10 < 80%",  "X10 < 90%")))

    dataset = merge(dataset, coverage[,.(symbol, protein_position, aff, over_10)], by=c("protein_position", "aff"), all.x=T)
    dataset = merge(dataset, unique(coverage[color == blackname,.(symbol, protein_position, black=over_10)]), by=c("protein_position"), all.x=T)
    dataset[,color := ifelse(!is.na(black), blackname, ifelse(over_10 < 0.8, "X10 < 80%", ifelse(over_10 < 0.9, "X10 < 90%", ifelse(aff==1, "Cases X10 > 80%", "Controls X10 > 80%"))))]
  }


  g1 = ggplot(dataset, aes(y = y, x = protein_position, color=color)) +

    geom_jitter(alpha=0.3, size=2) +

    stat_density(geom="line", data=dataset[aff==1], aes(y=1.5 + ..scaled..), size=1.5, color="#D8B70A") +
    stat_density(geom="line", data=dataset[aff==0], aes(y=0.5 + ..scaled..), size=1.5, color="#81A88D") +

    scale_color_manual(name="", values=COLORS) +

    guides(color=F) + labs(x = "Index in linear protein sequence", y = "") +

    theme_bw(15) + theme(legend.position="top")

  if(null1|null2){

    g1 = g1 +  geom_point(data=coverage, aes(color=color), size=3, shape=15) +
      geom_point(data=coverage[color=="X10 < 80%"], aes(color=color), size=3, shape=15) +
      geom_point(data=coverage[color=="X10 < 90%"], aes(color=color), size=3, shape=15) +
      geom_point(data=coverage[color==blackname], aes(color=color), size=3, shape=15)

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

  return(g1)


}
