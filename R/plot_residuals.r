#' Plot for chi-squared residuals
#'
#' This function allows you to plot the residuals of a chi-squared two-sample test to infer bins
#' that differ significantly from expected values.
#'
#' @author Adam Waring - adam.waring@@msdtc.ox.ac.uk
#' @return Returns an object of class ggplot.
#'
#' @param chisq output from chi-squared two-sample test.
#' @keywords RVAT, distribution, cluster, genetics, case-control, gene
#'
#' @details This method plots residuals as a 2 x m table where m is the number of bins in the chisq.test.
#'
#' Residuals are standardized by the standard deviation of the residuals.
#' Within each cell of the table is a circle whose color and size relate to the size of the residual.
#' Residuals that are significant at a nominal 0.05 (two-tailed) after Bonferonni correction are highlighted with an asterix.
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
#' chisq = BIN_test(case_residues, control_residues)
#' plot_residuals(chisq)

plot_residuals = function(chisq){

  contig = chisq$observed
  contig2 = chisq$expected

  resids = chisq$residuals / sd(chisq$residuals)

  n = ncol(resids)
  alpha = 0.025 / n
  threshold = abs(qnorm(alpha))
  signif = ifelse(as.vector(resids) < -threshold | as.vector(resids) > threshold, "*", "")

  bin = factor(rep(colnames(contig), 2), levels=colnames(contig))

  plot_table = data.frame(bin = bin,
                          Residuals = as.vector(t(resids)),
                          cohort=c(rep("Cases OE", n), rep("Controls OE", n)),
                          counts = paste0(c(paste0(contig[1,], "/", contig2[1,]), paste0(contig[2,], "/", contig2[2,])), signif))


  max_resid = max(abs(resids))
  limit = ifelse(max_resid > threshold, max_resid, threshold) + 1

  ggplot(plot_table, aes(x=bin, y=cohort)) + geom_point(aes(size=abs(Residuals), color=Residuals)) +
    geom_text(color="black", aes(label = counts)) +
    scale_color_gradient2(low = scales::muted("blue"), mid = "white", high = scales::muted("red"), limits=c(-limit, limit)) +
    scale_size_continuous(range = c(2, 10), guide="none", limits=c(0, limit/2)) +
    coord_fixed(ratio = 0.5) +
    xlab("\n Bins") + ylab("") +
    guides(colour = guide_colourbar(title="Standardized residuals", title.position="top", title.hjust = 0.5, barwidth = 20)) +
    labs(caption = paste0("\n* denotes two-tailed significance at an alpha of 0.05\n
                            Bonferroni corrected for ", n, " bins.\n
                            Alpha = ", signif(alpha, 2), ". Significant residual at +/- ", round(threshold, 2), "."), size=1) +
    theme(panel.background=element_blank(),
          axis.title = element_text(size=15),
          axis.text = element_text(size=10),
          axis.ticks = element_blank(),
          legend.position = "top",
          legend.title = element_text(size=10),
          plot.caption = element_text(hjust=0.9, size=8))
}
