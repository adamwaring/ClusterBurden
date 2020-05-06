#' Plot sequence features such as domains alongside distributions of case and control variant residues
#'
#' @author Adam Waring - adam.waring@@msdtc.ox.ac.uk
#' @return Returns a ggplot object
#'
#' @param pvals pvals data.table from ClusterBurden WES analysis
#' @param genename gene to plot
#' @param cases case data in format: data.table(aff, symbol, protein_position, ac)
#' @param controls control data in format: data.table(aff, symbol, protein_position, ac)
#' @param case_coverage optional coverage data for cases in format: data.table(symbol, protein_position, over_10)
#' @param control_coverage optional coverage data for controls in format: data.table(symbol, protein_position, over_10)
#' @param cov_threshold threshold at which to exclude a residue position from the analysis (choose 0 to keep all residues)
#' @param features which sequence features to plot e.g. c("Domain", "Region", "Structure")
#' @param SCALE scale of the plot
#'
#' @keywords RVAT, distribution, cluster, genetics, case-control, gene
#'
#' @export


plot_features = function(genename, pvals, cases, controls, case_coverage, control_coverage, cov_threshold=0.5,
                         features = c("Compositional-bias","Coiled-coil","Repeat","Zinc-Finger","Region","Motif","Domain","Topological domain","Intramembrane","Transmembrane","Cross-link","Disulphide","Glycosylation",
                                      "Lipidation","Peptide","Propeptide","Signal peptide","Transit peptide","Structure"),
                         SCALE=1){


  pos1 = 2.1 # case variants
  pos2 = 1.8 # control variants
  pos3 = 2.6 # case coverage
  pos4 = 1.3 # control coverage


  top_pvals = pvals[symbol==genename]

  dataset = rbind(cases, controls)[symbol%in%genename]
  dataset[,aff:=ifelse(aff==1, pos1, pos2)]


  case_prov = attributes(cases)$provenance
  control_prov = attributes(controls)$provenance

  case_auto = !is.null(case_prov) && case_prov == "auto"
  control_auto = !is.null(control_prov) && control_prov == "auto"

  if(case_auto) case_coverage = attributes(cases)$coverage
  if(control_auto) control_coverage = attributes(controls)$coverage


  coverage = NULL
  if(!is.null(case_coverage)){
    coverage = rbind(coverage, cbind(case_coverage[symbol%in%genename], aff=pos3))
  }
  if(!is.null(control_coverage)){
    coverage = rbind(coverage, cbind(control_coverage[symbol%in%genename], aff=pos4))
  }

  blackname = paste0("Blacklisted (X10 < ", cov_threshold, ")")
  dataset[,color := ifelse(aff==1.9, "Case X10 > 90%", "Control X10 > 90%")]
  coverage[,color := ifelse(over_10 < cov_threshold, blackname, ifelse(over_10 < 0.8, "X10 < 80%", ifelse(over_10 < 0.9, "X10 < 90%", ifelse(aff==pos3, "Case X10 > 90%", "Control X10 > 90%"))))]

  COLORS = setNames(c("#D8B70A", "#81A88D", "black", "darkred", "red"),
                    c("Case X10 > 90%", "Control X10 > 90%", blackname, "X10 < 80%",  "X10 < 90%"))



  data("seq_features")
  data("nresidues")

  seq_features = seq_features[symbol %in% genename & feature %in% features]

  seq_features[,color:=feature_name]
  nfeatures = length(unique(seq_features$color))

  # which number feature is it for features of the same type
  seq_features[,ifeature:=1:.N, by=.(symbol, feature, feature_name)]
  seq_features = seq_features[,.(protein_position=start:end), by=.(symbol, feature, feature_name, color, ifeature)]

  # how many lines in plot
  seq_features = seq_features[,nline:=sapply(feature, function(x) which(unique(feature)==x)), by=symbol]
  seq_features[,aff:=1 - (nline * 0.5)]

  findMidBigRange = function(x){

    if(length(x)==1) return(as.double(x))

    starts = x[1]
    stops = NULL
    for(i in 2:length(x)){
      if(x[i] != x[i-1] + 1){
        starts = c(starts, x[i])
        stops = c(stops, x[i-1])
      }
    }
    stops = c(stops, x[length(x)])
    x = data.table(starts, stops)
    x[which.max(stops-starts), mean(c(starts, stops))]
  }

  centres = seq_features[,.(x=findMidBigRange(protein_position)), by=.(symbol, feature, feature_name, color, aff)]
  centres = merge(centres, nresidues, by="symbol")
  centres[,loc:=(x/nresidues)-(nchar(feature_name)/2)*0.01]
  centres[,x:=ifelse(loc<0, ((x/nresidues)+abs(loc))*nresidues, x)]
  centres[order(x),aff:=ifelse(1:.N%%2==1, aff+0.15, aff-0.15), by=.(feature, symbol)]

  findMidRanges = function(x){

    if(length(x)==1) return(as.double(x))

    starts = x[1]
    stops = NULL
    for(i in 2:length(x)){
      if(x[i] != x[i-1] + 1){
        starts = c(starts, x[i])
        stops = c(stops, x[i-1])
      }
    }
    stops = c(stops, x[length(x)])
    x = data.table(starts, stops)
    apply(x, 1, mean)
  }

  findXindex = function(cov_centres){

    INDEX = NULL

    sample1 = function(x){
      x = which(x)
      if(length(x) == 1) return(x) else return(sample(x, size=1))
    }


    # cases
    if(cov_centres[aff==pos3, length(unique(color))]==1){

      index = cov_centres[,which(aff==pos3)]
      INDEX = c(INDEX, index)

    } else {

      min = round(optimise(function(i){
        set.seed(i)

        index = cov_centres[aff==pos3, sample1(cov_centres$color == color & cov_centres$aff == pos3), by=.(color)]$V1
        x = cov_centres[index, x]

        1/min(as.numeric(dist(x)))
      }, interval = 1:1000)$minimum)


      set.seed(min)
      INDEX = c(INDEX, cov_centres[aff==pos3, sample1(cov_centres$color == color & cov_centres$aff == pos3), by=.(color)]$V1)


    }

    if(cov_centres[aff!=pos3, length(unique(color))]==1){

      index = cov_centres[,which(aff!=pos3)]
      INDEX = c(INDEX, index)

    } else {

      min = round(optimise(function(i){
        set.seed(i)

        index = cov_centres[aff!=pos3,sample1(cov_centres$color == color & cov_centres$aff != pos3), by=.(color)]$V1
        x = cov_centres[index, x]

        1/min(as.numeric(dist(x)))
      }, interval = 1:1000)$minimum)


      set.seed(min)
      INDEX = c(INDEX, cov_centres[aff!=pos3,sample1(cov_centres$color == color & cov_centres$aff != pos3), by=.(color)]$V1)

    }

    return(INDEX)

  }


  cov_centres = coverage[,.(x=findMidRanges(protein_position)), by=.(symbol, color, aff)]
  cov_centres = cov_centres[findXindex(cov_centres)]
  cov_centres = merge(cov_centres, nresidues, by="symbol")
  cov_centres[,loc:=(x/nresidues)-(nchar(as.character(color))/2)*0.01]
  cov_centres[,x:=ifelse(loc<0, ((x/nresidues)+abs(loc))*nresidues, x)]
  cov_centres[order(x),y:=ifelse(1:.N%%2==1, aff+0.1, aff-0.1), by=aff]


  point_size = 3 * SCALE
  line_size = 2 * SCALE
  cov_size = 5 * SCALE
  text_size = 5 * SCALE
  theme_size = 15 * SCALE
  title_size = 15 * SCALE

  data("color_vector")
  feature_names = unique(seq_features$feature_name)
  COLORS = c(COLORS, setNames(color_vector[1:length(feature_names)], feature_names))

  ymin = min(seq_features$aff) - 0.25


  ggplot(dataset, aes(y = aff, x = protein_position, color=color)) +

    geom_jitter(data=dataset[aff==pos1], aes(size=ac * SCALE), alpha=0.2, height=0.1, color = "#D8B70A") +
    geom_jitter(data=dataset[aff==pos2], aes(size=ac * SCALE), alpha=0.2, height=0.1, color = "#81A88D") +

    stat_density(geom="line", data=dataset[aff==pos1], aes(y=pos1 - 0.2 + ..scaled..*(1/3)), size=line_size, color="#D8B70A") +
    stat_density(geom="line", data=dataset[aff==pos2], aes(y=pos2 - 0.2 + ..scaled..*(1/3)), size=line_size, color="#81A88D") +

    geom_point(data=coverage[aff==pos3], size=point_size, shape=15, color = "#D8B70A") +
    geom_point(data=coverage[aff==pos4], size=point_size, shape=15, color = "#81A88D") +

    geom_point(data=coverage[color=="X10 < 80%"], size=point_size, shape=15, color="darkred") +
    geom_point(data=coverage[color=="X10 < 90%"], size=point_size, shape=15, color="red") +
    geom_point(data=coverage[color==blackname], size=point_size, shape=15, color="black") +


    geom_text(data=cov_centres[aff==pos3], aes(x = x, y = y, label = color), size=text_size, color = "#D8B70A") +
    geom_text(data=cov_centres[aff==pos4], aes(x = x,  y = y, label = color), size=text_size, color = "#81A88D") +

    geom_text(data=cov_centres[color=="X10 < 80%"], aes(x = x,  y = y, label = color), size=cov_size, color="darkred") +
    geom_text(data=cov_centres[color=="X10 < 90%"], aes(x = x,  y = y, label = color), size=cov_size, color="red") +
    geom_text(data=cov_centres[color==blackname], aes(x = x,  y = y, label = color), size=cov_size, color="black") +

    #geom_text(data=cov_centres, aes(x = x, label = color), size=5) +

    geom_point(data=seq_features, size=point_size, shape=15) +
    geom_text(data=centres, aes(x = x, label=feature_name), size=text_size) +

    scale_color_manual(name="", values=COLORS) +

    scale_y_continuous(breaks = c(sort(unique(seq_features$aff)), pos4, pos2, pos1, pos3), limits = c(ymin, pos3+0.1), labels = c(seq_features[order(aff), as.character(unique(feature))], "Control coverage", "Control variants", "Case variants", "Case coverage")) +

    guides(size=F, color=F) + labs(x = "Index in linear protein sequence", y = "") +

    ggtitle(with(top_pvals, paste0(symbol, "        (BIN-test p < ", signif(BIN.test_pvalue, 3), ",    Burden p < ", signif(burden_pvalue, 3), ",    ncase_vars = ", ncases, ",    ncontrol_vars = ", ncontrols, ")"))) +

    theme_bw(theme_size) + theme(panel.grid = element_blank(),
                                 plot.title = element_text(size=title_size))


}
