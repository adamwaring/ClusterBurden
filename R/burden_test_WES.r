#' Calculate Fisher's-exact p-values across the whole-exome
#'
#' @author Adam Waring - adam.waring@@msdtc.ox.ac.uk
#' @return Returns a data.table: symbol, p-value, covstats (optional)
#'
#' @param cases case data in format: data.table(aff, symbol, protein_position, ac)
#' @param controls control data in format: data.table(aff, symbol, protein_position, ac)
#' @param case_coverage optional coverage data for cases in format: data.table(symbol, protein_position, over_10)
#' @param control_coverage optional coverage data for controls in format: data.table(symbol, protein_position, over_10)
#' @param cases_ss sample sizes for cases either a scalar for all genes or a data.table with two columns: symbol, ss
#' @param controls_ss sample sizes for controls either a scalar for all genes or a data.table with two columns: symbol, ss
#' @param cov_threshold threshold at which to exclude a residue position from the analysis (choose 0 to keep all residues)
#' @param alternative seearch for excess variant counts in cases only? (or "two-sided")
#' @param covstats include additional coverage information as a column in results?
#' @param messages print messages to the terminal (e.g. how many variants removed by coverage)
#'
#' @keywords RVAT, distribution, cluster, genetics, case-control, gene
#'
#' @details Calculates a Fisher's-exact test with coverage control across an exome set. Yates continuity correction used for 0 counts. Residue positions with mean 10X coverage less than "cov_threshold" excluded from the analysis.
#' @export
#' @examples
#' library(ClusterBurden)
#'
#' N = 10000 # dataset size
#' PL = 1000 # protein length
#' n_genes = 10
#' ss = 10000 # sample size
#' genes = paste0("gene", 1:n_genes) #gene names
#'
#' # example case dataset
#' cases = data.table(aff=1,
#'                    symbol = sample(genes, N, rep=T),
#'                    protein_position = sample(PL, N, rep=T),
#'                    ac = ceiling(rexp(N, 1.2)))
#'
#' # example control dataset
#' controls = data.table(aff=0,
#'                       symbol = sample(genes, N, rep=T),
#'                       protein_position = sample(PL, N, rep=T),
#'                       ac = ceiling(rexp(N, 1.2)))
#'
#' # provide cases and control sample sizes
#' cases_ss = controls_ss = ss
#'
#' # p-values without coverage control
#' pvals = burden_test_WES(cases, controls, cases_ss, controls_ss)
#'
#' # example coverage files
#' case_coverage = data.table(symbol=rep(genes, each=PL),
#'                            protein_position=rep(1:PL, n_genes),
#'                            over_10=1-rexp(n_genes*PL, 50))
#'
#' control_coverage = data.table(symbol=rep(genes, each=PL),
#'                               protein_position=rep(1:PL, n_genes),
#'                               over_10=1-rexp(n_genes*PL, 50))
#'
#' # p-values with coverage control
#' pvals = burden_test_WES(cases, controls, cases_ss, controls_ss, case_coverage, control_coverage)
#'
#' # example format for genes with different sample sizes
#' cases_ss = controls_ss = data.table(symbol=genes, ss=c(rep(ss, 9), ss/2))
#'
#' # p-values when genes have different sample sizes
#' pvals = burden_test_WES(cases, controls, cases_ss, controls_ss, case_coverage, control_coverage)


burden_test_WES = function(cases, controls, cases_ss=NULL, controls_ss=NULL, case_coverage=NULL, control_coverage=NULL, cov_threshold=0.5, alternative = "greater", covstats=F, messages=T){

  null1 = !is.null(case_coverage)
  null2 = !is.null(control_coverage)

  is1 = is.data.table(cases_ss)
  is2 = is.data.table(controls_ss)

  if(is1 & !setequal(colnames(cases_ss), c("symbol", "ss"))){
    stop("If a data.table  is provided for 'cases_ss' it must have two columns: symbol, ss. Where symbol is the genename and effective sample size for each case gene.")

  }
  if(is2 & !setequal(colnames(controls_ss), c("symbol", "ss"))){
    stop("If a data.table  is provided for 'controls_ss' it must have two columns: symbol, ss. Where symbol is the genename and effective sample size for each control gene.")

  }


  # if autocases or autocontrols automatically supply case and control coverage
  case_prov = attributes(cases)$provenance
  control_prov = attributes(controls)$provenance

  case_auto = !is.null(case_prov) && case_prov == "auto"
  control_auto = !is.null(control_prov) && control_prov == "auto"

  if(case_auto) case_coverage = attributes(cases)$coverage
  if(control_auto) control_coverage = attributes(controls)$coverage

  if(!case_auto & is.null(cases_ss)) stop("Must provide 'cases_ss' as either a scalar sample size for the whole group or a data.table with sample sizes per gene (symbol, ss)")
  if(!control_auto & is.null(controls_ss)) stop("Must provide 'controls_ss' as either a scalar sample size for the whole group or a data.table with sample sizes per gene (symbol, ss)")

  dataset = rbind(cases, controls, fill=T)


  if(nrow(dataset) == 0){
    stop("dataset has zero rows.")
  }


  newdata = remove_uncovered_residues(dataset, case_coverage, control_coverage, cov_threshold, messages)
  dataset = newdata[[1]]
  case_coverage = newdata[[2]]
  control_coverage = newdata[[3]]
  blacklist = newdata[[4]]


  # count variants per genes
  gene_counts = dataset[rep(1:.N, ceiling(ac^0.8)),.(ncases=sum(aff==1), ncontrols=sum(aff==0)), by=.(symbol)]
  gene_counts[ncases==0, ncases:=0.5]
  gene_counts[ncontrols==0, ncontrols:=0.5]


  # adjust case sample sizes
  if(null1){

    if(case_auto){

      mean_cov = case_coverage[,.(mean_cov=mean(over_10)), by=symbol]
      ss = attributes(cases)$ss
      mean_cov_ss = merge(mean_cov, ss, by="symbol")
      mean_cov_ss[,case_ss:=round(mean_cov*ss)]

      gene_counts = merge(gene_counts, mean_cov_ss[,.(symbol, case_ss)])

    } else {

      if(length(cases_ss) == 1){

        case_ss_by_gene = case_coverage[,.(case_ss=round(mean(over_10) * cases_ss)), by=symbol]
        gene_counts = merge(gene_counts, case_ss_by_gene)

      } else {

        mean_cov = case_coverage[,.(mean_cov=round(mean(over_10))), by=symbol]
        mean_cov = merge(mean_cov, cases_ss, by="symbol")
        adj_ss = mean_cov[,.(symbol, ss=mean_cov*ss)]

        gene_counts = merge(gene_counts, adj_ss[,.(symbol, case_ss=ss)])
        nna = gene_counts[,sum(is.na(case_ss))]
        if(nna>0) message(paste0("Removed ", nna, " genes with no sample size given in 'cases_ss'"))
        gene_counts = na.omit(gene_counts)

      }

    }


  } else {

    if(length(cases_ss) == 1){

      gene_counts[,case_ss:=cases_ss]

    } else {

      gene_counts = merge(gene_counts, cases_ss[,.(symbol, case_ss=ss)])
      nna = gene_counts[,sum(is.na(case_ss))]
      if(nna>0) message(paste0("Removed ", nna, " genes with no sample size given in 'cases_ss'"))
      gene_counts = na.omit(gene_counts)

    }

  }


  # adjust control sample sizes
  if(null2){

    if(control_auto){

      mean_cov = control_coverage[,.(mean_cov=round(mean(over_10))), by=symbol]
      ss = attributes(controls)$ss
      mean_cov_ss = merge(mean_cov, ss, by="symbol")
      mean_cov_ss[,control_ss:=mean_cov*ss]

      gene_counts = merge(gene_counts, mean_cov_ss[,.(symbol, control_ss)])

    } else {

      if(length(controls_ss) == 1){

        control_ss_by_gene = control_coverage[,.(control_ss=round(mean(over_10) * controls_ss)), by=symbol]
        gene_counts = merge(gene_counts, control_ss_by_gene)

      } else {

        mean_cov = control_coverage[,.(mean_cov=round(mean(over_10))), by=symbol]
        mean_cov = merge(mean_cov, controls_ss, by="symbol")
        adj_ss = mean_cov[,.(symbol, ss=mean_cov*ss)]

        gene_counts = merge(gene_counts, adj_ss[,.(symbol, control_ss=ss)])
        nna = gene_counts[,sum(is.na(control_ss))]
        if(nna>0) message(paste0("Removed ", nna, " genes with no sample size given in 'controls_ss'"))
        gene_counts = na.omit(gene_counts)

      }

    }

  } else {

    if(length(controls_ss) == 1){

      gene_counts[,control_ss:=controls_ss]

    } else {

      gene_counts = merge(gene_counts, controls_ss[,.(symbol, control_ss=ss)])
      nna = gene_counts[,sum(is.na(control_ss))]
      if(nna>0) message(paste0("Removed ", nna, " genes with no sample size given in 'cases_ss'"))
      gene_counts = na.omit(gene_counts)

    }
  }


  # calculate p-values
  gene_counts[,burden_pvalue:=mapply(function(n1, n2, ss1, ss2) burden_test(n1=n1, n2=n2, ss1=ss1, ss2=ss2, alternative=alternative),
                                     ncases, ncontrols, case_ss, control_ss)]


  if(covstats){

    gene_counts = merge(gene_counts, blacklist, all.x=T)

    gene_counts[,covstats:=paste0("adjSS_cases=", case_ss,
                                  ",adjSS_controls=", control_ss,
                                  ",rem=", ranges)]

    return(gene_counts[,!"ranges", with=F])

  } else {

    return(gene_counts)

  }

}

