#' Annotate a dataset with GnomAD frequencies
#' @author Adam Waring - adam.waring@@msdtc.ox.ac.uk
#' @return Returns a data.table
#'
#' @param dataset dataset to annotate with GnomAD freqs; must have either chr, pos, ref, alt or txconsq
#'
#' @keywords RVAT, distribution, cluster, genetics, case-control, gene
#' @export

annotate_with_gnomad_freqs = function(dataset){

  if(all(c("chr", "pos", "ref", "alt") %in% colnames(dataset))){
    method = "chrpos"
  } else if("txconsq" %in% colnames(dataset)){
    method = "txconsq"
  } else {
    stop("Must provide columns (chr, pos, ref, alt) or columns (txconsq) i.e. HGVSc transcript consequence")
  }


  data("econtrols")
  data("econtrols_nonpass")
  data("gcontrols")
  data("gcontrols_nonpass")

  econtrols = rbind(econtrols, econtrols_nonpass, fill=T)
  gcontrols = rbind(gcontrols, gcontrols_nonpass, fill=T)


  if(method == "chrpos"){

    data2 = merge(dataset, econtrols[,.(chr, pos, ref, alt, global, popmax, strict)], by = c("chr", "pos", "ref", "alt"), all.x=T)
    data3 = merge(data2, gcontrols[,.(chr, pos, ref, alt, af_g2)], by = c("chr", "pos", "ref", "alt"), all.x=T)

  } else if(method == "txconsq"){

    data2 = merge(dataset, econtrols[,.(txconsq, global, popmax, strict)], by = c("txconsq"), all.x=T)
    data3 = merge(data2, gcontrols[,.(txconsq, af_g2)], by = c("txconsq"), all.x=T)

  }

  freqs = c("af_g2", "global", "popmax", "strict")
  data3[,(freqs):=lapply(.SD, function(x) ifelse(is.na(x), 0, x)), .SDcols = freqs]

  return(data3)

}
