# library(data.table)
#
#
# # gnomad exomes v2
# e2_files = list.files("C:/Users/adamwar/Google Drive/dphil/RAW_DATA/gnomad_e2/", pattern="*.txt$", full.names = T)
# e2 = do.call(rbind, lapply(e2_files, fread))
# e2$dataset = "gnomad_e2"
#
#
# # gnomad genomes v2
# g2_files = list.files("C:/Users/adamwar/Google Drive/dphil/RAW_DATA/gnomad_g2/", pattern="*.txt$", full.names = T)[-24] # chromosome Y size = 0
# g2 = do.call(rbind, lapply(g2_files, fread))
# g2$dataset = "gnomad_g2"
#
#
# # clinvar
# cln_files = list.files("C:/Users/adamwar/Google Drive/dphil/RAW_DATA/clinvar_050220/", pattern="*.txt$", full.names = T)[-24]
# cln = do.call(rbind, lapply(cln_files, fread))
# cln$dataset = "clinvar"
#
#
# colnames(e2)[7:36] = paste0(colnames(e2)[7:36], "_e2")
# colnames(g2)[7:36] = paste0(colnames(g2)[7:36], "_g2")
#
#
# cln = merge(cln, e2[,c(1:4, 7:36)], by = c("CHROM", "POS", "REF", "ALT"), all.x=TRUE)
# cln = merge(cln, g2[,c(1:4, 7:36)], by = c("CHROM", "POS", "REF", "ALT"), all.x=TRUE)
#
#
# e2 = merge(e2, cln[,c(1:4, 7:9)], by = c("CHROM", "POS", "REF", "ALT"), all.x=TRUE)
# g2 = merge(g2, cln[,c(1:4, 7:9)], by = c("CHROM", "POS", "REF", "ALT"), all.x=TRUE)
#
#
# e2 = merge(e2, g2[,c(1:4, 7:36)], by = c("CHROM", "POS", "REF", "ALT"), all.x=TRUE)
# g2 = merge(g2, e2[,c(1:4, 7:36)], by = c("CHROM", "POS", "REF", "ALT"), all.x=TRUE)
#
#
#
# FEATURES = colnames(e2)[-c(52:80)]
# X = rbind(e2[,FEATURES, with=F], g2[,FEATURES, with=F], cln[,FEATURES, with=F])
#
#
# # for multi-position variants just take the first position affected
# X[, Protein_position := gsub("-.*", "", Protein_position)]
#
# # format txconsq and pconsq
# X[, c("txconsq", "pconsq"):= .(gsub(".*:", "", HGVSc), gsub(".*:", "", HGVSp))]
#
# # format clnsig
# clnsig = c("Benign", "Likely_benign", "Uncertain_significance", "Likely_pathogenic", "Pathogenic")
# X[grepl("[Cc]onflict", CLNSIG), CLNSIG:="Uncertain_significance"]
#
# format_cln = function(x){
#   if(is.na(x)) return(NA)
#
#   clnsigs = strsplit(gsub("/", ",", x), ",")[[1]]
#   scores = sapply(clnsigs, function(x) ifelse(x %in% clnsig, which(clnsig == x), NA))
#
#   agg_score = function(x) round(mean(c(rep(x, 4), 3), na.rm=T),0)
#   score = agg_score(scores)
#   ifelse(score %in% 1:5, clnsig[score], NA)
# }
#
# X[, CLNSIG2 := sapply(CLNSIG, format_cln)]
#
# # gnomAD freqs
# freq_cols = colnames(X)[grepl("^A[CF].*[eg]2", colnames(X))]
# X[, (freq_cols) := lapply(.SD, function(x) ifelse(is.na(x), 0, x)), .SDcols = freq_cols]
#
# # change column names
# colnames(X) = tolower(colnames(X))
# colnames(X) = gsub("\\+", "", gsub("-", "_", gsub("_score", "", colnames(X))))
#
# # numeric columns
# C=colnames(X)
# num_cols = c(C[grepl("position", C)], C[grepl(".*_[eg]2", C)])
# X[,(num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]
#
# # factor columns
# X[,clnsig2 := factor(clnsig2, levels = c("Benign", "Likely_benign", "Uncertain_significance", "Likely_pathogenic", "Pathogenic"))]
#
# # strict popmax
# freq_names = c(C[grepl("^af.*e2", C)], "af_g2")
# X$strict_popmax = X[,max(.SD, na.rm = T), by=seq_len(nrow(X)), .SDcols = freq_names][[2]]
#
# # add case-control column
# X[,cohort:=ifelse(grepl("gnomad", dataset), 0, 1)]
#
#
# canons = X[,names(which.max(table(feature))), by=symbol]
# X = X[feature%in%canons$V1]
#
# missense = c("missense_variant", "inframe_insertion", "inframe_deletion", "protein_altering_variant")
#
# X = X[consequence %in% missense]
#
# large_inframes = X[,which(consequence%in%missense[2:3] & abs(nchar(ref)-nchar(alt)>9))]
# X = X[-large_inframes]
#
# load("C:/Users/adamwar/Google Drive/dphil/RAW_DATA/gnomad_clinvar_missense_incl_passfilter.Rdata")
#
# X = X[,.(chrom, pos, ref, alt, filter, ac_e2, af_e2, af_popmax_e2, ac_g2, af_g2, txconsq, pconsq, protein_position, symbol, consequence, feature, clnsig, clnsig2, clndn, strict_popmax, dataset)]
# X[,n_res:=mapply(function(x, y) max(nchar(x)%/%3, nchar(y)%/%3), ref, alt)]
#
#
# econtrols = X[dataset=="gnomad_e2",.(filter, chr=chrom, pos, ref, alt, symbol, txconsq, pconsq, n_res, protein_position, ac=ac_e2, global=af_e2, popmax=af_popmax_e2, strict=strict_popmax)]
# save(econtrols, file="data/econtrols.Rdata", compress = "bzip2")
#
# gcontrols = X[dataset=="gnomad_g2",.(filter, chr=chrom, pos, ref, alt, symbol, txconsq, pconsq, n_res, protein_position, ac=ac_g2, af_g2, global=af_e2, popmax=af_popmax_e2, strict=strict_popmax)]
# save(gcontrols, file="data/gcontrols.Rdata", compress = "bzip2")

# load("data/econtrols.RData")
#
# econtrols_nonpass = econtrols[filter!="PASS"]
# econtrols_nonpass = econtrols_nonpass[,.(chr, pos, ref, alt, symbol, txconsq, global, popmax, strict)]
#
# econtrols = econtrols[filter=="PASS", .(chr, pos, ref, alt, symbol, txconsq, n_res, protein_position, ac, global, popmax, strict)]
#
# save(econtrols, file="data/econtrols.RData", compress = "bzip2")
# save(econtrols_nonpass, file="data/econtrols_nonpass.RData", compress = "bzip2")
#
#
# load("data/gcontrols.RData")
#
# gcontrols_nonpass = gcontrols[filter!="PASS"]
# gcontrols_nonpass = gcontrols_nonpass[,.(chr, pos, ref, alt, symbol, txconsq, af_g2, global, popmax, strict)]
#
# gcontrols = gcontrols[filter=="PASS", .(chr, pos, ref, alt, symbol, txconsq, n_res, protein_position, ac, af_g2, global, popmax, strict)]
#
# save(gcontrols, file="data/gcontrols.RData", compress = "bzip2")
# save(gcontrols_nonpass, file="data/gcontrols_nonpass.RData", compress = "bzip2")

# change functions to not include filter=="PASS"
# add documentation for new data
# change anno gnomad function to load both files
