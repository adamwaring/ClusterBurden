# library(data.table)
#
#
# gnomad exomes v2
e2_files = list.files("C:/Users/adamwar/Google Drive/dphil/RAW_DATA/gnomad_e2/", pattern="*.txt$", full.names = T)
e2 = do.call(rbind, lapply(e2_files, fread))
e2$dataset = "gnomad_e2"


# gnomad genomes v2
g2_files = list.files("C:/Users/adamwar/Google Drive/dphil/RAW_DATA/gnomad_g2/", pattern="*.txt$", full.names = T)[-24] # chromosome Y size = 0
g2 = do.call(rbind, lapply(g2_files, fread))
g2$dataset = "gnomad_g2"


# clinvar
cln_files = list.files("C:/Users/adamwar/Google Drive/dphil/RAW_DATA/clinvar_050220/", pattern="*.txt$", full.names = T)[-24]
cln = do.call(rbind, lapply(cln_files, fread))
cln$dataset = "clinvar"


colnames(e2)[7:36] = paste0(colnames(e2)[7:36], "_e2")
colnames(g2)[7:36] = paste0(colnames(g2)[7:36], "_g2")


cln = merge(cln, e2[,c(1:4, 7:36)], by = c("CHROM", "POS", "REF", "ALT"), all.x=TRUE)
cln = merge(cln, g2[,c(1:4, 7:36)], by = c("CHROM", "POS", "REF", "ALT"), all.x=TRUE)


e2 = merge(e2, cln[,c(1:4, 7:9)], by = c("CHROM", "POS", "REF", "ALT"), all.x=TRUE)
g2 = merge(g2, cln[,c(1:4, 7:9)], by = c("CHROM", "POS", "REF", "ALT"), all.x=TRUE)


e2 = merge(e2, g2[,c(1:4, 7:36)], by = c("CHROM", "POS", "REF", "ALT"), all.x=TRUE)
g2 = merge(g2, e2[,c(1:4, 7:36)], by = c("CHROM", "POS", "REF", "ALT"), all.x=TRUE)



FEATURES = colnames(e2)[-c(52:80)]
X = rbind(e2[,FEATURES, with=F], g2[,FEATURES, with=F], cln[,FEATURES, with=F])


# for multi-position variants just take the first position affected
X[, Protein_position := gsub("-.*", "", Protein_position)]

# format txconsq and pconsq
X[, c("txconsq", "pconsq"):= .(gsub(".*:", "", HGVSc), gsub(".*:", "", HGVSp))]

# format clnsig
clnsig = c("Benign", "Likely_benign", "Uncertain_significance", "Likely_pathogenic", "Pathogenic")
X[grepl("[Cc]onflict", CLNSIG), CLNSIG:="Uncertain_significance"]

format_cln = function(x){
  if(is.na(x)) return(NA)

  clnsigs = strsplit(gsub("/", ",", x), ",")[[1]]
  scores = sapply(clnsigs, function(x) ifelse(x %in% clnsig, which(clnsig == x), NA))

  agg_score = function(x) round(mean(c(rep(x, 4), 3), na.rm=T),0)
  score = agg_score(scores)
  ifelse(score %in% 1:5, clnsig[score], NA)
}

X[, CLNSIG2 := sapply(CLNSIG, format_cln)]

# gnomAD freqs
freq_cols = colnames(X)[grepl("^A[CF].*[eg]2", colnames(X))]
X[, (freq_cols) := lapply(.SD, function(x) ifelse(is.na(x), 0, x)), .SDcols = freq_cols]

# change column names
colnames(X) = tolower(colnames(X))
colnames(X) = gsub("\\+", "", gsub("-", "_", gsub("_score", "", colnames(X))))

# numeric columns
C=colnames(X)
num_cols = c(C[grepl("position", C)], C[grepl(".*_[eg]2", C)])
X[,(num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]

# factor columns
X[,clnsig2 := factor(clnsig2, levels = c("Benign", "Likely_benign", "Uncertain_significance", "Likely_pathogenic", "Pathogenic"))]

# strict popmax
freq_names = c(C[grepl("^af.*e2", C)], "af_g2")
X$strict_popmax = X[,max(.SD, na.rm = T), by=seq_len(nrow(X)), .SDcols = freq_names][[2]]



pops = c("afr","amr","asj","eas","fin","nfe","sas","oth")

poisCI = function(x, y, alpha){
  ifelse(is.na(x) | x==0, 0, qgamma(1-alpha, x)/y)
}


for(pop in pops){
  print(pop)
  X[,paste0("pci95_", pop):=.(poisCI(X[[paste0("ac_", pop, "_e2")]], X[[paste0("an_", pop, "_e2")]], 0.95))]
  X[,paste0("pci90_", pop):=.(poisCI(X[[paste0("ac_", pop, "_e2")]], X[[paste0("an_", pop, "_e2")]], 0.9))]
}

X[,"pci95_global":=.(poisCI(ac_e2, an_e2, 0.95))]
X[,"pci90_global":=.(poisCI(ac_e2, an_e2, 0.9))]

X[,"pci95_g2":=.(poisCI(ac_g2, an_g2, 0.95))]
X[,"pci90_g2":=.(poisCI(ac_g2, an_g2, 0.9))]

X[,pci95_popmax:=apply(.SD, 1, max), .SDcols=paste0("pci95_", pops[-c(3, 5, 8)])]
X[,pci95_strict:=apply(.SD, 1, max), .SDcols=paste0("pci95_", pops)]
X[,pci90_strict:=apply(.SD, 1, max), .SDcols=paste0("pci90_", pops)]
X[,pci90_popmax:=apply(.SD, 1, max), .SDcols=paste0("pci90_", pops[-c(3, 5, 8)])]




# add case-control column
X[,cohort:=ifelse(grepl("gnomad", dataset), 0, 1)]


canons = X[,names(which.max(table(feature))), by=symbol]
X = X[feature%in%canons$V1]

missense = c("missense_variant", "inframe_insertion", "inframe_deletion", "protein_altering_variant")

X = X[consequence %in% missense]

large_inframes = X[,which(consequence%in%missense[2:3] & abs(nchar(ref)-nchar(alt)>9))]
X = X[-large_inframes]

X[,n_res:=mapply(function(x, y) max(nchar(x)%/%3, nchar(y)%/%3), ref, alt)]


econtrols = X[dataset=="gnomad_e2"]
gcontrols = X[dataset=="gnomad_g2"]

econtrols_nonpass = econtrols[filter!="PASS"]
econtrols_nonpass = econtrols_nonpass[,.(chr=chrom, pos, ref, alt, symbol, txconsq, pci95_global, pci90_global, pci95_popmax, pci95_strict, pci90_popmax, pci90_strict)]
freq_cols = c("pci95_global", "pci90_global", "pci95_popmax", "pci95_strict", "pci90_popmax", "pci90_strict")
econtrols_nonpass[,(freq_cols):=lapply(.SD, signif, 3), .SDcols=freq_cols]

econtrols = econtrols[filter=="PASS", .(chr=chrom, pos, ref, alt, ac_e2, txconsq, pconsq, protein_position, symbol, consequence, pci95_g2, pci90_g2, pci95_global, pci90_global, pci95_popmax, pci95_strict, pci90_popmax, pci90_strict)]
freq_cols = c("pci95_g2", "pci90_g2", "pci95_global", "pci90_global", "pci95_popmax", "pci95_strict", "pci90_popmax", "pci90_strict")
econtrols[,(freq_cols):=lapply(.SD, signif, 3), .SDcols=freq_cols]

save(econtrols, file="data/econtrols.RData", compress = "bzip2")
save(econtrols_nonpass, file="data/econtrols_nonpass.RData", compress = "bzip2")

gcontrols_nonpass = gcontrols[filter!="PASS"]
gcontrols_nonpass = gcontrols_nonpass[,.(chr=chrom, pos, ref, alt, symbol, txconsq, pci95_g2, pci90_g2)]
freq_cols = c("pci95_g2", "pci90_g2")
gcontrols_nonpass[,(freq_cols):=lapply(.SD, signif, 3), .SDcols=freq_cols]


gcontrols = gcontrols[filter=="PASS", .(chr=chrom, pos, ref, alt, ac_g2, txconsq, pconsq, protein_position, symbol, consequence, pci95_g2, pci90_g2, pci95_global, pci90_global, pci95_popmax, pci95_strict, pci90_popmax, pci90_strict)]
freq_cols = c("pci95_g2", "pci90_g2", "pci95_global", "pci90_global", "pci95_popmax", "pci95_strict", "pci90_popmax", "pci90_strict")
gcontrols[,(freq_cols):=lapply(.SD, signif, 3), .SDcols=freq_cols]


save(gcontrols, file="data/gcontrols.RData", compress = "bzip2")
save(gcontrols_nonpass, file="data/gcontrols_nonpass.RData", compress = "bzip2")

# change functions to not include filter=="PASS"
# add documentation for new data
# change anno gnomad function to load both files

load("data/exome_cov.RData")
exome_cov[,over_10:=signif(over_10, 3)]
save(exome_cov, file="data/exome_cov.Rdata", compress = "bzip2")

load("data/genome_cov.RData")
genome_cov[,over_10:=signif(over_10, 3)]
save(genome_cov, file="data/genome_cov.Rdata", compress = "bzip2")
