#load("C:/Users/adamwar/Google Drive/dphil/RAW_DATA/gnomad_clinvar_missense_240320.Robj")


#
#
#
# x3 = X[dataset=="clinvar"]
#
# x3[,n_res:=mapply(function(x, y) max(nchar(x)%/%3, nchar(y)%/%3), ref, alt)]
#
# x3 = x3[n_res <= 5]
#
# x3 = x3[,.(symbol, n_res, protein_position, clnsig, clndn, global=af_e2, popmax=af_popmax_e2, strict=strict_popmax)]
#
# cases = x3
# save(cases, file="data/clinvar_cases.Rdata")

# library(adam)
# exome_cov1 = fread(pasteraw("coverage/gnomad.exomes.coverage.summary.reduced.txt"))
# exome_cov = format_coverage(exome_cov1[,.(chrom, pos, over_10)])
#
# genome_cov1 = rbindlist(lapply(list.files(pasteraw("coverage/"), full.names = T), fread))
# genome_cov = unique(genome_cov1)
#
# load("data/canon_txs.Rdata")
# genome_cov = genome_cov[tx%in%canon_txs$tx]
# #
# exome_cov = exome_cov[,.(symbol, protein_position=index, over_10)]
# genome_cov = genome_cov[,.(symbol, protein_position=index, over_10)]
#
# save(exome_cov, file="data/exome_cov.Rdata")
# save(genome_cov, file="data/genome_cov.Rdata")

# bioinfFolder = 'C:/Users/adamwar/Google Drive/dphil/RAW_DATA/bioinf_tables/'
# x = fread(file.path(bioinfFolder, "Homo_sapiens.GRCh37.87.gtf_AAseq.txt"))
#
# load("data/canon_txs.Rdata")
# x = x[V5%in%canon_txs$tx]
# gene_locations = setNames(x[,.(unique(V1), min(V2)), by=V14], c("symbol", "chr", "pos"))
#
# save(gene_locations, file="data/gene_locations.Rdata")

#
# load("C:/Users/adamwar/Google Drive/dphil/RAW_DATA/bioinf_tables/protein_feature_start_end.Rdata")
#
#
# long[,feature_name:=sub("\\d+$", "", feature_name)]
#
# long[feature_name=="none", feature_name:=""]
# long[feature=="Cross-link", feature_name:=""]
# long[grepl("Selectivity filter", feature_name), feature_name := "Selectivity filter"]
#
#
#
# library(wesanderson)
# color_vector = sample(unlist(wes_palettes))
#
# seq_features = long
#
# save(seq_features, file="data/seq_features.Rdata")
# save(color_vector, file="data/color_vector.Rdata")



# load("C:/Users/adamwar/Google Drive/dphil/RAW_DATA/bioinf_tables/protein_nresidues.Rdata")
#
# load("data/genome_cov.Rdata")
# load("data/exome_cov.Rdata")
#
# max_genomes = genome_cov[,.(nresidues2=max(protein_position)), by=symbol]
# max_exomes = exome_cov[,.(nresidues3=max(protein_position)), by=symbol]
#
# nresidues = nresidues[symbol!=""]
#
# x1 = merge(nresidues, max_genomes, all=T)
# x2 = merge(x1, max_exomes, all=T)
#
#
# x2[,nresidues4:=apply(x2[,2:4], 1, function(x) names(table(x))[which.max(table(x))])]
#
# nresidues = x2[,.(symbol, nresidues=as.numeric(nresidues4))]
# save(nresidues, file="data/nresidues.Rdata")
