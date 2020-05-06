#library(ClusterBurden)
library(ggrepel)
sapply(list.files("R/", pattern = "*.r$", full.names = T), source)

n_random_genes = 200
disease = "cardiomyopathy"

genes = find_genes(n_random_genes, disease)

controls = collect_gnomad_controls(genenames=genes)
cases = collect_example_cases(controls, disease)

binpval = BIN_test_WES(cases, controls, covstats = T)
# slightly inflated lambda(binpval[!is.na(binp), binp]) = 1.07


burdenpval = burden_test_WES(cases, controls, covstats=T)


# 2. add warnings
# 3. add more check of inputs going into these functions#
# 4. improve plot functions and add data files to data in package after reducing

pvals = merge(binpval, burdenpval, by="symbol", all = T)
pvals[,ClusterBurden:=mapply(function(x, y) combine_ps(x, y), BIN.test_pvalue, burden_pvalue)]

# manhattan plot
manhattan(pvals[cov_flag2 == "pass"], "ClusterBurden", 20)


plot_signif_distribs(pvals, cases, controls, n_genes = 8)


plot_features("STON2", pvals, cases, controls)

lambda = function(p) median(qchisq(p, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
lambda(binpval[!is.na(binp) & cov_flag2 != "vhigh", binp])
lambda(binpval[!is.na(binp) & cov_flag2 == "pass" & cov_flag1 == "pass" & pl_flag != "vhigh", binp])

# test out lambda after ceiling ^ 0.8 to allele counts

ceiling((1:19)^(0.8))
