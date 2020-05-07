# ClusterBurden

To install in R: 

> if(!"devtools" %in% installed.packages()[,"Package"]){\
>   install.packages("devtools")\
> }\
> library(devtools)\
> install_github("adamwaring/ClusterBurden")

Due to some large variant and coverage files, download may be slow. if interest in specific functions, feel free to copy directly from the R/ folder in this gitHub repo. 

## Rare variant association using both positional (variant amino acid residue clustering) and burden signals - with control for uneven coverage.

### Relevance of method:
This method is designed for case-control studies where the condition under investigation is caused by rare missense variants in protein-coding regions. This covers many Mendelian/rare diseases. When a simple burden test does not yield significant associations, consideration of clustering information can increase power to detect low penetrance undiscovered disease-genes. 

### In this package 

This package allows exome-wide scans searching for both burden and clustering signals as well as their combined significance. GnomAD controls (exome and genomes v2) are available to automatically generate a control set and their coverage files to adjust for uneven coverage. Once downloaded (some moderately large files present) exome-wide scans are quick. See the associated vignettes for more details about the methods, how to conduct exome-wide scans using automatic controls and formatting user data. 

### Functions overview

#### Statistical methods (all with coverage control)
* BIN-test: detect positional differences of missense variants between cases and controls
* burden_test: detect burden differences of missense variants between cases and controls 
* combine_ps: combine p-values from independent statistical tests
* BIN_test_WES: whole-exome scans for clustering signals
* burden_test_WES: whole-exome scans for burden signals 
* ClusterBurden_WES: whole-exome scans for burden, clustering and combined signals

#### Visualisation
* plot_distribs: plot variant and coverage distributions for a single gene
* plot_signif_distribs: plot variant and coverage distributions for top significant genes in WES analysis 
* manhhattan: manhattan plot for p-values from WES analysis 
* plot_features: plot variant positions alongside sequence features from uniprot 
* plot_residuals: plot standardized residuals from a chi-squared two-sample test

#### Formatting inputs
* collect_gnomad_controls: generate automatic controls based on GnomAD exomes or genomes v2
* remove_uncovered_residues: exclude residue positions from analysis which have poor coverage in cases or controls 
* get_residue_position: extract variant residue position from HGVSp protein consequence format
* annotate_with_gnomad_freqs: annotate a dataset with frequencies from GnomAD exomes and genomes for filtering 
* format_coverage: map a chr-pos coverage file into a gene-protein_index coverage file (canonical transcripts only)


