# ClusterBurden

## Rare variant association using both positional (variant amino acid residue clustering) and burden signals - with control for uneven coverage.

### Relevance of method:
This method is designed for case-control studies where the condition under investigation is caused by rare missense variants in protein-coding regions. This covers many Mendelian/rare diseases. When a simple burden test does not yield significant associations, consideration of clustering information can increase power to detect low penetrance undiscovered disease-genes. 

### In this package 

This package allows exome-wide scans searching for both burden and clustering signals as well as their combined significance. GnomAD controls (exome and genomes v2) are available to automatically generate a control set and their coverage files to adjust for uneven coverage. Once downloaded (some moderately large files present) exome-wide scans are quick. See the associated vignettes for more details about the methods, how to conduct exome-wide scans using automatic controls and formatting user data. 




