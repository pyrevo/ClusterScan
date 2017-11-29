# ClusterScan

ClusterScan is a tool which search for genomic clusters starting from a feature annotation. It allow the user to scan an annotation file (BED format) and get clusters coordinates in output. ClusterScan also need an additional two-columns file storing the feature names and the corresponding accession (including Gene Ontology, KEGG, Pfam accessions, etc). The user can also use a custom set of accessions, making ClusterScan very flexible.

## How the tool works:

ClusterScan is composed by two different algorithms that perform the search but it also offer filters to select the minimum number of features to validate a cluster (that can’t be a number lower than 2). The two algorithms are:

- **clusterdist**: scans the features using bedtools merge in order to find those features which are separated by a maximum distance in base pairs that can be selected by the user (p: --distance). Some studies based on gene families, for example in human and mouse, have estimated genes within 500 kb to be in cluster (Niimura et al. 2003; Tadepally et al. 2008).

- **clustermean**: divides the genome in sliding windows and calculates the mean number of features and the standard deviation for each accession. After that, clustermean searches for those windows containing a number of features higher or equal to the relation _mean+n*stdv_ in which *n* can be set by the user (p: --seed). Thus the seed parameter set the number of standard deviations to identify a window which serves as the beginning of the cluster. After this step, the algorithm similarly tries to extend the cluster in both the directions starting from the seed using the same relation _mean+n*stdv_ in which *n* can be set again by the user (p: --extension). Thus the extension parameter set the number of standard deviations to identify the window(s) which serve to extend the cluster. Lastly, clustermean trims the clusters in order to replace the cluster start/end represented by the first window start and the last window end respectively with the first gene in cluster start and the last gene in cluster end.

# Dependencies:
ClusterScan requires [Python](https://www.python.org/downloads/release/python-2714/) (v.2.7.x). [Bedtools](https://github.com/arq5x/bedtools2) (v.2.25.0+) and [R](https://www.r-project.org/) (v.3.0.0+) are needed to be in the user path. It also need some other Python libraries which can be easily installed via [pip](https://pip.pypa.io/en/stable/installing/) from [PyPI](https://pypi.python.org/pypi):

- [pybedtools](https://daler.github.io/pybedtools/)
- [pandas](https://pandas.pydata.org/)
- [rpy2](https://rpy2.readthedocs.io/en/version_2.8.x/)

Finally, in order to draw high quality clusters distributions for features in the top 10 clusters found (by number of features), it is also required to install the R library [ggplot2](http://ggplot2.org/) (v.2.0.0+).
