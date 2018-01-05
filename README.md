# ClusterScan
ClusterScan is a tool which search for genomic clusters starting from a feature annotation. It allow the user to scan an annotation file (BED format) and get clusters coordinates in output. ClusterScan also need an additional two-columns file storing the feature names and the corresponding categorical information (including Gene Ontology, KEGG, Pfam accessions, etc). The user can also define a custom set of categories, making ClusterScan very flexible.

## How the tool works:
ClusterScan is composed by two different algorithms that perform the search but it also offer filters to select the minimum number of features to validate a cluster (that can’t be a number lower than 2). The two algorithms are:

- **clusterdist**: scans the features using bedtools merge in order to find those features which are separated by a maximum distance in base pairs that can be selected by the user (p: --distance). Some studies based on gene families, for example in human and mouse, have estimated genes within 500 kb to be in cluster ([Niimura et al. 2003](https://www.ncbi.nlm.nih.gov/pubmed/14507991); [Tadepally et al. 2008](https://www.ncbi.nlm.nih.gov/pubmed/18559114)).

- **clustermean**: divides the genome in sliding windows and calculates the mean number of features and the standard deviation for each accession. After that, clustermean searches for those windows containing a number of features higher or equal to the relation _mean+n*stdv_ in which *n* can be set by the user (p: --seed). Thus, the seed parameter set the number of standard deviations to identify a window which serves as the beginning of the cluster. After this step, the algorithm similarly tries to extend the cluster in both the directions starting from the seed using the same relation _mean+n*stdv_ in which *n* can be set again by the user (p: --extension). Thus, the extension parameter set the number of standard deviations to identify the window(s) which serve to extend the cluster. Lastly, clustermean trims the clusters in order to replace the cluster start/end represented by the first window start and the last window end respectively, with the first gene in cluster start and the last gene in cluster end.

## Dependencies:
ClusterScan requires [Python](https://www.python.org/downloads/release/python-2714/) (v2.7.x). [Bedtools](https://github.com/arq5x/bedtools2) (v2.25.0+) and [R](https://www.r-project.org/) (v3.0.0+) are needed to be in the user path. It also need some other Python libraries which can be easily installed via [pip](https://pip.pypa.io/en/stable/installing/) from [PyPI](https://pypi.python.org/pypi):

- [pandas](https://pandas.pydata.org/) (v0.19.1)
- [pybedtools](https://daler.github.io/pybedtools/) (v0.7.8)
- [rpy2](https://rpy2.readthedocs.io/en/version_2.8.x/) (v2.8.3)

Finally, in order to draw high quality clusters distributions for features in the top 10 clusters found (by number of features), it is also required to install the R library [ggplot2](http://ggplot2.org/) (v2.0.0+). ClusterScan is tested on [Ubuntu](https://www.ubuntu.com/) (v14.04LTS+)

The program is distributed with a file named _setup.py_ that can be found in the main directory. This can be used to automatically install the python dependencies typyng:
```
python setup.py install
```

## Options:
ClusterScan provides different parameters in order to finely tune the cluster search. Some of them are algorithm-specific whereas other are in common between clusterdist and clustermean:

```
clusterdist:
  -d, --dist=<bp>	Maximum distance between features in bp [default: 500000].

clustermean:
  -w, --window=<bp>	Window size [default: 500000].
  -s, --slide=<bp>	Sliding size [default: 250000].
  -k, --seed=<n>	Number of standard deviations to identify a window which serves as the beginning of the cluster [default: 3].
  -e, --extension=<n>	Number of standard deviations to identify the window(s) which serve to extend the cluster [default: 2].

shared options:
  -o, --output PATH	Specify output path [default: ./].
  -a, --analysis NAME	Specify optional analysis name for output files.
  -n, --nf=<n>		Minimum number of features per cluster [default: 2].
  -c, --category LIST   Comma separated list of one or more specific categories to be analyzed [e.g. PF00001,PF00002].
  --info FILE           Specify optional file to describe categories.
  --singletons          Identify singletons after clusters and bystanders annotation.
```

An example of execution with clusterdist:
```
python clusterscan.py clusterdist my_genes.bed my_categories.txt -d 250000 -a analysis_01
```

An example of execution with clustermean:
```
python clusterscan.py clustermean my_genes.bed my_categories.txt --info my_descriptions.txt -w 10000 -s 5000
```

## Input/Output:
ClusterScan is composed by three files: _clusterscan.py_, _algos.py_ and _docopt.py_. The first represent the tool itself, the second stores the functions belonging to the two algorithms described above, and the latter is a dependence that provides an interface for the command-line menu. The three files need to be in the same folder.

ClusterScan requires two mandatory files and an optional file in input. The first mandatory file is represented by a six-field bed file of features annotation. It may follow the classical organization of the [BED file format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1). The second file is a two-column table containing the name of the features as first column and a string representing a categorical information used in the cluster search. This "category" can be represented for esample by Gene Ontology, KEGG, Pfam accessions etc. The user can also search for clusters based on a custom made series of categories describing them using a third two-column tab-delimited .txt file (category, description) that can be read by the program (p: --info). Moreover, the user can pass to the program a comma separated list of categories on which perform the analysis (p: --category). It is not uncommon that a feature is associated with different accessions at the same time and, in this case, the feature can take part of multiple different kind of clusters. The user can also choose to identify features which exist as _singletons_ (p: --singletons). Those represent features belonging to a certain category, for which at least one cluster was annotated, that are outside of any of the corresponding clusters.

ClusterScan gives six or seven files in output depending on whether the user has chosen or not to identify the singletons:

|File|Description|
|---|---|
|	clusters.tsv | stores the coordinates of all the clusters found. It contains an unique ID for each cluster; the category to which the cluster belongs; the chromosome/scaffold (chr) on which it resides; its start and end coordinates; the number of features within the cluster (n_features) and the number of features which overlap the cluster but belong to a different kind of category (n_bystanders). |
|		clusters.bed | a bed version of the previous file for bedtools/bedops compatibility. Fields respect the format and  are (from first to last column): chromosome/scaffold on which the cluster resides; its start/end coordinates using a 0-based start and 1-based end system of coordinates; the cluster ID; the number of features within the cluster; the strand; the category describing the cluster. |
|		summary.tsv | contains a category based summary of the cluster analysis with the category; the number of clusters found for each category (n_clusters); the total number of features found for clusters belonging to each category (n_ft); the number of bystanders (n_bs); the number of features and bystanders found in the clusters which contains the minimal and the maximal number of these features (max_ft and min_ft / max_bs and min_bs) for each category. |
|		features.tsv | a list of features found to be in overlap with clusters. It contains the chromosome/scaffold (chr) on which the feature resides; its start and end coordinates; the feature name; the strand; the ID of the cluster to which it overlaps; the category to which the cluster (and the feature) belongs. |
|		bystanders.tsv | a list of bystanders found to be in overlap with clusters. It contains the chromosome/scaffold (chr) on which the feature resides; its start and end coordinates; the feature name; the strand; the ID of the cluster to which it overlaps; the category to which the cluster (but no the feature) belongs. |
|		distribution.pdf | is an histogram which shows the distribution of features in the top-10 categories by number of features. |
|		singletons.tsv | a list of singletons found to be out from the annotated clusters for each category. It contains the chromosome/scaffold (chr) on which it resides; its start and end coordinates; the feature name; the strand; the category to which the feature belongs. |

The program assigns an unique ID for each cluster during the analysis, joining the “C” letter with a progressive number starting from 1. The user can trace back a specific cluster at any time in the different output files simply looking at the ID reported in the corresponding column.
