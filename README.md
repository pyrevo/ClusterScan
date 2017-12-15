# ClusterScan
ClusterScan is a tool which search for genomic clusters starting from a feature annotation. It allow the user to scan an annotation file (BED format) and get clusters coordinates in output. ClusterScan also need an additional two-columns file storing the feature names and the corresponding accession (including Gene Ontology, KEGG, Pfam accessions, etc). The user can also define a custom set of accessions, making ClusterScan very flexible.

## How the tool works:
ClusterScan is composed by two different algorithms that perform the search but it also offer filters to select the minimum number of features to validate a cluster (that can’t be a number lower than 2). The two algorithms are:

- **clusterdist**: scans the features using bedtools merge in order to find those features which are separated by a maximum distance in base pairs that can be selected by the user (p: --distance). Some studies based on gene families, for example in human and mouse, have estimated genes within 500 kb to be in cluster ([Niimura et al. 2003](https://www.ncbi.nlm.nih.gov/pubmed/14507991); [Tadepally et al. 2008](https://www.ncbi.nlm.nih.gov/pubmed/18559114)).

- **clustermean**: divides the genome in sliding windows and calculates the mean number of features and the standard deviation for each accession. After that, clustermean searches for those windows containing a number of features higher or equal to the relation _mean+n*stdv_ in which *n* can be set by the user (p: --seed). Thus the seed parameter set the number of standard deviations to identify a window which serves as the beginning of the cluster. After this step, the algorithm similarly tries to extend the cluster in both the directions starting from the seed using the same relation _mean+n*stdv_ in which *n* can be set again by the user (p: --extension). Thus the extension parameter set the number of standard deviations to identify the window(s) which serve to extend the cluster. Lastly, clustermean trims the clusters in order to replace the cluster start/end represented by the first window start and the last window end respectively with the first gene in cluster start and the last gene in cluster end.

## Dependencies:
ClusterScan requires [Python](https://www.python.org/downloads/release/python-2714/) (v2.7.x). [Bedtools](https://github.com/arq5x/bedtools2) (v2.25.0+) and [R](https://www.r-project.org/) (v3.0.0+) are needed to be in the user path. It also need some other Python libraries which can be easily installed via [pip](https://pip.pypa.io/en/stable/installing/) from [PyPI](https://pypi.python.org/pypi):

- [pandas](https://pandas.pydata.org/) (v0.19.1)
- [pybedtools](https://daler.github.io/pybedtools/) (v0.7.8)
- [rpy2](https://rpy2.readthedocs.io/en/version_2.8.x/) (v2.8.3)

Finally, in order to draw high quality clusters distributions for features in the top 10 clusters found (by number of features), it is also required to install the R library [ggplot2](http://ggplot2.org/) (v2.0.0+). ClusterScan is tested on [Ubuntu](https://www.ubuntu.com/) (v14.04LTS+)

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
  --info FILE		Specify optional file to describe accessions.
```

An example of execution:
```
python clusterscan.py clusterdist my_genes.bed my_accessions.txt -d 250000 -a analysis_01
```

## Input/Output:
ClusterScan is composed by two files: _clusterscan.py_ and _docopt.py_. The first represent the tool itself while the latter is a  dependence that provides an interface for the command-line menu. The two files need to be in the same folder.

ClusterScan requires two mandatory files and an optional file in input. The first mandatory file is represented by a six-field bed file of features annotation. It may follow the classical organization of the [BED file format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1). The second file is a two-column table containing the name of the features on the first column and an accession for which the cluster search need to be called. This accession can belong to Gene Ontology, KEGG, Pfam, etc. The user can also search for clusters based on a custom made series of accessions. It is not uncommon that a feature is associated with different accessions at the same time and in this case the feature can take part of multiple different kind of clusters. 

ClusterScan gives six files in output:

|File|Description|
|---|---|
|	clusters.csv | stores the coordinates of all the clusters found. It contains an ID for each cluster; the accession (ACC) for which the cluster was composed; the chromosome/scaffold (chr) on which it resides; its start and end coordinates; the number of features within the cluster (n_features) and the number of features which overlap the cluster but belong to a different kind of accession (n_bystanders). |
|		clusters.bed | a bed version of the previous file for bedtools/bedops compatibility. Fields respect the format and  are (from first to last column): chromosome/scaffold on which the cluster resides; its start/end coordinates using a 0-based start and 1-based end system of coordinates; the cluster ID; the number of features within the cluster; the strand; the cluster accession. |
|		summary.csv | contains a per-accession summary of the cluster analysis with the accession (ACC); the number of clusters found for that accession (n_clusters) the total number of features found for clusters belonging to that accession (n_ft); the number of bystanders (n_bs); the number of features and bystanders found in the clusters which contains the minimal and the maximal number of these features (max_ft and min_ft / max_bs and min_bs). |
|		features.csv | is a list of features found to be in overlap with clusters. It contains the chromosome/scaffold (chr) on which it resides; its start and end coordinates; the feature name; the score; the strand; an ID of the cluster on which it overlaps (cluster_ID); the accession of the cluster in which it overlaps (cluster_ACC). |
|		bystanders.csv | is a list of bystanders found to be in overlap with clusters. It contains exactly the same fields described for the file features.csv. |
|		distribution.pdf | is an histogram which shows the distribution of features in the per-accession top-10 clusters by number of features. |

## Disambiguating the meaning of ID and ACC acronyms:
In ClusterScan the word _ID_ refers to the cluster identifier. An ID for each cluster is assigned during the analysis joining the “C” letter with a progressive number starting from 1. Contrariwise the word _ACC_ refers to the database accessions used to classify the features. They can come from Gene Ontology, KEGG, Pfam, etc. but the user can also provide its own custom set of accessions and describe them using a third two-column tab-delimited .txt file (ACC, description) that can be read by the program (p: --info).
