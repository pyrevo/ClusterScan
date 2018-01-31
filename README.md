# ClusterScan
ClusterScan is a tool to search for genomic clusters starting from genomic feature locations and their annotations. It allow the user to scan an annotation file (BED format with locations of specific features such as gene, transcripts, regulatory regions or anything that can be mapped on a genome) and get clusters coordinates in output. In order to build annotated clusters, ClusterScan need an additional two-columns file, storing the feature ids and the corresponding categorical information (such as Gene Ontology classes, KEGG, Pfam accessions, etc). The user can also define a custom set of categories, which makes ClusterScan very flexible.

***********************

- [How the tool works](#how-the-tool-works)
- [Installation](#installation)
	- [Installing Bedtools](#installing-bedtools)
	- [Installing R and the ggplot2 library](#installing-r-and-the-ggplot2-library)
	- [Installing the required Python libraries](#installing-the-required-python-libraries)
	- [Installing ClusterScan](#installing-clusterscan)
- [Options](#options)
- [Input/Output](#inputoutput)
- [Other resources](#other-resources)

***********************

## How the tool works:
ClusterScan can use two different algorithms that perform the search and it also offer filters to select the minimum number of features to call a cluster (that can’t be a number lower than 2). The two algorithms are:

- **clusterdist**: scans the features using bedtools merge in order to find those features, for each category, which are separated by a maximum distance in base pairs that can be selected by the user (p: --distance). Some studies based on gene families, for example in human and mouse, have estimated genes within 500 Kb to be in cluster ([Niimura et al. 2003](https://www.ncbi.nlm.nih.gov/pubmed/14507991); [Tadepally et al. 2008](https://www.ncbi.nlm.nih.gov/pubmed/18559114)).

- **clustermean**: divides the genome in sliding windows and calculates, for each category, the mean number of features and the standard deviation for each accession. After that, clustermean searches for those windows showing a Z-score bigger than a given value, i.e. containing a number of features higher or equal to the relation _mean+n*stdv_ in which *n* can be set by the user (p: --seed). Thus, the seed parameter set the number of standard deviations to identify a window which serves as the beginning of the cluster. After this step, the algorithm similarly tries to extend the cluster in both the directions starting from the seed using the same relation _mean+n*stdv_ in which *n* can be set again by the user (p: --extension). Thus, the extension parameter set the number of standard deviations to identify the window(s) which serve to extend the cluster. Lastly, clustermean trims the clusters in order to replace the cluster start/end represented by the most 5' feature start and the most 3' feature end respectively.

## Installation:
ClusterScan requires [Python](https://www.python.org/downloads/release/python-2714/) (v2.7.x). [Bedtools](https://github.com/arq5x/bedtools2) (v2.25.0+) and [R](https://www.r-project.org/) (v3.0.0+) are needed to be in the user path. In order to draw high quality clusters distributions for features in the top 10 clusters (by number of features) found, it is also required to install the R library [ggplot2](http://ggplot2.org/) (v2.0.0+).

### Installing Bedtools
To install bedtools you need to download the package via GitHub and compile from source:
```
wget https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz
tar -zxvf bedtools-2.26.0.tar.gz
cd bedtools2
make
```
Once you have finished you can add bedtools to your path opening a new terminal and typing:
```
sudo nano ~/.bashrc
```
scroll until you reach the bottom of the page and add the following string:
```
export PATH=$PATH:/<path-to-bedtools>/bedtools2/bin
```
Save and exit. To test your bedtools installation you must open a new terminal and type:
```
bedtools -h
```
More informations about the bedtools installation can be found [here](http://bedtools.readthedocs.io/en/latest/content/installation.html).

### Installing R and the ggplot2 library
To install R simply type:
```
sudo apt install r-base r-base-dev
```
and enter the program with the command:
```
R
```
now you can install ggplotw with:
```
install.packages("ggplot2")
```
Finally you can test the ggplot installation with:
```
library(ggplot2)
```

### Installing the required Python libraries
ClusterScan needs few other Python libraries which can be easily installed via [pip](https://pip.pypa.io/en/stable/installing/) from [PyPI](https://pypi.python.org/pypi). As first step install python setuptools and pip on your distribution:
```
sudo apt install python-pip
sudo apt install python-setuptools
```
Then install the required dependencies via pip typing:
```
pip install pandas==0.19.1
pip install rpy2==2.8.3
pip install pybedtools==0.7.8
```
You can find more informations at: 
- [pandas](https://pandas.pydata.org/)
- [pybedtools](https://daler.github.io/pybedtools/)
- [rpy2](https://rpy2.readthedocs.io/en/version_2.8.x/)

### Installing ClusterScan
Now, you can download and extract ClusterScan in your preferred directory:
```
wget https://github.com/pyrevo/ClusterScan/archive/master.zip
unzip master.zip
```
Finally, you can add ClusterScan in your path by editing your bashrc file:
```
sudo nano ~/.bashrc
export PATH=$PATH:/<path-to-clusterscan>/ClusterScan-master
```
Save and exit. To test the clusterscan installation you must open a new terminal and type:
```
clusterscan.py -h
```
If you visualize the ClusterScan help page, you can proceed to run your first analysis by following the [ClusterScan Tutorial](https://github.com/pyrevo/ClusterScan/wiki/ClusterScan-Tutorial).

We have tested ClusterScan on [Ubuntu](https://www.ubuntu.com/) (v14.04LTS+). ClusterScan is distributed under the [GNU General Public License (GPL) Version 3](https://www.gnu.org/licenses/gpl-3.0.en.html).

## Options:
ClusterScan provides different options in order to finely tune the cluster search. Some of them are algorithm-specific whereas other are in common between clusterdist and clustermean:

```
Usage:
  clusterscan.py clusterdist FEATURES ANNOTATION [-o PATH] [-a NAME] [-c LIST] [--info FILE] [--singletons] [-n=<n>] [-d=<bp>]
  clusterscan.py clustermean FEATURES ANNOTATION [-o PATH] [-a NAME] [-c LIST] [--info FILE] [--singletons] [-n=<n>] [-w=<bp>] [-s=<bp>] [-k=<n>] [-e=<n>]
  clusterscan.py (-h | --help)
  clusterscan.py --version

clusterdist:
  -d, --dist=<bp>	Maximum distance between features in bp [default: 500000].

clustermean:
  -w, --window=<bp>	Window size [default: 500000].
  -s, --slide=<bp>	Sliding size [default: 250000].
  -k, --seed=<n>	Number of standard deviations to identify a window which serves as the beginning of the cluster [default: 3].
  -e, --extension=<n>	Number of standard deviations to identify the window(s) which serve to extend the cluster [default: 2].

shared options:
  -o, --output PATH	Specify output path [default: current directory].
  -a, --analysis NAME	Specify optional analysis name for output files.
  -n, --nf=<n>		Minimum number of features per cluster [default: 2].
  -c, --category LIST   Comma separated list of one or more specific categories to be analyzed [e.g. PF00001,PF00002].
                        Useful when you need to perform the analysis only for specific categories in the ANNOTATION file.
  --info FILE           Specify optional file to describe categories.
  --singletons          Identify singletons after clusters and bystanders annotation.
```

An example of execution with clusterdist:
```
clusterscan.py clusterdist my_genes.bed my_categories.txt -d 250000 -a analysis_01
```

An example of execution with clustermean:
```
clusterscan.py clustermean my_genes.bed my_categories.txt --info my_descriptions.txt -w 10000 -s 5000
```

## Input/Output:
ClusterScan is composed by three files: _clusterscan.py_, _algos.py_ and _docopt.py_. The first represent the tool itself, the second stores the functions belonging to the two algorithms described above, and the latter is a dependence that provides an interface for the command-line menu. The three files need to be in the same folder.

ClusterScan requires two mandatory files and a third optional file in input. The first mandatory file is represented by a six-field bed file (FEATURES) containing the location of the features for which to build cluster. It must follow the organization of the [BED file format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1). The second file is a two-column table (ANNOTATION) containing the name of the features as first column (must correspond to the name of the feature in the FEATURES file) and a string representing the categorical information used in the cluster search as second column. This "category" can be represented for esample by Gene Ontology classes, KEGG, Pfam accessions etc. The user can also search for clusters based on a custom made series of categories describing them using an optional additional two-column tab-delimited .txt file (category, description) that can be read by the program (p: --info). Moreover, the user can pass to the program a comma separated list of categories on which perform the analysis (p: --category) in the case in which he/she does not want to analyze all the categories in the ANNOTATION file. It is not uncommon that a feature is associated with different accessions at the same time and, in this case, the feature can be part of multiple different clusters. The user can also choose to identify features which exist as _singletons_ (p: --singletons). Those represent features belonging to a certain category, for which at least one cluster was annotated, that are outside of any of the corresponding clusters.

ClusterScan gives six or seven files in output depending on whether the user has chosen or not to identify the singletons:

|File|Description|
|---|---|
|	clusters.tsv | stores the coordinates of all the clusters found. It contains an unique ID for each cluster; the category to which the cluster belongs; the chromosome/scaffold (chr) on which it resides; its start and end coordinates; the number of features within the cluster (n_features) and the number of features which overlap the cluster but belong to a different category (n_bystanders). |
|		clusters.bed | a bed version of the previous file for compatibility with other tools for downstream analysis. Fields respect the bed format and are (from first to last column): chromosome/scaffold on which the cluster resides; its start/end coordinates using a 0-based start and 1-based end system of coordinates; the cluster ID; the number of features within the cluster; the strand; the category describing the cluster. |
|		summary.tsv | contains a category based summary of the cluster analysis with the category; the number of clusters found for each category (n_clusters); the total number of features found for clusters belonging to each category (n_ft); the number of bystanders (n_bs); the number of features and bystanders found in the clusters which contains the minimal and the maximal number of these features (max_ft and min_ft / max_bs and min_bs) for each category. |
|		features.tsv | a list of features found to be in overlap with clusters. It contains the chromosome/scaffold (chr) on which the feature resides; its start and end coordinates; the feature name; the strand; the ID of the cluster to which it overlaps; the category to which the cluster (and the feature) belongs. |
|		bystanders.tsv | a list of bystanders found to be in overlap with clusters. It contains the chromosome/scaffold (chr) on which the feature resides; its start and end coordinates; the feature name; the strand; the ID of the cluster to which it overlaps; the category to which the cluster (but no the feature) belongs. |
|		distribution.pdf | is an histogram which shows the distribution of features in the top-10 categories by number of features. |
|		singletons.tsv | a list of singletons found to be out from the annotated clusters for each category. It contains the chromosome/scaffold (chr) on which it resides; its start and end coordinates; the feature name; the strand; the category to which the feature belongs. |

The program assigns an unique ID for each cluster during the analysis, joining the “C” letter with a progressive number starting from 1. The user can trace back a specific cluster at any time in the different output files simply looking at the ID reported in the corresponding column.

## Other resources:
You can contact us by opening a new issue on GitHub or asking a new question on our dedicated [mailing list](https://groups.google.com/forum/#!forum/clusterscan-support). Take a look at our Wiki [here](https://github.com/pyrevo/ClusterScan/wiki).
