#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Description:
  ClusterScan, search for clusters of features in a given annotation.

Usage:
  clusterscan.py clusterdist FEATURES ANNOTATION [-o PATH] [-a NAME] [--info FILE] [-n=<n>] [-d=<bp>]
  clusterscan.py clustermean FEATURES ANNOTATION [-o PATH] [-a NAME] [--info FILE] [-n=<n>] [-w=<bp>] [-s=<bp>] [-k=<n>] [-e=<n>]
  clusterscan.py (-h | --help)
  clusterscan.py --version

Options:
  -h, --help                        Show this screen.
  -o, --output PATH                 Specify output path [default: ./].
  -a, --analysis NAME               Specify optional analysis name for output files.
  -n, --nf=<n>                      Minimum number of features per cluster [default: 2].
  -d, --dist=<bp>                   Maximum distance between features in bp [default: 500000].
  -w, --window=<bp>                 Window size [default: 500000].
  -s, --slide=<bp>                  Sliding size [default: 250000].
  -k, --seed=<n>                    Number of standard deviations to identify a window which serves as the beginning of the cluster [default: 3].
  -e, --extension=<n>               Number of standard deviations to identify the window(s) which serve to extend the cluster [default: 2].
  --info FILE                       Specify optional file to describe accessions.
  --version                         Show program version.
"""


import warnings

import pandas as pd
import pybedtools
from docopt import docopt
from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

from algos import *

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import rpy2.robjects.lib.ggplot2 as ggplot2


def input_tester(file_path):
    """Check for the presence of input files."""
    try:
        open(file_path)
    except IOError:
        raise SystemExit('Unable to open %s, file does not exist!' % (file_path.split('/')[-1]))
    else:
        pass


def options_tester(option, n, string):
    """Check if rules for parameters are respected."""
    if option < n:
        raise ValueError(string)
    else:
        pass

def rpy2_plotter(anno, clusters, name):
    """Plot genes distribution in clusters using ggplot2 from R."""
    pandas2ri.activate()
    grdevices = importr('grDevices')
    rprint = robjects.globalenv.get("print")

    anno = anno.sort_values(by="n_ft", ascending=False)
    anno = anno.head(n=10)
    accession = anno["ACC"].tolist()
    clusters = clusters[clusters["ACC"].isin(accession)]
    clusters = pandas2ri.py2ri(clusters)

    pp = ggplot2.ggplot(clusters) + ggplot2.aes_string(x="n_features") + ggplot2.geom_histogram(binwidth=1) + ggplot2.facet_wrap(robjects.Formula("~ACC"), ncol=5) + ggplot2.labs(x="Number of Features", y="Number of Clusters", title="Clusters distribution")

    grdevices.pdf(file=name, width=11.692, height=8.267)
    rprint(pp)
    grdevices.dev_off()


def main():
    # test for input files availability
    input_tester(arguments['FEATURES'])
    input_tester(arguments['ANNOTATION'])

    if arguments['--info'] is None:
        pass
    else:
        input_tester(arguments['--info'])

    # clusters can't contain less than 2 features
    error1 = "Minimum number of features per cluster must be a number higher than 1!"
    options_tester(int(arguments['--nf']), 2, error1)

    # window size can't be lower than sliding size
    error2 = "Sliding size can't be higher than window size!"
    options_tester(int(arguments['--window']), int(arguments['--slide']), error2)

    # window size can't be lower than sliding size
    error3 = "Seed or extension can't be a number lower than 1!"
    options_tester(int(arguments['--seed']), 1, error3)
    options_tester(int(arguments['--extension']), 1, error3)

    # build database
    feat = pd.read_table(arguments['FEATURES'], header=None, usecols=range(6))
    anno = pd.read_table(arguments['ANNOTATION'], header=None)

    feat.columns = ['chr', 'start', 'end', 'name', 'score', 'strand']
    anno.columns = ['name', 'ACC']
    #anno['ACC'] = anno['ACC'].fillna("Unknown")

    # pdtable stores genes annotation and corresponding accessions
    pdtable = pd.merge(feat, anno, on='name', how='outer')
    pdtable['ACC'] = pdtable['ACC'].fillna("Unknown")
    pdtable = pdtable[pd.notnull(pdtable['ACC'])]
    pdtable = pdtable[pd.notnull(pdtable['chr'].astype(str))]
    pdtable[['start', 'end']] = pdtable[['start', 'end']].astype(int)
    pdtable = pdtable.drop_duplicates(['name', 'ACC'])
    #movq print str(pdtable)
    all_features = pdtable
    pdtable = pdtable[pdtable['ACC'] != "Unknown"]

    # list unique accessions
    l = list(pdtable.ACC.unique())
    # inizialize empty table to be filled with clusters
    table = pd.DataFrame()

    # choose the algorithm
    if arguments['clusterdist'] is True:
        print "ClusterScan is running with clusterdist..."

        # movq: arguments should be changed
        table = do_clusterdist(l,pdtable, table, arguments)
        #print resTable
        # for each accession merge features getting clusters
        #for ACC in l:
        #    df = pdtable[pdtable.ACC == ACC]
        #    BEDtools_object = pybedtools.BedTool().from_dataframe(df).sort()

        #    try:
        #        merge = BEDtools_object.merge(d=int(arguments['--dist']), c=4, o="count_distinct")
        #    except:
        #        continue

        #    df = pd.read_table(merge.fn, header=None)
        #    df[4] = ACC
        #    table = table.append(df)
    else:
        print "ClusterScan is running with clustermean..."

        table = do_clustermean(l, pdtable, table, arguments)

        ## for each chromosome determine its length as last feature end
        #loc = list(pdtable.chr.unique())
        #chr_len = []

        #for chr in loc:
        #    df = pdtable[pdtable.chr == chr]
        #    chr_len.append((chr, 0, max(df.end)))

        #windows = []
        #window_maker(windows, chr_len, int(arguments['--window']), int(arguments['--slide']))
        #win_bed = pybedtools.BedTool(windows)

        ## for each accession compute clusters
        #for ACC in l:
        #    df = pdtable[pdtable.ACC == ACC]
        #    BEDtools_object = pybedtools.BedTool().from_dataframe(df)

        #    # intersect features to windows
        #    try:
        #        intersect_bed = win_bed.intersect(BEDtools_object, c=True)
        #    except:
        #        continue

        #    df = pd.read_table(intersect_bed.fn, header=None, dtype={0: str})
        #    df[4] = ACC

        #    # compute mean and stdv feature density per-window
        #    mean = df[3].mean()
        #    stdv = df[3].std()

        #    multi1 = mean + (int(arguments['--seed'])*stdv)
        #    multi2 = mean + (int(arguments['--extension'])*stdv)

        #    # extract seeds and try to extend them
        #    seed_list = df[df[3] >= multi1].index.tolist()
        #    extended_seed = []
        #    seed_extender(extended_seed, seed_list, intersect_bed, multi2)

        #    pre_clusters = pybedtools.BedTool(extended_seed)
        #    features_in_clusters = BEDtools_object.intersect(pre_clusters, wa=True)
        #    final_list = []
        #    cluster_composer(final_list, pre_clusters, features_in_clusters)

        #    final_clusters = pybedtools.BedTool(final_list)
        #    final_clusters = final_clusters.intersect(BEDtools_object, c=True)
        #    final_clusters = pd.read_table(final_clusters.fn, header=None)
        #    final_clusters[5] = ACC
        #    table = table.append(final_clusters)

    # generate cluster table and filter it
    table.columns = ["chr", "start", "end", "n_features", "ACC"]
    table = table.sort_values(["ACC", "chr"], ascending=[True, True])
    table = table[table["n_features"] >= int(arguments['--nf'])]
    table = table.sort_values(by=["ACC"], ascending=[True])
    # table['ID'] = range(1, len(table) + 1)
    table['ID'] = ["C"+str(i) for i in range(1, len(table) + 1)]

    if table.empty:
        print "ClusterScan didn't found any cluster!"
        exit()

    #print str(table)
    # generate output of clusters in BED format
    bed = table.copy()
    bed["strand"] = "+"
    bed = pybedtools.BedTool().from_dataframe(bed[[0, 1, 2, 5, 3, 6, 4]]).sort()

    # generate table of features by intersect feature with clusters
    all_features_bed = pybedtools.BedTool().from_dataframe(all_features)
    clusters = pybedtools.BedTool().from_dataframe(table)
    features = all_features_bed.intersect(clusters, wb=True)
    features = pd.read_table(features.fn, header=None, dtype={0: str, 7: str})
    cl_features = features[features[6] == features[11]]
    cl_features = cl_features[[0, 1, 2, 3, 4, 5, 12, 11]]
    cl_features.columns = ["chr", "start", "end", "name", "score",
                           "strand", "cluster_ID", "cluster_ACC"]

    # generate table of bystanders
    bystanders = features[features[6] != features[11]]

    # control for bystander = 0 (when program run with 1 accession)
    if bystanders.empty:
        table = table[[5, 4, 0, 1, 2, 3]]
        table["n_bystanders"] = 0
    else:
        bystanders = bystanders[[0, 1, 2, 3, 4, 5, 12, 11]]
        bystanders.columns = ["chr", "start", "end", "name", "score",
                              "strand", "cluster_ID", "cluster_ACC"]
        # prevent bystanders with 2+ different ACC to be counted twice
        bystanders = bystanders.drop_duplicates(['name', 'cluster_ID'])
        # prevent features with 2+ different ACC to be bystanders in theyr cluster
        bs_merge = pd.merge(bystanders, cl_features, how='outer', indicator=True)
        bystanders = bs_merge.ix[bs_merge._merge == 'left_only']
        bystanders = bystanders.drop(bystanders.columns[8], axis=1)
        # count bystanders number
        bs_count = bystanders.groupby('cluster_ID').count().reset_index()
        bs_count = bs_count[[0, 1]]
        bs_count.columns = ["ID", "n_bystanders"]
        table = table.merge(bs_count, on="ID", how='outer')
        table = table[[5, 4, 0, 1, 2, 3, 6]]
        table["n_bystanders"].fillna(0, inplace=True)
        table["n_bystanders"] = table["n_bystanders"].astype(int)

    # generate summary table
    summary = table.drop(table.columns[[0, 2, 3, 4]], axis=1)
    # calculate total number of clusters per-accession
    n_clusters = summary.groupby("ACC").count().reset_index()
    n_clusters = n_clusters.drop("n_bystanders", axis=1)
    n_clusters.columns = ["ACC", "n_clusters"]
    # calculate total number of features and bystanders per-accession
    n_ft_bs = summary.groupby("ACC").sum().reset_index()
    n_ft_bs.columns = ["ACC", "n_ft", "n_bs"]
    # calculate maximum number of features and bystander in cluster
    max_ft_bs = summary.groupby("ACC").max().reset_index()
    max_ft_bs.columns = ["ACC", "max_ft", "max_bs"]
    # calculate minimum number of features and bystander in cluster
    min_ft_bs = summary.groupby("ACC").min().reset_index()
    min_ft_bs.columns = ["ACC", "min_ft", "min_bs"]
    # add accession description if an info file is provided
    if arguments['--info'] is None:
        summary = n_clusters.merge(n_ft_bs, on='ACC').merge(max_ft_bs, on='ACC').merge(min_ft_bs, on='ACC')
    else:
        desc = pd.read_table(arguments['--info'], header=None)
        desc.columns = ["ACC", "DESC"]
        summary = n_clusters.merge(n_ft_bs, on='ACC').merge(max_ft_bs, on='ACC').merge(min_ft_bs, on='ACC').merge(desc, on='ACC')

    # assign file names and save tables as result
    if arguments['--analysis'] is None:
        feat_name = arguments['--output'] + 'features.csv'
        byst_name = arguments['--output'] + 'bystanders.csv'
        clus_name = arguments['--output'] + 'clusters.csv'
        summ_name = arguments['--output'] + 'summary.csv'
        bed_name = arguments['--output'] + 'clusters.bed'
        plot_name = arguments['--output'] + 'distribution.pdf'
    else:
        feat_name = arguments['--output'] + arguments['--analysis'] + '_features.csv'
        byst_name = arguments['--output'] + arguments['--analysis'] + '_bystanders.csv'
        clus_name = arguments['--output'] + arguments['--analysis'] + '_clusters.csv'
        summ_name = arguments['--output'] + arguments['--analysis'] + '_summary.csv'
        bed_name = arguments['--output'] + arguments['--analysis'] + '_clusters.bed'
        plot_name = arguments['--output'] + arguments['--analysis'] + '_distribution.pdf'

    table["start"] += 1
    cl_features.to_csv(feat_name, sep='\t', header=True, index=False)
    bystanders.to_csv(byst_name, sep='\t', header=True, index=False)
    table.to_csv(clus_name, sep='\t', header=True, index=False)
    summary.to_csv(summ_name, sep='\t', header=True, index=False)

    if arguments['--analysis'] is None:
        bed.saveas(bed_name, trackline='track name="%s" description="chr start end ID n_features strand ACC"' % (arguments['FEATURES']))
    else:
        bed.saveas(bed_name, trackline='track name="%s" description="chr start end ID n_features strand ACC"' % (arguments['--analysis']))

    # plot a duistribution for top 10 clusters (per n of features)
    rpy2_plotter(summary, table, plot_name)

    print "Done!"

# program execution
if __name__ == '__main__':
    arguments = docopt(__doc__, version='ClusterScan 0.1.0')
    # print arguments
    main()
