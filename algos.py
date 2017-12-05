import string
import pybedtools
import pandas as pd

def cluster_composer(pre_cluster_object, pre_cluster_intersection):
    final_list = []
    """Find real feature's positions."""
    tmp = []
    for item in pre_cluster_object:
        scaffold = str(item[0])
        for line in pre_cluster_intersection:
            line = str(line)
            if line.split()[0] == scaffold and (int(item[1]) <= int(line.split()[1]) <= int(item[2]) and int(item[1]) <= int(line.split()[2]) <= int(item[2])):
                tmp.append(int(line.split()[1]))
                tmp.append(int(line.split()[2]))
            else:
                continue
        final_list.append((scaffold, min(tmp), max(tmp)))
        tmp = []
    return final_list

def seed_extender(new_list, indexes, intersection, limit):
    """Extends culster's seeds."""
    for index in indexes:
        cluster_pos = []
        scaffold = intersection[index][0]
        right_step = index
        left_step = index
        while int(intersection[left_step][3]) >= limit and intersection[left_step][0] == scaffold:
            if int(intersection[left_step][1]) not in cluster_pos:
                cluster_pos.append(int(intersection[left_step][1]))
            if (left_step - 1) >= 0:
                left_step -= 1
            else:
                break
        while int(intersection[right_step][3]) >= limit and intersection[right_step][0] == scaffold:
            if int(intersection[right_step][2]) not in cluster_pos:
                cluster_pos.append(int(intersection[right_step][2]))
            if (right_step + 1) <= max(indexes):
                right_step += 1
            else:
                break
        cluster = (str(scaffold), min(cluster_pos), max(cluster_pos))
        if cluster not in new_list:
            new_list.append(cluster)
    return new_list

def do_clusterdist(accList, pdTbl, tbl, sargs):
    for ACC in accList:
        df = pdTbl[pdTbl.ACC == ACC]
        BEDtools_object = pybedtools.BedTool().from_dataframe(df).sort()

        try:
            merge = BEDtools_object.merge(d=int(sargs['--dist']), c=4, o="count_distinct")
        except Exception as e:
            continue

        df = pd.read_table(merge.fn, header=None)
        df[4] = ACC
        tbl = tbl.append(df)


    return tbl

def do_clustermean(accList, pdTbl, tbl, sargs):
    loc = list(pdTbl.chr.unique())
    chr_len = []

    for chr in loc:
        df = pdTbl[pdTbl.chr == chr]
        chr_len.append((chr, 0, max(df.end)))

    windows = []
    window_maker(windows, chr_len, int(sargs['--window']), int(sargs['--slide']))
    win_bed = pybedtools.BedTool(windows)

    # for each accession compute clusters
    for ACC in accList:
        df = pdTbl[pdTbl.ACC == ACC]
        BEDtools_object = pybedtools.BedTool().from_dataframe(df)

        # intersect features to windows
        try:
            intersect_bed = win_bed.intersect(BEDtools_object, c=True)
        except:
            continue

        df = pd.read_table(intersect_bed.fn, header=None, dtype={0: str})
        df[4] = ACC

        # compute mean and stdv feature density per-window
        mean = df[3].mean()
        stdv = df[3].std()

        multi1 = mean + (int(sargs['--seed'])*stdv)
        multi2 = mean + (int(sargs['--extension'])*stdv)

        # extract seeds and try to extend them
        seed_list = df[df[3] >= multi1].index.tolist()
        extended_seed = []
        seed_extender(extended_seed, seed_list, intersect_bed, multi2)

        pre_clusters = pybedtools.BedTool(extended_seed)
        features_in_clusters = BEDtools_object.intersect(pre_clusters, wa=True)
        
        final_list = cluster_composer(pre_clusters, features_in_clusters)

        try:
            final_clusters = pybedtools.BedTool(final_list)
            final_clusters = final_clusters.intersect(BEDtools_object, c=True)  
            final_clusters = pd.read_table(final_clusters.fn, header=None)
            final_clusters[5] = ACC
            tbl = tbl.append(final_clusters)
        except Exception as e:
            pass

    return tbl

def window_maker(list_name, filled_list, window_size, slide_size):
    """Make a bed file of sliding windows."""
    for scaffold, start, end in filled_list:
        width = window_size
        step = slide_size

        if width <= end:
            list_name.append((scaffold, start, width))
        else:
            list_name.append((scaffold, start, end))

        while width <= end:
            start += step
            width += step
            if width >= end:
                list_name.append((scaffold, start, end))
            else:
                list_name.append((scaffold, start, width))
    return list_name
