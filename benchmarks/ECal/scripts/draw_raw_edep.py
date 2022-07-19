import ROOT
import os
import pandas as pd
import argparse
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection


# helper function to truncate color map (for a better view from the rainbow colormap)
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


# read from RDataFrame and flatten a given collection, return pandas dataframe
def flatten_collection(rdf, collection, cols=None, event_colname='event'):
    if not cols:
        cols = [str(c) for c in rdf.GetColumnNames() if str(c).startswith('{}.'.format(collection))]
    else:
        cols = ['{}.{}'.format(collection, c) for c in cols]
    if not cols:
        print('cannot find any branch under collection {}'.format(collection))
        return pd.DataFrame()

    data = rdf.AsNumpy(cols)
    # flatten the data, add an event id to identify clusters from different events
    evns = []
    for i, vec in enumerate(data[cols[0]]):
        evns += [i]*vec.size()
    for n, vals in data.items():
        # make sure ints are not converted to floats
        typename = vals[0].__class__.__name__.lower()
        dtype = np.int64 if 'int' in typename or 'long' in typename else np.float32
        # type safe creation
        data[n] = np.asarray([v for vec in vals for v in vec], dtype=dtype)
    # build data frame
    dfp = pd.DataFrame({c: pd.Series(v) for c, v in data.items()})
    dfp.loc[:, event_colname] = evns
    return dfp


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Visualize the cluster from analysis')
    parser.add_argument('file', type=str, help='path to root file')
    parser.add_argument('branch', type=str, help='name of data branch')
    args = parser.parse_args()

    # read a specific event, and flatten data
    rdf = ROOT.RDataFrame("events", args.file)
    df = flatten_collection(rdf, args.branch)
    df.rename(columns={c: c.replace(args.branch + '.', '') for c in df.columns}, inplace=True)

    dfm = df.groupby(['event'])['energy'].sum()

    fig, ax = plt.subplots(figsize=(16, 9), dpi=160)
    ax.hist(dfm.values, bins=np.linspace(1.1, 1.5, 100), align='mid', ec='k')
    ax.set_xlabel('Edep (GeV)', fontsize=24)
    ax.tick_params(labelsize=24)
    ax.set_axisbelow(True)
    ax.grid(linestyle=':')

    fig.savefig('test2.png')

