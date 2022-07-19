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
    parser.add_argument('-e', dest='event', type=int, default=0, help='event number to plot')
    args = parser.parse_args()

    # read a specific event, and flatten data
    rdf = ROOT.RDataFrame("events", args.file)
    rdf = rdf.Range(args.event, args.event + 1)
    df = flatten_collection(rdf, args.branch)
    df.rename(columns={c: c.replace(args.branch + '.', '') for c in df.columns}, inplace=True)

    dfm = df.groupby(['position.x', 'position.y'])['energy'].sum()

    cmap = truncate_colormap(plt.get_cmap('jet'), 0.1, 0.9)
    energy_range = (1, 1000)

    # convert mm->cm, GeV->MeV
    patches, colors = [], []
    for px, py, ene in dfm.reset_index().values:
        if ene*1000. > energy_range[0]:
            w = np.log(ene*1000.)
            wmin, wmax = np.log(energy_range)
            patches.append(mpatches.RegularPolygon((px/10., py/10.), 6, 6.25*1.25, np.pi/2.))
            colors.append(cmap((w - wmin)/(wmax - wmin)))

    fig, axs = plt.subplots(1, 2, figsize=(13, 12), dpi=160, gridspec_kw={'wspace':0., 'width_ratios': [12, 1]})
    ax = axs[0]
    ax.add_collection(PatchCollection(patches, facecolor=colors))
    ax.add_patch(mpatches.Circle((0, 0), 90., facecolor='none', edgecolor='k', linewidth=2, ls=':'))
    ax.add_patch(mpatches.Circle((0, 0), 230., facecolor='none', edgecolor='k', linewidth=2, ls=':'))
    ax.set_xlim(-240, 240)
    ax.set_ylim(-240, 240)
    ax.set_ylabel('y (cm)', fontsize=24)
    ax.set_xlabel('x (cm)', fontsize=24)
    ax.tick_params(labelsize=22)
    sm = cm.ScalarMappable(norm=matplotlib.colors.LogNorm(vmin=energy_range[0], vmax=energy_range[1]), cmap=cmap)
    cb = plt.colorbar(sm, cax=axs[1], shrink=0.85, aspect=1.2*20)
    cb.ax.tick_params(labelsize=24)
    cb.ax.get_yaxis().labelpad = 10
    cb.ax.set_ylabel('Energy Deposit (MeV)', rotation=90, fontsize=24)
    fig.savefig('test.png')

