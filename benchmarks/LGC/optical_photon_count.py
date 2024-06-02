import os
import argparse
import awkward as ak
import uproot
import numpy as np
from matplotlib import pyplot as plt


mc_cols = [
    'MCParticles.generatorStatus',
    'MCParticles.momentum.x',
    'MCParticles.momentum.y',
    'MCParticles.momentum.z',
    'MCParticles.PDG'
]

lgc_cols = [
    'LightGasCherenkovHits.cellID',
    'LightGasCherenkovHits.time',
]

qe = 0.15

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Count optical photons on the PMT Arrays')
    parser.add_argument('sim_file', type=str,
                        help='path to simulation output (root file).')
    parser.add_argument('-o', '--output-dir', type=str,
                        default='.',
                        help='output directory.')
    parser.add_argument('-n', dest='nevents', type=int,
                        default=-1,
                        help='number of event to process.')
    parser.add_argument('-s', '--entry-start', type=int,
                        default=-1,
                        help='event number to start')
    args = parser.parse_args()

    # determine event range
    entry_start = None
    entry_stop = None
    if args.nevents > 0:
        entry_stop = args.nevents
    if args.entry_start > 0:
        entry_start = args.entry_start
        entry_stop += entry_start

    events = uproot.open(args.sim_file)['events']

    # generated particles
    mc_pars = events.arrays(mc_cols, entry_stop=entry_stop, library='ak')
    first_par_mask = mc_pars['MCParticles.generatorStatus'] == 1
    mc_pars = mc_pars[first_par_mask]

    momentum = np.sqrt(mc_pars['MCParticles.momentum.x']**2 + mc_pars['MCParticles.momentum.y']**2 + mc_pars['MCParticles.momentum.z']**2)
    phi = ak.flatten(np.arctan2(mc_pars['MCParticles.momentum.x'], mc_pars['MCParticles.momentum.y']))
    theta = ak.flatten(np.arccos(mc_pars['MCParticles.momentum.z']/momentum))

    # LGC readouts
    lgc_hits = events.arrays(lgc_cols, entry_stop=entry_stop, library='ak')
    raw_counts = ak.count(lgc_hits['LightGasCherenkovHits.cellID'], axis=-1)
    qe_counts = raw_counts*qe
    counts_err = np.sqrt(qe * (1. - qe) * raw_counts)

    # plots
    fig, ax = plt.subplots(figsize=(8, 6), dpi=160)
    ax.hist(qe_counts, bins=np.linspace(0, 60, 61), ec='k')
    fig.savefig(os.path.join(args.output_dir, 'oph_counts.png'))
