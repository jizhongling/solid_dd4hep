import os
import argparse
import awkward as ak
import uproot
import numpy as np
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
from scipy import interpolate
from matplotlib import pyplot as plt


qe_data = np.array([
    (1.907449123, 0.0000),
    (1.937253016, 0.0000),
    (1.968003063, 0.0000),
    (1.999745048, 0.0000),
    (2.032527754, 0.0000),
    (2.066403217, 0.0000),
    (2.101427,    0.0266),
    (2.1376585,   0.0345),
    (2.175161281, 0.0422),
    (2.214003446, 0.0511),
    (2.254258055, 0.0595),
    (2.296003574, 0.0698),
    (2.339324396, 0.0815),
    (2.384311404, 0.1025259891),
    (2.431062608, 0.1343072843),
    (2.47968386,  0.1626440114),
    (2.530289653, 0.1805672418),
    (2.583004021, 0.1942904219),
    (2.637961553, 0.2072868681),
    (2.695308543, 0.2238869711),
    (2.755204289, 0.241019145),
    (2.817822568, 0.2585520407),
    (2.883353326, 0.2726082555),
    (2.952004595, 0.2852769602),
    (3.024004707, 0.296623382),
    (3.099604825, 0.3061133496),
    (3.179081872, 0.3141449187),
    (3.262741921, 0.3210286163),
    (3.350924135, 0.3267333585),
    (3.444005361, 0.3309719776),
    (3.542405514, 0.3337052875),
    (3.646593912, 0.3344068971),
    (3.757096758, 0.3326316333),
    (3.874506031, 0.3301402224),
    (3.999490097, 0.3284878445),
    (4.132806433, 0.3284802486),
    (4.275317,    0.328480254),
    (4.428006893, 0.328480254),
    (4.592007148, 0.328480254),
    (4.768622808, 0.328480254),
    (4.95936772,  0.3284802348),
    (5.166008042, 0.3283322658),
    (5.390617087, 0.3278605334),
    (5.635645136, 0.3284802515),
    (5.90400919,  0.3284802474),
    (6.19920965,  0.3284802531),
    (6.525483842, 0.3284802539),
    (6.888010722, 0.328480254),
    (7.293187824, 0.328480254),
])

mc_cols = [
    'MCParticles.generatorStatus',
    'MCParticles.momentum.x',
    'MCParticles.momentum.y',
    'MCParticles.momentum.z',
    'MCParticles.PDG',
]

lgc_cols = [
    'LightGasCherenkovHits.cellID',
    'LightGasCherenkovHits.time',
    'LightGasCherenkovHits.EDep',
]


def ev2nm(ev):
    return 4.1357e-15*2.99792458e8*1e9 / ev


def nm2ev(nm):
    return 4.1357e-15*2.99792458e8*1e9 / nm


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
    parser.add_argument('--color-map', type=str,
                        default='viridis',
                        help='color map name in matplotlib')
    parser.add_argument('--vmax-theta', type=float,
                        default=None,
                        help='max value for 2d histogram with theta')
    parser.add_argument('--vmin-theta', type=float,
                        default=None,
                        help='min value for 2d histogram with theta')
    parser.add_argument('--vmax-phi', type=float,
                        default=None,
                        help='max value for 2d histogram with phi')
    parser.add_argument('--vmin-phi', type=float,
                        default=None,
                        help='min value for 2d histogram with phi')
    parser.add_argument('--bins-photon', type=str,
                        default="1,120,120",
                        help='optical photon bins in the format of \"<min>,<max>,<nbins>\"')
    parser.add_argument('--bins-theta', type=str,
                        default="5,20,151",
                        help='theta bins in the format of \"<min>,<max>,<nbins>\"')
    parser.add_argument('--bins-phi', type=str,
                        default="-180,180,361",
                        help='phi bins in the format of \"<min>,<max>,<nbins>\"')

    args = parser.parse_args()

    bins_npe = np.linspace(*eval(args.bins_photon))
    bins_theta = np.linspace(*eval(args.bins_theta))
    bins_phi = np.linspace(*eval(args.bins_phi))

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
    phi = ak.to_numpy(ak.flatten(np.arctan2(mc_pars['MCParticles.momentum.x'], mc_pars['MCParticles.momentum.y'])))
    theta = ak.to_numpy(ak.flatten(np.arccos(mc_pars['MCParticles.momentum.z']/momentum)))

    # LGC readouts
    lgc_hits = events.arrays(lgc_cols, entry_stop=entry_stop, library='ak')
    raw_counts = ak.to_numpy(ak.count(lgc_hits['LightGasCherenkovHits.cellID'], axis=-1))
    # interpolate qe data to get a curve, output 0 if extrapolating
    qe_curve = interpolate.interp1d(*qe_data.T, fill_value=(0., 0.), bounds_error=False)
    oph_energies = ak.to_numpy(ak.flatten(lgc_hits['LightGasCherenkovHits.EDep']*1e9))
    oph_qe = qe_curve(oph_energies) > np.random.uniform(0., 1., size=len(oph_energies))
    qe_counts = ak.to_numpy(ak.sum(ak.unflatten(oph_qe, ak.num(lgc_hits['LightGasCherenkovHits.EDep'])), axis=-1))

    # plots
    cmap = mpl.colormaps[args.color_map]

    # quantum efficiency curve
    fig, ax = plt.subplots(figsize=(8, 6), dpi=160)
    ax.plot(qe_data.T[0], qe_data.T[1]*100., '.')
    qe_points = np.linspace(1, 8, 701)
    ax.plot(qe_points, qe_curve(qe_points)*100, '--')
    ax.xaxis.set_minor_locator(mticker.MultipleLocator(0.5))
    ax.set_xlabel('Optical Photon Energy [eV]')
    ax.set_ylabel('Quantum Efficiency [%]')
    ax.set_xlim(1.5, 8)
    secax = ax.secondary_xaxis('top', functions=(ev2nm, nm2ev))
    secax.xaxis.set_major_locator(mticker.FixedLocator([160, 200, 300, 400, 500, 600, 800]))
    secax.xaxis.set_minor_locator(mticker.MultipleLocator(10))
    secax.set_xlabel('Wavelength [nm]')
    fig.savefig(os.path.join(args.output_dir, 'oph_quantum_efficiency.png'))

    # 1D distribution of optical photon counts
    fig, ax = plt.subplots(figsize=(8, 6), dpi=160)
    ax.hist(qe_counts, bins=bins_npe, ec='k')
    ax.set_ylabel('Event Counts')
    ax.set_xlabel('Number of Optical Photons Detected')
    fig.savefig(os.path.join(args.output_dir, 'oph_counts.png'))

    # 2D distribution of optical photon counts vs. theta angle
    fig, ax = plt.subplots(figsize=(8, 6), dpi=160, gridspec_kw=dict(right=0.98))
    h = ax.hist2d(theta/np.pi*180., qe_counts, bins=(bins_theta, bins_npe), cmap=cmap, vmin=args.vmin_theta, vmax=args.vmax_theta)
    cbar = fig.colorbar(h[3], ax=ax, pad=0.01)
    cbar.ax.set_ylabel('Event Counts')
    # ax.set_facecolor(cmap(0))
    # ax.set_axisbelow(True)
    ax.grid(ls=(0, (5, 15)), color='w', lw=0.5)
    ax.set_xlabel('Incident Particle Polar Angle ($^{\circ}$)')
    ax.set_ylabel('Number of Optical Photons Detected')
    fig.savefig(os.path.join(args.output_dir, 'oph_counts_theta.png'))

    # 2D distribution of optical photon counts vs. phi angle
    fig, ax = plt.subplots(figsize=(8, 6), dpi=160, gridspec_kw=dict(right=0.98))
    h = ax.hist2d(phi/np.pi*180., qe_counts, bins=(bins_phi, bins_npe), cmap=cmap, vmin=args.vmin_phi, vmax=args.vmax_phi)
    cbar = fig.colorbar(h[3], ax=ax, pad=0.01)
    cbar.ax.set_ylabel('Event Counts')
    # ax.set_facecolor(cmap(0))
    # ax.set_axisbelow(True)
    ax.grid(ls=(0, (5, 15)), color='w', lw=0.5)
    ax.set_xlabel('Incident Particle Azimuthal Angle ($^{\circ}$)')
    ax.set_ylabel('Number of Optical Photons Detected')
    fig.savefig(os.path.join(args.output_dir, 'oph_counts_phi.png'))

