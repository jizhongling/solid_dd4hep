#! /usr/local/bin/python3
'''
    A python script to facilitate the ML training benchmark for imaging Calorimeter:
    simulation -> reconstruction -> prepare tagging data for tensorflow
    Author: Chao Peng (ANL)
    Date: 06/20/2021
'''
import os
import sys
import subprocess
import argparse


default_compact = os.path.join(os.environ.get('DETECTOR_PATH',  os.environ.get('DETECTOR_PATH', '')),
        '{}.xml'.format(os.environ.get('JUGGLER_DETECTOR', '')))
parser = argparse.ArgumentParser()
parser.add_argument('-n', '--numberOfEvents', dest='nev', type=int, default=100, help='Number of events to process.')
parser.add_argument('-t', '--nametag', type=str, default='IMCAL_ML', help='Name tag for output files.')
parser.add_argument('--seed', type=int, default=-1, help='Random seed to scripts.')
parser.add_argument('--process', type=str, default='sim, rec', help='Processes to be executed (sim, rec).')
parser.add_argument('--nlayers', dest='nlayers', type=int, default=9, help='number of layers in ML data.')
parser.add_argument('--nhits', dest='nhits', type=int, default=20, help='number of hits in ML data.')
parser.add_argument('--particles', type=str, default='electron', help='Partcile names, separated by \",\".')
parser.add_argument('--angmin', type=float, default=12, help='Minimum momentum of particles.')
parser.add_argument('--angmax', type=float, default=15, help='Maximum momentum of particles.')
parser.add_argument('--pmin', type=float, default=4, help='Minimum momentum of particles.')
parser.add_argument('--pmax', type=float, default=10, help='Maximum momentum of particles.')
parser.add_argument('--compact', type=str, default=default_compact, help='Path to detector compact file.')
parser.add_argument('--combine-method', type=str, default='interlayer', help='Path to detector compact file.')
parser.add_argument('--physics-list', type=str, default='FTFP_BERT', help='Path to detector compact file.')

args = parser.parse_args()
kwargs = vars(args)

for mdir in ['gen_data', 'sim_data', 'rec_data']:
    os.makedirs(mdir, exist_ok=True)

gen_file = os.path.join('gen_data', '{nametag}_{pmin}_{pmax}.hepmc'.format(**kwargs))
sim_file = os.path.join('sim_data', '{nametag}_{pmin}_{pmax}.root'.format(**kwargs))
rec_file = os.path.join('rec_data', '{nametag}_{pmin}_{pmax}.root'.format(**kwargs))
tag_dir = os.path.join('tag_data', '{nametag}_{pmin}_{pmax}'.format(**kwargs))
os.makedirs(tag_dir, exist_ok=True)

procs = [p.strip() for p in args.process.split(',')]
sdir = os.path.dirname(os.path.realpath(__file__))

if 'sim' in procs:
    # generate particles
    gen_cmd = ['python', os.path.join(sdir, 'scripts', 'gen_particles.py'), gen_file,
            '-n', '{}'.format(args.nev),
            '-s', '{}'.format(args.seed),
            '--angmin', '{}'.format(args.angmin), '--angmax', '{}'.format(args.angmax),
            '--pmin', '{}'.format(args.pmin), '--pmax', '{}'.format(args.pmax),
            '--particles', args.particles]
    subprocess.run(gen_cmd)
    # simulation
    sim_cmd = ['npsim',
            '--part.minimalKineticEnergy', '1*TeV',
            '--numberOfEvents', '{}'.format(args.nev),
            '--runType', 'batch',
            # '--physics.list', args.physics_list,
            '--inputFiles', gen_file,
            '--outputFile', sim_file,
            '--compact', args.compact,
            '-v', 'WARNING']
    if args.seed > 0:
        sim_cmd += ['--random.seed', args.seed]
    return_code = subprocess.run(sim_cmd).returncode
    print(return_code)
    if return_code is not None and return_code < 0:
        print("ERROR running simulation!")
        exit(return_code)
    subprocess.run(['rootls', '-t', sim_file])


if 'rec' in procs:
    # export to environment variables (used to pass arguments to the option file)
    run_env = os.environ.copy()
    run_env.update({
        'JUGGLER_SIM_FILE': sim_file,
        'JUGGLER_REC_FILE': rec_file,
        'JUGGLER_COMPACT_PATH': args.compact,
        'JUGGLER_N_EVENTS': str(args.nev),
        'IMCAL_ML_IMG_NLAYERS': str(args.nlayers),
        'IMCAL_ML_NHITS': str(args.nhits),
        'IMCAL_ML_COMBINE': str(args.combine_method),
    })

    juggler_xenv = os.path.join(os.environ.get('JUGGLER_INTALL_PREFIX', '../local'), 'Juggler.xenv')
    rec_cmd = [
        # 'xenv', '-x', juggler_xenv,   # v35+ do not need xenv anymore
        'gaudirun.py', os.path.join(sdir, 'options', 'faec.py')
    ]
    return_code = subprocess.run(rec_cmd, env=run_env).returncode
    print(return_code)
    if return_code is not None and return_code < 0:
        print("ERROR running juggler (reconstruction)!")
        exit(return_code)
    process = subprocess.run(['rootls', '-t', rec_file])


