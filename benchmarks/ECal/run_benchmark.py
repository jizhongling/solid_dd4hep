#! /usr/local/bin/python3
'''
    A python script to facilitate the ML benchmarks for e/pi separation with the imaging calorimeter (single particles).
    This process follows the steps below:
    1. Simulation to generate training samples
    2. Study and apply E/p cut to reduce the training samples
    3. Train and test ML models with the "cleaned" (after E/p cut) samples
    4. Benchmark the performance

    Author: Chao Peng (ANL)
    Date: 11/11/2022
'''
import os
import sys
import json
import subprocess
import argparse


SDIR = os.path.dirname(os.path.realpath(__file__))
# {var} is from args
FILE_NAMES = dict(
    gen_script = os.path.join(SDIR, 'scripts', 'gen_particles.py'),
    rec_script = os.path.join(SDIR, 'options', 'faec.py'),

    sim_dir = os.path.join('{outdir}', '{run_type}', 'sim_data'),

    gen_file = os.path.join('{outdir}', '{run_type}', 'sim_data', '{ntag}_gen.hepmc'),
    sim_file = os.path.join('{outdir}', '{run_type}', 'sim_data', '{ntag}_sim.edm4hep.root'),
    rec_file = os.path.join('{outdir}', '{run_type}', 'sim_data', '{ntag}_rec.root'),
)
# default values for argument parser
DEFAULT_COMPACT = os.path.join(
        os.environ.get('DETECTOR_PATH', ''),
        '{}.xml'.format(os.environ.get('DETECTOR_CONFIG', ''))
        )
# defined steps
SCRIPT_STEPS = (
    'sim',      # step 1; simulation to generate samples
)


# simulation and reconstruction
def gen_sim_rec(**kwargs):
    # generate particles
    gen_cmd = [
        'python {gen_script} {gen_file}',
        '-n {nev}',
        '-s {seed}',
        '--angmin {angmin} --angmax {angmax}',
        '--pmin {pmin} --pmax {pmax}',
        '--particles {particles}',
        ]
    gen_cmd = ' '.join(gen_cmd).format(**kwargs).split(' ')
    subprocess.run(gen_cmd)

    # simulation
    sim_cmd = [
        'ddsim --runType batch --part.minimalKineticEnergy 1*TeV --filter.tracker edep0',
        '-v WARNING',
        '--numberOfEvents {nev}',
        # '--physics.list {physics_list}',
        '--inputFiles {gen_file}',
        '--outputFile {sim_file}',
        '--compact {compact}',
        ]
    if 'seed' in kwargs and kwargs['seed'] > 0:
        sim_cmd += ['--random.seed {seed}']
    sim_cmd = ' '.join(sim_cmd).format(**kwargs).split(' ')
    return_code = subprocess.run(sim_cmd).returncode
    print(return_code)
    if return_code is not None and return_code < 0:
        print("ERROR running simulation!")
        exit(return_code)
    subprocess.run(['rootls', '-t', kwargs['sim_file']])

    # reconstruction with juggler
    # export to environment variables (used to pass arguments to the option file)
    run_env = os.environ.copy()
    juggler_vars = [
        'JUGGLER_SIM_FILE {sim_file}',
        'JUGGLER_REC_FILE {rec_file}',
        'JUGGLER_COMPACT_PATH {compact}',
        'JUGGLER_N_EVENTS {nev}',
        ]
    lst = ' '.join(juggler_vars).format(**kwargs).split(' ')
    run_env.update({lst[i]: lst[i + 1] for i in range(0, len(lst), 2)})

    rec_cmd = 'gaudirun.py {rec_script}'.format(**kwargs).split(' ')
    print(rec_cmd)
    return_code = subprocess.run(rec_cmd, env=run_env).returncode
    print(return_code)
    if return_code is not None and return_code < 0:
        print("ERROR running juggler (reconstruction)!")
        exit(return_code)
    process = subprocess.run(['rootls', '-t', kwargs['rec_file']])



if __name__ == '__main__':

    # argument parser
    parser = argparse.ArgumentParser()

    parser.add_argument(
            '-n', '--n-events', type=int,
            dest='nev',
            default=100,
            help='number of events.'
            )
    parser.add_argument(
            '-o', '--outdir', type=str,
            dest='outdir',
            default='sim_output',
            help='output directory.'
            )
    parser.add_argument(
            '-r', '--run-type', type=str,
            dest='run_type',
            default='ecal',
            help='a name specify the run type.'
            )
    parser.add_argument(
            '-t', '--name-tag', type=str,
            dest='ntag',
            default='solid',
            help='a name tag for output files.'
            )
    parser.add_argument(
            '-c', '--compact', type=str,
            dest='compact',
            default=DEFAULT_COMPACT,
            help='path to detector compact file.'
            )
    parser.add_argument(
            '-s', '--seed', type=int,
            default=-1,
            help='random seed to child scripts (only pass it if > 0).'
            )
    parser.add_argument(
            '--batch-size', type=int,
            dest='batch',
            default=100000,
            help='batch size to process data.'
            )
    parser.add_argument(
            '--p-min', type=float,
            dest='pmin',
            default=5.0,
            help='minimum momentum of particles.'
            )
    parser.add_argument(
            '--p-max', type=float,
            dest='pmax',
            default=5.0,
            help='maximum momentum of particles.'
            )
    parser.add_argument(
            '--angle-min', type=float,
            dest='angmin',
            default=5,
            help='minimum scattering angle of particles.'
            )
    parser.add_argument(
            '--angle-max', type=float,
            dest='angmax',
            default=25,
            help='maximum scattering angle of particles.'
            )
    parser.add_argument(
            '--particles', type=str,
            default='electron',
            help='partcile names, separated by \",\".'
            )
    parser.add_argument(
            '--steps', type=str,
            default=', '.join(SCRIPT_STEPS),
            help='FOR DEV: choose the steps to be executed ({}).'.format(', '.join(SCRIPT_STEPS))
            )

    args = parser.parse_args()
    kwargs = vars(args)

    # prepare
    steps = [p.strip() for p in args.steps.split(',')]

    # make dirs, add paths to kwargs
    FILE_NAMES.update({key: val.format(**kwargs) for key, val in FILE_NAMES.items()})
    for key, val in FILE_NAMES.items():
        if key.endswith('_dir'):
            os.makedirs(val, exist_ok=True)
    kwargs.update(FILE_NAMES)

    # simulation for benchmark samples
    if SCRIPT_STEPS[0] in steps:
        gen_sim_rec(**kwargs)

    # save run information, combine runs with the same run_type
    run_data = {args.run_type: {args.ntag: kwargs}}
    try:
        with open(os.path.join(args.outdir, 'result.json'), 'r') as f:
            run_data = json.load(f)
            run_info = run_data.get(args.run_type, {})
            run_info.update({args.ntag: kwargs})
            run_data[args.run_type] = run_info
    except (FileNotFoundError, json.decoder.JSONDecodeError):
        pass

    with open(os.path.join(args.outdir, 'result.json'), 'w') as f:
        f.write(json.dumps(run_data, sort_keys=True, indent=4))

