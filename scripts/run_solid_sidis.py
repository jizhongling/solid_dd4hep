#!/usr/bin/env python
from __future__ import absolute_import, unicode_literals

import os
import time
import logging

import argparse
parser = argparse.ArgumentParser(
     prog='run_solid_sidis',
     description='''This runs the simulation  the but that is okay''',
     epilog='''
     This program should be run in the "compact" directory.
         ''')
parser.add_argument("-v","--verbose", help="increase output verbosity", type=int, default=0)
parser.add_argument("--compact", help="compact detector file",default="solid_sidis.xml")
parser.add_argument("--vis", help="vis true/false", action="store_true",default=False)
parser.add_argument("--ui", help="ui setting tcsh or qt; default=qt", type=str,default="qt",dest="ui")
parser.add_argument("-n","--n-events", help="number of events", type=int, default=0)
parser.add_argument("-s","--n-skip", help="number of events to skip", type=int, default=0)
parser.add_argument("-m","--macro", help="macro to execute", default="macro/vis.mac")
parser.add_argument("-o","--output", help="output file", default=None)
parser.add_argument("-i","--input", help="input data file", default=None)
parser.add_argument("--target-position", nargs=3, help="nominal location of the target in units of mm", default=[0.0,0.0,-3000.0],metavar=("X","Y","Z"))

args = parser.parse_args()

import DDG4
from DDG4 import OutputLevel as Output
from g4units import keV, GeV, mm, ns, MeV

def run():

    if not args.vis:
        os.environ['G4UI_USE_TCSH'] = "0"

    kernel = DDG4.Kernel()
    description = kernel.detectorDescription()
    kernel.loadGeometry(str("file:" + args.compact))
    DDG4.importConstants(description)

    geant4 = DDG4.Geant4(kernel)
    if args.vis > 0:
        geant4.printDetectors()

    n_events = args.n_events
    if args.vis:
        geant4.setupUI(typ=args.ui, vis=True,ui=True, macro=args.macro)
    else:
        geant4.setupUI(typ=args.ui, vis=False,ui=False,macro=False)
        if n_events <= 0:
            n_events = 1
    kernel.NumEvents = n_events
     
    geant4.setupTrackingField(stepper='ClassicalRK4', equation='Mag_UsualEqRhs')

    rndm = DDG4.Action(kernel, 'Geant4Random/Random')
    rndm.Seed = 854321
    rndm.initialize()
    rndm.showStatus()

    run1 = DDG4.RunAction(kernel, 'Geant4TestRunAction/RunInit')
    run1.Property_int = 12345
    run1.Property_double = -5e15 * keV
    run1.Property_string = 'Startrun: Hello_2'
    run1.enableUI()
    kernel.registerGlobalAction(run1)
    kernel.runAction().adopt(run1)

    prt = DDG4.EventAction(kernel, 'Geant4ParticlePrint/ParticlePrint')
    #prt.OutputLevel = Output.INFO
    #prt.OutputType = 3  # Print both: table and tree
    kernel.eventAction().adopt(prt)

    outputfile = args.output
    if outputfile is None:
        outputfile = 'data/solid_sidis_' + time.strftime('%Y-%m-%d_%H-%M')
    #rootoutput = geant4.setupROOTOutput('RootOutput', outputfile)
    #rootoutput.HandleMCTruth = True
    podio = DDG4.EventAction(kernel, 'Geant4Output2Podio/RootOutput', True)
    podio.HandleMCTruth = False
    podio.Control = True
    podio.Output = outputfile
    podio.enableUI()
    kernel.eventAction().adopt(podio)
    #--------------------------------

    gen = DDG4.GeneratorAction(kernel, "Geant4GeneratorActionInit/GenerationInit")
    gen.OutputLevel = 0
    gen.enableUI()
    kernel.generatorAction().adopt(gen)

    inputfile = args.input
    if inputfile is None:
        inputfile = "eg_data/solid.ep-2gluon.composite.run00001-1000.hepmc"
    gen = DDG4.GeneratorAction(kernel, "Geant4InputAction/hepmc3")
    gen.Mask = 0
    gen.Input = "HEPMC3FileReader|" + inputfile 
    gen.OutputLevel = 0  # generator_output_level
    gen.Sync = args.n_skip
    gen.enableUI()
    kernel.generatorAction().adopt(gen)

    gen = DDG4.GeneratorAction(kernel, "Geant4InteractionVertexSmear/SmearVert")
    gen.Mask = 0
    gen.Offset = (args.target_position[0]*mm, args.target_position[1]* mm, args.target_position[2]* mm, 0 * ns)
    gen.Sigma = (3 * mm, 3 * mm, 1 * mm, 0 * ns)
    kernel.generatorAction().adopt(gen)

    #gen = DDG4.GeneratorAction(kernel, "Geant4GeneratorWrapper/GPS")
    #gen.Uses = 'G4GeneralParticleSource'
    ##gen.OutputLevel = Output.WARNING
    ##gen.Mask = 1
    #gen.enableUI()
    #kernel.generatorAction().adopt(gen)

    gen = DDG4.GeneratorAction(kernel, "Geant4InteractionMerger/InteractionMerger")
    #gen.OutputLevel = 0  # generator_output_level
    gen.enableUI()
    kernel.generatorAction().adopt(gen)

    gen = DDG4.GeneratorAction(kernel, "Geant4PrimaryHandler/PrimaryHandler")
    #gen.OutputLevel = 0  # generator_output_level
    gen.enableUI()
    kernel.generatorAction().adopt(gen)

    part = DDG4.GeneratorAction(kernel, "Geant4ParticleHandler/ParticleHandler")
    # part.SaveProcesses = ['conv','Decay']
    #part.SaveProcesses = ['Decay']
    part.MinimalKineticEnergy = 100000 * GeV
    #part.OutputLevel = 5  # generator_output_level
    part.enableUI()
    kernel.generatorAction().adopt(part)

    #gen = DDG4.GeneratorAction(kernel, "Geant4InteractionVertexSmear/SmearPi+")
    #gen.Mask = 1
    #gen.Offset = (20 * mm, 10 * mm, 10 * mm, 0 * ns)
    #gen.Sigma = (4 * mm, 1 * mm, 1 * mm, 0 * ns)
    #kernel.generatorAction().adopt(gen)

    #gen = DDG4.GeneratorAction(kernel, "Geant4IsotropeGenerator/IsotropE-")
    #gen.Mask = 2
    #gen.Particle = 'e-'
    #gen.Energy = 25 * GeV
    #gen.Multiplicity = 3
    #gen.Distribution = 'uniform'
    #kernel.generatorAction().adopt(gen)
    #gen = DDG4.GeneratorAction(kernel, "Geant4InteractionVertexSmear/SmearE-")
    #gen.Mask = 2
    #gen.Offset = (-20 * mm, -10 * mm, -10 * mm, 0 * ns)
    #gen.Sigma = (12 * mm, 8 * mm, 8 * mm, 0 * ns)
    #kernel.generatorAction().adopt(gen)

    #gen = DDG4.GeneratorAction(kernel, "Geant4InteractionMerger/InteractionMerger")
    #gen.OutputLevel = 4  # generator_output_level
    #gen.enableUI()
    #kernel.generatorAction().adopt(gen)

    #gen = DDG4.GeneratorAction(kernel, "Geant4PrimaryHandler/PrimaryHandler")
    #gen.OutputLevel = 4  # generator_output_level
    #gen.enableUI()
    #kernel.generatorAction().adopt(gen)

    #part = DDG4.GeneratorAction(kernel, "Geant4ParticleHandler/ParticleHandler")
    #kernel.generatorAction().adopt(part)
    ## part.SaveProcesses = ['conv','Decay']
    #part.SaveProcesses = ['Decay']
    #part.MinimalKineticEnergy = 100 * MeV
    #part.OutputLevel = 5  # generator_output_level
    #part.enableUI()

    #user = DDG4.Action(kernel, "Geant4TCUserParticleHandler/UserParticleHandler")
    #user.TrackingVolume_Zmax = DDG4.EcalEndcap_zmin
    #user.TrackingVolume_Rmax = DDG4.EcalBarrel_rmin
    #user.enableUI()
    #part.adopt(user)

    #f1 = DDG4.Filter(kernel, 'GeantinoRejectFilter/GeantinoRejector')

    f2 = DDG4.Filter(kernel, 'ParticleRejectFilter/OpticalPhotonRejector')
    f2.particle = 'opticalphoton'

    f3 = DDG4.Filter(kernel, 'ParticleSelectFilter/OpticalPhotonSelector')
    f3.particle = 'opticalphoton'

    #f4 = DDG4.Filter(kernel, 'EnergyDepositMinimumCut')
    #f4.Cut = 10 * MeV

    #f4.enableUI()
    #kernel.registerGlobalFilter(f1)
    kernel.registerGlobalFilter(f2)
    kernel.registerGlobalFilter(f3)
    #kernel.registerGlobalFilter(f4)

    seq, act = geant4.setupDetector('LightGasCherenkov','PhotoMultiplierSDAction')
    act.adopt(f3)
    seq, act = geant4.setupDetector('HeavyGasCherenkov','PhotoMultiplierSDAction')
    act.adopt(f3)

    seq, act = geant4.setupTracker('GEMTracker_SIDIS')
    seq, act = geant4.setupCalorimeter('FAECPreShower')
    seq, act = geant4.setupCalorimeter('FAECShower')
    seq, act = geant4.setupCalorimeter('LAECPreShower')
    seq, act = geant4.setupCalorimeter('LAECShower')

    phys = geant4.setupPhysics('QGSP_BERT')
    geant4.addPhysics(str('Geant4PhysicsList/Myphysics'))

    ph = DDG4.PhysicsList(kernel, 'Geant4OpticalPhotonPhysics/OpticalPhotonPhys')
    ph.VerboseLevel = 0
    ph.enableUI()
    phys.adopt(ph)

    ph = DDG4.PhysicsList(kernel, 'Geant4CerenkovPhysics/CerenkovPhys')
    ph.MaxNumPhotonsPerStep = 10
    ph.MaxBetaChangePerStep = 10.0
    ph.TrackSecondariesFirst = True
    ph.VerboseLevel = 0
    ph.enableUI()
    phys.adopt(ph)

    ## Add special particle types from specialized physics constructor
    #part = geant4.addPhysics('Geant4ExtraParticles/ExtraParticles')
    #part.pdgfile = 'checkout/DDG4/examples/particle.tbl'

    # Add global range cut
    rg = geant4.addPhysics('Geant4DefaultRangeCut/GlobalRangeCut')
    rg.RangeCut = 0.7 * mm

    if args.verbose > 1 :
        phys.dump()

    #ui_action = dd4hep.sim.createAction(kernel, "Geant4UIManager/UI")
    #ui_action.HaveVIS = True
    #ui_action.HaveUI = True
    #ui_action.SessionType = qt
    #ui_action.SetupUI = macro
    #kernel.registerGlobalAction(ui_action)

    kernel.configure()
    kernel.initialize()

    # DDG4.setPrintLevel(Output.DEBUG)
    kernel.run()
    kernel.terminate()

if __name__ == "__main__":
    run()
