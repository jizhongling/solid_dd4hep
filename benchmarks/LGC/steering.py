"""
A steering file that enables optical photon counting for LGC
It sets up physics lists and particle filters for LGC, as well as setting a particle gun for SoLID configuration
   @author  C.Peng

[simulation]:
    ddsim --steeringFile benchmarks/LGC/steering.py --compactFile solid.xml
[visualization]:
    ddsim --steeringFile benchmarks/LGC/steering.py --compactFile solid.xml --runType qt --macro macro/vis.mac
"""
from __future__ import absolute_import, unicode_literals
import logging
import sys
import os

from DDSim.DD4hepSimulation import DD4hepSimulation


SIM = DD4hepSimulation()

# Ensure that Cerenkov and optical physics are always loaded
def setupCerenkov(kernel):
    from DDG4 import PhysicsList

    seq = kernel.physicsList()
    cerenkov = PhysicsList(kernel, "Geant4CerenkovPhysics/CerenkovPhys")
    cerenkov.MaxNumPhotonsPerStep = 10
    cerenkov.MaxBetaChangePerStep = 10.0
    cerenkov.TrackSecondariesFirst = False
    cerenkov.VerboseLevel = 0
    cerenkov.enableUI()
    seq.adopt(cerenkov)
    ph = PhysicsList(kernel, "Geant4OpticalPhotonPhysics/OpticalGammaPhys")
    ph.addParticleConstructor("G4OpticalPhoton")
    ph.VerboseLevel = 0
    ph.enableUI()
    seq.adopt(ph)
    return None

SIM.physics.setupUserPhysics(setupCerenkov)

# Allow energy depositions to 0 energy in trackers (which include optical detectors)
SIM.filter.tracker = "edep0"

# Some detectors are only sensitive to optical photons
SIM.filter.filters["opticalphotons"] = dict(
    name="ParticleSelectFilter/OpticalPhotonSelector",
    parameter={"particle": "opticalphoton"},
    )
SIM.filter.mapDetFilter["LightGasCherenkov"] = "opticalphotons"

# Use the optical tracker for the PFRICH
SIM.action.mapActions["LightGasCherenkov"] = "Geant4OpticalTrackerAction"

# Disable user tracker particle handler, so hits can be associated to photons
SIM.part.userParticleHandler = ""

# Particle gun settings: electrons with fixed energy and theta, varying phi
SIM.numberOfEvents = 100
SIM.enableGun = True
SIM.gun.position = (0., 0., "-300*cm")
SIM.gun.energy = "5*GeV"
SIM.gun.particle = "e-"
SIM.gun.thetaMin = "8.0*deg"
SIM.gun.thetaMax = "16.0*deg"
SIM.gun.distribution = "cos(theta)"

# Output file (assuming CWD)
SIM.outputFile = "lgc_hits.edm4hep.root"
# SIM.outputConfig.forceDD4HEP = True

