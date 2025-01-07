"""
A steering file that enables optical photon counting for LGC
It sets up physics lists and particle filters for LGC, as well as setting a particle gun for SoLID configuration
   @author  C.Peng

[simulation]:
    ddsim --steeringFile benchmarks/ECal/steering.py --compactFile solid.xml
[visualization]:
    ddsim --steeringFile benchmarks/ECal/steering.py --compactFile solid.xml --runType qt --macro macro/vis.mac
"""
from __future__ import absolute_import, unicode_literals
import logging
import sys
import os
from DDSim.DD4hepSimulation import DD4hepSimulation

SIM = DD4hepSimulation()

# Allow energy depositions to 0 energy in trackers (which include optical detectors)
SIM.filter.tracker = "edep0"

# Particle gun settings: electrons with fixed energy and theta, varying phi
SIM.numberOfEvents = 100
SIM.enableGun = True
SIM.gun.position = (0., 0., "-350*cm")
SIM.gun.energy = "5*GeV"
SIM.gun.particle = "e-"
SIM.gun.thetaMin = "10.0*deg"
SIM.gun.thetaMax = "10.0*deg"
SIM.gun.phiMin = "0.*deg"
SIM.gun.phiMax = "360.*deg"
SIM.gun.distribution = "cos(theta)"

# Output file (assuming CWD)
SIM.outputFile = "ecal_hits.edm4hep.root"
# SIM.outputConfig.forceDD4HEP = True

