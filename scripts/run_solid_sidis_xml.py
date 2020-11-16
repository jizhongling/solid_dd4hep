#!/usr/bin/env python
from __future__ import absolute_import, unicode_literals
import os
import time
import logging
import DDG4
from DDG4 import OutputLevel as Output
from g4units import keV, GeV, mm, ns, MeV

def run():
  #os.environ['G4UI_USE_TCSH'] = "1"
  kernel = DDG4.Kernel()
  description = kernel.detectorDescription()

  kernel.loadGeometry(str("file:" + "solid_sidis.xml"))
  kernel.loadXML(str("file:solid/sim/field.xml"))
  kernel.loadXML(str("file:solid/sim/sequences.xml"))
  kernel.loadXML(str("file:solid/sim/physics.xml"))

  kernel.configure()
  kernel.initialize()

  kernel.run()
  kernel.terminate()

if __name__ == "__main__":
  run()
