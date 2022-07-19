'''
    An example script to digitize/reconstruct/clustering endcap ecal hits
'''
from Gaudi.Configuration import *
import os
import ROOT

from Configurables import ApplicationMgr, EICDataSvc, PodioInput, PodioOutput, GeoSvc
from GaudiKernel.SystemOfUnits import MeV, GeV, mm, cm, mrad

detector_name = str(os.environ.get("JUGGLER_DETECTOR", "solid_sidis"))
detector_path = str(os.environ.get("DETECTOR_PATH", "."))
compact_path = str(os.environ.get("JUGGLER_COMPACT_PATH", "{}.xml".format(os.path.join(detector_path, detector_name))))

# input and output
input_sims = [f.strip() for f in str.split(os.environ["JUGGLER_SIM_FILE"], ",") if f.strip()]
output_rec = str(os.environ["JUGGLER_REC_FILE"])
n_events = int(os.environ["JUGGLER_N_EVENTS"])
print(input_sims)

# geometry service
geo_service = GeoSvc("GeoSvc", detectors=[compact_path], OutputLevel=INFO)
# data service
podioevent = EICDataSvc("EventDataSvc", inputs=input_sims)


# juggler components
from Configurables import Jug__Digi__CalorimeterHitDigi as CalHitDigi
from Configurables import Jug__Reco__CalorimeterHitReco as CalHitReco
from Configurables import Jug__Reco__CalorimeterHitsMerger as CalHitsMerger
from Configurables import Jug__Reco__CalorimeterIslandCluster as IslandCluster
from Configurables import Jug__Reco__ClusterRecoCoG as RecoCoG

# branches needed from simulation root file
sim_coll = [
    "MCParticles",
    "FAEC_ShHits",
]

# input and output
podin = PodioInput("PodioReader", collections=sim_coll)
podout = PodioOutput("out", filename=output_rec)


# FAEC Shower
ce_ecal_daq = dict(
        dynamicRangeADC=5.*GeV,
        capacityADC=32768,
        pedestalMean=400,
        pedestalSigma=3)

ce_ecal_digi = CalHitDigi("ce_ecal_digi",
        inputHitCollection="FAEC_ShHits",
        outputHitCollection="FAEC_ShHits_Digi",
        energyResolutions=[0., 0.02, 0.],
        **ce_ecal_daq)

ce_ecal_reco = CalHitReco("ce_ecal_reco",
        inputHitCollection=ce_ecal_digi.outputHitCollection,
        outputHitCollection="FAEC_ShHits_Reco",
        thresholdFactor=4,          # 4 sigma cut on pedestal sigma
        readoutClass="FAEC_ShHits",
        sectorField="system",
        **ce_ecal_daq)

# merge hits in different layer (projection to local x-y plane)
ce_ecal_merger = CalHitsMerger("ci_ecal_merger",
        # OutputLevel=DEBUG,
        inputHitCollection=ce_ecal_reco.outputHitCollection,
        outputHitCollection="FAEC_ShHits_RecoMod",
        fields=["layer", "slice"],
        fieldRefNumbers=[0, 0],
        readoutClass="FAEC_ShHits")

ce_ecal_cl = IslandCluster("ce_ecal_cl",
        # OutputLevel=DEBUG,
        inputHitCollection=ce_ecal_merger.outputHitCollection,
        outputProtoClusterCollection="FAEC_ShHits_ProtoClusters",
        splitCluster=False,
        minClusterHitEdep=1.0*MeV,  # discard low energy hits
        minClusterCenterEdep=30*MeV,
        sectorDist=5.0*cm,
        dimScaledLocalDistXY=[1.8, 1.8])          # hybrid calorimeter with different dimensions, using a scaled dist

ce_ecal_clreco = RecoCoG("ce_ecal_clreco",
        inputProtoClusterCollection=ce_ecal_cl.outputProtoClusterCollection,
        outputClusterCollection="FAEC_ShHits_Clusters",
        samplingFraction=0.998,      # this accounts for a small fraction of leakage
        logWeightBase=4.6)

podout.outputCommands = ['keep *']

ApplicationMgr(
    TopAlg = [podin,
              ce_ecal_digi, ce_ecal_reco, ce_ecal_merger,
              ce_ecal_cl, ce_ecal_clreco,
              podout],
    EvtSel = 'NONE',
    EvtMax = n_events,
    ExtSvc = [podioevent],
    OutputLevel=DEBUG
)
