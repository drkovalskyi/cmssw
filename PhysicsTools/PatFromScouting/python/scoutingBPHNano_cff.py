import FWCore.ParameterSet.Config as cms

# BPH NanoAOD content on top of scouting MiniAOD.
# Same producers and tables as the data branch of
# PhysicsTools/NanoAOD/custom_bph_cff.nanoAOD_customizeBPH, scheduled after the
# scouting NanoAOD tables (NANO:@ScoutBPH = @ScoutMini + this sequence).
#
# Inputs from scouting MiniAOD (PAT:@Scout): slimmedMuons, slimmedElectrons,
# packedPFCandidates, lostTracks, offlineSlimmedPrimaryVertices, offlineBeamSpot,
# slimmedKshortVertices and slimmedLambdaVertices.
#
# Scouting limitations:
# - No HLT trigger objects are stored, so BPHMuon_isTriggering and BPHMuon_fired_* are always 0.
# - Scouting tracks are stored only above 1 GeV, so tracks and V0 daughters below 1 GeV
#   have no track details and are not used.

# Standard NanoAOD gets the TransientTrackBuilder via PhysicsTools/NanoAOD/muons_cff,
# which the scouting NanoAOD does not load
from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *

from PhysicsTools.BPHNano.muons_cff import *
from PhysicsTools.BPHNano.MuMu_cff import *
from PhysicsTools.BPHNano.tracks_cff import *
from PhysicsTools.BPHNano.BToKLL_cff import *
from PhysicsTools.BPHNano.DiTrack_cff import *
from PhysicsTools.BPHNano.BToTrkTrkLL_cff import *
from PhysicsTools.BPHNano.V0_cff import *
from PhysicsTools.BPHNano.BToV0LL_cff import *
from PhysicsTools.BPHNano.BToV0TrkLL_cff import *

scoutingBPHNanoSequence = cms.Sequence(
    muonBPHSequence + muonBPHTables
    + MuMuSequence + MuMuTables
    + tracksBPHSequence + tracksBPHTables
    + BToKMuMuSequence + BToKMuMuTables
    + DiTrackSequence + BToTrkTrkMuMuSequence + BToTrkTrkMuMuTables
    + KshortToPiPiSequence + KshortToPiPiTables
    + BToKshortMuMuSequence + BToKshortMuMuTables
    + LambdaToProtonPiSequence + LambdaToProtonPiTables
    + LambdabToLambdaMuMuSequence + LambdabToLambdaMuMuTables
    + BToChargedKstarMuMuSequence + BToChargedKstarsMuMuTable
    + XibToXiMuMuSequence + XibToXiMuMuTable
)
