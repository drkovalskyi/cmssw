import FWCore.ParameterSet.Config as cms

simKBmtfStubs = cms.EDProducer("L1TMuonBarrelKalmanStubProducer",
    verbose = cms.int32(0),
    srcPhi = cms.InputTag("simTwinMuxDigis"),
    srcTheta = cms.InputTag("simDtTriggerPrimitiveDigis"),
    minPhiQuality = cms.int32(0),
    minThetaQuality = cms.int32(0),
    minBX = cms.int32(-3),                           
    maxBX = cms.int32(3),                           
    etaLUT_minus_2_1 = cms.vint32(-106,-102,-98,-94,-89,-84,-79),
    etaLUT_minus_2_2 = cms.vint32(-94,-90,-87,-82,-78,-73,-69),
    etaLUT_minus_2_3 = cms.vint32(-81,-78,-74,-71,-67,-63,-58),
    etaLUT_minus_1_1 = cms.vint32(-74,-68,-61,-55,-48,-41,-33),
    etaLUT_minus_1_2 = cms.vint32(-63,-58,-53,-47,-41,-35,-28),
    etaLUT_minus_1_3 = cms.vint32(-54,-49,-45,-40,-35,-29,-24),
    etaLUT_0_1 = cms.vint32(-25,-17,-9,0,9,17,25),
    etaLUT_0_2 = cms.vint32(-21,-15,-7,0,7,15,21),
    etaLUT_0_3 = cms.vint32(-18,-12,-6,0,6,12,18),
    etaLUT_plus_1_1 = cms.vint32(74,68,61,55,48,41,33),
    etaLUT_plus_1_2 = cms.vint32(63,58,53,47,41,35,28),
    etaLUT_plus_1_3 = cms.vint32(54,49,45,40,35,29,24),
    etaLUT_plus_2_1 = cms.vint32(106,102,98,94,89,84,79),
    etaLUT_plus_2_2 = cms.vint32(94,90,87,82,78,73,69),
    etaLUT_plus_2_3 = cms.vint32(81,78,74,71,67,63,58),
    etaCoarseLUT_minus_2 = cms.vint32(-89,-78,-67,-62),
    etaCoarseLUT_minus_1 = cms.vint32(-48,-41,-35,-34),
    etaCoarseLUT_0 = cms.vint32(0,0,0,0),
    etaCoarseLUT_plus_1= cms.vint32(48,41,35,34),
    etaCoarseLUT_plus_2= cms.vint32(89,78,67,62)
                           
)