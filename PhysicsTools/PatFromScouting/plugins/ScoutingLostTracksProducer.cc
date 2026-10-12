/**
 * ScoutingLostTracksProducer
 *
 * Builds the MiniAOD lostTracks collection from scouting tracks: every reco::Track
 * not matched to a packedPFCandidate becomes a pat::PackedCandidate with track details,
 * mirroring PATLostTracks. Also produces an Association from the reco::Tracks to the
 * lost track candidates, as needed by V0 slimming.
 *
 * Scouting tracks come from hltPFMuonMerging, which contains only high-purity tracks,
 * so all lost tracks are flagged as high purity.
 */

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

class ScoutingLostTracksProducer : public edm::global::EDProducer<> {
public:
  explicit ScoutingLostTracksProducer(const edm::ParameterSet&);
  ~ScoutingLostTracksProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  const edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  const edm::EDGetTokenT<edm::ValueMap<int>> trackVertexIndexToken_;
  const edm::EDGetTokenT<edm::Association<pat::PackedCandidateCollection>> trackToPFCandidateToken_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  const int covarianceVersion_;
  const int covarianceSchema_;
};

ScoutingLostTracksProducer::ScoutingLostTracksProducer(const edm::ParameterSet& iConfig)
    : trackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
      trackVertexIndexToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("trackVertexIndex"))),
      trackToPFCandidateToken_(consumes<edm::Association<pat::PackedCandidateCollection>>(
          iConfig.getParameter<edm::InputTag>("packedPFCandidates"))),
      vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
      covarianceVersion_(iConfig.getParameter<int>("covarianceVersion")),
      covarianceSchema_(iConfig.getParameter<int>("covarianceSchema")) {
  produces<pat::PackedCandidateCollection>();
  produces<edm::Association<pat::PackedCandidateCollection>>();
}

void ScoutingLostTracksProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  constexpr float pionMass = 0.13957018;

  auto output = std::make_unique<pat::PackedCandidateCollection>();

  auto tracksHandle = iEvent.getHandle(trackToken_);
  const auto& tracks = *tracksHandle;
  const auto& trackVertexIndex = iEvent.get(trackVertexIndexToken_);
  const auto& trackToPFCandidate = iEvent.get(trackToPFCandidateToken_);

  auto verticesHandle = iEvent.getHandle(vertexToken_);
  const auto& vertices = *verticesHandle;
  reco::VertexRefProd vertexRefProd(verticesHandle);

  std::vector<int> trackToCandidate(tracks.size(), -1);

  for (size_t iTk = 0; iTk < tracks.size(); ++iTk) {
    reco::TrackRef trackRef(tracksHandle, iTk);
    if (trackToPFCandidate[trackRef].isNonnull())
      continue;

    const auto& track = tracks[iTk];

    int vtxIdx = trackVertexIndex[trackRef];
    reco::VertexRef::key_type pvKey = 0;
    if (vtxIdx >= 0 && static_cast<size_t>(vtxIdx) < vertices.size()) {
      pvKey = static_cast<reco::VertexRef::key_type>(vtxIdx);
    }

    reco::Candidate::PolarLorentzVector p4(track.pt(), track.eta(), track.phi(), pionMass);
    pat::PackedCandidate cand(
        p4, track.vertex(), track.pt(), track.eta(), track.phi(), 211 * track.charge(), vertexRefProd, pvKey);

    // Scouting tracks carry no information on missing inner hits
    cand.setLostInnerHits(pat::PackedCandidate::noLostInnerHits);
    cand.setTrackHighPurity(true);
    cand.setTrackProperties(track, covarianceSchema_, covarianceVersion_);

    // Same PV association convention as for packedPFCandidates
    if (vtxIdx == 0) {
      cand.setAssociationQuality(pat::PackedCandidate::UsedInFitTight);
    } else if (vtxIdx > 0) {
      cand.setAssociationQuality(pat::PackedCandidate::CompatibilityDz);
    } else {
      cand.setAssociationQuality(pat::PackedCandidate::NotReconstructedPrimary);
    }

    trackToCandidate[iTk] = output->size();
    output->push_back(cand);
  }

  auto outputHandle = iEvent.put(std::move(output));

  auto trackToCandidateMap = std::make_unique<edm::Association<pat::PackedCandidateCollection>>(outputHandle);
  edm::Association<pat::PackedCandidateCollection>::Filler filler(*trackToCandidateMap);
  filler.insert(tracksHandle, trackToCandidate.begin(), trackToCandidate.end());
  filler.fill();
  iEvent.put(std::move(trackToCandidateMap));
}

void ScoutingLostTracksProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracks", edm::InputTag("scoutingTracks"))->setComment("Input reco::Track collection");
  desc.add<edm::InputTag>("trackVertexIndex", edm::InputTag("scoutingTracks", "vertexIndex"))
      ->setComment("Index of the primary vertex associated to each track");
  desc.add<edm::InputTag>("packedPFCandidates", edm::InputTag("packedPFCandidates"))
      ->setComment("Association from tracks to packedPFCandidates; matched tracks are not lost tracks");
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"))
      ->setComment("Input vertex collection for vertex references");
  desc.add<int>("covarianceVersion", 1)->setComment("Covariance parameterization version (0=Phase0, 1=Phase1)");
  desc.add<int>("covarianceSchema", 520)->setComment("Covariance packing schema");
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(ScoutingLostTracksProducer);
