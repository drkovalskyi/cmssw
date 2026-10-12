// -*- C++ -*-
//
// Package:    PhysicsTools/PatFromScouting
// Class:      PatFromScoutingMuonProducer
//
/**\class PatFromScoutingMuonProducer PatFromScoutingMuonProducer.cc PhysicsTools/PatFromScouting/plugins/PatFromScoutingMuonProducer.cc

 Description: Converts Run3ScoutingMuon to pat::Muon

 Implementation:
     Uses the pat::Muon(const Run3ScoutingMuon&) constructor which sets:
     - Kinematics (pt, eta, phi, mass, charge)
     - Embedded track with vertex
     - Isolation (trackIso, ecalIso, hcalIso -> isolationR03)
     Muons whose track matches a scouting PF muon candidate get the PFMuon type bit,
     which makes isPFMuon() and isLooseMuon() usable on scouting data.
*/
//
// Original Author:  Dmytro Kovalskyi
//         Created:  Thu, 05 Dec 2024 15:27:09 GMT
//
//

#include <cmath>
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"

class PatFromScoutingMuonProducer : public edm::stream::EDProducer<> {
public:
  explicit PatFromScoutingMuonProducer(const edm::ParameterSet&);
  ~PatFromScoutingMuonProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  const edm::EDGetTokenT<Run3ScoutingMuonCollection> muonToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingParticle>> pfCandidateToken_;
  const double maxPFMatchMetric2_;
};

PatFromScoutingMuonProducer::PatFromScoutingMuonProducer(const edm::ParameterSet& iConfig)
    : muonToken_(consumes<Run3ScoutingMuonCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      maxPFMatchMetric2_(std::pow(iConfig.getParameter<double>("maxPFMatchMetric"), 2)) {
  const auto& pfCandidateTag = iConfig.getParameter<edm::InputTag>("pfCandidates");
  if (!pfCandidateTag.label().empty())
    pfCandidateToken_ = consumes<std::vector<Run3ScoutingParticle>>(pfCandidateTag);
  produces<pat::MuonCollection>();
}

void PatFromScoutingMuonProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  auto patMuons = std::make_unique<pat::MuonCollection>();

  const auto& scoutingMuons = iEvent.get(muonToken_);

  // Track parameters of scouting PF muon candidates, used to flag PF muons
  struct PFMuonTrack {
    float pt, eta, phi;
    int charge;
    bool used;
  };
  std::vector<PFMuonTrack> pfMuonTracks;
  if (!pfCandidateToken_.isUninitialized()) {
    auto pfCandidates = iEvent.getHandle(pfCandidateToken_);
    if (pfCandidates.isValid()) {
      for (const auto& particle : *pfCandidates) {
        if (std::abs(particle.pdgId()) != 13)
          continue;
        bool relativeTrackVars = particle.relative_trk_vars();
        pfMuonTracks.push_back({relativeTrackVars ? particle.trk_pt() + particle.pt() : particle.trk_pt(),
                                relativeTrackVars ? particle.trk_eta() + particle.eta() : particle.trk_eta(),
                                relativeTrackVars ? particle.trk_phi() + particle.phi() : particle.trk_phi(),
                                particle.pdgId() > 0 ? -1 : 1,
                                false});
      }
    }
  }

  for (const auto& sMuon : scoutingMuons) {
    // Constructor now handles kinematics, track embedding, and isolation
    pat::Muon muon(sMuon);

    // Match the muon track to a PF muon track: dR^2 + (dpt/pt)^2 < maxPFMatchMetric^2
    PFMuonTrack* bestMatch = nullptr;
    float bestMetric2 = maxPFMatchMetric2_;
    for (auto& pfMuon : pfMuonTracks) {
      if (pfMuon.used || pfMuon.charge != sMuon.charge() || sMuon.trk_pt() <= 0)
        continue;
      float dEta = pfMuon.eta - sMuon.trk_eta();
      float dPhi = reco::deltaPhi(pfMuon.phi, sMuon.trk_phi());
      float dPtRel = (pfMuon.pt - sMuon.trk_pt()) / sMuon.trk_pt();
      float metric2 = dEta * dEta + dPhi * dPhi + dPtRel * dPtRel;
      if (metric2 < bestMetric2) {
        bestMetric2 = metric2;
        bestMatch = &pfMuon;
      }
    }
    if (bestMatch) {
      bestMatch->used = true;
      muon.setType(muon.type() | reco::Muon::PFMuon);
    }

    patMuons->push_back(muon);
  }

  iEvent.put(std::move(patMuons));
}

void PatFromScoutingMuonProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("hltScoutingMuonPacker"));
  desc.add<edm::InputTag>("pfCandidates", edm::InputTag("hltScoutingPFPacker"))
      ->setComment("Scouting PF candidates used to set the PFMuon type bit; empty label disables the matching");
  desc.add<double>("maxPFMatchMetric", 0.01)
      ->setComment("Maximum sqrt(dR^2 + (dpt/pt)^2) between the muon track and a PF muon track");
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(PatFromScoutingMuonProducer);
