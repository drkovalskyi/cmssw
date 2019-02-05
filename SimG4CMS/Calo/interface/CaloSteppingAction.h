#ifndef SimG4CMS_CaloSteppingAction_H
#define SimG4CMS_CaloSteppingAction_H

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"

#include "SimG4Core/Notification/interface/Observer.h"
#include "SimG4Core/Notification/interface/BeginOfJob.h"
#include "SimG4Core/Notification/interface/BeginOfRun.h"
#include "SimG4Core/Notification/interface/BeginOfEvent.h"
#include "SimG4Core/Notification/interface/EndOfEvent.h"
#include "SimG4Core/Watcher/interface/SimProducer.h"
#include "SimG4Core/Watcher/interface/SimWatcherFactory.h"

#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/SimHitMaker/interface/CaloSlaveSD.h"

#include "SimG4CMS/Calo/interface/CaloGVHit.h"
#include "SimG4CMS/Calo/interface/HcalNumberingScheme.h"

#include "Geometry/EcalCommonData/interface/EcalBarrelNumberingScheme.h"
#include "Geometry/EcalCommonData/interface/EcalBaseNumber.h"
#include "Geometry/EcalCommonData/interface/EcalEndcapNumberingScheme.h"
#include "Geometry/HcalCommonData/interface/HcalNumberingFromDDD.h"
#include "Geometry/HcalCommonData/interface/HcalDDDSimConstants.h"

#include "G4LogicalVolume.hh"
#include "G4Region.hh"
#include "G4Step.hh"
#include "G4UserSteppingAction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"
#include "G4Track.hh"

#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

class CaloSteppingAction : public SimProducer,
                           public Observer<const BeginOfJob *>, 
                           public Observer<const BeginOfRun *>, 
                           public Observer<const BeginOfEvent *>, 
                           public Observer<const EndOfEvent *>, 
                           public Observer<const G4Step *> {

public:
  CaloSteppingAction(const edm::ParameterSet &p);
  ~CaloSteppingAction() override;

  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  void fillHits(edm::PCaloHitContainer& cc, int type);
  // observer classes
  void update(const BeginOfJob * job)   override;
  void update(const BeginOfRun * run)   override;
  void update(const BeginOfEvent * evt) override;
  void update(const EndOfEvent * evt)   override;
  void update(const G4Step * step)      override;

  void NaNTrap(const G4Step*) const;
  uint32_t getDetIDEB(const G4Step * step) const;
  uint32_t getDetIDEE(const G4Step * step) const;
  uint32_t getDetIDHC(const G4Step * step) const;
  EcalBaseNumber getBaseNumber(const G4Step*) const; 
  void fillHit(const G4Step * step, uint32_t id, double dE, int flag);
  uint16_t getDepth(const G4Step* step, int type) const;
  G4ThreeVector setToLocal(const G4ThreeVector& global, 
			   const G4VTouchable* touch) const;
  double   curve_LY(const G4LogicalVolume* lv, double z);
  double   getBirkL3(const G4Step* aStep);
  double   getBirkHC(const G4Step* aStep);
  void     saveHits(int flag);

  static const int                      nSD_= 3;
  std::unique_ptr<EcalBarrelNumberingScheme> ebNumberingScheme_;
  std::unique_ptr<EcalEndcapNumberingScheme> eeNumberingScheme_;
  std::unique_ptr<HcalNumberingFromDDD>      hcNumbering_;
  std::unique_ptr<HcalNumberingScheme>       hcNumberingScheme_;
  std::unique_ptr<CaloSlaveSD>               slave_[nSD_];

  std::vector<std::string>              nameEBSD_, nameEESD_, nameHCSD_;
  std::vector<std::string>              nameHitC_;
  std::vector<const G4LogicalVolume*>   volEBSD_, volEESD_, volHCSD_;
  std::map<const G4LogicalVolume*,double> xtalMap_;
  int                                   count_;
  double                                slopeLY_, birkC1EC_, birkSlopeEC_;
  double                                birkCutEC_, birkC1HC_, birkC2HC_;
  double                                birkC3HC_;
  std::map<CaloHitID,CaloGVHit>         hitMap_[nSD_];
};

#endif
