#include "LYSimTrajectoryPoint.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4StepStatus.hh"

#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4AttDefStore.hh"

#include "G4UnitsTable.hh"

//#define G4ATTDEBUG
#ifdef G4ATTDEBUG
#include "G4AttCheck.hh"
#endif

G4Allocator<LYSimTrajectoryPoint> LYSimTrajPointAllocator;

LYSimTrajectoryPoint::LYSimTrajectoryPoint()
      : fTime(0.), fMomentum(0.,0.,0.),
        fStepStatus(fUndefined), fVolumeName("") { }

LYSimTrajectoryPoint::LYSimTrajectoryPoint(const G4Step* aStep)
    : G4TrajectoryPoint(aStep->GetPostStepPoint()->GetPosition())
{
      fTime = aStep->GetPostStepPoint()->GetGlobalTime();
      fMomentum = aStep->GetPostStepPoint()->GetMomentum();
      fStepStatus = aStep->GetPostStepPoint()->GetStepStatus();
      if (aStep->GetPostStepPoint()->GetPhysicalVolume())
      {
         fVolumeName = aStep->GetPostStepPoint()->
                              GetPhysicalVolume()->GetName();
      } else {
         fVolumeName = " ";
      }
}

LYSimTrajectoryPoint::LYSimTrajectoryPoint(const G4Track* aTrack)
    : G4TrajectoryPoint(aTrack->GetPosition())
{
      fTime = aTrack->GetGlobalTime();
      fMomentum = aTrack->GetMomentum();
      fStepStatus = fUndefined;
      fVolumeName = aTrack->GetVolume()->GetName();
}

LYSimTrajectoryPoint::LYSimTrajectoryPoint(const LYSimTrajectoryPoint &right)
    : G4TrajectoryPoint(right)
{
      fTime = right.fTime;
      fMomentum = right.fMomentum;
      fStepStatus = right.fStepStatus;
      fVolumeName = right.fVolumeName;
}

LYSimTrajectoryPoint::~LYSimTrajectoryPoint() { }

const std::map<G4String,G4AttDef>* LYSimTrajectoryPoint::GetAttDefs() const
{
  G4bool isNew;
  std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("TrajectoryPoint",isNew);
  if (isNew) {
    G4String Pos("Pos");
    (*store)[Pos] =
      G4AttDef(Pos, "Position", "Physics","G4BestUnit","G4ThreeVector");

    G4String Time("Time");
    (*store)[Time] =
      G4AttDef(Time, "Time", "Physics","G4BestUnit","G4double");

    G4String Momentum("Momentum");
    (*store)[Momentum] =
      G4AttDef(Momentum, "Momentum", "Physics","G4BestUnit","G4ThreeVector");

    G4String StepStatus("StepStatus");
    (*store)[StepStatus] =
      G4AttDef(StepStatus, "StepStatus", "Physics", "", "G4StepStatus");

    G4String VolumeName("VolumeName");
    (*store)[VolumeName] =
      G4AttDef(VolumeName, "VolumeName", "Physics", "", "G4String");

  }
  return store;
}

std::vector<G4AttValue>* LYSimTrajectoryPoint::CreateAttValues() const
{
  std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

  values->push_back(G4AttValue("Time",G4BestUnit(fTime,"Time"),""));

  values->push_back(G4AttValue("Momentum",G4BestUnit(fMomentum,"Momentum"),""));

  values->push_back(G4AttValue("StepStatus",fStepStatus,""));

  values->push_back(G4AttValue("VolumeName",fVolumeName,""));

#ifdef G4ATTDEBUG
  G4cout << G4AttCheck(values,GetAttDefs());
#endif

  return values;
}
