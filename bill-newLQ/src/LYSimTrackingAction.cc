#include "globals.hh"
#include "G4RunManager.hh"

#include "LYSimTrajectory.hh"

//#include "LYSimUserTrackInformation.hh"

#include "G4Track.hh"
#include "G4ParticleTypes.hh"
#include "G4TrackingManager.hh"

#include "LYSimTrackingAction.hh"

void LYSimTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  //Let this be up to the user via vis.mac
  //  fpTrackingManager->SetStoreTrajectory(true);

  //Use custom trajectory class
  fpTrackingManager->SetTrajectory(new LYSimTrajectory(aTrack));

/*  LYSimUserTrackInformation* trackInformation = new LYSimUserTrackInformation();

  if (aTrack->GetMomentumDirection().z()>0.0) {
     trackInformation->AddStatusFlag(right);
  } else {
     trackInformation->AddStatusFlag(left);
  }

  G4String PVName = aTrack->GetVolume()->GetName();

  if (PVName == "LYSimFiber" || PVName == "Clad1" || PVName == "Clad2")
     trackInformation->AddStatusFlag(InsideOfFiber);

  fpTrackingManager->SetUserTrackInformation(trackInformation);
*/
}

void LYSimTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  LYSimTrajectory* trajectory =
                      (LYSimTrajectory*)fpTrackingManager->GimmeTrajectory();

  if (aTrack->GetDefinition()==G4OpticalPhoton::OpticalPhotonDefinition())
  {
     if (aTrack->GetParentID()==0) trajectory->SetDrawTrajectory(true);
     else {
        const G4VProcess* creator = aTrack->GetCreatorProcess();
        if (creator && creator->GetProcessName()=="OpWLS")
        {
           trajectory->WLS();
           trajectory->SetDrawTrajectory(true);
        }
        if (creator && creator->GetProcessName()=="Scintillation")
        {
           trajectory->SetDrawTrajectory(false);
        }
     }
  }
  else //draw all other trajectories
    trajectory->SetDrawTrajectory(true);
}
