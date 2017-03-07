#include "Analysis.hh"

#include "LYSimRunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"


#include "LYSimDetectorConstruction.hh"

LYSimRunAction::LYSimRunAction(LYSimDetectorConstruction* ipDetectorConstruction)
{
//  runMessenger = new LYSimRunActionMessenger(this);
  pDetectorConstruction = ipDetectorConstruction;
}

LYSimRunAction::~LYSimRunAction()
{
//  delete runMessenger;
}

void LYSimRunAction::BeginOfRunAction(const G4Run* aRun)
{
  Analysis::GetInstance()->PrepareNewRun(aRun);
}

void LYSimRunAction::EndOfRunAction(const G4Run* aRun)
{
  Analysis::GetInstance()->EndOfRun(aRun);
}
