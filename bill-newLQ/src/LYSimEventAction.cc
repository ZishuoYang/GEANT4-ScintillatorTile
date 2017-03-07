// LYSimEventAction.cc


#include "LYSimEventAction.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4DigiManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "Analysis.hh"

LYSimEventAction::LYSimEventAction()
{
}


void LYSimEventAction::BeginOfEventAction(const G4Event* anEvent )
{
	if ( anEvent->GetEventID() % 100 == 0 )
	{
		G4cout<<"Starting Event: "<<anEvent->GetEventID()<<G4endl;
	}
	Analysis::GetInstance()->PrepareNewEvent(anEvent);
	//Retrieve the ID for the hit collection
	// if ( hitsCollID == -1 )
	// {
		// G4SDManager * SDman = G4SDManager::GetSDMpointer();
		// hitsCollID = SDman->GetCollectionID(hitsCollName);
	// }
}

void LYSimEventAction::EndOfEventAction(const G4Event* anEvent)
{
	Analysis::GetInstance()->EndOfEvent(anEvent);
}

