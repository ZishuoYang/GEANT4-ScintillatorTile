#include "Analysis.hh"
#include "AnalysisMessenger.hh"
#include "G4UnitsTable.hh"
#include "G4SDManager.hh"
#include "G4PhysicalConstants.hh"
#include "LYSimDetectorConstruction.hh"
#include "G4ios.hh"
#include <iostream>
#include <fstream>
using namespace std;

//ROOT Stuff
//#include "TProfile.h"
//#include "TFile.h"

Analysis* Analysis::singleton = 0;

//Constructor
Analysis::Analysis()
{
	fMessenger = new AnalysisMessenger(this);
}

Analysis::~Analysis()
{
	if(fMessenger) delete fMessenger;
}


void Analysis::PrepareNewEvent(const G4Event* /*anEvent*/)
{
}

void Analysis::EndOfEvent(const G4Event* anEvent)
{
	G4String hitCollName = "PMTHitsCollection";
	G4SDManager* SDman = G4SDManager::GetSDMpointer();	
	static G4int hitCollID = -1;
	if ( hitCollID < 0 )
	hitCollID = SDman->GetCollectionID( hitCollName );
	
	G4HCofThisEvent* hitsCollections = 0;
	hitsCollections = anEvent->GetHCofThisEvent();
	
	LYSimPMTHitsCollection* hits = 0;
	
	if ( hitsCollections )
	{
		hits = static_cast<LYSimPMTHitsCollection*> ( hitsCollections->GetHC(hitCollID) );
	}
	else
	{
		G4cerr << "hitsCollection not found" << G4endl;
		return;
	}
	LYSimPMTHit* hit = NULL;
	G4double EventEnergy = 0;
	G4int EventPhotonCount = 0;
	PhotonCount++;
	if(hits->entries()==1) //if there is an entry in hits collection
	{
		hit = (*hits)[0];
		EventEnergy = hit->GetEnergy();
		EventPhotonCount = hit->GetPhotonCount();
		HitCount++;
	}
	//G4cout	<< "Energy deposited in PMT this event is " << G4BestUnit(EventEnergy,"Energy") << G4endl
	//		<< "# of photon hits in PMT this event is " << EventPhotonCount << G4endl;
	//EnergyHist->Fill(EventEnergy/eV);
	//PhotonHitsHist->Fill(EventPhotonCount);
	
}


void Analysis::PrepareNewRun(const G4Run* /*aRun*/ )
{
	//Reset variables relative to the run
	PhotonCount = 0;
	HitCount = 0; //
}


void Analysis::EndOfRun(const G4Run* aRun)
{
	
	//Write to file
	//~ TFile* rootfile = TFile::Open(fROOTFileName.c_str(),"recreate");
	//~ EnergyHist->Write();
	//~ PhotonHitsHist->Write();
	//~ rootfile->Close();
	//~ 
	std::ofstream outputfile;
	outputfile.open(fOutputFileName.c_str(), ofstream::out | ofstream::app);
	G4cout	<< "Efficiency in this run is " << (G4double)HitCount/(G4double)PhotonCount << G4endl;
	if (outputfile.is_open())
	{
		outputfile << "#Attenuation Length (cm) \t Efficiency" << G4endl;
		outputfile << tileAbsLength << "\t"<< (G4double)HitCount/(G4double)PhotonCount << G4endl;
	}
	else
	{
		G4cout << "Output file not open" << G4endl;
		G4cout << "Beam spot position (mm) \t Efficiency" << G4endl;
		G4cout << tileAbsLength << "\t"<< (G4double)HitCount/(G4double)PhotonCount << G4endl;
	}
	outputfile.close();
}
