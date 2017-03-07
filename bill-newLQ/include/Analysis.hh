#ifndef ANALYSIS_HH_
#define ANALYSIS_HH_

#include "G4Event.hh"
#include "G4Run.hh"
//ROOT
//#include "TH1D.h"

#define NUMRUNS 10

class LYSimDetectorConstruction;
class AnalysisMessenger;

/*!
 * \brief Analysis class
 * This class contains the code to collect information from
 * the different UserActions.
 * The class is designed as a singleton.
 * To access it you need to use:
 * Analysis* analysis = Analysis::GetInstance()
 */
class Analysis {
public:
	//! Singleton pattern
	static Analysis* GetInstance() {
		if ( Analysis::singleton == NULL ) Analysis::singleton = new Analysis();
		return Analysis::singleton;
	}
	//! destructor
	virtual ~Analysis();
	//Set pointer to Detector Construction
	void SetDetector(LYSimDetectorConstruction* det){
		DetectorConstruction = det;
	};
	//Set output file name
	void SetOutputFile(std::string filename) {fOutputFileName = filename;};
	//Set ROOT file name
	void SetROOTFile(std::string filename) {fROOTFileName = filename;};
	//! Should be called at the beginning of an event
	void PrepareNewEvent(const G4Event* anEvent);
	//! Should be called at the end of an event
	void EndOfEvent(const G4Event* anEvent);
	//! Should be called at the beginning of a run
	void PrepareNewRun(const G4Run* aRun);
	//! Should be called at the end of a run
	void EndOfRun(const G4Run* aRun);
	
	// Called at the end of experiment
	void EndOfExperiment();
	
	// Increase the generated photon count
	void AddPhotonCount( G4int num ) { PhotonCount += num; }
	// Increase the detected photon hit count
	void AddHitCount( G4int num ) { HitCount += num; }
	
	// Set tile absorption length to record
	void SetTileAbsLength ( G4double length) { tileAbsLength = length;}
	
private:
	//! Private construtor: part of singleton pattern
	Analysis();
	//! Singleton static instance
	static Analysis* singleton;
	//Pointer to Analysis messenger
	AnalysisMessenger* fMessenger;
	//Pointer to DetectorConstruction class for access to detector properties
    LYSimDetectorConstruction* DetectorConstruction;   
    //Output file name
    std::string fOutputFileName;
    //ROOT file name
    std::string fROOTFileName;
   	//Total number of photons generated per run
	G4int PhotonCount;
	//Total number of photon hits detected per run
	G4int HitCount;
	//Array of Scintillator thickness for each run
	G4double ThicknessArray[NUMRUNS];
	//Array of (light collection) efficiency for each run
	G4double EfficiencyArray[NUMRUNS];
    //Histograms
    //TH1D* EnergyHist;
	//TH1D* PhotonHitsHist;
	
	G4double tileAbsLength;
};

#endif /* ANALYSIS_HH_ */
