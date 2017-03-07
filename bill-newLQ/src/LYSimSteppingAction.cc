//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id$
//
/// \file LYSimSteppingAction.cc
/// \brief Implementation of the LYSimSteppingAction class

#include "LYSimSteppingAction.hh"
#include "LYSimDetectorConstruction.hh"
#include "LYSimPMTSD.hh"

#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LYSimSteppingAction* LYSimSteppingAction::fgInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LYSimSteppingAction* LYSimSteppingAction::Instance()
{
	// Static acces function via G4RunManager

	return fgInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LYSimSteppingAction::LYSimSteppingAction()
: G4UserSteppingAction(),
fVolume(0),
PhotonHitCount(0)
{
	fgInstance = this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LYSimSteppingAction::~LYSimSteppingAction()
{
	fgInstance = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LYSimSteppingAction::UserSteppingAction(const G4Step* step)
{
	G4OpBoundaryProcessStatus boundaryStatus=Undefined;
	static G4OpBoundaryProcess* boundary=NULL;

	//find the boundary process only once
	if(!boundary){
		G4ProcessManager* pm
			= step->GetTrack()->GetDefinition()->GetProcessManager();
		G4int nprocesses = pm->GetProcessListLength();
		G4ProcessVector* pv = pm->GetProcessList();
		G4int i;
		for( i=0;i<nprocesses;i++){
			if((*pv)[i]->GetProcessName()=="OpBoundary"){
				boundary = (G4OpBoundaryProcess*)(*pv)[i];
				break;
			}
		}
	}

	if(boundary)
	{
		boundaryStatus=boundary->GetStatus();

		switch(boundaryStatus){
			case Detection: //Note, this assumes that the volume causing detection
						  //is the photocathode because it is the only one with
						  //non-zero efficiency
				{
					//Trigger sensitive detector manually since photon is
					//absorbed but status was Detection
					G4SDManager* SDman = G4SDManager::GetSDMpointer();
					G4String sdName="/LYSimPMT";
					LYSimPMTSD* pmtSD = (LYSimPMTSD*)SDman->FindSensitiveDetector(sdName);
					if(pmtSD)pmtSD->ProcessHits_constStep(step,NULL);
					break;
				}
			default:
				break;
		}
	}

	//kill tracks with length > 5000mm
	G4double tracklength = step->GetTrack()->GetTrackLength();
	if (tracklength > 5000.*mm)
	{
		G4cout << "Track length exceeded limit of 5000 mm" << G4endl;
		step->GetTrack()->SetTrackStatus(fStopAndKill);
	}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LYSimSteppingAction::Reset()
{
	PhotonHitCount = 0;
}

