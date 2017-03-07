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
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Analysis.hh"

#include "LYSimPrimaryGeneratorAction.hh"

#include "Randomize.hh"

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4RandomDirection.hh"

#include <CLHEP/Units/PhysicalConstants.h>


//#include "G4OpticalPhoton.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LYSimPrimaryGeneratorAction::LYSimPrimaryGeneratorAction(LYSimDetectorConstruction* det)
: PhotonEnergy(2.95*eV), GammaEnergy(660*keV), BetaEnergy(511*keV)
{
	fDetector = det;
	particleSource = new G4GeneralParticleSource();

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* particle = G4OpticalPhoton::OpticalPhotonDefinition();

	particleSource->SetParticleDefinition(particle);
	particleSource->SetParticleTime(0.0*ns);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LYSimPrimaryGeneratorAction::~LYSimPrimaryGeneratorAction()
{
	delete particleSource;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LYSimPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	if (particleSource->GetParticleDefinition()->GetParticleName() == "opticalphoton")
	{
		SetOptPhotonPolar();
	}
	
	particleSource->GeneratePrimaryVertex(anEvent);

	//Analysis
	//Analysis::GetInstance()->AddPhotonCount(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Photon polarization is randomized
void LYSimPrimaryGeneratorAction::SetOptPhotonPolar()
{
	G4double angle = G4UniformRand() * 360.0*deg;
	SetOptPhotonPolar(angle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void LYSimPrimaryGeneratorAction::SetOptPhotonPolar(G4double angle)
{
	if (particleSource->GetParticleDefinition()->GetParticleName() == "opticalphoton")
	{
		G4ThreeVector normal (1., 0., 0.);
		G4ThreeVector kphoton = particleSource->GetParticleMomentumDirection();
		G4ThreeVector product = normal.cross(kphoton); 
		G4double modul2       = product*product;

		G4ThreeVector e_perpend (0., 0., 1.);
		if (modul2 > 0.) e_perpend = (1./std::sqrt(modul2))*product; 
		G4ThreeVector e_paralle    = e_perpend.cross(kphoton);

		G4ThreeVector polar = std::cos(angle)*e_paralle + std::sin(angle)*e_perpend;
		particleSource->SetParticlePolarization(polar);
	}
					
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
