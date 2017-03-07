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
/// \file /LYSim/include/LYSimPMTHit.hh
/// \brief Definition of the LYSimPMTHit class
/// Stores hit information (Energy, detector number and position).
//
//
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef LYSimPMTHit_h
#define LYSimPMTHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

class G4VTouchable;

//--------------------------------------------------
// LYSimPMTHit Class
//--------------------------------------------------

class LYSimPMTHit : public G4VHit
{
public:

	LYSimPMTHit();
	~LYSimPMTHit();
	LYSimPMTHit(const LYSimPMTHit &right);
	
	const LYSimPMTHit& operator=(const LYSimPMTHit& right);
	G4int operator==(const LYSimPMTHit& right) const;

	inline void *operator new(size_t);
	inline void operator delete(void *aHit);

	void Draw();
	void Print();
	
	void SetEnergy(G4double energy) {fEnergy = energy;}
	void SetPhotonCount(G4int photonCount){fPhotonCount = photonCount;}
	
	void AddEnergy(G4double energy) {fEnergy += energy;}
	void IncPhotonCount(){fPhotonCount++;}
	
	G4double GetEnergy(){ return fEnergy; }
	G4int GetPhotonCount(){ return fPhotonCount;}

private:
	G4double		fEnergy; //Total photon energy deposited in PMT
	G4int			fPhotonCount; //Total number of photons detected by PMT


};

//--------------------------------------------------
// Type Definitions
//--------------------------------------------------

typedef G4THitsCollection<LYSimPMTHit> LYSimPMTHitsCollection;

extern G4Allocator<LYSimPMTHit> LYSimPMTHitAllocator;

//--------------------------------------------------
// Operator Overloads
//--------------------------------------------------

inline void* LYSimPMTHit::operator new(size_t)
{
	void *aHit;
	aHit = (void *) LYSimPMTHitAllocator.MallocSingle();
	return aHit;
}

inline void LYSimPMTHit::operator delete(void *aHit)
{
	LYSimPMTHitAllocator.FreeSingle((LYSimPMTHit*) aHit);
}

#endif
