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
/// \file LYSim/include/LYSimPMTSD.hh
/// \brief Definition of the LYSimPMTSD class
//
//
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef LYSimPMTSD_h
#define LYSimPMTSD_h 1

#include "LYSimPMTHit.hh"

#include "G4VSensitiveDetector.hh"

class G4Step;
class G4HCofThisEvent;

class LYSimPMTSD : public G4VSensitiveDetector
{
  public:

    LYSimPMTSD(G4String name);
    ~LYSimPMTSD();


    G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
    //A version of processHits that keeps aStep constant
    G4bool ProcessHits_constStep(const G4Step* aStep,
                                 G4TouchableHistory* ROhist);

    void Initialize(G4HCofThisEvent* HCE);
    void EndOfEvent(G4HCofThisEvent* HCE);

    void clear();
    void DrawAll();
    void PrintAll();

  private:

    LYSimPMTHitsCollection* fPMTHitsCollection;
 
};

#endif
