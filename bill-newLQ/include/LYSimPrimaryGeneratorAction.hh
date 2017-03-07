
#ifndef LYSimPrimaryGeneratorAction_h
#define LYSimPrimaryGeneratorAction_h 1

#include "LYSimDetectorConstruction.hh"

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4GeneralParticleSource;
class G4Event;
class LYSimDetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class LYSimPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    LYSimPrimaryGeneratorAction(LYSimDetectorConstruction*);
   ~LYSimPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event*);

    void SetOptPhotonPolar();
    void SetOptPhotonPolar(G4double);

  private:
  
    G4GeneralParticleSource* particleSource;
    LYSimDetectorConstruction* fDetector;
    G4double GammaEnergy;
    G4double BetaEnergy;
    G4double PhotonEnergy;
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*LYSimPrimaryGeneratorAction_h*/
