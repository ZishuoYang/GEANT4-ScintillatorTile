#include "LYSimScintillation.hh"

#include <G4VParticleChange.hh>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

//extern std::ofstream outFile;

LYSimScintillation::LYSimScintillation(const G4String &processName, G4ProcessType type)
  : G4Scintillation(processName, type)
{
}

LYSimScintillation::~LYSimScintillation()
{
}

G4VParticleChange* LYSimScintillation::PostStepDoIt(const G4Track& aTrack, const G4Step&  aStep)
{
  G4VParticleChange *result =  G4Scintillation::PostStepDoIt(aTrack, aStep);
  G4double depenergy = aStep.GetTotalEnergyDeposit();
  if (depenergy > 0.0) {
    G4ThreeVector pos = aStep.GetPreStepPoint()->GetPosition();

/*    if (outFile.is_open()) {
      outFile << "# scintillating: " << aStep.GetTotalEnergyDeposit()/keV << " keV of "
              << aStep.GetPreStepPoint()->GetKineticEnergy()/keV << " keV "
              << "from parent ID " << aTrack.GetTrackID() << " "
              << "deposited at (" << pos.x()/mm << " mm, " << pos.y()/mm << " mm, " << pos.z()/mm << " mm)"
              << " producing " << result->GetNumberOfSecondaries() << " optical photons"
              << std::endl;
    }
    else {
*/
      G4cout << "scintillating: " << aStep.GetTotalEnergyDeposit()/keV << " keV of "
             << aStep.GetPreStepPoint()->GetKineticEnergy()/keV << " keV "
             << "from parent ID " << aTrack.GetTrackID() << " "
             << "deposited at (" << pos.x()/mm << " mm, " << pos.y()/mm << " mm, " << pos.z()/mm << " mm)"
             << " producing " << result->GetNumberOfSecondaries() << " optical photons"
             << G4endl;
//    }
  }
  return result;
}
