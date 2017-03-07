#ifndef LYSimScintillation_h
#define LYSimScintillation_h 1

#include <G4Scintillation.hh>


class LYSimScintillation : public G4Scintillation
{
public:
    LYSimScintillation(const G4String &processName="Scintillation", G4ProcessType type=fElectromagnetic);
    virtual ~LYSimScintillation();

    G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step&  aStep);
};

#endif
