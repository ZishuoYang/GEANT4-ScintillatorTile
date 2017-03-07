#ifndef LYSimRunAction_h
#define LYSimRunAction_h 1

#include "globals.hh"

#include "G4UserRunAction.hh"

class G4Run;

//class LYSimRunActionMessenger;
class LYSimDetectorConstruction;

class LYSimRunAction : public G4UserRunAction
{
  public:

    LYSimRunAction(LYSimDetectorConstruction*);
    ~LYSimRunAction();

  public:

    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);

  private:
    LYSimDetectorConstruction* pDetectorConstruction;

};

#endif
