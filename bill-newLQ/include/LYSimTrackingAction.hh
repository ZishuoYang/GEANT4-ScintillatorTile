#ifndef LYSimTrackingAction_h
#define LYSimTrackingAction_h 1

#include "G4UserTrackingAction.hh"

class LYSimTrackingAction : public G4UserTrackingAction {

  public:

    LYSimTrackingAction() { };
    ~LYSimTrackingAction() { };

    void PreUserTrackingAction(const G4Track*);
    void PostUserTrackingAction(const G4Track*);

};

#endif
