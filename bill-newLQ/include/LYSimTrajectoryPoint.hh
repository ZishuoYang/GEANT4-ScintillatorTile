#ifndef LYSimTrajectoryPoint_h_seen
#define LYSimTrajectoryPoint_h_seen 1

#include "globals.hh"

#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4TrajectoryPoint.hh"

#include "G4StepStatus.hh"

class G4Track;
class G4Step;
class G4VProcess;

class LYSimTrajectoryPoint : public G4TrajectoryPoint {

//--------
  public: // without description
//--------

// Constructor/Destructor

    LYSimTrajectoryPoint();
    LYSimTrajectoryPoint(const G4Track* aTrack);
    LYSimTrajectoryPoint(const G4Step* aStep);
    LYSimTrajectoryPoint(const LYSimTrajectoryPoint &right);
    virtual ~LYSimTrajectoryPoint();

// Operators

    inline void *operator new(size_t);
    inline void operator delete(void *aTrajectoryPoint);
    inline int operator==(const LYSimTrajectoryPoint& right) const
    { return (this==&right); };

// Get/Set functions

    inline G4double GetTime() const { return fTime; };
    inline const G4ThreeVector GetMomentum() const { return fMomentum; };
    inline G4StepStatus GetStepStatus() const { return fStepStatus; };
    inline G4String GetVolumeName() const { return fVolumeName; };

// Get method for HEPRep style attributes

   virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
   virtual std::vector<G4AttValue>* CreateAttValues() const;

//---------
  private:
//---------

// Member data

    G4double fTime;
    G4ThreeVector fMomentum;
    G4StepStatus fStepStatus;
    G4String fVolumeName;

};

extern G4Allocator<LYSimTrajectoryPoint> LYSimTrajPointAllocator;

inline void* LYSimTrajectoryPoint::operator new(size_t)
{
    void *aTrajectoryPoint = (void *) LYSimTrajPointAllocator.MallocSingle();
    return aTrajectoryPoint;
}

inline void LYSimTrajectoryPoint::operator delete(void *aTrajectoryPoint)
{
    LYSimTrajPointAllocator.FreeSingle(
        (LYSimTrajectoryPoint *) aTrajectoryPoint);
}

#endif
