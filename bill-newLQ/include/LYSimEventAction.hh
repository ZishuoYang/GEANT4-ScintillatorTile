// LYSimEventAction.hh


#ifndef LYSimEVENTACTION_HH_
#define LYSimEVENTACTION_HH_


#include "G4UserEventAction.hh"
#include "G4String.hh"
class G4Event;

//User event action class. Prepares new event in analysis code at beginning of event.
class LYSimEventAction : public G4UserEventAction
{
public:
	//! Default constructor
	LYSimEventAction();
	//! Default destructor
	virtual ~LYSimEventAction() {};
	//! Beginning of event
	void BeginOfEventAction(const G4Event* anEvent);
	//! Digitize hits and store information
	void EndOfEventAction(const G4Event* anEvent);
private:
};

#endif /* TEST4EVENTACTION_HH_ */
