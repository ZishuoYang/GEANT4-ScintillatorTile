#ifndef AnalysisMessenger_h
#define AnalysisMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

#include "Analysis.hh"

class Analysis;

class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

class AnalysisMessenger : public G4UImessenger
{
  public:

    AnalysisMessenger(Analysis* );
    ~AnalysisMessenger();
 
    void SetNewValue(G4UIcommand*, G4String);

  private:

    Analysis*   analysis;
 
    G4UIdirectory*          analysisDir;

    G4UIcmdWithADouble*   	SetTileAbsLengthCmd;

};

#endif
