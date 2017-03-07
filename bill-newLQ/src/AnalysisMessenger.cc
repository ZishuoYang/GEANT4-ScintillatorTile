#include "AnalysisMessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

AnalysisMessenger::AnalysisMessenger(Analysis* instance)
 :analysis(instance)
{
  analysisDir = new G4UIdirectory("/analysis/");
  analysisDir->SetGuidance(" Analysis commands ");

  SetTileAbsLengthCmd = new G4UIcmdWithADouble("/analysis/setTileAbsLength",this);
  SetTileAbsLengthCmd->SetGuidance("Set tile absorption length (cm) to be recorded in output file.");
  SetTileAbsLengthCmd->AvailableForStates(G4State_Idle);

}

AnalysisMessenger::~AnalysisMessenger()
{
  delete analysisDir;
  delete SetTileAbsLengthCmd;
}

void AnalysisMessenger::SetNewValue(G4UIcommand* command,G4String val)
{
	if( command == SetTileAbsLengthCmd ) {
		analysis->SetTileAbsLength(G4UIcmdWithADouble::GetNewDoubleValue(val));
	}
}
