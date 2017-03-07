#ifndef LYSimDetectorMessenger_h
#define LYSimDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

#include "LYSimDetectorConstruction.hh"

class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;
class G4UIcmdWithoutParameter;

class LYSimDetectorMessenger : public G4UImessenger
{
  public:

    LYSimDetectorMessenger(LYSimDetectorConstruction* );
    ~LYSimDetectorMessenger();

    void SetNewValue(G4UIcommand*, G4String);

  private:

		LYSimDetectorConstruction*   Detector;

		G4UIdirectory*          detDir;

		G4UIcmdWithoutParameter*   	UpdateCmd;
		G4UIcmdWithABool*           SetQuartzToggleCmd;
		G4UIcmdWithABool*			SetFiberHoleCmd;
		G4UIcmdWithABool*			SetWrappingCmd;
		G4UIcmdWithABool*			SetFiberCmd;
		G4UIcmdWithABool*			SetWLSCmd;
		G4UIcmdWithABool*			SetShieldingCmd;
		G4UIcmdWithADouble*        	SetRefIndexCmd;
		G4UIcmdWithADoubleAndUnit*	SetScintThicknessCmd;
		G4UIcmdWithADoubleAndUnit*	SetScintSizeXYCmd;
		G4UIcmdWithADoubleAndUnit*	SetScintPMTGapThicknessCmd;
		G4UIcmdWithADoubleAndUnit*	SetAngle1Cmd;
		G4UIcmdWithADoubleAndUnit*	SetAngle2Cmd;
		G4UIcmdWithADoubleAndUnit*	SetDxCmd;
		G4UIcmdWithADoubleAndUnit*	SetDyCmd;
		G4UIcmdWithADoubleAndUnit*	SetDzCmd;
		G4UIcmdWithAnInteger* 			SetIetaCmd;
		G4UIcmdWithAnInteger* 			SetLayerNoCmd;
		G4UIcmdWithAnInteger* 			SetTileTypeCmd;
		G4UIcmdWithADoubleAndUnit*	SetTileAbsLengthCmd;

};

#endif
