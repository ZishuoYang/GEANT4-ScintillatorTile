#include "Analysis.hh"

//User initialization
#include "LYSimDetectorConstruction.hh"
#include "LYSimPhysicsList.hh"

//User action
#include "LYSimPrimaryGeneratorAction.hh"
#include "LYSimRunAction.hh"
#include "LYSimTrackingAction.hh"
#include "LYSimSteppingAction.hh"
#include "LYSimEventAction.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"


//Always include VI, UI
//#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
//#endif
//#ifdef G4UI_USE
#include "G4UIExecutive.hh"
//#endif

#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{

	// Construct the default run manager
	//
	G4RunManager * runManager = new G4RunManager;

	G4cout << "11111" << G4endl;
	// Set mandatory initialization classes
	//
	// Detector construction
	LYSimDetectorConstruction* detector= new LYSimDetectorConstruction();
	runManager->SetUserInitialization(detector);
	G4cout << "22222" << G4endl;
	// Physics list
	G4VUserPhysicsList* physics = new LYSimPhysicsList;
	runManager-> SetUserInitialization(physics);
	G4cout << "33333" << G4endl;

	//Construct Analysis class
	Analysis::GetInstance()->SetDetector(detector);
	
	// Set user action classes
	//    
	// Primary generator action
	runManager->SetUserAction(new LYSimPrimaryGeneratorAction(detector));
	G4cout << "44444" << G4endl;
	// Run action
	runManager->SetUserAction(new LYSimRunAction(detector));
	// Event action
	runManager->SetUserAction(new LYSimEventAction());
	// Tracking action
	runManager->SetUserAction(new LYSimTrackingAction());
	// Stepping action
	runManager->SetUserAction(new LYSimSteppingAction());
	
	// Initialize G4 kernel
	//
	runManager->Initialize();
	G4cout << "55555" << G4endl;
	//#ifdef G4VIS_USE
	// Initialize visualization
	G4VisManager* visManager = new G4VisExecutive("Quiet");
	// G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
	// G4VisManager* visManager = new G4VisExecutive("Quiet");
	visManager->Initialize();
	//#endif

	// Get the pointer to the User Interface manager
	G4UImanager* UImanager = G4UImanager::GetUIpointer();

	Analysis::GetInstance()->SetOutputFile("Analysis.txt");
	if (argc == 2) {
		G4cout<< "argv[1] is " << argv[1] <<G4endl;
		if(argv[1] == "-novis") {
			G4UIExecutive* ui = new G4UIExecutive(argc, argv);
			ui->SessionStart();
			UImanager->ApplyCommand("/control/execute init.mac");
			delete ui;
		}
		else {
			// batch mode
			G4String command = "/control/execute ";
			G4String fileName = argv[1];
			UImanager->ApplyCommand(command+fileName);
		}
	}
	else if (argc == 3) {
		// batch mode
		G4String command = "/control/execute ";
		G4String fileName = argv[1];
		std::string outFileName = argv[2];
		G4cout<<"outFileName is "<< outFileName <<G4endl; 
		Analysis::GetInstance()->SetOutputFile(outFileName);
		UImanager->ApplyCommand(command+fileName);
	}
	else {
		// interactive mode : define UI session
		//#ifdef G4UI_USE
		G4UIExecutive* ui = new G4UIExecutive(argc, argv, "tcsh");
		//#ifdef G4VIS_USE
		UImanager->ApplyCommand("/control/execute init_vis.mac"); 
		//#else
		//    UImanager->ApplyCommand("/control/execute init.mac"); 
		//#endif
		ui->SessionStart();
		delete ui;
		//#endif
	}

	// Job termination
	// Free the store: user actions, physics_list and detector_description are
	// owned and deleted by the run manager, so they should not be deleted 
	// in the main() program !

	//#ifdef G4VIS_USE
	delete visManager;
	//#endif
	delete runManager;

	return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
