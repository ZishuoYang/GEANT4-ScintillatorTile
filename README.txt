================================================================================
The following instructions are for the lxplus cluster at CERN.
IMPORTANT: Use lxplus5.cern.ch to access SL5 nodes.


Checkout instructions
----------------------
1. Create working area
mkdir scratcharea
cd scratcharea

2. Checkout source code and copy to working area
cvs co UserCode/YShin/LYSim-LO2
cp -r UserCode/YShin/LYSim-LO2 LYSim-LO2


Build and run instructions
-----------------------------
1. Create build directory
mkdir LYSim-LO2-build

2a. Check out CMSSW (Only required the first time), e.g.
cmsrel CMSSW_4_4_4

2b. cd to SCRAM based area and set up runtime environment, e.g.
cd CMSSW_4_4_4/src
cmsenv

3. Compile the source code
cd ~/scratcharea/LYSim-LO2-build
export G4BASE=/afs/cern.ch/sw/lcg/external/geant4
source $G4BASE/9.6/x86_64-slc5-gcc43-opt/CMake-setup.sh
cmake -DWITH_GEANT4_UIVIS=ON -DGeant4_DIR=/afs/cern.ch/sw/lcg/external/geant4/9.6/x86_64-slc5-gcc43-opt/lib64/Geant4-9.6.0/ $HOME/scratcharea/UserCode/YShin/LYSim-LO2
make

4a. The following commands executes the input macro files, and pipes the output to the corresponding output files
./LYSim LO2.in >& LO2.out &

4b. (optional) You can also manually execute visualization commands using the interactive mode, e.g.
./LYSim
Idle> /vis/viewer/set/viewpointThetaPhi 90 0
Idle> exit


See LYSimDetectorConstruction.cc for details on simulation setup, detector geometry and material properties

================================================================================

Description of files in LYSim project.

LYSim.cc
Main function of the LYSim program. For simulating light yield measurements.

LYSimPhysicsList
Activates required physics processes.

LYSimScintillation
Custom physics process. Inherits from G4Scintillation, with added output when scintillating.

LYSimPrimaryGeneratorAction
Initializes GeneralParticleSource. Sets random photon polarization. 
Energy and postion of generated photons are set in input file, as GPS commands,
e.g.
/gps/particle opticalphoton
/gps/pos/type Plane
/gps/pos/shape Circle
/gps/pos/centre 0. 0. 0. mm
/gps/pos/halfx 12.5 mm
/gps/pos/halfy 12.5 mm
/gps/ang/type iso
/gps/ene/type Arb
/gps/hist/type arb
/gps/hist/point 3.44e-6 0.00
/gps/hist/point 3.26e-6 0.06
/gps/hist/point 3.18e-6 0.28
/gps/hist/point 3.10e-6 0.72
/gps/hist/point 3.02e-6 1.40
/gps/hist/point 2.95e-6 2.00
/gps/hist/point 2.88e-6 2.20
/gps/hist/point 2.82e-6 2.06
/gps/hist/point 2.70e-6 1.48
/gps/hist/point 2.58e-6 0.94
/gps/hist/point 2.48e-6 0.60
/gps/hist/point 2.38e-6 0.40
/gps/hist/point 2.30e-6 0.30
/gps/hist/point 2.21e-6 0.20
/gps/hist/point 2.14e-6 0.10
/gps/hist/point 1.82e-6 0.00
/gps/hist/inter Lin

LYSimRunAction
Prepares new run and collects analysis information for a run in the Analysis code.
Calls Analysis::PrepareNewRun() and Analysis::EndOfRun() before and after each run.

LYSimEventAction
Prints event ID every 100 events.
Calls Analysis::PrepareNewEvent() and Analysis::EndOfEvent() before and after each event.

LYSimTrackingAction/LYSimTrajectory/LYSimTrajectoryPoint
Specifies track information to save. (Default should be good enough for now.)
Can also be used to control which tracks to draw and color of tracks. (Currently draws everything.)
