#include "LYSimDetectorConstruction.hh"
#include "LYSimDetectorMessenger.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4UnitsTable.hh"
#include "G4SubtractionSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Torus.hh"
#include "G4Trap.hh"
#include "G4RotationMatrix.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"
#include "G4SDManager.hh"
#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4RegionStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

using std::cos;
using std::sin;
using std::tan;
using std::atan;
using std::exp;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Initialize G4VUserDetecterConstruction class
LYSimPMTSD* LYSimDetectorConstruction::fPMTSD = NULL;

LYSimDetectorConstruction::LYSimDetectorConstruction()
: G4VUserDetectorConstruction()
{
	fdetectorMessenger = new LYSimDetectorMessenger(this);

	solidWorld = NULL;
	logicWorld = NULL;
	physWorld = NULL;

	fVacuum = fAir = fSiO2 = fPolystyrene = fPolycarbonate = fLYSO = fGaAs = fAluminum = NULL;

	fSCSN81 = fBC408 = fEJ309 = NULL;

	fScintPmtGapMat = NULL;

	fFiberCore = fFiberInnerCladding = fFiberOuterCladding = fLiquidFiberCore = NULL;

	fTyvekOpSurface = fIdealTyvekOpSurface = fUnifiedTyvekOpSurface = fUnifiedIdealTyvekOpSurface = NULL;
	fPolishedOpSurface = fIdealPolishedOpSurface = NULL;
	fMirrorOpSurface = fIdealMirrorOpSurface = NULL;

	RefIndex = 1.4;

	fUpdated = false;

	fQuartz = fSapphire = NULL;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LYSimDetectorConstruction::~LYSimDetectorConstruction()
{
	if (fdetectorMessenger) delete fdetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* LYSimDetectorConstruction::Construct()
{
	SetDefaults();
	DefineMaterials();
	DefineSurfaces();
	return ConstructDetector();
}

void LYSimDetectorConstruction::SetDefaults()
{
	//Set default parameters

	mirror_toggle = true;	//	turn on/off mirrors
	fiber_hole_toggle = false;
	wrapping_toggle = false;//	tyvek/Aluminum
	wls_toggle = true;
	shielding_toggle = false;

	fScintPmtGapMat = fWater;

	// Geometry parameters
	//
	scint_sizeXY = 100.*mm;
	scint_thickness = 4.0*mm;
	hole_radius = 0.65*mm;
	PMTwindow_sizeXY = 75.*mm; //scenario
	PMTwindow_thickness = 1.0*mm;
	Photocat_thickness = 0.3*mm; //arbirary thickness for Photocathode
	ScintPMT_gap = 0.5*mm; //rough estimate of air gap between scintillator and PMT face
	shielding_sizeXY = 2*scint_sizeXY;
	shielding_thickness = 3*mm;
	fib_radius = 0.30*mm;								//Y11 fiber radius = 0.45*mm
	fib_length = 12.*cm;	
	fibXY = 6.25*mm;		//not used
	quartz_radius = 0.35*mm;	//internal radius 0.65mm; external radius 1.00mm
	fibMinZ = -200.*mm, fibMaxZ = 200.*mm;

	Photocat_sizeXY = 10*fib_radius; //scenario
	FibPMT_gap = 0.5*mm;	//air or coupling between end of fiber and PMT window

	FibQuartz_gap = 0.2*mm;	//air gap between quartz and fiber. was 0.2*mm
//      air_gap = 0.*mm;	//air gap between quartz tube and scintillator, in the liquid case is 0.mm 
	angle1 = 0*degree;
	angle2 = 0*degree;
	Dx = 50*mm;	//50	//x-distance between two fibers  
	Dy = 100*mm;
	Dz = 3.7*mm;
	bendRadius = 30*mm;
	fibRadius = fib_radius;

//	distance = 1*mm + fibRadius; //not used
//	ieta = 29;	//not used
//	layerNo = 1;	//not used
//	tileAbsLength = 20.0*cm; //not used 
//	readoutCorner = 1;	//not used


	//world volume just needs to be big enough to accomodate everything
	world_sizeXY = 10*(scint_sizeXY);
	world_sizeZ  = 10*(fibMaxZ - fibMinZ);

}

G4VPhysicalVolume* LYSimDetectorConstruction::ConstructDetector() //Begin constructing detector in G4
{
	// Option to switch on/off checking of volumes overlaps
	//
	G4bool checkOverlaps = true;

	//
	// World
	//
	G4Box* solidWorld =
	new G4Box("World",                       //its name
	0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size

	G4LogicalVolume* logicWorld =
	new G4LogicalVolume(solidWorld,			//its solid
	fAir,  				     		//its material
	"World");					//its name
		

	G4VPhysicalVolume* physWorld =
	new G4PVPlacement(0,         //no rotation
	G4ThreeVector(),       //at (0,0,0)
	logicWorld,            //its logical volume
	"World",               //its name
	0,                     //its mother  volume
	false,                 //no boolean operation
	0,                     //copy number
	checkOverlaps);        //overlaps checking

	//logicWorld -> SetVisAttributes(new G4VisAttributes(white));
	logicWorld -> SetVisAttributes(G4VisAttributes::Invisible); //MakeInvisible

	///////////////////////////////////////////
	// Tile
	///////////////////////////////////////////


	G4ThreeVector tileCenter(0, 0, 0);

	
	/*	//Young's old code for megatile 
	G4VSolid* solidTile =
	ConstructTileSolid("TileTrap",
	angle1, 			//angle measured ccw from y axis for the side at -x
	angle2,				//angle measured ccw from y axis for the side at +x
	Dx,				//length along x of side at y=+Dy
	Dy,					//length along y
	Dz,					//length along z
	tileCenter				//coordinate of corner at -x, -y, -z with respect to origin
	);

	G4LogicalVolume* logicTile =
	new G4LogicalVolume(solidTile,
	fSCSN81,
	"Tile");

	G4ThreeVector TileOffset(0, 0, 0);
	//TileOffset -= tileCenter;
	G4VPhysicalVolume* physTile =
	new G4PVPlacement(0,
	TileOffset,
	logicTile,
	"Tile",
	logicWorld,
	false,
	0,
	checkOverlaps);

	minZposition = -0.5*scint_thickness;
	maxZposition = 0.5*scint_thickness;

	if(wrapping_toggle)
	{
		G4LogicalBorderSurface* TileTyvekSurface =
			new G4LogicalBorderSurface("TileTyvekSurface",
									   physTile,
									   physWorld,
									   fUnifiedIdealTyvekOpSurface);
	}
	else
	{
		G4LogicalBorderSurface* TileAirSurface =
			new G4LogicalBorderSurface("TileAirSurface",
									   physTile,
									   physWorld,
									   fPolishedOpSurface);
	}

	//Tile visualization attributes
	G4VisAttributes * TileVisAtt = new G4VisAttributes(G4Colour(1.,1.,0.));
	TileVisAtt->SetForceWireframe(true);
	TileVisAtt->SetVisibility(true);
	logicTile->SetVisAttributes(TileVisAtt);
	*/

	G4VSolid* solidTileBox = 
	new G4Box("TileBox", 
	0.5*scint_sizeXY, 
	0.5*scint_sizeXY, 
	0.5*scint_thickness); //create Tile box that is the size of scintillator tile
	
/*	G4VSolid* solidTileHole1 = 
	new G4Tubs("TileHole1",
	0.0,
	fib_radius + FibQuartz_gap + quartz_radius + air_gap,
	0.5*fib_length,
	0.,
	2.*pi);				//create hole for FiberGroove1
	
	G4RotationMatrix* rotHole1 = new G4RotationMatrix;
	rotHole1->rotateX(pi/2*rad); 


	G4double hole1_x = -0.5*Dx - fib_radius;
        G4double hole1_y = 0 ;
    	G4double hole1_z = 0;

	G4ThreeVector transHole1;
	transHole1 = G4ThreeVector(hole1_x, hole1_y, hole1_z);
	
	G4VSolid* solidTileHole2 = 
	new G4Tubs("TileHole2",
	0.0,
	fib_radius + FibQuartz_gap + quartz_radius + air_gap,
	0.5*fib_length,
	0.,
	2.*pi);				//create hole for FiberGroove2
	
	G4RotationMatrix* rotHole2 = new G4RotationMatrix;
	rotHole2->rotateX(pi/2*rad); 


	G4double hole2_x = 0.5*Dx - fib_radius;
        G4double hole2_y = 0 ;
    	G4double hole2_z = 0;

	G4ThreeVector transHole2;
	transHole2 = G4ThreeVector(hole2_x, hole2_y, hole2_z);
	G4VSolid* solidTile;		//declares pointer solidTile
	
	G4VSolid* UnionHole = 
	new G4UnionSolid("UnionHole", solidTileHole1, solidTileHole2, 0, G4ThreeVector(Dx+2*fib_radius, 0, 0)); 	//combine hole1 and hole2 with hole2 Dx+2*fib_radius away from hole1, to be rotated


	solidTile = new G4SubtractionSolid("Tile",
					 solidTileBox, 
					 UnionHole,
					 rotHole1,
					 transHole1); //if it has hole, solidTile = (solidTileBox - solidTileHole)
*/
	
	G4VSolid* solidTile;
	solidTile = solidTileBox;
	
	G4LogicalVolume* logicTile = //define the logic volume for tile, "logicTile"
	new G4LogicalVolume(solidTile,//reference to solid volume "solidTile"
			fEJ309,//define material to be EJ309
			"Tile"); //name the logic volume "Tile"
	
	G4ThreeVector TileOffSet(0,0,0);		//offset tile center

	G4VPhysicalVolume* physTile = 
	new G4PVPlacement(0,
			TileOffSet,
			logicTile,
			"Tile",
			logicWorld,
			false,
			0,
			checkOverlaps);//creates a new G4VPhysicalVolume: physTile
	
	minZposition = -0.5*scint_thickness;
	maxZposition = 0.5*scint_thickness;

	if(wrapping_toggle) //define Tile inner surface to be either Tyvek or polished optical surface (Aluminum)
	{
		G4LogicalBorderSurface* TileTyvekSurface =
			new G4LogicalBorderSurface("TileTyvekSurface",
							physTile,
							physWorld,
							fUnifiedTyvekOpSurface);
	}
	else
	{
		G4LogicalBorderSurface* TileAirSurface = 
			new G4LogicalBorderSurface("TileAirSurface",
							physTile,
							physWorld,
							fPolishedOpSurface);
	}
	
	//Tile Visualization setting:
	G4VisAttributes * TileVisAtt = new G4VisAttributes(G4Colour(1.,1.,0.));
	TileVisAtt->SetForceWireframe(true);
	TileVisAtt->SetVisibility(true);
	logicTile->SetVisAttributes(TileVisAtt);


	////////////////////////////////////////////
	// Fiber groove and fiber
	////////////////////////////////////////////

	///Fiber grooves and fibers in parallel rather than sigma configuration///

    //FiberGroove1
    
	G4Tubs* solidFiberGroove1 =
        new G4Tubs("FiberGroove1", 0., fib_radius + FibQuartz_gap + quartz_radius, 0.5*fib_length+0.5*mm, 0., 2.*pi);//Changed from 0.5*fib_length to 0.5*fib_length+0.5*mm in order to adjust for fiber mirror , so that fiber grooves accommodate mirror volumes. 

	G4LogicalVolume* logicFiberGroove1 =
	new G4LogicalVolume(solidFiberGroove1, fAir, "FiberGroove1");

	G4RotationMatrix* rotFiber1 = new G4RotationMatrix;
	rotFiber1->rotateX(pi/2*rad); //set rotation for fiber groove placement


	G4double fibgroove1_x = -0.5*Dx - fib_radius;
        G4double fibgroove1_y = 0 ;
    	G4double fibgroove1_z = 0;

    G4ThreeVector transFiberGroove1; //set location for fiber groove placement
	transFiberGroove1 = G4ThreeVector(fibgroove1_x, fibgroove1_y, fibgroove1_z);

	G4VPhysicalVolume* physFiberGroove1 =
	new G4PVPlacement(rotFiber1,
		 transFiberGroove1,
  		 logicFiberGroove1,
		"FiberGroove1",
		 logicTile,
		 false,
		 0,
    checkOverlaps);

    //FiberGroove2

	G4Tubs* solidFiberGroove2 =
            new G4Tubs("FiberGroove2", 0., fib_radius + FibQuartz_gap + quartz_radius, 0.5*fib_length+0.5*mm, 0., 2.*pi);

	G4LogicalVolume* logicFiberGroove2 =
	new G4LogicalVolume(solidFiberGroove2, fAir, "FiberGroove2");

	G4RotationMatrix* rotFiber2 = new G4RotationMatrix;
	rotFiber2->rotateX(pi/2*rad);
	rotFiber2->rotateY(angle2);  //angle2

	G4double fibgroove2_x = 0.5*Dx + fib_radius;
	G4double fibgroove2_y = 0;	//original: 0.5*Dy - 0.5*fib_length
	G4double fibgroove2_z = 0;

	G4ThreeVector transFiberGroove2;
	transFiberGroove2 = G4ThreeVector(fibgroove2_x, fibgroove2_y, fibgroove2_z);

	G4VPhysicalVolume* physFiberGroove2 =
	new G4PVPlacement(rotFiber2,
   			transFiberGroove2,
			logicFiberGroove2,
			"FiberGroove2",
    			logicTile,
			false,
			0,
			checkOverlaps);

    //Fiber Groove visualization
	G4VisAttributes * FibGrooveVisAtt = new G4VisAttributes(G4Colour(0.7,0.2,0.2));
	FibGrooveVisAtt->SetForceWireframe(true);
	FibGrooveVisAtt->SetVisibility(false);
	logicFiberGroove1->SetVisAttributes(FibGrooveVisAtt);
	logicFiberGroove2->SetVisAttributes(FibGrooveVisAtt);

	//Fiber1

	G4Tubs* solidFiberCore1 =
	new G4Tubs("FiberCore1", 0., 0.97*fib_radius, 0.5*fib_length, 0., 2.*pi);

	G4LogicalVolume* logicFiberCore1 =
	new G4LogicalVolume(solidFiberCore1,
	fFiberCore,	//or fLiquidFiberCore
	"FiberCore1");

    	G4ThreeVector transFiber1;	
	transFiber1 = G4ThreeVector(0, 0, -0.5*mm);	

	G4VPhysicalVolume* physFiberCore1 =
	new G4PVPlacement(0,
		transFiber1,		//G4ThreeVector()
		logicFiberCore1,
		"FiberCore1",
		logicFiberGroove1,
		false,
		0,
		checkOverlaps);
		
		//Cladding surrounding core, usually 3% of fiber radius
	G4Tubs* solidFiberCladding1 =
	new G4Tubs("FiberCladding1", 0.97*fib_radius, fib_radius, 0.5*fib_length, 0., 2.*pi);

	G4LogicalVolume* logicFiberCladding1 =
	new G4LogicalVolume(solidFiberCladding1,
	fQuartz,	//fFiberOuterCladding
	"FiberCladding1");

	G4VPhysicalVolume* physFiberCladding1 =
	new G4PVPlacement(0,
		transFiber1,		//G4ThreeVector()
		logicFiberCladding1,
		"FiberCladding1",
		logicFiberGroove1,
		false,
		0,
		checkOverlaps);

    	//Fiber2

    G4Tubs* solidFiberCore2 =
	new G4Tubs("FiberCore2", 0., 0.97*fib_radius, 0.5*fib_length, 0., 2.*pi);

	G4LogicalVolume* logicFiberCore2 =
	new G4LogicalVolume(solidFiberCore2,
	fFiberCore,	//fLiquidFiberCore
	"FiberCore2");

    G4ThreeVector transFiber2;
	transFiber2 = G4ThreeVector(0, 0, -0.5*mm);

	G4VPhysicalVolume* physFiberCore2 =
	new G4PVPlacement(0,
  		transFiber2,
 		logicFiberCore2,
		"FiberCore2",
		logicFiberGroove2,
		false,
		0,
		checkOverlaps);

	G4Tubs* solidFiberCladding2 =
	new G4Tubs("FiberCladding2", 0.97*fib_radius, fib_radius, 0.5*fib_length, 0., 2.*pi);

	G4LogicalVolume* logicFiberCladding2 =
	new G4LogicalVolume(solidFiberCladding2,
	fQuartz,	//fFiberOuterCladding
	"FiberCladding2");

	G4VPhysicalVolume* physFiberCladding2 =
	new G4PVPlacement(0,
		transFiber2,		
		logicFiberCladding2,
		"FiberCladding2",
		logicFiberGroove2,
		false,
		0,
		checkOverlaps);
    	//Fiber Visualization

	G4VisAttributes * FibCoreVisAtt = new G4VisAttributes(G4Colour(1.0,0.8,0.0));
	FibCoreVisAtt->SetForceWireframe(true);
	FibCoreVisAtt->SetVisibility(true); 
	logicFiberCore1->SetVisAttributes(FibCoreVisAtt);
	logicFiberCore2->SetVisAttributes(FibCoreVisAtt);

	G4VisAttributes * FibCladdingVisAtt = new G4VisAttributes(G4Colour(1.0,0.6,0.1));
	FibCladdingVisAtt->SetForceWireframe(true);
	FibCladdingVisAtt->SetVisibility(true); 
	logicFiberCladding1->SetVisAttributes(FibCladdingVisAtt);
	logicFiberCladding2->SetVisAttributes(FibCladdingVisAtt);

    //Air Gap between fiber and quartz
/*  	  G4Tubs* solidAirGap1;
  	  G4Tubs* solidAirGap2;

        solidAirGap1 =
            new G4Tubs("AirGap1", fib_radius, fib_radius + FibQuartz_gap, 0.5*fib_length+ 0.5*mm, 0., 2.*pi); 	//was quartz_radius, quartz_radius + air_gap

        solidAirGap2 =
            new G4Tubs("AirGap2", fib_radius, fib_radius + FibQuartz_gap, 0.5*fib_length+ 0.5*mm, 0., 2.*pi);	//as above


	G4LogicalVolume* logicAirGap1 =
            new G4LogicalVolume(solidAirGap1,
            fAir,
            "AirGap1");

            G4VPhysicalVolume* physAirGap1 =
            new G4PVPlacement(0,
            G4ThreeVector(),
            logicAirGap1,
            "AirGap1",
            logicFiberGroove1,
            false,
            0,
            checkOverlaps);

	G4LogicalVolume* logicAirGap2 =
            new G4LogicalVolume(solidAirGap2,
            fAir,
            "AirGap2");

            G4VPhysicalVolume* physAirGap2 =
            new G4PVPlacement(0,
            G4ThreeVector(),
            logicAirGap2,
            "AirGap2",
            logicFiberGroove2,
            false,
            0,
            checkOverlaps);

	G4VisAttributes * AirGapVisAtt = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
	AirGapVisAtt->SetForceWireframe(true);
	AirGapVisAtt->SetVisibility(false);
	logicAirGap1->SetVisAttributes(AirGapVisAtt);
	logicAirGap2->SetVisAttributes(AirGapVisAtt);
*/
    //Quartz Coating

            G4Tubs* solidQuartz1 =
            new G4Tubs("Quartz1", fib_radius + FibQuartz_gap, fib_radius + FibQuartz_gap + quartz_radius, 0.5*fib_length+0.5*mm, 0., 2.*pi);//again changed to 0.5*fib_length+0.5*mm to make room for mirrors

            G4LogicalVolume* logicQuartz1 =
            new G4LogicalVolume(solidQuartz1,
            fQuartz,	//fQuartz or fSapphire
            "Quartz1");

            G4VPhysicalVolume* physQuartz1 =
            new G4PVPlacement(0,
            G4ThreeVector(),
            logicQuartz1,
            "Quartz1",
            logicFiberGroove1,
            false,
            0,
            checkOverlaps);

            G4Tubs* solidQuartz2 =
            new G4Tubs("Quartz2", fib_radius + FibQuartz_gap, fib_radius + FibQuartz_gap + quartz_radius, 0.5*fib_length+0.5*mm, 0., 2.*pi);

            G4LogicalVolume* logicQuartz2 =
            new G4LogicalVolume(solidQuartz2,
            fQuartz,
            "Quartz2");

            G4VPhysicalVolume* physQuartz2 =
            new G4PVPlacement(0,
            G4ThreeVector(),
            logicQuartz2,
            "Quartz2",
            logicFiberGroove2,
            false,
            0,
            checkOverlaps);

            G4VisAttributes * QuartzVisAtt = new G4VisAttributes(G4Colour(0.2,0.2,1.0));
            QuartzVisAtt->SetForceWireframe(true);
            QuartzVisAtt->SetVisibility(false);
            logicQuartz1->SetVisAttributes(QuartzVisAtt);
            logicQuartz2->SetVisAttributes(QuartzVisAtt);

/*	//Air Gap between Fiber and Quartz	//duplicated. volume between fiber and quartz already defined as air. --Zishuo
	
            G4Tubs* solidFibQuartzGap1 =
            new G4Tubs("FibQuartzGap1", fib_radius, fib_radius + FibQuartz_gap, 0.5*fib_length+1.0*mm, 0., 2.*pi);

            G4LogicalVolume* logicFibQuartzGap1 =
            new G4LogicalVolume(solidFibQuartzGap1,
            fAir,
            "FibQuartzGap1");

            G4VPhysicalVolume* physFibQuartzGap1 =
            new G4PVPlacement(0,
            G4ThreeVector(),
            logicFibQuartzGap1,
            "FibQuartzGap1",
            logicFiberGroove1,
            false,
            0,
            checkOverlaps);

            G4Tubs* solidFibQuartzGap2 =
            new G4Tubs("FibQuartzGap2", fib_radius, fib_radius + FibQuartz_gap, 0.5*fib_length+1.0*mm, 0., 2.*pi);

            G4LogicalVolume* logicFibQuartzGap2 =
            new G4LogicalVolume(solidFibQuartzGap2,
            fAir,
            "FibQuartzGap2");

            G4VPhysicalVolume* physFibQuartzGap2 =
            new G4PVPlacement(0,
            G4ThreeVector(),
            logicFibQuartzGap2,
            "FibQuartzGap2",
            logicFiberGroove2,
            false,
            0,
            checkOverlaps);

            G4VisAttributes * FibQuartzGapVisAtt = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
            FibQuartzGapVisAtt->SetForceWireframe(true);
            FibQuartzGapVisAtt->SetVisibility(false);
            logicFibQuartzGap1->SetVisAttributes(FibQuartzGapVisAtt);
            logicFibQuartzGap2->SetVisAttributes(FibQuartzGapVisAtt);
*/


	///////////////////////////////////////////
	// Mirror
	///////////////////////////////////////////


	if(mirror_toggle){       
 	
	//Mirror corresponding to fiber 1
        G4Tubs* solidMirror1;

        solidMirror1 =
        new G4Tubs("Mirror1",
        0.,
        fib_radius,
        0.5*Photocat_thickness,
        0.,
        2.*pi);
    
    
	G4LogicalVolume* logicMirror1 =
	new G4LogicalVolume(solidMirror1,
	fAluminum,
	"Mirror1");

	G4RotationMatrix* rotMirror1 = rotFiber1;

	G4double mirror1_x = fibgroove1_x;
	G4double mirror1_y = fibgroove1_y + 0.5*fib_length + 0.5*Photocat_thickness - 0.5*mm;
	G4double mirror1_z = fibgroove1_z;

	G4ThreeVector transMirror1;
	transMirror1 = G4ThreeVector(mirror1_x, mirror1_y, mirror1_z);

	G4VPhysicalVolume* physMirror1 =
	new G4PVPlacement(rotMirror1,
	transMirror1,
	logicMirror1,
	"Mirror1",
	logicWorld,
	false,
	0,
	checkOverlaps);

	//Mirror corresponding to fiber 2

	G4Tubs* solidMirror2;

        solidMirror2 =
        new G4Tubs("Mirror2",
        0.,
        fib_radius,
        0.5*Photocat_thickness,
        0.,
        2.*pi);

	G4LogicalVolume* logicMirror2 =
	new G4LogicalVolume(solidMirror2,
	fAluminum,
	"Mirror2");

	G4RotationMatrix* rotMirror2 = rotFiber2;

	G4double mirror2_x = fibgroove2_x;
	G4double mirror2_y = fibgroove2_y+(0.5*fib_length+0.5*Photocat_thickness) - 0.5*mm;	
	G4double mirror2_z = fibgroove2_z;

	G4ThreeVector transMirror2;
	transMirror2 = G4ThreeVector(mirror2_x, mirror2_y, mirror2_z);

	G4VPhysicalVolume* physMirror2 =
	new G4PVPlacement(rotMirror2,
	transMirror2,
	logicMirror2,
	"Mirror2",
	logicWorld,
	false,
	0,
	checkOverlaps);

	//define surface using SkinSurface
/*	G4LogicalSkinSurface* MirrorSurface1 = 
			new G4LogicalSkinSurface("MirrirSurface1",
									logicMirror1,
									fMirrorOpSurface);
	
	
	G4LogicalSkinSurface* MirrorSurface2 = 
			new G4LogicalSkinSurface("MirrirSurface2",
									logicMirror2,
									fMirrorOpSurface);
*/
	
	//define surface using BorderSurface
	G4LogicalBorderSurface* MirrorSurface1 =
			new G4LogicalBorderSurface("MirrorSurface1",
									   physFiberGroove1,
									   physMirror1,
									   fMirrorOpSurface);

	
//	G4LogicalBorderSurface* MirrorSurface11 =
//			new G4LogicalBorderSurface("MirrorSurface11",
//									   physMirror1,
//									   physFiberGroove1,
//									   fMirrorOpSurface);

	G4LogicalBorderSurface* MirrorSurface2 =
			new G4LogicalBorderSurface("MirrorSurface2",
									   physFiberGroove2,
									   physMirror2,
									   fMirrorOpSurface);

	
//	G4LogicalBorderSurface* MirrorSurface22 =
//			new G4LogicalBorderSurface("MirrorSurface22",
//									   physMirror2,
//									   physFiberGroove2,
//									   fMirrorOpSurface);


	// Instantiation of a set of visualization attributes with grey colour
	G4VisAttributes * MirrorVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
	// Set the forced wireframe style
	MirrorVisAtt->SetForceWireframe(true);
	MirrorVisAtt->SetVisibility(true);
 	logicMirror1->SetVisAttributes(MirrorVisAtt);
 	logicMirror2->SetVisAttributes(MirrorVisAtt);
	
	}

	///////////////////////////////////////////
	// PMT
	///////////////////////////////////////////

    //Photocathode corresponding to fiber 1

	G4Tubs* solidPhotocat1 =
	new G4Tubs("Photocathode",
	0.,
	0.5*Photocat_sizeXY,
	0.5*Photocat_thickness,
	0.,
	2.*pi);

	G4LogicalVolume* logicPhotocat1 =
	new G4LogicalVolume(solidPhotocat1,
	fGaAs, 
	"Photocathode");

	G4double Photocat1_x = fibgroove1_x;
	G4double Photocat1_y = fibgroove1_y - (0.5*fib_length + 0.5*Photocat_thickness + 1*mm);//0.5*mm
	G4double Photocat1_z = fibgroove1_z;

	G4RotationMatrix* rotPhotocat1 = rotFiber1;

	G4ThreeVector transPhotocat1;
	transPhotocat1 = G4ThreeVector(Photocat1_x, Photocat1_y, Photocat1_z);

	G4VPhysicalVolume* physPhotocat1 =
	new G4PVPlacement(rotPhotocat1,
	transPhotocat1,
	logicPhotocat1,
	"Photocathode1",
	logicWorld,
	false,
	0,
	checkOverlaps);

	//Photocathode corresponding to fiber 2

	G4Tubs* solidPhotocat2 =
	new G4Tubs("Photocathode",
	0.,
	0.5*Photocat_sizeXY,
	0.5*Photocat_thickness,
	0.,
	2.*pi);

	G4LogicalVolume* logicPhotocat2 =
	new G4LogicalVolume(solidPhotocat2,
	fGaAs, 
	"Photocathode");

	G4RotationMatrix* rotPhotocat2 = rotFiber2;

	G4double Photocat2_x = fibgroove2_x;
	G4double Photocat2_y = fibgroove2_y-(0.5*fib_length+0.5*Photocat_thickness + 1.0*mm);//0.5*mm
	G4double Photocat2_z = fibgroove2_z;

	G4ThreeVector transPhotocat2;
	transPhotocat2 = G4ThreeVector(Photocat2_x, Photocat2_y, Photocat2_z);

	G4VPhysicalVolume* physPhotocat2 =
	new G4PVPlacement(rotPhotocat2,
	transPhotocat2,
	logicPhotocat2,
	"Photocathode2",
	logicWorld,
	false,
	0,
	checkOverlaps);


	G4OpticalSurface* PMTOpSurface = new G4OpticalSurface("PMT_Surface");
	G4LogicalSkinSurface* PMTSurface1 = new G4LogicalSkinSurface("name",logicPhotocat1,PMTOpSurface);
	G4LogicalSkinSurface* PMTSurface2 = new G4LogicalSkinSurface("name",logicPhotocat2,PMTOpSurface);

	PMTOpSurface -> SetType(dielectric_metal);
	PMTOpSurface -> SetModel(unified);

	G4MaterialPropertiesTable *OpSurfaceProperty = new G4MaterialPropertiesTable();

	//Hamamatsu R7600-00-4M
	const G4int numentries = 11;
	G4double energies[numentries] =
	{1.0*eV,1.5*eV,2.0*eV,2.5*eV,3.0*eV,3.5*eV,4.0*eV,4.5*eV,5.0*eV,5.5*eV,6.0*eV};
	G4double reflectivity[numentries] =
	{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
	G4double perfectefficiency[numentries] =
	{1.0, 1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
	G4double specsheetefficiency[numentries] =
	{0.01,0.01,0.02,0.18,0.23,0.25,0.20,0.10,0.09,0.07,0.05};

	OpSurfaceProperty -> AddProperty("REFLECTIVITY",energies,reflectivity,numentries);
	OpSurfaceProperty -> AddProperty("EFFICIENCY",energies,perfectefficiency,numentries);	//change from perfectefficiency to specsheetefficiency

	PMTOpSurface -> SetMaterialPropertiesTable(OpSurfaceProperty);
	if(!fPMTSD)
	{
		fPMTSD = new LYSimPMTSD("/LYSimPMT");
		G4SDManager* sdman = G4SDManager::GetSDMpointer();
		sdman->AddNewDetector(fPMTSD);
	}

	//Photocathode visualization attributes
	G4VisAttributes * PhotocatVisAtt = new G4VisAttributes(G4Colour(0.7,0.7,0.7));
	PhotocatVisAtt->SetForceWireframe(false);
	PhotocatVisAtt->SetVisibility(true);
	logicPhotocat1->SetVisAttributes(PhotocatVisAtt);
	logicPhotocat2->SetVisAttributes(PhotocatVisAtt);

	//
	//always return the physical World
	//
	return physWorld;
}
/*
G4VSolid* LYSimDetectorConstruction::ConstructTileSolid (	const G4String& name,
															G4double angle1,
															G4double angle2,
															G4double Dx2,
															G4double Dy,
															G4double Dz,
															G4ThreeVector& center)
{
	G4double x2offset = -Dy*tan(angle1);
	G4double Dx1 = Dx2 + x2offset + Dy*tan(angle2);

	G4ThreeVector centerOfGravity(0,0,0);
	corners[0] = G4ThreeVector(0.,				0.,		0.);
	corners[1] = G4ThreeVector(Dx1,				0.,		0.);
	corners[2] = G4ThreeVector(x2offset,		Dy,		0.);
	corners[3] = G4ThreeVector(Dx2+x2offset,	Dy,		0.);
	corners[4] = G4ThreeVector(0.,				0.,		Dz);
	corners[5] = G4ThreeVector(Dx1,				0.,		Dz);
	corners[6] = G4ThreeVector(x2offset,		Dy,		Dz);
	corners[7] = G4ThreeVector(Dx2+x2offset,	Dy,		Dz);

	for(int i = 0; i < 8; i++){
		centerOfGravity += corners[i];
	}
	centerOfGravity /= 8;

	for(int i = 0; i < 8; i++){
		corners[i] -= centerOfGravity;
	}

	center = centerOfGravity;

	G4VSolid* solidTile =
	new G4Trap(name, corners);

	return solidTile;
}
*/



void LYSimDetectorConstruction::DefineMaterials()
{
	// Get nist material manager
	G4NistManager* nist = G4NistManager::Instance();
	// Get elements
	G4Element *C = nist->FindOrBuildElement("C");
	G4Element *H = nist->FindOrBuildElement("H");
	G4Element *Si = nist->FindOrBuildElement("Si");
	G4Element *O = nist->FindOrBuildElement("O");
	G4Element *Sb = nist->FindOrBuildElement("Sb");
	G4Element *Rb = nist->FindOrBuildElement("Rb");
	G4Element *Cs = nist->FindOrBuildElement("Cs");
	G4Element *Lu = nist->FindOrBuildElement("Lu");
	G4Element *Y = nist->FindOrBuildElement("Y");
	G4Element *Ce = nist->FindOrBuildElement("Ce");
	G4Element *La = nist->FindOrBuildElement("La");
	G4Element *Br = nist->FindOrBuildElement("Br");
	G4Element *Na = nist->FindOrBuildElement("Na");
	G4Element *I = nist->FindOrBuildElement("I");
	G4Element *Tl = nist->FindOrBuildElement("Tl");
	G4Element *Al = nist->FindOrBuildElement("Al");
	//Making quartz for capillary fibers

    const G4int NUMENTRIES = 3;
    G4double density,      // density
    a,                   // atomic mass
    z;                   // atomic number
    G4String name,         // name
    symbol;              // symbol
    G4int ncomponents,     // n components
    iz,                  // number of protons
    in;                  // number of nuceons
    G4double abundance,    // abundance
    temperature,         // temperature
    pressure;            // pressure

	G4Material* quartz = new G4Material("quartz", 2.200*g/cm3, 2, kStateSolid);
    quartz->AddElement(Si, 1);
    quartz->AddElement(O , 2);

    G4double quartz_PP[NUMENTRIES]   = { 1.93*eV, 2.44*eV, 3.06*eV}; // lambda range 4 ri
    G4double quartz_RIND[NUMENTRIES] = { 1.46, 1.46, 1.46};     // ref index 	{ 1.54, 1.55, 1.56 }
    G4double quartz_ABSL[NUMENTRIES] = { 10.0*m, 10.0*m, 10.0*m };// atten length
    G4MaterialPropertiesTable *quartz_mt = new G4MaterialPropertiesTable();
    quartz_mt->AddProperty("RINDEX", quartz_PP, quartz_RIND, NUMENTRIES);
    quartz_mt->AddProperty("ABSLENGTH", quartz_PP, quartz_ABSL, NUMENTRIES);
    quartz->SetMaterialPropertiesTable(quartz_mt);
    fQuartz = quartz;
	
	//Make Sapphire (Al2O3) for fiber tube

	const G4int nEntries = 2;
	G4Material* sapphire = new G4Material("sapphire", 3.980*g/cm3, 2, kStateSolid);
	sapphire->AddElement(Al, 2);
	sapphire->AddElement(O, 3);
	G4double PhotonEnergy[nEntries] = {2*eV, 6*eV};	//1st order approx, const. across spectrum
	G4double Rindex[nEntries] = {1.768, 1.768};	//refractive index 1.768
	G4double AbsLength[nEntries] = {10.0*m, 10.0*m}; //inaccurate Abs length--to change later
	G4MaterialPropertiesTable *sapphireMPT = new G4MaterialPropertiesTable();
	sapphireMPT->AddProperty("RINDEX", PhotonEnergy, Rindex, nEntries);
	sapphireMPT->AddProperty("ABSLENGTH",PhotonEnergy, AbsLength, nEntries);
	sapphire->SetMaterialPropertiesTable(sapphireMPT);
	fSapphire = sapphire;	//to add fSapphire	



    ///////////////////////////////////////////////////
	fVacuum = nist->FindOrBuildMaterial("G4_VACUUM");
	fAir = nist->FindOrBuildMaterial("G4_AIR");
	fSiO2 = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
	fPolystyrene = nist->FindOrBuildMaterial("G4_POLYSTYRENE");
	fPolycarbonate = nist->FindOrBuildMaterial("G4_POLYCARBONATE");
	fGaAs = nist->FindOrBuildMaterial("G4_GALLIUM_ARSENIDE");
	fPyrex = nist->FindOrBuildMaterial("G4_Pyrex_Glass");
	fWater = nist->FindOrBuildMaterial("G4_WATER");
	fAluminum = nist->FindOrBuildMaterial("G4_Al"); 
	{
		//Air definition
		/*
		const G4int nEntries = 50;

		//Photon energy corresponds to ~620nm, ~611nm, ... ~357nm wavelengths
		G4double PhotonEnergy[nEntries] =
	   {2.00*eV,2.03*eV,2.06*eV,2.09*eV,2.12*eV,
		2.15*eV,2.18*eV,2.21*eV,2.24*eV,2.27*eV,
		2.30*eV,2.33*eV,2.36*eV,2.39*eV,2.42*eV,
		2.45*eV,2.48*eV,2.51*eV,2.54*eV,2.57*eV,
		2.60*eV,2.63*eV,2.66*eV,2.69*eV,2.72*eV,
		2.75*eV,2.78*eV,2.81*eV,2.84*eV,2.87*eV,
		2.90*eV,2.93*eV,2.96*eV,2.99*eV,3.02*eV,
		3.05*eV,3.08*eV,3.11*eV,3.14*eV,3.17*eV,
		3.20*eV,3.23*eV,3.26*eV,3.29*eV,3.32*eV,
		3.35*eV,3.38*eV,3.41*eV,3.44*eV,3.47*eV};

		G4double RefractiveIndexAir[nEntries] =
	  { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
		1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
		1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
		1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
		1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00};
		*/
		const G4int nEntries = 2;

		G4double PhotonEnergy[nEntries] =
		{1.0*eV, 6.0*eV};
		G4double RefractiveIndexAir[nEntries] =
		{1.0, 1.0};

		G4MaterialPropertiesTable* MPTAir = new G4MaterialPropertiesTable();
		MPTAir->AddProperty("RINDEX", PhotonEnergy, RefractiveIndexAir, nEntries);

		fAir->SetMaterialPropertiesTable(MPTAir);
	}

	{
		//Silicon Dioxide (SiO2)
		const G4int nEntries = 50;
		//Photon energy corresponds to ~620nm, ~611nm, ... ~357nm wavelengths
		G4double PhotonEnergy[nEntries] =
	   {2.00*eV,2.03*eV,2.06*eV,2.09*eV,2.12*eV,
		2.15*eV,2.18*eV,2.21*eV,2.24*eV,2.27*eV,
		2.30*eV,2.33*eV,2.36*eV,2.39*eV,2.42*eV,
		2.45*eV,2.48*eV,2.51*eV,2.54*eV,2.57*eV,
		2.60*eV,2.63*eV,2.66*eV,2.69*eV,2.72*eV,
		2.75*eV,2.78*eV,2.81*eV,2.84*eV,2.87*eV,
		2.90*eV,2.93*eV,2.96*eV,2.99*eV,3.02*eV,
		3.05*eV,3.08*eV,3.11*eV,3.14*eV,3.17*eV,
		3.20*eV,3.23*eV,3.26*eV,3.29*eV,3.32*eV,
		3.35*eV,3.38*eV,3.41*eV,3.44*eV,3.47*eV};

		G4double RefractiveIndexFiber[nEntries] =
	  { 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60,
		1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60,
		1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60,
		1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60,
		1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60, 1.60};

		  G4double AbsWLSFiber[nEntries] =
	  { 5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,
		5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,
		5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,1.10*m,
		1.10*m,1.10*m,1.10*m,1.10*m,1.10*m,1.10*m, 1.*mm, 1.*mm, 1.*mm, 1.*mm,
		1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm};

		  G4double EmissionWLSFiber[nEntries] =
	  { 0.05, 0.10, 0.30, 0.50, 0.75, 1.00, 1.50, 1.85, 2.30, 2.75,
		3.25, 3.80, 4.50, 5.20, 6.00, 7.00, 8.50, 9.50, 11.1, 12.4,
		12.9, 13.0, 12.8, 12.3, 11.1, 11.0, 12.0, 11.0, 17.0, 16.9,
		15.0, 9.00, 2.50, 1.00, 0.05, 0.00, 0.00, 0.00, 0.00, 0.00,
		0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};

		G4MaterialPropertiesTable* MPTFiber = new G4MaterialPropertiesTable();
		MPTFiber->AddProperty("RINDEX", PhotonEnergy, RefractiveIndexFiber, nEntries)
		->SetSpline(true);
		if (wls_toggle)
		{
			MPTFiber->AddProperty("WLSABSLENGTH",PhotonEnergy,AbsWLSFiber,nEntries)
			->SetSpline(true);
			MPTFiber->AddProperty("WLSCOMPONENT",PhotonEnergy,EmissionWLSFiber,nEntries)
			->SetSpline(true);
			MPTFiber->AddConstProperty("WLSTIMECONSTANT", 0.5*ns);
		}

		fSiO2->SetMaterialPropertiesTable(MPTFiber);
	}

	{
		//Polystyrene
		const G4int nEntries = 50;
		//Photon energy corresponds to ~620nm, ~611nm, ... ~357nm wavelengths
		G4double PhotonEnergy[nEntries] =
	   {2.00*eV,2.03*eV,2.06*eV,2.09*eV,2.12*eV,
		2.15*eV,2.18*eV,2.21*eV,2.24*eV,2.27*eV,
		2.30*eV,2.33*eV,2.36*eV,2.39*eV,2.42*eV,
		2.45*eV,2.48*eV,2.51*eV,2.54*eV,2.57*eV,
		2.60*eV,2.63*eV,2.66*eV,2.69*eV,2.72*eV,
		2.75*eV,2.78*eV,2.81*eV,2.84*eV,2.87*eV,
		2.90*eV,2.93*eV,2.96*eV,2.99*eV,3.02*eV,
		3.05*eV,3.08*eV,3.11*eV,3.14*eV,3.17*eV,
		3.20*eV,3.23*eV,3.26*eV,3.29*eV,3.32*eV,
		3.35*eV,3.38*eV,3.41*eV,3.44*eV,3.47*eV};

		G4double RefractiveIndexPS[nEntries] =
	  { 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50,
		1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50,
		1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50,
		1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50,
		1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50};

		G4double AbsPS[nEntries] =
		{2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,
		2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,
		2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,
		2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,
		2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm,2.*cm};

		G4double ScintilFast[nEntries] =
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
		1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

		// Add entries into properties table
		G4MaterialPropertiesTable* MPTPolystyrene = new G4MaterialPropertiesTable();
		MPTPolystyrene->AddProperty("RINDEX",PhotonEnergy,RefractiveIndexPS,nEntries);
		MPTPolystyrene->AddProperty("ABSLENGTH",PhotonEnergy,AbsPS,nEntries);
		MPTPolystyrene->AddProperty("FASTCOMPONENT",PhotonEnergy, ScintilFast,nEntries);
		MPTPolystyrene->AddConstProperty("SCINTILLATIONYIELD",10./keV);
		MPTPolystyrene->AddConstProperty("RESOLUTIONSCALE",1.0);
		MPTPolystyrene->AddConstProperty("FASTTIMECONSTANT", 10.*ns);

		fPolystyrene->SetMaterialPropertiesTable(MPTPolystyrene);

		// Set the Birks Constant for the Polystyrene scintillator
		fPolystyrene->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
	}

	{
		//EJ-309 (scintillation not built.--Bill) 
		
		fEJ309 = nist->FindOrBuildMaterial("G4_NAPHTHALENE");	//get NIST data for EJ309 base 
		
		const G4int nEntries = 50;
		G4double PhotonEnergy[nEntries] = 		//array for energy spectrum
		{2.07*eV,2.10*eV,2.13*eV,2.16*eV,2.19*eV,
		 2.22*eV, 2.25*eV, 2.28*eV, 2.31*eV, 2.34*eV,
		 2.37*eV, 2.40*eV, 2.43*eV, 2.46*eV, 2.49*eV,
		 2.52*eV, 2.55*eV, 2.58*eV, 2.61*eV, 2.64*eV, 
		 2.67*eV, 2.70*eV, 2.73*eV, 2.76*eV, 2.79*eV,
		 2.82*eV, 2.85*eV, 2.88*eV, 2.91*eV, 2.94*eV,
		 2.97*eV, 3.00*eV, 3.03*eV, 3.06*eV, 3.09*eV,
		 3.12*eV, 3.15*eV, 3.18*eV, 3.21*eV, 3.24*eV,
		 3.27*eV, 3.30*eV, 3.33*eV, 3.36*eV, 3.39*eV,
		 3.42*eV, 3.45*eV, 3.48*eV, 3.51*eV, 3.54*eV};
	
		 G4double RefractiveIndex[nEntries] = 		//array for refractive index
		{1.57, 1.57, 1.57, 1.57, 1.57,
		 1.57, 1.57, 1.57, 1.57, 1.57,
		 1.57, 1.57, 1.57, 1.57, 1.57,
		 1.57, 1.57, 1.57, 1.57, 1.57,
		 1.57, 1.57, 1.57, 1.57, 1.57,
		 1.57, 1.57, 1.57, 1.57, 1.57,
		 1.57, 1.57, 1.57, 1.57, 1.57,
		 1.57, 1.57, 1.57, 1.57, 1.57,
		 1.57, 1.57, 1.57, 1.57, 1.57,
		 1.57, 1.57, 1.57, 1.57, 1.57};
		 
		 G4double AbsLength[nEntries] = 			//array for absorption length (bulk)
		{2*m, 2*m, 2*m, 2*m, 2*m,
		 2*m, 2*m, 2*m, 2*m, 2*m,
		 2*m, 2*m, 2*m, 2*m, 2*m,
		 2*m, 2*m, 2*m, 2*m, 2*m,
		 2*m, 2*m, 2*m, 2*m, 2*m,
		 2*m, 2*m, 2*m, 2*m, 2*m,
		 2*m, 2*m, 2*m, 2*m, 2*m,
		 2*m, 2*m, 2*m, 2*m, 2*m,
		 2*m, 2*m, 2*m, 2*m, 2*m,
		 2*m, 2*m, 2*m, 2*m, 2*m};
		 
/*		 G4double ScintilFast[nEntries] = 		//array for fast scintillation component
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
		 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

		 G4double ScintilSlow[nEntries] = 		//array for slow 
		{};
*/
		 G4cout << "Tile filled with: EJ309" << G4endl;
		 G4MaterialPropertiesTable* MPT = new G4MaterialPropertiesTable();
		 MPT->AddProperty("RINDEX",PhotonEnergy,RefractiveIndex,nEntries);	//add refractive index to properties table
		 MPT->AddProperty("ABSLENGTH", PhotonEnergy, AbsLength, nEntries);	//add ABS length to table
		 
		 fEJ309->SetMaterialPropertiesTable(MPT);		//Set MPT to fEJ309
	}




	{
		//Pyrex glass
		const G4int nEntries = 2;
		G4double PhotonEnergy[nEntries] =
		{1.0*eV, 6.0*eV};
		G4double RefractiveIndexPyrex[nEntries] =
		{1.52, 1.52};
		G4double AbsLengthPyrex[nEntries] =
		{10*m, 10*m};

		G4MaterialPropertiesTable* MPTPyrex = new G4MaterialPropertiesTable();
		MPTPyrex->AddProperty("RINDEX",PhotonEnergy,RefractiveIndexPyrex,nEntries);
		MPTPyrex->AddProperty("ABSLENGTH",PhotonEnergy,AbsLengthPyrex,nEntries);

		fPyrex->SetMaterialPropertiesTable(MPTPyrex);
	}

	{
		//Water (Coupling between Scintillator and PMT window)
		const G4int nEntries = 2;
		G4double PhotonEnergy[nEntries] =
		{1.0*eV, 6.0*eV};
		G4double RefractiveIndexWater[nEntries];
		for (int i = 0; i <nEntries; i++)
		{
			RefractiveIndexWater[i] = RefIndex;
		}
		 G4double AbsLengthWater[nEntries] =
		{10*m, 10*m};

		G4MaterialPropertiesTable* MPTWater = new G4MaterialPropertiesTable();
		MPTWater->AddProperty("RINDEX",PhotonEnergy,RefractiveIndexWater,nEntries);
		MPTWater->AddProperty("ABSLENGTH",PhotonEnergy,AbsLengthWater,nEntries);

		fWater->SetMaterialPropertiesTable(MPTWater);

		//G4cout << "Refractive Index of Scintillator-PMT gap set to " << RefIndex << G4endl;
	}

	/*
	const G4int numentrieslysolal = 55;
	G4double lysoenergieslal[numentrieslysolal] =
  {	3.542*eV, 3.493*eV, 3.444*eV, 3.397*eV, 3.351*eV,
	3.306*eV, 3.263*eV, 3.220*eV, 3.179*eV, 3.139*eV,
	3.100*eV, 3.061*eV, 3.024*eV, 2.988*eV, 2.952*eV,
	2.917*eV, 2.883*eV, 2.850*eV, 2.818*eV, 2.786*eV,
	2.755*eV, 2.725*eV, 2.695*eV, 2.666*eV, 2.638*eV,
	2.610*eV, 2.583*eV, 2.556*eV, 2.530*eV, 2.505*eV,
	2.480*eV, 2.455*eV, 2.431*eV, 2.407*eV, 2.384*eV,
	2.362*eV, 2.339*eV, 2.317*eV, 2.296*eV, 2.275*eV,
	2.254*eV, 2.234*eV, 2.214*eV, 2.194*eV, 2.175*eV,
	2.156*eV, 2.138*eV, 2.119*eV, 2.101*eV, 2.084*eV,
	2.066*eV, 2.049*eV, 2.033*eV, 2.016*eV, 2.000*eV };
	G4double lysolal[numentrieslysolal] =
	{0.001*cm, 0.0022387211*cm, 0.0050118723*cm, 0.0112201845*cm, 0.0251188643*cm, 0.0562341325*cm, 0.1258925412*cm, 0.2818382931*cm, 0.6309573445*cm, 1.4125375446*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm, 50*cm}; //scenario
	*/
	/*
	//Pre-05/01/2013
	const G4int numentrieslysolal = 10;
	G4double lysoenergieslal[numentrieslysolal] =
	{
		3.26*eV, 3.18*eV, 3.10*eV, 3.02*eV, 2.95*eV,
		2.82*eV, 2.48*eV, 2.21*eV, 1.94*eV, 1.55*eV,
	};
	G4double lysolal[numentrieslysolal] =
	{
		0.*cm, 4.*cm, 20.*cm, 60.*cm, 75.*cm,
		150.*cm, 330.*cm, 550.*cm, 640.*cm, 640.*cm
	};
	*/



	{
		//LYSO
		//------------------------------
		// common LYSO
		G4Material* LYSOtemplate = new G4Material("LYSOtemplate", 7.1*g/cm3, 5, kStateSolid);
		LYSOtemplate->AddElement(Lu, 71.43*perCent);
		LYSOtemplate->AddElement(Y, 4.03*perCent);
		LYSOtemplate->AddElement(Si, 6.37*perCent);
		LYSOtemplate->AddElement(O, 18.14*perCent);
		LYSOtemplate->AddElement(Ce, 0.02*perCent); // cooke2000

		// LYSO
		if(!fLYSO){fLYSO = new G4Material("LYSO", LYSOtemplate->GetDensity(), LYSOtemplate, kStateSolid);}

		//LYSO Sellmeier fit from Sasha
		G4double n0 = 1.43923e0;
		G4double n1 = 3.67622e-1;
		G4double lambda1 = 2.95130e2;

		const G4int numentriesrefindex = 35;
		//wavelength array in nm
		G4double lysowavelength[numentriesrefindex] =
		{
			360, 370, 380, 390, 400,
			410, 420, 430, 440, 450,
			460, 470, 480, 490, 500,
			510, 520, 530, 540, 550,
			560, 570, 580, 590, 600,
			610, 620, 630, 640, 650,
			660, 670, 680, 690, 700
		};

		G4double lysoenergies[numentriesrefindex];
		for (int i = 0; i<numentriesrefindex; i++)
		{
			lysoenergies[i] = 1239.842 / lysowavelength[i] * eV;
		}

		G4double lysorefindex[numentriesrefindex];
		for (int i = 0; i<numentriesrefindex; i++)
		{
			lysorefindex[i] =
			sqrt(1 + pow(n0, 2.0) + pow(n1, 2.0) / (1 - pow(lambda1, 2.0) / pow(lysowavelength[i], 2.0)));
		}

		G4double lysoconstrefindex[numentriesrefindex];
		for (int i = 0; i<numentriesrefindex; i++)
		{
			lysoconstrefindex[i] = 1.82;
		}

		const G4int numentrieslal = 2;
		G4double energieslal[numentrieslal] =
		{1.0*eV, 6.0*eV};
		G4double lal[numentrieslal] =
		{200*cm, 200*cm};	//scenario


		//Light Absorption Length
		//From 1mm sample transmission data
		const G4int numentrieslysolal = 10;
		G4double lysoenergieslal[numentrieslysolal] =
		{3.351*eV, 3.263*eV, 3.179*eV, 3.100*eV, 3.024*eV, 2.952*eV, 2.883*eV, 2.695*eV, 2.384*eV, 2.066*eV};
		G4double lysolal[numentrieslysolal] =
		{0.025*cm, 0.1*cm, 1*cm, 4*cm, 6*cm, 7*cm, 8*cm, 9*cm, 10*cm, 12*cm};


		//Scintillation emission spectrum (fast component)
		//Gamma-ray emission
		const G4int numentriesemissiongamma = 16;
		G4double lysoenergiesemissiongamma[numentriesemissiongamma] =
		{
			3.44*eV, 3.26*eV, 3.18*eV, 3.10*eV, 3.02*eV,
			2.95*eV, 2.88*eV, 2.82*eV, 2.70*eV, 2.58*eV,
			2.48*eV, 2.38*eV, 2.30*eV, 2.21*eV, 2.14*eV,
			1.82*eV
		};
		G4double lysoemissiongamma[numentriesemissiongamma] =
		{
			0.00, 0.06, 0.28, 0.72, 1.40,
			2.00, 2.20, 2.06, 1.48, 0.94,
			0.60, 0.40, 0.30, 0.20, 0.10,
			0.00
		};


		//Photoluminescence (theta = 10 deg)
		const G4int numentrieslysoemission = 21;
		G4double lysoenergiesemission[numentrieslysoemission] =
		{
			3.54*eV, 3.35*eV, 3.26*eV, 3.18*eV, 3.13*eV,
			3.10*eV, 3.02*eV, 2.95*eV, 2.88*eV, 2.82*eV,
			2.76*eV, 2.70*eV, 2.64*eV, 2.58*eV, 2.53*eV,
			2.48*eV, 2.43*eV, 2.38*eV, 2.30*eV, 2.21*eV,
			2.00*eV
		};
		G4double lysoemission[numentrieslysoemission] =
		{
			0, 0.26, 1.26, 2.14, 2.2,
			2.16, 2.04, 1.9, 1.64, 1.3,
			0.9, 0.62, 0.38, 0.26, 0.14,
			0.1, 0.08, 0.06, 0.04, 0.02,
			0
		};

		G4MaterialPropertiesTable* lysoprop = new G4MaterialPropertiesTable();
		lysoprop->AddProperty("FASTCOMPONENT", lysoenergiesemission, lysoemission, numentrieslysoemission);
		lysoprop->AddProperty("RINDEX",        lysoenergies, lysorefindex, numentriesrefindex);
		lysoprop->AddProperty("ABSLENGTH",     lysoenergieslal, lysolal,  numentrieslysolal); //scenario
		lysoprop->AddConstProperty("SCINTILLATIONYIELD",32./keV); // saint-gobain
		lysoprop->AddConstProperty("RESOLUTIONSCALE",1.0);
		lysoprop->AddConstProperty("FASTTIMECONSTANT",41.0*ns); // saint-gobain
		lysoprop->AddConstProperty("YIELDRATIO",1.0);
		fLYSO->SetMaterialPropertiesTable(lysoprop);
	}


	{
		if(!fFiberCore)
		{
			fFiberCore = new G4Material("FiberCorePS", 1.05*g/cm3, 2, kStateSolid);
			fFiberCore->AddElement(C, 85.71*perCent);
			fFiberCore->AddElement(H, 14.28*perCent);
		}
		//fFiberCore = nist->FindOrBuildMaterial("G4_POLYSTYRENE");


		//Fiber material definition
		
		const G4int NUMENTRIES = 2;
		G4double PhotonEnergyFiberCore[NUMENTRIES] = {1.0*eV, 6.0*eV};
		G4double RefractiveIndexFiberCore[NUMENTRIES] = {1.59, 1.59};
		G4double AbsLengthFiberCore[NUMENTRIES] =
		{2.5*m,2.5*m};	//was 2.5*m, 2.5*m

		const G4int NUMENTRIES1 = 91;
		G4double PhotonEnergy_WLS_ABS_FiberCore[NUMENTRIES1] =
	  {	1.776*eV, 1.794*eV, 1.807*eV, 1.821*eV, 1.832*eV,
		1.851*eV, 1.868*eV, 1.876*eV, 1.887*eV, 1.890*eV,
		1.902*eV, 1.908*eV, 1.917*eV, 1.926*eV, 1.932*eV,
		1.941*eV, 1.947*eV, 1.959*eV, 1.969*eV, 1.981*eV,
		1.994*eV, 2.004*eV, 2.010*eV, 2.020*eV, 2.027*eV,
		2.030*eV, 2.040*eV, 2.047*eV, 2.054*eV, 2.061*eV,
		2.071*eV, 2.081*eV, 2.088*eV, 2.099*eV, 2.110*eV,
		2.135*eV, 2.154*eV, 2.161*eV, 2.177*eV, 2.212*eV,
		2.244*eV, 2.273*eV, 2.285*eV, 2.302*eV, 2.311*eV,
		2.320*eV, 2.333*eV, 2.359*eV, 2.377*eV, 2.391*eV,
		2.410*eV, 2.424*eV, 2.443*eV, 2.458*eV, 2.478*eV,
		2.513*eV, 2.572*eV, 2.594*eV, 2.616*eV, 2.632*eV,
		2.649*eV, 2.666*eV, 2.678*eV, 2.695*eV, 2.701*eV,
		2.731*eV, 2.749*eV, 2.768*eV, 2.792*eV, 2.811*eV,
		2.824*eV, 2.831*eV, 2.850*eV, 2.877*eV, 2.904*eV,
		2.910*eV, 2.931*eV, 2.952*eV, 2.980*eV, 3.017*eV,
		3.046*eV, 3.069*eV, 3.092*eV, 3.123*eV, 3.155*eV,
		3.212*eV, 3.271*eV, 3.315*eV, 3.378*eV, 3.454*eV,
		3.522*eV};
		G4double WLS_ABSLENGTH_FiberCore[NUMENTRIES1] =
	  {	71.2971*cm, 117.49*cm, 146.611*cm, 181.757*cm, 211.883*cm,
		224.937*cm, 207.866*cm, 204.854*cm, 188.787*cm, 174.728*cm,
		155.649*cm, 139.582*cm, 128.536*cm, 131.548*cm, 141.59*cm,
		152.636*cm, 167.699*cm, 185.774*cm, 198.828*cm, 204.854*cm,
		200.837*cm, 187.782*cm, 165.69*cm, 123.515*cm, 85.3556*cm,
		67.2803*cm, 61.2552*cm, 63.2636*cm, 69.2887*cm, 86.3598*cm,
		111.464*cm, 139.582*cm, 156.653*cm, 168.703*cm, 178.745*cm,
		177.741*cm, 166.695*cm, 160.669*cm, 152.636*cm, 144.603*cm,
		136.569*cm, 129.54*cm, 119.498*cm, 108.452*cm, 99.4142*cm,
		88.3682*cm, 82.3431*cm, 84.3515*cm, 81.3389*cm, 74.3096*cm,
		65.272*cm, 56.2343*cm, 42.1757*cm, 31.1297*cm, 22.0921*cm,
		11.046*cm, 1.64583*cm, 0.51974*cm, 0.214673*cm, 0.121914*cm,
		0.0742481*cm, 0.0539618*cm, 0.0416667*cm, 0.0337031*cm, 0.0298338*cm,
		0.0277388*cm, 0.029216*cm, 0.0309561*cm, 0.0321661*cm, 0.0317524*cm,
		0.0301988*cm, 0.0278955*cm, 0.0261243*cm, 0.025*cm, 0.0261936*cm,
		0.0282951*cm, 0.0321661*cm, 0.0347711*cm, 0.0387255*cm, 0.0404713*cm,
		0.0418432*cm, 0.046801*cm, 0.0536685*cm, 0.0671769*cm, 0.0822918*cm,
		0.109722*cm, 0.147388*cm, 0.205729*cm, 0.308594*cm, 0.448864*cm,
		0.548611*cm};
		const G4int NUMENTRIES2 = 42;
		G4double PhotonEnergy_WLS_Em_FiberCore[NUMENTRIES2] =
	  {	1.993*eV, 2.029*eV, 2.070*eV, 2.109*eV, 2.153*eV,
		2.187*eV, 2.222*eV, 2.246*eV, 2.271*eV, 2.305*eV,
		2.331*eV, 2.353*eV, 2.366*eV, 2.384*eV, 2.394*eV,
		2.407*eV, 2.417*eV, 2.431*eV, 2.445*eV, 2.460*eV,
		2.475*eV, 2.490*eV, 2.510*eV, 2.520*eV, 2.535*eV,
		2.546*eV, 2.562*eV, 2.572*eV, 2.583*eV, 2.594*eV,
		2.605*eV, 2.616*eV, 2.627*eV, 2.644*eV, 2.661*eV,
		2.666*eV, 2.678*eV, 2.689*eV, 2.701*eV, 2.719*eV,
		2.749*eV, 2.780*eV };
		G4double WLS_Emission_FiberCore[NUMENTRIES2] =
	  {	0.00505051, 0.012626, 0.0252525, 0.035353, 0.0555556,
		0.0782828, 0.126263, 0.164141, 0.222222, 0.270202,
		0.315657, 0.373737, 0.444444, 0.515152, 0.580808,
		0.65404, 0.719697, 0.762626, 0.792929, 0.777778,
		0.747475, 0.70202, 0.686869, 0.69697, 0.739899,
		0.787879, 0.858586, 0.919192, 0.969697, 1,
		0.984848, 0.924242, 0.815657, 0.64899, 0.517677,
		0.39899, 0.287879, 0.186869, 0.103535, 0.0530303,
		0.0151515, 0 };
		G4MaterialPropertiesTable* MPTFiberCore = new G4MaterialPropertiesTable();
		MPTFiberCore->AddProperty("RINDEX",PhotonEnergyFiberCore,RefractiveIndexFiberCore,NUMENTRIES);
		MPTFiberCore->AddProperty("ABSLENGTH",PhotonEnergyFiberCore,AbsLengthFiberCore,NUMENTRIES);
		MPTFiberCore->AddProperty("WLSABSLENGTH",PhotonEnergy_WLS_ABS_FiberCore,WLS_ABSLENGTH_FiberCore,NUMENTRIES1);
		MPTFiberCore->AddProperty("WLSCOMPONENT",PhotonEnergy_WLS_Em_FiberCore,WLS_Emission_FiberCore,NUMENTRIES2);
		MPTFiberCore->AddConstProperty("WLSTIMECONSTANT", 0.5*ns);
		fFiberCore->SetMaterialPropertiesTable(MPTFiberCore);
	}


	{
		if (!fFiberInnerCladding)
		{
			fFiberInnerCladding = new G4Material("FiberInnerCladding", 1.05*g/cm3, 2, kStateSolid);
			fFiberInnerCladding->AddElement(C, 85.71*perCent);
			fFiberInnerCladding->AddElement(H, 14.28*perCent);
		}
		const G4int NUMENTRIES3 = 2;
		G4double PhotonEnergyFiberInnerCladding[NUMENTRIES3] = {1.0*eV, 6.0*eV};
		G4double RefractiveIndexFiberInnerCladding[NUMENTRIES3] = {1.49, 1.49};

		G4MaterialPropertiesTable* MPTFiberInnerCladding = new G4MaterialPropertiesTable();
		MPTFiberInnerCladding->AddProperty("RINDEX",PhotonEnergyFiberInnerCladding,RefractiveIndexFiberInnerCladding,NUMENTRIES3);
		fFiberInnerCladding->SetMaterialPropertiesTable(MPTFiberInnerCladding);
	}

	{
		if (!fFiberOuterCladding)
		{
			fFiberOuterCladding = new G4Material("FiberOuterCladding", 1.05*g/cm3, 2, kStateSolid);
			fFiberOuterCladding->AddElement(C, 85.71*perCent);
			fFiberOuterCladding->AddElement(H, 14.28*perCent);
		}
		const G4int NUMENTRIES4 = 2;
		G4double PhotonEnergyFiberOuterCladding[NUMENTRIES4] = {1.0*eV, 6.0*eV};
		G4double RefractiveIndexFiberOuterCladding[NUMENTRIES4] = {1.42, 1.42};

		G4MaterialPropertiesTable* MPTFiberOuterCladding = new G4MaterialPropertiesTable();
		MPTFiberOuterCladding->AddProperty("RINDEX",PhotonEnergyFiberOuterCladding,RefractiveIndexFiberOuterCladding,NUMENTRIES4);
		fFiberOuterCladding->SetMaterialPropertiesTable(MPTFiberOuterCladding);
	}

	{
		//SCSN81
		fSCSN81 = nist->FindOrBuildMaterial("G4_POLYSTYRENE");
		const G4int nEntries = 2;
		G4double PhotonEnergy[nEntries] = {1.0*eV, 6.0*eV};
		G4double RefractiveIndex[nEntries] = {1.59, 1.59};
		G4double absLength = GetTileAbsLength();
	//	G4cout << "Tile abs length set to " << G4BestUnit(absLength, "Length") << G4endl;
		G4double AbsLength[nEntries] = {absLength, absLength};
		// Add entries into properties table
		G4MaterialPropertiesTable* MPT = new G4MaterialPropertiesTable();
		MPT->AddProperty("RINDEX",PhotonEnergy,RefractiveIndex,nEntries);
		MPT->AddProperty("ABSLENGTH",PhotonEnergy,AbsLength,nEntries);
		fSCSN81->SetMaterialPropertiesTable(MPT);
	}
	
	{
		if(!fLiquidFiberCore)	//Build liquid fibercore simulating EJ309
		{
			fLiquidFiberCore = new G4Material("LiquidFiberCorePS", 0.959*g/cm3, 2, kStateSolid);
			fLiquidFiberCore->AddElement(C, 44.4*perCent);
			fLiquidFiberCore->AddElement(H, 55.6*perCent);
		}
		//fLiquidFiberCore = nist->FindOrBuildMaterial("G4_POLYSTYRENE");
		
		const G4int NUMENTRIES = 2;
		G4double PhotonEnergyFiberCore[NUMENTRIES] = {1.0*eV, 6.0*eV};
		G4double RefractiveIndexFiberCore[NUMENTRIES] = {1.57, 1.57};
		G4double AbsLengthFiberCore[NUMENTRIES] =
		{2.5*m, 2.5*m};	//was 2.5*m, 2.5*m

		const G4int NUMENTRIES1 = 91;
		G4double PhotonEnergy_WLS_ABS_FiberCore[NUMENTRIES1] =
	  {	1.776*eV, 1.794*eV, 1.807*eV, 1.821*eV, 1.832*eV,
		1.851*eV, 1.868*eV, 1.876*eV, 1.887*eV, 1.890*eV,
		1.902*eV, 1.908*eV, 1.917*eV, 1.926*eV, 1.932*eV,
		1.941*eV, 1.947*eV, 1.959*eV, 1.969*eV, 1.981*eV,
		1.994*eV, 2.004*eV, 2.010*eV, 2.020*eV, 2.027*eV,
		2.030*eV, 2.040*eV, 2.047*eV, 2.054*eV, 2.061*eV,
		2.071*eV, 2.081*eV, 2.088*eV, 2.099*eV, 2.110*eV,
		2.135*eV, 2.154*eV, 2.161*eV, 2.177*eV, 2.212*eV,
		2.244*eV, 2.273*eV, 2.285*eV, 2.302*eV, 2.311*eV,
		2.320*eV, 2.333*eV, 2.359*eV, 2.377*eV, 2.391*eV,
		2.410*eV, 2.424*eV, 2.443*eV, 2.458*eV, 2.478*eV,
		2.513*eV, 2.572*eV, 2.594*eV, 2.616*eV, 2.632*eV,
		2.649*eV, 2.666*eV, 2.678*eV, 2.695*eV, 2.701*eV,
		2.731*eV, 2.749*eV, 2.768*eV, 2.792*eV, 2.811*eV,
		2.824*eV, 2.831*eV, 2.850*eV, 2.877*eV, 2.904*eV,
		2.910*eV, 2.931*eV, 2.952*eV, 2.980*eV, 3.017*eV,
		3.046*eV, 3.069*eV, 3.092*eV, 3.123*eV, 3.155*eV,
		3.212*eV, 3.271*eV, 3.315*eV, 3.378*eV, 3.454*eV,
		3.522*eV};
		G4double WLS_ABSLENGTH_FiberCore[NUMENTRIES1] =
	  {	71.2971*cm, 117.49*cm, 146.611*cm, 181.757*cm, 211.883*cm,
		224.937*cm, 207.866*cm, 204.854*cm, 188.787*cm, 174.728*cm,
		155.649*cm, 139.582*cm, 128.536*cm, 131.548*cm, 141.59*cm,
		152.636*cm, 167.699*cm, 185.774*cm, 198.828*cm, 204.854*cm,
		200.837*cm, 187.782*cm, 165.69*cm, 123.515*cm, 85.3556*cm,
		67.2803*cm, 61.2552*cm, 63.2636*cm, 69.2887*cm, 86.3598*cm,
		111.464*cm, 139.582*cm, 156.653*cm, 168.703*cm, 178.745*cm,
		177.741*cm, 166.695*cm, 160.669*cm, 152.636*cm, 144.603*cm,
		136.569*cm, 129.54*cm, 119.498*cm, 108.452*cm, 99.4142*cm,
		88.3682*cm, 82.3431*cm, 84.3515*cm, 81.3389*cm, 74.3096*cm,
		65.272*cm, 56.2343*cm, 42.1757*cm, 31.1297*cm, 22.0921*cm,
		11.046*cm, 1.64583*cm, 0.51974*cm, 0.214673*cm, 0.121914*cm,
		0.0742481*cm, 0.0539618*cm, 0.0416667*cm, 0.0337031*cm, 0.0298338*cm,
		0.0277388*cm, 0.029216*cm, 0.0309561*cm, 0.0321661*cm, 0.0317524*cm,
		0.0301988*cm, 0.0278955*cm, 0.0261243*cm, 0.025*cm, 0.0261936*cm,
		0.0282951*cm, 0.0321661*cm, 0.0347711*cm, 0.0387255*cm, 0.0404713*cm,
		0.0418432*cm, 0.046801*cm, 0.0536685*cm, 0.0671769*cm, 0.0822918*cm,
		0.109722*cm, 0.147388*cm, 0.205729*cm, 0.308594*cm, 0.448864*cm,
		0.548611*cm};
		const G4int NUMENTRIES2 = 42;
		G4double PhotonEnergy_WLS_Em_FiberCore[NUMENTRIES2] =
	  {	1.993*eV, 2.029*eV, 2.070*eV, 2.109*eV, 2.153*eV,
		2.187*eV, 2.222*eV, 2.246*eV, 2.271*eV, 2.305*eV,
		2.331*eV, 2.353*eV, 2.366*eV, 2.384*eV, 2.394*eV,
		2.407*eV, 2.417*eV, 2.431*eV, 2.445*eV, 2.460*eV,
		2.475*eV, 2.490*eV, 2.510*eV, 2.520*eV, 2.535*eV,
		2.546*eV, 2.562*eV, 2.572*eV, 2.583*eV, 2.594*eV,
		2.605*eV, 2.616*eV, 2.627*eV, 2.644*eV, 2.661*eV,
		2.666*eV, 2.678*eV, 2.689*eV, 2.701*eV, 2.719*eV,
		2.749*eV, 2.780*eV };
		G4double WLS_Emission_FiberCore[NUMENTRIES2] =
	  {	0.00505051, 0.012626, 0.0252525, 0.035353, 0.0555556,
		0.0782828, 0.126263, 0.164141, 0.222222, 0.270202,
		0.315657, 0.373737, 0.444444, 0.515152, 0.580808,
		0.65404, 0.719697, 0.762626, 0.792929, 0.777778,
		0.747475, 0.70202, 0.686869, 0.69697, 0.739899,
		0.787879, 0.858586, 0.919192, 0.969697, 1,
		0.984848, 0.924242, 0.815657, 0.64899, 0.517677,
		0.39899, 0.287879, 0.186869, 0.103535, 0.0530303,
		0.0151515, 0 };
		G4MaterialPropertiesTable* MPTLiquidFiberCore = new G4MaterialPropertiesTable();
		MPTLiquidFiberCore->AddProperty("RINDEX",PhotonEnergyFiberCore,RefractiveIndexFiberCore,NUMENTRIES);
		MPTLiquidFiberCore->AddProperty("ABSLENGTH",PhotonEnergyFiberCore,AbsLengthFiberCore,NUMENTRIES);
		MPTLiquidFiberCore->AddProperty("WLSABSLENGTH",PhotonEnergy_WLS_ABS_FiberCore,WLS_ABSLENGTH_FiberCore,NUMENTRIES1);
		MPTLiquidFiberCore->AddProperty("WLSCOMPONENT",PhotonEnergy_WLS_Em_FiberCore,WLS_Emission_FiberCore,NUMENTRIES2);
		MPTLiquidFiberCore->AddConstProperty("WLSTIMECONSTANT", 0.5*ns);
		fLiquidFiberCore->SetMaterialPropertiesTable(MPTLiquidFiberCore);
	}
}


void LYSimDetectorConstruction::DefineSurfaces()
{
	{
		//////////////////////////////////
		//Realistic Crystal-Tyvek surface
		//////////////////////////////////
		fTyvekOpSurface = new G4OpticalSurface("TyvekOpSurface");
		fTyvekOpSurface->SetType(dielectric_LUT);
		fTyvekOpSurface->SetModel(LUT);
		fTyvekOpSurface->SetFinish(polishedtyvekair);

		const G4int num = 2;
		G4double Ephoton[num] = {1.5*eV, 8.0*eV};
		G4double Reflectivity[num] = {0.979, 0.979};

		G4MaterialPropertiesTable *MPT = new G4MaterialPropertiesTable();

		MPT->AddProperty("REFLECTIVITY", Ephoton, Reflectivity, num);

		fTyvekOpSurface->SetMaterialPropertiesTable(MPT);
	}

	{
		//////////////////////////////////
		//Ideal Crystal-Tyvek surface
		//////////////////////////////////
		fIdealTyvekOpSurface = new G4OpticalSurface("IdealTyvekOpSurface");
		fIdealTyvekOpSurface->SetType(dielectric_LUT);
		fIdealTyvekOpSurface->SetModel(LUT);
		fIdealTyvekOpSurface->SetFinish(polishedtyvekair);
	}

	{
		//////////////////////////////////
		//Unified Tyvek surface
		//////////////////////////////////
		fUnifiedTyvekOpSurface = new G4OpticalSurface("UnifiedTyvekOpSurface");
		fUnifiedTyvekOpSurface->SetType(dielectric_dielectric);
		fUnifiedTyvekOpSurface->SetModel(unified);
		fUnifiedTyvekOpSurface->SetFinish(groundbackpainted);
		fUnifiedTyvekOpSurface->SetSigmaAlpha(1.3*degree);

		const G4int NUM = 2;
		G4double pp[NUM] = {2.0*eV, 6.0*eV};
		G4double specularlobe[NUM] = {1.0, 1.0};
		G4double specularspike[NUM] = {0., 0.};
		G4double backscatter[NUM] = {0., 0.};
		G4double rindex[NUM] = {1.0, 1.0};
		G4double reflectivity[NUM] = {0.979, 0.979};
		G4double efficiency[NUM] = {0.0, 0.0};

		G4MaterialPropertiesTable* SMPT = new G4MaterialPropertiesTable();

		SMPT->AddProperty("RINDEX",pp,rindex,NUM);
		SMPT->AddProperty("SPECULARLOBECONSTANT",pp,specularlobe,NUM);
		SMPT->AddProperty("SPECULARSPIKECONSTANT",pp,specularspike,NUM);
		SMPT->AddProperty("BACKSCATTERCONSTANT",pp,backscatter,NUM);
		SMPT->AddProperty("REFLECTIVITY",pp,reflectivity,NUM);
		SMPT->AddProperty("EFFICIENCY",pp,efficiency,NUM);

		fUnifiedTyvekOpSurface -> SetMaterialPropertiesTable(SMPT);

	}


	{
		//////////////////////////////////
		//Unified Ideal Tyvek surface
		//////////////////////////////////
		fUnifiedIdealTyvekOpSurface = new G4OpticalSurface("UnifiedIdealTyvekOpSurface");
		fUnifiedIdealTyvekOpSurface->SetType(dielectric_dielectric);
		fUnifiedIdealTyvekOpSurface->SetModel(unified);
		fUnifiedIdealTyvekOpSurface->SetFinish(groundbackpainted);
		fUnifiedIdealTyvekOpSurface->SetSigmaAlpha(1.3*degree);

		const G4int NUM = 2;
		G4double pp[NUM] = {2.0*eV, 6.0*eV};
		G4double specularlobe[NUM] = {1.0, 1.0};
		G4double specularspike[NUM] = {0., 0.};
		G4double backscatter[NUM] = {0., 0.};
		G4double rindex[NUM] = {1.0, 1.0};
		G4double reflectivity[NUM] = {1.0, 1.0};
		G4double efficiency[NUM] = {0.0, 0.0};

		G4MaterialPropertiesTable* SMPT = new G4MaterialPropertiesTable();

		SMPT->AddProperty("RINDEX",pp,rindex,NUM);
		SMPT->AddProperty("SPECULARLOBECONSTANT",pp,specularlobe,NUM);
		SMPT->AddProperty("SPECULARSPIKECONSTANT",pp,specularspike,NUM);
		SMPT->AddProperty("BACKSCATTERCONSTANT",pp,backscatter,NUM);
		SMPT->AddProperty("REFLECTIVITY",pp,reflectivity,NUM);
		SMPT->AddProperty("EFFICIENCY",pp,efficiency,NUM);

		fUnifiedIdealTyvekOpSurface -> SetMaterialPropertiesTable(SMPT);

	}

	{
		//////////////////////////////////
		//Realistic polished surface
		//////////////////////////////////
		fPolishedOpSurface = new G4OpticalSurface("PolishedOpSurface");//, unified);
		fPolishedOpSurface->SetType(dielectric_metal);
		fPolishedOpSurface->SetModel(unified);
		fPolishedOpSurface->SetFinish(ground); //ground	// necessary even for polished surfaces to enable UNIFIED code
		fPolishedOpSurface->SetSigmaAlpha(1.3 * degree); // Janecek2010

		const G4int NUM = 2;
		G4double pp[NUM] = {2.0*eV, 6.0*eV};
		G4double specularlobe[NUM] = {1.0, 1.0};
		G4double specularspike[NUM] = {0., 0.};
		G4double backscatter[NUM] = {0., 0.};
		G4double reflectivity[NUM] = {0.99, 0.99};	//was 0.9
		G4double efficiency[NUM] = {0.0, 0.0};		
		
		G4MaterialPropertiesTable *PolishedOpSurfaceProperty = new G4MaterialPropertiesTable();
		
		PolishedOpSurfaceProperty->AddProperty("SPECULARLOBECONSTANT",pp,specularlobe,NUM);
		PolishedOpSurfaceProperty->AddProperty("SPECULARSPIKECONSTANT",pp,specularspike,NUM);
		PolishedOpSurfaceProperty->AddProperty("BACKSCATTERCONSTANT",pp,backscatter,NUM);
		PolishedOpSurfaceProperty->AddProperty("REFLECTIVITY",pp,reflectivity,NUM);	
		PolishedOpSurfaceProperty->AddProperty("EFFICIENCY",pp,efficiency,NUM);

		fPolishedOpSurface->SetMaterialPropertiesTable(PolishedOpSurfaceProperty);
	}

	{
		//////////////////////////////////
		//Ideal polished surface
		//////////////////////////////////
		fIdealPolishedOpSurface = new G4OpticalSurface("IdealOpSurface");
		fIdealPolishedOpSurface->SetType(dielectric_dielectric);
		fIdealPolishedOpSurface->SetModel(glisur);
		fIdealPolishedOpSurface->SetFinish(polished);
	}

	{
		//////////////////////////////////
		//Mirror surface
		//////////////////////////////////

		fMirrorOpSurface = new G4OpticalSurface("MirrorOpSurface");
		fMirrorOpSurface->SetType(dielectric_metal);
		fMirrorOpSurface->SetFinish(polished);	
		fMirrorOpSurface->SetModel(unified);

/*		G4OpticalSurface* fMirrorOpSurface = new G4OpticalSurface("MirrorSurface",
		glisur,
		ground,
		dielectric_metal,
		1.0);//
*/
		G4MaterialPropertiesTable *MirrorOpSurfaceProperty = new G4MaterialPropertiesTable();
		const G4int NUM = 2;
		G4double pp[NUM] = {1.0*eV, 6.0*eV};
		G4double reflectivity[NUM] = {1.0, 1.0};//1.0
		G4double efficiency[NUM] = {0.0, 0.0};//

		MirrorOpSurfaceProperty->AddProperty("REFLECTIVITY",pp,reflectivity,NUM);
		MirrorOpSurfaceProperty->AddProperty("EFFICIENCY",pp,efficiency,NUM);//

		fMirrorOpSurface->SetMaterialPropertiesTable(MirrorOpSurfaceProperty);
	}

	{
		//////////////////////////////////
		//Ideal mirror surface
		//////////////////////////////////
		fIdealMirrorOpSurface = new G4OpticalSurface("MirrorOpSurface");
		fIdealMirrorOpSurface->SetType(dielectric_metal);
		fIdealMirrorOpSurface->SetFinish(polished);
		fIdealMirrorOpSurface->SetModel(unified);

		G4MaterialPropertiesTable *IdealMirrorOpSurfaceProperty = new G4MaterialPropertiesTable();
		const G4int NUM = 2;
		G4double pp[NUM] = {1.0*eV, 6.0*eV};
		G4double reflectivity[NUM] = {1.0, 1.0};	
		IdealMirrorOpSurfaceProperty->AddProperty("REFLECTIVITY",pp,reflectivity,NUM);
		fIdealMirrorOpSurface->SetMaterialPropertiesTable(IdealMirrorOpSurfaceProperty);
	}

	{
		//////////////////////////////////
		//Absorbing surface
		//////////////////////////////////
		fAbsorbingOpSurface = new G4OpticalSurface("AbsorbingOpSurface");
		fAbsorbingOpSurface->SetType(dielectric_dielectric);
		fAbsorbingOpSurface->SetFinish(groundfrontpainted);
		fAbsorbingOpSurface->SetModel(unified);

		G4MaterialPropertiesTable *AbsorbingOpSurfaceProperty = new G4MaterialPropertiesTable();
		const G4int NUM = 2;
		G4double pp[NUM] = {1.0*eV, 6.0*eV};
		G4double reflectivity[NUM] = {0.0, 0.0};
		AbsorbingOpSurfaceProperty->AddProperty("REFLECTIVITY",pp,reflectivity,NUM);
		fAbsorbingOpSurface->SetMaterialPropertiesTable(AbsorbingOpSurfaceProperty);
	}
}

void LYSimDetectorConstruction::SetTileType(G4int type)
{
	if (type == 1) {
		angle1 = 0 * degree;
		angle2 = 10 * degree;
	}
	else if (type == 2) {
		angle1 = 0 * degree;
		angle2 = 5 * degree;
	}
	else if (type == 3) {
		angle1 = 5 * degree;
		angle2 = 10 * degree;
	}
	else {
		G4cout << "Tile Type error: " << type << "is not a valid value." << G4endl;
	}
}

void LYSimDetectorConstruction::UpdateGeometry()
{

  // clean-up previous geometry
  G4GeometryManager::GetInstance()->OpenGeometry();

  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  G4LogicalSkinSurface::CleanSurfaceTable();
  G4LogicalBorderSurface::CleanSurfaceTable();
  G4SurfaceProperty::CleanSurfacePropertyTable();
  DefineMaterials();
  DefineSurfaces();
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructDetector());

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  //G4RunManager::GetRunManager()->PhysicsHasBeenModified();

//  G4RegionStore::GetInstance()->UpdateMaterialList(physWorld);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
