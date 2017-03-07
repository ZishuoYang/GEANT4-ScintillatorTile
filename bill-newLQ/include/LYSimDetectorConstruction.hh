#ifndef LYSimDetectorConstruction_h
#define LYSimDetectorConstruction_h 1

#include "LYSimPMTSD.hh"

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class LYSimDetectorMessenger;
class G4Box;
class G4LogicalVolume;
class G4Material;
class G4OpticalSurface;
class LYSimPMTSD;

class LYSimDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    LYSimDetectorConstruction();
    virtual ~LYSimDetectorConstruction();

    virtual G4VPhysicalVolume* Construct();

    void UpdateGeometry();

    	//Toggle quartz on or off
	void    SetQuartzToggle     (G4bool toggle) {quartz_toggle = toggle;}
	G4bool  GetQuartzToggle     ()  {return quartz_toggle;}

	//Toggle mirror on or off
	void	SetMirrorToggle     (G4bool toggle) {mirror_toggle = toggle;}
	G4bool  GetMirrorToggle     ()  {return mirror_toggle;}

	//Toggle hole on or off
	void	SetFiberHoleToggle	(G4bool toggle)	{fiber_hole_toggle = toggle;}
	G4bool	GetFiberHoleToggle	()	{return fiber_hole_toggle;}

	//Toggle Tyvek on or off
	void	SetWrappingToggle	(G4bool toggle)	{wrapping_toggle = toggle;}
	G4bool	GetWrappingToggle	()	{return wrapping_toggle;}

	//Toggle fiber on or off
	void	SetFiberToggle	(G4bool toggle)	{fiber_toggle = toggle;}
	G4bool	GetFiberToggle	()	{return fiber_toggle;}

	//Toggle WLS on or off in fiber
	void	SetWLSToggle	(G4bool toggle)	{wls_toggle = toggle;}
	G4bool	GetWLSToggle	()	{return wls_toggle;}

	//Toggle shielding on or off
	void	SetShieldingToggle	(G4bool toggle)	{shielding_toggle = toggle;}
	G4bool	GetShieldingToggle	()	{return shielding_toggle;}

	//Set ref index of fiber
	void		SetRefIndex	(G4double iRefIndex)	{RefIndex = iRefIndex;}
	G4double	GetRefIndex	()	{return RefIndex;}

	void		SetScintThickness	(G4double iscint_thickness)	{scint_thickness = iscint_thickness;}
	G4double	GetScintThickness()	{return scint_thickness;}

	void		SetScintSizeXY	(G4double iscint_sizeXY)	{scint_sizeXY = iscint_sizeXY;}
	G4double	GetScintSizeXY()	{return scint_sizeXY;}

	void		SetScintPMTGapThickness	(G4double ithickness)	{ScintPMT_gap = ithickness;}
	G4double	GetScintPMTGapThickness()	{return ScintPMT_gap;}

	void		SetTileAbsLength	(G4double iAbsLength)	{tileAbsLength = iAbsLength;}
	G4double	GetTileAbsLength()	{return tileAbsLength;}

	void		SetAngle1	(G4double iAngle1)	{angle1 = iAngle1;}
	G4double	GetAngle1() {return angle1;}

	void		SetAngle2	(G4double iAngle2)	{angle2 = iAngle2;}
	G4double	GetAngle2() {return angle2;}

	void		SetDx	(G4double iDx)	{Dx = iDx;}
	G4double	GetDx() {return Dx;}

	void		SetDy	(G4double iDy)	{Dy = iDy;}
	G4double	GetDy() {return Dy;}

	void		SetDz	(G4double iDz)	{Dz = iDz;}
	G4double	GetDz() {return Dz;}

	void		SetBendRadius	(G4double iBendRadius)	{bendRadius = iBendRadius;}
	G4double	GetBendRadius() {return bendRadius;}

	void		SetDistance	(G4double iDistance)	{distance = iDistance;}
	G4double	GetDistance() {return distance;}

	void		SetFibRadius	(G4double iFibRadius)	{fibRadius = iFibRadius;}
	G4double	GetFibRadius() {return fibRadius;}

	void    SetIeta	(G4int integer)	{ieta = integer;}
	G4int   GetIeta() {return ieta;}

	void		SetLayerNo	(G4int integer)	{layerNo = integer;}
	G4int   GetLayerNo() {return layerNo;}

	G4double	GetMinZ() {return minZposition;}
	G4double	GetMaxZ() {return maxZposition;}

	//Tile type 1, 2, 3 based on phi segmentation
	void		SetTileType(G4int);


  private:
	//Subfunctions for cleaner code
	void DefineMaterials();
	void DefineSurfaces();
	void SetDefaults(); //*-*Doesn't seem to work right now
	G4VPhysicalVolume* ConstructDetector();
	G4VSolid* ConstructTileSolid(const	G4String& name,
										G4double angle1, //angle measured ccw from y axis for the side at -x
										G4double angle2, //angle measured ccw from y axis for the side at +x
										G4double Dx,	//length along x of side at y=+Dy
										G4double Dy,	//length along y
										G4double Dz,	//length along z
										G4ThreeVector& centerCoord); //coordinate of center of gravity w.r.t. corner 0 at -x, -y
	G4VSolid* ConstructFiberSolid(const	G4String& name,
										G4double radius,
										G4double bendRadius,
										G4double distance,		//distance from tile edge to fiber center axis
										G4double angle1,
										G4double angle2,
										G4int readoutCorner);	//index of corner with fiber readout

	//Pointer to detector messenger class
	LYSimDetectorMessenger* fdetectorMessenger;

	G4bool fUpdated;

	//Contains the minimum and maximum Z-position of the previously defined volume
	//Use this variable to define volume successively along the z-axis
	//(max for +z direction, min for -z direction)
	//Set value after defining a physical volume
	G4double minZposition;
	G4double maxZposition;

	//Pointers to main volumes
	G4Box* solidWorld;
    G4LogicalVolume* logicWorld;
    G4VPhysicalVolume* physWorld;

	//Pointers to materials
	G4Material* fVacuum;
	G4Material* fAir;
	G4Material* fSiO2;
	G4Material* fPolystyrene;
	G4Material* fPolycarbonate;
	G4Material* fLiquidFiberCore;
	G4Material* fFiberCore;
	G4Material* fFiberInnerCladding;
	G4Material* fFiberOuterCladding;
	G4Material* fLYSO;
	G4Material* fGaAs;
	G4Material* fPyrex;
	G4Material* fWater;
	G4Material* fSCSN81;
	G4Material* fBC408;
	G4Material* fEJ309;
	G4Material* fAluminum;
	
	G4Material* fScintPmtGapMat;

	G4Material* fQuartz;
	G4Material* fSapphire;

	//Pointers to surfaces
	G4OpticalSurface* fTyvekOpSurface;
	G4OpticalSurface* fIdealTyvekOpSurface;
	G4OpticalSurface* fUnifiedTyvekOpSurface;
	G4OpticalSurface* fUnifiedIdealTyvekOpSurface;
	G4OpticalSurface* fPolishedOpSurface;
	G4OpticalSurface* fIdealPolishedOpSurface;
	G4OpticalSurface* fMirrorOpSurface;
	G4OpticalSurface* fIdealMirrorOpSurface;
	G4OpticalSurface* fAbsorbingOpSurface;

	//Pointers for access to Sensitive Detector
	static LYSimPMTSD*	fPMTSD;

	//Geometry parameters
	//Important parameters marked with ////
	G4bool quartz_toggle;
	G4bool mirror_toggle;
	G4bool fiber_hole_toggle;		//hole for fiber in scintillator
	G4bool wrapping_toggle;			//(Tyvek) wrapping around scintillator
	G4bool fiber_toggle;			////
	G4bool wls_toggle;				////
	G4bool shielding_toggle;		////
	G4double RefIndex;				//Refractive index of fiber or scint-PMT gap
	G4double world_sizeXY;
	G4double world_sizeZ;
	G4double scint_sizeXY;
	G4double scint_thickness;		////
	G4double hole_radius;			////
	G4double PMTwindow_sizeXY;		////
	G4double Photocat_sizeXY;		////
	G4double PMTwindow_thickness;	////
	G4double Photocat_thickness;
	G4double ScintPMT_gap;			//thickness of gap between scintillator and PMT
	G4double shielding_sizeXY;
	G4double shielding_thickness;	////
	G4double fib_radius;			////
	G4double fib_length;
	G4double fibXY;					////
	G4double fibMinZ, fibMaxZ;		////
	G4double FibPMT_gap;			////
	G4double air_gap;				//thickness of air gap between fiber and tile
	G4double tileAbsLength;
	G4double quartz_radius;
	G4double FibQuartz_gap;

	//Tile geometry parameters
	G4double angle1;
	G4double angle2;
	G4double Dx;
	G4double Dy;
	G4double Dz;
	G4double bendRadius;
	G4double distance;
	G4double fibRadius;
	G4int ieta;
	G4int layerNo;

	G4ThreeVector corners[8];		//stores the corner coordinates of the tile
	G4ThreeVector readout0, readout1;	//location of readout at corner 0/1
	G4ThreeVector mirror0, mirror1;	//location of mirrors at corner 0/1
	G4int readoutCorner;			//determines location of fiber readout (0:corner 0, 1:corner 1)

};

#endif

