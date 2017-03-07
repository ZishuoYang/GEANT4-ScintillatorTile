#include "LYSimTrajectory.hh"
#include "LYSimTrajectoryPoint.hh"

#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4AttDefStore.hh"

#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"

#include "G4Polyline.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"
#include "G4Polymarker.hh"

//#define G4ATTDEBUG
#ifdef G4ATTDEBUG
#include "G4AttCheck.hh"
#endif

G4Allocator<LYSimTrajectory> LYSimTrajectoryAllocator;

LYSimTrajectory::LYSimTrajectory()
    : fpPointsContainer(0), fTrackID(0), fParentID(0),
      PDGCharge(0.0), PDGEncoding(0), ParticleName(""),
      initialMomentum(G4ThreeVector())
{
    wls         = false;
    drawit      = false;
    forceNoDraw = false;
    forceDraw   = false;

    particleDefinition = NULL;
}

LYSimTrajectory::LYSimTrajectory(const G4Track* aTrack)
{
    particleDefinition = aTrack->GetDefinition();
    ParticleName = particleDefinition->GetParticleName();
    PDGCharge = particleDefinition->GetPDGCharge();
    PDGEncoding = particleDefinition->GetPDGEncoding();
    fTrackID = aTrack->GetTrackID();
    fParentID = aTrack->GetParentID();
    initialMomentum = aTrack->GetMomentum();
    fpPointsContainer = new LYSimTrajectoryPointContainer();
    // Following is for the first trajectory point
    fpPointsContainer->push_back(new LYSimTrajectoryPoint(aTrack));

    wls         = false;
    drawit      = false;
    forceNoDraw = false;
    forceDraw   = false;
}

LYSimTrajectory::LYSimTrajectory(LYSimTrajectory & right) : G4VTrajectory()
{
    particleDefinition=right.particleDefinition;
    ParticleName = right.ParticleName;
    PDGCharge = right.PDGCharge;
    PDGEncoding = right.PDGEncoding;
    fTrackID = right.fTrackID;
    fParentID = right.fParentID;
    initialMomentum = right.initialMomentum;
    fpPointsContainer = new LYSimTrajectoryPointContainer();

    for(size_t i=0;i<right.fpPointsContainer->size();++i) {
        LYSimTrajectoryPoint* rightPoint
            = (LYSimTrajectoryPoint*)((*(right.fpPointsContainer))[i]);
        fpPointsContainer->push_back(new LYSimTrajectoryPoint(*rightPoint));
    }

    wls         = right.wls;
    drawit      = right.drawit;
    forceNoDraw = right.forceNoDraw;
    forceDraw   = right.forceDraw;
}

LYSimTrajectory::~LYSimTrajectory()
{
    for(size_t i=0;i<fpPointsContainer->size();++i){
        delete  (*fpPointsContainer)[i];
    }
    fpPointsContainer->clear();

    delete fpPointsContainer;
}

void LYSimTrajectory::ShowTrajectory(std::ostream& os) const
{
    // Invoke the default implementation in G4VTrajectory...
    G4VTrajectory::ShowTrajectory(os);
    // ... or override with your own code here.
}

void LYSimTrajectory::DrawTrajectory() const
{
  // Invoke the default implementation in G4VTrajectory...
  G4VTrajectory::DrawTrajectory();
  // ... or override with your own code here.
}

void LYSimTrajectory::DrawTrajectory(G4int i_mode) const
{
    // Invoke the default implementation in G4VTrajectory...
    // G4VTrajectory::DrawTrajectory(i_mode);
    // ... or override with your own code here.

    //Taken from G4VTrajectory and modified to select colours based on particle
    //type and to selectively eliminate drawing of certain trajectories.

    if (!forceDraw && (!drawit || forceNoDraw)) return;

    // If i_mode>=0, draws a trajectory as a polyline and, if i_mode!=0,
    // adds markers - yellow circles for step points and magenta squares
    // for auxiliary points, if any - whose screen size in pixels is
    // given by std::abs(i_mode)/1000.  E.g: i_mode = 5000 gives easily
    // visible markers.

    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if (!pVVisManager) return;

    const G4double markerSize = std::abs(i_mode)/1000;
    G4bool lineRequired (i_mode >= 0);
    G4bool markersRequired (markerSize > 0.);

    G4Polyline trajectoryLine;
    G4Polymarker stepPoints;
    G4Polymarker auxiliaryPoints;

    for (G4int i = 0; i < GetPointEntries() ; i++) {
      G4VTrajectoryPoint* aTrajectoryPoint = GetPoint(i);
      const std::vector<G4ThreeVector>* auxiliaries
        = aTrajectoryPoint->GetAuxiliaryPoints();
      if (auxiliaries) {
        for (size_t iAux = 0; iAux < auxiliaries->size(); ++iAux) {
          const G4ThreeVector pos((*auxiliaries)[iAux]);
          if (lineRequired) {
            trajectoryLine.push_back(pos);
          }
          if (markersRequired) {
            auxiliaryPoints.push_back(pos);
          }
        }
      }
      const G4ThreeVector pos(aTrajectoryPoint->GetPosition());
      if (lineRequired) {
        trajectoryLine.push_back(pos);
      }
      if (markersRequired) {
        stepPoints.push_back(pos);
      }
    }

    if (lineRequired) {
      G4Colour colour;

      if(particleDefinition==G4OpticalPhoton::OpticalPhotonDefinition()){
        if(wls) //WLS photons are red
          colour = G4Colour(1.,0.,0.);
        else{ //Scintillation and Cerenkov photons are green
          colour = G4Colour(0.,1.,0.);
        }
      }
      else //All other particles are blue
        colour = G4Colour(0.,0.,1.);

      G4VisAttributes trajectoryLineAttribs(colour);
      trajectoryLine.SetVisAttributes(&trajectoryLineAttribs);
      pVVisManager->Draw(trajectoryLine);
    }
    if (markersRequired) {
      auxiliaryPoints.SetMarkerType(G4Polymarker::squares);
      auxiliaryPoints.SetScreenSize(markerSize);
      auxiliaryPoints.SetFillStyle(G4VMarker::filled);
      G4VisAttributes auxiliaryPointsAttribs(G4Colour(0.,1.,1.));  // Magenta
      auxiliaryPoints.SetVisAttributes(&auxiliaryPointsAttribs);
      pVVisManager->Draw(auxiliaryPoints);

      stepPoints.SetMarkerType(G4Polymarker::circles);
      stepPoints.SetScreenSize(markerSize);
      stepPoints.SetFillStyle(G4VMarker::filled);
      G4VisAttributes stepPointsAttribs(G4Colour(1.,1.,0.));  // Yellow.
      stepPoints.SetVisAttributes(&stepPointsAttribs);
      pVVisManager->Draw(stepPoints);
    }
}

void LYSimTrajectory::AppendStep(const G4Step* aStep)
{
    fpPointsContainer->push_back(new LYSimTrajectoryPoint(aStep));
}

G4ParticleDefinition* LYSimTrajectory::GetParticleDefinition()
{
    return (G4ParticleTable::GetParticleTable()->FindParticle(ParticleName));
}

void LYSimTrajectory::MergeTrajectory(G4VTrajectory* secondTrajectory)
{
    if(!secondTrajectory) return;

    LYSimTrajectory* second = (LYSimTrajectory*)secondTrajectory;
    G4int ent = second->GetPointEntries();
    // initial point of the second trajectory should not be merged
    for(G4int i=1; i<ent; ++i) {
        fpPointsContainer->push_back((*(second->fpPointsContainer))[i]);
    }
    delete (*second->fpPointsContainer)[0];
    second->fpPointsContainer->clear();
}

const std::map<G4String,G4AttDef>* LYSimTrajectory::GetAttDefs() const
{
    G4bool isNew;
    std::map<G4String,G4AttDef>* store
        = G4AttDefStore::GetInstance("Trajectory",isNew);

    if (isNew) {

      G4String ID("ID");
      (*store)[ID] = G4AttDef(ID,"Track ID","Bookkeeping","","G4int");

      G4String PID("PID");
      (*store)[PID] = G4AttDef(PID,"Parent ID","Bookkeeping","","G4int");

      G4String PN("PN");
      (*store)[PN] = G4AttDef(PN,"Particle Name","Physics","","G4String");

      G4String Ch("Ch");
      (*store)[Ch] = G4AttDef(Ch,"Charge","Physics","e+","G4double");

      G4String PDG("PDG");
      (*store)[PDG] = G4AttDef(PDG,"PDG Encoding","Physics","","G4int");

      G4String IMom("IMom");
      (*store)[IMom] = G4AttDef(IMom,
                       "Momentum of track at start of trajectory",
                       "Physics","G4BestUnit","G4ThreeVector");

      G4String IMag("IMag");
      (*store)[IMag] = G4AttDef(IMag,
                       "Magnitude of momentum of track at start of trajectory",
                       "Physics","G4BestUnit","G4double");

        G4String NTP("NTP");
        (*store)[NTP] = G4AttDef(NTP,"No. of points","Bookkeeping","","G4int");

    }
    return store;
}

std::vector<G4AttValue>* LYSimTrajectory::CreateAttValues() const
{
  std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

  values->push_back
    (G4AttValue("ID",G4UIcommand::ConvertToString(fTrackID),""));

  values->push_back
    (G4AttValue("PID",G4UIcommand::ConvertToString(fParentID),""));

  values->push_back(G4AttValue("PN",ParticleName,""));

  values->push_back
    (G4AttValue("Ch",G4UIcommand::ConvertToString(PDGCharge),""));

  values->push_back
    (G4AttValue("PDG",G4UIcommand::ConvertToString(PDGEncoding),""));

  values->push_back
    (G4AttValue("IMom",G4BestUnit(initialMomentum,"Energy"),""));

  values->push_back
    (G4AttValue("IMag",G4BestUnit(initialMomentum.mag(),"Energy"),""));

  values->push_back
    (G4AttValue("NTP",G4UIcommand::ConvertToString(GetPointEntries()),""));

#ifdef G4ATTDEBUG
  G4cout << G4AttCheck(values,GetAttDefs());
#endif
    return values;
}
