#include "RecoilScintHit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "SimulationManager.h"


G4Allocator<RecoilScintHit> RecoilScintHitAllocator;

RecoilScintHit::RecoilScintHit()
{;}
//______________________________________________________________________________

RecoilScintHit::~RecoilScintHit()
{;}
//______________________________________________________________________________

G4int RecoilScintHit::operator==(const RecoilScintHit &right) const
{
  return (this==&right) ? 1 : 0;
}
//______________________________________________________________________________

void RecoilScintHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(fPosition);
    circle.SetScreenSize(10.04);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(0.,0.,1.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}
//______________________________________________________________________________

void RecoilScintHit::Print()
{;}
//______________________________________________________________________________


