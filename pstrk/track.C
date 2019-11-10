/// Toy macro to simulate and display particles flying in Mu2e PS
/// Modified from tutorials/eve/track.C
/// Needs to be run in compiled mode.
/// root
///   .L track.C+
///   track();
///   newTrack(px, py, pz, charge, refresh);  // newTrack(0.0, 0.18, -0.025, 1, 1)
///
/// \author Zhengyun You

#include "TEveTrackPropagator.h"
#include "TEveTrack.h"
#include "TEveVSDStructs.h"
#include "TEveManager.h"
#include "TEveViewer.h"
#include "TSystem.h"
#include "TGLViewer.h"
#include "TMath.h"

#include "TEveViewer.h"
#include "TEvePointSet.h"

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoNode.h"
#include "TEveGeoNode.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

const double kX0=2804.0, kY0=-1200.0, kZ0=-9929.0; // mm
const int knX=89, knY=97, knZ=281;
const double kdX=25.0, kdY=25.0, kdZ=25.0; // mm
const double kX1=kX0+(knX-1)*kdX, kY1=kY0+(knY-1)*kdY, kZ1=kZ0+(knZ-1)*kdZ;
TVector3 kCenter(3909.0, 0.0, -6164.5); // mm

TEveTrackList *g_list = 0;
TEveTrackPropagator* g_prop = 0;
Bool_t isRungeKutta = kTRUE;

class GappedField : public TEveMagField
{
public:
   GappedField():TEveMagField(){}
   ~GappedField(){};
   using   TEveMagField::GetField;

   virtual TEveVectorD GetFieldD(Double_t /*x*/, Double_t /*y*/, Double_t z) const
   {
      if (TMath::Abs(z) < 300) return TEveVectorD(0, 0, -4);
      if (TMath::Abs(z) < 600) return TEveVectorD(0, 0, 0);
      return TEveVectorD(0, 0, 4);
   }
};

//==============================================================================

class CmsMagField: public TEveMagField
{
   bool m_magnetIsOn;
   bool m_reverse;
   bool m_simpleModel;

public:
   CmsMagField():
      m_magnetIsOn(true),
      m_reverse(false),
      m_simpleModel(true){}

   virtual ~CmsMagField(){}
   virtual Double_t   GetMaxFieldMagD() const { return m_magnetIsOn ? 3.8 : 0.0; }
   void               setMagnetState( bool state )
   {
      if (state != m_magnetIsOn)
      {
         if ( state )
            std::cout << "Magnet state is changed to ON" << std::endl;
         else
            std::cout << "Magnet state is changed to OFF" << std::endl;
      }
      m_magnetIsOn = state;
   }

   bool isMagnetOn() const               { return m_magnetIsOn;}
   void setReverseState(bool state)      { m_reverse = state; }
   bool isReverse() const                { return m_reverse;}
   void setSimpleModel(bool simpleModel) { m_simpleModel = simpleModel; }
   bool isSimpleModel() const            { return m_simpleModel;}

   using   TEveMagField::GetField;

   virtual TEveVectorD GetFieldD(Double_t x, Double_t y, Double_t z) const
   {
      double R = sqrt(x*x+y*y);
      double field = m_reverse?-GetMaxFieldMag():GetMaxFieldMag();
      //barrel
      if ( TMath::Abs(z)<724 )
      {
         //inside solenoid
         if ( R < 300) return TEveVectorD(0,0,field);
         // outside solinoid
         if ( m_simpleModel ||
              ( R>461.0 && R<490.5 ) ||
              ( R>534.5 && R<597.5 ) ||
              ( R>637.0 && R<700.0 ) )
            return TEveVectorD(0,0,-field/3.8*1.2);

      } else {
         // endcaps
         if (m_simpleModel)
         {
            if ( R < 50 ) return TEveVectorD(0,0,field);
            if ( z > 0 )
               return TEveVectorD(x/R*field/3.8*2.0, y/R*field/3.8*2.0, 0);
            else
               return TEveVectorD(-x/R*field/3.8*2.0, -y/R*field/3.8*2.0, 0);
         }
         // proper model
         if ( ( TMath::Abs(z)>724 && TMath::Abs(z)<786  ) ||
              ( TMath::Abs(z)>850 && TMath::Abs(z)<910  ) ||
              ( TMath::Abs(z)>975 && TMath::Abs(z)<1003 ) )
         {
            if ( z > 0 )
               return TEveVectorD(x/R*field/3.8*2.0, y/R*field/3.8*2.0, 0);
            else
               return TEveVectorD(-x/R*field/3.8*2.0, -y/R*field/3.8*2.0, 0);
         }
      }
      return TEveVectorD(0,0,0);
   }
};

class Mu2eMagField: public TEveMagField
{
   bool m_magnetIsOn;
   bool m_reverse;
   bool m_simpleModel;

   vector<TVector3> gbPos;
   vector<TVector3> gbField;

public:
   Mu2eMagField():
      m_magnetIsOn(true),
      m_reverse(false),
      m_simpleModel(true){}

   virtual ~Mu2eMagField(){}
   virtual Double_t   GetMaxFieldMagD() const { return m_magnetIsOn ? 3.8 : 0.0; }
   void               setMagnetState( bool state )
   {
      if (state != m_magnetIsOn)
      {
         if ( state )
            std::cout << "Magnet state is changed to ON" << std::endl;
         else
            std::cout << "Magnet state is changed to OFF" << std::endl;
      }
      m_magnetIsOn = state;
   }

   bool isMagnetOn() const               { return m_magnetIsOn;}
   void setReverseState(bool state)      { m_reverse = state; }
   bool isReverse() const                { return m_reverse;}
   void setSimpleModel(bool simpleModel) { m_simpleModel = simpleModel; }
   bool isSimpleModel() const            { return m_simpleModel;}
 
   void readMap();
   TVector3 getBField(TVector3 pos) const;

   using   TEveMagField::GetField;

   virtual TEveVectorD GetFieldD(Double_t x, Double_t y, Double_t z) const
   {
      double R = sqrt(x*x+y*y);
      double field = m_reverse?-GetMaxFieldMag():GetMaxFieldMag();

      TVector3 v=getBField(TVector3(x*10.0,y*10.0,z*10.0));
      if (m_reverse) v *= -1.0;

      return TEveVectorD(v.x(), v.y(), v.z());
   }
};

void Mu2eMagField::readMap()
{
    cout << "Initialize PS BFieldMap..." << endl;
    ifstream f("../BFieldMaps/Mau13/PSMap.txt");
    string line;
    int nLine = 0;
    double x, y, z;
    double bx, by, bz;

    while (!f.eof()) { // && nLine < 100000) {
        getline(f, line);
        istringstream iss(line);

        if (iss >> x >> y >> z >> bx >> by >> bz) {
            TVector3 pos(x, y, z); 
            pos = pos - kCenter; 
            TVector3 b(bx, by, bz);
            //cout << x << " " << y << " " << z << endl;
            gbPos.push_back(pos);
            gbField.push_back(b);
        }
        nLine++;
        if (nLine%100000 == 0) cout << "Read line " << nLine << "..." << endl;
    }
}

TVector3 Mu2eMagField::getBField(TVector3 pos) const
{
    //pos.Print();
    // in range judgement
    if (pos.x() < kX0-kCenter.x() || pos.x() > kX1-kCenter.x() || 
        pos.y() < kY0-kCenter.y() || pos.y() > kY1-kCenter.y() ||
        pos.z() < kZ0-kCenter.z() || pos.z() > kZ1-kCenter.z())
       return TVector3(0,0,0);
    //
    TVector3 bP[8];
    TVector3 bF[8];
    int ix0 = int((pos.x()-(kX0-kCenter.x()))/kdX);
    int iy0 = int((pos.y()-(kY0-kCenter.y()))/kdY);
    int iz0 = int((pos.z()-(kZ0-kCenter.z()))/kdZ);
    //cout << "ix0 " << ix0 << " iy0 " << iy0 << " iz0 " << iz0 << endl;

    int iEntry = (ix0*knY+iy0)*knZ+iz0;
    //cout << "iEntry " << iEntry << endl;
    //cout << "gpPos size " << gbPos.size() << endl;

    bP[0] = gbPos[iEntry];       // -x,-y,-z
    //bP[0].Print();
    bP[1] = gbPos[iEntry+1];     // -x,-y, z
    bP[2] = gbPos[iEntry+knZ];   // -x, y,-z
    bP[3] = gbPos[iEntry+knZ+1]; // -x, y, z
    bP[4] = gbPos[iEntry+knY*knZ];   //  x,-y,-z
    bP[5] = gbPos[iEntry+knY*knZ+1]; //  x,-y, z
    bP[6] = gbPos[iEntry+(knY+1)*knZ];   //  x, y,-z
    bP[7] = gbPos[iEntry+(knY+1)*knZ+1]; //  x, y,-z

    bF[0] = gbField[iEntry];       // -x,-y,-z
    bF[1] = gbField[iEntry+1];     // -x,-y, z
    bF[2] = gbField[iEntry+knZ];   // -x, y,-z
    bF[3] = gbField[iEntry+knZ+1]; // -x, y, z
    bF[4] = gbField[iEntry+knY*knZ];   //  x,-y,-z
    bF[5] = gbField[iEntry+knY*knZ+1]; //  x,-y, z
    bF[6] = gbField[iEntry+(knY+1)*knZ];   //  x, y,-z
    bF[7] = gbField[iEntry+(knY+1)*knZ+1]; //  x, y,-z

    TVector3 b(0, 0, 0);
    double xd = (pos.x()-bP[0].x())/kdX;
    double yd = (pos.y()-bP[0].y())/kdY;
    double zd = (pos.z()-bP[0].z())/kdZ;

    for (int i = 0; i < 8; i++) {
        //cout << bP[i].x() << " " << bP[i].y() << " " << bP[i].z() << " "  << bF[i].x() << " " << bF[i].y() << " " << bF[i].z()  << endl;
        double dist = (pos-bP[i]).Mag();
        if (dist < 1e-3) return bF[i];

        double wt = 0;
        if (i==0) wt = (1-xd)*(1-yd)*(1-zd);
        else if (i==1) wt = (1-xd)*(1-yd)*zd;
        else if (i==2) wt = (1-xd)*yd*(1-zd);
        else if (i==3) wt = (1-xd)*yd*zd;
        else if (i==4) wt = xd*(1-yd)*(1-zd);
        else if (i==5) wt = xd*(1-yd)*zd;
        else if (i==6) wt = xd*yd*(1-zd);
        else if (i==7) wt = xd*yd*zd;

        TVector3 bWeight = bF[i]*wt;
        b += bWeight;
    }  
    //b.Print();
    //cout << "pos(" << pos.x() << ", " << pos.y() << ", " << pos.z() << ") B(" << b.x() << ", " << b.y() << ", " << b.z() << ")" << endl;

    return b;
}

//==============================================================================
vector<TGeoNode*> findNode(TGeoVolume* vol, TString name)
{
  vector<TGeoNode*> vec;
  int nChild = vol->GetNdaughters();
  for (int i = 0; i < nChild; i++) {
    TGeoNode* node = vol->GetNode(i);
    if (TString(node->GetName()).Contains(name)) vec.push_back(node);
  }

  return vec;
}

void initDet()
{
  TGeoManager::Import("mu2e_gdml_root5.root");
  TGeoVolume* volHallAir  = gGeoManager->GetVolume("HallAir");
  cout << volHallAir->GetName() << endl;
  TGeoVolume* volPSVacuum = gGeoManager->GetVolume("PSVacuum");
  cout << volPSVacuum->GetName() << endl;
  //volPSVacuum->GetShape()->Inspect();  // Tube:R=75.375,dZ=249.648

  bool kDrawDet = 1;
  TVector3 posTarget(0,0,0);
  if (kDrawDet) {
    // Draw detector, draw PSVacuum only
    vector<TGeoNode*> vecNodePSVacuum = findNode(volHallAir, TString("PSVacuum"));
    cout << "vecNodePSVacuum size " << vecNodePSVacuum.size() << endl;
    TGeoNode* node0 = vecNodePSVacuum[0];
    if (!node0) cout << "node0 not found" << endl;
    cout << "node0 " << node0->GetName() << endl;
    //cout << "Translation in its mother:" << endl;
    //node0->GetMatrix()->Print();

    // Find the PSShield and set it invisible or transparent to see inside
    TGeoVolume* volPSVacuum = node0->GetVolume();
    vector<TGeoNode*> vecNodePSShieldShell = findNode(volPSVacuum, TString("PSShieldShell"));
    cout << "vecNodePSShieldShell size " << vecNodePSShieldShell.size() << endl;
    for (int iNode = 0; iNode < int(vecNodePSShieldShell.size()); iNode++) {
      TGeoNode* nodeShell = vecNodePSShieldShell[iNode];
      //nodeShell->SetVisibility(0);
      nodeShell->GetVolume()->SetTransparency(92);
    }

    // Coordinate origin is the center of PSVacuum, get the relative positon of target center
    TGeoNode* nodeTarget = node0->GetVolume()->GetNode(0);
    cout << nodeTarget->GetName() << endl;
    TGeoMatrix* matrixTarget = nodeTarget->GetMatrix();
    matrixTarget->Print();
    const double *trans = matrixTarget->GetTranslation();
    cout << "PSTarget in PSVacuum Translation: (" << trans[0] << ", " << trans[1] << ", " << trans[2] << ")" << endl; // (0, 0, 33.5975) 
    posTarget = TVector3(trans[0], trans[1], trans[2]);

    // The center of node PSVacuum is not PSTarget center, add a TOP as mother to shift world coordinate center to PSTarget center
    TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
    TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
    TGeoVolume *top = gGeoManager->MakeBox("TOP", Vacuum, 80.0, 80.0, 300.0);
    TGeoTranslation *trPSVacuumInMother = new TGeoTranslation(-1.0*trans[0], -1.0*trans[1], -1.0*trans[2]);
    top->AddNode(volPSVacuum, 1, trPSVacuumInMother);
    TGeoNodeOffset *topNode = new TGeoNodeOffset(top, 1, 0.0);

    // Add the top node to eve
    TEveGeoTopNode* eveTopNode = new TEveGeoTopNode(gGeoManager, topNode); // node0
    gEve->AddGlobalElement(eveTopNode);
  }
}

//==============================================================================
TEveTrack* make_track(TEveTrackPropagator* prop, Int_t sign)
{
  // Make track with given propagator.
  // Add to math-marks to test fit.

  TEveRecTrackD *rc = new TEveRecTrackD();
  //rc->fV.Set(0.028558, -0.000918, 3.691919);
  //rc->fP.Set(0.767095, -2.400006, -0.313103);
  rc->fV.Set(0.0, 0.0, 0.0);
  rc->fP.Set(1, 0, 0.1);
  rc->fSign = sign;

  TEveTrack* track = new TEveTrack(rc, prop);
  track->SetName(Form("Charge %d", sign));
/*  // daughter 0
  TEvePathMarkD* pm1 = new TEvePathMarkD(TEvePathMarkD::kDaughter);
  pm1->fV.Set(1.479084, -4.370661, 3.119761);
  track->AddPathMark(*pm1);
  // daughter 1
  TEvePathMarkD* pm2 = new TEvePathMarkD(TEvePathMarkD::kDaughter);
  pm2->fV.Set(57.72345, -89.77011, -9.783746);
  track->AddPathMark(*pm2);
*/
  return track;
}

void newTrack(double px=0.0, double py=0.18, double pz=-0.025, int charge=1, int refresh=0)
{
  if (refresh) g_list->RemoveElements();

  TEveRecTrackD *rc = new TEveRecTrackD();
  rc->fV.Set(0, 0, 0);
  rc->fP.Set(px, py, pz);
  rc->fSign = -1*charge;
  TEveTrack *track = new TEveTrack(rc, g_prop);

  track->SetRnrPoints(kTRUE);
  track->SetMarkerStyle(4);

  if (isRungeKutta)
     g_list->SetLineColor(kMagenta);
  else
     g_list->SetLineColor(kCyan);
  track->SetLineColor(g_list->GetLineColor());

  g_list->AddElement(track);
  track->MakeTrack();

  TEveViewer *ev = gEve->GetDefaultViewer();
  TGLViewer  *gv = ev->GetGLViewer();
  gv->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, 0);

  gEve->Redraw3D(kTRUE);
  gSystem->ProcessEvents();

  //gv->CurrentCamera().RotateRad(-0.5, 1.4);
  gv->RequestDraw();
}

void track(Int_t mode = 4)
{
  gSystem->IgnoreSignal(kSigSegmentationViolation, true);
  TEveManager::Create();
  initDet();

  g_list = new TEveTrackList();
  gEve->AddElement(g_list);
  TEveTrackPropagator* prop = g_prop = g_list->GetPropagator();
  prop->SetFitDaughters(kFALSE);
  prop->SetMaxR(10000);
  prop->SetMaxZ(20000);
  prop->SetMaxOrbs(1.0); // circles

  if (isRungeKutta)
  {
    prop->SetStepper(TEveTrackPropagator::kRungeKutta);
    g_list->SetName("RK Propagator");
  }
  else
  {
    g_list->SetName("Heix Propagator");
  }

   TEveTrack *track = 0;
   switch (mode)
   {
      case 0:
      {
         // B = 0 no difference btween signed and charge particles
         prop->SetMagField(0.);
         g_list->SetElementName(Form("%s, zeroB", g_list->GetElementName()));
         track = make_track(prop, 1);
         break;
      }
      case 1:
      {
         // constant B field, const angle
         prop->SetMagFieldObj(new TEveMagFieldConst(0., 0., 1.0));
         g_list->SetElementName(Form("%s, constB", g_list->GetElementName()));
         track = make_track(prop, 1);
         break;
      }
      case 2:
      {
         // variable B field, sign change at  R = 200 cm
         prop->SetMagFieldObj(new TEveMagFieldDuo(200, -4.4, 2));
         g_list->SetElementName(Form("%s, duoB", g_list->GetElementName()));
         track = make_track(prop, 1);
         break;
      }
      case 3:
      {
         // gapped field
         prop->SetMagFieldObj(new GappedField());
         g_list->SetElementName(Form("%s, gappedB", g_list->GetElementName()));

         TEveRecTrackD *rc = new TEveRecTrackD();
         rc->fV.Set(0.028558, -0.000918, 3.691919);
         rc->fP.Set(0.767095, -0.400006, 2.313103);
         rc->fSign = 1;
         track = new TEveTrack(rc, prop);

         TEvePointSet* marker = new TEvePointSet(2);
         marker->SetElementName("B field break points");
         marker->SetPoint(0, 0., 0., 300.f);
         marker->SetPoint(1, 0., 0., 600.f);
         marker->SetMarkerColor(3);
         gEve->AddElement(marker);
         break;
      }
      case 4:
      {
         // Magnetic field of Mu2e.
         Mu2eMagField* mf = new Mu2eMagField;
         mf->setReverseState(false);
         mf->readMap();
         //TVector3 pos(2804, -1200, -9891.5);
         //TVector3 pos(3904, 0, -6164.5);
         TVector3 pos(0,0,0);
         TVector3 b = mf->getBField( pos );

         prop->SetMagFieldObj(mf);
         prop->SetMaxR(100);
         prop->SetMaxZ(300);
         prop->SetMaxOrbs(20.0); // circles
         prop->SetRnrDaughters(kTRUE);
         prop->SetRnrDecay(kTRUE);
         prop->RefPMAtt().SetMarkerStyle(4);
         g_list->SetElementName(Form("%s, Mu2e field", g_list->GetElementName()));

         break;
      }
   };

   newTrack();
}

