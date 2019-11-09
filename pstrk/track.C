/// \file
/// \ingroup tutorial_eve
/// Demonstrates usage of TEveTrackPRopagator with different magnetic
/// field configurations.
/// Needs to be run in compiled mode.
/// root
/// ~~~{.cpp}
///   .L track.C+
///   track(3, kTRUE)
/// ~~~
/// void track(Int_t mode = 5, Bool_t isRungeKutta = kTRUE)
/// Modes are
///  0. B = 0, no difference btween signed and charge particles;
///  1. constant B field (along z, but could have arbitrary direction);
///  2. variable B field, sign change at  R = 200 cm;
///  3. magnetic field with a zero-field region;
///  4. CMS magnetic field - simple track;
///  5. CMS magnetic field - track with different path-marks.
///  6. Concpetual ILC detector, problematic track
///
/// \image html eve_track.png
/// \macro_code
///
/// \author Alja Mrak-Tadel


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

TEveTrackPropagator* g_prop = 0;

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
    ifstream f("BFieldMaps/Mau13/PSMap.txt");
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
        if (nLine%100000 == 0) cout << "line " << nLine << endl;
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
    cout << "pos(" << pos.x() << ", " << pos.y() << ", " << pos.z() << ") B(" << b.x() << ", " << b.y() << ", " << b.z() << ")" << endl;

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

//==============================================================================

//______________________________________________________________________________
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


void track(Int_t mode = 4, Bool_t isRungeKutta = kTRUE)
{
   //readMap();
   //TVector3 pos(2804, -1200, -9891.5);
   //TVector3 b = getBField( pos );

   gSystem->IgnoreSignal(kSigSegmentationViolation, true);
   TEveManager::Create();

  TGeoManager::Import("mu2e_gdml.root");
  TGeoVolume* volHallAir  = gGeoManager->GetVolume("HallAir");
  cout << volHallAir->GetName() << endl;
  TGeoVolume* volPSVacuum = gGeoManager->GetVolume("PSVacuum");
  cout << volPSVacuum->GetName() << endl;

  bool kDrawDet = 1;
  TVector3 posTarget(0,0,0);
  if (kDrawDet) {
    // Draw detector, draw PSVacuum only
    vector<TGeoNode*> vecNodePSVacuum = findNode(volHallAir, TString("PSVacuum"));
    cout << "vecNodePSVacuum size " << vecNodePSVacuum.size() << endl;
    TGeoNode* node0 = vecNodePSVacuum[0];
    if (!node0) cout << "node0 not found" << endl;
    cout << "node0 " << node0->GetName() << endl;

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
    //cout << "Translation: (" << trans[0] << ", " << trans[1] << ", " << trans[2] << ")" << endl;
    posTarget = TVector3(trans[0], trans[1], trans[2]);

    TEveGeoTopNode* eveNode0 = new TEveGeoTopNode(gGeoManager, node0);
    gEve->AddGlobalElement(eveNode0);
  }

   TEveTrackList *list = new TEveTrackList();
   TEveTrackPropagator* prop = g_prop = list->GetPropagator();
   prop->SetFitDaughters(kFALSE);
   prop->SetMaxR(10000);
   prop->SetMaxZ(20000);
   prop->SetMaxOrbs(1.0); // circles

   if (isRungeKutta)
   {
      prop->SetStepper(TEveTrackPropagator::kRungeKutta);
      list->SetName("RK Propagator");
   }
   else
   {
      list->SetName("Heix Propagator");
   }

   TEveTrack *track = 0;
   switch (mode)
   {
      case 0:
      {
         // B = 0 no difference btween signed and charge particles
         prop->SetMagField(0.);
         list->SetElementName(Form("%s, zeroB", list->GetElementName()));
         track = make_track(prop, 1);
         break;
      }

      case 1:
      {
         // constant B field, const angle
         prop->SetMagFieldObj(new TEveMagFieldConst(0., 0., 1.0));
         list->SetElementName(Form("%s, constB", list->GetElementName()));
         track = make_track(prop, 1);
         break;
      }

      case 2:
      {
         // variable B field, sign change at  R = 200 cm
         prop->SetMagFieldObj(new TEveMagFieldDuo(200, -4.4, 2));
         list->SetElementName(Form("%s, duoB", list->GetElementName()));
         track = make_track(prop, 1);
         break;
      }

      case 3:
      {
         // gapped field
         prop->SetMagFieldObj(new GappedField());
         list->SetElementName(Form("%s, gappedB", list->GetElementName()));


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
         prop->SetMaxR(500);
         prop->SetMaxZ(300);
         prop->SetMaxOrbs(20.0); // circles
         prop->SetRnrDaughters(kTRUE);
         prop->SetRnrDecay(kTRUE);
         prop->RefPMAtt().SetMarkerStyle(4);
         list->SetElementName(Form("%s, Mu2e field", list->GetElementName()));


         TEveRecTrackD *rc = new TEveRecTrackD();
         //rc->fV.Set(390.4, 0.0, -616.45);
         rc->fV.Set(0, 0, 30);
         rc->fP.Set(-0.17, 0.0, -0.05);
         rc->fSign = -1;
         track = new TEveTrack(rc, prop);

         track->SetRnrPoints(kTRUE);
         track->SetMarkerStyle(4);

         break;
      }

      case 5:
      {
         // Magnetic field of CMS I.
         CmsMagField* mf = new CmsMagField;
         mf->setReverseState(true);
         mf->setSimpleModel(false);

         prop->SetMagFieldObj(mf);
         prop->SetMaxR(1000);
         prop->SetMaxZ(1000);
         prop->SetRnrReferences(kTRUE);
         prop->SetRnrDaughters(kTRUE);
         prop->SetRnrDecay(kTRUE);
         prop->RefPMAtt().SetMarkerStyle(4);
         list->SetElementName(Form("%s, CMS field", list->GetElementName()));

         TEveRecTrackD *rc = new TEveRecTrackD();
         rc->fV.Set(-16.426592, 16.403185, -19.782692);
         rc->fP.Set(3.631100, 3.643450, 0.682254);
         rc->fSign = -1;
         track = new TEveTrack(rc, prop);

         track->AddPathMark(TEvePathMarkD(TEvePathMarkD::kReference,
                  TEveVectorD(-1.642659e+01, 1.640318e+01, -1.978269e+01),
                  TEveVectorD(3.631100, 3.643450, 0.682254)));
         track->AddPathMark(TEvePathMarkD(TEvePathMarkD::kReference,
                  TEveVectorD(-1.859987e+00, 3.172243e+01, -1.697866e+01),
                  TEveVectorD(3.456056, 3.809894, 0.682254)));
         track->AddPathMark(TEvePathMarkD(TEvePathMarkD::kReference,
                  TEveVectorD(4.847579e+01, 9.871711e+01, -5.835719e+00),
                  TEveVectorD(2.711614, 4.409945, 0.687656)));
         track->AddPathMark(TEvePathMarkD(TEvePathMarkD::kDaughter,
                  TEveVectorD(1.342045e+02, 4.203950e+02, 3.846268e+01)));
         track->AddPathMark(TEvePathMarkD(TEvePathMarkD::kDaughter,
                  TEveVectorD(1.483827e+02, 5.124750e+02, 5.064311e+01)));
         track->AddPathMark(TEvePathMarkD(TEvePathMarkD::kDaughter,
                  TEveVectorD(1.674676e+02, 6.167731e+02, 6.517403e+01)));
         track->AddPathMark(TEvePathMarkD(TEvePathMarkD::kDecay,
                  TEveVectorD(1.884976e+02, 7.202000e+02, 7.919290e+01)));

         track->SetRnrPoints(kTRUE);
         track->SetMarkerStyle(4);

         break;
      }

      case 6:
      {
         // Problematic track from Druid
         prop->SetMagFieldObj(new TEveMagFieldDuo(350, -3.5, 2.0));
         prop->SetMaxR(1000);
         prop->SetMaxZ(1000);
         prop->SetRnrReferences(kTRUE);
         prop->SetRnrDaughters(kTRUE);
         prop->SetRnrDecay(kTRUE);
         prop->RefPMAtt().SetMarkerStyle(4);
         list->SetElementName(Form("%s, Some ILC Detector field",
                                   list->GetElementName()));

         TEveRecTrackD *rc = new TEveRecTrackD();
         rc->fV.Set(57.1068, 31.2401, -7.07629);
         rc->fP.Set(4.82895, 2.35083, -0.611757);
         rc->fSign = 1;
         track = new TEveTrack(rc, prop);

         track->AddPathMark(TEvePathMarkD(TEvePathMarkD::kDaughter,
                  TEveVectorD(1.692235e+02, 7.047929e+01, -2.064785e+01)));
         track->AddPathMark(TEvePathMarkD(TEvePathMarkD::kDaughter,
                  TEveVectorD(5.806180e+02, 6.990633e+01, -6.450000e+01)));
         track->AddPathMark(TEvePathMarkD(TEvePathMarkD::kDecay,
                  TEveVectorD(6.527213e+02, 1.473249e+02, -8.348498e+01)));

         track->SetRnrPoints(kTRUE);
         track->SetMarkerStyle(4);

         break;
      }
   };

   if (isRungeKutta)
      list->SetLineColor(kMagenta);
   else
      list->SetLineColor(kCyan);

   track->SetLineColor(list->GetLineColor());

   gEve->AddElement(list);
   list->AddElement(track);

   prop->PrintMagField(1,1,1);
   cout << "nPoints " << prop->GetCurrentPoint() << endl;

   track->MakeTrack();

   TEveViewer *ev = gEve->GetDefaultViewer();
   TGLViewer  *gv = ev->GetGLViewer();
   gv->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, 0);

   gEve->Redraw3D(kTRUE);
   gSystem->ProcessEvents();

   //gv->CurrentCamera().RotateRad(-0.5, 1.4);
   gv->RequestDraw();
}

