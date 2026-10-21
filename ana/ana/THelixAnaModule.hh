///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef Stntuple_ana_THelixAnaModule_hh
#define Stntuple_ana_THelixAnaModule_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TStnTimeClusterBlock.hh"
#include "Stntuple/obj/TStnHelixBlock.hh"
#include "Stntuple/obj/TStnTrackSeedBlock.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "Stntuple/obj/TStrawHitBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"
#include "Stntuple/obj/TStepPointMCBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/geom/TStnCrystal.hh"
#include "Stntuple/geom/TDiskCalorimeter.hh"

namespace  stntuple {
class THelixAnaModule: public TStnModule {
public:

  struct TrackPar_t {
    int     fNHPl;
    int     fNEPl;
    int     fNDPl;
    float   fDpF ;    // tracker-only resolution
    float   fDp0 ;
    float   fDp2 ;
    float   fDpFSt;
    double  fDioWt;

    double  fXc;      // X,Y coordinates of the track helix center
    double  fYc;
    double  fPhic;      // phi angle defined by XXc and Yx

    double  fEcl;
    double  fEp;
    double  fDx;
    double  fDy;
    double  fDz;
    double  fDt;
    double  fDu;			// rotated residuals
    double  fDv;
    double  fChi2Match;
    double  fPath;
  };

  struct EventHist_t {
    TH1F*      fRunNumber;
    TH1F*      fEventNumber;
    TH1F*      fSectionNumber; //subrun
    TH1F*    fRv;			// MC truth information
    TH1F*    fZv;
    TH1F*    fEleMom;
    TH1F*    fEleCosTh;
    TH1F*    fNTimeClusters[2];
    TH1F*    fNHelices[10];
    TH1F*    fNTracks[10];
    TH1F*    fNHelPos;		       // N(helices) with helicity > 0
    TH1F*    fNHelNeg;		       // N(helices) with helicity < 0
    TH1F*    fNGenp;                   // N(particles in GENP block)
  };

  struct TimeClusterHist_t {
    TH1F*    fNHits;
    TH1F*    fNComboHits;
    TH1F*    fT0;
    TH1F*    fT0Zoom;        // 500-900 ns
    TH1F*    fClusterTime;   // calorimeter cluster time, 500-900 ns
    TH1F*    fClusterEnergy;
  
  };


  struct HelixHist_t {
    TH1F*    fNHits;
    TH1F*    fHelicity;
    TH1F*    fClusterTime;
    TH1F*    fClusterEnergy;
    TH1F*    fRadius;    // fabs(1/omega)
    TH1F*    fMom;
    TH1F*    fPt;
    TH1F*    fLambda;   //dz/dphi
    TH1F*    fT0;
    TH1F*    fT0Zoom;        // 500-900 ns
    TH1F*    fDeltaT0;       // T0_i - T0_j (i<j) over helix pairs in one event
    TH1F*    fClusterTimeZoom; // calorimeter cluster time, 500-900 ns
    TH1F*    fT0Err;
    TH1F*    fD0;
    TH1F*    fBestAlg;
    TH1F*    fAlgMask;
    TH1F*    fChi2XY;
    TH1F*    fChi2ZPhi;
  };

  struct GenpHist_t {
    TH1F*    fPdgCode[2];		// same distribution in different scale
    TH1F*    fGenID;			// 
    TH1F*    fZ0;			// 
    TH1F*    fT0;			// 
    TH1F*    fR0;			// 
    TH1F*    fP;			// 
    TH1F*    fCosTh;			// 
  };
					// histograms for the simulated CE
  struct SimpHist_t {
    TH1F*    fPdgCode;
    TH1F*    fMomTargetEnd;
    TH1F*    fMomTrackerFront;
    TH1F*    fNStrawHits;
    TH1F*    fNSimp;
    TH1F*    fEffNSimp[10];
    TH1F*    fTStart;        // StartPos().T(), 500-900 ns
    TH1F*    fTEnd;          // EndPos().T(), 500-900 ns
    TH1F*    fMomPt;         // pT
    TH1F*    fMom;           // total momentum
    TH1F*    fDtStart;       // signed dT of StartPos().T() between the 2 pbar-origin particles
    TH1F*    fDtDirect;      // signed dT of the direct pbar daughters' start times for the 2 particles
  };

 struct TrackHist_t {
    TH1F*    fP[3];     // total momentum, 3 hists with different binning
    TH1F*    fP0;
    TH1F*    fP2;
    TH1F*    fPt;
    TH1D*    fPDio;                     // momentum dist weighted with the DIO weight
    TH1F*    fFitMomErr;
    TH1F*    fPFront;
    TH1F*    fDpFront;
    TH1F*    fDpFront0;
    TH1F*    fDpFront2;
    TH2F*    fDpFVsZ1;
    TH1F*    fPStOut;
    TH1F*    fDpFSt;      // P(TT_Hollow) - P(ST_Out)
    TH1F*    fCosTh;
    TH1F*    fChi2;
    TH1F*    fNDof;
    TH1F*    fChi2Dof;
    TH1F*    fNActive;
    TH1F*    fT0;
    TH1F*    fT0Err;
    TH1F*    fQ;
    TH1F*    fFitCons[2];   // fit consistency (0 to 1)
    TH1F*    fD0;
    TH1F*    fZ0;
    TH1F*    fTanDip;
    TH1F*    fResid;
    TH1F*    fAlgMask;

    TH1F*    fXc;
    TH1F*    fYc;
    TH1F*    fPhic;
          // matching histograms
    TH1F*    fNClusters;
    TH1F*    fVaneID;
    TH1F*    fXCal;
    TH1F*    fYCal;
    TH1F*    fZCal;
    TH1F*    fXTrk;
    TH1F*    fYTrk;
    TH1F*    fZTrk;
    TH1F*    fRTrk;
    TH1F*    fDt;     // track-cluster residuals
    TH1F*    fChi2Match;
    TH1F*    fDt_eMinus;
    TH1F*    fDt_ePlus;
    TH1F*    fDt_muMinus;
    TH1F*    fDt_muPlus;
    TH1F*    fDx;
    TH1F*    fDy;
    TH1F*    fDz;
    TH1F*    fDu;
    TH1F*    fDv;
    TH2F*    fDvVsDu;
    TH1F*    fPath;
    TH2F*    fDuVsPath;
    TH2F*    fDucVsPath;
    TH2F*    fDvVsPath;
    TH2F*    fDvcVsPath;
    TH2F*    fDtVsPath;
    TH2F*    fDuVsTDip;
    TH2F*    fDvVsTDip;
    TH1F*    fZ1;
    TH1F*    fECl;
    TH1F*    fEClEKin;
    TH1F*    fEp;
    TH2F*    fEpVsPath;
    TH1F*    fEp_eMinus;
    TH1F*    fEp_ePlus;
    TH1F*    fEp_muMinus;
    TH1F*    fEp_muPlus;
    TH2F*    fNHVsStation;
    TH2F*    fNHVsNSt;

    TH1F*    fRSlope;
    TH1F*    fXSlope;
          // likelihoods
    TH2F*    fEpVsDt;
    TH1F*    fEleLogLHCal;
    TH1F*    fMuoLogLHCal;
    TH1F*    fLogLHRCal;
    TH1F*    fLogLHRDeDx;
    TH1F*    fLogLHRXs;
    TH1F*    fLogLHRTrk;
    TH1F*    fLogLHR;
       // MC truth
    TH1F*    fPdgCode;                  // PDG code of the particle produced most hits
    TH1F*    fFrGH;     // fraction of hits produced by the particle

    TH2F*    fNEPlVsNHPl;
    TH2F*    fNDPlVsNHPl;
    TH2F*    fChi2dVsNDPl;
    TH2F*    fDpFVsNDPl;

    TH1F*    fFrE1;
    TH1F*    fFrE2;
  };


  //-----------------------------------------------------------------------------
  // helix MC truth matching for events with fEffNSimp[1] == 2 (needs TStnHelix v8)
  // owner of a helix = dominant MC particle McSimID(0), purity = McPurity()
  //-----------------------------------------------------------------------------
  struct PbarHelixHist_t {
    TH1F*    fPurityPbar;    // purity of helices owned by one of the 2 pbar-origin particles
    TH1F*    fPurityOther;   // purity of the other helices
    TH1F*    fSharePbar;     // share of the owner in helices owned by a pbar-origin particle
    TH1F*    fTimeDiff[3];   // timediff T0_i - T0_j (signed, i < j): [0] S-pure, [1] S-mixed, [2] B
    TH1F*    fPairCat;       // pair category: 0 S-pure, 1 S-mixed, 2 B, 3 D
                             // stage-0 helix variables, all helices (no cut):
                             // [0] owned by a pbar-origin particle, [1] other
    TH1F*    fNComboHits[2];
    TH1F*    fP[2];
    TH1F*    fD0[2];
    TH1F*    fTanDip[2];
    TH1F*    fPropDir[2];
    TH1F*    fChi2XY[2];
    TH1F*    fChi2ZPhi[2];
    TH1F*    fT0Err[2];
    TH1F*    fSelTruth;      // owners among the kept helices: 0 both pbar particles, 1 one, 2 none, 3 no kept helix
    TH1F*    fSelT0Span;     // max - min T0 of the selected helices
  };

  struct TrackEffHist_t {
    TH1F*    fPtMc;			// denominator
    TH1F*    fPtReco;			// numerator
  };

  enum { kNEventHistSets       =  10 };
  enum { kNHelixHistSets       =  10 };
  enum { kNTimeClusterHistSets =  10};
  enum { kNGenpHistSets        =  10 };
  enum { kNSimpHistSets        =  10 };
  enum { kNTrackHistSets       =  10 };

  struct Hist_t {
    EventHist_t*       fEvent       [kNEventHistSets];
    TimeClusterHist_t* fTimeCluster [kNTimeClusterHistSets];
    HelixHist_t*       fHelix       [kNHelixHistSets];
    GenpHist_t*        fGenp        [kNGenpHistSets];
    SimpHist_t*        fSimp        [kNSimpHistSets];
    TrackHist_t*     fTrack         [kNTrackHistSets];
    PbarHelixHist_t*   fPbarHelix;
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
					// 0: TPR, 1: CPR

  TStnTimeClusterBlock*  fTimeClusterBlock[2];
  TStnHelixBlock*        fHelixBlock;                 // 
  TStrawHitBlock*        fStrawHitBlock;
  TGenpBlock*            fGenpBlock;
  TSimpBlock*            fSimpBlock;
  TStepPointMCBlock*     fVDetBlock;          // virtual detector step points (compressDigiMCs:virtualdetector)
  TStnTrackBlock*        fTrackBlock;

  TString                fHelixBlockName;

					// additional track parameters (assume ntracks < 20)
  TrackPar_t             fTrackPar[20];
					// histograms filled
  Hist_t                 fHist;
					// cut values
  double            fPtMin;

  TGenParticle*     fParticle;		// electron or muon
  int               fPdgCode;		// determines which one
  int               fProcessCode;      

  TSimParticle*     fSimp;
  double            fEleE;		// electron energy

  TStnTrack*        fTrack;
  int               fNTracks[10];
  
  int               fNClusters;
  int               fNTimeClusters[2];

  int               fNHelices[10];	// 0: TPR, 1: CPR, 2: merged
  int               fNHelPos ;
  int               fNHelNeg ;

  int               fNStrawHits;
  int               fNGenp;		// N(generated particles)

  int               fNSimp; // N simparticles
  int               fEffNSimp[10];// Effective N simparticles
  double            fMomPt;   // particle pT
  double            fMomTot;  // particle total momentum

  int               fNHyp;
  int               fBestHyp[10];
  int               fFillDioHist;

  double            fMinT0;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  THelixAnaModule(const char* name="Stntuple_HelixAna", const char* title="Stntuple_HelixAna");
  ~THelixAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*  GetHist        () { return &fHist;        }
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  void     SetPdgCode    (int Code ) { fPdgCode     = Code ; }
  void     SetProcessCode(int Code ) { fProcessCode = Code ; }
  void     SetHelixBlockName(const char* Name) { fHelixBlockName = Name; }
//-----------------------------------------------------------------------------
// overloaded methods of TStnModule
//-----------------------------------------------------------------------------
  int     BeginJob();
  int     BeginRun();
  int     Event   (int ientry);
  int     EndJob  ();
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  void    BookGenpHistograms        (GenpHist_t*    Hist, const char* Folder);
  void    BookSimpHistograms        (SimpHist_t*    Hist, const char* Folder);
  void    BookEventHistograms       (EventHist_t*   Hist, const char* Folder);
  void    BookTrackHistograms       (TrackHist_t*   Hist, const char* Folder);

  void    BookTimeClusterHistograms (TimeClusterHist_t* Hist, const char* Folder);
  void    BookHelixHistograms       (HelixHist_t*       Hist, const char* Folder);
  void    BookPbarHelixHistograms   (PbarHelixHist_t*   Hist, const char* Folder);

  void    FillEventHistograms       (EventHist_t* Hist);
  void    FillGenpHistograms        (GenpHist_t*  Hist, TGenParticle* Genp   );
  void    FillSimpHistograms        (SimpHist_t*  Hist, TSimParticle* Simp   );
  void    FillTrackHistograms       (TrackHist_t*   Hist, TStnTrack*  Trk    );

  void    FillHelixHistograms       (HelixHist_t*       Hist, TStnHelix*       Helix      );
  void    FillTimeClusterHistograms (TimeClusterHist_t* Hist, TStnTimeCluster* TimeCluster);

  void    BookHistograms();
  void    FillHistograms();
                                        // stage-0 helix quality selection
  bool    IsGoodHelix(TStnHelix* Helix);


  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(THelixAnaModule,0)
};
}
#endif
