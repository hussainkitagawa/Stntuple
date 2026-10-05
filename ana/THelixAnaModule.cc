//////////////////////////////////////////////////////////////////////////////
// use of tmp:
//
// use of debug bits: bits 0-2 are reserved
//  0  : all events
//  1  : passed events
//  2  : rejected events
//  3  : N(CalHelixFinder helices hel>0) > 0
//  4  : N(CalHelixFinder helices hel<0) > 0
//  5  : N(helices) > 1
//
///////////////////////////////////////////////////////////////////////////////
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib> // std::strtod
#include "TGraph.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/geom/TDisk.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
//Track files
#include "Offline/RecoDataProducts/inc/RobustHelix.hh"
//#include "Offline/RecoDataProducts/inc/KalSegment.hh"
//#include "Offline/RecoDataProducts/inc/KalSeed.hh"
//#include "Offline/GeometryService/inc/GeometryService.hh"
//#include "Offline/GeometryService/inc/GeomHandle.hh"
//#include "Offline/BTrkData/inc/TrkStrawHit.hh"

#include <cmath>
#include <algorithm>

#include "art/Framework/Principal/Handle.h"

#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"

#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"


#include "TDatabasePDG.h"
//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------
#include "Stntuple/ana/THelixAnaModule.hh"

ClassImp(stntuple::THelixAnaModule)

    namespace stntuple {
  //-----------------------------------------------------------------------------
  THelixAnaModule::THelixAnaModule(const char* name, const char* title) : TStnModule(name, title) {
    fPtMin = 1.;

    fMinT0 = 0; // do not cut on time by default
                //-----------------------------------------------------------------------------
                // MC truth: define which MC particle to consider as signal
                //-----------------------------------------------------------------------------
    fPdgCode = 11;
    fProcessCode = 2; // conversionGun, 28:StoppedParticleReactionGun
                      //-----------------------------------------------------------------------------
                      // default names of the helix block branches correspond to settings of StntupleMakerDe
                      // configuration in Stntuple/fcl/prolog.fcl
                      //-----------------------------------------------------------------------------
    fHelixBlockName = "HelixBlock";
  }

  //-----------------------------------------------------------------------------
  THelixAnaModule::~THelixAnaModule() {
  }


  //-----------------------------------------------------------------------------
  void THelixAnaModule::BookTimeClusterHistograms(TimeClusterHist_t * Hist, const char* Folder) {

    HBook1F(Hist->fNHits, "nhits", Form("%s: # of straw hits", Folder), 75, 0, 150, Folder);
    HBook1F(Hist->fNComboHits, "ncombohits", Form("%s: # of combo hits", Folder), 150, 0, 150, Folder);
    HBook1F(Hist->fT0, "t0", Form("%s: TimeCluster; t_{0}[ns]", Folder), 400, 0, 2000, Folder);
    HBook1F(Hist->fT0Zoom, "t0_zoom", Form("%s: TimeCluster; t_{0}[ns]", Folder), 640, 400, 1200, Folder);
    HBook1F(Hist->fClusterTime, "clusterTime", Form("%s: cluster time; t_{cluster}[ns]", Folder), 400, 500, 900, Folder);
    HBook1F(Hist->fClusterEnergy, "clusterE", Form("%s: cluster energy; E [MeV]      ", Folder), 400, 0, 200, Folder);
  }

  //-----------------------------------------------------------------------------
  void THelixAnaModule::BookHelixHistograms(HelixHist_t * Hist, const char* Folder) {

    HBook1F(Hist->fNHits, "nhits", Form("%s: # of straw hits", Folder), 150, 0, 150, Folder);
    HBook1F(Hist->fHelicity, "hel", Form("%s: Helicity", Folder), 10, -5, 5, Folder);
    HBook1F(Hist->fClusterTime, "clusterTime", Form("%s: cluster time; t_{cluster}[ns]", Folder), 400, 0, 2000, Folder);
    HBook1F(Hist->fClusterEnergy, "clusterE", Form("%s: cluster energy; E [MeV]      ", Folder), 200, 0, 200, Folder);
    HBook1F(Hist->fRadius, "radius", Form("%s: helix radius; r [mm]", Folder), 250, 0, 500, Folder);
    HBook1F(Hist->fMom, "p", Form("%s: momentum; p [MeV/c]", Folder), 250, 0, 500, Folder);
    HBook1F(Hist->fPt, "pT", Form("%s: pT; pT [MeV/c]", Folder), 250, 0, 500, Folder);
    HBook1F(Hist->fLambda, "lambda", Form("%s: lambda; #lambda", Folder), 300, -1500, 1500, Folder);
    HBook1F(Hist->fT0, "t0", Form("%s: t0; t0[ns]", Folder), 100, 0, 2000, Folder);
    HBook1F(Hist->fT0Zoom, "t0_zoom", Form("%s: t0; t0[ns]", Folder), 80, 400, 1200, Folder);
    HBook1F(Hist->fClusterTimeZoom, "clusterTime_zoom", Form("%s: cluster time; t_{cluster}[ns]", Folder), 400, 500, 900, Folder);
    HBook1F(Hist->fDeltaT0, "deltat0", Form("%s: #Deltat0 between helix pairs; t0_{i}-t0_{j} [ns]", Folder), 1000, -1000, 1000, Folder);
    HBook1F(Hist->fT0Err, "t0err", Form("%s: t0err; t0err [ns]", Folder), 100, 0, 10, Folder);
    HBook1F(Hist->fD0, "d0", Form("%s: D0; d0 [mm]", Folder), 1600, -400, 400, Folder);
    HBook1F(Hist->fAlgMask, "algmask", Form("%s: Algorithm Mask", Folder), 10, 0, 10, Folder);
    HBook1F(Hist->fBestAlg, "bestalg", Form("%s: Best Algorithm", Folder), 10, 0, 10, Folder);
    HBook1F(Hist->fChi2XY, "chi2xy", Form("%s: Chi2(XY)/DOF", Folder), 100, 0, 10, Folder);
    HBook1F(Hist->fChi2ZPhi, "chi2zphi", Form("%s: Chi2(ZPhi)/DOF", Folder), 100, 0, 10, Folder);
    //HBook1F(Hist->fChi2ZPhi      ,"chi2zphi"   ,Form("%s: Chi2(ZPhi)/DOF"               ,Folder),  500,     0,   50,Folder);
  }

  //-----------------------------------------------------------------------------
  void THelixAnaModule::BookGenpHistograms(GenpHist_t * Hist, const char* Folder) {
    //   char name [200];
    //   char title[200];

    HBook1F(Hist->fP, "p", Form("%s: Momentum", Folder), 1000, 0, 200, Folder);
    HBook1F(Hist->fPdgCode[0], "pdg_code_0", Form("%s: PDG Code[0]", Folder), 200, -100, 100, Folder);
    HBook1F(Hist->fPdgCode[1], "pdg_code_1", Form("%s: PDG Code[1]", Folder), 500, -2500, 2500, Folder);
    HBook1F(Hist->fGenID, "gen_id", Form("%s: Generator ID", Folder), 100, 0, 100, Folder);
    HBook1F(Hist->fZ0, "z0", Form("%s: Z0", Folder), 500, 5400, 6400, Folder);
    HBook1F(Hist->fT0, "t0", Form("%s: T0", Folder), 200, 0, 2000, Folder);
    HBook1F(Hist->fR0, "r", Form("%s: R0", Folder), 100, 0, 100, Folder);
    HBook1F(Hist->fCosTh, "cos_th", Form("%s: Cos(Theta)", Folder), 200, -1., 1., Folder);
  }


  //-----------------------------------------------------------------------------
  void THelixAnaModule::BookEventHistograms(EventHist_t * Hist, const char* Folder) {
    //  char name [200];
    //  char title[200];

    HBook1F(Hist->fEventNumber, "evtnum", Form("%s: Event Number", Folder), 1000, 0, 1e4, Folder);
    HBook1F(Hist->fRunNumber, "runnum", Form("%s: Run   Number", Folder), 1000, 0, 1e6, Folder);
    HBook1F(Hist->fSectionNumber, "subrunnum", Form("%s: SubRun   Number", Folder), 1000, 0, 1e6, Folder);
    HBook1F(Hist->fEleCosTh, "ce_costh", Form("%s: Conversion Electron Cos(Theta)", Folder), 100, -1, 1, Folder);
    HBook1F(Hist->fEleMom, "ce_mom", Form("%s: Conversion Electron Momentum", Folder), 1000, 0, 200, Folder);
    HBook1F(Hist->fRv, "rv", Form("%s: R(Vertex)", Folder), 100, 0, 1000, Folder);
    HBook1F(Hist->fZv, "zv", Form("%s: Z(Vertex)", Folder), 300, 0, 15000, Folder);
    HBook1F(Hist->fNGenp, "ngenp", Form("%s: N(Gen Particles)", Folder), 500, 0, 500, Folder);

    HBook1F(Hist->fNHelices[0], "nhel_0", Form("%s: N(helices)", Folder), 10, 0, 10, Folder);
    HBook1F(Hist->fNHelices[1], "nhel_1", Form("%s: N(helices)", Folder), 10, 0, 10, Folder);
    // HBook1F(Hist->fNHelices[1],"nhel_1"  ,Form("%s: N(helices)[2]"                   ,Folder), 10,0, 10,Folder);
    // HBook1F(Hist->fNHelices[2],"nhel_2"  ,Form("%s: N(helices)[3]"                   ,Folder), 10,0, 10,Folder);

    HBook1F(Hist->fNTracks[0], "ntrk_0", Form("%s: N(tracks)", Folder), 5, 0, 5, Folder);
    HBook1F(Hist->fNTracks[1], "ntrk_1", Form("%s: N(tracks)", Folder), 5, 0, 5, Folder);

    HBook1F(Hist->fNHelPos, "nhp_0", Form("%s: N(helices) hel > 0", Folder), 10, 0, 10, Folder);
    // HBook1F(Hist->fNHelPos [1],"nhp_1"   ,Form("%s: N(helices)[2] hel > 0"           ,Folder), 10,0, 10,Folder);
    // HBook1F(Hist->fNHelPos [2],"nhp_2"   ,Form("%s: N(helices)[3] hel > 0"           ,Folder), 10,0, 10,Folder);

    HBook1F(Hist->fNHelNeg, "nhn_0", Form("%s: N(helices) hel < 0", Folder), 10, 0, 10, Folder);
    // HBook1F(Hist->fNHelNeg [1],"nhn_1"   ,Form("%s: N(helices)[2] hel < 0"           ,Folder), 10,0, 10,Folder);
    // HBook1F(Hist->fNHelNeg [2],"nhn_2"   ,Form("%s: N(helices)[3] hel < 0"           ,Folder), 10,0, 10,Folder);

    HBook1F(Hist->fNTimeClusters[0], "ntcl_0", Form("%s: N(time clusters)[0]", Folder), 100, 0, 100, Folder);
    HBook1F(Hist->fNTimeClusters[1], "ntcl_1", Form("%s: N(time clusters)[1]", Folder), 100, 0, 100, Folder);
  }

  //-----------------------------------------------------------------------------
  void THelixAnaModule::BookSimpHistograms(SimpHist_t * Hist, const char* Folder) {
    //  char name [200];
    //  char title[200];

    HBook1F(Hist->fPdgCode, "pdg", Form("%s: PDG code", Folder), 200, -100, 100, Folder);
    //HBook1F(Hist->fMom           ,"p"          ,Form("%s: momentum; p [MeV/c]"          ,Folder),  250,    0,  500,Folder);
    HBook1F(Hist->fMomTargetEnd, "ptarg", Form("%s: Mom after Stopping Target", Folder), 250, 0, 500, Folder);
    HBook1F(Hist->fMomTrackerFront, "ptrack", Form("%s: Mom before Tracker", Folder), 250, 0, 500, Folder);
    HBook1F(Hist->fNStrawHits, "nsth", Form("%s: n straw hits", Folder), 200, 0, 200, Folder);
    HBook1F(Hist->fNSimp, "nsimp", Form("%s: N(sim particles)", Folder), 500, 0, 500, Folder);
    HBook1F(Hist->fEffNSimp[0], "effnsimp0", Form("%s: N(effective sim particles)", Folder), 10, 0, 10, Folder);
    HBook1F(Hist->fEffNSimp[1], "effnsimp1", Form("%s: N(sim particles)", Folder), 10, 0, 10, Folder);
    HBook1F(Hist->fTStart, "tstart", Form("%s: sim particle start time; t_{start} [ns]", Folder), 200, 0, 2000, Folder);
    HBook1F(Hist->fTEnd, "tend", Form("%s: sim particle end time; t_{end} [ns]", Folder), 200, 0, 2000, Folder);
    HBook1F(Hist->fMomPt, "momPt", Form("%s: pT; pT [MeV/c]", Folder), 250, 0, 500, Folder);
    HBook1F(Hist->fMom, "mom", Form("%s: momentum; p [MeV/c]", Folder), 250, 0, 500, Folder);
    HBook1F(Hist->fDtStart, "dtstart", Form("%s: #Deltat_{start} of 2 pbar-origin particles; t_{j} - t_{i} [ns]", Folder), 2000, -100, 100, Folder);
    HBook1F(Hist->fDtDirect, "dtdirect", Form("%s: #Deltat of direct pbar daughters of 2 pbar-origin particles; t_{j} - t_{i} [ns]", Folder), 400, -0.002, 0.002, Folder);
  }

  //-----------------------------------------------------------------------------
  void THelixAnaModule::BookTrackHistograms(TrackHist_t * Hist, const char* Folder) {

    HBook1F(Hist->fP[0], "p", Form("%s: Track P(Z1)", Folder), 250, 0, 500., Folder);
    HBook1F(Hist->fP[1], "p_1", Form("%s: Track P(total)[1]", Folder), 100, 104.5, 105.5, Folder);
    HBook1F(Hist->fP[2], "p_2", Form("%s: Track P(total)[1]", Folder), 1000, 0, 500., Folder);
    HBook1F(Hist->fP0, "p0", Form("%s: Track P(Z0)", Folder), 250, 0, 500., Folder);
    HBook1F(Hist->fP2, "p2", Form("%s: Track P(z=-1540)", Folder), 250, 0, 500., Folder);
    HBook1D(Hist->fPDio, "pdio", Form("%s: Track P(DIO WT)", Folder), 400, 90, 110., Folder);
    Hist->fPDio->Sumw2(kTRUE);

    HBook1F(Hist->fFitMomErr, "momerr", Form("%s: Track FitMomError", Folder), 200, 0, 1., Folder);
    HBook1F(Hist->fPFront, "pf", Form("%s: Track P(front)   ", Folder), 400, 90, 110., Folder);
    HBook1F(Hist->fDpFront, "dpf", Form("%s: Track P-P(front) ", Folder), 200, -5., 5., Folder);
    HBook1F(Hist->fDpFront0, "dp0f", Form("%s: Track P0-P(front)", Folder), 200, -5., 5., Folder);
    HBook1F(Hist->fDpFront2, "dp2f", Form("%s: Track P2-P(front)", Folder), 200, -5., 5., Folder);
    HBook1F(Hist->fPStOut, "pstout", Form("%s: Track P(ST_Out)  ", Folder), 400, 90, 110., Folder);
    HBook1F(Hist->fDpFSt, "dpfst", Form("%s: Track Pf-Psto", Folder), 200, -5, 5., Folder);
    HBook2F(Hist->fDpFVsZ1, "dpf_vs_z1", Form("%s: Track DPF Vs Z1", Folder), 200, -2000., 0, 200, -5., 5, Folder);

    HBook1F(Hist->fPt, "pt", Form("%s: Track Pt", Folder), 250, 0, 500., Folder);
    HBook1F(Hist->fCosTh, "costh", Form("%s: Track cos(theta)", Folder), 100, -1, 1, Folder);
    HBook1F(Hist->fChi2, "chi2", Form("%s: Track chi2 total", Folder), 200, 0, 200, Folder);
    HBook1F(Hist->fNDof, "ndof", Form("%s: Number of DOF", Folder), 200, 0, 200, Folder);
    HBook1F(Hist->fChi2Dof, "chi2d", Form("%s: track chi2/N(dof)", Folder), 100, 0, 20, Folder);
    HBook1F(Hist->fNActive, "nactv", Form("%s: N(active)", Folder), 200, 0, 200, Folder);
    HBook1F(Hist->fT0, "t0", Form("%s: track T0", Folder), 200, 0, 2000, Folder);
    HBook1F(Hist->fT0Err, "t0err", Form("%s: track T0Err", Folder), 100, 0, 10, Folder);
    HBook1F(Hist->fQ, "q", Form("%s: track Q", Folder), 4, -2, 2, Folder);
    HBook1F(Hist->fFitCons[0], "fcon", Form("%s: track fit cons [0]", Folder), 200, 0, 1, Folder);
    HBook1F(Hist->fFitCons[1], "fcon1", Form("%s: track fit cons [1]", Folder), 1000, 0, 0.1, Folder);
    HBook1F(Hist->fD0, "d0", Form("%s: track D0      ", Folder), 200, -200, 200, Folder);
    HBook1F(Hist->fZ0, "z0", Form("%s: track Z0      ", Folder), 200, -2000, 2000, Folder);
    HBook1F(Hist->fTanDip, "tdip", Form("%s: track tan(dip)", Folder), 160, -4, 4, Folder);
    HBook1F(Hist->fResid, "resid", Form("%s: hit residuals", Folder), 500, -0.5, 0.5, Folder);
    HBook1F(Hist->fAlgMask, "alg", Form("%s: algorithm mask", Folder), 10, 0, 10, Folder);

    HBook1F(Hist->fXc, "xc", Form("%s: track Xc      ", Folder), 100, -1000, 1000, Folder);
    HBook1F(Hist->fYc, "yc", Form("%s: track Yc      ", Folder), 100, -1000, 1000, Folder);
    HBook1F(Hist->fPhic, "phic", Form("%s: track Phic    ", Folder), 130, 0, 6.5, Folder);

    HBook1F(Hist->fDt, "dt", Form("%s: track delta(T)", Folder), 200, -20, 20, Folder);
    HBook1F(Hist->fChi2Match, "chi2tcm", Form("%s: chi2(t-c match)", Folder), 250, 0, 250, Folder);

    HBook1F(Hist->fDx, "dx", Form("%s: track delta(X)", Folder), 200, -500, 500, Folder);
    HBook1F(Hist->fDy, "dy", Form("%s: track delta(Y)", Folder), 200, -500, 500, Folder);
    HBook1F(Hist->fDz, "dz", Form("%s: track delta(Z)", Folder), 200, -250, 250, Folder);
    HBook1F(Hist->fDz, "dz", Form("%s: track delta(Z)", Folder), 200, -250, 250, Folder);
    HBook1F(Hist->fDu, "du", Form("%s: track-cluster DU)", Folder), 250, -250, 250, Folder);
    HBook1F(Hist->fDv, "dv", Form("%s: track-cluster DV)", Folder), 200, -100, 100, Folder);
    HBook2F(Hist->fDvVsDu, "dv_vs_du", Form("%s: Track Dv Vs Du", Folder), 100, -250, 250, 100, -100., 100, Folder);
    HBook1F(Hist->fPath, "path", Form("%s: track sdisk", Folder), 50, 0, 500, Folder);
    HBook2F(Hist->fDuVsPath, "du_vs_path", Form("%s: Track Du Vs Path", Folder), 50, 0, 500, 200, -200., 200., Folder);
    HBook2F(Hist->fDucVsPath, "duc_vs_path", Form("%s: T-C Duc Vs Path", Folder), 50, 0, 500, 200, -200., 200., Folder);
    HBook2F(Hist->fDvVsPath, "dv_vs_path", Form("%s: T-C  Dv Vs Path", Folder), 50, 0, 500, 200, -200., 200., Folder);
    HBook2F(Hist->fDvcVsPath, "dvc_vs_path", Form("%s: T-C Dvc Vs Path", Folder), 50, 0, 500, 200, -200., 200., Folder);
    HBook2F(Hist->fDtVsPath, "dt_vs_path", Form("%s: T-C DT Vs Path", Folder), 50, 0, 500, 100, -5., 5., Folder);
    HBook2F(Hist->fDuVsTDip, "du_vs_tdip", Form("%s: Track Du Vs TDip", Folder), 100, 0.5, 1.5, 200, -200., 200., Folder);
    HBook2F(Hist->fDvVsTDip, "dv_vs_tdip", Form("%s: Track Dv Vs TDip", Folder), 100, 0.5, 1.5, 200, -200., 200., Folder);

    HBook1F(Hist->fZ1, "z1", Form("%s: track Z1      ", Folder), 200, -2000, 2000, Folder);
    HBook1F(Hist->fNClusters, "ncl", Form("%s: track N(clusters)", Folder), 10, 0, 10, Folder);
    HBook1F(Hist->fVaneID, "vid", Form("%s: track vane ID", Folder), 10, -5, 5, Folder);
    HBook1F(Hist->fXCal, "xcal", Form("%s: track XCal", Folder), 200, -1000, 1000, Folder);
    HBook1F(Hist->fYCal, "ycal", Form("%s: track YCal", Folder), 200, -1000, 1000, Folder);
    HBook1F(Hist->fZCal, "zcal", Form("%s: track ZCal", Folder), 200, 1500, 3500, Folder);
    HBook1F(Hist->fXTrk, "xtrk", Form("%s: track XTrk", Folder), 200, -1000, 1000, Folder);
    HBook1F(Hist->fYTrk, "ytrk", Form("%s: track YTrk", Folder), 200, -1000, 1000, Folder);
    HBook1F(Hist->fRTrk, "rtrk", Form("%s: track RTrk", Folder), 200, -1000, 1000, Folder);
    HBook1F(Hist->fZTrk, "ztrk", Form("%s: track ZTrk", Folder), 200, -1000, 1000, Folder);
    HBook1F(Hist->fECl, "ecl", Form("%s: cluster E", Folder), 300, 0, 150, Folder);
    HBook1F(Hist->fEClEKin, "ecl_ekin", Form("%s: cluster E/Ekin(mu)", Folder), 200, 0, 2, Folder);
    HBook1F(Hist->fEp, "ep", Form("%s: track E/P", Folder), 300, 0, 1.5, Folder);
    HBook2F(Hist->fEpVsPath, "ep_vs_path", Form("%s: E/P Vs Path", Folder), 50, 0, 500, 150, 0., 1.5, Folder);
    HBook2F(Hist->fNHVsStation, "nh_vs_st", Form("%s: N(hits) Vs Station", Folder), 40, 0, 40, 10, -0.5, 9.5, Folder);
    HBook2F(Hist->fNHVsNSt, "nh_vs_nst", Form("%s: N(hits) Vs NSt", Folder), 10, -0.5, 9.5, 40, -0.5, 39.5, Folder);

    HBook1F(Hist->fRSlope, "rslope", Form("%s: Res Slope", Folder), 200, -20, 20, Folder);
    HBook1F(Hist->fXSlope, "xslope", Form("%s: Res/Sig Slope", Folder), 200, -20, 20, Folder);

    HBook2F(Hist->fEpVsDt, "ep_vs_dt", Form("%s: E/P vs Dt", Folder), 200, -10, 10, 150, 0., 1.5, Folder);
    HBook1F(Hist->fEleLogLHCal, "ele_llh_c", Form("%s: ELE Log(LH) Cal", Folder), 200, -100, 0, Folder);
    HBook1F(Hist->fMuoLogLHCal, "muo_llh_c", Form("%s: MUO Log(LH) Cal", Folder), 200, -100, 0, Folder);
    HBook1F(Hist->fLogLHRCal, "llhr_cal", Form("%s: LogLH(e/m) Cal", Folder), 200, -100, 100, Folder);
    HBook1F(Hist->fLogLHRDeDx, "llhr_dedx", Form("%s: LogLH(e/m) De/Dx", Folder), 200, -20, 20, Folder);
    HBook1F(Hist->fLogLHRXs, "llhr_xs", Form("%s: LogLH(e/m) XSlope", Folder), 200, -20, 20, Folder);
    HBook1F(Hist->fLogLHRTrk, "llhr_trk", Form("%s: LogLH(e/m) Trk", Folder), 200, -20, 20, Folder);
    HBook1F(Hist->fLogLHR, "llhr", Form("%s: LogLH(e/m)", Folder), 200, -100, 100, Folder);

    HBook1F(Hist->fPdgCode, "pdg", Form("%s: track PDG code", Folder), 500, -250, 250, Folder);
    HBook1F(Hist->fFrGH, "fgh", Form("%s: Fraction Goog Hits", Folder), 100, 0, 1, Folder);
    HBook2F(Hist->fNEPlVsNHPl, "nep_vs_nhp", Form("%s: Track NEXP vs NHit", Folder), 100, 0, 100, 100, 0., 100, Folder);
    HBook2F(Hist->fNDPlVsNHPl, "ndp_vs_nhp", Form("%s: Track NDIF vs NHit", Folder), 100, 0, 100, 100, 0., 100, Folder);
    HBook2F(Hist->fChi2dVsNDPl, "chi2d_vs_ndp", Form("%s: Track Chi2/Dof vs NDP", Folder), 30, 0, 30, 100, 0., 10, Folder);
    HBook2F(Hist->fDpFVsNDPl, "dpf_vs_ndp", Form("%s: Track DpF vs NDP", Folder), 30, 0, 30, 100, -5, 5, Folder);

    HBook1F(Hist->fFrE1, "fre1", Form("%s: E1/Etot", Folder), 200, 0, 1, Folder);
    HBook1F(Hist->fFrE2, "fre2", Form("%s: (E1+E2)/Etot", Folder), 200, 0, 1, Folder);
  }


  //_____________________________________________________________________________
  void THelixAnaModule::BookHistograms() {

    //  char name [200];
    //  char title[200];

    TFolder* fol;
    TFolder* hist_folder;
    char folder_name[200];

    DeleteHistograms();
    hist_folder = (TFolder*)GetFolder()->FindObject("Hist");
    //--------------------------------------------------------------------------------
    // book timecluster histograms
    //--------------------------------------------------------------------------------
    int book_timecluster_histset[kNTimeClusterHistSets];
    for (int i = 0; i < kNTimeClusterHistSets; ++i) book_timecluster_histset[i] = 0;

    book_timecluster_histset[0] = 1; // all time clusters
    book_timecluster_histset[1] = 1; // time clusters with NHits >= 10

    for (int i = 0; i < kNTimeClusterHistSets; i++) {
      if (book_timecluster_histset[i] != 0) {
        sprintf(folder_name, "timecluster_%i", i);
        fol = (TFolder*)hist_folder->FindObject(folder_name);
        if (!fol) fol = hist_folder->AddFolder(folder_name, folder_name);
        fHist.fTimeCluster[i] = new TimeClusterHist_t;
        BookTimeClusterHistograms(fHist.fTimeCluster[i], Form("Hist/%s", folder_name));
      }
    }
    //--------------------------------------------------------------------------------
    // book helix histograms
    //--------------------------------------------------------------------------------
    int book_helix_histset[kNHelixHistSets];
    for (int i = 0; i < kNHelixHistSets; ++i) book_helix_histset[i] = 0;

    book_helix_histset[0] = 1; // all helices
    book_helix_histset[1] = 1; // helices in events with effective N(SimParticle) = 2
    book_helix_histset[2] = 1; // helices in events with effective N(SimParticle) = 2, chi2xy/chi2zphi <= 5, N(helices) <= 3
                               //-----------------------------------------------------------------------------
                               // the sets below are not filled at the moment - the corresponding calls in
                               // FillHistograms() are commented out. Re-enable both sides together.
                               //-----------------------------------------------------------------------------
                               //book_helix_histset[  3] = 1;   // helices with |lambda| > 1500

    for (int i = 0; i < kNHelixHistSets; i++) {
      if (book_helix_histset[i] != 0) {
        sprintf(folder_name, "helix_%i", i);
        fol = (TFolder*)hist_folder->FindObject(folder_name);
        if (!fol) fol = hist_folder->AddFolder(folder_name, folder_name);
        fHist.fHelix[i] = new HelixHist_t;
        BookHelixHistograms(fHist.fHelix[i], Form("Hist/%s", folder_name));
      }
    }
    //-----------------------------------------------------------------------------
    // book event histograms
    //-----------------------------------------------------------------------------
    int book_event_histset[kNEventHistSets];
    for (int i = 0; i < kNEventHistSets; i++) book_event_histset[i] = 0;

    book_event_histset[0] = 1; // all events
    book_event_histset[1] = 1; // events with N>0 reconstructed helices
    book_event_histset[2] = 1; // events with N=0 reconstructed helices
    book_event_histset[3] = 1; // events with N>0 reconstructed tracks
    book_event_histset[4] = 1; // events with effective N(SimParticle) = 2  (was 5)
    book_event_histset[5] = 1; // all events, filled after the nSimp==2 block (was 6)

    for (int i = 0; i < kNEventHistSets; i++) {
      if (book_event_histset[i] != 0) {
        sprintf(folder_name, "evt_%i", i);
        fol = (TFolder*)hist_folder->FindObject(folder_name);
        if (!fol) fol = hist_folder->AddFolder(folder_name, folder_name);
        fHist.fEvent[i] = new EventHist_t;
        BookEventHistograms(fHist.fEvent[i], Form("Hist/%s", folder_name));
      }
    }
    //-----------------------------------------------------------------------------
    // book simp histograms
    //-----------------------------------------------------------------------------
    int book_simp_histset[kNSimpHistSets];
    for (int i = 0; i < kNSimpHistSets; i++) book_simp_histset[i] = 0;

    book_simp_histset[0] = 1; // fSimp - the "signal" SimParticle
    book_simp_histset[1] = 1; // pbar origin particles                (was 4)
    book_simp_histset[2] = 1; // fSimp, filled in the pbar loop       (was 5)
    book_simp_histset[3] = 1; // fEffNSimp[1] == 2: |dT| of production times (fDtStart only)

    for (int i = 0; i < kNSimpHistSets; i++) {
      if (book_simp_histset[i] != 0) {
        sprintf(folder_name, "sim_%i", i);
        fol = (TFolder*)hist_folder->FindObject(folder_name);
        if (!fol) fol = hist_folder->AddFolder(folder_name, folder_name);
        fHist.fSimp[i] = new SimpHist_t;
        BookSimpHistograms(fHist.fSimp[i], Form("Hist/%s", folder_name));
      }
    }
    //-----------------------------------------------------------------------------
    // book track histograms
    //-----------------------------------------------------------------------------
    int book_track_histset[kNTrackHistSets];
    for (int i = 0; i < kNTrackHistSets; i++) book_track_histset[i] = 0;

    book_track_histset[0] = 1; // all tracks
    book_track_histset[1] = 1; // tracks  NActive >= 15
    book_track_histset[2] = 1; // tracks  NActive >= 30
    book_track_histset[3] = 1; // tracks  TanDip  >  0.8
    book_track_histset[4] = 1; // tracks  TanDip  >  0.8 and NActive >= 30
    book_track_histset[5] = 1; // tracks  nSimp == 1, NActive >= 20, P >= 80 MeV/c
    book_track_histset[6] = 1; // tracks  effective nSimp == 2, NActive >= 15, P >= 60 MeV/c

    for (int i = 0; i < kNTrackHistSets; i++) {
      if (book_track_histset[i] != 0) {
        sprintf(folder_name, "trk_%i", i);
        fol = (TFolder*)hist_folder->FindObject(folder_name);
        if (!fol) fol = hist_folder->AddFolder(folder_name, folder_name);
        fHist.fTrack[i] = new TrackHist_t;
        BookTrackHistograms(fHist.fTrack[i], Form("Hist/%s", folder_name));
      }
    }
    //-----------------------------------------------------------------------------
    // book Genp histograms
    //-----------------------------------------------------------------------------
    int book_genp_histset[kNGenpHistSets];
    for (int i = 0; i < kNGenpHistSets; i++) book_genp_histset[i] = 0;

    book_genp_histset[0] = 1; // all particles

    for (int i = 0; i < kNGenpHistSets; i++) {
      if (book_genp_histset[i] != 0) {
        sprintf(folder_name, "gen_%i", i);
        fol = (TFolder*)hist_folder->FindObject(folder_name);
        if (!fol) fol = hist_folder->AddFolder(folder_name, folder_name);
        fHist.fGenp[i] = new GenpHist_t;
        BookGenpHistograms(fHist.fGenp[i], Form("Hist/%s", folder_name));
      }
    }
  }


  //-----------------------------------------------------------------------------
  // need MC truth branch
  //-----------------------------------------------------------------------------
  void THelixAnaModule::FillEventHistograms(EventHist_t * Hist) {

    int event_number = GetHeaderBlock()->EventNumber();
    int run_number = GetHeaderBlock()->RunNumber();
    int subrun_number = GetHeaderBlock()->SectionNumber();

    double cos_th(-2.), xv(-1.e6), yv(-1.e6), rv(-1.e6), zv(-1.e6), p(-1.);
    TLorentzVector mom;

    if (fParticle) {
      fParticle->Momentum(mom);
      p = mom.P();
      cos_th = mom.Pz() / p;
      xv = fParticle->Vx() + 3904.;
      yv = fParticle->Vy();
      rv = sqrt(xv * xv + yv * yv);
      zv = fParticle->Vz();
    }

    Hist->fEventNumber->Fill(event_number);
    Hist->fRunNumber->Fill(run_number);
    Hist->fSectionNumber->Fill(subrun_number);
    Hist->fEleMom->Fill(p);
    Hist->fEleCosTh->Fill(cos_th);
    Hist->fRv->Fill(rv);
    Hist->fZv->Fill(zv);


    Hist->fNHelices[0]->Fill(fNHelices[0]);
    Hist->fNHelices[1]->Fill(fNHelices[1]);
    Hist->fNTracks[0]->Fill(fNTracks[0]);
    Hist->fNTracks[1]->Fill(fNTracks[1]);
    //for (int i=0; i<3; i++) {
    // Hist->fNHelPos [i]->Fill(fNHelPos [i]);
    // Hist->fNHelNeg [i]->Fill(fNHelNeg [i]);
    //}

    Hist->fNTimeClusters[0]->Fill(fNTimeClusters[0]);
    Hist->fNTimeClusters[1]->Fill(fNTimeClusters[1]);
  }

  //--------------------------------------------------------------------------------
  // function to fill Helix block
  //--------------------------------------------------------------------------------
  void THelixAnaModule::FillTimeClusterHistograms(TimeClusterHist_t * Hist, TStnTimeCluster * TimeCluster) {

    int nhits = TimeCluster->NHits();
    int ncombohits = TimeCluster->NComboHits();

    double time = TimeCluster->T0();
    double clusterE = TimeCluster->ClusterEnergy();

    Hist->fNHits->Fill(nhits);
    Hist->fNComboHits->Fill(ncombohits);
    Hist->fT0->Fill(time);
    Hist->fT0Zoom->Fill(time);
    Hist->fClusterTime->Fill(TimeCluster->ClusterTime());
    Hist->fClusterEnergy->Fill(clusterE);
  }

  //--------------------------------------------------------------------------------
  // function to fill Helix block
  //--------------------------------------------------------------------------------
  void THelixAnaModule::FillHelixHistograms(HelixHist_t * Hist, TStnHelix * Helix) {

    int nhits = Helix->NHits();
    double clusterT = Helix->ClusterTime();
    double clusterE = Helix->ClusterEnergy();

    double radius = Helix->Radius();

    double lambda = Helix->Lambda();
    double tanDip = lambda / radius;
    // assuming B=1T
    double mm2MeV = 3 / 10.;
    double pT = radius * mm2MeV;
    double p = pT / std::cos(std::atan(tanDip));


    Hist->fHelicity->Fill(Helix->Helicity());
    Hist->fNHits->Fill(nhits);
    Hist->fClusterTime->Fill(clusterT);
    Hist->fClusterEnergy->Fill(clusterE);

    Hist->fRadius->Fill(radius);
    Hist->fMom->Fill(p);
    Hist->fPt->Fill(pT);
    Hist->fLambda->Fill(lambda);

    Hist->fBestAlg->Fill(Helix->BestAlg());
    Hist->fAlgMask->Fill(Helix->AlgMask());
    Hist->fD0->Fill(Helix->D0());
    Hist->fT0->Fill(Helix->T0());
    Hist->fT0Zoom->Fill(Helix->T0());
    Hist->fClusterTimeZoom->Fill(clusterT);
    Hist->fT0Err->Fill(Helix->T0Err());
    Hist->fChi2XY->Fill(Helix->Chi2XY());
    Hist->fChi2ZPhi->Fill(Helix->Chi2ZPhi());
  }

  //-----------------------------------------------------------------------------
  void THelixAnaModule::FillGenpHistograms(GenpHist_t * Hist, TGenParticle * Genp) {
    int gen_id;
    float p, cos_th, z0, t0, r0, x0, y0;

    TLorentzVector mom, v;

    Genp->Momentum(mom);
    //  Genp->ProductionVertex(v);

    p = mom.P();
    cos_th = mom.CosTheta();

    x0 = Genp->Vx() + 3904.;
    y0 = Genp->Vy();

    z0 = Genp->Vz();
    t0 = Genp->T();
    r0 = sqrt(x0 * x0 + y0 * y0);
    gen_id = Genp->GetStatusCode();

    Hist->fPdgCode[0]->Fill(Genp->GetPdgCode());
    Hist->fPdgCode[1]->Fill(Genp->GetPdgCode());
    Hist->fGenID->Fill(gen_id);
    Hist->fZ0->Fill(z0);
    Hist->fT0->Fill(t0);
    Hist->fR0->Fill(r0);
    Hist->fP->Fill(p);
    Hist->fCosTh->Fill(cos_th);
  }

  //-----------------------------------------------------------------------------
  void THelixAnaModule::FillSimpHistograms(SimpHist_t * Hist, TSimParticle * Simp) {

    Hist->fPdgCode->Fill(Simp->fPdgCode);
    Hist->fMomTargetEnd->Fill(Simp->fMomTargetEnd);
    Hist->fMomTrackerFront->Fill(Simp->fMomTrackerFront);
    Hist->fNStrawHits->Fill(Simp->NStrawHits());
    Hist->fNSimp->Fill(fNSimp);
    Hist->fEffNSimp[0]->Fill(fEffNSimp[0]);
    Hist->fEffNSimp[1]->Fill(fEffNSimp[1]);
    Hist->fTStart->Fill(Simp->StartPos()->T());
    Hist->fTEnd->Fill(Simp->EndPos()->T());
    Hist->fMomPt->Fill(fMomPt);
    Hist->fMom->Fill(fMomTot);
  }

  //-----------------------------------------------------------------------------
  // for DIO : ultimately, one would need to renormalize the distribution
  //-----------------------------------------------------------------------------
  void THelixAnaModule::FillTrackHistograms(TrackHist_t * Hist, TStnTrack * Track) {

    TLorentzVector mom;
    double r;
    int itrk;
    TrackPar_t* tp;
    // pointer to local track parameters
    itrk = Track->Number();
    tp = fTrackPar + itrk;

    Hist->fP[0]->Fill(Track->fP);
    Hist->fP[1]->Fill(Track->fP);
    Hist->fP[2]->Fill(Track->fP);
    Hist->fP0->Fill(Track->fP0);
    Hist->fP2->Fill(Track->fP2);

    Hist->fPDio->Fill(Track->fP, tp->fDioWt);

    Hist->fFitMomErr->Fill(Track->fFitMomErr);

    Hist->fPt->Fill(Track->fPt);
    Hist->fPFront->Fill(Track->fPFront);
    Hist->fPStOut->Fill(Track->fPStOut);
    // dp: Tracker-only resolution
    Hist->fDpFront->Fill(tp->fDpF);
    Hist->fDpFront0->Fill(tp->fDp0);
    Hist->fDpFront2->Fill(tp->fDp2);
    Hist->fDpFSt->Fill(tp->fDpFSt);
    Hist->fDpFVsZ1->Fill(Track->fZ1, tp->fDpF);

    Hist->fCosTh->Fill(Track->Momentum()->CosTheta());
    Hist->fChi2->Fill(Track->fChi2);
    Hist->fNDof->Fill(Track->NActive() - 5.);
    //Hist->fChi2Dof->Fill(Track->fChi2/(Track->NActive()-5.));
    Hist->fChi2Dof->Fill(Track->Chi2Dof());
    Hist->fNActive->Fill(Track->NActive());
    Hist->fT0->Fill(Track->fT0);
    Hist->fT0Err->Fill(Track->fT0Err);
    Hist->fQ->Fill(Track->Charge());
    Hist->fFitCons[0]->Fill(Track->fFitCons);
    Hist->fFitCons[1]->Fill(Track->fFitCons);

    Hist->fD0->Fill(Track->fD0);
    Hist->fZ0->Fill(Track->fZ0);
    Hist->fTanDip->Fill(Track->fTanDip);
    Hist->fAlgMask->Fill(Track->AlgMask());

    Hist->fXc->Fill(tp->fXc);
    Hist->fYc->Fill(tp->fYc);
    Hist->fPhic->Fill(tp->fPhic);
    //-----------------------------------------------------------------------------
    // track-cluster matching part:
    // - for residuals, determine intersection with the most energetic cluster
    // - for track -only parameters use intersection with lowest trajectory length
    //-----------------------------------------------------------------------------
    TStnTrack::InterData_t* vt = Track->fVMinS; // track-only

    if (vt) {
      Hist->fVaneID->Fill(vt->fID);
      Hist->fXTrk->Fill(vt->fXTrk);
      Hist->fYTrk->Fill(vt->fYTrk);

      r = sqrt(vt->fXTrk * vt->fXTrk + vt->fYTrk * vt->fYTrk);
      Hist->fRTrk->Fill(r);

      Hist->fZTrk->Fill(vt->fZTrk);
    } else {
      //-----------------------------------------------------------------------------
      // fill histograms with numbers easy to recognize as dummy
      //-----------------------------------------------------------------------------
      Hist->fVaneID->Fill(-1.);
      Hist->fXTrk->Fill(999.);
      Hist->fYTrk->Fill(999.);
      Hist->fRTrk->Fill(999.);
      Hist->fZTrk->Fill(-1.);
    }
    //-----------------------------------------------------------------------------
    // there is an inconsistency in the SIMP block filling - in Mu2e offline
    // the particle momentumis is kept in MeV/c, while the PDG mass  -in GeV/c^2..
    // thus the energy is screwed up... kludge around
    // assign muon mass
    //-----------------------------------------------------------------------------
    double ekin(-1.);
    if (fParticle) {
      double p, m;
      //    p    = fSimp->fStartMom.P();
      p = Track->fP;
      m = 105.658; // in MeV
      ekin = sqrt(p * p + m * m) - m;
    }
    Hist->fECl->Fill(tp->fEcl);
    Hist->fEClEKin->Fill(tp->fEcl / ekin);
    Hist->fEp->Fill(tp->fEp);
    Hist->fEpVsPath->Fill(tp->fPath, tp->fEp);

    Hist->fDx->Fill(tp->fDx);
    Hist->fDy->Fill(tp->fDy);
    Hist->fDz->Fill(tp->fDz);

    Hist->fDt->Fill(tp->fDt);
    Hist->fChi2Match->Fill(tp->fChi2Match);

    Hist->fDu->Fill(tp->fDu);
    Hist->fDv->Fill(tp->fDv);
    Hist->fDvVsDu->Fill(tp->fDu, tp->fDv);

    Hist->fPath->Fill(tp->fPath);
    Hist->fDuVsPath->Fill(tp->fPath, tp->fDu);

    //  double cu[4] = { -60.3049, -0.749111,    0.00522242,  -7.52018e-06};
    double cu[4] = {-59.5174, -0.541226, 0.00414309, -5.84989e-06};
    double cv[4] = {6.44161, 0.0722353, -0.000653084, 1.14054e-06};

    double x = tp->fPath;
    double corr_u = cu[0] + cu[1] * x + cu[2] * x * x + cu[3] * x * x * x;
    double corr_v = cv[0] + cv[1] * x + cv[2] * x * x + cv[3] * x * x * x;

    //  double duc = tp->fDu-0.34*(tp->fPath-350.);
    double duc = tp->fDu - corr_u;
    double dvc = tp->fDv - corr_v;

    Hist->fDucVsPath->Fill(tp->fPath, duc);
    Hist->fDvcVsPath->Fill(tp->fPath, dvc);

    Hist->fDvVsPath->Fill(tp->fPath, tp->fDv);
    Hist->fDtVsPath->Fill(tp->fPath, tp->fDt);

    Hist->fDuVsTDip->Fill(Track->fTanDip, tp->fDu);
    Hist->fDvVsTDip->Fill(Track->fTanDip, tp->fDv);

    Hist->fZ1->Fill(Track->fZ1);

    int ncl = Track->NClusters();
    Hist->fNClusters->Fill(ncl);
    Hist->fRSlope->Fill(Track->RSlope());
    Hist->fXSlope->Fill(Track->XSlope());

    double llhr_dedx, llhr_xs, llhr_cal, llhr_trk, llhr;

    Hist->fEleLogLHCal->Fill(Track->EleLogLHCal());
    Hist->fMuoLogLHCal->Fill(Track->MuoLogLHCal());

    llhr_cal = Track->LogLHRCal();
    Hist->fLogLHRCal->Fill(llhr_cal);

    llhr_dedx = Track->LogLHRDeDx();
    llhr_xs = Track->LogLHRXs();
    llhr_trk = Track->LogLHRTrk();
    llhr = llhr_cal + llhr_trk;

    Hist->fEpVsDt->Fill(tp->fDt, tp->fEp);
    Hist->fLogLHRDeDx->Fill(llhr_dedx);
    Hist->fLogLHRXs->Fill(llhr_xs);
    Hist->fLogLHRTrk->Fill(llhr_trk);
    Hist->fLogLHR->Fill(llhr);

    Hist->fPdgCode->Fill(Track->fPdgCode);
    Hist->fFrGH->Fill(Track->fNGoodMcHits / (Track->NActive() + 1.e-5));

    Hist->fNEPlVsNHPl->Fill(tp->fNEPl, tp->fNHPl);
    Hist->fNDPlVsNHPl->Fill(tp->fNDPl, tp->fNHPl);
    Hist->fChi2dVsNDPl->Fill(tp->fNDPl, Track->Chi2Dof());
    Hist->fDpFVsNDPl->Fill(tp->fNDPl, tp->fDpF);
  }


  //-----------------------------------------------------------------------------
  // register data blocks and book histograms
  //-----------------------------------------------------------------------------
  int THelixAnaModule::BeginJob() {
    //-----------------------------------------------------------------------------
    // register data blocks, 'HelixBlock' - OR of TPR and CPR
    //-----------------------------------------------------------------------------
    RegisterDataBlock(fHelixBlockName.Data(), "TStnHelixBlock", &fHelixBlock);

    RegisterDataBlock("TimeClusterBlock", "TStnTimeClusterBlock", &fTimeClusterBlock[0]);
    RegisterDataBlock("PhClusterBlock", "TStnTimeClusterBlock", &fTimeClusterBlock[1]);

    RegisterDataBlock("GenpBlock", "TGenpBlock", &fGenpBlock);
    RegisterDataBlock("SimpBlock", "TSimpBlock", &fSimpBlock);
    RegisterDataBlock("TrackBlock", "TStnTrackBlock", &fTrackBlock);
    //-----------------------------------------------------------------------------
    // book histograms
    //-----------------------------------------------------------------------------
    BookHistograms();

    return 0;
  }


  //_____________________________________________________________________________
  int THelixAnaModule::BeginRun() {
    int rn = GetHeaderBlock()->RunNumber();
    TStntuple::Init(rn);
    return 0;
  }


  //_____________________________________________________________________________
  void THelixAnaModule::FillHistograms() {

    int event_number = GetHeaderBlock()->EventNumber();
    int run_number = GetHeaderBlock()->RunNumber();
    int subrun_number = GetHeaderBlock()->SectionNumber();
    std::cout << "======================================" << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "          FillHistograms              " << std::endl;
    std::cout << "run/subrun/eventNumber = " << run_number << "/" << subrun_number << "/" << event_number << std::endl;
    std::cout << run_number << "_" << subrun_number << "_" << event_number << std::endl;
    std::cout << "======================================" << std::endl;
    std::cout << "======================================" << std::endl;


    //-----------------------------------------------------------------------------
    // global helper functions
    //-----------------------------------------------------------------------------
    int nSimp = 0;
    fEffNSimp[0] = -1;
    fMomPt = 0;
    fMomTot = 0;


    //---------------------------------------------
    // helper: same kinematic definition everywhere
    //---------------------------------------------
    auto GetSimpKinematics = [](TSimParticle* simp,
                                double& px, double& py, double& pz,
                                double& p, double& pt,
                                const TLorentzVector*& pos0,
                                const TLorentzVector*& mom0) -> bool {
      if (!simp) return false;

      pos0 = simp->StartPos();
      mom0 = simp->StartMom();

      if (!pos0 || !mom0) return false;

      TVector3 mom0_vec = mom0->Vect();

      px = mom0_vec.X();
      py = mom0_vec.Y();
      pz = mom0_vec.Z();

      p = mom0_vec.Mag();
      pt = mom0_vec.Pt();

      return true;
    };

    fEffNSimp[0] = nSimp;

    //---------------------------------------------
    // helper:
    // Store useful info for particles whose ancestry traces back
    // to the primary pbar annihilation
    //---------------------------------------------

    struct PbarOriginParticle_t {
      int i;
      TSimParticle* simp;

      int pdg_code;
      int simp_id;
      int parent_id;
      int creationCode;
      int nhits;

      float ptrk;

      double px;
      double py;
      double pz;
      double p;
      double pt;

      const TLorentzVector* pos0;
      const TLorentzVector* mom0;

      // ancestry information
      bool fromPrimaryPbar;
      int ancestryDepth;

      int directDaughterID;
      int directDaughterPDG;
      int directDaughterCreationCode;

      double directDaughter_dR;
      double directDaughter_dT;
      bool directDaughterSameVertex;
    };

    std::vector<PbarOriginParticle_t> pbarOriginParticles;


    //------------------------------------------------------------
    // Analysis cuts
    //------------------------------------------------------------
    //const int    min_nhit_pbar = 1;     // require tracker hits
    const int min_nhit_pbar = 20;   // require tracker hits
    const double min_p_pbar = 80.0; // MeV/c, tune to 60 or 70

    const double vertexTolerance = 1.0e-2; // mm
    const double timeTolerance = 1.0e-3;   // same unit as StartPos.T()


    std::cout << "======================================" << std::endl;
    std::cout << "pbar_origin_particle_check_nested_loop" << std::endl;
    std::cout << "======================================" << std::endl;


    //------------------------------------------------------------
    // Step 1: find primary pbar
    //------------------------------------------------------------
    int primaryPbarID = -1;
    int primaryPbarIndex = -1;

    TVector3 primaryPbarVtx;
    double primaryPbarTime = 0.0;

    bool foundPrimaryPbar = false;

    for (int i = 0; i < fNSimp; ++i) {
      TSimParticle* s = fSimpBlock->Particle(i);
      if (!s) continue;

      int sid = s->GetUniqueID();
      int pid = s->ParentID();
      int pdg = s->PDGCode();

      const TLorentzVector* pos = s->StartPos();
      if (!pos) continue;

      if (pdg == -2212 && pid < 0) {
        primaryPbarID = sid;
        primaryPbarIndex = i;

        primaryPbarVtx.SetXYZ(pos->X(), pos->Y(), pos->Z());
        primaryPbarTime = pos->T();

        foundPrimaryPbar = true;

        std::cout << "[primary pbar found] "
                  << "index=" << i
                  << " simID=" << sid
                  << " parentID=" << pid
                  << " pdg=" << pdg
                  << " creationCode=" << s->fCreationCode
                  << " startPos=("
                  << pos->X() << ", "
                  << pos->Y() << ", "
                  << pos->Z() << ", "
                  << pos->T() << ")"
                  << std::endl;

        break;
      }
    }


    if (!foundPrimaryPbar) {
      std::cout << "[WARNING] No primary pbar found in this event." << std::endl;
    } else {

      //------------------------------------------------------------
      // Step 2:
      // Loop over all MC particles.
      // First require tracker hits and momentum.
      // Then trace parent -> grandparent -> ... back to pbar.
      //------------------------------------------------------------
      for (int i = 0; i < fNSimp; ++i) {

        TSimParticle* s = fSimpBlock->Particle(i);
        if (!s) continue;

        int sid = s->GetUniqueID();
        int pid = s->ParentID();
        int pdg = s->PDGCode();
        int nhits = s->NStrawHits();

        // skip primary pbar itself
        if (sid == primaryPbarID) continue;

        //--------------------------------------------------------
        // Condition 1: this MC particle must make tracker hits
        //--------------------------------------------------------
        if (nhits < min_nhit_pbar) continue;

        //--------------------------------------------------------
        // Get kinematics of this tracker-hit particle
        //--------------------------------------------------------
        const TLorentzVector* pos0 = nullptr;
        const TLorentzVector* mom0 = nullptr;

        double px = 0.0;
        double py = 0.0;
        double pz = 0.0;
        double p = 0.0;
        double pt = 0.0;

        bool goodKin = GetSimpKinematics(s, px, py, pz, p, pt, pos0, mom0);

        if (!goodKin) {
          std::cout << "[WARNING] tracker-hit particle has null StartPos or StartMom, skipping. "
                    << "index=" << i
                    << " simID=" << sid
                    << std::endl;
          continue;
        }

        //--------------------------------------------------------
        // Condition 2: momentum threshold
        //--------------------------------------------------------
        if (p < min_p_pbar) continue;


        //--------------------------------------------------------
        // Step 3:
        // Trace ancestry backward using nested loops.
        //
        // current particle -> parent -> grandparent -> ...
        //--------------------------------------------------------
        bool fromPrimaryPbar = false;
        //bool fromPrimaryPbar = true;

        int ancestryDepth = -1;

        int directDaughterID = -1;
        int directDaughterPDG = 0;
        int directDaughterCreationCode = -999;

        double directDaughter_dR = 1.0e30;
        double directDaughter_dT = 1.0e30;
        bool directDaughterSameVertex = false;
        //bool directDaughterSameVertex = true;

        TSimParticle* cur = s;

        for (int depth = 0; depth < 50000; ++depth) {

          if (!cur) break;

          int curID = cur->GetUniqueID();
          int curParentID = cur->ParentID();

          //------------------------------------------------------
          // If cur's parent is the primary pbar, then cur is the
          // direct daughter from pbar annihilation.
          //------------------------------------------------------
          if (curParentID == primaryPbarID) {

            fromPrimaryPbar = true;
            ancestryDepth = depth + 1;

            directDaughterID = curID;
            directDaughterPDG = cur->PDGCode();
            directDaughterCreationCode = cur->fCreationCode;

            const TLorentzVector* dpos = cur->StartPos();

            if (dpos) {
              TVector3 directDaughterVtx(dpos->X(), dpos->Y(), dpos->Z());

              directDaughter_dR = (directDaughterVtx - primaryPbarVtx).Mag();
              directDaughter_dT = dpos->T() - primaryPbarTime;

              directDaughterSameVertex =
                  (directDaughter_dR < vertexTolerance &&
                   std::abs(directDaughter_dT) < timeTolerance);
            }

            break;
          }

          //------------------------------------------------------
          // Search the whole SimParticle list to find cur's parent.
          //------------------------------------------------------
          TSimParticle* parent = nullptr;

          for (int j = 0; j < fNSimp; ++j) {

            TSimParticle* candidateParent = fSimpBlock->Particle(j);
            if (!candidateParent) continue;

            if ((int)candidateParent->GetUniqueID() == curParentID) {
              parent = candidateParent;
              break;
            }
          }

          //------------------------------------------------------
          // If parent was not found, ancestry chain stops here.
          //------------------------------------------------------
          if (!parent) {
            break;
          }

          //------------------------------------------------------
          // Move one generation upward.
          //------------------------------------------------------
          cur = parent;
        }


        //--------------------------------------------------------
        // Step 4:
        // Require that ancestry reaches pbar annihilation.
        //--------------------------------------------------------
        if (!fromPrimaryPbar) continue;

        //--------------------------------------------------------
        // Step 5:
        // Require that the direct daughter of pbar was produced
        // at the same pbar annihilation vertex.
        //
        // Important:
        // We do NOT require the final tracker-hit particle itself
        // to start at the pbar vertex.
        //--------------------------------------------------------
        if (!directDaughterSameVertex) continue;


        //--------------------------------------------------------
        // Store selected pbar-origin tracker-hit particle
        //--------------------------------------------------------
        PbarOriginParticle_t d;

        d.i = i;
        d.simp = s;

        d.pdg_code = pdg;
        d.simp_id = sid;
        d.parent_id = pid;
        d.creationCode = s->fCreationCode;
        d.nhits = nhits;

        d.ptrk = s->fMomTrackerFront;

        d.px = px;
        d.py = py;
        d.pz = pz;
        d.p = p;
        d.pt = pt;

        d.pos0 = pos0;
        d.mom0 = mom0;

        d.fromPrimaryPbar = fromPrimaryPbar;
        d.ancestryDepth = ancestryDepth;

        d.directDaughterID = directDaughterID;
        d.directDaughterPDG = directDaughterPDG;
        d.directDaughterCreationCode = directDaughterCreationCode;

        d.directDaughter_dR = directDaughter_dR;
        d.directDaughter_dT = directDaughter_dT;
        d.directDaughterSameVertex = directDaughterSameVertex;

        pbarOriginParticles.push_back(d);


        std::cout << "[pbar-origin tracker particle] "
                  << "index=" << d.i
                  << " simID=" << d.simp_id
                  << " parentID=" << d.parent_id
                  << " pdg=" << d.pdg_code
                  << " creationCode=" << d.creationCode
                  << " nhits=" << d.nhits
                  << " p=" << d.p
                  << " pt=" << d.pt
                  << " ancestryDepth=" << d.ancestryDepth
                  << " directDaughterID=" << d.directDaughterID
                  << " directDaughterPDG=" << d.directDaughterPDG
                  << " directDaughterCreationCode=" << d.directDaughterCreationCode
                  << " directDaughter_dR=" << d.directDaughter_dR
                  << " directDaughter_dT=" << d.directDaughter_dT
                  << " directDaughterSameVertex=" << d.directDaughterSameVertex
                  << " startPos=("
                  << d.pos0->X() << ", "
                  << d.pos0->Y() << ", "
                  << d.pos0->Z() << ", "
                  << d.pos0->T() << ")"
                  << std::endl;
      }


      //------------------------------------------------------------
      // Step 6: summary
      //------------------------------------------------------------
      std::cout << "--------------------------------------" << std::endl;
      std::cout << "Number of pbar-origin tracker particles = "
                << pbarOriginParticles.size() << std::endl;

      std::cout << "Cuts: nhits >= " << min_nhit_pbar
                << ", p >= " << min_p_pbar << " MeV/c"
                << std::endl;

      std::cout << "Primary pbar index=" << primaryPbarIndex
                << " simID=" << primaryPbarID
                << std::endl;

      std::cout << "run/subrun/eventNumber = "
                << run_number << "/"
                << subrun_number << "/"
                << event_number
                << std::endl;

      std::cout << "--------------------------------------" << std::endl;
    }

    fEffNSimp[1] = 0;
    for (const auto& d : pbarOriginParticles) {
      TSimParticle* simp = d.simp;
      if (!simp) continue;

      if (d.nhits < 20) continue;
      fMomPt = d.pt;
      fMomTot = d.p;
      fEffNSimp[1]++;
      FillSimpHistograms(fHist.fSimp[1], simp);
    }
    //if(fEffNSimp[1] != 0) FillSimpHistograms(fHist.fSimp[2], fSimp);
    if (fEffNSimp[1] == 2) {
      for (const auto& d : pbarOriginParticles) {

        TSimParticle* s = d.simp;
        if (!s) continue;

        int sid = s->GetUniqueID();
        int nhits = s->NStrawHits();

        // skip primary pbar itself
        if (sid == primaryPbarID) continue;

        //--------------------------------------------------------
        // Condition 1: this MC particle must make tracker hits
        //--------------------------------------------------------
        if (nhits < min_nhit_pbar) continue;

        //--------------------------------------------------------
        // Get kinematics of this tracker-hit particle
        //--------------------------------------------------------
        const TLorentzVector* pos0 = nullptr;
        const TLorentzVector* mom0 = nullptr;

        double px = 0.0;
        double py = 0.0;
        double pz = 0.0;
        double p = 0.0;
        double pt = 0.0;

        bool goodKin = GetSimpKinematics(s, px, py, pz, p, pt, pos0, mom0);

        if (!goodKin) {
          std::cout << "[WARNING] tracker-hit particle has null StartPos or StartMom, skipping. "
                    << "index=" << d.i
                    << " simID=" << sid
                    << std::endl;
          continue;
        }

        //--------------------------------------------------------
        // Condition 2: momentum threshold
        //--------------------------------------------------------
        if (p < min_p_pbar) continue;
        fMomPt = pt;
        fMomTot = p;
        FillSimpHistograms(fHist.fSimp[2], s);
      }

      //--------------------------------------------------------
      // Production-time difference between the 2 pbar-origin
      // particles (nhits >= 20, same pbar annihilation vertex).
      // Collect StartPos().T() and fill signed dT = t_j - t_i for every pair
      // (i < j, in SimpBlock order).
      //
      // Same for the direct pbar daughters each particle descends from:
      // directDaughter_dT = T(direct daughter) - T(pbar), so the pairwise
      // difference equals the difference of the direct daughters' start times.
      //--------------------------------------------------------
      std::vector<double> tStart;
      std::vector<double> tDirect;
      for (const auto& d : pbarOriginParticles) {
        if (!d.simp) continue;
        if (d.nhits < 20) continue;
        if (!d.pos0) continue;
        tStart.push_back(d.pos0->T());
        tDirect.push_back(d.directDaughter_dT);
      }

      for (size_t i = 0; i < tStart.size(); ++i) {
        for (size_t j = i + 1; j < tStart.size(); ++j) {
          double dt = tStart[j] - tStart[i];
          fHist.fSimp[3]->fDtStart->Fill(dt);

          double dt_direct = tDirect[j] - tDirect[i];
          fHist.fSimp[3]->fDtDirect->Fill(dt_direct);
        }
      }
    }


    //-----------------------------------------------------------------------------
    // fill Simp histograms
    //-----------------------------------------------------------------------------
    if (fSimp) FillSimpHistograms(fHist.fSimp[0], fSimp);


    //--------------------------------------------------------------------------------
    // fill TimeCluster histograms
    //--------------------------------------------------------------------------------
    TStnTimeCluster* tCluster = nullptr;

    for (int i = 0; i < fNTimeClusters[0]; ++i) {
      tCluster = fTimeClusterBlock[0]->TimeCluster(i);
      if (!tCluster) continue;

      FillTimeClusterHistograms(fHist.fTimeCluster[0], tCluster);

      int nhits = tCluster->NHits();
      if (nhits >= 10) FillTimeClusterHistograms(fHist.fTimeCluster[1], tCluster);
    }

    //--------------------------------------------------------------------------------
    // fill helix histograms
    //--------------------------------------------------------------------------------
    TStnHelix* helix;
    double mm2MeV{3 / 10.};

    int hist_offset = 0;

    for (int i = 0; i < fNHelices[0]; ++i) {

      helix = fHelixBlock->Helix(i);


      // int         nhits    = helix->NHits();
      double radius = helix->Radius();
      double lambda = helix->Lambda();
      double tanDip = lambda / radius;
      // [[maybe_unused]]: these are kept for the diagnostics below, which are
      // commented out at the moment (pre-existing -Werror=unused-variable)
      [[maybe_unused]] double p = radius * mm2MeV / cos(atan(tanDip));
      [[maybe_unused]] double pT = radius * mm2MeV;
      [[maybe_unused]] int nhits = helix->NHits();

      [[maybe_unused]] double chi2xy = helix->Chi2XY();
      [[maybe_unused]] double chi2zphi = helix->Chi2ZPhi();
      //double p        = pT / std::cos(std::atan(tanDip));
      //double rcent = helix->RCent();
      //double fcent = polyAtan2(yC, xC)
      //double fcent = helix->FCent();
      //double fz0 = helix->FZ0();
      //double helix_xc = rcent*cos(fcent);
      //double helix_yc = rcent*sin(fcent);

      FillHelixHistograms(fHist.fHelix[hist_offset + 0], helix);
      if (fEffNSimp[1] == 2 and chi2xy <= 5.0 and chi2zphi <= 5.0 and fNHelices[0] <= 2) FillHelixHistograms(fHist.fHelix[hist_offset + 1], helix);
      if (fEffNSimp[1] == 2 and chi2xy <= 5.0 and chi2zphi <= 5.0 and fNHelices[0] >= 3) FillHelixHistograms(fHist.fHelix[hist_offset + 2], helix);
      //if (abs(lambda) <   100.) FillHelixHistograms(fHist.fHelix[hist_offset+2], helix);
      //if (abs(lambda) >  1500.) FillHelixHistograms(fHist.fHelix[hist_offset+3], helix);

      /*if (nSimp == 2 and chi2zphi > 50 and fNTracks[0] == 0){
        std::cout<<"fNSimp_2_chi2zphi_5"<<std::endl;
        std::cout<<"chi2xy/chi2zphi = "<<chi2xy<<"/"<<chi2zphi<<std::endl;
        std::cout<<"run/subrun/eventNumber = "<<run_number<<"/"<<subrun_number<<"/"<<event_number<<std::endl;
      }*/
    }

    //--------------------------------------------------------------------------------
    // delta T0 = T0_i - T0_j between all pairs of helices in the event (i<j)
    //--------------------------------------------------------------------------------
    const double deltaT0DumpThre = -400.0; // ns, dump events having a helix pair with deltaT0 below this
    bool deltaT0BelowThre = false;
    double deltaT0Min = 0.0;
    int deltaT0MinI = -1;
    int deltaT0MinJ = -1;

    static std::ofstream g_deltaT0Out;
    static bool g_deltaT0Opened = false;
    if (!g_deltaT0Opened) {
      g_deltaT0Out.open("/exp/mu2e/app/users/kitagawa/20251126/deltaT0.txt", std::ios::out | std::ios::app);
      if (!g_deltaT0Out) std::cerr << "[WARN] Can't open deltaT0.txt for append.\n";
      else g_deltaT0Out << "# run subrun event NHelices idx_i idx_j T0_i[ns] T0_j[ns] deltaT0[ns]\n";
      g_deltaT0Opened = true;
    }

    for (int i = 0; i < fNHelices[0]; ++i) {
      TStnHelix* helix_i = fHelixBlock->Helix(i);
      for (int j = i + 1; j < fNHelices[0]; ++j) {
        TStnHelix* helix_j = fHelixBlock->Helix(j);

        double deltaT0 = helix_i->T0() - helix_j->T0();

        if (fEffNSimp[1] == 2 and helix_i->Chi2XY() <= 5.0 and helix_i->Chi2ZPhi() <= 5.0 and helix_j->Chi2XY() <= 5.0 and helix_j->Chi2ZPhi() <= 5.0) {
          fHist.fHelix[hist_offset + 0]->fDeltaT0->Fill(deltaT0);

          if (deltaT0 < deltaT0DumpThre and (not deltaT0BelowThre or deltaT0 < deltaT0Min)) {
            deltaT0BelowThre = true;
            deltaT0Min = deltaT0;
            deltaT0MinI = i;
            deltaT0MinJ = j;
          }

          if (g_deltaT0Out.is_open()) {
            g_deltaT0Out << run_number << " " << subrun_number << " " << event_number << " " << fNHelices[0] << " " << i << " " << j
                         << std::fixed << std::setprecision(3) << " " << helix_i->T0() << " " << helix_j->T0() << " " << deltaT0 << "\n";
            g_deltaT0Out.flush();
          }
        }
      }
    }

    //--------------------------------------------------------------------------------
    // fill track histograms
    //--------------------------------------------------------------------------------
    int nhits_track_thre = 20;    //threshold of track StrawHits
    double mom_track_thre = 80.0; //threshold of track momentum
    TStnTrack* trk;
    //  TrackPar_t*  tp;
    for (int i = 0; i < fNTracks[0]; ++i) {
      trk = fTrackBlock->Track(i);

      FillTrackHistograms(fHist.fTrack[0], trk);
      if (trk->NActive() >= 20) FillTrackHistograms(fHist.fTrack[1], trk);
      if (fEffNSimp[1] == 2) FillTrackHistograms(fHist.fTrack[2], trk);
      if (fEffNSimp[1] == 2 and trk->NActive() >= 20) FillTrackHistograms(fHist.fTrack[3], trk);
      if (fEffNSimp[1] == 2 and trk->NActive() >= 20 and trk->P() >= mom_track_thre) FillTrackHistograms(fHist.fTrack[4], trk);
      if (fEffNSimp[1] == 1) FillTrackHistograms(fHist.fTrack[5], trk);
      if (nSimp == 1 and trk->NActive() >= nhits_track_thre and trk->P() >= 80.0) FillTrackHistograms(fHist.fTrack[6], trk);
    }

    //--------------------------------------------------------------------------------
    // fill event histograms
    //--------------------------------------------------------------------------------
    FillEventHistograms(fHist.fEvent[0]);

    if (fNHelices[0] > 0) FillEventHistograms(fHist.fEvent[1]);
    else FillEventHistograms(fHist.fEvent[2]);
    if (fNTracks[0] > 0) FillEventHistograms(fHist.fEvent[3]);

    if (fEffNSimp[1] == 2) {
      int nhelix = 0;
      for (int i = 0; i < fNHelices[0]; ++i) {
        auto helix = fHelixBlock->Helix(i);
        double chi2xy = helix->Chi2XY();
        double chi2zphi = helix->Chi2ZPhi();
        //int    nhits    = helix->NHits();
        //double     mm2MeV{3/10.};
        //double      radius   = helix->Radius();
        //double      lambda   = helix->Lambda();
        //double      tanDip   = lambda/radius;
        //double      pT       = radius*mm2MeV;
        //if(nhits < 10) continue;
        if (chi2xy > 5.0) continue;
        if (chi2zphi > 5.0) continue;
        nhelix++;
      }
      fNHelices[1] = nhelix;

      int ntrack = 0;
      //TStnTrack*   trk;
      for (int i = 0; i < fNTracks[0]; ++i) {
        auto trk = fTrackBlock->Track(i);
        if (trk->NActive() < nhits_track_thre) continue;
        if (trk->P() < mom_track_thre) continue;
        ntrack++;
      }
      fNTracks[1] = ntrack;
      FillEventHistograms(fHist.fEvent[4]);
    }
    FillEventHistograms(fHist.fEvent[5]);
    //  first_entry = 0;


    //-----------------------------------------------------------------------------
    // fill GENP histograms
    //-----------------------------------------------------------------------------
    //GEN_0: all particles
    TGenParticle* genp;
    for (int i = 0; i < fNGenp; i++) {
      genp = fGenpBlock->Particle(i);
      FillGenpHistograms(fHist.fGenp[0], genp);
    }

    //--------------------------------------------------------------------------------
    // Dump txt file
    //--------------------------------------------------------------------------------
#include <fstream>
#include <iomanip>
#include <sstream>

    // open once, reuse
    static std::ofstream g_selOut;
    static bool g_filesOpened = false;

    auto ensure_out_open = []() {
      if (!g_filesOpened) {
        // Open the original events file
        g_selOut.open("/exp/mu2e/app/users/kitagawa/20251126/selected_events.txt",
                      std::ios::out | std::ios::app);
        if (!g_selOut) {
          std::cerr << "[WARN] Can't open selected_events.txt for append.\n";
        }

        g_filesOpened = true;
      }
    };

    // =========================================================
    // Skip events with >= 2 tracks in fNTracks[1]
    // =========================================================
    if (fEffNSimp[1] == 2 and fNHelices[1] == 2 and fNTracks[1] <= 1) {

      // =========================================================
      // Prepare event header (run/sub/event + counts)
      // =========================================================
      ensure_out_open();

      std::ostringstream tag;
      tag << run_number << ":"
          << std::setw(2) << std::setfill('0') << subrun_number << ":"
          << std::setw(3) << std::setfill('0') << event_number;

      std::ostringstream header;
      header << "\n============================================================\n";
      header << "Event: " << tag.str() << "\n";
      header << "Tracks (fNTracks[1]):  " << fNTracks[1] << "\n";
      header << "Helices (fNHelices[1]): " << fNHelices[1] << "\n";
      header << "============================================================\n\n";

      if (g_selOut.is_open()) {
        g_selOut << header.str();
        g_selOut.flush();
      }
      std::cout << header.str();


      // =========================================================
      // Dump TRACK INFO only
      // =========================================================

      std::ostringstream trackTitle;
      trackTitle << "--- TRACK ---\n";
      if (g_selOut.is_open()) g_selOut << trackTitle.str();
      std::cout << trackTitle.str();

      std::ostringstream trackHeader;
      trackHeader << std::right
                  << std::setw(8) << "NActive"
                  << std::setw(10) << "Chi2Dof"
                  << std::setw(9) << "P"
                  << std::setw(9) << "Pt"
                  << "\n";

      if (g_selOut.is_open()) g_selOut << trackHeader.str();
      std::cout << trackHeader.str();


      for (int i = 0; i < fNTracks[0]; ++i) {

        auto trk = fTrackBlock->Track(i);

        std::ostringstream line;
        line << std::right << std::fixed
             << std::setw(8) << trk->NActive()
             << std::setw(10) << std::setprecision(3) << trk->Chi2Dof()
             << std::setw(9) << std::setprecision(2) << trk->P()
             << std::setw(9) << std::setprecision(2) << trk->Pt()
             << "\n";

        if (g_selOut.is_open()) {
          g_selOut << line.str();
          g_selOut.flush();
        }
        std::cout << line.str();
      }
      std::cout << "\n";
      if (g_selOut.is_open()) g_selOut << "\n";

      // =========================================================
      // Dump HELIX INFO
      // =========================================================

      std::ostringstream helixTitle;
      helixTitle << "--- HELIX ---\n";

      std::ostringstream helixHeader;
      helixHeader << std::right
                  << std::setw(9) << "Chi2XY"
                  << std::setw(10) << "Chi2ZPhi"
                  << std::setw(7) << "NHits"
                  << std::setw(10) << "Radius"
                  << std::setw(10) << "Lambda"
                  << std::setw(9) << "TanDip"
                  << std::setw(9) << "pT"
                  << "\n";
      bool headerPrinted = false;


      for (int i = 0; i < fNHelices[0]; ++i) {

        auto helix = fHelixBlock->Helix(i);

        double chi2xy = helix->Chi2XY();
        double chi2zphi = helix->Chi2ZPhi();
        int nhits = helix->NHits();

        double mm2MeV = 3.0 / 10.0;
        double radius = helix->Radius();
        double lambda = helix->Lambda();
        double tanDip = lambda / radius;
        double pT = radius * mm2MeV;

        // Apply your quality cuts
        if (chi2xy > 5.0) continue;
        if (chi2zphi > 5.0) continue;

        if (!headerPrinted) {
          if (g_selOut.is_open()) {
            g_selOut << helixTitle.str();
            g_selOut << helixHeader.str();
          }
          std::cout << helixTitle.str() << helixHeader.str();
          headerPrinted = true;
        }

        std::ostringstream line;
        line << std::right << std::fixed
             << std::setw(9) << std::setprecision(3) << chi2xy
             << std::setw(10) << std::setprecision(3) << chi2zphi
             << std::setw(7) << nhits
             << std::setw(10) << std::setprecision(3) << radius
             << std::setw(10) << std::setprecision(3) << lambda
             << std::setw(9) << std::setprecision(3) << tanDip
             << std::setw(9) << std::setprecision(3) << pT
             << "\n";

        if (g_selOut.is_open()) {
          g_selOut << line.str();
          g_selOut.flush();
        }
        std::cout << line.str();
      }

      // =========================================================
      // Dump SIM PARTICLE INFO (simID and NHits)
      // =========================================================

      std::ostringstream simTitle;
      simTitle << "--- SIM PARTICLE ---\n";
      if (g_selOut.is_open()) g_selOut << simTitle.str();
      std::cout << simTitle.str();

      std::ostringstream simHeader;
      simHeader << std::right
                << std::setw(10) << "SimID"
                << std::setw(9) << "NHits"
                << "\n";

      if (g_selOut.is_open()) g_selOut << simHeader.str();
      std::cout << simHeader.str();

      for (const auto& d : pbarOriginParticles) {

        std::ostringstream line;
        line << std::right
             << std::setw(10) << d.simp_id
             << std::setw(9) << d.nhits
             << "\n";

        if (g_selOut.is_open()) {
          g_selOut << line.str();
          g_selOut.flush();
        }
        std::cout << line.str();
      }
      std::cout << "\n";
      if (g_selOut.is_open()) g_selOut << "\n";
    }

    // =========================================================
    // Dump events with fEffNSimp[1] == 2 and (quality-cut) N(helices) >= 3
    // =========================================================

    static std::ofstream g_selOut3;
    static bool g_filesOpened3 = false;

    auto ensure_out_open3 = []() {
      if (!g_filesOpened3) {
        g_selOut3.open("/exp/mu2e/app/users/kitagawa/20251126/selected_events_nhelix_ge3.txt",
                       std::ios::out | std::ios::app);
        if (!g_selOut3) {
          std::cerr << "[WARN] Can't open selected_events_nhelix_ge3.txt for append.\n";
        }
        g_filesOpened3 = true;
      }
    };

    auto dump_event_info = [&](std::ofstream& out, const std::string& extraHeader) {

      std::ostringstream tag3;
      tag3 << run_number << ":"
           << std::setw(2) << std::setfill('0') << subrun_number << ":"
           << std::setw(3) << std::setfill('0') << event_number;

      std::ostringstream header3;
      header3 << "\n============================================================\n";
      header3 << "Event: " << tag3.str() << "\n";
      header3 << "Run: " << run_number << "  Subrun: " << subrun_number << "  Event: " << event_number << "\n";
      header3 << "Tracks (fNTracks[0]):                         " << fNTracks[0] << "\n";
      header3 << "Helices (fNHelices[0]):                       " << fNHelices[0] << "\n";
      header3 << "Helices passing quality cuts (fNHelices[1]):  " << fNHelices[1] << "\n";
      header3 << extraHeader;
      header3 << "============================================================\n\n";

      if (out.is_open()) {
        out << header3.str();
        out.flush();
      }
      std::cout << header3.str();

      // ---------------------------------------------------------
      // TimeCluster info: [0] TimeClusterBlock (TZClusterFinder),
      // [1] PhiClusterBlock (PhiClusterFinder).
      // T0 = t-vs-z line-fit intercept at z=0 (tracker center) [ns].
      // SimID/PDG/NHSim = MC particle making most straw hits in the cluster.
      // ---------------------------------------------------------
      const char* tclNames3[2] = {"TimeClusterBlock", "PhiClusterBlock"};
      for (int ib = 0; ib < 2; ++ib) {
        std::ostringstream tclTitle3;
        tclTitle3 << "--- TIME CLUSTER: " << tclNames3[ib]
                  << " (N=" << fNTimeClusters[ib] << ") ---\n";

        std::ostringstream tclHeader3;
        tclHeader3 << std::right
                   << std::setw(4) << "Idx"
                   << std::setw(7) << "NHits"
                   << std::setw(7) << "NCH"
                   << std::setw(10) << "T0"
                   << std::setw(9) << "T0Err"
                   << std::setw(9) << "PosX"
                   << std::setw(9) << "PosY"
                   << std::setw(9) << "PosZ"
                   << std::setw(7) << "HelIdx"
                   << std::setw(8) << "SimID"
                   << std::setw(8) << "PDG"
                   << std::setw(7) << "NHSim"
                   << std::setw(9) << "McMom"
                   << std::setw(11) << "ClusterT"
                   << "\n";

        if (out.is_open()) out << tclTitle3.str() << tclHeader3.str();
        std::cout << tclTitle3.str() << tclHeader3.str();

        for (int i = 0; i < fNTimeClusters[ib]; ++i) {
          TStnTimeCluster* tc = fTimeClusterBlock[ib]->TimeCluster(i);
          if (!tc) continue;

          std::ostringstream line;
          line << std::right << std::fixed
               << std::setw(4) << i
               << std::setw(7) << tc->NHits()
               << std::setw(7) << tc->NComboHits()
               << std::setw(10) << std::setprecision(3) << tc->T0()
               << std::setw(9) << std::setprecision(3) << tc->T0Err()
               << std::setw(9) << std::setprecision(1) << tc->posX()
               << std::setw(9) << std::setprecision(1) << tc->posY()
               << std::setw(9) << std::setprecision(1) << tc->posZ()
               << std::setw(7) << tc->HelixSeedIndex()
               << std::setw(8) << tc->fSimID
               << std::setw(8) << tc->fPdgID
               << std::setw(7) << tc->fNHitsSimID
               << std::setw(9) << std::setprecision(2) << tc->fMcMom
               << std::setw(11) << std::setprecision(3) << tc->ClusterTime()
               << "\n";

          if (out.is_open()) {
            out << line.str();
            out.flush();
          }
          std::cout << line.str();
        }
        std::cout << "\n";
        if (out.is_open()) out << "\n";
      }

      // ---------------------------------------------------------
      // MC particle info: particles with >= 10 straw hits
      // "Sel" flags the (<=2) pbar-origin particles with >= 20 hits
      // that fEffNSimp[1] was computed from.
      // ---------------------------------------------------------
      std::vector<int> selectedPbarIDs3;
      for (const auto& d : pbarOriginParticles) {
        if (d.nhits >= 20) selectedPbarIDs3.push_back(d.simp_id);
      }

      std::ostringstream simTitle3;
      simTitle3 << "--- MC PARTICLES (NStrawHits >= 10) ---\n";
      simTitle3 << "    TCl0/TCl1: time clusters of TimeClusterBlock/PhiClusterBlock in which this particle\n";
      simTitle3 << "    is the dominant MC particle, as idx(n) with n = N(straw hits) from this particle.\n";
      simTitle3 << "    Hits of a particle in clusters where it is NOT dominant are not stored in the stntuple.\n";

      std::ostringstream simHeader3;
      simHeader3 << std::right
                 << std::setw(10) << "SimID"
                 << std::setw(10) << "ParentID"
                 << std::setw(8) << "PDG"
                 << std::setw(8) << "NHits"
                 << std::setw(10) << "P"
                 << std::setw(10) << "Pt"
                 << std::setw(10) << "Px"
                 << std::setw(10) << "Py"
                 << std::setw(10) << "Pz"
                 << std::setw(12) << "Tstart"
                 << std::setw(12) << "Tend"
                 << std::setw(5) << "Sel"
                 << std::setw(18) << "TCl0"
                 << std::setw(18) << "TCl1"
                 << "\n";

      if (out.is_open()) out << simTitle3.str() << simHeader3.str();
      std::cout << simTitle3.str() << simHeader3.str();

      for (int i = 0; i < fNSimp; ++i) {
        TSimParticle* s = fSimpBlock->Particle(i);
        if (!s) continue;

        int nhits = s->NStrawHits();
        if (nhits < 10) continue;

        const TLorentzVector* mom0 = s->StartMom();
        double px = 0, py = 0, pz = 0, p = 0, pt = 0;
        if (mom0) {
          TVector3 mv = mom0->Vect();
          px = mv.X();
          py = mv.Y();
          pz = mv.Z();
          p = mv.Mag();
          pt = mv.Pt();
        }

        const TLorentzVector* pos0 = s->StartPos();
        const TLorentzVector* pos1 = s->EndPos();
        double tStart = pos0 ? pos0->T() : 0.0;
        double tEnd = pos1 ? pos1->T() : 0.0;

        int sid = s->GetUniqueID();
        bool sel = std::find(selectedPbarIDs3.begin(), selectedPbarIDs3.end(), sid) != selectedPbarIDs3.end();

        // time clusters in which this particle is the dominant MC particle
        std::string tclList[2];
        for (int ib = 0; ib < 2; ++ib) {
          for (int itc = 0; itc < fNTimeClusters[ib]; ++itc) {
            TStnTimeCluster* tc = fTimeClusterBlock[ib]->TimeCluster(itc);
            if (!tc) continue;
            if (tc->fSimID != sid) continue;
            if (!tclList[ib].empty()) tclList[ib] += ",";
            tclList[ib] += std::to_string(itc) + "(" + std::to_string(tc->fNHitsSimID) + ")";
          }
          if (tclList[ib].empty()) tclList[ib] = "-";
        }

        std::ostringstream line;
        line << std::right << std::fixed
             << std::setw(10) << sid
             << std::setw(10) << s->ParentID()
             << std::setw(8) << s->PDGCode()
             << std::setw(8) << nhits
             << std::setw(10) << std::setprecision(2) << p
             << std::setw(10) << std::setprecision(2) << pt
             << std::setw(10) << std::setprecision(2) << px
             << std::setw(10) << std::setprecision(2) << py
             << std::setw(10) << std::setprecision(2) << pz
             << std::setw(12) << std::setprecision(3) << tStart
             << std::setw(12) << std::setprecision(3) << tEnd
             << std::setw(5) << (sel ? "Y" : "")
             << std::setw(18) << tclList[0]
             << std::setw(18) << tclList[1]
             << "\n";

        if (out.is_open()) {
          out << line.str();
          out.flush();
        }
        std::cout << line.str();
      }
      std::cout << "\n";
      if (out.is_open()) out << "\n";

      // ---------------------------------------------------------
      // HELIX info: all reconstructed helices
      // ---------------------------------------------------------
      std::ostringstream helixTitle3;
      helixTitle3 << "--- HELIX (N=" << fNHelices[0] << ") ---\n";

      std::ostringstream helixHeader3;
      helixHeader3 << std::right
                   << std::setw(9) << "Chi2XY"
                   << std::setw(10) << "Chi2ZPhi"
                   << std::setw(7) << "NHits"
                   << std::setw(10) << "Radius"
                   << std::setw(10) << "Lambda"
                   << std::setw(9) << "TanDip"
                   << std::setw(9) << "P"
                   << std::setw(9) << "pT"
                   << std::setw(5) << "Hel"
                   << std::setw(10) << "T0"
                   << std::setw(10) << "T0Err"
                   << std::setw(11) << "ClusterT"
                   << std::setw(7) << "TclIdx"
                   << std::setw(8) << "SimID1"
                   << std::setw(8) << "PDG1"
                   << std::setw(6) << "NCH1"
                   << std::setw(11) << "Tstart1"
                   << "\n";

      if (out.is_open()) out << helixTitle3.str() << helixHeader3.str();
      std::cout << helixTitle3.str() << helixHeader3.str();

      for (int i = 0; i < fNHelices[0]; ++i) {
        auto helix = fHelixBlock->Helix(i);

        double chi2xy = helix->Chi2XY();
        double chi2zphi = helix->Chi2ZPhi();
        int nhits = helix->NHits();
        double radius = helix->Radius();
        double lambda = helix->Lambda();
        double tanDip = (radius != 0) ? lambda / radius : 0.0;
        double p = helix->P();
        double pT = helix->Pt();
        int hel = helix->Helicity();
        double t0 = helix->T0();
        double t0Err = helix->T0Err();
        double clusterT = helix->ClusterTime();

        std::ostringstream line;
        line << std::right << std::fixed
             << std::setw(9) << std::setprecision(3) << chi2xy
             << std::setw(10) << std::setprecision(3) << chi2zphi
             << std::setw(7) << nhits
             << std::setw(10) << std::setprecision(3) << radius
             << std::setw(10) << std::setprecision(3) << lambda
             << std::setw(9) << std::setprecision(3) << tanDip
             << std::setw(9) << std::setprecision(3) << p
             << std::setw(9) << std::setprecision(3) << pT
             << std::setw(5) << hel
             << std::setw(10) << std::setprecision(3) << t0
             << std::setw(10) << std::setprecision(3) << t0Err
             << std::setw(11) << std::setprecision(3) << clusterT
             << std::setw(7) << helix->TimeClusterIndex()
             << std::setw(8) << helix->fSimpID1
             << std::setw(8) << helix->PDG1()
             << std::setw(6) << helix->ComboHitsFrom1()
             << std::setw(11) << std::setprecision(3) << helix->SimpOrigin1().T()
             << "\n";

        if (out.is_open()) {
          out << line.str();
          out.flush();
        }
        std::cout << line.str();
      }
      std::cout << "\n";
      if (out.is_open()) out << "\n";

      // ---------------------------------------------------------
      // TRACK info: all reconstructed tracks
      // ---------------------------------------------------------
      std::ostringstream trackTitle3;
      trackTitle3 << "--- TRACK (N=" << fNTracks[0] << ") ---\n";

      std::ostringstream trackHeader3;
      trackHeader3 << std::right
                   << std::setw(8) << "NActive"
                   << std::setw(10) << "Chi2Dof"
                   << std::setw(9) << "P"
                   << std::setw(9) << "Pt"
                   << std::setw(9) << "D0"
                   << std::setw(10) << "T0"
                   << std::setw(10) << "T0Err"
                   << "\n";

      if (out.is_open()) out << trackTitle3.str() << trackHeader3.str();
      std::cout << trackTitle3.str() << trackHeader3.str();

      for (int i = 0; i < fNTracks[0]; ++i) {
        auto trk = fTrackBlock->Track(i);

        std::ostringstream line;
        line << std::right << std::fixed
             << std::setw(8) << trk->NActive()
             << std::setw(10) << std::setprecision(3) << trk->Chi2Dof()
             << std::setw(9) << std::setprecision(2) << trk->P()
             << std::setw(9) << std::setprecision(2) << trk->Pt()
             << std::setw(9) << std::setprecision(2) << trk->D0()
             << std::setw(10) << std::setprecision(3) << trk->T0()
             << std::setw(10) << std::setprecision(3) << trk->T0Err()
             << "\n";

        if (out.is_open()) {
          out << line.str();
          out.flush();
        }
        std::cout << line.str();
      }
      std::cout << "\n";
      if (out.is_open()) out << "\n";
    };

    if (fEffNSimp[1] == 2 and fNHelices[1] >= 3) {
      ensure_out_open3();
      dump_event_info(g_selOut3, "");
    }

    // =========================================================
    // Dump events having a helix pair with deltaT0 < deltaT0DumpThre
    // (same selection as the deltaT0 histogram: fEffNSimp[1] == 2 and
    // both helices passing Chi2XY/Chi2ZPhi <= 5). Same content as above.
    // =========================================================
    static std::ofstream g_selOutDt;
    static bool g_filesOpenedDt = false;

    if (deltaT0BelowThre) {
      if (!g_filesOpenedDt) {
        g_selOutDt.open("/exp/mu2e/app/users/kitagawa/20251126/selected_events_deltat0_below_m400.txt",
                        std::ios::out | std::ios::app);
        if (!g_selOutDt) {
          std::cerr << "[WARN] Can't open selected_events_deltat0_below_m400.txt for append.\n";
        }
        g_filesOpenedDt = true;
      }

      std::ostringstream extra;
      extra << std::fixed << std::setprecision(3);
      extra << "Min deltaT0 = T0_i - T0_j:                     " << deltaT0Min
            << " ns (helix i=" << deltaT0MinI << ", j=" << deltaT0MinJ << ")\n";
      extra << "Selection: deltaT0 < " << deltaT0DumpThre << " ns\n";

      dump_event_info(g_selOutDt, extra.str());
    }
  }


  //-----------------------------------------------------------------------------
  // 2014-04-30: it looks that reading the straw hits takes a lot of time -
  //              turn off by default by commenting it out
  //-----------------------------------------------------------------------------
  int THelixAnaModule::Event(int ientry) {

    //  TEmuLogLH::PidData_t  dat;
    //  TStnTrack*            track;
    //  int                   id_word;
    TLorentzVector mom;

    //  TDiskCalorimeter::GeomData_t disk_geom;

    fHelixBlock->GetEntry(ientry);

    fGenpBlock->GetEntry(ientry);
    fSimpBlock->GetEntry(ientry);
    fTrackBlock->GetEntry(ientry);

    fTimeClusterBlock[0]->GetEntry(ientry);
    fTimeClusterBlock[1]->GetEntry(ientry);
    //-----------------------------------------------------------------------------
    // assume electron in the first particle, otherwise the logic will need to
    // be changed
    //-----------------------------------------------------------------------------
    fNGenp = fGenpBlock->NParticles();

    TGenParticle* genp;
    int pdg_code, process_code;

    fParticle = NULL;
    for (int i = fNGenp - 1; i >= 0; i--) {
      genp = fGenpBlock->Particle(i);
      pdg_code = genp->GetPdgCode();
      process_code = genp->GetStatusCode();
      if ((abs(pdg_code) == fPdgCode) && (process_code == fProcessCode)) {
        fParticle = genp;
        break;
      }
    }
    // may want to revisit the definition of fSimp
    fSimp = fSimpBlock->Particle(0);
    fNSimp = fSimpBlock->NParticles();
    fNTracks[0] = fTrackBlock->NTracks();

    double p = 0;

    if (fParticle) {
      fParticle->Momentum(mom);
      // this is a kludge, to be removed at the next
      // ntupling
      p = mom.P();
      fEleE = sqrt(p * p + 0.511 * 0.511);
    }

    fNTimeClusters[0] = fTimeClusterBlock[0]->NTimeClusters();
    fNTimeClusters[1] = fTimeClusterBlock[1]->NTimeClusters();

    fNHelices[0] = fHelixBlock->NHelices();

    fNHelPos = 0;
    fNHelNeg = 0;
    for (int i = 0; i < fNHelices[0]; i++) {
      TStnHelix* h = fHelixBlock->Helix(i);
      if (h->Helicity() > 0) fNHelPos++;
      else if (h->Helicity() < 0) fNHelNeg++;
      else {
        printf("ERROR 001 in THelixAnaModule: helix with helicity=0\n");
      }
    }

    FillHistograms();

    Debug();

    return 0;
  }

  //-----------------------------------------------------------------------------
  void THelixAnaModule::Debug() {

    if (GetDebugBit(3) == 1) {
      if (fNHelPos > 0) {
        GetHeaderBlock()->Print(Form("N(CalHelixFinder helices hel > 0) = %2i", fNHelPos));
      }
    }

    if (GetDebugBit(4) == 1) {
      if (fNHelNeg > 0) {
        GetHeaderBlock()->Print(Form("N(CalHelixFinder helices hel < 0) = %2i", fNHelNeg));
      }
    }

    if (GetDebugBit(5) == 1) {
      if (fNHelNeg + fNHelPos > 1) {
        GetHeaderBlock()->Print(Form("N(helices = 0) = %2i", fNHelNeg + fNHelPos));
      }
    }
  }

  //_____________________________________________________________________________
  int THelixAnaModule::EndJob() {
    return 0;
  }

  //_____________________________________________________________________________
  void THelixAnaModule::Test001() {

    // mu2e::HexMap* hmap      = new mu2e::HexMap();

    // mu2e::HexLK hex_index(0,0);

    // for (int i=0; i<40; i++) {
    //   hex_index = hmap->lk(i);
    //   printf(" i,l,k = %5i %5i %5i\n",i,hex_index._l,hex_index._k);
    // }
  }
}
