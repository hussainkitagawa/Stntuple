//
// Read the track 
//
// $Id:   $
// $Author:$
// $Date: $
//
// Contact person Pavel Murat
//
#ifndef Stntuple_obj_TStnHelix_hh
#define Stntuple_obj_TStnHelix_hh

// storable objects (data products)
// 

// C++ includes.
#include <iostream>

#include "TString.h"
#include "TFolder.h"
#include "TFile.h"
#include "TBuffer.h"
#include "TLorentzVector.h"

//namespace murat {

namespace mu2e {
  class HelixSeed;
  class StrawHit;
  class CaloCluster;
}

class TStnHelix : public TObject {

  enum {
    kNFreeIntsV2   = 10,		// V2
    kNFreeFloatsV2 = 10,  		// V2

    kNFreeIntsV3   =  4,	        // v3: added the indices of the two SimParticles contributing 
    kNFreeFloatsV3 = 10,                //     the most and their fraction of hits within the helix, 
   			                //     also added the p and pT of the SimParticles.

    kNFreeIntsV4   =  3,	        // v4: added helicity
    kNFreeFloatsV4 = 10,                //

    kNFreeIntsV5   =  3,	        // v5: added number of loops
    kNFreeFloatsV5 =  9,                //

    kNFreeIntsV6   =  3,	        // v6: added TZSlope, TZSlope error, TZChi2NDof, hitRatio (expected/collected)
    kNFreeFloatsV6 =  5,		// 

    kNFreeInts     =  1,	        // v7: added propDir, sim 1 ID
    kNFreeFloats   =  5			// 
  };

public:
  enum {
    kMaxMcParticles = 50		// v8: max N(MC particles) stored per helix
  };
//-----------------------------------------------------------------------------
// vectors, added in V3
//-----------------------------------------------------------------------------
  TLorentzVector            fSimpMom1;  
  TLorentzVector            fSimpOrigin1;
  TLorentzVector            fSimpMom2;  
  TLorentzVector            fSimpOrigin2;
//-----------------------------------------------------------------------------
// integers
//-----------------------------------------------------------------------------
  int                       fNHits;	     
  int                       fAlgorithmID;     // bit-packed : (alg_mask << 16 ) | best
  int                       fTimeClusterIndex;
  int                       fTrackSeedIndex; 
  int                       fNComboHits;     
  int                       fSimpPDG1;          // added in v3
  int                       fSimpPDGM1;         // added in v3
  int                       fSimpId1Hits;       // added in v3
  int                       fSimpPDG2;          // added in v3
  int                       fSimpPDGM2;         // added in v3
  int                       fSimpId2Hits;       // added in v3
  int                       fHelicity;          // added in v4
  int                       fPropDir;           // added in v7: -1,0,1 stand for upstream, ambigous and downstream respectively
  int                       fSimpID1;           // added in v7: ID for sim with the most hits
  int                       fInt[kNFreeInts];   // provision for future I/O expansion
//-----------------------------------------------------------------------------
// floats
//-----------------------------------------------------------------------------
  float                     fT0;    	  
  float                     fT0Err; 	  
			                  
  float                     fRCent;  // radius of circle center
  float                     fFCent;  // azimuth of circle center
  float                     fRadius; // transverse radius of the helix (mm).  Always positive
  float                     fLambda; // dz/dphi (mm/radian)
  float                     fFZ0;    // azimuth (phi) at the center z position (radians)
  float                     fD0;     // impact paramter 
  float                     fChi2XYNDof;  
  float                     fChi2PhiZNDof;
			                  
  float                     fClusterTime;   
  float			    fClusterEnergy; 
  float			    fClusterX;      
  float			    fClusterY;      
  float			    fClusterZ;      
  float                     fNLoops;
  float                     fTZSlope;
  float                     fTZSlopeError;
  float                     fChi2TZNDof;
  float                     fHitRatio; // expected hits in the tracker/hits_collected
  float                     fFloat[kNFreeFloats]; // provision for future I/O expansion
//-----------------------------------------------------------------------------
// transients
//-----------------------------------------------------------------------------
  const mu2e::HelixSeed*    fHelix;  //!
  int                       fNumber; //! sequential number in the list, set by TStnHelixBlock::Streamer
//-----------------------------------------------------------------------------
// added in v8: all MC particles contributing (non-outlier) combo hits to the
// helix, sorted by N(combo hits) in descending order, at most
// kMaxMcParticles stored. Written explicitly by the Streamer (variable length
// on disk); must stay after fHelix, the fixed part is streamed by pointer arithmetic.
//-----------------------------------------------------------------------------
  int                       fNMcParticles;              // N(distinct MC particles), -1: not available (v1-v7)
  int                       fMcSimID[kMaxMcParticles];  // SimParticle ID
  int                       fMcPdgID[kMaxMcParticles];  // PDG code
  int                       fMcNHits[kMaxMcParticles];  // N(combo hits) assigned to the particle (majority of straw hits)
  float                     fMcShare[kMaxMcParticles];  // share: fMcNHits / sum of fMcNHits of this particle over all helices of the block
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
  TStnHelix(int Number = -1);
  ~TStnHelix();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int     NHits           () { return fNHits;       }
  int     Number          () { return fNumber;      }
  int     Helicity        () { return fHelicity;    }
  int     NComboHits      () { return fNComboHits;  }
  int     AlgorithmID     () { return fAlgorithmID; }
  int     AlgMask         () { return (fAlgorithmID >> 16) & 0xffff; }
  int     BestAlg         () { return fAlgorithmID & 0xffff; }
  int     TimeClusterIndex() { return fTimeClusterIndex; }
  int     TrackSeedIndex  () { return fTrackSeedIndex;   }
  int     PDG1            () { return fSimpPDG1; }
  int     PDGMother1      () { return fSimpPDGM1; }
  int     ComboHitsFrom1  () { return fSimpId1Hits; }
  int     PDG2            () { return fSimpPDG2; }
  int     PDGMother2      () { return fSimpPDGM2; }
  int     ComboHitsFrom2  () { return fSimpId2Hits; }
  int     PropDir         () { return fPropDir;}
  int     SimpID1         () { return fSimpID1;}

  int     NMcParticles      () const { return fNMcParticles; }
  int     NStoredMcParticles() const {
    return (fNMcParticles < 0) ? 0 : ((fNMcParticles < kMaxMcParticles) ? fNMcParticles : kMaxMcParticles);
  }
  int     McSimID       (int I) const { return fMcSimID[I]; }
  int     McPdgID       (int I) const { return fMcPdgID[I]; }
  int     McNHits       (int I) const { return fMcNHits[I]; }
  float   McShare       (int I) const { return fMcShare[I]; }
                                        // purity: fraction of MC-matched combo hits from the dominant
                                        // particle (index 0); -1 if not available
  float   McPurity      () const;
                                        // N(combo hits) from SimID in this helix:
                                        // 0 if not found, -1 if not available (v1-v7)
  int     NHitsFromSimID(int SimID) const;

  float   T0         () { return  fT0;     }
  float   T0Err      () { return  fT0Err;  }

  float   RCent      () { return  fRCent;  }
  float   FCent      () { return  fFCent;  }
  float   Radius     () { return  fRadius; }
  float   Lambda     () { return  fLambda; }
  float   FZ0        () { return  fFZ0;    }
  float   CenterX    () { return  fRCent*cos(fFCent);}
  float   CenterY    () { return  fRCent*sin(fFCent);}
  float   D0         () { return  fD0;     }

  float   Pz         () { return  fLambda*3./10;}//assumes 1T magnetic field!
  float   Pt         () { return  fRadius*3./10;}//assumes 1T magnetic field!
  float   P          () { return  sqrt(fRadius*fRadius + fLambda*fLambda)*3./10;}//assumes 1T magnetic field!

  float   Chi2XY     () { return  fChi2XYNDof;}
  float   Chi2ZPhi   () { return  fChi2PhiZNDof;}

  float   ClusterTime   () { return fClusterTime;  }
  float   ClusterEnergy () { return fClusterEnergy;}
  float   ClusterX      () { return fClusterX;     }
  float   ClusterY      () { return fClusterY;     }
  float   ClusterZ      () { return fClusterZ;     }

  float   NLoops        () { return fNLoops;       }

  float   TZSlope       () { return fTZSlope;      }
  float   TZSlopeError  () { return fTZSlopeError; }
  float   TZSlopeSig    () { return std::fabs(fTZSlope/fTZSlopeError);      }
  float   Chi2TZNDof    () { return fChi2TZNDof;   }
  float   HitRatio      () { return fHitRatio;     }
  
  TLorentzVector  SimpMom1     () { return fSimpMom1; }
  TLorentzVector  SimpOrigin1  () { return fSimpOrigin1; }
  TLorentzVector  SimpMom2     () { return fSimpMom2; }
  TLorentzVector  SimpOrigin2  () { return fSimpOrigin2; }

//----------------------------------------------------------------------------
// setters
//----------------------------------------------------------------------------
  void    SetTimeClusterIndex(int I) { fTimeClusterIndex = I; }
  void    SetTrackSeedIndex  (int I) { fTrackSeedIndex   = I; }
  void    SetNumber          (int I) { fNumber           = I; }
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void Clear(Option_t* Opt = "") ;
  virtual void Print(Option_t* Opt = "") const ;
//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
  void ReadV1(TBuffer& R__b);
  void ReadV2(TBuffer& R__b);
  void ReadV3(TBuffer& R__b);   // 2018-12-05 P.M.
  void ReadV4(TBuffer& R__b);   // 2019-02-27 G.P.
  void ReadV5(TBuffer& R__b);   // 2024-03-07 G.P.
  void ReadV6(TBuffer& R__b);   // 2024-10-17 G.P.
  void ReadV7(TBuffer& R__b);

  ClassDef(TStnHelix,8);
};

#endif
