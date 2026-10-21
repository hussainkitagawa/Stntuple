//
#include <iostream>

#include "TMatrix.h"
#include "obj/TStnHelix.hh"

// #include "RecoDataProducts/inc/Helix.hh"

namespace {
  //  const double BF = 1.4116 ;  // CDF case
  const double BF = 1.0 ;         // MU2E case
}

//namespace murat {

ClassImp(TStnHelix)
//-----------------------------------------------------------------------------
void TStnHelix::ReadV1(TBuffer &R__b) {

  struct TStnHelixDataV1_t {
//-----------------------------------------------------------------------------
// integers
//-----------------------------------------------------------------------------
    int                       fNHits;
    int                       fAlgorithmID;     // bit-packed : (alg_mask << 16 ) | best
    int                       fTimeClusterIndex;
    int                       fTrackSeedIndex;
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
    float		      fClusterEnergy; 
    float                     fClusterX;      
    float		      fClusterY;      
    float		      fClusterZ;      
//-----------------------------------------------------------------------------
// transients
//-----------------------------------------------------------------------------
    const mu2e::HelixSeed*    fHelix;  //!
 };

  TStnHelixDataV1_t data;

  int            nwi, nwf;
  
  nwi      = ((int*  ) &data.fT0            ) - &data.fNHits;
  nwf      = ((float*) &data.fHelix         ) - &data.fT0;
    
  R__b.ReadFastArray(&data.fNHits ,nwi);
  R__b.ReadFastArray(&data.fT0    ,nwf);

  fSimpMom1.SetPxPyPzE(0.,0.,0.,0.);
  fSimpOrigin1.SetXYZT(0.,0.,0.,0.);

  fSimpMom2.SetPxPyPzE(0.,0.,0.,0.);
  fSimpOrigin2.SetXYZT(0.,0.,0.,0.);


  fNHits             = data.fNHits;	     
  fAlgorithmID       = data.fAlgorithmID;    
  fTimeClusterIndex  = data.fTimeClusterIndex;
  fTrackSeedIndex    = data.fTrackSeedIndex; 
  fSimpPDG1          = -1;         // added in v3
  fSimpPDGM1         = -1;         // added in v3
  fSimpId1Hits       = -1;         // added in v3
  fSimpPDG2          = -1;         // added in v3
  fSimpPDGM2         = -1;         // added in v3
  fSimpId2Hits       = -1;         // added in v3
  fHelicity          =  0;         // added in V4
  fPropDir           =  0;         // added in V7
  fSimpID1           = -1;         // added in v7


  fT0                = data.fT0;    	  
  fT0Err             = data.fT0Err; 	  
                
  fRCent             = data.fRCent;
  fFCent             = data.fFCent;
  fRadius            = data.fRadius;
  fLambda            = data.fLambda;
  fFZ0               = data.fFZ0;   
  fD0                = data.fD0;    
  fChi2XYNDof        = data.fChi2XYNDof;  
  fChi2PhiZNDof      = data.fChi2PhiZNDof;
                
  fClusterTime       = data.fClusterTime; 
  fClusterEnergy     = data.fClusterEnergy; 
  fClusterX          = data.fClusterX;      
  fClusterY          = data.fClusterY;    
  fClusterZ          = data.fClusterZ;    
}

//-----------------------------------------------------------------------------
void TStnHelix::ReadV2(TBuffer &R__b) {

  struct TStnHelixDataV2_t {
//-----------------------------------------------------------------------------
// integers
//-----------------------------------------------------------------------------
    int                       fNHits;
    int                       fAlgorithmID;     // bit-packed : (alg_mask << 16 ) | best
    int                       fTimeClusterIndex;
    int                       fTrackSeedIndex;
    int                       fInt[kNFreeIntsV2];
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
    float		      fClusterEnergy; 
    float                     fClusterX;      
    float		      fClusterY;      
    float		      fClusterZ;      
    float                     fFloat[kNFreeFloatsV2]; // provision for future I/O expansion
//-----------------------------------------------------------------------------
// transients
//-----------------------------------------------------------------------------
    const mu2e::HelixSeed*    fHelix;  //!
 };

  TStnHelixDataV2_t data;

  int            nwi, nwf;
  
  nwi      = ((int*  ) &data.fT0            ) - &data.fNHits;
  nwf      = ((float*) &data.fHelix         ) - &data.fT0;
    
  R__b.ReadFastArray(&data.fNHits ,nwi);
  R__b.ReadFastArray(&data.fT0    ,nwf);

  fSimpMom1.SetPxPyPzE(0.,0.,0.,0.);
  fSimpOrigin1.SetXYZT(0.,0.,0.,0.);

  fSimpMom2.SetPxPyPzE(0.,0.,0.,0.);
  fSimpOrigin2.SetXYZT(0.,0.,0.,0.);

  fNHits             = data.fNHits;	     
  fAlgorithmID       = data.fAlgorithmID;    
  fTimeClusterIndex  = data.fTimeClusterIndex;
  fTrackSeedIndex    = data.fTrackSeedIndex; 
  fSimpPDG1          = -1;         // added in v3
  fSimpPDGM1         = -1;         // added in v3
  fSimpId1Hits       = -1;         // added in v3
  fSimpPDG2          = -1;         // added in v3
  fSimpPDGM2         = -1;         // added in v3
  fSimpId2Hits       = -1;         // added in v3
  fHelicity          =  0;         // added in V4
  fPropDir           =  0;         // added in V7
  fSimpID1           = -1;         // added in v7

  fT0                = data.fT0;    	  
  fT0Err             = data.fT0Err; 	  
                
  fRCent             = data.fRCent;
  fFCent             = data.fFCent;
  fRadius            = data.fRadius;
  fLambda            = data.fLambda;
  fFZ0               = data.fFZ0;   
  fD0                = data.fD0;    
  fChi2XYNDof        = data.fChi2XYNDof;  
  fChi2PhiZNDof      = data.fChi2PhiZNDof;
                
  fClusterTime       = data.fClusterTime; 
  fClusterEnergy     = data.fClusterEnergy; 
  fClusterX          = data.fClusterX;      
  fClusterY          = data.fClusterY;    
  fClusterZ          = data.fClusterZ;    
  
}

//-----------------------------------------------------------------------------
void TStnHelix::ReadV3(TBuffer &R__b) {

  struct TStnHelixDataV3_t {
//-----------------------------------------------------------------------------
// vectors, added in V3, I believe
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
    int                       fInt[kNFreeIntsV3];
//-----------------------------------------------------------------------------
// floats
//-----------------------------------------------------------------------------
    float                     fT0;    
    float                     fT0Err; 

    float                     fRCent;   // radius of circle center
    float                     fFCent;   // azimuth of circle center
    float                     fRadius;  // transverse radius of the helix (mm).  Always positive
    float                     fLambda;  // dz/dphi (mm/radian)
    float                     fFZ0;     // azimuth (phi) at the center z position (radians)
    float                     fD0;      // impact paramter 
    float                     fChi2XYNDof;
    float                     fChi2PhiZNDof;

    float                     fClusterTime;   
    float		      fClusterEnergy; 
    float                     fClusterX;      
    float		      fClusterY;      
    float		      fClusterZ;      
    float                     fFloat[kNFreeFloatsV3]; // provision for future I/O expansion
//-----------------------------------------------------------------------------
// transients
//-----------------------------------------------------------------------------
    const mu2e::HelixSeed*    fHelix;  //!
 };

  TStnHelixDataV3_t data;

  int            nwi, nwf;
  
  nwi      = ((int*  ) &data.fT0            ) - &data.fNHits;
  nwf      = ((float*) &data.fHelix         ) - &data.fT0;
    
  fSimpMom1.Streamer(R__b);
  fSimpOrigin1.Streamer(R__b);
  fSimpMom2.Streamer(R__b);
  fSimpOrigin2.Streamer(R__b);

  R__b.ReadFastArray(&data.fNHits ,nwi);
  R__b.ReadFastArray(&data.fT0    ,nwf);

  fNHits             = data.fNHits;	     
  fAlgorithmID       = data.fAlgorithmID;    
  fTimeClusterIndex  = data.fTimeClusterIndex;
  fTrackSeedIndex    = data.fTrackSeedIndex;
  fNComboHits        = data.fNComboHits ;
  fSimpPDG1          = data.fSimpPDG1   ;
  fSimpPDGM1         = data.fSimpPDGM1  ;
  fSimpId1Hits       = data.fSimpId1Hits;
  fSimpPDG2          = data.fSimpPDG2   ;
  fSimpPDGM2         = data.fSimpPDGM2  ;
  fSimpId2Hits       = data.fSimpId2Hits;
  fHelicity          =  0               ;         // added in V4
  fPropDir           =  0;         // added in V7
  fSimpID1           = -1               ;         // added in v7

  fT0                = data.fT0;    	  
  fT0Err             = data.fT0Err; 	  
                
  fRCent             = data.fRCent;
  fFCent             = data.fFCent;
  fRadius            = data.fRadius;
  fLambda            = data.fLambda;
  fFZ0               = data.fFZ0;   
  fD0                = data.fD0;    
  fChi2XYNDof        = data.fChi2XYNDof;  
  fChi2PhiZNDof      = data.fChi2PhiZNDof;
                
  fClusterTime       = data.fClusterTime; 
  fClusterEnergy     = data.fClusterEnergy; 
  fClusterX          = data.fClusterX;      
  fClusterY          = data.fClusterY;    
  fClusterZ          = data.fClusterZ;    
}

//-----------------------------------------------------------------------------
void TStnHelix::ReadV4(TBuffer &R__b) {

  struct TStnHelixDataV4_t {
//-----------------------------------------------------------------------------
// vectors, added in V3, I believe
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
    int                       fInt[kNFreeIntsV4];
//-----------------------------------------------------------------------------
// floats
//-----------------------------------------------------------------------------
    float                     fT0;    
    float                     fT0Err; 

    float                     fRCent;   // radius of circle center
    float                     fFCent;   // azimuth of circle center
    float                     fRadius;  // transverse radius of the helix (mm).  Always positive
    float                     fLambda;  // dz/dphi (mm/radian)
    float                     fFZ0;     // azimuth (phi) at the center z position (radians)
    float                     fD0;      // impact paramter 
    float                     fChi2XYNDof;
    float                     fChi2PhiZNDof;

    float                     fClusterTime;   
    float		      fClusterEnergy; 
    float                     fClusterX;      
    float		      fClusterY;      
    float		      fClusterZ;      
    float                     fFloat[kNFreeFloatsV4]; // provision for future I/O expansion
//-----------------------------------------------------------------------------
// transients
//-----------------------------------------------------------------------------
    const mu2e::HelixSeed*    fHelix;  //!
 };

  TStnHelixDataV4_t data;

  int            nwi, nwf;
  
  nwi      = ((int*  ) &data.fT0            ) - &data.fNHits;
  nwf      = ((float*) &data.fHelix         ) - &data.fT0;
    
  fSimpMom1.Streamer(R__b);
  fSimpOrigin1.Streamer(R__b);
  fSimpMom2.Streamer(R__b);
  fSimpOrigin2.Streamer(R__b);

  R__b.ReadFastArray(&data.fNHits ,nwi);
  R__b.ReadFastArray(&data.fT0    ,nwf);

  fNHits             = data.fNHits;	     
  fAlgorithmID       = data.fAlgorithmID;    
  fTimeClusterIndex  = data.fTimeClusterIndex;
  fTrackSeedIndex    = data.fTrackSeedIndex;
  fNComboHits        = data.fNComboHits ;
  fSimpPDG1          = data.fSimpPDG1   ;
  fSimpPDGM1         = data.fSimpPDGM1  ;
  fSimpId1Hits       = data.fSimpId1Hits;
  fSimpPDG2          = data.fSimpPDG2   ;
  fSimpPDGM2         = data.fSimpPDGM2  ;
  fSimpId2Hits       = data.fSimpId2Hits;
  fHelicity          = data.fHelicity   ;         

  fT0                = data.fT0;    	  
  fT0Err             = data.fT0Err; 	  
                
  fRCent             = data.fRCent;
  fFCent             = data.fFCent;
  fRadius            = data.fRadius;
  fLambda            = data.fLambda;
  fFZ0               = data.fFZ0;   
  fD0                = data.fD0;    
  fChi2XYNDof        = data.fChi2XYNDof;  
  fChi2PhiZNDof      = data.fChi2PhiZNDof;
                
  fClusterTime       = data.fClusterTime; 
  fClusterEnergy     = data.fClusterEnergy; 
  fClusterX          = data.fClusterX;      
  fClusterY          = data.fClusterY;    
  fClusterZ          = data.fClusterZ;    
  fNLoops            = 0.;                  //added in V5
  fPropDir           =  0;                  // added in V7
  fSimpID1           = -1;                  // added in v7

}


//--------------------------------------------------------------------------------
void TStnHelix::ReadV5(TBuffer &R__b) {

  struct TStnHelixDataV5_t {
//-----------------------------------------------------------------------------
// vectors, added in V3, I believe
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
    int                       fInt[kNFreeIntsV5];
//-----------------------------------------------------------------------------
// floats
//-----------------------------------------------------------------------------
    float                     fT0;    
    float                     fT0Err; 

    float                     fRCent;   // radius of circle center
    float                     fFCent;   // azimuth of circle center
    float                     fRadius;  // transverse radius of the helix (mm).  Always positive
    float                     fLambda;  // dz/dphi (mm/radian)
    float                     fFZ0;     // azimuth (phi) at the center z position (radians)
    float                     fD0;      // impact paramter 
    float                     fChi2XYNDof;
    float                     fChi2PhiZNDof;

    float                     fClusterTime;   
    float		      fClusterEnergy; 
    float                     fClusterX;      
    float		      fClusterY;      
    float		      fClusterZ;      
    float                     fNLoops;
    float                     fFloat[kNFreeFloatsV5]; // provision for future I/O expansion
//-----------------------------------------------------------------------------
// transients
//-----------------------------------------------------------------------------
    const mu2e::HelixSeed*    fHelix;  //!
 };

  TStnHelixDataV5_t data;

  int            nwi, nwf;
  
  nwi      = ((int*  ) &data.fT0            ) - &data.fNHits;
  nwf      = ((float*) &data.fHelix         ) - &data.fT0;
    
  fSimpMom1.Streamer(R__b);
  fSimpOrigin1.Streamer(R__b);
  fSimpMom2.Streamer(R__b);
  fSimpOrigin2.Streamer(R__b);

  R__b.ReadFastArray(&data.fNHits ,nwi);
  R__b.ReadFastArray(&data.fT0    ,nwf);

  fNHits             = data.fNHits;	     
  fAlgorithmID       = data.fAlgorithmID;    
  fTimeClusterIndex  = data.fTimeClusterIndex;
  fTrackSeedIndex    = data.fTrackSeedIndex;
  fNComboHits        = data.fNComboHits ;
  fSimpPDG1          = data.fSimpPDG1   ;
  fSimpPDGM1         = data.fSimpPDGM1  ;
  fSimpId1Hits       = data.fSimpId1Hits;
  fSimpPDG2          = data.fSimpPDG2   ;
  fSimpPDGM2         = data.fSimpPDGM2  ;
  fSimpId2Hits       = data.fSimpId2Hits;
  fHelicity          = data.fHelicity   ;         
  fNLoops            = data.fNLoops;
  fPropDir           =  0;                  // added in V7
  fSimpID1           = -1;                  // added in v7

  fT0                = data.fT0;    	  
  fT0Err             = data.fT0Err; 	  
                
  fRCent             = data.fRCent;
  fFCent             = data.fFCent;
  fRadius            = data.fRadius;
  fLambda            = data.fLambda;
  fFZ0               = data.fFZ0;   
  fD0                = data.fD0;    
  fChi2XYNDof        = data.fChi2XYNDof;  
  fChi2PhiZNDof      = data.fChi2PhiZNDof;
                
  fClusterTime       = data.fClusterTime; 
  fClusterEnergy     = data.fClusterEnergy; 
  fClusterX          = data.fClusterX;      
  fClusterY          = data.fClusterY;    
  fClusterZ          = data.fClusterZ;    

  fTZSlope           = 0; //added in V6
  fTZSlopeError      = 0; //added in V6
  fChi2TZNDof        = 0; //added in V6
  fHitRatio          = 0; //added in V6
}


//--------------------------------------------------------------------------------
void TStnHelix::ReadV6(TBuffer &R__b) {

  struct TStnHelixDataV6_t {
//-----------------------------------------------------------------------------
// vectors, added in V3, I believe
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
    int                       fInt[kNFreeIntsV6];
//-----------------------------------------------------------------------------
// floats
//-----------------------------------------------------------------------------
    float                     fT0;    
    float                     fT0Err; 

    float                     fRCent;   // radius of circle center
    float                     fFCent;   // azimuth of circle center
    float                     fRadius;  // transverse radius of the helix (mm).  Always positive
    float                     fLambda;  // dz/dphi (mm/radian)
    float                     fFZ0;     // azimuth (phi) at the center z position (radians)
    float                     fD0;      // impact paramter 
    float                     fChi2XYNDof;
    float                     fChi2PhiZNDof;

    float                     fClusterTime;   
    float		      fClusterEnergy; 
    float                     fClusterX;      
    float		      fClusterY;      
    float		      fClusterZ;      
    float                     fNLoops;
    float                     fTZSlope;
    float                     fTZSlopeError;
    float                     fChi2TZNDof;
    float                     fHitRatio; 
    float                     fFloat[kNFreeFloatsV6]; // provision for future I/O expansion
//-----------------------------------------------------------------------------
// transients
//-----------------------------------------------------------------------------
    const mu2e::HelixSeed*    fHelix;  //!
 };

  TStnHelixDataV6_t data;

  int            nwi, nwf;
  
  nwi      = ((int*  ) &data.fT0            ) - &data.fNHits;
  nwf      = ((float*) &data.fHelix         ) - &data.fT0;
    
  fSimpMom1.Streamer(R__b);
  fSimpOrigin1.Streamer(R__b);
  fSimpMom2.Streamer(R__b);
  fSimpOrigin2.Streamer(R__b);

  R__b.ReadFastArray(&data.fNHits ,nwi);
  R__b.ReadFastArray(&data.fT0    ,nwf);

  fNHits             = data.fNHits;	     
  fAlgorithmID       = data.fAlgorithmID;    
  fTimeClusterIndex  = data.fTimeClusterIndex;
  fTrackSeedIndex    = data.fTrackSeedIndex;
  fNComboHits        = data.fNComboHits ;
  fSimpPDG1          = data.fSimpPDG1   ;
  fSimpPDGM1         = data.fSimpPDGM1  ;
  fSimpId1Hits       = data.fSimpId1Hits;
  fSimpPDG2          = data.fSimpPDG2   ;
  fSimpPDGM2         = data.fSimpPDGM2  ;
  fSimpId2Hits       = data.fSimpId2Hits;
  fHelicity          = data.fHelicity   ;         
  fNLoops            = data.fNLoops;
  fPropDir           =  0;                  // added in V7
  fSimpID1           = -1;                  // added in v7

  fT0                = data.fT0;    	  
  fT0Err             = data.fT0Err; 	  
                
  fRCent             = data.fRCent;
  fFCent             = data.fFCent;
  fRadius            = data.fRadius;
  fLambda            = data.fLambda;
  fFZ0               = data.fFZ0;   
  fD0                = data.fD0;    
  fChi2XYNDof        = data.fChi2XYNDof;  
  fChi2PhiZNDof      = data.fChi2PhiZNDof;
                
  fClusterTime       = data.fClusterTime; 
  fClusterEnergy     = data.fClusterEnergy; 
  fClusterX          = data.fClusterX;      
  fClusterY          = data.fClusterY;    
  fClusterZ          = data.fClusterZ;    

  fTZSlope           = data.fTZSlope; 
  fTZSlopeError      = data.fTZSlopeError; 
  fChi2TZNDof        = data.fChi2TZNDof;
  fHitRatio          = data.fHitRatio;

  fPropDir           = 0; //added in v7
}



//-----------------------------------------------------------------------------
// v7: same fixed part as v8, no MC particle list
//-----------------------------------------------------------------------------
void TStnHelix::ReadV7(TBuffer &R__b) {
  int nwi = ((int*  ) &fT0            ) - &fNHits;
  int nwf = ((float*) &fHelix         ) - &fT0;

  fSimpMom1.Streamer(R__b);
  fSimpOrigin1.Streamer(R__b);
  fSimpMom2.Streamer(R__b);
  fSimpOrigin2.Streamer(R__b);

  R__b.ReadFastArray(&fNHits, nwi);
  R__b.ReadFastArray(&fT0   , nwf);
}

//-----------------------------------------------------------------------------
void TStnHelix::Streamer(TBuffer& R__b) {
  int   nwi, nwf;
  
  nwi      = ((int*  ) &fT0            ) - &fNHits;
  nwf      = ((float*) &fHelix         ) - &fT0;

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); 

    fNMcParticles = -1;                 // added in v8, not available in older versions

    if      (R__v == 1) ReadV1 (R__b);
    else if (R__v == 2) ReadV2 (R__b);
    else if (R__v == 3) ReadV3 (R__b);
    else if (R__v == 4) ReadV4 (R__b);
    else if (R__v == 5) ReadV5 (R__b);
    else if (R__v == 6) ReadV6 (R__b);
    else if (R__v == 7) ReadV7 (R__b);
    else {
      fSimpMom1.Streamer(R__b);
      fSimpOrigin1.Streamer(R__b);
      fSimpMom2.Streamer(R__b);
      fSimpOrigin2.Streamer(R__b);
					// current version: V8
      R__b.ReadFastArray(&fNHits, nwi);
      R__b.ReadFastArray(&fT0   , nwf);

      R__b >> fNMcParticles;
      int ns = NStoredMcParticles();
      R__b.ReadFastArray(fMcSimID, ns);
      R__b.ReadFastArray(fMcPdgID, ns);
      R__b.ReadFastArray(fMcNHits, ns);
      R__b.ReadFastArray(fMcShare, ns);
    }
  }
  else {
    R__b.WriteVersion(TStnHelix::IsA());
    fSimpMom1.Streamer(R__b);
    fSimpOrigin1.Streamer(R__b);
    fSimpMom2.Streamer(R__b);
    fSimpOrigin2.Streamer(R__b);

    R__b.WriteFastArray(&fNHits, nwi);
    R__b.WriteFastArray(&fT0   , nwf);

    R__b << fNMcParticles;
    int ns = NStoredMcParticles();
    R__b.WriteFastArray(fMcSimID, ns);
    R__b.WriteFastArray(fMcPdgID, ns);
    R__b.WriteFastArray(fMcNHits, ns);
    R__b.WriteFastArray(fMcShare, ns);
  }
}

//-----------------------------------------------------------------------------
TStnHelix::TStnHelix(int Number) {

  fSimpMom1.SetPxPyPzE(1e6,1e6,1e6,1e6);
  fSimpOrigin1.SetXYZT(0.,0.,0.,0.);
  
  fSimpMom2.SetPxPyPzE(1e6,1e6,1e6,1e6);
  fSimpOrigin2.SetXYZT(0.,0.,0.,0.);
  
  fNumber           = Number;
  fNHits            =  0;
  fAlgorithmID      = -1;
  fTimeClusterIndex = -1;
  fTrackSeedIndex   = -1;
  fNComboHits       = -1;

  fSimpPDG1     = 0;    
  fSimpPDGM1    = 0;    
  fSimpId1Hits  = 0;
  fSimpPDG2     = 0; 
  fSimpPDGM2    = 0; 
  fSimpId2Hits  = 0;
  fHelicity     = 0;
  fPropDir      = 0;
  fSimpID1      =-1;

  for (int i=0; i<kNFreeInts; i++) fInt[i] = 0;

  fT0     = 0;    
  fT0Err  = 0; 
  fRCent  = 0;    
  fFCent  = 0;  
  fRadius = 0; 
  fLambda = 0;    
  fFZ0    = 0;
  fD0     = 0;

  fChi2XYNDof   = 0;
  fChi2PhiZNDof = 0;

  fClusterTime   = 0;  
  fClusterEnergy = 0;
  fClusterX      = 0;     
  fClusterY      = 0;     
  fClusterZ      = 0;     

  fNLoops        = 0;
  fTZSlope       = 0;
  fTZSlopeError  = 0;
  fChi2TZNDof    = 0;
  fHitRatio      = 0; 

  for (int i=0; i<kNFreeFloats; i++) fFloat[i] = 0;

  fNMcParticles  = -1;
}

//_____________________________________________________________________________
TStnHelix::~TStnHelix() {
}

//-----------------------------------------------------------------------------
int TStnHelix::NHitsFromSimID(int SimID) const {
  if (fNMcParticles < 0) return -1;

  int ns = NStoredMcParticles();
  for (int i=0; i<ns; i++) {
    if (fMcSimID[i] == SimID) return fMcNHits[i];
  }
  return 0;
}

//-----------------------------------------------------------------------------
float TStnHelix::McPurity() const {
  int ns = NStoredMcParticles();
  if (ns <= 0) return -1;

  int nsum = 0;
  for (int i=0; i<ns; i++) nsum += fMcNHits[i];
  return (nsum > 0) ? float(fMcNHits[0])/nsum : -1;
}


//-----------------------------------------------------------------------------
void TStnHelix::Clear(Option_t* opt) {
  Error("Print", "Not implemented yet");
}

//-----------------------------------------------------------------------------
void TStnHelix::Print(Option_t* Option) const {
  TString opt(Option);

  opt.ToLower();
					// "non-const *this" for printing purposes
  TStnHelix* h = (TStnHelix*) this;

  if ((opt == "") || (opt.Index("banner") >= 0)) {
//-----------------------------------------------------------------------------
// print banner
//-----------------------------------------------------------------------------
    printf("--------------------------------------------------------------------------------------------------------");
    printf("---------------------------------------\n");
    printf(" i hel  nh nch tci tsi    AlgID      Radius      X0        Y0        D0      Lambda      p         pT   ");
    printf("    T0        T0Err    Chi2XY  Chi2ZPhi\n");
    printf("--------------------------------------------------------------------------------------------------------");
    printf("---------------------------------------\n");
  }

  if ((opt == "") || (opt.Index("data") >= 0)) {
    printf("%2i %3i %3i %3i %3i %3i 0x%08x %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f",
	   h->Number(), h->Helicity(), h->NHits(), h->NComboHits(), h->TimeClusterIndex(), h->TrackSeedIndex(),
	   h->AlgorithmID(),
	   h->Radius(),h->CenterX(),h->CenterY(),h->D0(),h->Lambda(),h->P(),h->Pt(),h->T0(),h->T0Err(),
	   h->Chi2XY(),h->Chi2ZPhi());
    printf("\n");
  }
}

// } // end namespace




