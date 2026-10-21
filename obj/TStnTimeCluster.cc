//
#include <iostream>

#include "TMatrix.h"
#include "obj/TStnTimeCluster.hh"

// #include "RecoDataProducts/inc/TimeCluster.hh"

namespace {
  //  const double BF = 1.4116 ;  // CDF case
  const double BF = 1.0 ;         // MU2E case
}

//namespace murat {

ClassImp(TStnTimeCluster)

//-----------------------------------------------------------------------------
void TStnTimeCluster::ReadV1(TBuffer &R__b) {

  struct TStnTimeClusterDataV1_t {
    int                       fNHits;
    int                       fHelixSeedIndex;

    float                     fT0;    
    float                     fT0Err; 

    float                     fPosX;
    float                     fPosY;
    float                     fPosZ;
 
    float                     fClusterTime;   
    float		      fClusterEnergy; 
    float		      fClusterX;      
    float		      fClusterY;      
    float		      fClusterZ;      
    
    const mu2e::TimeCluster*  fTimeCluster;  //!
  };

  TStnTimeClusterDataV1_t data;

  int            nwi, nwf;
  
  nwi      = ((int*  ) &data.fT0            ) - &data.fNHits;
  nwf      = ((float*) &data.fTimeCluster   ) - &data.fT0;
    
  R__b.ReadFastArray(&data.fNHits ,nwi);
  R__b.ReadFastArray(&data.fT0    ,nwf);

  fNHits             = data.fNHits;	     
  fHelixSeedIndex    = data.fHelixSeedIndex;

  fT0                = data.fT0;    	  
  fT0Err             = data.fT0Err; 	  
                
  fPosX              = data.fPosX; 
  fPosY              = data.fPosY; 
  fPosZ              = data.fPosZ; 

  fClusterTime       = data.fClusterTime; 
  fClusterEnergy     = data.fClusterEnergy; 
  fClusterX          = data.fClusterX;      
  fClusterY          = data.fClusterY;    
  fClusterZ          = data.fClusterZ;    

  fNComboHits        = -1;		// added in V2

  fSimID             = -1;		// added in V3
  fPdgID             = -9999999;	// added in V3
  fNHitsSimID        = -1;		// added in V3
  fMcMom             = -1.;		// added in V3

  fNMcParticles      = -1;		// added in V4
}


//-----------------------------------------------------------------------------
void TStnTimeCluster::ReadV2(TBuffer &R__b) {

  struct TStnTimeClusterDataV2_t {
    int                       fNHits;
    int                       fHelixSeedIndex;
    int                       fNComboHits;        // added in V2: 
    int                       fInt[kNFreeIntsV2]; // added in V2: provision for future I/O expansion

    float                     fT0;    
    float                     fT0Err; 

    float                     fPosX;
    float                     fPosY;
    float                     fPosZ;
 
    float                     fClusterTime;   
    float		      fClusterEnergy; 
    float		      fClusterX;      
    float		      fClusterY;      
    float		      fClusterZ;      
    float                     fFloat[kNFreeFloatsV2]; // added in V2: provision for future I/O expansion
    
    const mu2e::TimeCluster*  fTimeCluster;  //!
  };

  TStnTimeClusterDataV2_t data;

  int            nwi, nwf;
  
  nwi      = ((int*  ) &data.fT0            ) - &data.fNHits;
  nwf      = ((float*) &data.fTimeCluster   ) - &data.fT0;
    
  R__b.ReadFastArray(&data.fNHits ,nwi);
  R__b.ReadFastArray(&data.fT0    ,nwf);

  fNHits             = data.fNHits;	     
  fHelixSeedIndex    = data.fHelixSeedIndex;
  fNComboHits        = data.fNComboHits;

  fT0                = data.fT0;    	  
  fT0Err             = data.fT0Err; 	  
                
  fPosX              = data.fPosX; 
  fPosY              = data.fPosY; 
  fPosZ              = data.fPosZ; 

  fClusterTime       = data.fClusterTime; 
  fClusterEnergy     = data.fClusterEnergy; 
  fClusterX          = data.fClusterX;      
  fClusterY          = data.fClusterY;    
  fClusterZ          = data.fClusterZ;

  fSimID             = -1;		// added in V3
  fPdgID             = -999999;		// added in V3
  fNHitsSimID        = -1;		// added in V3
  fMcMom             = -1.;		// added in V3

  fNMcParticles      = -1;		// added in V4
}

//-----------------------------------------------------------------------------
// V3: same fixed part as V4, no MC particle list
//-----------------------------------------------------------------------------
void TStnTimeCluster::ReadV3(TBuffer &R__b) {

  int nwi = ((int*  ) &fT0            ) - &fNHits;
  int nwf = ((float*) &fTimeCluster   ) - &fT0;

  R__b.ReadFastArray(&fNHits, nwi);
  R__b.ReadFastArray(&fT0   , nwf);

  fNMcParticles      = -1;		// added in V4
}

//-----------------------------------------------------------------------------
void TStnTimeCluster::Streamer(TBuffer& R__b) {
  
  int nwi = ((int*  ) &fT0            ) - &fNHits;
  int nwf = ((float*) &fTimeCluster   ) - &fT0;

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); 
    
    if      (R__v == 1) ReadV1(R__b);
    else if (R__v == 2) ReadV2(R__b);
    else if (R__v == 3) ReadV3(R__b);
    else {
                                        // current version: V4
      R__b.ReadFastArray(&fNHits, nwi);
      R__b.ReadFastArray(&fT0   , nwf);

      R__b >> fNMcParticles;
      int ns = NStoredMcParticles();
      R__b.ReadFastArray(fMcSimID, ns);
      R__b.ReadFastArray(fMcPdgID, ns);
      R__b.ReadFastArray(fMcNHits, ns);
    }
  }
  else {
    R__b.WriteVersion(TStnTimeCluster::IsA());

    R__b.WriteFastArray(&fNHits, nwi);
    R__b.WriteFastArray(&fT0   , nwf);

    R__b << fNMcParticles;
    int ns = NStoredMcParticles();
    R__b.WriteFastArray(fMcSimID, ns);
    R__b.WriteFastArray(fMcPdgID, ns);
    R__b.WriteFastArray(fMcNHits, ns);
  }
}


//-----------------------------------------------------------------------------
TStnTimeCluster::TStnTimeCluster(int i) {

  fNHits         = 0;

  fNComboHits    = 0;
  fT0            = 0;    
  fT0Err         = 0; 

  fPosX          = 0;
  fPosY          = 0;
  fPosZ          = 0;

  fClusterTime   = 0;  
  fClusterEnergy = 0;
  fClusterX      = 0;     
  fClusterY      = 0;     
  fClusterZ      = 0;     

  fSimID         = -1;
  fPdgID         = -999999;
  fNHitsSimID    = -1;
  fMcMom         = -1.;

  fNMcParticles  = -1;
}

//_____________________________________________________________________________
TStnTimeCluster::~TStnTimeCluster() {
}

//-----------------------------------------------------------------------------
int TStnTimeCluster::NHitsFromSimID(int SimID) const {
  if (fNMcParticles < 0) return -1;

  int ns = NStoredMcParticles();
  for (int i=0; i<ns; i++) {
    if (fMcSimID[i] == SimID) return fMcNHits[i];
  }
  return 0;
}


//-----------------------------------------------------------------------------
void TStnTimeCluster::Clear(Option_t* opt) {
  Error("Print", "Not implemented yet");
}

//-----------------------------------------------------------------------------
void TStnTimeCluster::Print(Option_t* Option) const {
  TString opt(Option);

  opt.ToLower();
					// "non-const *this" for printing purposes
  //  TStnTimeCluster* cl = (TStnTimeCluster*) this;

  if ((opt == "") || (opt.Index("banner") >= 0)) {
//-----------------------------------------------------------------------------
// print banner
//-----------------------------------------------------------------------------
    printf("-------------------------------------------------------\n");
    printf("  i  nch nsh   T0  T0Err \n");
    printf("-------------------------------------------------------\n");
  }

  if ((opt == "") || (opt.Index("data") >= 0)) {
    printf("%4i %3i %3i %8.3f %8.3f ",
	   -1,fNComboHits,fNHits,fT0, fT0Err);
    printf("\n");
  }
}

// } // end namespace




