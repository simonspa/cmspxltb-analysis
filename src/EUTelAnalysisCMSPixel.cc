// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#if defined USE_GEAR

// EUTelescope includes:
#include "EUTelAnalysisCMSPixel.h"
#include "EUTELESCOPE.h"
#include "EUTelSparseDataImpl.h"
#include "EUTelBaseSparsePixel.h"
#include "EUTelSimpleSparsePixel.h"
#include "EUTelExceptions.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelPStream.h" // process streams redi::ipstream

// AIDA histogram package (on top of ROOT):

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
#include <AIDA/IProfile2D.h>
#include <AIDA/ITree.h>
#endif

// GBL:
#include "include/GblTrajectory.h"
#include "include/MilleBinary.h"


// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/ProcessorMgr.h"
#include "marlin/Global.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// LCIO includes <.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>

// system includes <>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sstream>
#include <memory>
#include <string.h>
#include <map>
#include <cstdlib>
#include <limits>

// ROOT includes ".h"
#include <TMath.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TVector3.h>
#include <TRotation.h>
#include "TH1D.h"

using namespace std;
//FIXME stupid:
#include "CMSpixelClust.h"
using namespace lcio;
using namespace marlin;
using namespace eutelescope;
using namespace gbl;

// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

AIDA::IHistogram1D *t1Histo, *t10Histo, *t100Histo, *t300Histo, *t600Histo, *t1000Histo, *t1800Histo, *t3600Histo;
AIDA::IHistogram1D *dtHisto, *dtmsHisto, *logdtHisto, *logdtcmsHisto;
AIDA::IProfile1D *dtfvstau, *tfvstau;
AIDA::IHistogram1D * dtfHisto;
AIDA::IProfile1D * dtfvst, *dtfvsdt;
AIDA::IHistogram1D * tfHisto;
AIDA::IProfile1D * tfvst, * tfvst1, * tfvst10, * tfvst100, * tfvst300;

AIDA::IHistogram1D * nAllHitHisto, * hits0Histo, * hits1Histo, * hits2Histo, * hits3Histo, * hits4Histo, * hits5Histo;

AIDA::IHistogram1D * dutnclusHisto, * dutcolHisto, * dutrowHisto, * dutnpxHisto, * dutadcHisto;
AIDA::IHistogram1D * refnclusHisto, * refcolHisto, * refrowHisto, * refnpxHisto, * refadcHisto;

AIDA::IHistogram1D * cmsdtHisto;
AIDA::IHistogram1D * dutrefddtHisto;
AIDA::IHistogram1D * sysrtHisto;
AIDA::IHistogram1D * sysrdtHisto;
AIDA::IHistogram1D * dutddtnsHisto;
AIDA::IHistogram1D * refddtnsHisto;
AIDA::IHistogram1D * dutddtusHisto;
AIDA::IHistogram1D * dutddtmsHisto;
AIDA::IProfile1D * dutddtvst;
AIDA::IProfile1D * dutddtvsdt;
AIDA::IProfile1D * ddtvst;
AIDA::IProfile1D * ddtvstms;
AIDA::IProfile1D * ddtvsdt;
AIDA::IHistogram1D * gapdtHisto;


// Correlation plots for telescope planes
AIDA::IHistogram1D * dx01Histo, * dy01Histo, * du01Histo, * dx02Histo, * dx03Histo, * dx04Histo, * dx05Histo, * dx12Histo, * dy12Histo, * du12Histo, * dx23Histo, * dy23Histo, * du23Histo, * dx34Histo, * dy34Histo, * du34Histo, * dx45Histo, * dy45Histo, * du45Histo;


// triplets 0-1-2:

AIDA::IHistogram1D * da02Histo;
AIDA::IHistogram1D * db02Histo;

AIDA::IProfile2D * dzcvsxy;
AIDA::IProfile2D * z3vsxy;

AIDA::IHistogram1D * tridxHisto;
AIDA::IHistogram1D * tridyHisto;

AIDA::IProfile1D * tridxvsx;
AIDA::IProfile1D * tridxvsy;
AIDA::IProfile1D * tridxvstx;
AIDA::IProfile1D * tridxvsty;
AIDA::IProfile1D * tridyvsx;
AIDA::IProfile1D * tridyvsy;
AIDA::IProfile1D * tridyvstx;
AIDA::IProfile1D * tridyvsty;

AIDA::IHistogram1D * tridx1Histo, * tridy1Histo, * tridx3Histo, * tridy3Histo, * tridx3bHisto, * tridy3bHisto, * tridx4Histo, * tridy4Histo, * tridx4bHisto, * tridy4bHisto, * tridx5Histo, * tridy5Histo, * tridx5bHisto, * tridy5bHisto, * trixHisto, * triyHisto, * tritxHisto, * trityHisto;
AIDA::IHistogram2D * trixyHisto;

AIDA::IHistogram1D * trixdutHisto; // at DUT
AIDA::IHistogram1D * triydutHisto;
AIDA::IHistogram2D * trixydutHisto;

AIDA::IHistogram2D * cmsxxHisto;
AIDA::IHistogram2D * cmsyyHisto;

AIDA::IHistogram1D * cmspxqHisto, * cmssxaHisto, * cmsdyaHisto, * cmsdxaHisto, * cmssyaHisto, * cmsdx4Histo, * cmsdy4Histo, * cmsdxHisto, * cmsdyHisto, * cmsdxfHisto, * cmsdyfHisto, * cmsdxfcHisto, * cmsdyfcHisto, * cmsdyfc1Histo, * cmsdyfc2Histo, * cmsdyfc3Histo, * cmsdyq0Histo, * cmsdyq1Histo, * cmsdyq2Histo, * cmsdyeta0Histo, * cmsdyeta1Histo, * cmsdxfctHisto, * cmsdyfctHisto, * cmsdyfcntHisto, * cmsdxfctqHisto, * cmsdyfctqHisto, * cmsdyfcntqHisto, * cmsdxfctq1Histo, * cmsdyfctq1Histo, * cmsdyfcntq1Histo, * cmsdyfctq1lHisto, * cmsdyfctq1rHisto, * cmsdxfctq2Histo, * cmsdyfctq2Histo, * cmsdxfctq3Histo, * cmsdyfctq3Histo, * cmsdyfctqdotHisto, * cmsdyfctq3dHisto, * cmscolHisto, * cmsrowHisto, * cmsqHisto, * cmsq0Histo, * trixlkHisto, * triylkHisto;

AIDA::IHistogram2D * trixylkHisto;
AIDA::IProfile1D * cmsdxvsx;
AIDA::IProfile1D * cmsdyvsx;
AIDA::IProfile1D * cmsdxvsy;
AIDA::IProfile1D * cmsdyvsy;
AIDA::IProfile1D * cmsdxvstx;
AIDA::IProfile1D * cmsdyvsty;
AIDA::IHistogram2D * cmsdyvsxHisto;

AIDA::IHistogram1D * cmsnpxHisto;
AIDA::IHistogram1D * cmsnpx0Histo;
AIDA::IHistogram1D * cmsnpx1Histo;
AIDA::IHistogram1D * cmsnpx2Histo;
AIDA::IHistogram1D * cmsnpx3Histo;
AIDA::IHistogram1D * cmsncolHisto;
AIDA::IHistogram1D * cmsnrowHisto;
AIDA::IHistogram1D * cmsnrowqHisto;
AIDA::IProfile1D * cmsnrowvst1;
AIDA::IHistogram1D * cmsetaHisto;
AIDA::IHistogram1D * cmsqfHisto;
AIDA::IHistogram1D * cmsq0fHisto;
AIDA::IHistogram1D * cmsqf0Histo;
AIDA::IHistogram1D * cmsqf1Histo;
AIDA::IHistogram1D * cmsqf2Histo;
AIDA::IHistogram1D * cmsqf3Histo;

AIDA::IProfile1D * cmsdyvsxm;
AIDA::IProfile1D * cmsdyvsym;
AIDA::IHistogram2D * cmspixvsxmym;

AIDA::IProfile1D * cmsqvsx;
AIDA::IProfile1D * cmsqvsy;
AIDA::IProfile1D * cmsqvsxm;
AIDA::IProfile1D * cmsqvsym;
AIDA::IProfile2D * cmsqvsxmym;
AIDA::IProfile1D * cmsqvsddt;
AIDA::IProfile1D * cmsqvst1;
AIDA::IProfile1D * cmsqvst2;
AIDA::IProfile1D * cmsqvst3;
AIDA::IProfile1D * cmsqvst4;

AIDA::IProfile1D * cmsrmsxvsx;
AIDA::IProfile1D * cmsrmsyvsx;
AIDA::IProfile1D * cmsrmsxvsy;
AIDA::IProfile1D * cmsrmsyvsy;
AIDA::IProfile1D * cmsrmsxvsxm;
AIDA::IProfile1D * cmsrmsyvsxm;
AIDA::IProfile1D * cmsncolvsxm;
AIDA::IProfile1D * cmsnrowvsxm;
AIDA::IProfile1D * cmsrmsxvsym;
AIDA::IProfile1D * cmsrmsyvsym;
AIDA::IProfile1D * cmsrmsyvsym3;
AIDA::IProfile1D * cmsrmsyvsym6;
AIDA::IProfile1D * cmsrmsyvst;
AIDA::IProfile1D * cmsrmsyvsddt;
AIDA::IProfile1D * cmsrmsxvsq;
AIDA::IProfile1D * cmsrmsyvsq;

AIDA::IProfile1D * cmsdyvseta;
AIDA::IProfile1D * cmsrmsyvseta;

AIDA::IProfile1D * cmspMoyalvsq;
AIDA::IHistogram1D * cmspMoyalHisto;
AIDA::IProfile1D * cmsrmsyvsp;

AIDA::IProfile2D * cmsnpxvsxmym;
AIDA::IProfile1D * cmsncolvsym;
AIDA::IProfile1D * cmsnrowvsym;
AIDA::IProfile1D * cmsetavsym;
AIDA::IProfile1D * cmsetavsym3;
AIDA::IProfile1D * cmsetavsym2;
AIDA::IHistogram1D * cmsym1Histo;
AIDA::IHistogram1D * cmsym2Histo;

AIDA::IHistogram2D * effxyHisto;
AIDA::IProfile2D * effvsxy;
AIDA::IProfile1D * effvsx;
AIDA::IProfile1D * effvsxg;
AIDA::IProfile1D * effvsy;
AIDA::IProfile1D * eff300;
AIDA::IProfile1D * eff600;
AIDA::IProfile1D * eff1200;
AIDA::IProfile1D * eff1800;
AIDA::IProfile1D * eff3600;
AIDA::IProfile1D * effvsddt;
AIDA::IProfile2D * effvsxmym;
AIDA::IProfile1D * effd600;
AIDA::IProfile1D * effn600;
AIDA::IProfile1D * effm600;

AIDA::IProfile2D * rffvsxy;
AIDA::IProfile1D * rffvsx;

AIDA::IHistogram1D * ntriHisto;
AIDA::IProfile1D * lkAvst;

// triplet eff w.r.t. CMS:

AIDA::IHistogram1D * cmsxeHisto;
AIDA::IHistogram1D * cmsyeHisto;

AIDA::IHistogram1D * cmsdxeHisto;
AIDA::IHistogram1D * cmsdyeHisto;

AIDA::IHistogram1D * cmsnmHisto;
AIDA::IProfile2D * trieffvsxy;

// driplets 3-4-5:

AIDA::IHistogram1D * dx35Histo;
AIDA::IHistogram1D * dy35Histo;

AIDA::IHistogram1D * dridxHisto;
AIDA::IHistogram1D * dridyHisto;
AIDA::IHistogram1D * drixHisto;
AIDA::IHistogram1D * driyHisto;
AIDA::IHistogram2D * drixyHisto;
AIDA::IHistogram1D * dritxHisto;
AIDA::IHistogram1D * drityHisto;

AIDA::IProfile1D * dridxvsx;
AIDA::IProfile1D * dridxvsy;
AIDA::IProfile1D * dridxvstx;
AIDA::IProfile1D * dridxvsty;
AIDA::IProfile1D * dridyvsx;
AIDA::IProfile1D * dridyvsy;
AIDA::IProfile1D * dridyvstx;
AIDA::IProfile1D * dridyvsty;

AIDA::IHistogram1D * drixrefHisto; // at REF
AIDA::IHistogram1D * driyrefHisto;
AIDA::IHistogram2D * drixyrefHisto;

AIDA::IHistogram1D * drixlkHisto;
AIDA::IHistogram1D * driylkHisto;
AIDA::IHistogram2D * drixylkHisto;
AIDA::IHistogram2D * refpixvsxmym;

AIDA::IHistogram1D * refqHisto;
AIDA::IProfile2D * refqvsxmym;

AIDA::IHistogram2D * refxxHisto; //REF vs driplet
AIDA::IHistogram2D * refyyHisto;

AIDA::IHistogram1D * refsxaHisto;
AIDA::IHistogram1D * refdyaHisto;
AIDA::IHistogram1D * refsxHisto;
AIDA::IHistogram1D * refdyHisto;
AIDA::IHistogram1D * refsxcHisto;
AIDA::IHistogram1D * refdycHisto;

AIDA::IProfile1D * refdyvsx;
AIDA::IProfile1D * refdyvsy;
AIDA::IProfile1D * refdyvsty;

AIDA::IHistogram1D * reflkcolHisto;
AIDA::IHistogram1D * reflkrowHisto;

AIDA::IHistogram1D * bacsxaHisto;
AIDA::IHistogram1D * bacdyaHisto;
AIDA::IHistogram1D * bacsxcHisto;
AIDA::IHistogram1D * bacdycHisto;
AIDA::IHistogram1D * bacsxcqHisto;
AIDA::IHistogram1D * bacdycqHisto;

AIDA::IHistogram1D * ndriHisto;
AIDA::IHistogram1D * ndrirefHisto;
AIDA::IProfile1D * lkBvst;

AIDA::IHistogram1D * nsixHisto;
AIDA::IHistogram1D * sixkxHisto; //driplet-triplet
AIDA::IHistogram1D * sixkyHisto;
AIDA::IHistogram1D * sixdxHisto;
AIDA::IHistogram1D * sixdyHisto;
AIDA::IHistogram1D * sixdxcHisto;
AIDA::IHistogram1D * sixdycHisto;

AIDA::IHistogram1D * sixkxcHisto;
AIDA::IHistogram1D * sixkycHisto;
AIDA::IHistogram1D * sixxHisto;
AIDA::IHistogram1D * sixyHisto;
AIDA::IHistogram2D * sixxyHisto;
AIDA::IHistogram2D * sixxycHisto;
AIDA::IProfile2D * kinkvsxy;

AIDA::IHistogram1D * sixx0Histo;
AIDA::IHistogram1D * sixy0Histo;
AIDA::IHistogram1D * sixx1Histo;
AIDA::IHistogram1D * sixy1Histo;
AIDA::IHistogram1D * sixx2Histo;
AIDA::IHistogram1D * sixy2Histo;
AIDA::IHistogram1D * sixx3Histo;
AIDA::IHistogram1D * sixy3Histo;
AIDA::IHistogram1D * sixx4Histo;
AIDA::IHistogram1D * sixy4Histo;
AIDA::IHistogram1D * sixx5Histo;
AIDA::IHistogram1D * sixy5Histo;

AIDA::IHistogram2D * sixxylkHisto;

AIDA::IHistogram1D * derxtiltHisto;
AIDA::IHistogram1D * derytiltHisto;
AIDA::IHistogram1D * derxturnHisto;
AIDA::IHistogram1D * deryturnHisto;

AIDA::IHistogram1D * selxHisto;
AIDA::IHistogram1D * selyHisto;
AIDA::IHistogram1D * selaxHisto;
AIDA::IHistogram1D * selayHisto;
AIDA::IHistogram1D * seldxHisto;
AIDA::IHistogram1D * seldyHisto;
AIDA::IHistogram1D * selkxHisto;
AIDA::IHistogram1D * selkyHisto;

AIDA::IHistogram1D * seldx1Histo;
AIDA::IHistogram1D * seldy1Histo;
AIDA::IHistogram1D * seldx3Histo;
AIDA::IHistogram1D * seldy3Histo;
AIDA::IHistogram1D * seldx4Histo;
AIDA::IHistogram1D * seldy4Histo;
AIDA::IHistogram1D * seldx5Histo;
AIDA::IHistogram1D * seldy5Histo;
AIDA::IHistogram1D * seldx6Histo;
AIDA::IHistogram1D * seldy6Histo;

AIDA::IHistogram1D * gblndfHisto;
AIDA::IHistogram1D * gblchi2aHisto;
AIDA::IHistogram1D * gblchi2bHisto;
AIDA::IHistogram1D * gblprbHisto;

AIDA::IHistogram1D * badxHisto;
AIDA::IHistogram1D * badyHisto;
AIDA::IHistogram1D * badaxHisto;
AIDA::IHistogram1D * badayHisto;
AIDA::IHistogram1D * baddxHisto;
AIDA::IHistogram1D * baddyHisto;
AIDA::IHistogram1D * badkxHisto;
AIDA::IHistogram1D * badkyHisto;

AIDA::IHistogram1D * baddx1Histo;
AIDA::IHistogram1D * baddy1Histo;
AIDA::IHistogram1D * baddx3Histo;
AIDA::IHistogram1D * baddy3Histo;
AIDA::IHistogram1D * baddx4Histo;
AIDA::IHistogram1D * baddy4Histo;
AIDA::IHistogram1D * baddx5Histo;
AIDA::IHistogram1D * baddy5Histo;
AIDA::IHistogram1D * baddx6Histo;
AIDA::IHistogram1D * baddy6Histo;

AIDA::IHistogram1D * goodx1Histo;
AIDA::IHistogram1D * goody1Histo;
AIDA::IHistogram1D * goodx6Histo;
AIDA::IHistogram1D * goody6Histo;

AIDA::IHistogram1D * gblax0Histo;
AIDA::IHistogram1D * gbldx0Histo;
AIDA::IHistogram1D * gblrx0Histo;
AIDA::IHistogram1D * gblpx0Histo;
AIDA::IHistogram1D * gblqx0Histo;

AIDA::IHistogram1D * gblax1Histo;
AIDA::IHistogram1D * gbldx1Histo;
AIDA::IHistogram1D * gblrx1Histo;
AIDA::IHistogram1D * gblpx1Histo;
AIDA::IHistogram1D * gblqx1Histo;
AIDA::IHistogram1D * gblsx1Histo;
AIDA::IHistogram1D * gbltx1Histo;

AIDA::IHistogram1D * gblax2Histo;
AIDA::IHistogram1D * gbldx2Histo;
AIDA::IHistogram1D * gblrx2Histo;
AIDA::IHistogram1D * gblpx2Histo;
AIDA::IHistogram1D * gblqx2Histo;

AIDA::IHistogram1D * gblax3Histo;
AIDA::IHistogram1D * gbldx3Histo;
AIDA::IHistogram1D * gblrx3Histo;
AIDA::IHistogram1D * gblpx3Histo;
AIDA::IHistogram1D * gblqx3Histo;

AIDA::IHistogram1D * gblax4Histo;
AIDA::IHistogram1D * gbldx4Histo;
AIDA::IHistogram1D * gblrx4Histo;
AIDA::IHistogram1D * gblpx4Histo;
AIDA::IHistogram1D * gblqx4Histo;

AIDA::IHistogram1D * gblax5Histo;
AIDA::IHistogram1D * gbldx5Histo;
AIDA::IHistogram1D * gblrx5Histo;
AIDA::IHistogram1D * gblpx5Histo;
AIDA::IHistogram1D * gblqx5Histo;

AIDA::IHistogram1D * gblax6Histo;
AIDA::IHistogram1D * gbldx6Histo;
AIDA::IHistogram1D * gbldy6Histo;
AIDA::IHistogram1D * gblrx6Histo;
AIDA::IHistogram1D * gblry6Histo;
AIDA::IHistogram1D * gblpx6Histo;
AIDA::IHistogram1D * gblpy6Histo;
AIDA::IHistogram1D * gblqx6Histo;
AIDA::IHistogram1D * gblsx6Histo;
AIDA::IHistogram1D * gbltx6Histo;

AIDA::IHistogram1D * gblkx1Histo;
AIDA::IHistogram1D * gblkx2Histo;
AIDA::IHistogram1D * gblkx3Histo;
AIDA::IHistogram1D * gblkx4Histo;
AIDA::IHistogram1D * gblkx5Histo;
AIDA::IHistogram1D * gblkx6Histo;

AIDA::IHistogram1D * sixzx3Histo;
AIDA::IHistogram1D * sixzy3Histo;
AIDA::IHistogram1D * sixzx2Histo;
AIDA::IHistogram1D * sixzy2Histo;
AIDA::IHistogram1D * sixzx1Histo;
AIDA::IHistogram1D * sixzy1Histo;

AIDA::IHistogram1D * sixkxzyHisto;
AIDA::IHistogram1D * sixkyzxHisto;
AIDA::IHistogram1D * sixkxzxHisto;
AIDA::IHistogram1D * sixkyzyHisto;

#endif

int dutnsync = 0;
int refnsync = 0;
int nclst = 0;
int nmtch = 0;
int naldut = 0;
int ngbl = 0;
bool ldb = 0;

MilleBinary * mille; // for producing MillePede-II binary file

double turn = 0;
double tilt = 0;
double DUTz = 40;
double DUTalignx = 0;
double DUTaligny = 0;
double DUTrot = 0;


EUTelAnalysisCMSPixel::EUTelAnalysisCMSPixel() : Processor("EUTelAnalysisCMSPixel"), _isFirstEvent(0), _nSkip(0),_nSkipTelescope(0), _nSkipRef(0), _gTLU(0),_nEvtBeg(0), _nRefBeg(0), _DUT_data(""), _REF_data(""), _DUT_gain(""), _REF_gain(""), _DUT_run(0), _REF_run(0), _DUT_weibull(0), _REF_weibull(0), _DUT_chip(0), _REF_chip(0), _DUTalignx(0), _DUTaligny(0), _DUTz(0), _DUTrot(0), _DUTtilt(0), _DUTturn(0), _REFalignx(0), _REFaligny(0), _REFz(0), _REFrot(0), _eBeam(0), _nRun(0), _nEvt(0), _leff_val(0), _nTelPlanes(0) {

  // modify processor description
  _description = "Analysis for CMS PSI46 Pixel Detectors as DUT in AIDA telescopes ";


  // processor parameters
  registerInputCollection( LCIO::TRACKERHIT,
                           "InputCollectionTelescope" ,
                           "Name of the input TrackerHit collection of the telescope",
                           _inputCollectionTelescope,
                           std::string("alignedhits") );
  registerInputCollection( LCIO::TRACKERHIT,
                           "InputCollectionDUT" ,
                           "Name of the input TrackerHit collection of the DUT",
                           _inputCollectionDUT,
                           std::string("duthits") );
  registerOutputCollection(LCIO::TRACK,"InputTrackCollectionTelescope",
                           "Name of the input Track collection of the telescope",
                           _inputTrackCollection, string( "fittracks"));

  registerProcessorParameter( "Ebeam",
                              "Beam energy [GeV]",
                             _eBeam, static_cast < double >( 4.0));

  registerOptionalParameter( "nSkipEventDUT",
                              "Skip n CMS DUT events at begin of file",
                              _nSkip, static_cast < int >(0) );
  registerOptionalParameter( "nSkipEventREF",
                              "Skip n CMS REF events at begin of file",
                              _nSkipRef, static_cast < int >(0) );
  registerOptionalParameter( "nSkipEventTel",
                              "Skip n Telescope events at begin of file",
                              _nSkipTelescope, static_cast < int >(0) );

  registerProcessorParameter( "gTLU",
                              "DESY-II Umlauftakt",
                             _gTLU, static_cast < double >(0.384678));

  registerProcessorParameter("CMS_data_path",
                             "Path to the native CMS data files",
                             _CMS_data_path, string("/dev/null"));
  registerProcessorParameter("CMS_gain_path",
                             "Path to the CMS gain calibration files",
                             _CMS_gain_path, string("/dev/null"));


  registerProcessorParameter("DUT_run",
                             "Runnumber of the DUT associated with this telescope run",
                             _DUT_run, static_cast< int > (0) );
  registerProcessorParameter("REF_run",
                             "Runnumber of the REF associated with this telescope run",
                             _REF_run, static_cast< int > (0) );

  registerProcessorParameter("DUT_chip",
                             "Serial number of the DUT used in this run",
                             _DUT_chip, static_cast< int > (0) );
  registerProcessorParameter("DUT_gain",
                             "CMS DUT gain file to be used",
                             _DUT_gain, string("/dev/null"));
  registerProcessorParameter("REF_chip",
                             "Serial number of the REF used in this run",
                             _REF_chip, static_cast< int > (0) );
  registerProcessorParameter("DUT_weibull",
                             "Activate DUT Weibull calibration",
                             _DUT_weibull, static_cast< bool >(true));
  registerProcessorParameter("REF_gain",
                             "CMS REF gain file to be used",
                             _REF_gain, string("/dev/null"));
  registerProcessorParameter("REF_weibull",
                             "Activate REF Weibull calibration",
                             _REF_weibull, static_cast< bool >(true));

  // Alignment stuff:

  registerProcessorParameter( "DUT_align_x",
                              "DUT alignment in x",
                             _DUTalignx, static_cast < double >(0));
  registerProcessorParameter( "DUT_align_y",
                              "DUT alignment in y",
                             _DUTaligny, static_cast < double >(0));
  registerProcessorParameter( "DUT_rot",
                              "DUT rotation in xy",
                             _DUTrot, static_cast < double >(0));
  registerProcessorParameter( "DUT_tilt",
                              "DUT tilt angle",
                             _DUTtilt, static_cast < double >(0));
  registerProcessorParameter( "DUT_turn",
                              "DUT turn angle",
                             _DUTturn, static_cast < double >(0));
  registerProcessorParameter( "DUT_pos_z",
                              "DUT position z in mm after the telescope plane",
                             _DUTz, static_cast < double >(0));

  registerProcessorParameter( "REF_align_x",
                              "REF alignment in x",
                             _REFalignx, static_cast < double >(0));
  registerProcessorParameter( "REF_align_y",
                              "REF alignment in y",
                             _REFaligny, static_cast < double >(0));
  registerProcessorParameter( "REF_rot",
                              "REF rotation in xy",
                             _REFrot, static_cast < double >(0));
  registerProcessorParameter( "REF_pos_z",
                              "REF position z in mm after the last telescope plane",
                             _REFz, static_cast < double >(0));

  // Stuff only needed for the printout of the updated runlist line:
  registerOptionalParameter( "DUT_board",
                              "PSI testboard ID of the DUT chip",
			     _DUT_board, string("0") );
  registerOptionalParameter( "REF_board",
                              "PSI testboard ID of the REF chip",
			     _REF_board, string("0") );
  registerOptionalParameter( "date",
                              "Date when the test beam was taken",
			     _DATE_run, string("1970-01-01") );
  registerOptionalParameter( "gearfile",
                              "Again, the gear fiel of the used telescope configuration",
			     _gearfile, string("none") );
  
  
}


void EUTelAnalysisCMSPixel::init() {

  // usually a good idea to
  printParameters();

  _nRun = 0;
  _nEvt = 0;

  _isFirstEvent = true;

  streamlog_out(MESSAGE0) << "Beam energy " << _eBeam << " GeV" <<  endl;

  // Read geometry information from GEAR

  streamlog_out( MESSAGE0 ) << "Reading telescope geometry description from GEAR " << endl;

  _siPlanesParameters  = const_cast<gear::SiPlanesParameters* >( &(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<gear::SiPlanesLayerLayout*>(  &(_siPlanesParameters->getSiPlanesLayerLayout() ));


  // Take all layers defined in GEAR geometry
  _nTelPlanes = _siPlanesLayerLayout->getNLayers();

  _planeSort = new int[_nTelPlanes];
  _planePosition = new double[_nTelPlanes]; // z pos
  _planeID         = new int[_nTelPlanes];
  _planeThickness  = new double[_nTelPlanes];
  _planeX0         = new double[_nTelPlanes];
  _planeResolution = new double[_nTelPlanes];

  for(int i = 0; i < _nTelPlanes; i++) {
    
      _planeID[i]=_siPlanesLayerLayout->getID(i);
      _planePosition[i]=_siPlanesLayerLayout->getLayerPositionZ(i);
      _planeThickness[i]=_siPlanesLayerLayout->getLayerThickness(i);
      _planeX0[i]=_siPlanesLayerLayout->getLayerRadLength(i);
      _planeResolution[i] = _siPlanesLayerLayout->getSensitiveResolution(i);
  }

  streamlog_out( MESSAGE2 ) <<  "Telescope configuration with " << _nTelPlanes << " planes" << endl;

  for( int ipl = 0; ipl < _nTelPlanes; ipl++) {
    stringstream ss;
      ss << "  ID = " << _planeID[ipl]
         << "  at Z [mm] = " << _planePosition[ipl]
         << " dZ [um] = " << _planeThickness[ipl]*1000.;

        ss << "  Res [um] = " << _planeResolution[ipl]*1000.;

	streamlog_out( MESSAGE2 ) <<  ss.str() << endl;
    
  }

  // Book histograms:

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  bookHistos();
#endif

  mille = new gbl::MilleBinary( "mille.bin" );

}//init

//------------------------------------------------------------------------------
void EUTelAnalysisCMSPixel::processRunHeader( LCRunHeader* runHeader) {

  auto_ptr<EUTelRunHeaderImpl> eutelHeader( new EUTelRunHeaderImpl( runHeader ) );
  eutelHeader->addProcessor( type() );

  _nRun++;

  // Decode and print out Run Header information - just a check

  streamlog_out( MESSAGE2 )  << "Processing run header " << _nRun
                             << ", run nr " << runHeader->getRunNumber() << endl;
  _TEL_run = runHeader->getRunNumber();

  const std::string detectorName = runHeader->getDetectorName();
  const std::string detectorDescription = runHeader->getDescription();

  streamlog_out( MESSAGE0 ) << detectorName << " : " << detectorDescription << endl;

} // processRunHeader


TMatrixD Jac5( double ds ) {
  /* for GBL:
     Jacobian for straight line track
     track = q/p, x', y', x, y
             0,   1,  2,  3, 4
  */
  TMatrixD jac(5, 5);
  jac.UnitMatrix();
  jac[3][1] = ds; // x = x0 + xp * ds
  jac[4][2] = ds; // y = y0 + yp * ds
  return jac;
}

void EUTelAnalysisCMSPixel::processEvent( LCEvent * event ) {

  EUTelEventImpl * euEvent = static_cast<EUTelEventImpl*>( event );

  if( euEvent->getEventType() == kEORE ) {
    streamlog_out( DEBUG ) <<  "EORE found: nothing else to do." << endl;
    return;
  }

  bool debug = false;

  // ################## TIMING CALCULATIONS ###########################
  static long prevtime = event->getTimeStamp(); //init
  static long time0 = 0;
  if( _isFirstEvent ) time0 = event->getTimeStamp();
  static long time1 = time0;
  if( _nEvt == 1 ) time1 = event->getTimeStamp();

  // phase from TLU time between events:

  // TLU time stamp: count 384 MHz clock (384 = 48*8)
  // CMS timestamp: count 40 MHz clock, take as reference (HP pulser)
  // dutddt: TLU/CMS = 1.00007
  // => TLU/1.00007 = CMS
  // dutddtns: <TLU/1.00007344 - CMS> = 0

  double gTLU = _gTLU; // [GHz]

  // from dutddtvsdt: gtlu_new = gtlu_old * ( 1 + slope )

  //DP double fTLU = gTLU * 1E9; // [Hz] w.r.t. 40 MHz
  double fTLU = 375.056 * 1.024 * 1E6; // [Hz] = 384.05734 MHz, from Umlauftakt, see below


  // Informational message about event number and timing:
  if( _nEvt < 10 || _nEvt % 100 == 0 ) {

    streamlog_out( MESSAGE2 ) << "Processing event "
                              << setw( 7) << setiosflags(ios::right)
			      << event->getEventNumber()
			      << " in run " << setw( 4) << setiosflags(ios::right)
			      << event->getRunNumber()
			      << ", time " << " = " << fixed << setprecision(1)
			      << ( event->getTimeStamp() - time0 )/fTLU << " s"
			      << ", dt = " << int( (event->getTimeStamp() - prevtime)/gTLU/1E3 ) << " us"
			      << ", rate "
			      << (_nEvt+1) / ( event->getTimeStamp() - time0 + 0.1 ) * fTLU
			      << " Hz"
                              << endl;
  }

  //Filling event timing histograms for the telescope:
  t1Histo->fill( ( event->getTimeStamp() - time0 ) / fTLU ); //event time
  t10Histo->fill( ( event->getTimeStamp() - time0 ) / fTLU ); //event time
  t100Histo->fill( ( event->getTimeStamp() - time0 ) / fTLU ); //event time
  t300Histo->fill( ( event->getTimeStamp() - time0 ) / fTLU ); //event time
  t600Histo->fill( ( event->getTimeStamp() - time0 ) / fTLU ); //event time
  t1000Histo->fill( ( event->getTimeStamp() - time0 ) / fTLU ); //event time
  t1800Histo->fill( ( event->getTimeStamp() - time0 ) / fTLU ); //event time
  t3600Histo->fill( ( event->getTimeStamp() - time0 ) / fTLU ); //event time


  if( _nEvt > 0 ) {

    dtHisto->fill( ( event->getTimeStamp() - prevtime ) / 384E0 ); //us
    dtmsHisto->fill( ( event->getTimeStamp() - prevtime ) / 384E3 ); //ms
    logdtHisto->fill( std::log( ( event->getTimeStamp() - prevtime ) / 384E3 ) / std::log(10.0) );

    // scan for Umlauftakt:

    for( double tau = 372; tau < 378; tau += 0.01 ) {
      double ganz;
      double frac = std::modf( (event->getTimeStamp() - prevtime) / tau, &ganz );
      dtfvstau->fill( tau, abs( frac - 0.5 ) ); // peak at 375.055 for run 5513
      frac = std::modf( (event->getTimeStamp() - time1) / tau, &ganz ); // absolute phase
      tfvstau->fill( tau, abs( frac - 0.5 ) );
    }

    double tau = 375.06106; // Umlauftakt [TLU ticks] for run 5359
    if( event->getRunNumber() >= 5476 ) tau = 375.06109;
    if( event->getRunNumber() >= 5509 ) tau = 375.06104;

    double ganz;
    double frac = std::modf( (event->getTimeStamp() - prevtime) / tau, &ganz );
    if( frac > 0.5 ) frac = frac - 1; // complement for peak around zero
    dtfHisto->fill( frac ); // phase [-0.5,0.5]
    dtfvsdt->fill( ( event->getTimeStamp() - prevtime ) / gTLU/1E3, frac ); // new = old*(1+slope)
    dtfvst->fill( ( event->getTimeStamp() - time0 ) / fTLU, frac ); // phase [-0.5,0.5]

    // absolute phase:

    frac = std::modf( (event->getTimeStamp() - time1) / tau, &ganz ); // ganz = number of turns
    if( frac > 0.5 ) frac = frac - 1; // complement
    tfHisto->fill( frac ); // phase [-0.5,0.5]
    tfvst->fill( ganz, frac ); // phase [-0.5,0.5]
    tfvst1->fill( ganz, frac ); // phase [-0.5,0.5]
    tfvst10->fill( ganz, frac ); // phase [-0.5,0.5]
    tfvst100->fill( ganz, frac ); // phase [-0.5,0.5]
    tfvst300->fill( ganz, frac ); // phase [-0.5,0.5]

    // Detection of long gaps between two telescope events:
    if( ( event->getTimeStamp() - prevtime ) / fTLU > 10 ) {
      cout << "long gap event "
	   << setw( 6) << setiosflags(ios::right)
	   << event->getEventNumber()
	   << " in run " << setw( 6) << setiosflags(ios::right)
	   << event->getRunNumber()
	   << ", time " << setw(10) << setprecision(9)
	   << ( event->getTimeStamp() - time0 ) /fTLU << " s"
	   << ", dt = " << int( (event->getTimeStamp() - prevtime)/384E0 ) << " us"
	   << endl;
    }// Gap detection

  }//dt plots


  //----------------------------------------------------------------------------
  // CMS pixel:

  // Open the CMS pixel data stream:
  if( _isFirstEvent ) {
    std::string dpath = _CMS_data_path + "/" + ZeroPadNumber(_DUT_run,6) + "/mtb.bin";
    std::string rpath = _CMS_data_path + "/" + ZeroPadNumber(_REF_run,6) + "/mtb.bin";
    std::string dgain = _CMS_gain_path + "/chip" + ZeroPadNumber(_DUT_chip,3) + "/" + _DUT_gain;
    std::string rgain = _CMS_gain_path + "/chip" + ZeroPadNumber(_REF_chip,3) + "/" + _REF_gain;
    _dutReader.open(dpath, _DUT_chip, dgain, _DUT_weibull );
    _refReader.open(rpath, _REF_chip, rgain, _REF_weibull );
  }//first event

  vector <cluster> ClustDUT;
  vector <cluster> ClustREF;

  int NpixDUT = 0;
  int NpixREF = 0;

  static long dutTimeStamp = 0;
  static long refTimeStamp = 0;

  static int nmiss = 0; // for DUT
  static int kmiss = 0; // for REF

  // Fixing starting events, sometimes skip one:
  int nEvtBeg = _nSkipTelescope;
  int nRefBeg = _nSkipRef;

  if( _dutReader.isOpen() && _nEvt >= nEvtBeg ) {
    if( nmiss < 3 ) { // read CMS event:
      ClustDUT = _dutReader.getCMSclust( dutTimeStamp, NpixDUT );
    }
    else { // skip one CMS
      if(  event->getRunNumber() == 6174 && event->getEventNumber() < 112000 ) { // read 2 CMS
	ClustDUT = _dutReader.getCMSclust( dutTimeStamp, NpixDUT );
	ClustDUT = _dutReader.getCMSclust( dutTimeStamp, NpixDUT );
      }
      if(  event->getRunNumber() == 6349 && event->getEventNumber() < 3000 ) { // read 2 CMS
	ClustDUT = _dutReader.getCMSclust( dutTimeStamp, NpixDUT );
	ClustDUT = _dutReader.getCMSclust( dutTimeStamp, NpixDUT );
      }
      if(  event->getRunNumber() == 6352 && event->getEventNumber() < 38000 ) { // read 2 CMS
	ClustDUT = _dutReader.getCMSclust( dutTimeStamp, NpixDUT );
	ClustDUT = _dutReader.getCMSclust( dutTimeStamp, NpixDUT );
      }
      if(  event->getRunNumber() == 6361 && event->getEventNumber() < 29000 ) { // read 2 CMS
	ClustDUT = _dutReader.getCMSclust( dutTimeStamp, NpixDUT );
	ClustDUT = _dutReader.getCMSclust( dutTimeStamp, NpixDUT );
      }
      cout << "==>> event " << event->getEventNumber()
	   << ": nmiss = " << nmiss
	   << ": skip DUT" << endl;
      dutnsync++;
      nmiss = 0; // reset
    }
  }

  if( _refReader.isOpen() && _nEvt >= nRefBeg ) { // read REF event:

    if( kmiss < 3 ) { // read CMS event:
      ClustREF = _refReader.getCMSclust( refTimeStamp, NpixREF );
    }
    else {
      if(  event->getRunNumber() == 5636 &&
	   ( event->getEventNumber() < 66000 || 
	     ( event->getEventNumber() > 166000 && event->getEventNumber() < 255000 ) ) ) {
	ClustREF = _refReader.getCMSclust( refTimeStamp, NpixREF );
	ClustREF = _refReader.getCMSclust( refTimeStamp, NpixREF );
      }
      if(  event->getRunNumber() == 6174 && event->getEventNumber() < 112000 ) { // read 2 CMS
	ClustREF = _refReader.getCMSclust( refTimeStamp, NpixREF );
	ClustREF = _refReader.getCMSclust( refTimeStamp, NpixREF );
      }
      if(  event->getRunNumber() == 6349 && event->getEventNumber() < 3000 ) { // read 2 CMS
	ClustREF = _refReader.getCMSclust( refTimeStamp, NpixREF );
	ClustREF = _refReader.getCMSclust( refTimeStamp, NpixREF );
      }
      if(  event->getRunNumber() == 6352 && event->getEventNumber() < 38000 ) { // read 2 CMS
	ClustREF = _refReader.getCMSclust( refTimeStamp, NpixREF );
	ClustREF = _refReader.getCMSclust( refTimeStamp, NpixREF );
      }
      if(  event->getRunNumber() == 6361 && event->getEventNumber() < 29000 ) { // read 2 CMS
	ClustREF = _refReader.getCMSclust( refTimeStamp, NpixREF );
	ClustREF = _refReader.getCMSclust( refTimeStamp, NpixREF );
      }
      cout << "==>> event " << event->getEventNumber() << ": Re-sync REF" << endl;
      refnsync++;
      kmiss = 0; // reset
    }
  }

  static int nskip = 0;

  if( _isFirstEvent ) {

    nskip = _nSkip; 

    cout << "skip first " << nskip << " DUT events" << endl;

    if( _dutReader.isOpen() ) {
      for( int iskip = 0; iskip < nskip; ++iskip ){
	ClustDUT = _dutReader.getCMSclust( dutTimeStamp, NpixDUT ); // skip first DUT events
      }
    }

    if( _refReader.isOpen() ) {

      int nrefskip = 0;
      cout << "skip first " << nrefskip << " REF events" << endl;

      for( int iskip = 0; iskip < nrefskip; ++iskip ) {
	ClustREF = _refReader.getCMSclust( refTimeStamp, NpixREF ); // skip first REF events
      }

    }

  } // first

  static long dutTimeStamp0 = dutTimeStamp; // init
  static long prevdutTimeStamp = dutTimeStamp; // init
  static long olddutTimeStamp = dutTimeStamp; // init
  static long prevrefTimeStamp = refTimeStamp; // init
  //static long oldrefTimeStamp = refTimeStamp; // init

  static long prevgapTimeStamp = event->getTimeStamp(); // init

  streamlog_out(DEBUG) << setw(6) << setiosflags(ios::right) << event->getEventNumber() << ". DUT pixels " << NpixDUT;
  if( NpixDUT > 0 ) {
    streamlog_out(DEBUG) << ", clusters at";
    for( vector<cluster>::iterator c = ClustDUT.begin(); c != ClustDUT.end(); c++ ){
      streamlog_out(DEBUG) << "  (" << c->col << ", " << c->row << ")";
    }
  }
  streamlog_out(DEBUG) << endl;


  long dutdt = dutTimeStamp - prevdutTimeStamp; // [clocks]
  if( event->getRunNumber() > 4600 && nskip > 0 ) dutdt = prevdutTimeStamp - olddutTimeStamp;
  long refdt = refTimeStamp - prevrefTimeStamp; // [clocks]

  long dutrefddt = dutdt - refdt; // [clocks]
  static long prevdutrefddt = dutrefddt;

  double dutddtns = ( event->getTimeStamp() - prevtime )/gTLU - dutdt / 0.040; // +-40 ns
  double refddtns = ( event->getTimeStamp() - prevtime )/gTLU - refdt / 0.040; // REF

  if( _dutReader.isOpen() ) {

    dutnclusHisto->fill( ClustDUT.size() );

    if( dutTimeStamp > dutTimeStamp0 && event->getTimeStamp() > time0 ) {

      cmsdtHisto->fill( dutdt / 40E0 ); // [us]

      dutrefddtHisto->fill( dutrefddt ); // [clocks]

      // time between events, TLU - DUT:

      dutddtnsHisto->fill( dutddtns ); // [ns]
      refddtnsHisto->fill( refddtns ); // [ns]
      dutddtusHisto->fill( dutddtns/1E3 ); // [us]. peak at 25 us in run 5348
      dutddtmsHisto->fill( dutddtns/1E6 ); // [ms]

      dutddtvst->fill( ( event->getTimeStamp() - time0 ) / fTLU, dutddtns );
      dutddtvsdt->fill( ( event->getTimeStamp() - prevtime ) / gTLU/1E3, dutddtns );

      if( abs( dutddtns ) < 99 ) {
	ddtvst->fill( ( event->getTimeStamp() - time0 ) / fTLU, dutddtns );
	ddtvstms->fill( ( event->getTimeStamp() - time0 ) / fTLU, dutddtns ); // no trend
	ddtvsdt->fill( ( event->getTimeStamp() - prevtime ) / gTLU/1E3, dutddtns );
      }

      if( abs( dutddtns ) > 99 ) {
	long dtgap = event->getTimeStamp() - prevgapTimeStamp;
	gapdtHisto->fill( dtgap / fTLU ); // [s]
	prevgapTimeStamp = event->getTimeStamp();
      }

      if( nskip == 0 ) {
	sysrtHisto->fill( ( event->getTimeStamp() - time0 ) / fTLU /
			  ( dutTimeStamp - dutTimeStamp0 ) * 40E6 ); // mean = 1.00007
	sysrdtHisto->fill( ( event->getTimeStamp() - prevtime) / fTLU /
			   ( dutTimeStamp - prevdutTimeStamp ) * 40E6 );
      }
      else {
	sysrtHisto->fill( ( event->getTimeStamp() - time0 ) / fTLU /
			  ( prevdutTimeStamp - dutTimeStamp0 ) * 40E6 );
	sysrdtHisto->fill( ( event->getTimeStamp() - prevtime) / fTLU /
			   ( prevdutTimeStamp - olddutTimeStamp ) * 40E6 );
      }

      //FIXME run dependend variable leftover. Possible to make more general?
      if( event->getRunNumber() > 4600 && abs( dutddtns ) > 99 ) { // digi DUT missed trigger
	nmiss++;
      }
      else {
	nmiss = 0; // reset
      }


    } // > time0

  } // DUT open

  if( ( event->getTimeStamp() - prevtime ) / fTLU > 0.1 ) { // [s] long gap
    nmiss = 0; // reset
    kmiss = 0;
  }

  if( NpixDUT > 0 ) {

    for( vector<cluster>::iterator c = ClustDUT.begin(); c != ClustDUT.end(); c++ ){

      nclst++;
      dutcolHisto->fill( c->col );
      dutrowHisto->fill( c->row );
      dutnpxHisto->fill( c->size );
      dutadcHisto->fill( c->charge );

    }//DUT clusters

    if( _nEvt > 0 ) {
      logdtcmsHisto->fill( std::log( ( event->getTimeStamp() - prevtime ) / gTLU/1E6 ) / std::log(10.0) );
    }

  }//pix

  if( _refReader.isOpen() ) {
    refnclusHisto->fill( ClustREF.size() );

    if( abs( refddtns ) > 99 ) { // REF missed trigger
      kmiss++;
    }
    else {
      kmiss = 0;
    }
  }

  if( NpixREF > 0 ) {

    for( vector<cluster>::iterator c = ClustREF.begin(); c != ClustREF.end(); c++ ){

      refcolHisto->fill( c->col );
      refrowHisto->fill( c->row );
      refnpxHisto->fill( c->size );
      refadcHisto->fill( c->charge );

    }//REF clusters

  }//pixREF

  olddutTimeStamp = prevdutTimeStamp;
  prevdutTimeStamp = dutTimeStamp;
  prevrefTimeStamp = refTimeStamp;

  _nEvt++;

  prevtime = event->getTimeStamp(); // remember

  //----------------------------------------------------------------------------

  if( _isFirstEvent ) {
 
    // apply all GEAR/alignment offsets to get corrected X,Y,Z position of the
    // sensor center

    for( int iplane = 0; iplane < _siPlanesLayerLayout->getNLayers(); iplane++ ) {

      map< unsigned int , double > _planeCenter;
      map< unsigned int , double > _planeNormal;
 
      int sensorID = _siPlanesLayerLayout->getID(iplane);

      // start filling the map with Gear values:

      _planeCenter[ 0 ] =  _siPlanesLayerLayout->getLayerPositionX(iplane); // X
      _planeCenter[ 1 ] =  _siPlanesLayerLayout->getLayerPositionY(iplane); // Y
      _planeCenter[ 2 ] =  _siPlanesLayerLayout->getLayerPositionZ(iplane); // Z

      _planeNormal[ 0 ] =  0.; // X
      _planeNormal[ 1 ] =  0.; // Y
      _planeNormal[ 2 ] =  1.; // Z

      TVector3  _normalTVec( _planeNormal[0], _planeNormal[1], _planeNormal[2]); 

      // do initial rotation( from GEAR)
      try{
 
	double gRotation[3] = { 0., 0., 0.}; // not rotated
                         
	gRotation[0] = _siPlanesLayerLayout->getLayerRotationXY(iplane); // Euler alpha
	gRotation[1] = _siPlanesLayerLayout->getLayerRotationZX(iplane); // Euler alpha
	gRotation[2] = _siPlanesLayerLayout->getLayerRotationZY(iplane); // Euler alpha
                          
	// input angles are in DEGREEs !!!
	// translate into radians
	gRotation[0] =  gRotation[0]*3.1415926/180.; // 
	gRotation[1] =  gRotation[1]*3.1415926/180.; //
	gRotation[2] =  gRotation[2]*3.1415926/180.; //

	TRotation	r;
	r.RotateX( gRotation[2] );                          
	r.RotateY( gRotation[1] );                          
	r.RotateZ( gRotation[0] );                          

	_normalTVec.Transform( r );
                                  
	_planeNormal[0] = _normalTVec[0];
	_planeNormal[1] = _normalTVec[1];
	_planeNormal[2] = _normalTVec[2];
      }
      catch(...) {
	printf(" no sensor rotation is given in the GEAR steering file, assume NONE \n" );
      }
 
      printf("sensorID: %5d    %8.3f %8.3f %8.3f \n", sensorID, _planeCenter[ 0 ], _planeCenter[ 1 ], _planeCenter[ 2 ] );
             
      _siPlaneCenter[ sensorID]=  _planeCenter;
      _siPlaneNormal[ sensorID]=  _planeNormal;
    }

    if( isFirstEvent() ) _isFirstEvent = false;

  }//isFirstEvent

  //----------------------------------------------------------------------------
  // check input collection (aligned hits):

  LCCollection* col;
  try {
    col = event->getCollection( _inputCollectionTelescope );
  }
  catch( lcio::DataNotAvailableException& e) {
    streamlog_out( DEBUG ) << "Not able to get collection "
			   << _inputCollectionTelescope
			   << "\nfrom event " << event->getEventNumber()
			   << " in run " << event->getRunNumber()  << endl;

    return;
  }

  //----------------------------------------------------------------------------
  // Copy hits to local table
  // Assign hits to sensor planes

  int nHit = col->getNumberOfElements();  //  from event collection

  if( debug ) {
    streamlog_out( TESTFITTERMESSAGE )  << "Total of " << nHit
					<< " tracker hits in input collection " << endl;
  }

  nAllHitHisto->fill(nHit);

  //----------------------------------------------------------------------------

  double * hitX  = new double[nHit];
  double * hitEx = new double[nHit];
  double * hitY  = new double[nHit];
  double * hitEy = new double[nHit];
  double * hitZ  = new double[nHit];
  int    * hitPlane = new int[nHit];

  IntVec * planeHitID = new IntVec[_nTelPlanes];

  // ####################### Loop over telescope hits ##############################
  // This basically extracts all hits from the LCIO collection and
  // writes them into a array (hitX, hitY, hitZ) and calculates the
  // uncertainty.

  int nGoodHit = 0;

  for( int ihit = 0; ihit < nHit; ihit++ ) {

    TrackerHit * meshit = dynamic_cast<TrackerHit*>( col->getElementAt(ihit) );

    // Telescope hit position:
    const double * pos = meshit->getPosition();

    hitZ[ihit] = pos[2];

    // We have to find Plane ID of the hit by looking at the Z position

    double distMin = 99; // [mm]
    hitPlane[ihit] = -1;

    for( int ipl = 0; ipl < _nTelPlanes; ipl++ ) {

      double dist =  hitZ[ihit] - _planePosition[ipl];
      if( dist < 0 ) dist = -dist; // always positive defined !

      if( dist < distMin ) {
        hitPlane[ihit] = ipl;
        distMin        = dist;
      }
    }//planes

    hitX[ihit] = pos[0];
    hitY[ihit] = pos[1];

    // Add hit to hit list for given plane - to be used in track selection:
    planeHitID[hitPlane[ihit]].push_back(ihit);
    nGoodHit++;

    // Position uncertainty. Use nominal resolution if not properly defined:
    const EVENT::FloatVec cov = meshit->getCovMatrix();

    if( cov.at(0) > 0. ) {
      hitEx[ihit]=sqrt(cov.at(0));
    }
    else {
      hitEx[ihit]=_planeResolution[hitPlane[ihit]];
    }
    if( cov.at(2) > 0. ) {
      hitEy[ihit]=sqrt(cov.at(2));
    }
    else {
      hitEy[ihit]=_planeResolution[hitPlane[ihit]];
    }

    streamlog_out(DEBUG) << "Hit " << ihit
			    << "   X = " << hitX[ihit] << " +/- " << hitEx[ihit]
			    << "   Y = " << hitY[ihit] << " +/- " << hitEy[ihit]
			    << "   Z = " << hitZ[ihit] << "( plane " << hitPlane[ihit] << ")" << endl;

  }//loop over hits

  streamlog_out(DEBUG) << "event " << event->getEventNumber()
			 << ", hits " << nHit
			 << ", good " << nGoodHit << endl;
  
  hits0Histo->fill( planeHitID[0].size() );
  hits1Histo->fill( planeHitID[1].size() );
  hits2Histo->fill( planeHitID[2].size() );
  hits3Histo->fill( planeHitID[3].size() );
  hits4Histo->fill( planeHitID[4].size() );
  hits5Histo->fill( planeHitID[5].size() );

  //###########################################################################


  //####################### loop over telescope hit pairs #####################
  // This only fills some correlation histograms.

  for( int ihit = 0; ihit < nHit; ihit++ ) {

    int ipl = hitPlane[ihit];

    for( int jhit = 0; jhit < nHit; jhit++ ) {

      int jpl = hitPlane[jhit];

      double dx = hitX[jhit] - hitX[ihit];
      double dy = hitY[jhit] - hitY[ihit];

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

      if( ipl == 0 ) {

	if( jpl == 1 ) {
	  dx01Histo->fill( dx );
	  dy01Histo->fill( dy );
	  if( abs(dy) < 1 ) du01Histo->fill( dx );
	}
	if( jpl == 2 ) {
	  dx02Histo->fill( dx );
	}
	if( jpl == 3 ) {
	  dx03Histo->fill( dx );
	}
	if( jpl == 4 ) {
	  dx04Histo->fill( dx );
	}
	if( jpl == 5 ) {
	  dx05Histo->fill( dx );
	}

      }//ipl==0

      if( ipl == 1 ) {

	if( jpl == 2 ) {
	  dx12Histo->fill( dx );
	  dy12Histo->fill( dy );
	  if( abs(dy) < 1 ) du12Histo->fill( dx );
	}
      }

      if( ipl == 2 ) {

	if( jpl == 3 ) {
	  dx23Histo->fill( dx );
	  dy23Histo->fill( dy );
	  if( abs(dy) < 1 ) du23Histo->fill( dx );
	}
      }

      if( ipl == 3 ) {

	if( jpl == 4 ) {
	  dx34Histo->fill( dx );
	  dy34Histo->fill( dy );
	  if( abs(dy) < 1 ) du34Histo->fill( dx );
	}
      }

      if( ipl == 4 ) {

	if( jpl == 5 ) {
	  dx45Histo->fill( dx );
	  dy45Histo->fill( dy );
	  if( abs(dy) < 1 ) du45Histo->fill( dx );
	}
      }

#endif

    }//loop over telescope hits

  }// loop over telescope hits

  //##########################################################################


  //----------------------------------------------------------------------------
  // DUT alignment: relative to _planePosition[2]
  // REF alignment: relative to _planePosition[5]

  bool lDUT = 1; // DUT present
  DUTalignx = _DUTalignx;
  DUTaligny = _DUTaligny;
  DUTz = _DUTz + _planePosition[2]; // default

  DUTrot = _DUTrot;
  tilt = _DUTtilt;
  turn = _DUTturn;

  double wt = atan(1.0) / 45.0; // pi/180 deg
  double DUTX0 = 0.07; // init

  double REFz = _REFz + _planePosition[5]; // timing reference plane
  double REFalignx =  _REFalignx;
  double REFaligny = _REFaligny;
  double REFrot = _REFrot;

  DUTX0 = DUTX0 / cos( tilt * wt ); // for flat prob in run 158

  double pitchcol = 0.150; // [mm]
  double pitchrow = 0.100; // [mm]
  double pitchrowt = 0.100 * cos( tilt * wt ); // [mm] in telescope system

  // normal vector on DUT surface:
  // N = ( 0, 0, -1 ) on DUT
  // transform into tele system:
  // tilt alpha around x
  // turn omega around y

  double Nx =-sin( turn*wt )*cos( tilt*wt );
  double Ny = sin( tilt*wt );
  double Nz =-cos( turn*wt )*cos( tilt*wt );

  double co = cos( turn*wt );
  double so = sin( turn*wt );
  double ca = cos( tilt*wt );
  double sa = sin( tilt*wt );
  double cf = cos( DUTrot );
  double sf = sin( DUTrot );

  double dNxdo =-co*ca;
  double dNxda = so*sa;
  double dNydo = 0;
  double dNyda = ca;
  double dNzdo = so*ca;
  double dNzda = co*sa;

  double path = sqrt( 1 + tan(turn*wt)*tan(turn*wt) + tan(tilt*wt)*tan(tilt*wt) );
  double norm = cos( turn*wt ) * cos( tilt*wt ); // length of Nz. Landau 21.6

  //norm = 1/path; // test. Landau 22.4

  //FIXME very hacky way of determining the coordinate transformation mode:
  bool transform = 0; // old projection style with manual align
  if( abs(tilt) > 2 && abs(turn) > 2 )
    transform = 1; // new transformation style with 6-D MillePede align


  //################### 3-4-5 downstream ("driplets") ###########################
  // This runs over all hits in plane 3 and 5 of the telescope and
  // tries to match triplets by comparing with the middle plane (plane
  // 4). if a triplet is found it is extrapolated to both REF and DUT planes.

  double driCut = 0.1; // [mm]

  int ndri = 0;
  int ndriref = 0;
  int ndrilk = 0;
  double xmB[99];
  double ymB[99];
  double zmB[99];
  double sxB[99];
  double syB[99];
  bool lkB[99];
  int hts[6][99]; // 6 planes

  for( int ihit = 0; ihit < nHit; ihit++ ) {

    if( hitPlane[ihit] != 3 ) continue; // want plane 3

    for( int jhit = 0; jhit < nHit; jhit++ ) {

      if( hitPlane[jhit] != 5 ) continue; //want last plane 5

      double dx35 = hitX[jhit] - hitX[ihit];
      double dy35 = hitY[jhit] - hitY[ihit];

      if( abs(dy35) < 1 ) dx35Histo->fill( dx35 );
      if( abs(dx35) < 1 ) dy35Histo->fill( dy35 );

      double dz35 = hitZ[jhit] - hitZ[ihit];

      if( abs(dy35) > 0.010 * dz35 ) continue; // [mm]
      if( abs(dx35) > 0.010 * dz35 ) continue; // cuts on track angle:+-10 mrad

      double avx = 0.5 * ( hitX[ihit] + hitX[jhit] );
      double avy = 0.5 * ( hitY[ihit] + hitY[jhit] );
      double avz = 0.5 * ( hitZ[jhit] + hitZ[ihit] ); // mid z

      double tx = ( hitX[jhit] - hitX[ihit] ) / dz35; // angle theta x
      double ty = ( hitY[jhit] - hitY[ihit] ) / dz35; // angle theta y

      // extrapolate to REF timing plane: [mm]

      // telescope coordinates (right handed, look along beam):
      //   x_tele = inward (towards DESY)
      //   y_tele = down
      //   z = beam
      // REF coordinates (mounted upside down in old black socket carrier):
      //   col = -x = outward
      //   row =  y = down
      // DP coordinates (right handed, look at beam):
      //   x_DP = outward (away from DESY)
      //   y_DP = up
      //   z = beam
      //   x_DP = -x_tele
      //   y_DP = -y_tele, together = 180 deg rot

      double zB = REFz - avz; // z REF from mid of driplet
      double xB = avx + tx * zB; // driplet impact point on REF
      double yB = avy + ty * zB; // tele coord

      // transform into REF system:

      double xBm = -xB + REFalignx; // invert x, shift to mid
      double yBm =  yB + REFaligny; // shift to mid
      double xBr = xBm - DUTrot*yBm; // rotate CMS pixel
      double yBr = yBm + DUTrot*xBm;
      //double xBt = xBr + 0.075; // shift by half bin?
      double xBt = xBr;
      double yBt = yBr; // no tilt for REF

      double xmod = fmod( 20 + xBt, 0.3 ) * 1E3; // [0,300] um, 2 pixel wide
      double ymod = fmod( 20 + yBt, 0.2 ) * 1E3; // [0,200] um

      // middle plane k for driplets:

      for( int khit = 0; khit < nHit; khit++ ) {

	if( hitPlane[khit] != 4 ) continue; //want middle plane 4

	//interpolate to plane k:

	double dz = hitZ[khit] - avz;
	double xk = avx + tx * dz; // driplet at k
	double yk = avy + ty * dz;

	double dx = hitX[khit] - xk;
	double dy = hitY[khit] - yk;

	if( abs(dy) < driCut ) dridxHisto->fill( dx*1E3 ); //sig = 7 um at 5 GeV
	if( abs(dx) < driCut ) dridyHisto->fill( dy*1E3 );

	// alignment check:

	if( abs(dy) < driCut ) {
	  dridxvsx->fill( avx, dx*1E3 ); // check for rot
	  dridxvsy->fill( avy, dx*1E3 );
	  dridxvstx->fill( tx*1E3, dx*1E3 ); // check for z shift
	  dridxvsty->fill( ty*1E3, dx*1E3 );
	}
	if( abs(dx) < driCut ) {
	  dridyvsx->fill( avx, dy*1E3 );
	  dridyvsy->fill( avy, dy*1E3 );
	  dridyvstx->fill( tx*1E3, dy*1E3 );
	  dridyvsty->fill( ty*1E3, dy*1E3 );
	}

	// driplet found:

	if( abs(dx) < driCut && abs(dy) < driCut ){

	  if( ndri < 99 ) {
	    xmB[ndri] = avx;
	    ymB[ndri] = avy;
	    zmB[ndri] = avz;
	    sxB[ndri] = tx;
	    syB[ndri] = ty;
	    lkB[ndri] = 0; //no link to CMS yet
	    hts[3][ndri] = ihit;
	    hts[5][ndri] = jhit;
	    hts[4][ndri] = khit;
	  }

	  ndri++;

	  drixHisto->fill( -hitX[khit] ); // -x = x_DP = out
	  driyHisto->fill( -hitY[khit] ); // -y = y_DP =  up
	  drixyHisto->fill( -hitX[khit], -hitY[khit] );
	  dritxHisto->fill( tx*1E3 );
	  drityHisto->fill( ty*1E3 );

	  drixrefHisto->fill( -xB ); // -xB = x_DP = out
	  driyrefHisto->fill( -yB ); // -yB = y_DP = up
	  drixyrefHisto->fill( -xB, -yB );

	  // correlate with REF:

	  if( NpixREF > 0 ){

	    ndriref++;

	    for( vector<cluster>::iterator c = ClustREF.begin(); c != ClustREF.end(); c++ ){

	      refxxHisto->fill( c->col, xB ); // anti-correlation
	      refyyHisto->fill( c->row, yB ); // correlation

	      double refx = ( c->col - 26 ) * 0.15; // -3.9..3.9 mm
	      double refy = ( c->row - 40 ) * 0.10; // -4..4 mm

	      refsxaHisto->fill( refx + xB );
	      refdyaHisto->fill( refy - yB );

	      double refsx = refx + REFrot*refy - REFalignx + xB; // residual x
	      double refdy = refy - REFrot*refx - REFaligny - yB; // residual y

	      refsxHisto->fill( refsx );
	      refdyHisto->fill( refdy );

	      if( abs( refdy ) < 1 ) refsxcHisto->fill( refsx );
	      if( abs( refsx ) < 1 ) refdycHisto->fill( refdy );

	      if( abs( refdy ) < 1 && abs( refsx ) < 1 ) {

		drixlkHisto->fill( -xB ); // -xB = x_DP = out
		driylkHisto->fill( -yB ); // -yB = y_DP = up
		drixylkHisto->fill( -xB, -yB );
		refpixvsxmym->fill( xmod, ymod ); // occupancy map

		refqHisto->fill( c->charge );
		refqvsxmym->fill( xmod, ymod, c->charge ); // cluster charge profile

		refdyvsx->fill( refx, refdy*1E3 ); // dy vs x: rot?
		refdyvsy->fill( refy, refdy*1E3 ); // dy vs y: tilt?
		refdyvsty->fill( ty*1E3, refdy*1E3 ); // dy vs tety: z shift?

	      }

	      if( abs( refdy ) < 0.3 && abs( refsx ) < 0.3 ) {
		if( ndri < 99 ) lkB[ndri-1] = 1; // link to CMS
		reflkcolHisto->fill( c->col );
		reflkrowHisto->fill( c->row );
		ndrilk++;
	      }

	    }//loop over CMS clusters

	  }//have CMS cluster

	}//have telescope driplet

	// correlate 3-4-5 driplet with DUT:

	if( NpixDUT > 0 && abs(dx) < driCut && abs(dy) < driCut ){

	  //extrapolate to DUT plane: [mm]

	  double zs = DUTz - avz; // z DUT from mid of driplet
	  double xs = avx + tx * zs; // driplet impact point on DUT
	  double ys = avy + ty * zs;

	  for( vector<cluster>::iterator c = ClustDUT.begin(); c != ClustDUT.end(); c++ ){

	    double bacx = ( c->col - 26 ) * pitchcol; // -3.9..3.9 mm
	    double bacy = ( c->row - 40 ) * pitchrowt; // -4..4 mm

	    bacsxaHisto->fill( bacx + xs );
	    bacdyaHisto->fill( bacy - ys );

	    double bacsx = bacx + DUTrot*bacy - DUTalignx + xs; // residual x
	    double bacdy = bacy - DUTrot*bacx - DUTaligny - ys; // residual y

	    if( _dutReader.getChip() >= 100 ) { // xdb and dig mounted upright
	      bacsx = bacx + DUTrot*bacy - DUTalignx - xs; // residual x
	      bacdy = bacy - DUTrot*bacx - DUTaligny + ys; // residual y
	    }

	    if( abs( bacdy ) < 1 ) bacsxcHisto->fill( bacsx*1E3 );
	    if( abs( bacsx ) < 1 ) bacdycHisto->fill( bacdy*1E3 );

	    if( c->charge * norm > 18 ){
	      if( abs( bacdy ) < 1 ) bacsxcqHisto->fill( bacsx*1E3 );
	      if( abs( bacsx ) < 1 ) bacdycqHisto->fill( bacdy*1E3 ); // 18 um @ 5 GeV
	    }

	  }//loop over BAC clusters

	}//have telescope driplet

      }//loop over hits

    }//loop over hits

  }// loop over hits

  ndriHisto->fill( ndri );
  if( NpixDUT > 0 ) ndrirefHisto->fill( ndriref );
  lkBvst->fill( (event->getTimeStamp()-time0)/fTLU, ndrilk );

  // ##########################################################################


  //##########################  0-1-2 upstream triplets #######################
  // here again a telescope hit triplet is formed, from planes 0 and 2
  // and the correlation to plane 1. the found triplets are
  // extrapolated and matched to the DUT plane and exhaustive
  // histograms are written.
  // Furthermore the hits of a triplet matched to the DUT are
  // collected and fed into GBL for a track fit.

  double triCut = 0.1; // [mm]

  int ntri = 0;
  int ntrilk = 0;
  double xmA[99];
  double ymA[99];
  double zmA[99];
  double sxA[99];
  double syA[99];
  bool lkA[99];
  double DUTsxA[99];
  double DUTdyA[99];

  for( int ihit = 0; ihit < nHit; ihit++ ) { //hit in 0

    if( hitPlane[ihit] > 0 ) continue; //want plane 0

    for( int jhit = 0; jhit < nHit; jhit++ ) { // hit in 2

      if( hitPlane[jhit] != 2 ) continue; //want plane 2

      double dx02 = hitX[jhit] - hitX[ihit];
      double dy02 = hitY[jhit] - hitY[ihit];

      if( abs(dy02) < 1 ) da02Histo->fill( dx02 );
      if( abs(dx02) < 1 ) db02Histo->fill( dy02 );

      double dz02 = hitZ[jhit] - hitZ[ihit]; // lever arm from 0 to 2

      if( abs(dy02) > 0.005 * dz02 ) continue; // [mm]
      if( abs(dx02) > 0.005 * dz02 ) continue; // cuts on track angle:+-5 mrad

      double avx = 0.5 * ( hitX[ihit] + hitX[jhit] ); // average
      double avy = 0.5 * ( hitY[ihit] + hitY[jhit] ); // = mid point
      double avz = 0.5 * ( hitZ[jhit] + hitZ[ihit] ); // mid z

      double tx = ( hitX[jhit] - hitX[ihit] ) / dz02; // angle theta x
      double ty = ( hitY[jhit] - hitY[ihit] ) / dz02; // angle theta x

      //extrapolate to CMS: [mm]

      // telescope coordinates (right handed, look along beam):
      //   x_tele = inward (towards DESY)
      //   y_tele = down
      //   z = beam
      // REF coordinates (mounted upside down in old black socket carrier):
      //   col = -x = outward
      //   row =  y = down
      // DP coordinates (right handed, look at beam):
      //   x_DP = outward (away from DESY)
      //   y_DP = up
      //   z = beam
      //   x_DP = -x_tele
      //   y_DP = -y_tele, together = 180 deg rot around z

      double zA = DUTz - avz; // z CMS from mid of triplet
      double xA = avx + tx * zA; // triplet impact point on CMS
      double yA = avy + ty * zA;
      xA = xA / ( 1.0 - tx * tan(turn*wt) ); // turn correction
      yA = yA / ( 1.0 - ty * tan(tilt*wt) ); // tilt correction

      // intersect inclined track with tilted DUT plane:

      double zc = (Nz*zA - Ny*avy - Nx*avx) / (Nx*tx + Ny*ty + Nz); // from avz
      double yc = avy + ty * zc;
      double xc = avx + tx * zc;

      double dxcdz = tx*Nz / (Nx*tx + Ny*ty + Nz);
      double dycdz = ty*Nz / (Nx*tx + Ny*ty + Nz);

      double dzc = zc + avz - DUTz; // from DUT z0 [-8,8] mm

      double dzcda = ( dNzda*zA - dNyda*avy - dNxda*avx ) / (Nx*tx + Ny*ty + Nz) - zc*( dNxda*tx + dNyda*ty + dNzda ) / (Nx*tx + Ny*ty + Nz);
      double dzcdo = ( dNzdo*zA - dNydo*avy - dNxdo*avx ) / (Nx*tx + Ny*ty + Nz) - zc*( dNxdo*tx + dNydo*ty + dNzdo ) / (Nx*tx + Ny*ty + Nz);

      double ddzcdz = Nz / (Nx*tx + Ny*ty + Nz) - 1; // ddzc/dDUTz

      dzcvsxy->fill( xc, yc, dzc ); // tilted plane

      // transform into DUT system: (passive).
      // large rotations don't commute: careful with order

      double x1 = co*xc - so*dzc; // turn
      double y1 = yc;
      double z1 = so*xc + co*dzc;

      double x2 = x1;
      double y2 = ca*y1 + sa*z1; // tilt
      double z2 =-sa*y1 + ca*z1; // should be zero (in DUT plane). is zero

      double x3 = cf*x2 + sf*y2; // rot
      double y3 =-sf*x2 + cf*y2;
      double z3 = z2; // should be zero (in DUT plane). is zero

      z3vsxy->fill( xc, yc, z3 ); // should be zero. is zero

      double upsign = -1; // analog psi46: upside down
      if( _dutReader.getChip() >= 100 )
	upsign = 1; // xdb and dig mounted upright
      double x4 = upsign*x3 + DUTalignx; // shift to mid
      double y4 =-upsign*y3 + DUTaligny; // invert y, shift to mid

      //cout << "UPSIGN " << upsign << endl;
      //cout << "X4 " << x4 << endl;
      //cout << "DUTalign " << DUTalignx << endl;
      //printf("MONITOR: x4 = %f*%f + %f = %f\n",upsign,x3,DUTalignx,x4);

      // for tilted: track at DUT:

      double xAm =-xA + DUTalignx; // invert x, shift to mid
      double yAm = yA + DUTaligny; // shift to mid
      if( _dutReader.getChip() >= 100 ) { // xdb and dig mounted upright
	xAm = xA + DUTalignx; // shift to mid
	yAm =-yA + DUTaligny; // invert y, shift to mid
      }
      double xAr = xAm - DUTrot*yAm; // rotate CMS pixel
      double yAr = yAm + DUTrot*xAm;
      double xAt = xAr / cos(turn*wt); // turn in x
      double yAt = yAr / cos(tilt*wt); // tilt in y

      if( transform ) { // new style
	xAt = x4;
	yAt = y4;
      }

      // reduce to 2x2 pixel region:

      double xmod = fmod( 9.075 + xAt, 0.3 ) * 1E3; // [0,300] um, 2 pixel wide
      double ymod = fmod( 9.050 + yAt, 0.2 ) * 1E3; // [0,200] um
      double ymd3 = fmod( 9.050 + yAt, 0.3 ) * 1E3; // [0,300] um
      double ymd6 = fmod( 9.050 + yAt, 0.6 ) * 1E3; // [0,600] um

      // 0 deg:
      bool ldot = 1; // bias dot, from cmsqvsxmym
      if( xmod < 105 ) ldot = 0; // dot at x = 125
      if( xmod > 195 ) ldot = 0; // and at x = 175
      if( tilt < 6 ) {
	if( ymod <  55 ) ldot = 0; // dot at y =  75
	if( ymod > 195 ) ldot = 0; // dot at y = 175
	if( ymod >  95 && ymod < 155 ) ldot = 0; // band between dots
      }

      bool lcore = 1; // pixel core, 2x2 pixel region
      if( xmod <  20 ) lcore = 0; // outer edge, see cmsncolvsxm
      if( xmod > 280 ) lcore = 0; // outer edge
      if( ymod <  20 ) lcore = 0; // outer edge, see cmsnrowvsym
      if( ymod > 180 ) lcore = 0; // outer edge
      if( xmod > 130 && xmod < 170 ) lcore = 0; // inner edge
      if( ymod >  80 && ymod < 120 ) lcore = 0; // inner edge

      // middle plane k for triplets:

      for( int khit = 0; khit < nHit; khit++ ) { //hit in 1

	if( hitPlane[khit] != 1 ) continue; //want middle plane

	//interpolate to plane k:

	double dz = hitZ[khit] - avz;
	double xk = avx + tx * dz; // driplet at k
	double yk = avy + ty * dz;

	double dx = hitX[khit] - xk;
	double dy = hitY[khit] - yk;

	if( abs(dy) < triCut ) tridxHisto->fill( dx*1E3 );
	if( abs(dx) < triCut ) tridyHisto->fill( dy*1E3 );

	if( abs(dy) < triCut ) tridx1Histo->fill( dx*1E0 );
	if( abs(dx) < triCut ) tridy1Histo->fill( dy*1E0 );

	// alignment check:

	if( abs(dy) < triCut ) {
	  tridxvsx->fill( avx, dx*1E3 ); // check for rot
	  tridxvsy->fill( avy, dx*1E3 );
	  tridxvstx->fill( tx*1E3, dx*1E3 ); // check for z shift
	  tridxvsty->fill( ty*1E3, dx*1E3 );
	}
	if( abs(dx) < triCut ) {
	  tridyvsx->fill( avx, dy*1E3 );
	  tridyvsy->fill( avy, dy*1E3 );
	  tridyvstx->fill( tx*1E3, dy*1E3 );
	  tridyvsty->fill( ty*1E3, dy*1E3 );
	}

	// triplet found:

	bool ltri = 0;

	if( abs(dx) < triCut && abs(dy) < triCut ) {

	  if( ntri < 99 ) {
	    ltri = 1;
	    xmA[ntri] = avx;
	    ymA[ntri] = avy;
	    zmA[ntri] = avz;
	    sxA[ntri] = tx;
	    syA[ntri] = ty;
	    lkA[ntri] = 0; //no link to CMS yet
	    hts[0][ntri] = ihit;
	    hts[2][ntri] = jhit;
	    hts[1][ntri] = khit;
	  }

	  ntri++;

	  trixHisto->fill( -hitX[khit] ); // -x = x_DP = out
	  triyHisto->fill( -hitY[khit] ); // -y = y_DP =  up
	  trixyHisto->fill( -hitX[khit], -hitY[khit] );
	  tritxHisto->fill( tx*1E3 );
	  trityHisto->fill( ty*1E3 );

	  trixdutHisto->fill( -xA ); // -xA = x_DP = out
	  triydutHisto->fill( -yA ); // -yA = y_DP = up
	  trixydutHisto->fill( -xA, -yA );

	  // scale hit to CMS pixel:
	  // sizeX="21.2"  sizeY="10.6"
	  // npixelX="1152"  npixelY="576" 

	  streamlog_out(DEBUG) << "  x " << int( ( hitX[jhit] + 10.6 ) / 0.15 )
		 << ", y " << int( ( hitY[jhit] +  5.3 ) / 0.10 );
	  
	  // extrapolate to planes 3,4,5: resolution study

	  for( int lhit = 0; lhit < nHit; lhit++ ) {

	    if( hitPlane[lhit] <= 2 ) continue; //want 3,4, or 5

	    double dz = hitZ[lhit] - avz;
	    double xs = avx + tx * dz; // triplet impact point on CMS
	    double ys = avy + ty * dz;
	    double dx = hitX[lhit] - xs;
	    double dy = hitY[lhit] - ys;

	    if( hitPlane[lhit] == 3 ) {
	      tridx3Histo->fill( dx*1E0 );
	      tridy3Histo->fill( dy*1E0 ); // 65 um at 4.7 GeV with CMS
	      tridx3bHisto->fill( dx*1E3 ); // finer binning
	      tridy3bHisto->fill( dy*1E3 ); // 
	    }
	    else if( hitPlane[lhit] == 4 ) {
	      tridx4Histo->fill( dx*1E0 );
	      tridy4Histo->fill( dy*1E0 ); //174 um at 4.7 GeV
	      tridx4bHisto->fill( dx*1E3 ); // finer binning
	      tridy4bHisto->fill( dy*1E3 ); // 
	    }
	    else if( hitPlane[lhit] == 5 ) {
	      tridx5Histo->fill( dx*1E0 );
	      tridy5Histo->fill( dy*1E0 ); //273 um at 4.7 GeV
	      tridx5bHisto->fill( dx*1E3 ); // finer binning
	      tridy5bHisto->fill( dy*1E3 ); // 
	    }

	  }//lhit in 3,4,5

	  //----------------------------------------------------------------------
	  // correlate with CMS DUT:
	  if( NpixDUT > 0 ) {

	    // CMS pixel clusters:

	    for( vector<cluster>::iterator c = ClustDUT.begin(); c != ClustDUT.end(); c++ ){

	      cmsxxHisto->fill( c->col, xA ); // anti-correlation x-x
	      cmsyyHisto->fill( c->row, yA ); // correlation y-y

	      // pix in clus:

	      int colmin = 99;
	      int colmax = -1;
	      int rowmin = 99;
	      int rowmax = -1;

	      double qcol[52];
	      for( int icol = 0; icol < 80; ++icol ) qcol[icol] = 0;

	      double qrow[80];
	      for( int irow = 0; irow < 80; ++irow ) qrow[irow] = 0;

	      for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); px++ ){

		cmspxqHisto->fill( px->q );
		qcol[px->col] += fabs(px->q); // project cluster onto cols
		qrow[px->row] += fabs(px->q); // project cluster onto rows
		if( px->col < colmin ) colmin = px->col;
		if( px->col > colmax ) colmax = px->col;
		if( px->row < rowmin ) rowmin = px->row;
		if( px->row > rowmax ) rowmax = px->row;

	      }//pix

	      bool fiducial = 1;
	      if( rowmin ==  0 ) fiducial = 0;
	      if( rowmax == 79 ) fiducial = 0;
	      if( colmin ==  0 ) fiducial = 0;
	      if( colmax == 51 ) fiducial = 0;
	      if( _dutReader.getChip() == 7 ) {
		if( colmax == 25 ) fiducial = 0; // 26 dead
		if( colmin == 28 ) fiducial = 0; // 27 dead
		if( colmax == 35 ) fiducial = 0; // 36 dead
		if( colmin == 38 ) fiducial = 0; // 37 dead
	      }

	      int ncol = colmax - colmin + 1;
	      int nrow = rowmax - rowmin + 1;

	      double a1 = 0; //1st
	      double a2 = 0; //2nd
	      int i1 = 99;
	      int i2 = 99;

	      // eta-algo in rows:

	      for( int irow = 0; irow < 80; ++irow ) {
		if( qrow[irow] > a1 ) {
		  a2 = a1;
		  a1 = qrow[irow];
		  i2 = i1;
		  i1 = irow;
		}
		else if( qrow[irow] > a2 ) {
		  a2 = qrow[irow];
		  i2 = irow;
		}
	      }
	      double a12 = a1 + a2;
	      double eta = 0;
	      if( a12 > 1 ) eta = ( a1 - a2 ) / a12;
	      if( i1 > i2 ) eta = -eta;

	      // head-tail in cols:

	      if( ncol > 2 ) {
		double colmid = colmin + 0.5*(ncol-1);
		double qht  = qcol[colmin] + qcol[colmax];
		double asy = ( qcol[colmax] - qcol[colmin] ) / qht; // -1..1
		double range = 0.5; // range 0..1. 0.5 is best at 46 deg turn
		c->col = colmid + range * 0.5*(ncol-1) * asy;
		//printf("MONITOR: c->col = %f, c->row - %f\n",c->col,c->row);
	      }

	      bool lq = 0;
	      double Q0 = c->charge * norm; // cluster charge normalized to vertical incidence
	      if( Q0 > 18 &&  Q0 < 35 ) lq = 1;

	      // DUT - triplet:

	      double cmsx = ( c->col - 26 ) * pitchcol; // -3.9..3.9 mm
	      double cmsy = ( c->row - 40 ) * pitchrow; // -4..4 mm
	      double cmsyt = ( c->row - 40 ) * pitchrowt; // -4..4 mm

	      //printf("MONITOR: cmsx = %f, cmsy - %f\n",cmsx,cmsy);

	      //FIXME check this carefully, "new style" is true 3D
	      //transformation while "old style is just rotation somehow
	      if( transform ) { // new style

		cmssxaHisto->fill( cmsx + x3 ); // rot, tilt and turn but no shift
		cmsdxaHisto->fill( cmsx - x3 );

		cmssyaHisto->fill( cmsy + y3 );
		cmsdyaHisto->fill( cmsy - y3 );
	      }
	      else {

		cmssxaHisto->fill( cmsx + xA ); // rot and tilt or turn but no shift
		cmsdxaHisto->fill( cmsx - xA );

		cmssyaHisto->fill( cmsy + yA );
		cmsdyaHisto->fill( cmsy - yA );
	      }

	      double dx4 = cmsx - x4;
	      double dy4 = cmsy - y4;

	      //printf("MONITOR: x4 = %f, y4 = %f\n",x4,y4);
	      //printf("MONITOR: cmsdx = %f, cmsdy = %f\n",dx4,dy4);

	      cmsdx4Histo->fill( dx4 );
	      cmsdy4Histo->fill( dy4 );

	      double cmssx = cmsx + DUTrot*cmsy - DUTalignx + xA; // residual x
	      double cmsdy = cmsyt - DUTrot*cmsx - DUTaligny - yA; // residual y

	      if( _dutReader.getChip() >= 100 ) { // xdb and dig mounted upright
		cmssx = cmsx + DUTrot*cmsy - DUTalignx - xA; // residual x
		cmsdy = cmsyt - DUTrot*cmsx - DUTaligny + yA; // residual y
	      }

	      double cmsdx = dx4;
	      if( transform ) cmsdy = dy4; // new style

	      if( _dutReader.getChip() >= 200 ) { // even/odd col effect for dig
		int iodd = floor( fmod( c->col, 2 ) );
		if( iodd ) // odd
		  cmsdy = cmsdy + 1.5E-3;
		else
		  cmsdy = cmsdy - 1.5E-3;
	      }

	      cmsdxHisto->fill( cmsdx*1E3 );
	      cmsdyHisto->fill( cmsdy*1E3 );

	      if( fiducial ) {

		cmsdxfHisto->fill( cmsdx*1E3 );
		cmsdyfHisto->fill( cmsdy*1E3 );

		if( abs( cmsdy ) < 0.10 ) cmsdxfcHisto->fill( cmsdx*1E3 );
		if( abs( cmsdx ) < 0.15 ) cmsdyfcHisto->fill( cmsdy*1E3 );

	      }//CMS fiducial

	      // accumulate cuts for y:

	      if( fiducial && abs( cmsdx ) < 0.15 ) {

		if(      nrow == 1 )
		  cmsdyfc1Histo->fill( cmsdy*1E3 ); // 3972: 7.7
		else if( nrow == 2 )
		  cmsdyfc2Histo->fill( cmsdy*1E3 ); // 3972: 9.5
		else
		  cmsdyfc3Histo->fill( cmsdy*1E3 ); // 3872: 68

		if(      Q0 < 18 )
		  cmsdyq0Histo->fill( cmsdy*1E3 );

		else if( Q0 < 40 ){

		  cmsdyq1Histo->fill( cmsdy*1E3 );
		  if( eta < 0 ) {
		    cmsdyeta0Histo->fill( cmsdy*1E3 );
		  }
		  else {
		    cmsdyeta1Histo->fill( cmsdy*1E3 );
		  }
		}
		else
		  cmsdyq2Histo->fill( cmsdy*1E3 );

		if( abs( ty-0.000 ) < 0.002 &&  abs( tx-0.000 ) < 0.002 ) {
		  cmsdyfctHisto->fill( cmsdy*1E3 );
		  if( nrow <= 2 ) cmsdyfcntHisto->fill( cmsdy*1E3 );
		}

		if( abs( ty-0.000 ) < 0.002 &&
		    abs( tx-0.000 ) < 0.002 &&
		    Q0 > 18 ) {

		  cmsdyfctqHisto->fill( cmsdy*1E3 ); // 7.8 um @ 4 GeV, 19 deg
		  if( nrow <= 2 ) cmsdyfcntqHisto->fill( cmsdy*1E3 );

		  if( Q0 < 40 ) {
		    cmsdyfctq1Histo->fill( cmsdy*1E3 ); // 7.8 um @ 4 GeV, 19 deg, more Gaussian
		    if( nrow <= 2 ) cmsdyfcntq1Histo->fill( cmsdy*1E3 );
		    if( c->col < 26 )
		      cmsdyfctq1lHisto->fill( cmsdy*1E3 ); // xdb
		    else
		      cmsdyfctq1rHisto->fill( cmsdy*1E3 ); // xdb
		  }
		  if( Q0 < 35 ) {
		    cmsdyfctq2Histo->fill( cmsdy*1E3 ); // inserted 26.12.2012
		  }
		  if( Q0 < 30 ) {
		    cmsdyfctq3Histo->fill( cmsdy*1E3 ); // was fctq2. 7.4 um @ 4 GeV, 19 deg
		    if( ldot ) 
		      cmsdyfctqdotHisto->fill( cmsdy*1E3 ); // 8.1 um in run 5234
		    else
		      cmsdyfctq3dHisto->fill( cmsdy*1E3 ); // 7.2 um in run 5234
		  }
		}

	      } // CMS fiducial for y

	      // accumulate cuts for x:

	      if( fiducial && abs( cmsdy ) < 0.100 ) {

		if( abs( ty-0.000 ) < 0.002 &&  abs( tx-0.000 ) < 0.002 ){
		  cmsdxfctHisto->fill( cmsdx*1E3 );
		}

		if( abs( ty-0.000 ) < 0.002 &&
		    abs( tx-0.000 ) < 0.002 &&
		    Q0 > 18 ) {

		  cmsdxfctqHisto->fill( cmsdx*1E3 );

		  if( Q0 < 40 ) {
		    cmsdxfctq1Histo->fill( cmsdx*1E3 );
		  }

		  if( Q0 < 35 ) {
		    cmsdxfctq2Histo->fill( cmsdx*1E3 );
		  }

		  if( Q0 < 30 ) {
		    cmsdxfctq3Histo->fill( cmsdx*1E3 );
		  }

		}

	      } // CMS fiducial for x

	      // match CMS cluster and telescope track:
	      //FIXME this looks like cuts on the DUT hit position
	      //w.r.t. to the "track":
	      //printf("MONITOR: abs(cmsdx) = %f && abs(cmsdy) = %f\n",abs(cmsdx),abs(cmsdy));
	      if( abs( cmsdx ) < 0.3 && abs( cmsdy ) < 0.3 ){

		nmtch++;

		cmscolHisto->fill( c->col );
		cmsrowHisto->fill( c->row );

		cmsqHisto->fill( c->charge );
		cmsq0Histo->fill( Q0 );

		trixlkHisto->fill( -xA ); // -xA = x_DP = out
		triylkHisto->fill( -yA ); // -yA = y_DP = up
		trixylkHisto->fill( -xA, -yA );

		if( lq ) {
		  cmsdxvsx->fill( xAt, cmsdx*1E3 ); // dx vs x: turn?
		  cmsdxvsy->fill( yAt, cmsdx*1E3 ); // dx vs y: rot?
		  cmsdyvsx->fill( xAt, cmsdy*1E3 ); // dy vs x: rot?
		  cmsdyvsy->fill( yAt, cmsdy*1E3 ); // dy vs y: tilt?
		  cmsdxvstx->fill( tx*1E3, cmsdx*1E3 ); // dx vs tetx: z shift?
		  cmsdyvsty->fill( ty*1E3, cmsdy*1E3 ); // dy vs tety: z shift?
		  cmsdyvsxHisto->fill( xAt, cmsdy*1E3 ); // 
		}

		if( fiducial ) {

		  cmsnpxHisto->fill( c->size );
		  cmsncolHisto->fill( ncol );
		  cmsnrowHisto->fill( nrow );
		  cmsnrowvst1->fill( (event->getTimeStamp()-time0)/fTLU, nrow ); // cluster rows vs time

		  if( lq ) cmsnrowqHisto->fill( nrow ); // no 3-rows anymore

		  if( nrow == 2 ) cmsetaHisto->fill( eta );

		  cmsqfHisto->fill( c->charge );
		  cmsq0fHisto->fill( Q0 );

		  if( ldot ) { // sensor bias dot
		    cmsnpx0Histo->fill( c->size );
		    cmsqf0Histo->fill( c->charge );
		  }
		  else {
		    cmsnpx1Histo->fill( c->size );
		    cmsqf1Histo->fill( c->charge );
		    if( lcore ) { // pixel core
		      cmsnpx2Histo->fill( c->size );
		      cmsqf2Histo->fill( c->charge );
		    }
		    else {
		      cmsnpx3Histo->fill( c->size );
		      cmsqf3Histo->fill( c->charge );
		    }//core
		  }//dot

		  cmsqvsx->fill( xAt, c->charge ); // cluster charge profile
		  cmsqvsy->fill( yAt, c->charge ); // cluster charge profile
		  cmsqvsxm->fill( xmod, c->charge ); //q within pixel
		  cmsqvsym->fill( ymod, c->charge ); //q within pixel
		  cmsqvsxmym->fill( xmod, ymod, c->charge ); // cluster charge profile

		  cmsqvsddt->fill( dutddtns, c->charge ); // cluster charge vs phase
		  cmsqvst1->fill( (event->getTimeStamp()-time0)/fTLU, c->charge ); // cluster charge vs time
		  cmsqvst2->fill( (event->getTimeStamp()-time0)/fTLU, c->charge ); // cluster charge vs time
		  cmsqvst3->fill( (event->getTimeStamp()-time0)/fTLU, c->charge ); // cluster charge vs time
		  cmsqvst4->fill( (event->getTimeStamp()-time0)/fTLU, c->charge ); // cluster charge vs time

		  cmsrmsxvsq->fill( Q0, abs(cmsdx)*1E3 ); //resolution vs charge
		  cmsrmsyvsq->fill( Q0, abs(cmsdy)*1E3 ); //resolution vs charge

		  double pMoyal = exp( -exp( -( Q0 - 28 ) / 3.3 ) ); // fitMoyal
		  cmspMoyalvsq->fill( Q0, pMoyal );
		  cmspMoyalHisto->fill( pMoyal );
		  cmsrmsyvsp->fill( pMoyal, abs(cmsdy)*1E3 ); // resolution vs charge

		  if( lq ) {

		    cmsdyvsxm->fill( xmod, cmsdy*1E3 );
		    cmsdyvsym->fill( ymod, cmsdy*1E3 );

		    cmspixvsxmym->fill( xmod, ymod ); // occupancy map

		    cmsrmsxvsx->fill( xAt, abs(cmsdx)*1E3 ); //resolution across cols
		    cmsrmsyvsx->fill( xAt, abs(cmsdy)*1E3 ); //resolution across cols
		    cmsrmsxvsy->fill( yAt, abs(cmsdx)*1E3 ); //resolution across rows
		    cmsrmsyvsy->fill( yAt, abs(cmsdy)*1E3 ); //resolution across rows
		    cmsrmsxvsxm->fill( xmod, abs(cmsdx)*1E3 ); //resolution within pixel
		    cmsrmsyvsxm->fill( xmod, abs(cmsdy)*1E3 ); //resolution within pixel
		    cmsncolvsxm->fill( xmod, ncol );
		    cmsnrowvsxm->fill( xmod, nrow );
		    if( !ldot ) {
		      cmsrmsxvsym->fill( ymod, abs(cmsdx)*1E3 ); //resolution within pixel
		      cmsrmsyvsym->fill( ymod, abs(cmsdy)*1E3 ); //resolution within pixel
		      cmsrmsyvsym3->fill( ymd3, abs(cmsdy)*1E3 ); //resolution within pixel
		      cmsrmsyvsym6->fill( ymd6, abs(cmsdy)*1E3 ); //resolution within pixel
		    }
		    cmsrmsyvst->fill( (event->getTimeStamp()-time0)/fTLU, abs(cmsdy)*1E3 ); //resolution vs time
		    cmsrmsyvsddt->fill( dutddtns, abs(cmsdy)*1E3 ); // flat

		    if( nrow <= 2 ) {
		      cmsdyvseta->fill( eta, cmsdy*1E3 );
		      cmsrmsyvseta->fill( eta, abs(cmsdy)*1E3 );
		    }

		    cmsnpxvsxmym->fill( xmod, ymod, c->size ); // cluster size map
		    cmsncolvsym->fill( ymod, ncol ); // within pixel
		    cmsnrowvsym->fill( ymod, nrow ); // within pixel

		    cmsetavsym->fill( ymod, eta ); // eta within pixel
		    if( nrow == 2 ) cmsetavsym2->fill( ymod, eta ); // eta within pixel
		    cmsetavsym3->fill( ymd3, eta ); // eta within pixel

		    if(      nrow == 1 ) 
		      cmsym1Histo->fill( ymod ); //where are 1-rows?
		    else if( nrow == 2 ) 
		      cmsym2Histo->fill( ymod ); //where are 2-rows?

		  } // q Landau peak

		} // fiducial

	      } // CMS - triplet match

	      if( ltri && abs( cmsdx ) < 0.5 && abs( cmsdy ) < 0.5 ){ // link to CMS
		lkA[ntri-1] = 1; // ntri was already incremented
		DUTsxA[ntri-1] = cmsdx;
		DUTdyA[ntri-1] = cmsdy;
		ntrilk++;

	      }

	      // alignment:

	      if( fiducial && ltri && abs( cmsdx ) < 0.2 && abs( cmsdy ) < 0.15 ) {

		// GBL track fit for alignment
		vector<GblPoint> traj_points;
		// GblTrajectory traj( false ); // curvature = false

		// build up trajectory:

		vector<double> sPoint;

		// plane 0:

		double s = 0;

		TMatrixD proL2m(2,2);
		proL2m.UnitMatrix();

		TVectorD meas(2);

		double res = 3.5E-3; // [mm] Anemone telescope intrinsic resolution
		//res = 4.5E-3; // EUDET

		TVectorD measPrec(2); // precision = 1/resolution^2
		measPrec[0] = 1.0 / res / res;
		measPrec[1] = 1.0 / res / res;

		// scatter:
		TVectorD scat(2);
		scat.Zero(); //mean is zero

		double p = _eBeam; // beam momentum
		double X0Si = 65e-3 / 94; // Si + Kapton
		double tetSi = 0.0136 * sqrt(X0Si) / p * ( 1 + 0.038*std::log(X0Si) );

		TVectorD wscatSi(2);
		wscatSi[0] = 1.0 / ( tetSi * tetSi ); //weight
		wscatSi[1] = 1.0 / ( tetSi * tetSi );

		vector<unsigned int> ilab; // 0-2 = telescope, 3 = DUT

		// plane 0-2:

		double rx[3];
		double ry[3];
		double zprev = _planePosition[0];
		double step = 0;

		int iA = ntri-1; // ntri already incremented

		TMatrixD jacPointToPoint( 5, 5 );

		for( int ipl = 0; ipl < 3; ++ipl ){

		  int lhit = hts[ipl][iA];

		  step = _planePosition[ipl] - zprev;
		  zprev = _planePosition[ipl];

		  jacPointToPoint = Jac5( step );

		  GblPoint * point = new GblPoint( jacPointToPoint );
		  s += step;

		  double dz = hitZ[lhit] - zmA[iA];
		  double xs = xmA[iA] + sxA[iA] * dz; // Ax at plane
		  double ys = ymA[iA] + syA[iA] * dz; // Ay at plane

		  rx[ipl] = hitX[lhit] - xs;
		  ry[ipl] = hitY[lhit] - ys;

		  meas[0] = rx[ipl];
		  meas[1] = ry[ipl];

		  point->addMeasurement( proL2m, meas, measPrec );

		  point->addScatterer( scat, wscatSi );

		  traj_points.push_back(*point);
		  sPoint.push_back( s );

		  delete point;

		} // loop over planes

		// add DUT:

		step = DUTz - zprev;

		jacPointToPoint = Jac5( step );
		GblPoint * point = new GblPoint( jacPointToPoint );
		s += step;

		double tetDUT = 0.0136 * sqrt(DUTX0) / p * ( 1 + 0.038*std::log(DUTX0) );

		TVectorD wscatDUT(2);
		wscatDUT[0] = 1.0 / ( tetDUT * tetDUT ); //weight
		wscatDUT[1] = 1.0 / ( tetDUT * tetDUT );

		point->addScatterer( scat, wscatDUT );

		// DUT measurement:

		naldut++;

		meas[0] = DUTsxA[iA]; // cmsdx
		meas[1] = DUTdyA[iA]; // cmsdy

		double resx = 48E-3; // [mm] CMS col resolution
		if( abs( turn ) > 11 ) resx = 22E-3;
		double resy = 12E-3; // [mm] CMS row resolution at 20 deg tilt
		resy = 22E-3; // same weight as x

		TVectorD measWeight(2);
		measWeight[0] = 1.0 / resx / resx; // weight = 1/resolution^2
		measWeight[1] = 1.0 / resy / resy;

		// rot*tilt*turn:

		proL2m[0][0] = cf*co+sf*sa*so;
		proL2m[0][1] = sf*ca;
		proL2m[1][0] =-sf*co+cf*sa*so;
		proL2m[1][1] = cf*ca;

		// inverse:

		point->addMeasurement( proL2m, meas, measWeight );

		// DUT alignment: dx, dy, drot, dtilt, dturn, (dz)

		int nAlignmentParametersForDUT = 6;

		TMatrixD DUTDer( 2, nAlignmentParametersForDUT ); // alignment derivatives
		// residuals in DUT plane
		// transformed track impact point coordinates
		// residx = col - x4;
		// residy = row - y4;
		// alignment = linearized correction around initial values
		// may need to iterate

		DUTDer[0][0] =-1.0; // dresidx/ddeltax
		DUTDer[1][0] = 0.0;

		DUTDer[0][1] = 0.0;
		DUTDer[1][1] =-1.0; // dresidy/ddeltay

		DUTDer[0][2] =-upsign*y2; // dresidx/drot, linearized
		DUTDer[1][2] =-upsign*x2; // dresidy/drot

		DUTDer[0][3] =-upsign*( cf*(-so*dzcda)+sf*(-sa*y1 + ca*z1 + sa*co*dzcda)); // dresidx/dtilt a
		DUTDer[1][3] = upsign*(-sf*(-so*dzcda)+cf*(-sa*y1 + ca*z1 + sa*co*dzcda)); // dresidy/dtilt a

		derxtiltHisto->fill( DUTDer[0][3] );
		derytiltHisto->fill( DUTDer[1][3] );

		DUTDer[0][4] =-upsign*( cf*(-so*xc - co*dzc - so*dzcdo) + sf*sa*(co*xc-so*dzc+co*dzcdo)); // dresidx/dturn o
		DUTDer[1][4] = upsign*(-sf*(-so*xc - co*dzc - so*dzcdo) + cf*sa*(co*xc-so*dzc+co*dzcdo)); // dresidy/dturn o

		derxturnHisto->fill( DUTDer[0][4] );
		deryturnHisto->fill( DUTDer[1][4] );

		// dz:

		DUTDer[0][5] =-upsign*( cf*(co*dxcdz-so*ddzcdz) + sf*(ca*dycdz+sa*(so*dxcdz+co*ddzcdz))); // dresidx/dz
		DUTDer[1][5] = upsign*(-sf*(co*dxcdz-so*ddzcdz) + cf*(ca*dycdz+sa*(so*dxcdz+co*ddzcdz))); // dresidy/dz

		// global labels for Pede:

		std::vector<int> DUTLabels( nAlignmentParametersForDUT );
		DUTLabels[0] = 1; // dx
		DUTLabels[1] = 2; // dy
		DUTLabels[2] = 3; // drot
		DUTLabels[3] = 4; // dtilt
		DUTLabels[4] = 5; // dturn
		DUTLabels[5] = 6; // dz

		point->addGlobals( DUTLabels, DUTDer ); // for MillePede alignment
		traj_points.push_back(*point);
		//unsigned int iLabel = traj.addPoint(*point);

		// ilab[3] = iLabel;

		delete point;

		double Chi2;
		int Ndf;
		double lostWeight;

		GblTrajectory traj(traj_points, false ); // curvature = false
		traj.getLabels(ilab);

		traj.fit( Chi2, Ndf, lostWeight );

		traj.milleOut( *mille ); // write to mille.bin

	      } // linked triplet

	    } // loop over CMS clusters

	  }// have some CMS hit

	} // have telescope triplet

      }//loop over plane 1 hits k

    }//loop over plane 2 hits j

  }// loop over plane 0 hits i


  ntriHisto->fill( ntri );
  lkAvst->fill( (event->getTimeStamp()-time0)/fTLU, ntrilk );


  // ##########################################################################


  //----------------------------------------------------------------------------
  // telescope triplet efficiency w.r.t. CMS pixel:

  if( NpixDUT > 0 ){

    // CMS pixel clusters:

    for( vector<cluster>::iterator c = ClustDUT.begin(); c != ClustDUT.end(); c++ ){

      double cmsx = ( c->col - 26 ) * pitchcol; // -3.9..3.9 mm
      double cmsy = ( c->row - 40 ) * pitchrowt; // -4..4 mm

      double cmsxe = cmsx + DUTrot*cmsy - DUTalignx; // aligned w.r.t. telescope
      double cmsye = cmsy - DUTrot*cmsx - DUTaligny; // 

      cmsxeHisto->fill( cmsxe );
      cmsyeHisto->fill( cmsye );

      // triplets:

      int  nm = 0;
      bool im = 0;

      for( int i = 0; i < ntri && i < 99; ++i ) {

	double xA = xmA[i] + sxA[i] * ( DUTz - zmA[i] ); // A at DUT
	double yA = ymA[i] + syA[i] * ( DUTz - zmA[i] );

	double sx = cmsxe + xA;
	double dy = cmsye - yA;

	if( _dutReader.getChip() >= 100 ) { // xdb and dig mounted upright
	  sx = cmsxe - xA;
	  dy = cmsye + yA;
	}

	cmsdxeHisto->fill( sx*1E3 );
	cmsdyeHisto->fill( dy*1E3 );

	if( abs(sx) < 0.300 && abs(dy) < 0.200 ) {
	  nm++;
	  im = 1;
	}

      }//triplets

      cmsnmHisto->fill( nm );
      trieffvsxy->fill( cmsxe, cmsye, im ); //efficiency profile

    }//clusters

  }//have CMS data

  //----------------------------------------------------------------------------
  // six: triplets A and driplets B
  // matching and GBL fit
  // kinks: triplets A vs driplets B
  // scattering point = DUT:

  double sixCut = 0.1; // [mm]

  int nsix = 0;

  for( int iA = 0; iA < ntri && iA < 99; ++iA ) { // i = A = upstream

    double avx = xmA[iA];
    double avy = ymA[iA];
    double avz = zmA[iA];
    double tx = sxA[iA];
    double ty = syA[iA];

    double zA = DUTz - avz; // z CMS from mid of triplet
    double xA = avx + tx * zA; // triplet impact point on CMS
    double yA = avy + ty * zA;

    xA = xA / ( 1.0 - tx * tan(turn*wt) ); // turn correction
    yA = yA / ( 1.0 - syA[iA] * tan( tilt*wt ) ); // tilt correction

    // transform into DUT system:

    double xAm =-xA + DUTalignx; // invert x, shift to mid
    double yAm = yA + DUTaligny; // shift to mid
    if( _dutReader.getChip() >= 100 ) { // xdb and dig mounted upright
      xAm = xA + DUTalignx; // shift to mid
      yAm =-yA + DUTaligny; // invert y, shift to mid
    }
    double xAr = xAm - DUTrot*yAm; // rotate CMS pixel
    double yAr = yAm + DUTrot*xAm;
    double xAt = xAr / cos(turn*wt); // turn in x
    double yAt = yAr / cos(tilt*wt); // tilt in y

    // intersect inclined track with tilted DUT plane:

    double zc = (Nz*zA - Ny*avy - Nx*avx) / (Nx*tx + Ny*ty + Nz); // from avz
    double yc = avy + ty * zc;
    double xc = avx + tx * zc;
    double dzc = zc + avz - DUTz; // from DUT z0 [-8,8] mm

    // transform into DUT system: (passive).
    // large rotations don't commute: careful with order

    double x1 = co*xc - so*dzc; // turn
    double y1 = yc;
    double z1 = so*xc + co*dzc;

    double x2 = x1;
    double y2 = ca*y1 + sa*z1; // tilt
    double z2 =-sa*y1 + ca*z1; // should be zero (in DUT plane). is zero

    double x3 = cf*x2 + sf*y2; // rot
    double y3 =-sf*x2 + cf*y2;

    double upsign = -1; // analog psi46: upside down
    if( _dutReader.getChip() >= 100 )
      upsign = 1; // xdb and dig mounted upright
    double x4 = upsign*x3 + DUTalignx; // shift to mid
    double y4 =-upsign*y3 + DUTaligny; // invert y, shift to mid

    if( transform ) { // new style
      xAt = x4;
      yAt = y4;
    }

    double xmod = fmod( 9.075 + xAt, 0.3 ) * 1E3; // [0,300] um, 2 pixel wide
    double ymod = fmod( 9.050 + yAt, 0.2 ) * 1E3; // [0,200] um

    for( int jB = 0; jB < ndri && jB < 99; ++jB ) { // j = B = downstream

      double kx = sxB[jB] - sxA[iA]; //kink
      double ky = syB[jB] - syA[iA];
      double k2 = kx*kx + ky*ky;
      double xB = xmB[jB] + sxB[jB] * ( DUTz - zmB[jB] ); //B at DUTz
      double yB = ymB[jB] + syB[jB] * ( DUTz - zmB[jB] );
      xB = xB / ( 1.0 - sxB[jB] * tan( turn*wt ) ); // tilt correction
      yB = yB / ( 1.0 - syB[jB] * tan( tilt*wt ) ); // tilt correction

      double dx = xB - xA; // driplet - triplet
      double dy = yB - yA;

      double xR = xmB[jB] + sxB[jB] * ( REFz - zmB[jB] ); //B at REFz
      double yR = ymB[jB] + syB[jB] * ( REFz - zmB[jB] );

      sixkxHisto->fill( kx*1E3 );
      sixkyHisto->fill( ky*1E3 );
      sixdxHisto->fill( dx );
      sixdyHisto->fill( dy );
      if( abs(dy) < 0.4 ) sixdxcHisto->fill( dx*1E3 ); // sig = 17 um at 5 GeV
      if( abs(dx) < 0.4 ) sixdycHisto->fill( dy*1E3 );

      double probchi = 0;

      // match driplet and triplet:

      if( abs(dx) < sixCut && abs(dy) < sixCut ) {

	nsix++;

	sixkxcHisto->fill( kx*1E3 );
	sixkycHisto->fill( ky*1E3 );
	sixxHisto->fill( -xA ); // -xA = x_DP = out
	sixyHisto->fill( -yA ); // -yA = y_DP = up
	sixxyHisto->fill( -xA, -yA ); // DP: x_out, y_up

	// large kink map:

	if( abs( kx ) > 0.002 || abs( ky ) > 0.002 ) sixxycHisto->fill( -xA, -yA );

	kinkvsxy->fill( -xA, -yA, k2*1E6 ); //<kink^2> [mrad^2]

	// GBL with triplet A as seed:

	std::vector<GblPoint> traj_points;
	//GblTrajectory traj( false ); // curvature = false

	// build up trajectory:

	vector<double> sPoint;

	// plane 0:

	double s = 0;

	TMatrixD proL2m(2,2);
	proL2m.UnitMatrix();

	TVectorD meas(2);

	double res = 3.5E-3; // [mm] Anemone telescope intrinsic resolution
	//res = 4.5E-3; // EUDET

	TVectorD measPrec(2); // precision = 1/resolution^2
	measPrec[0] = 1.0 / res / res;
	measPrec[1] = 1.0 / res / res;

	// scatter:
	TVectorD scat(2);
	scat.Zero(); //mean is zero

	double p = _eBeam; // beam momentum
	double X0Si = 65e-3 / 94; // Si + Kapton
	double tetSi = 0.0136 * sqrt(X0Si) / p * ( 1 + 0.038*std::log(X0Si) );

	TVectorD wscatSi(2);
	wscatSi[0] = 1.0 / ( tetSi * tetSi ); //weight
	wscatSi[1] = 1.0 / ( tetSi * tetSi );

	TMatrixD alDer( 2, 3 ); // alignment derivatives
	alDer[0][0] = 1.0; // dx/dx
	alDer[1][0] = 0.0; // dy/dx
	alDer[0][1] = 0.0; // dx/dy
	alDer[1][1] = 1.0; // dy/dy

	std::vector<unsigned int> ilab; // 0-5 = telescope, 6 = DUT, 7 = REF

	// plane 0-5:

	double rx[6];
	double ry[6];
	double zprev = _planePosition[0];

	for( int ipl = 0; ipl < 6; ++ipl ){

	  int lhit;
	  if( ipl < 3 )
	    lhit = hts[ipl][iA];
	  else
	    lhit = hts[ipl][jB];

	  double step = _planePosition[ipl] - zprev;
	  zprev = _planePosition[ipl];

	  TMatrixD jacPointToPoint( 5, 5 );
	  jacPointToPoint = Jac5( step );

	  GblPoint * point = new GblPoint( jacPointToPoint );
	  s += step;

	  double dz = hitZ[lhit] - zmA[iA];
	  double xs = xmA[iA] + sxA[iA] * dz; // Ax at plane
	  double ys = ymA[iA] + syA[iA] * dz; // Ay at plane

	  rx[ipl] = hitX[lhit] - xs;
	  ry[ipl] = hitY[lhit] - ys;

	  if( ipl == 0 ) {
	    sixx0Histo->fill( -hitX[lhit] );
	    sixy0Histo->fill( -hitY[lhit] );
	  }
	  if( ipl == 1 ) {
	    sixx1Histo->fill( -hitX[lhit] );
	    sixy1Histo->fill( -hitY[lhit] );
	  }
	  if( ipl == 2 ) {
	    sixx2Histo->fill( -hitX[lhit] );
	    sixy2Histo->fill( -hitY[lhit] );
	  }
	  if( ipl == 3 ) {
	    sixx3Histo->fill( -hitX[lhit] );
	    sixy3Histo->fill( -hitY[lhit] );
	  }
	  if( ipl == 4 ) {
	    sixx4Histo->fill( -hitX[lhit] );
	    sixy4Histo->fill( -hitY[lhit] );
	  }
	  if( ipl == 5 ) {
	    sixx5Histo->fill( -hitX[lhit] );
	    sixy5Histo->fill( -hitY[lhit] );
	  }

	  meas[0] = rx[ipl];
	  meas[1] = ry[ipl];

	  point->addMeasurement( proL2m, meas, measPrec );

	  point->addScatterer( scat, wscatSi );

	  std::vector<int> globalLabels(3);
	  globalLabels[0] = 10 + ipl; // dx
	  globalLabels[1] = 20 + ipl; // dy
	  globalLabels[2] = 40 + ipl; // rot
	  alDer[0][2] = -ys; // dx/rot
	  alDer[1][2] =  xs; // dy/rot

	  //DP 2013 point->addGlobals( globalLabels, alDer ); // for MillePede alignment
	  traj_points.push_back(*point);
	  //unsigned int iLabel = traj.addPoint(*point);
	  //ilab[ipl] = iLabel;
	  sPoint.push_back( s );

	  delete point;

	  if( lDUT && ipl == 2 ) { // insert DUT

	    step = DUTz - zprev;
	    zprev = DUTz;

	    jacPointToPoint = Jac5( step );
	    point = new GblPoint( jacPointToPoint );
	    s += step;

	    double tetDUT = 0.0136 * sqrt(DUTX0) / p * ( 1 + 0.038*std::log(DUTX0) );

	    TVectorD wscatDUT(2);
	    wscatDUT[0] = 1.0 / ( tetDUT * tetDUT ); //weight
	    wscatDUT[1] = 1.0 / ( tetDUT * tetDUT );

	    point->addScatterer( scat, wscatDUT );

	    // DUT measurement:

	    if( lkA[iA] ) {

	      meas[0] = DUTsxA[iA]; // cmsdx
	      meas[1] = DUTdyA[iA]; // cmsdy

	      double resx = 48E-3; // [mm] CMS col resolution
	      if( abs( turn ) > 11 ) resx = 22E-3;
	      double resy =  8E-3; // [mm] CMS row resolution at 20 deg tilt

	      TVectorD measWeight(2);
	      measWeight[0] = 1.0 / resx / resx; // weight = 1/resolution^2
	      measWeight[1] = 1.0 / resy / resy;

	      point->addMeasurement( proL2m, meas, measWeight );

	    }//lkA linked hit in DUT

	    traj_points.push_back(*point);
	    //iLabel = traj.addPoint(*point);

	    //ilab[6] = iLabel;

	    delete point;

	  }//DUT present

	} // loop over planes

	// REF:

	// monitor what we put into GBL:

	selxHisto->fill( -xA ); // triplet at DUT
	selyHisto->fill( -yA );
	selaxHisto->fill( sxA[iA]*1E3 );
	selayHisto->fill( syA[iA]*1E3 );
	seldxHisto->fill( dx*1E3 ); // triplet-driplet match
	seldyHisto->fill( dy*1E3 );
	selkxHisto->fill( kx*1E3 ); // triplet-driplet kink
	selkyHisto->fill( ky*1E3 );

	seldx1Histo->fill( rx[1]*1E3 ); // triplet interpol
	seldy1Histo->fill( ry[1]*1E3 );
	seldx3Histo->fill( rx[3]*1E3 ); // triplet extrapol
	seldy3Histo->fill( ry[3]*1E3 );
	seldx4Histo->fill( rx[4]*1E3 );
	seldy4Histo->fill( ry[4]*1E3 );
	seldx5Histo->fill( rx[5]*1E3 );
	seldy5Histo->fill( ry[5]*1E3 );
	if( lkA[iA] ) {
	  seldx6Histo->fill( DUTsxA[iA]*1E3 );
	  seldy6Histo->fill( DUTdyA[iA]*1E3 );
	}

	double Chi2;
	int Ndf;
	double lostWeight;

	GblTrajectory traj(traj_points, false ); // curvature = false
	traj.getLabels(ilab);
	traj.fit( Chi2, Ndf, lostWeight );

	ngbl++;


	gblndfHisto->fill( Ndf );
	if( Ndf == 8 ) 
	  gblchi2aHisto->fill( Chi2 );
	else
	  gblchi2bHisto->fill( Chi2 );
	probchi = TMath::Prob( Chi2, Ndf );
	gblprbHisto->fill( probchi );

	// bad fits:

	if( probchi < 0.01 ) {

	  badxHisto->fill( -xA ); // triplet at DUT
	  badyHisto->fill( -yA );
	  badaxHisto->fill( sxA[iA]*1E3 );
	  badayHisto->fill( syA[iA]*1E3 );
	  baddxHisto->fill( dx*1E3 ); // triplet-driplet match
	  baddyHisto->fill( dy*1E3 );
	  badkxHisto->fill( kx*1E3 ); // triplet-driplet kink
	  badkyHisto->fill( ky*1E3 );

	  baddx1Histo->fill( rx[1]*1E3 ); // triplet interpol
	  baddy1Histo->fill( ry[1]*1E3 );
	  baddx3Histo->fill( rx[3]*1E3 ); // triplet extrapol
	  baddy3Histo->fill( ry[3]*1E3 );
	  baddx4Histo->fill( rx[4]*1E3 );
	  baddy4Histo->fill( ry[4]*1E3 );
	  baddx5Histo->fill( rx[5]*1E3 );
	  baddy5Histo->fill( ry[5]*1E3 );
	  if( lkA[iA] ) {
	    baddx6Histo->fill( DUTsxA[iA]*1E3 );
	    baddy6Histo->fill( DUTdyA[iA]*1E3 );
	  }

	}// bad fit

	else {

	  goodx1Histo->fill( rx[1]*1E3 ); // triplet interpol
	  goody1Histo->fill( ry[1]*1E3 );
	  if( lkA[iA] ) {
	    goodx6Histo->fill( DUTsxA[iA]*1E3 );
	    goody6Histo->fill( DUTdyA[iA]*1E3 );
	  }

	} // OK fit

	// look at fit:

	TVectorD aCorrection(5);
	TMatrixDSym aCovariance(5);

	double ax[8];
	double ay[8];
	unsigned int k = 0;

	// at plane 0:

	int ipos = ilab[0];
	traj.getResults( ipos, aCorrection, aCovariance );

	unsigned int ndim = 2;
	TVectorD aResiduals(ndim);
	TVectorD aMeasErrors(ndim);
	TVectorD aResErrors(ndim);
	TVectorD aDownWeights(ndim);

	traj.getMeasResults( static_cast<unsigned int>(ipos), ndim, aResiduals, aMeasErrors, aResErrors, aDownWeights );

	TVectorD aKinks(ndim);
	TVectorD aKinkErrors(ndim);
	TVectorD kResErrors(ndim);
	TVectorD kDownWeights(ndim);
	traj.getScatResults( static_cast<unsigned int>(ipos), ndim, aKinks, aKinkErrors, kResErrors, kDownWeights );

	//track = q/p, x', y', x, y
	//        0,   1,  2,  3, 4

	gblax0Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
	gbldx0Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
	gblrx0Histo->fill( ( rx[0] - aCorrection[3] ) * 1E3 ); // residual x [um]
	gblpx0Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
	gblqx0Histo->fill( aKinks[0]*1E3 ); // kink
	ax[k] = aCorrection[1]; // angle correction at plane, for kinks
	ay[k] = aCorrection[2]; // angle correction at plane, for kinks
	k++;

	ipos = ilab[1];
	traj.getResults( ipos, aCorrection, aCovariance );
	traj.getMeasResults(static_cast<unsigned int>(ipos), ndim, aResiduals, aMeasErrors, aResErrors, aDownWeights );
	traj.getScatResults(static_cast<unsigned int>(ipos), ndim, aKinks, aKinkErrors, kResErrors, kDownWeights );
	gblax1Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
	gbldx1Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
	gblrx1Histo->fill( ( rx[1] - aCorrection[3] ) * 1E3 ); // residual x [um]
	gblpx1Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
	gblqx1Histo->fill( aKinks[0]*1E3 ); // kink
	gblsx1Histo->fill( aKinks[0]/aKinkErrors[0] ); // x kink pull
	gbltx1Histo->fill( aKinks[0]/kResErrors[0] ); // x kink pull
	ax[k] = aCorrection[1]; // angle correction at plane, for kinks
	ay[k] = aCorrection[2]; // angle correction at plane, for kinks
	k++;

	ipos = ilab[2];
	traj.getResults( ipos, aCorrection, aCovariance );
	traj.getMeasResults(static_cast<unsigned int>(ipos), ndim, aResiduals, aMeasErrors, aResErrors, aDownWeights );
	traj.getScatResults(static_cast<unsigned int>(ipos), ndim, aKinks, aKinkErrors, kResErrors, kDownWeights );
	gblax2Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
	gbldx2Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
	gblrx2Histo->fill( ( rx[2] - aCorrection[3] ) * 1E3 ); // residual x [um]
	gblpx2Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
	gblqx2Histo->fill( aKinks[0]*1E3 ); // kink
	ax[k] = aCorrection[1]; // angle correction at plane, for kinks
	ay[k] = aCorrection[2]; // angle correction at plane, for kinks
	k++;

	if( lDUT ) {
	  ipos = ilab[6]; // 6 = DUT
	  traj.getResults( ipos, aCorrection, aCovariance );
	  traj.getScatResults(static_cast<unsigned int>(ipos), ndim, aKinks, aKinkErrors, kResErrors, kDownWeights );
	  gblax6Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
	  gbldx6Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
	  gbldy6Histo->fill( aCorrection[4]*1E3 ); // shift y [um]
	  gblqx6Histo->fill( aKinks[0]*1E3 ); // x kink
	  gblsx6Histo->fill( aKinks[0]/aKinkErrors[0] ); // x kink pull
	  gbltx6Histo->fill( aKinks[0]/kResErrors[0] ); // x kink pull
	  if( lkA[iA] ) {
	    traj.getMeasResults(static_cast<unsigned int>(ipos), ndim, aResiduals, aMeasErrors, aResErrors, aDownWeights );
	    gblrx6Histo->fill( ( DUTsxA[iA] - aCorrection[3] ) * 1E3 ); // residual x [um]
	    gblry6Histo->fill( ( DUTdyA[iA] - aCorrection[4] ) * 1E3 ); // residual y [um]
	    gblpx6Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
	    gblpy6Histo->fill( aResiduals[1] / aResErrors[1] ); // pull
	  }
	  ax[k] = aCorrection[1]; // angle correction at plane, for kinks
	  ay[k] = aCorrection[2]; // angle correction at plane, for kinks
	  k++;
	}//DUT

	ipos = ilab[3];
	traj.getResults( ipos, aCorrection, aCovariance );
	traj.getMeasResults(static_cast<unsigned int>(ipos), ndim, aResiduals, aMeasErrors, aResErrors, aDownWeights );
	traj.getScatResults(static_cast<unsigned int>(ipos), ndim, aKinks, aKinkErrors, kResErrors, kDownWeights );
	gblax3Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
	gbldx3Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
	gblrx3Histo->fill( ( rx[3] - aCorrection[3] ) * 1E3 ); // residual x [um]
	gblpx3Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
	gblqx3Histo->fill( aKinks[0]*1E3 ); // kink
	ax[k] = aCorrection[1]; // angle correction at plane, for kinks
	ay[k] = aCorrection[2]; // angle correction at plane, for kinks
	k++;

	ipos = ilab[4];
	traj.getResults( ipos, aCorrection, aCovariance );
	traj.getMeasResults(static_cast<unsigned int>(ipos), ndim, aResiduals, aMeasErrors, aResErrors, aDownWeights );
	traj.getScatResults(static_cast<unsigned int>(ipos), ndim, aKinks, aKinkErrors, kResErrors, kDownWeights );
	gblax4Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
	gbldx4Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
	gblrx4Histo->fill( ( rx[4] - aCorrection[3] ) * 1E3 ); // residual x [um]
	gblpx4Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
	gblqx4Histo->fill( aKinks[0]*1E3 ); // kink
	ax[k] = aCorrection[1]; // angle correction at plane, for kinks
	ay[k] = aCorrection[2]; // angle correction at plane, for kinks
	k++;

	ipos = ilab[5];
	traj.getResults( ipos, aCorrection, aCovariance );
	traj.getMeasResults(static_cast<unsigned int>(ipos), ndim, aResiduals, aMeasErrors, aResErrors, aDownWeights );
	traj.getScatResults(static_cast<unsigned int>(ipos), ndim, aKinks, aKinkErrors, kResErrors, kDownWeights );
	gblax5Histo->fill( aCorrection[1]*1E3 ); // angle x [mrad]
	gbldx5Histo->fill( aCorrection[3]*1E3 ); // shift x [um]
	gblrx5Histo->fill( ( rx[5] - aCorrection[3] ) * 1E3 ); // residual x [um]
	gblpx5Histo->fill( aResiduals[0] / aResErrors[0] ); // pull
	gblqx5Histo->fill( aKinks[0]*1E3 ); // kink
	ax[k] = aCorrection[1]; // angle correction at plane, for kinks
	ay[k] = aCorrection[2]; // angle correction at plane, for kinks
	k++;

	// kinks: 1,2 = tele, 3 = DUT, 4,5 = tele

	gblkx1Histo->fill( (ax[1] - ax[0])*1E3 ); // kink at 1 [mrad]
	gblkx2Histo->fill( (ax[2] - ax[1])*1E3 ); // kink at 2 [mrad]
	gblkx3Histo->fill( (ax[3] - ax[2])*1E3 ); // kink at 3 [mrad]
	gblkx4Histo->fill( (ax[4] - ax[3])*1E3 ); // kink at 4 [mrad]
	gblkx5Histo->fill( (ax[5] - ax[4])*1E3 ); // kink at 5 [mrad]
	gblkx6Histo->fill( (ax[6] - ax[5])*1E3 ); // kink at 6 [mrad]

      } // driplet-triplet match

      //------------------------------------------------------------------------
      // CMS DUT efficiency:

      if( abs(dx) < 0.1 && abs(dy) < 0.1 && lkB[jB] ) { // six with link

	sixxylkHisto->fill( -xA, -yA ); // six-tracks with REF link at CMS DUT

	bool nm = 0;
	if( lkA[iA] ) nm = 1;

	double fidxmax =  12;
	double fidxmin = -12;
	double fidymax =   5;
	double fidymin =  -6;

	if(  event->getRunNumber() < 188 ) { // Feb 2012
	  fidymin = -3;
	  fidymax =  3;
	  fidxmin = -3;
	  fidxmax =  4;
	}
	if(  event->getRunNumber() > 2842 ) { // Apr 2012 TB 22
	  fidymin = -2.7;
	  fidymax =  3.5;
	  fidxmin = -3.8;
	  fidxmax =  2.5;
	}
	if(  event->getRunNumber() > 2877 ) { // Apr 2012 TB 22
	  fidymin = -3.8;
	  fidymax =  3.8;
	  fidxmin = -4.0;
	  fidxmax =  3.0;
	}
	if(  event->getRunNumber() >= 3048 ) { // Apr 2012 TB 21
	  fidymin = -3.8;
	  fidymax =  3.8;
	  fidxmin = -3.8;
	  fidxmax =  3.7;
	}
	if(  event->getRunNumber() >= 3702 ) { // Jul 2012 xdb TB 21
	  fidymin = -3.8;
	  fidymax =  3.8;
	  fidxmin = -3.8;
	  fidxmax =  3.8;
	}

	if( abs(turn) > 2 ) {
	  fidxmin = -3.0 - ( yAt + 4 ) / 10;
	  fidxmax =  3.8 - ( yAt + 4 ) / 10;
	}

	bool leff = 1;

	//FIXME parameter csv
	//use _leff_val: if((event->getTimeStamp()-time0)/fTLU > _leff_val) leff = 0;
	if(  event->getRunNumber() == 5511 && (event->getTimeStamp()-time0)/fTLU > 305 ) leff = 0;

	if( leff ) effxyHisto->fill( xAt, yAt );
	if( leff ) effvsxy->fill( xAt, yAt, nm ); // CMS DUT efficiency profile

	if( leff && yAt > fidymin && yAt < fidymax ) {
	  effvsx->fill( xAt, nm ); // CMS DUT efficiency profile
	  if( probchi > 0.01 ) effvsxg->fill( xAt, nm );
	}

	if( leff && xAt > fidxmin && xAt < fidxmax ) {
	  effvsy->fill( yAt, nm ); // CMS DUT efficiency profile
	}

	if( xAt > fidxmin && xAt < fidxmax && yAt > fidymin && yAt < fidymax ) {
	  eff300->fill( (event->getTimeStamp()-time0)/fTLU, nm );
	  eff600->fill( (event->getTimeStamp()-time0)/fTLU, nm );
	  eff1200->fill( (event->getTimeStamp()-time0)/fTLU, nm );
	  eff1800->fill( (event->getTimeStamp()-time0)/fTLU, nm );
	  eff3600->fill( (event->getTimeStamp()-time0)/fTLU, nm );
	  if( leff ) effvsddt->fill( dutddtns, nm );
	  if( leff ) effvsxmym->fill( xmod, ymod, nm ); // CMS DUT efficiency profile
	  if( prevdutrefddt == 0 && dutrefddt == 0 )
	    effd600->fill( (event->getTimeStamp()-time0)/fTLU, nm );
	  if( ndri < 3 ) effn600->fill( (event->getTimeStamp()-time0)/fTLU, nm );
	  if( ndri < 2 ) effm600->fill( (event->getTimeStamp()-time0)/fTLU, nm );

	}

      }// triplet-driplet match

      //------------------------------------------------------------------------
      // eff(REF) with DUT as timing plane:

      if( abs(dx) < 0.1 && abs(dy) < 0.1 && lkA[iA] ) {

	bool nm = 0;
	if( lkB[jB] ) nm = 1;

	rffvsxy->fill( -xR, -yR, nm ); // CMS REF efficiency profile

	if( abs( yR ) < 3 ) {
	  rffvsx->fill( -xR, nm ); // CMS REF efficiency profile
	}

      }// triplet-driplet match

      //------------------------------------------------------------------------
      // intersect point in z:

      if( abs(dy) < 0.1 ) { // no cut on dx
	if( abs( kx ) > 0.003 ) { // 
	  double zx = ( xmA[iA] - sxA[iA]*zmA[iA] - xmB[jB] + sxB[jB]*zmB[jB] ) / kx; //z intersect
	  sixzx3Histo->fill( zx - _planePosition[2] );
	}
	if( abs( kx ) > 0.002 ) { // 
	  double zx = ( xmA[iA] - sxA[iA]*zmA[iA] - xmB[jB] + sxB[jB]*zmB[jB] ) / kx; //z intersect
	  sixzx2Histo->fill( zx - _planePosition[2] );
	}
      }

      if( abs(dx) < 0.1 ) { // no cut on dy
	if( abs( ky ) > 0.003 ) { // 
	  double zy = ( ymA[iA] - syA[iA]*zmA[iA] - ymB[jB] + syB[jB]*zmB[jB] ) / ky;
	  sixzy3Histo->fill( zy - _planePosition[2] );
	}
	if( abs( ky ) > 0.002 ) { // 
	  double zy = ( ymA[iA] - syA[iA]*zmA[iA] - ymB[jB] + syB[jB]*zmB[jB] ) / ky;
	  sixzy2Histo->fill( zy - _planePosition[2] );
	}
      }

      //------------------------------------------------------------------------
      // z intersect:

      if( abs(dx) < 0.2 && abs(dy) < 0.2 ) { // looser cut allows more z range

	if( abs( kx ) > 0.001 ) {
	  double zx = ( xmA[iA] - sxA[iA]*zmA[iA] - xmB[jB] + sxB[jB]*zmB[jB] ) / kx; //z intersect
	  sixzx1Histo->fill( zx - _planePosition[2] );
	}

	if( abs( ky ) > 0.001 ) {
	  double zy = ( ymA[iA] - syA[iA]*zmA[iA] - ymB[jB] + syB[jB]*zmB[jB] ) / ky;
	  sixzy1Histo->fill( zy - _planePosition[2] );
	}

	// measure scattering angle x after cuts in y:
	// cut on ky creates bias in kx

	if( abs( ky ) > 0.001 ) {
	  double zy = ( ymA[iA] - syA[iA]*zmA[iA] - ymB[jB] + syB[jB]*zmB[jB] ) / ky;
	  if( abs( zy - DUTz ) < 30 ) {
	    sixkxzyHisto->fill( kx*1E3 );
	  }
	}

	if( abs( kx ) > 0.001 ) {
	  double zx = ( xmA[iA] - sxA[iA]*zmA[iA] - xmB[jB] + sxB[jB]*zmB[jB] ) / kx;
	  if( abs( zx - DUTz ) < 30 ) {
	    sixkyzxHisto->fill( ky*1E3 );
	  }
	}

	// kx at DUT:

	if( abs( kx ) > 0.001 ) {
	  double zx = ( xmA[iA] - sxA[iA]*zmA[iA] - xmB[jB] + sxB[jB]*zmB[jB] ) / kx;
	  if( abs( zx - DUTz ) < 30 ) {
	    sixkxzxHisto->fill( kx*1E3 ); // plot with gap, fittp0g.C("sixkxzx")
	  }
	}

	// ky at DUT:

	if( abs( ky ) > 0.001 ) {
	  double zy = ( ymA[iA] - syA[iA]*zmA[iA] - ymB[jB] + syB[jB]*zmB[jB] ) / ky;
	  if( abs( zy - DUTz ) < 30 ) {
	    sixkyzyHisto->fill( ky*1E3 ); // plot with gap
	  }
	}

      }//match

      //------------------------------------------------------------------------

    }//driplets B j

  }//triplets A i

  nsixHisto->fill( nsix );

  prevdutrefddt = dutrefddt;

  // Clear all working arrays

  delete [] planeHitID;
  delete [] hitPlane;
  delete [] hitZ;
  delete [] hitEy;
  delete [] hitY;
  delete [] hitEx;
  delete [] hitX;

  return;

}


//------------------------------------------------------------------------------
void EUTelAnalysisCMSPixel::check( LCEvent * /* evt */  ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


//------------------------------------------------------------------------------
void EUTelAnalysisCMSPixel::end(){

  // Print the summary:

  streamlog_out( MESSAGE ) << endl
			   << "Processed events:    "
			   << setw(10) << setiosflags(ios::right)
			   << _nEvt << resetiosflags(ios::right) << endl
			   << "DUT re-syncs " << dutnsync << endl
			   << "DUT clusters " << nclst << endl
			   << "DUT matches  " << nmtch << endl
			   << "REF re-syncs " << refnsync << endl
                           << endl;

  cout << "Histo: " << cmssyaHisto->title() << endl;
  cout << "Entries " << cmssyaHisto->entries() << endl;
  cout << "Maximum " << cmssyaHisto->maxBinHeight() << endl;

  // Clean memory:

  delete [] _planeSort;
  delete [] _planeID;
  delete [] _planePosition;
  delete [] _planeThickness;
  delete [] _planeX0;
  delete [] _planeResolution;

  // Pede:

  // close the output file:
  delete mille;

  streamlog_out( MESSAGE4 ) << endl << "Generating the steering file for the pede program..." << endl;

  streamlog_out( MESSAGE2 ) << "have " << ngbl << " GBL tracks" << endl;
  streamlog_out( MESSAGE2 ) << "have " << naldut << " DUT links" << endl;

  bool ldut = naldut > 99;

  ofstream steerFile;
  steerFile.open( "steerPede.txt" );

  if( steerFile.is_open() && ldut ) {

    steerFile << "! generated by EUTelTestFitter" << endl;
    steerFile << "Cfiles" << endl;
    steerFile << "mille.bin" << endl;
    steerFile << endl;

    steerFile << "Parameter" << endl;
    steerFile << 1 << "  0.0  0.0" << endl; // dx
    steerFile << 2 << "  0.0  0.0" << endl; // dy
    steerFile << 3 << "  0.0  0.0" << endl; // rot
    steerFile << 4 << "  0.0  0.0" << endl; // dtilt
    steerFile << 5 << "  0.0  0.0" << endl; // dturn
    steerFile << 6 << "  0.0  0.0" << endl; // dz
    steerFile << endl;
    steerFile << "! chiscut 5.0 2.5" << endl;
    steerFile << "outlierdownweighting 4" << endl;
    steerFile << endl;
    steerFile << "method inversion 10  0.1" << endl;
    steerFile << endl;
    steerFile << "! histprint" << endl;
    steerFile << endl;
    steerFile << "end" << endl;

    steerFile.close();

    streamlog_out( MESSAGE2 ) << "Pede steer file written." << endl;

    // before starting pede, let's check if it is in the path
    bool isPedeInPath = true;

    // create a new process
    redi::ipstream which("which pede");

    // wait for the process to finish
    which.close();

    // get the status
    // if it is 255 then the program wasn't found in the path
    isPedeInPath = !( which.rdbuf()->status() == 255 );

    //isPedeInPath = 0; // DP initially no pede

    if( !isPedeInPath ) {
      streamlog_out( ERROR ) << "Pede cannot be executed because not found in the path" << endl;
    }
    else {

      std::string command = "pede  steerPede.txt";

      streamlog_out( MESSAGE2 ) << endl;
      streamlog_out( MESSAGE2 ) << "Starting pede..." << endl;
      streamlog_out( MESSAGE2 ) << command.c_str() << endl;

      redi::ipstream pede( command.c_str() );
      string output;
      while ( getline( pede, output ) ) {
	streamlog_out( MESSAGE2 ) << output << endl;
      }

      // wait for the pede execution to finish
      pede.close();

      // check the exit value of pede
      if( pede.rdbuf()->status() == 0 )
	streamlog_out( MESSAGE2 ) << "Pede successfully finished" << endl;

      // reading back the millepede.res file:

      string millepedeResFileName = "millepede.res";

      streamlog_out( MESSAGE2 ) << "Reading back the " << millepedeResFileName << endl;

      // open the millepede ASCII output file
      ifstream millepede( millepedeResFileName.c_str() );

      if( millepede.bad() || !millepede.is_open() ) {
	streamlog_out( ERROR4 )
	  << "Error opening the " << millepedeResFileName << endl
	  << "The alignment file cannot be saved" << endl;
      }
      else {
	vector<double > tokens;
	stringstream tokenizer;
	string line;
	double buffer;

	// get the first line and throw it away since it is a comment!

	getline( millepede, line );
	std::cout << "line: " <<  line  << std::endl;

	unsigned int numpars = 6; // DUT dx, dy, drot, dtilt, dturn, dz

	map< unsigned int, double > alpar; // map = associative array

	for( unsigned int ipar = 0; ipar < numpars; ++ipar ) {

	  if( millepede.eof() ) break;

	  getline( millepede, line );

	  if( line.empty() ) continue;

	  tokens.clear();
	  tokenizer.clear();
	  tokenizer.str( line );

	  while( tokenizer >> buffer ) {
	    tokens.push_back( buffer );
	  }

	  int lpar = (int) ( tokens[0] + 0.5 ); // par label

	  bool isFixed = ( tokens.size() == 3 );

	  if( isFixed )
	    cout << "Parameter " << lpar
		 << " is at " << tokens[1]
		 << " (fixed)"  << endl;
	  else
	    cout << "Parameter " << lpar
		 << " is at " << tokens[1]
		 << " +/- " << tokens[4] << endl;

	  alpar[lpar] = tokens[1];

	}//loop param

	if( ldut ){
	  cout << "DUT alignment corrections:" << endl;
	  cout << "dx    " << alpar[1]*1E3 << " um" << endl;
	  cout << "dy    " << alpar[2]*1E3 << " um" << endl;
	  cout << "drot  " << alpar[3]*1E3 << " mrad" << endl;
	  cout << "dtilt " << alpar[4]*180/3.141592654 << " deg" << endl;
	  cout << "dturn " << alpar[5]*180/3.141592654 << " deg" << endl;
	  cout << "dz    " << alpar[6]*1E3 << " um" << endl;
	  cout << "DUT alignment:" << endl;
	  cout << "   DUTalignx = " << DUTalignx-alpar[1] << ";" << endl;
	  cout << "   DUTaligny = " << DUTaligny-alpar[2] << ";" << endl;
	  cout << "   DUTrot = " << DUTrot-alpar[3] << ";" << endl;
	  cout << "   tilt = " << tilt-alpar[4]*180/3.141592654 << ";" << endl;
	  cout << "   turn = " << turn-alpar[5]*180/3.141592654 << ";" << endl;
	  cout << "   DUTz = " << DUTz-alpar[6] - _planePosition[2] << " + _planePosition[2];" << endl;
	  cout << "Copy this to your runlist.csv file in order to update the DUT alignment parameters:" << endl
	  // RunNumber,Date,GearFile,BeamEnergy,dut_run,dut_chip,dut_board,dut_gain,ref_run,ref_chip,ref_board,ref_gain,dut_align_x,dut_align_y,dut_pos_z,dut_tilt,dut_turn,dut_rot,ref_align_x,ref_align_y,ref_pos_z,ref_rot,nskipdut,nskipref,nskiptel
	  << _TEL_run << "," << _DATE_run << "," << _gearfile << "," <<_eBeam << "," << _DUT_run << "," << _DUT_chip << "," << _DUT_board << "," << _DUT_gain << "," << _REF_run << "," << _REF_chip << "," << _REF_board << "," << _REF_gain << "," << DUTalignx-alpar[1] << "," << DUTaligny-alpar[2] << "," << DUTz-alpar[6] - _planePosition[2] << "," << tilt-alpar[4]*180/3.141592654 << "," << turn-alpar[5]*180/3.141592654 << "," << DUTrot-alpar[3] << "," << _REFalignx << "," << _REFaligny << "," << _REFz << "," << _REFrot << "," << _nSkip << "," << _nSkipRef << "," << _nSkipTelescope << endl;
	}
      }//millepede OK

      millepede.close();

    }//PedeInPath

  }// pede steer file open

  else {
    streamlog_out( ERROR2 ) << "Could not open steering file." << endl;
  }

} // end end


//------------------------------------------------------------------------------
void EUTelAnalysisCMSPixel::bookHistos()
{

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)


  // event time:
  AIDAProcessor::tree(this)->mkdir("Timing");

  t1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/t1", 125, 0, 1 );
  t1Histo->setTitle( "telescope event time;telescope event time [s];events/s" );

  t10Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/t10", 250, 0, 10 );
  t10Histo->setTitle( "telescope event time;telescope event time [s];events/s" );

  t100Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/t100", 100, 0, 100 );
  t100Histo->setTitle( "telescope event time;telescope event time [s];events/s" );

  t300Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/t300", 350, 0, 350 );
  t300Histo->setTitle( "telescope event time;telescope event time [s];events/s" );

  t600Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/t600", 700, 0, 700 );
  t600Histo->setTitle( "telescope event time;telescope event time [s];events/s" );

  t1000Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/t1000", 100, 0, 1000 );
  t1000Histo->setTitle( "telescope event time;telescope event time [s];events/10s" );

  t1800Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/t1800", 200, 0, 2000 );
  t1800Histo->setTitle( "telescope event time;telescope event time [s];events/10s" );

  t3600Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/t3600", 400, 0, 4000 );
  t3600Histo->setTitle( "telescope event time;telescope event time [s];events/10s" );

  dtHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/dt", 500, 0, 5000 );
  dtHisto->setTitle( "telescope time between events;telescope time between events [#mus];events" );

  dtmsHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/dtms", 500, 0, 100 );
  dtmsHisto->setTitle( "telescope time between events;telescope time between events [ms];events" );

  logdtHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/logdt", 500, -1, 4 );
  logdtHisto->setTitle( "telescope time between events;telescope time between events log_{10}(#Deltat [ms]);events" );

  logdtcmsHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/logdtcms", 500, -1, 4 );
  logdtcmsHisto->setTitle( "telescope time between events with DUT;telescope time between events log_{10}(#Deltat [ms]);events" );

  dtfvstau = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Timing/dtfvstau", 600, 372, 378, 0, 1 );
  dtfvstau->setTitle( "dt/tau fractional part;tau [ticks];dt/tau [fractional part]" );

  tfvstau = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Timing/tfvstau", 600, 372, 378, 0, 1 );
  tfvstau->setTitle( "t/tau fractional part;tau [ticks];t/tau [fractional part]" );

  dtfHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/dtf", 500, -0.5, 0.5 );
  dtfHisto->setTitle( "phase of telescope time between events;phase of time between events [turns];telescope events" );

  dtfvst = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "dtfvst", 300, 0, 300, -1, 1 );
  dtfvst->setTitle( "dt phase vs t;run time [s];<TLU phase> [turns];" );

  dtfvsdt = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Timing/dtfvsdt", 100, 0, 10000, -1, 1 );
  dtfvsdt->setTitle( "dt phase vs dt;delta time [us];<TLU phase> [turns];" );

  tfHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/tf", 100, -0.5, 0.5 );
  tfHisto->setTitle( "phase of telescope event time;phase of event time [turns];telescope events" );

  tfvst = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Timing/tfvst", 300, 0, 300E6, -1, 1 );
  tfvst->setTitle( "t phase vs t;turns;<TLU phase> [turns];" );

  tfvst1 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Timing/tfvst1", 100, 0, 1E6, -1, 1 );
  tfvst1->setTitle( "t phase vs t;turns;<TLU phase> [turns];" );

  tfvst10 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Timing/tfvst10", 1000, 0, 10E6, -1, 1 );
  tfvst10->setTitle( "t phase vs t;turns;<TLU phase> [turns];" );

  tfvst100 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Timing/tfvst100", 1000, 0, 100E6, -1, 1 );
  tfvst100->setTitle( "t phase vs t;turns;<TLU phase> [turns];" );

  tfvst300 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Timing/tfvst300", 1000, 0, 300E6, -1, 1 );
  tfvst300->setTitle( "t phase vs t;turns;<TLU phase> [turns];" );



  // DUT clusters:

  dutnclusHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "dutnclus", 11, -0.5, 11.5 );
  dutnclusHisto->setTitle( "DUT clusters/event;DUT clusters/event;events" );

  dutcolHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "dutcol", 52, -0.5, 51.5 );
  dutcolHisto->setTitle( "DUT column;DUT cluster col;DUT clusters" );

  dutrowHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "dutrow", 80, -0.5, 79.5 );
  dutrowHisto->setTitle( "DUT row;DUT cluster row;DUT clusters" );

  dutnpxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "dutnpx", 11, -0.5, 10.5 );
  dutnpxHisto->setTitle( "DUT cluster size;DUT pixel per cluster;DUT clusters" );

  dutadcHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "dutadc", 100, 0, 100 );
  dutadcHisto->setTitle( "DUT cluster charge;DUT cluster charge [ke];DUT clusters" );


  refnclusHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "refnclus", 11, -0.5, 11.5 );
  refnclusHisto->setTitle( "REF clusters/event;REF clusters/event;events" );

  refcolHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "refcol", 52, -0.5, 51.5 );
  refcolHisto->setTitle( "REF column;REF cluster col;REF clusters" );

  refrowHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "refrow", 80, -0.5, 79.5 );
  refrowHisto->setTitle( "REF row;REF cluster row;REF clusters" );

  refnpxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "refnpx", 11, -0.5, 10.5 );
  refnpxHisto->setTitle( "REF cluster size;REF pixel per cluster;REF clusters" );

  refadcHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "refadc", 100, 0, 100 );
  refadcHisto->setTitle( "REF cluster charge;REF cluster charge [ke];REF clusters" );

  // telescope hits per plane:
  AIDAProcessor::tree(this)->mkdir("Telescope");

  nAllHitHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/nallhit", 201, -0.5, 200.5 );
  nAllHitHisto->setTitle( "Telescope hits/event;telescope hits;events" );

  hits0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/hits0", 101, -0.5, 100.5 );
  hits0Histo->setTitle( "hits in plane 0;hits in plane 0;events" );

  hits1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/hits1", 101, -0.5, 100.5 );
  hits1Histo->setTitle( "hits in plane 1;hits in plane 1;events" );

  hits2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/hits2", 101, -0.5, 100.5 );
  hits2Histo->setTitle( "hits in plane 2;hits in plane 2;events" );

  hits3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/hits3", 101, -0.5, 100.5 );
  hits3Histo->setTitle( "hits in plane 3;hits in plane 3;events" );

  hits4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/hits4", 101, -0.5, 100.5 );
  hits4Histo->setTitle( "hits in plane 4;hits in plane 4;events" );

  hits5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/hits5", 101, -0.5, 100.5 );
  hits5Histo->setTitle( "hits in plane 5;hits in plane 5;events" );

  // telescope dx:

  dx01Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx01", 100, -1, 1 );
  dx01Histo->setTitle( "x1-x0;x_{1}-x_{0} [mm];hit pairs" );

  dy01Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dy01", 100, -1, 1 );
  dy01Histo->setTitle( "y1-y0;y_{1}-y_{0} [mm];hit pairs" );

  du01Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/du01", 100, -1, 1 );
  du01Histo->setTitle( "x1-x0, |dy| < 1;x_{1}-x_{0} [mm];hit pairs" );

  dx02Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx02", 100, -1, 1 );
  dx02Histo->setTitle( "x2-x0;x_{2}-x_{0} [mm];hit pairs" );

  dx03Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx03", 100, -1, 1 );
  dx03Histo->setTitle( "x3-x0;x_{3}-x_{0} [mm];hit pairs" );

  dx04Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx04", 100, -1, 1 );
  dx04Histo->setTitle( "x4-x0;x_{4}-x_{0} [mm];hit pairs" );

  dx05Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx05", 100, -1, 1 );
  dx05Histo->setTitle( "x5-x0;x_{5}-x_{0} [mm];hit pairs" );

  dx12Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx12", 100, -1, 1 );
  dx12Histo->setTitle( "x2-x1;x_{2}-x_{1} [mm];hit pairs" );

  dy12Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dy12", 100, -1, 1 );
  dy12Histo->setTitle( "y2-y1;y_{2}-y_{1} [mm];hit pairs" );

  du12Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/du12", 100, -1, 1 );
  du12Histo->setTitle( "x2-x1, |dy| < 1;x_{2}-x_{1} [mm];hit pairs" );

  dx23Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx23", 100, -1, 1 );
  dx23Histo->setTitle( "x3-x2;x_{3}-x_{2} [mm];hit pairs" );

  dy23Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dy23", 100, -1, 1 );
  dy23Histo->setTitle( "y3-y2;y_{3}-y_{2} [mm];hit pairs" );

  du23Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/du23", 100, -1, 1 );
  du23Histo->setTitle( "x3-x2, |dy| < 1;x_{3}-x_{2} [mm];hit pairs" );

  dx34Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx34", 100, -1, 1 );
  dx34Histo->setTitle( "x4-x3;x_{4}-x_{3} [mm];hit pairs" );

  dy34Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dy34", 100, -1, 1 );
  dy34Histo->setTitle( "y4-y3;y_{4}-y_{3} [mm];hit pairs" );

  du34Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/du34", 100, -1, 1 );
  du34Histo->setTitle( "x4-x3, |dy| < 1;x_{4}-x_{3} [mm];hit pairs" );

  dx45Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dx45", 100, -1, 1 );
  dx45Histo->setTitle( "x5-x4;x_{5}-x_{4} [mm];hit pairs" );

  dy45Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/dy45", 100, -1, 1 );
  dy45Histo->setTitle( "y5-y4;y_{5}-y_{4} [mm];hit pairs" );

  du45Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Telescope/du45", 100, -1, 1 );
  du45Histo->setTitle( "x5-x4, |dy| < 1;x_{5}-x_{4} [mm];hit pairs" );

  // triplets:
  AIDAProcessor::tree(this)->mkdir("Upstream");

  double dz = _planePosition[2] - _planePosition[0];
  double dx = 0.005 * dz;

  da02Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/da02", 100, -dx, dx );
  da02Histo->setTitle( "x2-x0, |dy| < 2 mm;x_{2}-x_{0} [mm];hit pairs" );

  db02Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/db02", 100, -dx, dx );
  db02Histo->setTitle( "y2-y0, |dx| < 2 mm;y_{2}-y_{0} [mm];hit pairs" );

  dzcvsxy = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "Upstream/dzcvsxy", 120, -12, 12, 60, -6, 6, -999, 999 );
  dzcvsxy->setTitle( "DUT plane;telescope track x_{DUT} [mm];telescope track y_{DUT} [mm];z_{DUT} [mm]" );

  z3vsxy = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "Upstream/z3vsxy", 120, -12, 12, 60, -6, 6, -999, 999 );
  z3vsxy->setTitle( "DUT plane;telescope track x_{DUT} [mm];telescope track y_{DUT} [mm];z_{DUT} [mm]" );

  tridxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx", 100, -100, 100 );
  tridxHisto->setTitle( "triplet dx;x_{1}-x_{m} [#mum];telescope triplets" );

  tridyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy", 100, -100, 100 );
  tridyHisto->setTitle( "triplet dy;y_{1}-y_{m} [#mum];telescope triplets" );

  tridx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx1", 100, -1, 1 );
  tridx1Histo->setTitle( "triplet dx;x_{1}-x_{t} [mm];telescope triplets" );

  tridy1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy1", 100, -1, 1 );
  tridy1Histo->setTitle( "triplet dy;y_{1}-y_{t} [mm];telescope triplets" );

  tridxvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridxvsx", 100, -10, 10, -100, 100 );
  tridxvsx->setTitle( "triplet x resid vs x;x [mm];triplet <#Deltax> [#mum]" );

  tridxvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridxvsy", 50, -5, 5, -100, 100 );
  tridxvsy->setTitle( "triplet x resid vs y;y [mm];triplet <#Deltax> [#mum]" );

  tridxvstx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridxvstx", 80, -2, 2, -100, 100 );
  tridxvstx->setTitle( "triplet x resid vs tx;t_{x} [mrad];triplet <#Deltax> [#mum]" );

  tridxvsty = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridxvsty", 80, -2, 2, -100, 100 );
  tridxvsty->setTitle( "triplet x resid vs ty;t_{y} [mrad];triplet <#Deltax> [#mum]" );

  tridyvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridyvsx", 100, -10, 10, -100, 100 );
  tridyvsx->setTitle( "triplet y resid vs x;x [mm];triplet <#Deltay> [#mum]" );

  tridyvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridyvsy", 50, -5, 5, -100, 100 );
  tridyvsy->setTitle( "triplet y resid vs y;y [mm];triplet <#Deltay> [#mum]" );

  tridyvstx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridyvstx", 80, -2, 2, -100, 100 );
  tridyvstx->setTitle( "triplet y resid vs tx;t_{x} [mrad];triplet <#Deltay> [#mum]" );

  tridyvsty = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Upstream/tridyvsty", 80, -2, 2, -100, 100 );
  tridyvsty->setTitle( "triplet y resid vs ty;t_{y} [mrad];triplet <#Deltay> [#mum]" );

  tridx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx3", 100, -1, 1 );
  tridx3Histo->setTitle( "triplet dx;x_{3}-x_{t} [mm];telescope triplets" );

  tridy3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy3", 100, -1, 1 );
  tridy3Histo->setTitle( "triplet dy;y_{3}-y_{t} [mm];telescope triplets" );

  tridx3bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx3b", 100, -250, 250 );
  tridx3bHisto->setTitle( "triplet dx;x_{3}-x_{t} [um];telescope triplets" );

  tridy3bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy3b", 100, -250, 250 );
  tridy3bHisto->setTitle( "triplet dy;y_{3}-y_{t} [um];telescope triplets" );

  tridx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx4", 100, -2, 2 );
  tridx4Histo->setTitle( "triplet dx;x_{4}-x_{t} [mm];telescope triplets" );

  tridy4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy4", 100, -2, 2 );
  tridy4Histo->setTitle( "triplet dy;y_{4}-y_{t} [mm];telescope triplets" );

  tridx4bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx4b", 100, -400, 400 );
  tridx4bHisto->setTitle( "triplet dx;x_{4}-x_{t} [um];telescope triplets" );

  tridy4bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy4b", 100, -400, 400 );
  tridy4bHisto->setTitle( "triplet dy;y_{4}-y_{t} [um];telescope triplets" );

  tridx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx5", 100, -3, 3 );
  tridx5Histo->setTitle( "triplet dx;x_{5}-x_{t} [mm];telescope triplets" );

  tridy5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy5", 100, -3, 3 );
  tridy5Histo->setTitle( "triplet dy;y_{5}-y_{t} [mm];telescope triplets" );

  tridx5bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridx5b", 100, -1000, 1000 );
  tridx5bHisto->setTitle( "triplet dx;x_{5}-x_{t} [um];telescope triplets" );

  tridy5bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tridy5b", 100, -1000, 1000 );
  tridy5bHisto->setTitle( "triplet dy;y_{5}-y_{t} [um];telescope triplets" );

  trixHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/trix", 240, -12, 12 );
  trixHisto->setTitle( "triplet x1;x1_{out} [mm];telescope triplets" );

  triyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/triy", 120, -6, 6 );
  triyHisto->setTitle( "triplet y1;y1_{up} [mm];telescope triplets" );

  trixyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "Upstream/trixy", 240, -12, 12, 120, -6, 6 );
  trixyHisto->setTitle( "triplet y1 vs x1;x1_{out} [mm];y1_{up} [mm];telescope triplets" );

  tritxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/tritx", 100, -10, 10 );
  tritxHisto->setTitle( "triplet slope x;#theta_{x} [mrad];telescope triplets" );

  trityHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/trity", 100, -10, 10 );
  trityHisto->setTitle( "triplet slope y;#theta_{y} [mrad];telescope triplets" );

  trixdutHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/trixdut", 240, -12, 12 );
  trixdutHisto->setTitle( "triplet at DUT;triplet x_{out} at DUT [mm];telescope triplets" );

  triydutHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Upstream/triydut", 120, -6, 6 );
  triydutHisto->setTitle( "triplet at DUT;triplet y_{up} at DUT [mm];telescope triplets" );

  trixydutHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "Upstream/trixydut", 240, -12, 12, 120, -6, 6 );
  trixydutHisto->setTitle( "triplet at DUT;triplet x_{out} at DUT [mm];triplet y_{up} at DUT [mm];telescope triplets" );

  // DUT pixel vs triplets:

  cmsxxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "cmsxx", 52, -0.5, 51.5, 110, -11, 11 );
  cmsxxHisto->setTitle( "x correlation;DUT cluster col;telescope triplet x [mm];clusters" );

  cmsyyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "cmsyy", 80, -0.5, 79.5, 55, -5.5, 5.5 );
  cmsyyHisto->setTitle( "y correlation;DUT cluster row;telescope triplet y [mm];clusters" );

  cmspxqHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmspxq", 100, 0, 25 );
  cmspxqHisto->setTitle( "DUT pixel charge linked;DUT pixel charge [ke];DUT linked pixels" );

  cmssxaHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmssxa", 440, -11, 11 );
  cmssxaHisto->setTitle( "Pixel + Telescope x;cluster + triplet #Sigmax [mm];clusters" );

  cmsdyaHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdya", 220, -5.5, 5.5 );
  cmsdyaHisto->setTitle( "Pixel - telescope y;cluster - triplet #Deltay [mm];clusters" );

  cmsdxaHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdxa", 440, -11, 11 );
  cmsdxaHisto->setTitle( "Pixel - Telescope x;cluster - triplet #Deltax [mm];clusters" );

  cmssyaHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmssya", 220, -5.5, 5.5 );
  cmssyaHisto->setTitle( "Pixel + telescope y;cluster + triplet #Sigmay [mm];clusters" );

  cmsdx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdx4", 440, -11, 11 );
  cmsdx4Histo->setTitle( "Pixel - Telescope x;cluster + triplet #Sigmax [mm];clusters" );

  cmsdy4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdy4", 220, -5.5, 5.5 );
  cmsdy4Histo->setTitle( "Pixel - telescope y;cluster - triplet #Deltay [mm];clusters" );

  cmsdxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdx", 200, -500, 500 );
  cmsdxHisto->setTitle( "Pixel + Telescope x;cluster + triplet #Sigmax [#mum];clusters" );

  cmsdyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdy", 500, -500, 500 );
  cmsdyHisto->setTitle( "Pixel - telescope y;cluster - triplet #Deltay [#mum];clusters" );

  cmsdxfHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdxf", 200, -500, 500 );
  cmsdxfHisto->setTitle( "fiducial Pixel - telescope x;fiducial cluster - triplet #Deltax [#mum];fiducial clusters" );

  cmsdyfHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyf", 500, -500, 500 );
  cmsdyfHisto->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdxfcHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdxfc", 200, -500, 500 );
  cmsdxfcHisto->setTitle( "fiducial Pixel - telescope x;fiducial cluster - triplet #Deltax [#mum];fiducial clusters" );

  cmsdyfcHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfc", 500, -500, 500 );
  cmsdyfcHisto->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdyfc1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfc1", 500, -500, 500 );
  cmsdyfc1Histo->setTitle( "Pixel - telescope y;1-row cluster - triplet #Deltay [#mum];fiducial 1-row clusters" );

  cmsdyfc2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfc2", 500, -500, 500 );
  cmsdyfc2Histo->setTitle( "Pixel - telescope y;2-row cluster - triplet #Deltay [#mum];fiducial 2-row clusters" );

  cmsdyfc3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfc3", 500, -500, 500 );
  cmsdyfc3Histo->setTitle( "Pixel - telescope y;3-row cluster - triplet #Deltay [#mum];fiducial 3-row clusters" );

  cmsdyq0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyq0", 500, -500, 500 );
  cmsdyq0Histo->setTitle( "Pixel - telescope y, Q < 18;cluster - triplet #Deltay [#mum];low dE/dx clusters" );

  cmsdyq1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyq1", 500, -500, 500 );
  cmsdyq1Histo->setTitle( "Pixel - telescope y, 18 < Q < 40;cluster - triplet #Deltay [#mum];med dE/dx clusters" );

  cmsdyeta0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyeta0", 500, -500, 500 );
  cmsdyeta0Histo->setTitle( "Pixel - telescope y, neg eta;cluster - triplet #Deltay [#mum];neg eta clusters" );

  cmsdyeta1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyeta1", 500, -500, 500 );
  cmsdyeta1Histo->setTitle( "Pixel - telescope y, pos eta;cluster - triplet #Deltay [#mum];pos eta clusters" );

  cmsdyq2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyq2", 500, -500, 500 );
  cmsdyq2Histo->setTitle( "Pixel - telescope y, Q > 40;cluster - triplet #Deltay [#mum];high dE/dx clusters" );

  cmsdxfctHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdxfct", 500, -500, 500 );
  cmsdxfctHisto->setTitle( "fiducial Pixel - telescope x;fiducial cluster - triplet #Deltax [#mum];fiducial clusters" );

  cmsdyfctHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfct", 500, -500, 500 );
  cmsdyfctHisto->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdyfcntHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfcnt", 500, -500, 500 );
  cmsdyfcntHisto->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial 1-,2-row clusters" );

  cmsdxfctqHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdxfctq", 500, -500, 500 );
  cmsdxfctqHisto->setTitle( "fiducial Pixel - telescope x;fiducial cluster - triplet #Deltax [#mum];fiducial clusters" );

  cmsdyfctqHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfctq", 500, -500, 500 );
  cmsdyfctqHisto->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdyfcntqHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfcntq", 500, -500, 500 );
  cmsdyfcntqHisto->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial 1-,2-row clusters" );

  cmsdxfctq1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdxfctq1", 500, -500, 500 );
  cmsdxfctq1Histo->setTitle( "fiducial Pixel - telescope x;fiducial cluster - triplet #Deltax [#mum];fiducial clusters" );

  cmsdyfctq1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfctq1", 500, -500, 500 );
  cmsdyfctq1Histo->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdyfcntq1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfcntq1", 500, -500, 500 );
  cmsdyfcntq1Histo->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial 1-,2-row clusters" );

  cmsdyfctq1lHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfctq1l", 500, -500, 500 );
  cmsdyfctq1lHisto->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdyfctq1rHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfctq1r", 500, -500, 500 );
  cmsdyfctq1rHisto->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdxfctq2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdxfctq2", 500, -500, 500 );
  cmsdxfctq2Histo->setTitle( "fiducial Pixel - telescope x;fiducial cluster - triplet #Deltax [#mum];fiducial clusters" );

  cmsdyfctq2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfctq2", 500, -500, 500 );
  cmsdyfctq2Histo->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdxfctq3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdxfctq3", 500, -500, 500 );
  cmsdxfctq3Histo->setTitle( "fiducial Pixel - telescope x;fiducial cluster - triplet #Deltax [#mum];fiducial clusters" );

  cmsdyfctq3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfctq3", 500, -500, 500 );
  cmsdyfctq3Histo->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdyfctqdotHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfctqdot", 500, -500, 500 );
  cmsdyfctqdotHisto->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );

  cmsdyfctq3dHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdyfctq3d", 500, -500, 500 );
  cmsdyfctq3dHisto->setTitle( "fiducial Pixel - telescope y;fiducial cluster - triplet #Deltay [#mum];fiducial clusters" );


  cmscolHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmscol", 52, -0.5, 51.5 );
  cmscolHisto->setTitle( "DUT linked column;DUT linked cluster col;DUT linked clusters" );

  cmsrowHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsrow", 80, -0.5, 79.5 );
  cmsrowHisto->setTitle( "DUT linked row;DUT linked cluster row;DUT linked clusters" );

  cmsqHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsq", 100, 0, 100 );
  cmsqHisto->setTitle( "DUT cluster charge linked;DUT cluster charge [ke];DUT linked clusters" );

  cmsq0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsq0", 100, 0, 100 );
  cmsq0Histo->setTitle( "DUT cluster charge linked;normal DUT cluster charge [ke];DUT linked clusters" );

  trixlkHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "trixlk", 240, -12, 12 );
  trixlkHisto->setTitle( "DUT - triplet link;triplet x_{out} [mm];linked clusters" );

  triylkHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "triylk", 120, -6, 6 );
  triylkHisto->setTitle( "DUT - triplet link;triplet y_{up} [mm];linked clusters" );

  trixylkHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "trixylk", 240, -12, 12, 120, -6, 6 );
  trixylkHisto->setTitle( "DUT - triplet link;triplet x_{out} [mm];triplet y_{up} [mm];linked clusters" );


  cmsdxvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsdxvsx", 52, -26*0.15, 26*0.15, -150, 150 ); // 16.2.2013
  cmsdxvsx->setTitle( "DUT x resid vs x;x [mm];<DUT cluster - telescope triplet #Deltax> [#mum]" );

  cmsdyvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsdyvsx", 52, -26*0.15, 26*0.15, -150, 150 );
  cmsdyvsx->setTitle( "DUT y resid vs x;x [mm];<DUT cluster - telescope triplet #Deltay> [#mum]" );

  cmsdxvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsdxvsy", 80, -4, 4, -150, 150 );
  cmsdxvsy->setTitle( "DUT x resid vs y;y [mm];<DUT cluster - telescope triplet #Deltax> [#mum]" );

  cmsdyvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsdyvsy", 80, -4, 4, -150, 150 );
  cmsdyvsy->setTitle( "DUT y resid vs y;y [mm];<DUT cluster - telescope triplet #Deltay> [#mum]" );

  cmsdxvstx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsdxvstx", 80, -2, 2, -150, 150 );
  cmsdxvstx->setTitle( "DUT x resid vs tet x;#theta_{x} [mrad];<DUT cluster - telescope triplet #Deltax> [#mum]" );

  cmsdyvsty = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsdyvsty", 80, -2, 2, -150, 150 );
  cmsdyvsty->setTitle( "DUT y resid vs tet y;#theta_{y} [mrad];<DUT cluster - telescope triplet #Deltay> [#mum]" );

  cmsdyvsxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "cmsdyvsxh", 52, -26*0.15, 26*0.15, 50, -250, 250 );
  cmsdyvsxHisto->setTitle( "DUT resolution;x [mm];#Deltay [#mum];clusters" );

  cmsnpxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsnpx", 11, -0.5, 10.5 );
  cmsnpxHisto->setTitle( "linked CMS cluster size;CMS pixel per linked cluster;linked CMS clusters" );

  cmsnpx0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsnpx0", 11, -0.5, 10.5 );
  cmsnpx0Histo->setTitle( "linked CMS cluster size, bias dot;CMS pixel per linked cluster;linked CMS clusters in bias dot" );

  cmsnpx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsnpx1", 11, -0.5, 10.5 );
  cmsnpx1Histo->setTitle( "linked CMS cluster size, no bias dot;CMS pixel per linked cluster;linked CMS clusters no bias dot" );

  cmsnpx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsnpx2", 11, -0.5, 10.5 );
  cmsnpx2Histo->setTitle( "linked CMS cluster size, no bias dot, core;CMS pixel per linked cluster;linked CMS clusters no bias dot, core" );

  cmsnpx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsnpx3", 11, -0.5, 10.5 );
  cmsnpx3Histo->setTitle( "linked CMS cluster size, no bias dot, edge;CMS pixel per linked cluster;linked CMS clusters no bias dot, edge" );

  cmsncolHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsncol", 11, -0.5, 10.5 );
  cmsncolHisto->setTitle( "DUT cluster col size;columns per cluster;DUT clusters" );

  cmsnrowHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsnrow", 11, -0.5, 10.5 );
  cmsnrowHisto->setTitle( "DUT cluster row size;rows per cluster;DUT clusters" );

  cmsnrowqHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsnrowq", 11, -0.5, 10.5 );
  cmsnrowqHisto->setTitle( "DUT Landau peak cluster row size;rows per cluster in Landau peak ;DUT Landau peak clusters" );

  cmsnrowvst1 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsnrowvst1", 70, 0, 350, 0, 5.5 );
  cmsnrowvst1->setTitle( "DUT cluster rows vs time;run time [s];<cluster rows>" );

  cmsetaHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmseta", 100, -1, 1 );
  cmsetaHisto->setTitle( "DUT 2-row eta;DUT cluster 2-row eta;DUT 2-row clusters" );

  cmsqfHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsqf", 100, 0, 100 );
  cmsqfHisto->setTitle( "DUT cluster charge linked fiducial;DUT cluster charge [ke];DUT linked fiducial clusters" );

  cmsq0fHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsq0f", 100, 0, 100 );
  cmsq0fHisto->setTitle( "DUT cluster charge linked fiducial;normal DUT cluster charge [ke];DUT linked fiducial clusters" );

  cmsqf0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsqf0", 100, 0, 100 );
  cmsqf0Histo->setTitle( "DUT cluster charge linked fiducial bias dot;DUT cluster charge [ke];DUT linked fiducial clusters bias dot" );

  cmsqf1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsqf1", 100, 0, 100 );
  cmsqf1Histo->setTitle( "DUT cluster charge linked fiducial no dot;DUT cluster charge [ke];DUT linked fiducial clusters no dot" );

  cmsqf2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsqf2", 100, 0, 100 );
  cmsqf2Histo->setTitle( "DUT cluster charge linked fiducial no dot core;DUT cluster charge [ke];DUT linked fiducial clusters no dot core" );

  cmsqf3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsqf3", 100, 0, 100 );
  cmsqf3Histo->setTitle( "DUT cluster charge linked fiducial no dot edge;DUT cluster charge [ke];DUT linked fiducial clusters no dot edge" );

  cmsdyvsxm = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsdyvsxm", 60, 0, 300, -150, 150 );
  cmsdyvsxm->setTitle( "DUT y resid vs xmod;telescope x mod 300 [#mum];<DUT cluster - telescope triplet #Deltay> [#mum]" );

  cmsdyvsym = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsdyvsym", 40, 0, 200, -150, 150 );
  cmsdyvsym->setTitle( "DUT y resid vs ymod;telescope y mod 200 [#mum];<DUT cluster - telescope triplet #Deltay> [#mum]" );

  cmspixvsxmym = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "cmspixvsxmym", 60, 0, 300, 40, 0, 200 );
  cmspixvsxmym->setTitle( "DUT pixel occupancy;x mod 300 #mum;y mod 200 #mum;clusters" );

  cmsqvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvsx", 380, -3.8, 3.8, 0, 100 );
  cmsqvsx->setTitle( "DUT q vs x;telescope x at DUT [mm];<q> [ke]" );

  cmsqvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvsy", 380, -3.8, 3.8, 0, 100 );
  cmsqvsy->setTitle( "DUT q vs y;telescope y at DUT [mm];<q> [ke]" );

  cmsqvsxm = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvsxm", 60, 0, 300, 0, 100 );
  cmsqvsxm->setTitle( "DUT q vs xmod;telescope x_{DUT} mod 300 [#mum];<q> [ke]" );

  cmsqvsym = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvsym", 40, 0, 200, 0, 100 );
  cmsqvsym->setTitle( "DUT q vs ymod;telescope y_{DUT} mod 200 [#mum];<q> [ke]" );

  cmsqvsxmym = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "cmsqvsxmym", 60, 0, 300, 40, 0, 200, 0, 250 );
  cmsqvsxmym->setTitle( "DUT cluster charge map;x_{track} mod 300 #mum;y_{track} mod 200 #mum;<cluster charge> [ke]" );

  cmsqvsddt = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvsddt", 40, -40, 40, 0, 250 );
  cmsqvsddt->setTitle( "DUT cluster charge vs phase;TLU-DUT #delta#Deltat [ns];<cluster charge> [ke]" );

  cmsqvst1 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvst1", 70, 0, 350, 0, 250 );
  cmsqvst1->setTitle( "DUT cluster charge vs time;run time [s];<cluster charge> [ke]" );

  cmsqvst2 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvst2", 140, 0, 700, 0, 250 );
  cmsqvst2->setTitle( "DUT cluster charge vs time;run time [s];<cluster charge> [ke]" );

  cmsqvst3 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvst3", 200, 0, 2000, 0, 250 );
  cmsqvst3->setTitle( "DUT cluster charge vs time;run time [s];<cluster charge> [ke]" );

  cmsqvst4 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsqvst4", 400, 0, 4000, 0, 250 );
  cmsqvst4->setTitle( "DUT cluster charge vs time;run time [s];<cluster charge> [ke]" );

  cmsrmsxvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsxvsx", 152, -3.8, 3.8, 0, 100 );
  cmsrmsxvsx->setTitle( "DUT x resolution vs x;telescope x [mm];RMS(#Deltax) [#mum]" );

  cmsrmsyvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsyvsx", 152, -3.8, 3.8, 0, 100 );
  cmsrmsyvsx->setTitle( "DUT y resolution vs x;telescope x [mm];RMS(#Deltay) [#mum]" );

  cmsrmsxvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsxvsy", 381, -3.81, 3.81, 0, 100 ); // shift by 10 um
  cmsrmsxvsy->setTitle( "DUT x resolution vs y;telescope y [mm];RMS(#Deltax) [#mum]" );

  cmsrmsyvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsyvsy", 381, -3.81, 3.81, 0, 100 ); // shift by 10 um
  cmsrmsyvsy->setTitle( "DUT y resolution vs y;telescope y [mm];RMS(#Deltay) [#mum]" );

  cmsrmsxvsxm = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsxvsxm", 60, 0, 300, 0, 100 );
  cmsrmsxvsxm->setTitle( "DUT x resolution;telescope x mod 300 [#mum];RMS(#Deltax) [#mum]" );

  cmsrmsyvsxm = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsyvsxm", 60, 0, 300, 0, 100 );
  cmsrmsyvsxm->setTitle( "DUT y resolution vs x mod 300;telescope x mod 300 [#mum];RMS(#Deltay) [#mum]" );

  cmsncolvsxm = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsncolvsxm", 60, 0, 300, 0, 10 );
  cmsncolvsxm->setTitle( "DUT cluster cols vs x mod 300;telescope x mod 300 [#mum];<cols/cluster>" );

  cmsnrowvsxm = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsnrowvsxm", 60, 0, 300, 0, 10 );
  cmsnrowvsxm->setTitle( "DUT cluster rows vs x mod 300;telescope x mod 300 [#mum];<rows/cluster>" );

  cmsrmsxvsym = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsxvsym", 40, 0, 200, 0, 100 );
  cmsrmsxvsym->setTitle( "DUT x resolution vs y mod 200;telescope y mod 200 [#mum];RMS(#Deltax) [#mum]" );

  cmsrmsyvsym = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsyvsym", 40, 0, 200, 0, 100 );
  cmsrmsyvsym->setTitle( "DUT y resolution vs y mod 200;telescope y mod 200 [#mum];RMS(#Deltay) [#mum]" );

  cmsrmsyvsym3 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsyvsym3", 60, 0, 300, 0, 100 );
  cmsrmsyvsym3->setTitle( "DUT y resolution vs y mod 300;telescope y mod 300 [#mum];RMS(#Deltay) [#mum]" );

  cmsrmsyvsym6 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsyvsym6", 120, 0, 600, 0, 100 );
  cmsrmsyvsym6->setTitle( "DUT y resolution vs y mod 600;telescope y mod 600 [#mum];RMS(#Deltay) [#mum]" );

  cmsrmsyvst = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsyvst", 60, 0, 3600, 0, 100 );
  cmsrmsyvst->setTitle( "DUT y resolution;time [s];RMS(#Deltay) [#mum]" );

  cmsrmsyvsddt = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsyvsddt", 50, -25, 25, 0, 100 );
  cmsrmsyvsddt->setTitle( "DUT y resolution vs ddt;TLU-DUT #delta#Deltat [ns];RMS(#Deltay) [#mum]" );

  cmsrmsxvsq = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsxvsq", 70, 10, 80, 0, 100 );
  cmsrmsxvsq->setTitle( "DUT x resolution vs q;cluster charge [ke];RMS(#Deltax) [#mum]" );

  cmsrmsyvsq = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsyvsq", 70, 10, 80, 0, 100 );
  cmsrmsyvsq->setTitle( "DUT y resolution vs q;cluster charge [ke];RMS(#Deltay) [#mum]" );

  cmsdyvseta = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsdyvseta", 100, -1, 1, -200, 200 );
  cmsdyvseta->setTitle( "DUT dy vs eta;charge sharing eta;<#Deltay> [#mum]" );

  cmsrmsyvseta = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsyvseta", 100, -1, 1, 0, 100 );
  cmsrmsyvseta->setTitle( "DUT y resolution vs eta;charge sharing eta;RMS(#Deltay) [#mum]" );

  cmspMoyalvsq = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmspMoyalvsq", 100, 0, 100, 0, 1 );
  cmspMoyalvsq->setTitle( "cluster Moyal prob vs q;cluster charge [ke];<Moyal prob>" );

  cmspMoyalHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmspMoyal", 100, 0, 1 );
  cmspMoyalHisto->setTitle( "cluster Moyal prob;Moyal prob;DUT linked clusters" );

  cmsrmsyvsp = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsrmsyvsp", 100, 0, 1, 0, 100 );
  cmsrmsyvsp->setTitle( "DUT y resolution vs prob;cluster charge [Moyal prob];RMS(#Deltay) [#mum]" );

  cmsnpxvsxmym = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "cmsnpxvsxmym", 60, 0, 300, 40, 0, 200, 0, 4.5 );
  cmsnpxvsxmym->setTitle( "DUT cluster size vs xm-ym;telescope x mod 300 [#mum];telescope y mod 200 [#mum];<pixels/cluster>" );

  cmsncolvsym = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsncolvsym", 40, 0, 200, 0, 4.5 );
  cmsncolvsym->setTitle( "DUT cols vs xmod;telescope x mod 300 [#mum];<cols/cluster>" );

  cmsnrowvsym = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsnrowvsym", 40, 0, 200, 0, 4.5 );
  cmsnrowvsym->setTitle( "DUT rows vs ymod;telescope y mod 200 [#mum];<rows/cluster>" );

  cmsetavsym = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsetavsym", 40, 0, 200, -1.5, 1.5 );
  cmsetavsym->setTitle( "DUT eta vs ymod;telescope y mod 200 [#mum];<eta>" );

  cmsetavsym2 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsetavsym2", 40, 0, 200, -1.5, 1.5 );
  cmsetavsym2->setTitle( "DUT eta 2-row vs ymod;telescope y mod 200 [#mum];<eta> 2-row" );

  cmsetavsym3 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "cmsetavsym3", 60, 0, 300, -1.5, 1.5 );
  cmsetavsym3->setTitle( "DUT eta vs ymod;telescope y mod 200 [#mum];<eta>" );

  cmsym1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsym1", 40, 0, 200 );
  cmsym1Histo->setTitle( "ymod at DUT 1-row clus;telescope y mod 200 [#mum];1-row clusters" );

  cmsym2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsym2", 40, 0, 200 );
  cmsym2Histo->setTitle( "ymod at DUT 2-row clus;telescope y mod 200 [#mum];2-row clusters" );


  effxyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "effxy", 60, -4.5, 4.5, 90, -4.5, 4.5 );
  effxyHisto->setTitle( "six with REF link;triplet x_{DUT} [mm];triplet y_{DUT} [mm];tracks" );

  effvsxy = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "effvsxy", 60, -4.5, 4.5, 90, -4.5, 4.5, -1, 2 );
  effvsxy->setTitle( "DUT efficiency;telescope track x_{DUT} [mm];telescope track y_{DUT} [mm];CMS DUT efficiency" );

  effvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "effvsx", 60, -4.5, 4.5, -1, 2 );
  effvsx->setTitle( "DUT efficiency;telescope track x_{DUT} [mm];CMS DUT efficiency" );

  effvsxg = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "effvsxg", 60, -4.5, 4.5, -1, 2 );
  effvsxg->setTitle( "DUT efficiency;telescope track x_{DUT} [mm];CMS DUT efficiency" );

  effvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "effvsy", 90, -4.5, 4.5, -1, 2 );
  effvsy->setTitle( "DUT efficiency;telescope track y_{DUT} [mm];<CMS DUT efficiency>" );

  eff300 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "eff300",  70, 0, 350, -1, 2 );
  eff300->setTitle( "DUT efficiency;time [s];<CMS DUT efficiency>" );

  eff600 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "eff600",  70, 0, 700, -1, 2 );
  eff600->setTitle( "DUT efficiency;time [s];<CMS DUT efficiency>" );

  effd600 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "effd600",  70, 0, 700, -1, 2 );
  effd600->setTitle( "DUT efficiency DUT-REF ddt=0;time [s];<CMS DUT efficiency>" );

  effn600 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "effn600",  70, 0, 700, -1, 2 );
  effn600->setTitle( "DUT efficiency ndri < 3;time [s];<CMS DUT efficiency>" );

  effm600 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "effm600",  70, 0, 700, -1, 2 );
  effm600->setTitle( "DUT efficiency ndri = 1;time [s];<CMS DUT efficiency>" );

  eff1200 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "eff1200", 130, 0, 1300, -1, 2 );
  eff1200->setTitle( "DUT efficiency;time [s];<CMS DUT efficiency>" );

  eff1800 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "eff1800", 200, 0, 2000, -1, 2 );
  eff1800->setTitle( "DUT efficiency;time [s];<CMS DUT efficiency>" );

  eff3600 = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "eff3600", 400, 0, 4000, -1, 2 );
  eff3600->setTitle( "DUT efficiency;time [s];<CMS DUT efficiency>" );

  effvsddt = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "effvsddt", 40, -40, 40, -1, 2 );
  effvsddt->setTitle( "DUT efficiency;TLU-DUT #delta#Deltat [ns];<CMS DUT efficiency>" );

  effvsxmym = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "effvsxmym", 60, 0, 300, 40, 0, 200, -1, 2 );
  effvsxmym->setTitle( "DUT efficiency;telescope x' mod 300 [um];telescope y' mod 200 [um];efficiency" );


  rffvsxy = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "rffvsxy", 120, -12, 12, 60, -6, 6, -1, 2 );
  rffvsxy->setTitle( "REF efficiency;telescope x [mm];telescope y [mm];DUT REF efficiency" );

  rffvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "rffvsx", 100, -2, 8, -1, 2 );
  rffvsx->setTitle( "REF efficiency;telescope x [mm];DUT REF efficiency" );


  lkAvst = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "lkAvst", 360, 0, 3600, -0.5, 2.5 );
  lkAvst->setTitle( "DUT links/event;time [s];<tri-DUT links/event>" );

  ntriHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "ntri", 31, -0.5, 30.5 );
  ntriHisto->setTitle( "telescope triplets;0-1-2 triplets;events" );

  // triplet efficiency w.r.t. DUT:

  cmsxeHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsxe", 240, -12, 12 );
  cmsxeHisto->setTitle( "DUT hits;aligned x [mm];clusters" );

  cmsyeHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsye", 240, -6, 6 );
  cmsyeHisto->setTitle( "DUT hits;aligned y [mm];clusters" );

  cmsdxeHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdxe", 200, -500, 500 );
  cmsdxeHisto->setTitle( "Pixel + Triplet x;cluster + triplet #Sigmax [#mum];clusters" );

  cmsdyeHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsdye", 500, -500, 500 );
  cmsdyeHisto->setTitle( "Pixel - triplet y;cluster - triplet #Deltay [#mum];clusters" );

  cmsnmHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "cmsnm", 11, -0.5, 10.5 );
  cmsnmHisto->setTitle( "triplets linked to clusters;linked triplets;clusters" );

  trieffvsxy = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "trieffvsxy", 120, -12, 12, 60, -6, 6, -1, 2 );
  trieffvsxy->setTitle( "triplet efficiency;DUT cluster x [mm];DUT cluster y [mm];efficiency" );

  // driplets:
  AIDAProcessor::tree(this)->mkdir("Downstream");

  dz = _planePosition[5] - _planePosition[3];
  dx = 0.01 * dz;

  dx35Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/dx35", 100, -dx, dx );
  dx35Histo->setTitle( "x5-x3, |dy| < 2 mm;x_{5}-x_{3} [mm];hit pairs" );

  dy35Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/dy35", 100, -dx, dx );
  dy35Histo->setTitle( "y5-y3, |dx| < 2 mm;y_{5}-y_{3} [mm];hit pairs" );

  dridxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/dridx", 100, -100, 100 );
  dridxHisto->setTitle( "driplet dx;x_{4}-x_{m} [#mum];driplets" );

  dridyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/dridy", 100, -100, 100 );
  dridyHisto->setTitle( "driplet dy;y_{4}-y_{m} [#mum];driplets" );

  dridxvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridxvsx", 100, -10, 10, -100, 100 );
  dridxvsx->setTitle( "driplet x resid vs x;x [mm];driplet <#Deltax> [#mum]" );

  dridxvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridxvsy", 50, -5, 5, -100, 100 );
  dridxvsy->setTitle( "driplet x resid vs y;y [mm];driplet <#Deltax> [#mum]" );

  dridxvstx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridxvstx", 80, -2, 2, -100, 100 );
  dridxvstx->setTitle( "driplet x resid vs tx;t_{x} [mrad];driplet <#Deltax> [#mum]" );

  dridxvsty = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridxvsty", 80, -2, 2, -100, 100 );
  dridxvsty->setTitle( "driplet x resid vs ty;t_{y} [mrad];driplet <#Deltax> [#mum]" );

  dridyvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridyvsx", 100, -10, 10, -100, 100 );
  dridyvsx->setTitle( "driplet y resid vs x;x [mm];driplet <#Deltay> [#mum]" );

  dridyvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridyvsy", 50, -5, 5, -100, 100 );
  dridyvsy->setTitle( "driplet y resid vs y;y [mm];driplet <#Deltay> [#mum]" );

  dridyvstx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridyvstx", 80, -2, 2, -100, 100 );
  dridyvstx->setTitle( "driplet y resid vs tx;t_{x} [mrad];driplet <#Deltay> [#mum]" );

  dridyvsty = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Downstream/dridyvsty", 80, -2, 2, -100, 100 );
  dridyvsty->setTitle( "driplet y resid vs ty;t_{y} [mrad];driplet <#Deltay> [#mum]" );

  drixHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/drix", 240, -12, 12 );
  drixHisto->setTitle( "driplet x4;x4_{out} [mm];driplets" );

  driyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/driy", 120, -6, 6 );
  driyHisto->setTitle( "driplet y4;y4_{up} [mm];driplets" );

  drixyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "Downstream/drixy", 240, -12, 12, 120, -6, 6 );
  drixyHisto->setTitle( "driplet y4 vs x4;x4_{out} [mm];y4_{up} [mm];driplets" );

  dritxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/dritx", 100, -10, 10 );
  dritxHisto->setTitle( "driplet slope x;#theta_{x} [mrad];driplets" );

  drityHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/drity", 100, -10, 10 );
  drityHisto->setTitle( "driplet slope y;#theta_{y} [mrad];driplets" );




  drixrefHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "drixref", 240, -12, 12 );
  drixrefHisto->setTitle( "driplet at REF;driplet x_{out} at REF [mm];driplets" );

  driyrefHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "driyref", 120, -6, 6 );
  driyrefHisto->setTitle( "driplet at REF;driplet y_{up} at REF [mm];driplets" );

  drixyrefHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "drixyref", 240, -12, 12, 120, -6, 6 );
  drixyrefHisto->setTitle( "driplet at REF;driplet x_{out} at REF [mm];driplet y_{up} at REF [mm];driplets" );


  drixlkHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "drixlk", 240, -12, 12 );
  drixlkHisto->setTitle( "driplet with REF link;x_{REF} [mm];linked driplets" );

  driylkHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "driylk", 120, -6, 6 );
  driylkHisto->setTitle( "driplet with REF link;y_{REF} [mm];linked driplets" );

  drixylkHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "drixylk", 240, -12, 12, 120, -6, 6 );
  drixylkHisto->setTitle( "driplet with REF link;driplet x_{out} at REF [mm];driplet y_{up} at REF [mm];linked driplets" );

  refpixvsxmym = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "refpixvsxmym", 60, 0, 300, 40, 0, 200 );
  refpixvsxmym->setTitle( "REF pixel occupancy;x mod 300 #mum;y mod 200 #mum;clusters" );

  refqHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "refq", 100, 0, 100 );
  refqHisto->setTitle( "REF cluster charge linked;REF cluster charge [ke];REF linked clusters" );

  refqvsxmym = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "refqvsxmym", 60, 0, 300, 40, 0, 200, 0, 100 );
  refqvsxmym->setTitle( "REF cluster charge map;x mod 300 #mum;y mod 200 #mum;<cluster charge> [ke]" );


  refxxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "refxx", 52, -0.5, 51.5, 110, -11, 11 );
  refxxHisto->setTitle( "x correlation;REF pixel cluster col;driplet x [mm];REF clusters" );

  refyyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "refyy", 80, -0.5, 79.5, 55, -5.5, 5.5 );
  refyyHisto->setTitle( "y correlation;REF pixel cluster row;driplet y [mm];REF clusters" );

  refsxaHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "refsxa", 200, -10, 10 );
  refsxaHisto->setTitle( "REF + driplet x;cluster + driplet #Sigmax [mm];REF clusters" );

  refdyaHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "refdya", 200, -10, 10 );
  refdyaHisto->setTitle( "REF - driplet y;cluster - driplet #Deltay [mm];REF clusters" );

  refsxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "refsx", 100, -5, 5 );
  refsxHisto->setTitle( "REF + driplet x;cluster + driplet #Sigmax [mm];REF clusters" );

  refdyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "refdy", 100, -5, 5 );
  refdyHisto->setTitle( "REF - driplet y;cluster - driplet #Deltay [mm];REF clusters" );

  refsxcHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "refsxc", 100, -1, 1 );
  refsxcHisto->setTitle( "REF + driplet x;cluster + driplet #Sigmax [mm];REF clusters" );

  refdycHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "refdyc", 100, -1, 1 );
  refdycHisto->setTitle( "REF - driplet y;cluster - driplet #Deltay [mm];REF clusters" );

  refdyvsx = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "refdyvsx", 52, -26*0.15, 26*0.15, -150, 150 );
  refdyvsx->setTitle( "REF y resid vs x;x [mm];<REF cluster - driplet #Deltay> [#mum]" );

  refdyvsy = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "refdyvsy", 80, -4, 4, -150, 150 );
  refdyvsy->setTitle( "REF y resid vs y;y [mm];<REF cluster - driplet #Deltay> [#mum]" );

  refdyvsty = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "refdyvsty", 80, -2, 2, -150, 150 );
  refdyvsty->setTitle( "REF y resid vs tet y;#theta_{y} [mrad];<REF cluster - driplet #Deltay> [#mum]" );

  reflkcolHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "reflkcol", 52, -0.5, 51.5 );
  reflkcolHisto->setTitle( "REF linked column;REF linked cluster col;REF linked clusters" );

  reflkrowHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "reflkrow", 80, -0.5, 79.5 );
  reflkrowHisto->setTitle( "REF linked row;REF linked cluster row;REF linked clusters" );


  bacsxaHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "bacsxa", 100, -10, 10 );
  bacsxaHisto->setTitle( "DUT + driplet x;DUT cluster + driplet #Sigmax [mm];DUT clusters" );

  bacdyaHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "bacdya", 100, -10, 10 );
  bacdyaHisto->setTitle( "DUT - driplet y;DUT cluster - driplet #Deltay [mm];DUT clusters" );

  bacsxcHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "bacsxc", 200, -500, 500 );
  bacsxcHisto->setTitle( "DUT + driplet x;DUT cluster + driplet #Sigmax [mm];DUT clusters" );

  bacdycHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "bacdyc", 200, -500, 500 );
  bacdycHisto->setTitle( "DUT - driplet y;DUT cluster - driplet #Deltay [mm];DUT clusters" );

  bacsxcqHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "bacsxcq", 200, -500, 500 );
  bacsxcqHisto->setTitle( "DUT + driplet x;DUT cluster + driplet #Sigmax [mm];DUT clusters" );

  bacdycqHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "bacdycq", 100, -500, 500 );
  bacdycqHisto->setTitle( "DUT - driplet y;DUT cluster - driplet #Deltay [mm];DUT clusters" );


  ndriHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Downstream/ndri", 31, -0.5, 30.5 );
  ndriHisto->setTitle( "telescope driplets;3-4-5 driplets;events" );

  ndrirefHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "ndriref", 31, -0.5, 30.5 );
  ndrirefHisto->setTitle( "telescope driplets and DUT hit in event;3-4-5 driplets and DUT;events" );

  lkBvst = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "lkBvst", 360, 0, 3600, -0.5, 2.5 );
  lkBvst->setTitle( "REF links/event;time [s];<dri-REF links/event>" );

  //driplets-triplets
  // Tracks
  AIDAProcessor::tree(this)->mkdir("Tracks");

  nsixHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/nsix", 21, -0.5, 20.5 );
  nsixHisto->setTitle( "telescope six-plane-tracks;six-plane-tracks;events" );

  sixkxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixkx", 100, -10, 10 );
  sixkxHisto->setTitle( "kink x;kink x [mrad];triplet-driplet pairs" );

  sixkyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixky", 100, -10, 10 );
  sixkyHisto->setTitle( "kink y;kink y [mrad];triplet-driplet pairs" );

  sixdxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixdx", 100, -1, 1 );
  sixdxHisto->setTitle( "six match x;match x [mm];triplet-driplet pairs" );

  sixdyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixdy", 100, -1, 1 );
  sixdyHisto->setTitle( "six match y;match y [mm];triplet-driplet pairs" );

  sixdxcHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixdxc", 100, -250, 250 );
  sixdxcHisto->setTitle( "six match x;track #Deltax[#mum];triplet-driplet pairs" );

  sixdycHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixdyc", 100, -250, 250 );
  sixdycHisto->setTitle( "six match y;track #Deltay[#mum];triplet-driplet pairs" );

  sixkxcHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixkxc", 100, -10, 10 );
  sixkxcHisto->setTitle( "kink x, x-y matched;kink x [mrad];triplet-driplet pairs" );

  sixkycHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixkyc", 100, -10, 10 );
  sixkycHisto->setTitle( "kink y, x-y matched;kink y [mrad];triplet-driplet pairs" );

  sixxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixx", 240, -12, 12 );
  sixxHisto->setTitle( "six x at DUT;six x_{out} at DUT [mm];six-plane tracks" );

  sixyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixy", 120, -6, 6 );
  sixyHisto->setTitle( "six y at DUT;six y_{up} at DUT [mm];six-plane tracks" );

  sixxyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "Tracks/sixxy", 240, -12, 12, 120, -6, 6 );
  sixxyHisto->setTitle( "six at DUT;six x_{out} at DUT [mm];six y_{up} at DUT [mm];six-plane tracks" );

  sixxycHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "Tracks/sixxyc", 240, -12, 12, 120, -6, 6 );
  sixxycHisto->setTitle( "six large kink;six x_{out} at DUT [mm];six y_{up} at DUT [mm];large kink tracks" );

  sixxylkHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram2D( "Tracks/sixxylk", 240, -12, 12, 120, -6, 6 );
  sixxylkHisto->setTitle( "six with REF link at DUT;six x_{out} at DUT [mm];six y_{up} at DUT [mm];six-plane tracks with REF link" );

  kinkvsxy = AIDAProcessor::histogramFactory(this)->
    createProfile2D( "Tracks/kinkvsxy", 120, -12, 12, 60, -6, 6, 0, 100 );
  kinkvsxy->setTitle( "kink;six x_{out} at DUT [mm];six y_{up} at DUT [mm];<kink^{2}> [mrad^{2}]" );

  sixx0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixx0", 240, -12, 12 );
  sixx0Histo->setTitle( "six x at 0;six x_{out} at 0 [mm];six-plane tracks" );

  sixy0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixy0", 120, -6, 6 );
  sixy0Histo->setTitle( "six y at 0;six y_{up} at 0 [mm];six-plane tracks" );

  sixx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixx1", 240, -12, 12 );
  sixx1Histo->setTitle( "six x at 1;six x_{out} at 1 [mm];six-plane tracks" );

  sixy1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixy1", 120, -6, 6 );
  sixy1Histo->setTitle( "six y at 1;six y_{up} at 1 [mm];six-plane tracks" );

  sixx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixx2", 240, -12, 12 );
  sixx2Histo->setTitle( "six x at 2;six x_{out} at 2 [mm];six-plane tracks" );

  sixy2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixy2", 120, -6, 6 );
  sixy2Histo->setTitle( "six y at 2;six y_{up} at 2 [mm];six-plane tracks" );

  sixx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixx3", 240, -12, 12 );
  sixx3Histo->setTitle( "six x at 3;six x_{out} at 3 [mm];six-plane tracks" );

  sixy3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixy3", 120, -6, 6 );
  sixy3Histo->setTitle( "six y at 3;six y_{up} at 3 [mm];six-plane tracks" );

  sixx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixx4", 240, -12, 12 );
  sixx4Histo->setTitle( "six x at 4;six x_{out} at 4 [mm];six-plane tracks" );

  sixy4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixy4", 120, -6, 6 );
  sixy4Histo->setTitle( "six y at 4;six y_{up} at 4 [mm];six-plane tracks" );

  sixx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixx5", 240, -12, 12 );
  sixx5Histo->setTitle( "six x at 5;six x_{out} at 5 [mm];six-plane tracks" );

  sixy5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixy5", 120, -6, 6 );
  sixy5Histo->setTitle( "six y at 5;six y_{up} at 5 [mm];six-plane tracks" );

  // GBL:
  AIDAProcessor::tree(this)->mkdir("GBL");

  derxtiltHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/derxtilt", 100, -0.1, 0.1 );
  derxtiltHisto->setTitle( "ddx/dtilt;ddx/dtilt [mm/rad];align hits" );

  derytiltHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/derytilt", 100, -10, 10 );
  derytiltHisto->setTitle( "ddy/dtilt;ddy/dtilt [mm/rad];align hits" );

  derxturnHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/derxturn", 100, -10, 10 );
  derxturnHisto->setTitle( "ddx/dturn;ddx/dturn [mm/rad];align hits" );

  deryturnHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/deryturn", 100, -1, 1 );
  deryturnHisto->setTitle( "ddy/dturn;ddy/dturn [mm/rad];align hits" );

  selxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/selx", 240, -12, 12 );
  selxHisto->setTitle( "x at DUT, sel GBL;six x_{out} at DUT [mm];selected tracks" );

  selyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/sely", 120, -6, 6 );
  selyHisto->setTitle( "y at DUT, sel GBL;six y_{up} at DUT [mm];selected tracks" );

  selaxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/selax", 100, -5, 5 );
  selaxHisto->setTitle( "track angle x, sel GBL;x angle [mrad];tracks" );

  selayHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/selay", 100, -5, 5 );
  selayHisto->setTitle( "track angle y, sel GBL;y angle [mrad];tracks" );

  seldxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldx", 100, -150, 150 );
  seldxHisto->setTitle( "track match x, sel GBL;#Deltax [#mum];tracks" );

  seldyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldy", 100, -150, 150 );
  seldyHisto->setTitle( "track match y, sel GBL;#Deltay [#mum];tracks" );

  selkxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/selkx", 100, -10, 10 );
  selkxHisto->setTitle( "kink x, sel GBL;kink x [mrad];tracks" );

  selkyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/selky", 100, -10, 10 );
  selkyHisto->setTitle( "kink y, sel GBL;kink y [mrad];tracks" );

  seldx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldx1", 100, -100, 100 );
  seldx1Histo->setTitle( "triplet resid x at 1, sel GBL;#Deltax [#mum];tracks" );

  seldy1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldy1", 100, -100, 100 );
  seldy1Histo->setTitle( "triplet resid y at 1, sel GBL;#Deltay [#mum];tracks" );

  seldx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldx3", 100, -1000, 1000 );
  seldx3Histo->setTitle( "triplet resid x at 3, sel GBL;#Deltax [#mum];tracks" );

  seldy3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldy3", 100, -1000, 1000 );
  seldy3Histo->setTitle( "triplet resid y at 3, sel GBL;#Deltay [#mum];tracks" );

  seldx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldx4", 100, -1500, 1500 );
  seldx4Histo->setTitle( "triplet resid x at 4, sel GBL;#Deltax [#mum];tracks" );

  seldy4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldy4", 100, -1500, 1500 );
  seldy4Histo->setTitle( "triplet resid y at 4, sel GBL;#Deltay [#mum];tracks" );

  seldx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldx5", 100, -3000, 3000 );
  seldx5Histo->setTitle( "triplet resid x at 5, sel GBL;#Deltax [#mum];tracks" );

  seldy5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldy5", 100, -3000, 3000 );
  seldy5Histo->setTitle( "triplet resid y at 5, sel GBL;#Deltay [#mum];tracks" );

  seldx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldx6", 100, -500, 500 );
  seldx6Histo->setTitle( "triplet resid x at DUT, sel GBL;#Deltax [#mum];tracks" );

  seldy6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/seldy6", 100, -500, 500 );
  seldy6Histo->setTitle( "triplet resid y at DUT, sel GBL;#Deltay [#mum];tracks" );

  gblndfHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblndf", 16, -0.5, 15.5 );
  gblndfHisto->setTitle( "GBL fit NDF;GBL NDF;tracks" );

  gblchi2aHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblchi2a", 100, 0, 100 );
  gblchi2aHisto->setTitle( "GBL fit chi2, DoF 8;GBL chi2;tracks" );

  gblchi2bHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblchi2b", 100, 0, 100 );
  gblchi2bHisto->setTitle( "GBL fit chi2;GBL chi2;tracks" );

  gblprbHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblprb", 100, 0, 1 );
  gblprbHisto->setTitle( "GBL fit probability;GBL fit probability;tracks" );

  // bad fits:

  badxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/badx", 240, -12, 12 );
  badxHisto->setTitle( "x at DUT, bad GBL;six x_{out} at DUT [mm];bad tracks" );

  badyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/bady", 120, -6, 6 );
  badyHisto->setTitle( "y at DUT, bad GBL;six y_{up} at DUT [mm];bad tracks" );

  badaxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/badax", 100, -5, 5 );
  badaxHisto->setTitle( "track angle x, bad GBL;x angle [mrad];tracks" );

  badayHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baday", 100, -5, 5 );
  badayHisto->setTitle( "track angle y, bad GBL;y angle [mrad];tracks" );

  baddxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddx", 100, -150, 150 );
  baddxHisto->setTitle( "track match x, bad GBL;#Deltax [#mum];tracks" );

  baddyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddy", 100, -150, 150 );
  baddyHisto->setTitle( "track match y, bad GBL;#Deltay [#mum];tracks" );

  badkxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/badkx", 100, -10, 10 );
  badkxHisto->setTitle( "kink x, bad GBL;kink x [mrad];tracks" );

  badkyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/badky", 100, -10, 10 );
  badkyHisto->setTitle( "kink y, bad GBL;kink y [mrad];tracks" );

  baddx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddx1", 100, -100, 100 );
  baddx1Histo->setTitle( "triplet resid x at 1, bad GBL;#Deltax [#mum];tracks" );

  baddy1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddy1", 100, -100, 100 );
  baddy1Histo->setTitle( "triplet resid y at 1, bad GBL;#Deltay [#mum];tracks" );

  baddx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddx3", 100, -1000, 1000 );
  baddx3Histo->setTitle( "triplet resid x at 3, bad GBL;#Deltax [#mum];tracks" );

  baddy3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddy3", 100, -1000, 1000 );
  baddy3Histo->setTitle( "triplet resid y at 3, bad GBL;#Deltay [#mum];tracks" );

  baddx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddx4", 100, -1500, 1500 );
  baddx4Histo->setTitle( "triplet resid x at 4, bad GBL;#Deltax [#mum];tracks" );

  baddy4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddy4", 100, -1500, 1500 );
  baddy4Histo->setTitle( "triplet resid y at 4, bad GBL;#Deltay [#mum];tracks" );

  baddx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddx5", 100, -3000, 3000 );
  baddx5Histo->setTitle( "triplet resid x at 5, bad GBL;#Deltax [#mum];tracks" );

  baddy5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddy5", 100, -3000, 3000 );
  baddy5Histo->setTitle( "triplet resid y at 5, bad GBL;#Deltay [#mum];tracks" );

  baddx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddx6", 100, -250, 250 );
  baddx6Histo->setTitle( "triplet resid x at DUT, bad GBL;#Deltax [#mum];tracks" );

  baddy6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/baddy6", 100, -250, 250 );
  baddy6Histo->setTitle( "triplet resid y at DUT, bad GBL;#Deltay [#mum];tracks" );

  // good fits:

  goodx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/goodx1", 100, -100, 100 );
  goodx1Histo->setTitle( "triplet resid x at 1, good GBL;#Deltax [#mum];tracks" );

  goody1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/goody1", 100, -100, 100 );
  goody1Histo->setTitle( "triplet resid y at 1, good GBL;#Deltay [#mum];tracks" );

  goodx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/goodx6", 100, -250, 250 );
  goodx6Histo->setTitle( "triplet resid x at 6, good GBL;#Deltax [#mum];tracks" );

  goody6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/goody6", 100, -250, 250 );
  goody6Histo->setTitle( "triplet resid y at 6, good GBL;#Deltay [#mum];tracks" );

  // look at fit:

  gblax0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax0", 100, -1, 1 );
  gblax0Histo->setTitle( "GBL x angle at plane 0;x angle at plane 0 [mrad];tracks" );

  gbldx0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx0", 100, -10, 10 );
  gbldx0Histo->setTitle( "GBL x shift at plane 0;x shift at plane 0 [#mum];tracks" );

  gblrx0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx0", 100, -25, 25 );
  gblrx0Histo->setTitle( "GBL x resid at plane 0;x resid at plane 0 [#mum];tracks" );

  gblpx0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx0", 100, -10, 10 );
  gblpx0Histo->setTitle( "GBL x pull at plane 0;x pull at plane 0 [#sigma];tracks" );

  gblqx0Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx0", 100, -1, 1 );
  gblqx0Histo->setTitle( "GBL x kink at plane 0;x kink at plane 0 [mrad];tracks" );


  gblax1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax1", 100, -1, 1 );
  gblax1Histo->setTitle( "GBL x angle at plane 1;x angle at plane 1 [mrad];tracks" );

  gbldx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx1", 100, -100, 100 );
  gbldx1Histo->setTitle( "GBL x shift at plane 1;x shift at plane 1 [#mum];tracks" );

  gblrx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx1", 100, -25, 25 );
  gblrx1Histo->setTitle( "GBL x resid at plane 1;x resid at plane 1 [#mum];tracks" );

  gblpx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx1", 100, -10, 10 );
  gblpx1Histo->setTitle( "GBL x pull at plane 1;x pull at plane 1 [#sigma];tracks" );

  gblqx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx1", 100, -1, 1 );
  gblqx1Histo->setTitle( "GBL x kink at plane 1;x kink at plane 1 [mrad];tracks" );

  gblsx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblsx1", 100, -10, 10 );
  gblsx1Histo->setTitle( "GBL x kink at plane 1/error;x kink at plane 1/error;tracks" );

  gbltx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbltx1", 100, -10, 10 );
  gbltx1Histo->setTitle( "GBL x kink pull at plane 1;x kink pull at plane 1;tracks" );


  gblax2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax2", 100, -1, 1 );
  gblax2Histo->setTitle( "GBL x angle at plane 2;x angle at plane 2 [mrad];tracks" );

  gbldx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx2", 100, -20, 20 );
  gbldx2Histo->setTitle( "GBL x shift at plane 2;x shift at plane 2 [#mum];tracks" );

  gblrx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx2", 100, -25, 25 );
  gblrx2Histo->setTitle( "GBL x resid at plane 2;x resid at plane 2 [#mum];tracks" );

  gblpx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx2", 100, -10, 10 );
  gblpx2Histo->setTitle( "GBL x pull at plane 2;x pull at plane 2 [#sigma];tracks" );

  gblqx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx2", 100, -1, 1 );
  gblqx2Histo->setTitle( "GBL x kink at plane 2;x kink at plane 2 [mrad];tracks" );


  gblax3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax3", 100, -1, 1 );
  gblax3Histo->setTitle( "GBL x angle at plane 3;x angle at plane 3 [mrad];tracks" );

  gbldx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx3", 100, -250, 250 );
  gbldx3Histo->setTitle( "GBL x shift at plane 3;x shift at plane 3 [#mum];tracks" );

  gblrx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx3", 100, -25, 25 );
  gblrx3Histo->setTitle( "GBL x resid at plane 3;x resid at plane 3 [#mum];tracks" );

  gblpx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx3", 100, -10, 10 );
  gblpx3Histo->setTitle( "GBL x pull at plane 3;x pull at plane 3 [#sigma];tracks" );

  gblqx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx3", 100, -1, 1 );
  gblqx3Histo->setTitle( "GBL x kink at plane 3;x kink at plane 3 [mrad];tracks" );


  gblax4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax4", 100, -1, 1 );
  gblax4Histo->setTitle( "GBL x angle at plane 4;x angle at plane 4 [mrad];tracks" );

  gbldx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx4", 100, -500, 500 );
  gbldx4Histo->setTitle( "GBL x shift at plane 4;x shift at plane 4 [#mum];tracks" );

  gblrx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx4", 100, -25, 25 );
  gblrx4Histo->setTitle( "GBL x resid at plane 4;x resid at plane 4 [#mum];tracks" );

  gblpx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx4", 100, -10, 10 );
  gblpx4Histo->setTitle( "GBL x pull at plane 4;x pull at plane 4 [#sigma];tracks" );

  gblqx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx4", 100, -1, 1 );
  gblqx4Histo->setTitle( "GBL x kink at plane 4;x kink at plane 4 [mrad];tracks" );


  gblax5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax5", 100, -1, 1 );
  gblax5Histo->setTitle( "GBL x angle at plane 5;x angle at plane 5 [mrad];tracks" );

  gbldx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx5", 100, -1000, 1000 );
  gbldx5Histo->setTitle( "GBL x shift at plane 5;x shift at plane 5 [#mum];tracks" );

  gblrx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx5", 100, -25, 25 );
  gblrx5Histo->setTitle( "GBL x resid at plane 5;x resid at plane 5 [#mum];tracks" );

  gblpx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx5", 100, -10, 10 );
  gblpx5Histo->setTitle( "GBL x pull at plane 5;x pull at plane 5 [#sigma];tracks" );

  gblqx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx5", 100, -1, 1 );
  gblqx5Histo->setTitle( "GBL x kink at plane 5;x kink at plane 5 [mrad];tracks" );


  gblax6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblax6", 100, -1, 1 );
  gblax6Histo->setTitle( "GBL x angle at DUT;x angle at DUT [mrad];tracks" );

  gbldx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldx6", 100, -1000, 1000 );
  gbldx6Histo->setTitle( "GBL x shift at DUT;x shift at DUT [#mum];tracks" );

  gbldy6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbldy6", 100, -1000, 1000 );
  gbldy6Histo->setTitle( "GBL y shift at DUT;y shift at DUT [#mum];tracks" );

  gblrx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblrx6", 100, -250, 250 );
  gblrx6Histo->setTitle( "GBL x resid at DUT;x resid at DUT [#mum];tracks" );

  gblry6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblry6", 100, -100, 100 );
  gblry6Histo->setTitle( "GBL y resid at DUT;y resid at DUT [#mum];tracks" );

  gblpx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpx6", 100, -10, 10 );
  gblpx6Histo->setTitle( "GBL x pull at DUT;x pull at DUT [#sigma];tracks" );

  gblpy6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblpy6", 100, -10, 10 );
  gblpy6Histo->setTitle( "GBL y pull at DUT;y pull at DUT [#sigma];tracks" );

  gblqx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblqx6", 100, -10, 10 );
  gblqx6Histo->setTitle( "GBL x kink at DUT;x kink at DUT [mrad];tracks" );

  gblsx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblsx6", 100, -10, 10 );
  gblsx6Histo->setTitle( "GBL x kink at DUT/error;x kink at DUT/error;tracks" );

  gbltx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gbltx6", 100, -10, 10 );
  gbltx6Histo->setTitle( "GBL x kink pull at DUT;x kink pull at DUT;tracks" );


  gblkx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkx1", 100, -1, 1 );
  gblkx1Histo->setTitle( "GBL kink angle at plane 1;plane 1 kink [mrad];tracks" );

  gblkx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkx2", 100, -1, 1 );
  gblkx2Histo->setTitle( "GBL kink angle at plane 2;plane 2 kink [mrad];tracks" );

  gblkx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkx3", 100, -10, 10 );
  gblkx3Histo->setTitle( "GBL kink angle at plane 3;plane 3 kink [mrad];tracks" );

  gblkx4Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkx4", 100, -1, 1 );
  gblkx4Histo->setTitle( "GBL kink angle at plane 4;plane 4 kink [mrad];tracks" );

  gblkx5Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkx5", 100, -1, 1 );
  gblkx5Histo->setTitle( "GBL kink angle at plane 5;plane 5 kink [mrad];tracks" );

  gblkx6Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "GBL/gblkx6", 100, -1, 1 );
  gblkx6Histo->setTitle( "GBL kink angle at plane 6;plane 6 kink [mrad];tracks" );

  // z intersect:

  sixzx3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixzx3", 100, -50, 250 );
  sixzx3Histo->setTitle( "intersect z-x, kink > 3 mrad;intersect z(x) - z_{2} [mm];tracks" );

  sixzy3Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixzy3", 100, -50, 250 );
  sixzy3Histo->setTitle( "intersect z-y, kink > 3 mrad;intersect z(y) - z_{2} [mm];tracks" );

  sixzx2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixzx2", 100, -50, 250 );
  sixzx2Histo->setTitle( "intersect z-x, kink > 2 mrad;intersect z(x) - z_{2} [mm];tracks" );

  sixzy2Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixzy2", 100, -50, 250 );
  sixzy2Histo->setTitle( "intersect z-y, kink > 2 mrad;intersect z(y) - z_{2} [mm];tracks" );

  sixzx1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixzx1", 100, -50, 250 );
  sixzx1Histo->setTitle( "intersect z-x, kink > 1 mrad;intersect z(x) - z_{2} [mm];tracks" );

  sixzy1Histo = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixzy1", 100, -50, 250 );
  sixzy1Histo->setTitle( "intersect z-y, kink > 1 mrad;intersect z(y) - z_{2} [mm];tracks" );

  sixkxzyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixkxzy", 200, -20, 20 );
  sixkxzyHisto->setTitle( "kink x at DUT;kink x [mrad];tracks" );

  sixkyzxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixkyzx", 200, -20, 20 );
  sixkyzxHisto->setTitle( "kink y at DUT;kink y [mrad];tracks" );

  sixkxzxHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixkxzx", 200, -20, 20 );
  sixkxzxHisto->setTitle( "kink x at DUT;kink x [mrad];tracks" );

  sixkyzyHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Tracks/sixkyzy", 200, -20, 20 );
  sixkyzyHisto->setTitle( "kink y at DUT;kink y [mrad];tracks" );

  // system times:

  cmsdtHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/cmsdt", 500, 0, 5000 );
  cmsdtHisto->setTitle( "DUT time between events;DUT time between events [us];events" );

  dutrefddtHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/dutrefddt", 21, -10.5, 10.5 );
  dutrefddtHisto->setTitle( "DUT - REF #delta#Deltat;DUT - REF #delta#Deltat[clocks];events" );

  sysrtHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/sysrt", 200, 0.9998, 1.0002 );
  sysrtHisto->setTitle( "TLU time / DUT time;event time ratio;events" );

  sysrdtHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/sysrdt", 200, 0.999, 1.001 );
  sysrdtHisto->setTitle( "TLU time / DUT time;time between events ratio;events" );

  dutddtnsHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/dutddtns", 100, -100, 100 );
  dutddtnsHisto->setTitle( "TLU - DUT time;TLU - DUT #delta#Deltat [ns];events" );

  refddtnsHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/refddtns", 100, -100, 100 );
  refddtnsHisto->setTitle( "TLU - REF time;TLU - REF #delta#Deltat [ns];events" );

  dutddtusHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/dutddtus", 100, -100, 100 );
  dutddtusHisto->setTitle( "TLU - DUT time;TLU - DUT #delta#Deltat [us];events" );

  dutddtmsHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/dutddtms", 200, -1, 1 );
  dutddtmsHisto->setTitle( "TLU - DUT time;TLU - DUT #delta#Deltat [ms];events" );

  dutddtvst = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Timing/dutddtvst", 100, 0, 1000, -1E9, 1E9 );
  dutddtvst->setTitle( "TLU - DUT time;time [s];<#delta#Deltat> [ns]" );

  dutddtvsdt = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Timing/dutddtvsdt", 100, 0, 10000, -1E9, 1E9 );
  dutddtvsdt->setTitle( "TLU - DUT time;time between events [us];<#delta#Deltat> [ns]" );

  ddtvst = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Timing/ddtvst", 3000, 0, 300, -99, 99 );
  ddtvst->setTitle( "TLU - DUT time;time [s];<#delta#Deltat> [ns]" );

  ddtvstms = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Timing/ddtvstms", 4000, 8, 12, -99, 99 ); // ms bins
  ddtvstms->setTitle( "TLU - DUT time;time [s];<#delta#Deltat> [ns]" );

  ddtvsdt = AIDAProcessor::histogramFactory(this)->
    createProfile1D( "Timing/ddtvsdt", 200, 0, 2000, -99, 99 );
  ddtvsdt->setTitle( "TLU - DUT time;time between events [us];<#delta#Deltat> [ns]" );

  gapdtHisto = AIDAProcessor::histogramFactory(this)->
    createHistogram1D( "Timing/gapdt", 100, 0, 1 );
  gapdtHisto->setTitle( "time between long clock gaps;time between long clock gaps [s];long clock gaps" );

  // List all booked histogram - check of histogram map filling
  streamlog_out( DEBUG ) <<  _aidaHistoMap.size() << " histograms booked" << endl;


  map<string, AIDA::IBaseHistogram *>::iterator mapIter;
  for(mapIter = _aidaHistoMap.begin(); mapIter != _aidaHistoMap.end(); mapIter++ )
    streamlog_out( DEBUG ) <<  mapIter->first << " : " << ( mapIter->second)->title()  << endl;

  streamlog_out( DEBUG ) << "Histogram booking completed \n\n" << endl;

#else

  streamlog_out( MESSAGE4 ) << "No histogram produced because Marlin doesn't use AIDA" << endl;

#endif

  return;
}

std::string EUTelAnalysisCMSPixel::ZeroPadNumber(int num, int len)
{
    std::ostringstream ss;
    ss << std::setw( len ) << std::setfill( '0' ) << num;
    return ss.str();
}

#endif
