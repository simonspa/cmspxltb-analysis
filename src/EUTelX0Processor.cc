// Version: $Id: EUTelX0Processor.cc 2662 2013-06-02 14:31:25Z hamnett $
#include "EUTelX0Processor.h"
#include <cmath>
#include <cstdlib>
using namespace marlin;
using namespace eutelescope;
using namespace std;

std::string EUTelX0Processor::_histoResidualX = "ResidualX";
std::string EUTelX0Processor::_histoResidualXPlane1 = "ResidualXPlane1";
std::string EUTelX0Processor::_histoResidualXPlane2 = "ResidualXPlane2";
std::string EUTelX0Processor::_histoResidualXPlane3 = "ResidualXPlane3";
std::string EUTelX0Processor::_histoResidualXPlane4 = "ResidualXPlane4";
std::string EUTelX0Processor::_histoResidualYPlane1 = "ResidualYPlane1";
std::string EUTelX0Processor::_histoResidualYPlane2 = "ResidualYPlane2";
std::string EUTelX0Processor::_histoResidualYPlane3 = "ResidualYPlane3";
std::string EUTelX0Processor::_histoResidualYPlane4 = "ResidualYPlane4";
std::string EUTelX0Processor::_histoResidualY = "ResidualY";
std::string EUTelX0Processor::_histoResidualXZ = "ResidualXZ";
std::string EUTelX0Processor::_histoResidualYZ = "ResidualYZ";
std::string EUTelX0Processor::_histoResidualXY = "ResidualYZ";

EUTelX0Processor::EUTelX0Processor()
  :Processor("EUTelX0Processor"),
  _trackColName(""),
  _cutValue1(0.0),
  _cutValue2(0.0),
  _debug(false),
  _debugCount(0),
  _eventNumber(0),
  _finalEvent(false),
  _histoData(),
  _histoFile(""),
  _histoThing(),
  _hitInfo(),
  _inputHitCollectionVec(NULL), 
  _inputTrackCollectionVec(NULL), 
  _inputHitColName(""),
  _inputHitCollectionName(""),
  _inputTrackColName(""),
  _maxRecords(0),
  _projectedHits(),
  _referenceHitCollectionName(""),
  _referenceHitVec(NULL),
  _residual(),
  _residualAngle(),
  _residualCut(0.0),
  _residualProfile(),
  _runNumber(0),
  _trackCollectionName(""),
  nobins(0),
  nobinsangle(0),
  minbin(0),
  maxbin(0),
  minbinangle(0),
  maxbinangle(0),
  minbinalpha(0),
  maxbinalpha(0),
  binsx(0),
  minbinsx(0),//(mm)
  maxbinsx(0),
  binsy(0),
  minbinsy(0),
  maxbinsy(0),
  X0ProcessorDirectory(NULL),
  SinglePointResidualXPlane0(NULL),
  SinglePointResidualXPlane1(NULL),
  SinglePointResidualXPlane2(NULL),
  SinglePointResidualXPlane3(NULL),
  SinglePointResidualXPlane4(NULL),
  SinglePointResidualXPlane5(NULL),
  SinglePointResidualYPlane0(NULL),
  SinglePointResidualYPlane1(NULL),
  SinglePointResidualYPlane2(NULL),
  SinglePointResidualYPlane3(NULL),
  SinglePointResidualYPlane4(NULL),
  SinglePointResidualYPlane5(NULL),
  ThreePointResidualXPlane1(NULL),
  ThreePointResidualXPlane2(NULL),
  ThreePointResidualXPlane3(NULL),
  ThreePointResidualXPlane4(NULL),
  ThreePointResidualYPlane1(NULL),
  ThreePointResidualYPlane2(NULL),
  ThreePointResidualYPlane3(NULL),
  ThreePointResidualYPlane4(NULL),
  AngleXForwardPlane0(NULL),
  AngleXForwardPlane1(NULL),
  AngleXForwardPlane2(NULL),
  AngleXForwardPlane3(NULL),
  AngleXForwardPlane4(NULL),
  AngleYForwardPlane0(NULL),
  AngleYForwardPlane1(NULL),
  AngleYForwardPlane2(NULL),
  AngleYForwardPlane3(NULL),
  AngleYForwardPlane4(NULL),
  ScatteringAngleXPlane1(NULL),
  ScatteringAngleXPlane2(NULL),
  ScatteringAngleXPlane3(NULL),
  ScatteringAngleXPlane4(NULL),
  ScatteringAngleYPlane1(NULL),
  ScatteringAngleYPlane2(NULL),
  ScatteringAngleYPlane3(NULL),
  ScatteringAngleYPlane4(NULL),
  KinkAnglePlane1(NULL),
  KinkAnglePlane2(NULL),
  KinkAnglePlane3(NULL),
  KinkAnglePlane4(NULL),
  ScatteringAngleXPlane1Map(NULL),
  ScatteringAngleXPlane2Map(NULL),
  ScatteringAngleXPlane3Map(NULL),
  ScatteringAngleXPlane4Map(NULL),
  ScatteringAngleYPlane1Map(NULL),
  ScatteringAngleYPlane2Map(NULL),
  ScatteringAngleYPlane3Map(NULL),
  ScatteringAngleYPlane4Map(NULL),
  RadiationLengthPlane1Map(NULL),
  RadiationLengthPlane2Map(NULL),
  RadiationLengthPlane3Map(NULL),
  RadiationLengthPlane4Map(NULL),
  ScatteringAngleXMapData(),
  ScatteringAngleYMapData()
{
  streamlog_out(DEBUG1) << "Constructing the EUTelX0Processor, setting all values to zero or NULL" << std::endl;
  _inputHitCollectionVec = new LCCollectionVec(LCIO::TRACKERHIT);//Used to store the values of the hit events
  _inputTrackCollectionVec = new LCCollectionVec(LCIO::TRACK);//Used to store the values of the hit events
  registerInputCollection(LCIO::TRACKERHIT,"InputTrackCollectionName",
                           "Collection name for corrected particle positions",
                           _trackColName, string ("alignedHit"));

  
  registerInputCollection(LCIO::TRACK,"OutputTrackCollectionName",
                           "Collection name for fitted tracks",
                           _trackCollectionName, string ("AlignedTrack"));

  registerOptionalParameter("ReferenceCollection","This is the name of the reference it collection (init at 0,0,0)", _referenceHitCollectionName, static_cast< string > ( "referenceHit" ) );//Necessary for working out which layer the particle is detected in
  registerProcessorParameter("ResidualCutValue","Used to determine cuts in the system, measured in XXX", _residualCut, static_cast< double > (50000.0));
  registerProcessorParameter("MaxRecords","Will be used to determine the final event if the final event must come before EOF", _maxRecords, static_cast< int > (0));
  registerProcessorParameter("HistoFile","Will be used to add the gaussian to the kink angle at the end of the run", _histoFile, static_cast< std::string > (""));
}

void EUTelX0Processor::init()
{
  streamlog_out(DEBUG5) << "Running EUTelX0Processor::init()" << std::endl;
  nobins = 10000;
  nobinsangle = 1000;//Number of bins in the histograms
  minbin = -0.1;
  maxbin = 0.1;//Maximum and minimum bin values
  minbinangle = -0.01;
  maxbinangle = 0.01;
  minbinalpha = -0.01;
  maxbinalpha = 0.01;
  std::vector<double> empty;  
  binsx = 22;
  minbinsx = -11;//(mm)
  maxbinsx = 11;
  binsy = 12;
  minbinsy = -6;
  maxbinsy = 6;
  double binsizex = (maxbinsx-minbinsx)/binsx;
  double binsizey = (maxbinsy-minbinsy)/binsy;

  X0ProcessorDirectory = new TDirectory("X0Processor","X0Processor");
  X0ProcessorDirectory->cd();

  SinglePointResidualXPlane1 = new TH1D("SinglePointResidualXPlane1",
                                         "Single Point Residual in X on Plane 1;\\Delta X (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["SinglePointResidualXPlane1"] = SinglePointResidualXPlane1;
 
  SinglePointResidualXPlane2 = new TH1D("SinglePointResidualXPlane2",
                                         "Single Point Residual in X on Plane 2;\\Delta X (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["SinglePointResidualXPlane2"] = SinglePointResidualXPlane2;
  
  SinglePointResidualXPlane3 = new TH1D("SinglePointResidualXPlane3",
                                         "Single Point Residual in X on Plane 3;\\Delta X (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["SinglePointResidualXPlane3"] = SinglePointResidualXPlane3;
  
  SinglePointResidualXPlane4 = new TH1D("SinglePointResidualXPlane4",
                                         "Single Point Residual in X on Plane 4;\\Delta X (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["SinglePointResidualXPlane4"] = SinglePointResidualXPlane4;
  
  SinglePointResidualXPlane5 = new TH1D("SinglePointResidualXPlane5",
                                         "Single Point Residual in X on Plane 5;\\Delta X (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["SinglePointResidualXPlane5"] = SinglePointResidualXPlane5;
 
  SinglePointResidualYPlane0 = new TH1D("SinglePointResidualYPlane0",
                                         "Single Point Residual in Y on Plane 0;\\Delta Y (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["SinglePointResidualYPlane0"] = SinglePointResidualYPlane0;
  
  SinglePointResidualYPlane1 = new TH1D("SinglePointResidualYPlane1",
                                         "Single Point Residual in Y on Plane 1;\\Delta Y (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["SinglePointResidualYPlane1"] = SinglePointResidualYPlane1;
  
  SinglePointResidualYPlane2 = new TH1D("SinglePointResidualYPlane2",
                                         "Single Point Residual in Y on Plane 2;\\Delta Y (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["SinglePointResidualYPlane2"] = SinglePointResidualYPlane2;
  
  SinglePointResidualYPlane3 = new TH1D("SinglePointResidualYPlane3",
                                         "Single Point Residual in Y on Plane 3;\\Delta Y (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["SinglePointResidualYPlane3"] = SinglePointResidualYPlane3;
  
  SinglePointResidualYPlane4 = new TH1D("SinglePointResidualYPlane4",
                                         "Single Point Residual in Y on Plane 4;\\Delta Y (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["SinglePointResidualYPlane4"] = SinglePointResidualYPlane4;
  
  SinglePointResidualYPlane5 = new TH1D("SinglePointResidualYPlane5",
                                         "Single Point Residual in Y on Plane 5;\\Delta Y (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["SinglePointResidualYPlane5"] = SinglePointResidualYPlane5;
 
  ThreePointResidualXPlane1 = new TH1D("ThreePointResidualXPlane1",
                                         "Three Point Residual in X on Plane 1;\\Delta X (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["ThreePointResidualXPlane1"] = ThreePointResidualXPlane1;
  
  ThreePointResidualXPlane2 = new TH1D("ThreePointResidualXPlane2",
                                         "Three Point Residual in X on Plane 2;\\Delta X (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["ThreePointResidualXPlane2"] = ThreePointResidualXPlane2;
  
  ThreePointResidualXPlane3 = new TH1D("ThreePointResidualXPlane3",
                                         "Three Point Residual in X on Plane 3;\\Delta X (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["ThreePointResidualXPlane3"] = ThreePointResidualXPlane3;
  
  ThreePointResidualXPlane4 = new TH1D("ThreePointResidualXPlane4",
                                         "Three Point Residual in X on Plane 4;\\Delta X (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["ThreePointResidualXPlane4"] = ThreePointResidualXPlane4;
  
  ThreePointResidualYPlane1 = new TH1D("ThreePointResidualYPlane1",
                                         "Three Point Residual in Y on Plane 1;\\Delta Y (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["ThreePointResidualYPlane1"] = ThreePointResidualYPlane1;
  
  ThreePointResidualYPlane2 = new TH1D("ThreePointResidualYPlane2",
                                         "Three Point Residual in Y on Plane 2;\\Delta Y (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["ThreePointResidualYPlane2"] = ThreePointResidualYPlane2;
  
  ThreePointResidualYPlane3 = new TH1D("ThreePointResidualYPlane3",
                                         "Three Point Residual in Y on Plane 3;\\Delta Y (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["ThreePointResidualYPlane3"] = ThreePointResidualYPlane3;
  
  ThreePointResidualYPlane4 = new TH1D("ThreePointResidualYPlane4",
                                         "Three Point Residual in Y on Plane 4;\\Delta Y (mm);Count",
					 nobins,minbin,maxbin);
  _histoThing["ThreePointResidualYPlane4"] = ThreePointResidualYPlane4;
  
  AngleXForwardPlane0 = new TH1D("AngleXForwardPlane0",
                                 "Angle of Tracks in X Direction Relative to the Z Axis on Plane 0;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleXForwardPlane0"] = AngleXForwardPlane0;
 
  AngleXForwardPlane1 = new TH1D("AngleXForwardPlane1",
                                 "Angle of Tracks in X Direction Relative to the Z Axis on Plane 1;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleXForwardPlane1"] = AngleXForwardPlane1;
 
  AngleXForwardPlane2 = new TH1D("AngleXForwardPlane2",
                                 "Angle of Tracks in X Direction Relative to the Z Axis on Plane 2;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleXForwardPlane2"] = AngleXForwardPlane2;
 
  AngleXForwardPlane3 = new TH1D("AngleXForwardPlane3",
                                 "Angle of Tracks in X Direction Relative to the Z Axis on Plane 3;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleXForwardPlane3"] = AngleXForwardPlane3;
 
  AngleXForwardPlane4 = new TH1D("AngleXForwardPlane4",
                                 "Angle of Tracks in X Direction Relative to the Z Axis on Plane 4;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleXForwardPlane4"] = AngleXForwardPlane4;
  
  AngleYForwardPlane0 = new TH1D("AngleYForwardPlane0",
                                 "Angle of Tracks in Y Direction Relative to the Z Axis on Plane 0;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleYForwardPlane0"] = AngleYForwardPlane0;
 
  AngleYForwardPlane1 = new TH1D("AngleYForwardPlane1",
                                 "Angle of Tracks in Y Direction Relative to the Z Axis on Plane 1;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleYForwardPlane1"] = AngleYForwardPlane1;
 
  AngleYForwardPlane2 = new TH1D("AngleYForwardPlane2",
                                 "Angle of Tracks in Y Direction Relative to the Z Axis on Plane 2;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleYForwardPlane2"] = AngleYForwardPlane2;
 
  AngleYForwardPlane3 = new TH1D("AngleYForwardPlane3",
                                 "Angle of Tracks in Y Direction Relative to the Z Axis on Plane 3;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleYForwardPlane3"] = AngleYForwardPlane3;
 
  AngleYForwardPlane4 = new TH1D("AngleYForwardPlane4",
                                 "Angle of Tracks in Y Direction Relative to the Z Axis on Plane 4;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["AngleYForwardPlane4"] = AngleYForwardPlane4;
  
  ScatteringAngleXPlane1 = new TH1D("ScatteringAngleXPlane1",
                                    "Scattering Angle in X Direction Relative to the Z Axis on Plane 1;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["ScatteringAngleXPlane1"] = ScatteringAngleXPlane1;
 
  ScatteringAngleXPlane2 = new TH1D("ScatteringAngleXPlane2",
                                    "Scattering Angle in X Direction Relative to the Z Axis on Plane 2;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["ScatteringAngleXPlane2"] = ScatteringAngleXPlane2;
 
  ScatteringAngleXPlane3 = new TH1D("ScatteringAngleXPlane3",
                                    "Scattering Angle in X Direction Relative to the Z Axis on Plane 3;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["ScatteringAngleXPlane3"] = ScatteringAngleXPlane3;
 
  ScatteringAngleXPlane4 = new TH1D("ScatteringAngleXPlane4",
                                    "Scattering Angle in X Direction Relative to the Z Axis on Plane 4;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["ScatteringAngleXPlane4"] = ScatteringAngleXPlane4;
  
  ScatteringAngleYPlane1 = new TH1D("ScatteringAngleYPlane1",
                                    "Scattering Angle in Y Direction Relative to the Z Axis on Plane 1;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["ScatteringAngleYPlane1"] = ScatteringAngleYPlane1;
 
  ScatteringAngleYPlane2 = new TH1D("ScatteringAngleYPlane2",
                                    "Scattering Angle in Y Direction Relative to the Z Axis on Plane 2;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["ScatteringAngleYPlane2"] = ScatteringAngleYPlane2;
 
  ScatteringAngleYPlane3 = new TH1D("ScatteringAngleYPlane3",
                                    "Scattering Angle in Y Direction Relative to the Z Axis on Plane 3;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["ScatteringAngleYPlane3"] = ScatteringAngleYPlane3;
 
  ScatteringAngleYPlane4 = new TH1D("ScatteringAngleYPlane4",
                                    "Scattering Angle in Y Direction Relative to the Z Axis on Plane 4;\\theta_x (rads);Count",
				 nobinsangle,minbinangle,maxbinangle);
  _histoThing["ScatteringAngleYPlane4"] = ScatteringAngleYPlane4;
     
  KinkAnglePlane1 = new TH2D("KinkAnglePlane1",
                             "Kink Angle on Plane 1;\\theta_x (rads);\\theta_y (rads)",
			     nobinsangle,minbinalpha,maxbinalpha,nobinsangle,minbinalpha,maxbinalpha);
  _histoThing["KinkAnglePlane1"] = KinkAnglePlane1;
  
  KinkAnglePlane2 = new TH2D("KinkAnglePlane2",
                             "Kink Angle on Plane 2;\\theta_x (rads);\\theta_y (rads)",
			     nobinsangle,minbinalpha,maxbinalpha,nobinsangle,minbinalpha,maxbinalpha);
  _histoThing["KinkAnglePlane2"] = KinkAnglePlane2;
  
  KinkAnglePlane3 = new TH2D("KinkAnglePlane3",
                             "Kink Angle on Plane 3;\\theta_x (rads);\\theta_y (rads)",
			     nobinsangle,minbinalpha,maxbinalpha,nobinsangle,minbinalpha,maxbinalpha);
  _histoThing["KinkAnglePlane3"] = KinkAnglePlane3;
  
  KinkAnglePlane4 = new TH2D("KinkAnglePlane4",
                             "Kink Angle on Plane 4;\\theta_x (rads);\\theta_y (rads)",
			     nobinsangle,minbinalpha,maxbinalpha,nobinsangle,minbinalpha,maxbinalpha);
  _histoThing["KinkAnglePlane4"] = KinkAnglePlane4;
  
  ScatteringAngleXPlane1Map = new TH2D("ScatteringAngleXPlane1Map",
                                       "Scattering Angle in X on Plane 1;x (mm);y(mm)"
				       ,binsx,minbinsx,maxbinsx,binsy,minbinsy,maxbinsy);
  _histoThing["ScatteringAngleXPlane1Map"] = ScatteringAngleXPlane1Map;
  
  ScatteringAngleXPlane2Map = new TH2D("ScatteringAngleXPlane2Map",
                                       "Scattering Angle in X on Plane 2;x (mm);y(mm)"
				       ,binsx,minbinsx,maxbinsx,binsy,minbinsy,maxbinsy);
  _histoThing["ScatteringAngleXPlane2Map"] = ScatteringAngleXPlane2Map;
  
  ScatteringAngleXPlane3Map = new TH2D("ScatteringAngleXPlane3Map",
                                       "Scattering Angle in X on Plane 3;x (mm);y(mm)"
				       ,binsx,minbinsx,maxbinsx,binsy,minbinsy,maxbinsy);
  _histoThing["ScatteringAngleXPlane3Map"] = ScatteringAngleXPlane3Map;
  
  ScatteringAngleXPlane4Map = new TH2D("ScatteringAngleXPlane4Map",
                                       "Scattering Angle in X on Plane 4;x (mm);y(mm)"
				       ,binsx,minbinsx,maxbinsx,binsy,minbinsy,maxbinsy);
  _histoThing["ScatteringAngleXPlane4Map"] = ScatteringAngleXPlane4Map;
  
  ScatteringAngleYPlane1Map = new TH2D("ScatteringAngleYPlane1Map",
                                       "Scattering Angle in Y on Plane 1;x (mm);y(mm)"
				       ,binsx,minbinsx,maxbinsx,binsy,minbinsy,maxbinsy);
  _histoThing["ScatteringAngleYPlane1Map"] = ScatteringAngleYPlane1Map;
  
  ScatteringAngleYPlane2Map = new TH2D("ScatteringAngleYPlane2Map",
                                       "Scattering Angle in Y on Plane 2;x (mm);y(mm)"
				       ,binsx,minbinsx,maxbinsx,binsy,minbinsy,maxbinsy);
  _histoThing["ScatteringAngleYPlane2Map"] = ScatteringAngleYPlane2Map;
  
  ScatteringAngleYPlane3Map = new TH2D("ScatteringAngleYPlane3Map",
                                       "Scattering Angle in Y on Plane 3;x (mm);y(mm)"
				       ,binsx,minbinsx,maxbinsx,binsy,minbinsy,maxbinsy);
  _histoThing["ScatteringAngleYPlane3Map"] = ScatteringAngleYPlane3Map;
  
  ScatteringAngleYPlane4Map = new TH2D("ScatteringAngleYPlane4Map",
                                       "Scattering Angle in Y on Plane 4;x (mm);y(mm)"
				       ,binsx,minbinsx,maxbinsx,binsy,minbinsy,maxbinsy);
  _histoThing["ScatteringAngleYPlane4Map"] = ScatteringAngleYPlane4Map;
  
  RadiationLengthPlane1Map = new TH2D("RadiationLengthPlane1Map",
                                       "Radiation Length on Plane 1;x (mm);y(mm)"
				       ,binsx,minbinsx,maxbinsx,binsy,minbinsy,maxbinsy);
  _histoThing["RadiationLengthPlane1Map"] = RadiationLengthPlane1Map;
  
  RadiationLengthPlane2Map = new TH2D("RadiationLengthPlane2Map",
                                       "Radiation Length on Plane 2;x (mm);y(mm)"
				       ,binsx,minbinsx,maxbinsx,binsy,minbinsy,maxbinsy);
  _histoThing["RadiationLengthPlane2Map"] = RadiationLengthPlane2Map;
  
  RadiationLengthPlane3Map = new TH2D("RadiationLengthPlane3Map",
                                       "Radiation Length on Plane 3;x (mm);y(mm)"
				       ,binsx,minbinsx,maxbinsx,binsy,minbinsy,maxbinsy);
  _histoThing["RadiationLengthPlane3Map"] = RadiationLengthPlane3Map;
  
  RadiationLengthPlane4Map = new TH2D("RadiationLengthPlane4Map",
                                       "Radiation Length on Plane 4;x (mm);y(mm)"
				       ,binsx,minbinsx,maxbinsx,binsy,minbinsy,maxbinsy);
  _histoThing["RadiationLengthPlane4Map"] = RadiationLengthPlane4Map;
}

void EUTelX0Processor::processRunHeader(LCRunHeader *run)
{
  streamlog_out(DEBUG1) << "Running EUTelX0Processor::processRunHeader(LCRunHeader *run) with run = " << run << std::endl;
  _eventNumber = 0;
  _runNumber++;
}

void EUTelX0Processor::printTrackParameters( Track* eventtrack ){
    streamlog_out(DEBUG0) << "Type: " << eventtrack->getType() << std::endl << std::endl
    << "***TRACK PARAMETERS***" << std::endl << "(Note: In EUTelTestFitter the track parameters are not yet instantiated, all currently set to zero)" << std::endl
    << "D0: " << eventtrack->getD0() << std::endl
    << "Omega: " << eventtrack->getOmega() << std::endl
    << "Phi: " << eventtrack->getPhi() << std::endl
    << "Z0: " << eventtrack->getZ0() << std::endl
    << "Tan Lambda: " << eventtrack->getTanLambda() << std::endl 
    << "Covariance Matrix of Track Parameters: (This is a square matrix with parameter order: D0, Phi, Omega, Z0, TanLambda)" << std::endl;
    std::vector< float > covmatrix = eventtrack->getCovMatrix();
    streamlog_out(DEBUG0) << covmatrix[0] << ", " << covmatrix[1] << ", " << covmatrix[3] << ", " << covmatrix[6] << ", " << covmatrix[10] << std::endl
                          << covmatrix[1] << ", " << covmatrix[2] << ", " << covmatrix[4] << ", " << covmatrix[7] << ", " << covmatrix[11] << std::endl
                          << covmatrix[3] << ", " << covmatrix[4] << ", " << covmatrix[5] << ", " << covmatrix[8] << ", " << covmatrix[12] << std::endl
                          << covmatrix[6] << ", " << covmatrix[7] << ", " << covmatrix[8] << ", " << covmatrix[9] << ", " << covmatrix[13] << std::endl
                          << covmatrix[10] << ", " << covmatrix[11] << ", " << covmatrix[12] << ", " << covmatrix[13] << ", " << covmatrix[14] << std::endl;    
    streamlog_out(DEBUG0) << std::endl
    << "***ADDITIONAL INFORMATION***" << std::endl
    << "Reference Point: " << eventtrack->getReferencePoint()[0] << ", " <<  eventtrack->getReferencePoint()[1] << ", " << eventtrack->getReferencePoint()[2] << std::endl
    << "Chi2: " << eventtrack->getChi2() << std::endl
    << "Degrees of Freedom: " << eventtrack->getNdf() << std::endl
    << "dE/dx: " << eventtrack->getdEdx() << std::endl
    << "dE/dx error: " << eventtrack->getdEdxError() << std::endl << std::endl;
  streamlog_out(DEBUG0) << "***TRACK HITS***" << std::endl;
  std::vector< TrackerHit* > trackhits = eventtrack->getTrackerHits();
  int tracknumber(0);
  for(std::vector< TrackerHit* >::iterator it = trackhits.begin(); it != trackhits.end(); ++it){
    streamlog_out(DEBUG0) << "Track Hit " << tracknumber << ":" << std::endl;
    tracknumber++;
    printHitParameters(*it);
  }
}

void EUTelX0Processor::printHitParameters( TrackerHit* trackhit ){
  streamlog_out(DEBUG0) << "Type: " << trackhit->getType() << std::endl
    << "Time of hit (in nanoseconds): " << trackhit->getTime() << std::endl
    << "dE/dx (in GeV): " << trackhit->getEDep() << std::endl
    << "Position (x, y, z and measured in mm): " << trackhit->getPosition()[0] << ", "
    << trackhit->getPosition()[1] << ", " << trackhit->getPosition()[2] << std::endl
    << "Covariance of the position (x, y, z): " << std::endl;
  std::vector< float > hitcov = trackhit->getCovMatrix(); //Size 6
  streamlog_out(DEBUG0) << hitcov[0] << ", " << hitcov[1] << ", " << hitcov[3] << std::endl
                        << hitcov[1] << ", " << hitcov[2] << ", " << hitcov[4] << std::endl
                        << hitcov[3] << ", " << hitcov[4] << ", " << hitcov[5] << std::endl << std::endl;
}

void EUTelX0Processor::kinkEstimate(Track* track){
  //This function works out an angle based on a straight line fitted from plane 0 to 2 and plane 5 to 3
  //It will also store all other angles in histograms too
  streamlog_out(DEBUG0) << "Running function kinkEstimate(Track* " << &track << ")" << std::endl;
  
  //First we extract the relevant hits from the track
  try{
    std::vector< TVector3* > hits = getHitsFromTrack(track);
  } catch(...)
  { 
    streamlog_out(ERROR3) << "Could not get hits from track" << endl;
    return;
  }
  std::vector< TVector3* > hits = getHitsFromTrack(track);

  streamlog_out(DEBUG3) << "Successfully got hits from track" << std::endl;
  //Then we find the position where those lines would have hit the DUT
  //THIS IS TO BE DECIDED IF IT IS NEEDED LATER
  const double dutposition(hits[2]->z());
  streamlog_out(DEBUG3) << "For now DUT position is set to the position of the 2nd plane, that is at " << dutposition << endl;
  //Then we work out the angles of these lines with respect to XZ and YZ, plot results in histograms
  const size_t hitsize = hits.size();
  if(hitsize <= 1)
  {
    streamlog_out(DEBUG5) << "This track has only one hit, aborting track" << endl;
    return;
  }
  std::vector< double > scatterx, scattery;
  streamlog_out(DEBUG3) << "hitsize = " << hitsize << std::endl;
  for(size_t i = 0; i < hitsize-1; ++i){
    streamlog_out(DEBUG3) << "Running now in the first for loop in element " << i << std::endl;
    const double x0 = hits[i]->x();
    const double y0 = hits[i]->y();
    const double z0 = hits[i]->z();
    const double x1 = hits[i+1]->x();
    const double y1 = hits[i+1]->y();
    const double z1 = hits[i+1]->z();
    const double deltaz = z1-z0;
    streamlog_out(DEBUG3) << x0 << endl << y0 << endl << z0 << endl << x1 << endl << y1 << endl << z1 << endl;
    const double frontanglex = atan2(x1-x0,deltaz);
    const double frontangley = atan2(y1-y0,deltaz);
    scatterx.push_back(frontanglex);
    scattery.push_back(frontangley);
    std::stringstream xforward,yforward,xbackward,ybackward;
    xforward << "AngleXForwardPlane" << i;
    yforward << "AngleYForwardPlane" << i;
    try{
      dynamic_cast< TH1D* >(_histoThing[xforward.str().c_str()])->Fill(frontanglex);
    } catch(std::bad_cast &bc){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << xforward.str().c_str() << endl;
    }
    try{
      dynamic_cast< TH1D* >(_histoThing[yforward.str().c_str()])->Fill(frontangley);
    } catch(std::bad_cast &bc){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << yforward.str().c_str() << endl;
    }
  }
 
  const size_t scatterxsize = scatterx.size();
  streamlog_out(DEBUG3) << "scatterxsize = " << scatterxsize << std::endl;
  for(size_t i = 0; i < scatterxsize-1; ++i){
    streamlog_out(DEBUG3) << "Running now in the second for loop in element " << i << std::endl;
    const double scatteringanglex = scatterx[i+1] - scatterx[i];
    const double scatteringangley = scattery[i+1] - scattery[i];
    std::stringstream ssscatterx,ssscattery,ssscatterxy;
    ssscatterx << "ScatteringAngleXPlane" << i+1;
    ssscattery << "ScatteringAngleYPlane" << i+1;
    ssscatterxy << "KinkAnglePlane" << i+1;
    if(i == 1)
    {
      double binsizex = (maxbinsx - minbinsx)/binsx;
      double binsizey = (maxbinsx - minbinsx)/binsx;
      double x = hits[i+1]->x();
      double y = hits[i+1]->y();
      double newminbinx(minbinsx);
      double newminbiny(minbinsy);
      int actualbinx(0);
      int actualbiny(0);
      bool notfoundx = true;
      while(notfoundx)
      {
        if(x >= newminbinx && x < newminbinx + binsizex)
        {
          notfoundx = false;
          bool notfoundy = true;
          while(notfoundy)
          {
            if(y >= newminbiny && y < newminbiny + binsizey)
	    {
	      //track is in bin actualbinx, actualbiny
	      notfoundy = false;
              pair< int, int > position(actualbinx,actualbiny);
              ScatteringAngleXMapData[position].push_back(scatteringanglex);
	    } else if(newminbiny <= maxbinsy)
    	    {
	      newminbiny += binsizey;
	      ++actualbiny;
	    } else
	    {
	      streamlog_out(ERROR5) << "Somehow there is a hit beyond the scope of the sensor" << endl;
	    }
	  }
        }
        else if(newminbinx <= maxbinsx)
        {
          newminbinx += binsizex;
	  ++actualbinx;
        } else
	{
	  streamlog_out(ERROR5) << "Somehow there is a hit beyond the scope of the sensor" << endl;
	}
      }
    }
    try{
      dynamic_cast< TH1D* >(_histoThing[ssscatterx.str().c_str()])->Fill(scatteringanglex);
    } catch(std::bad_cast &bc){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << ssscatterx.str().c_str() << endl;
    }
    try{
      dynamic_cast< TH1D* >(_histoThing[ssscattery.str().c_str()])->Fill(scatteringangley);
    } catch(std::bad_cast &bc){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << ssscattery.str().c_str() << endl;
    }
    try{
      dynamic_cast< TH2D* >(_histoThing[ssscatterxy.str().c_str()])->Fill(scatteringanglex,scatteringangley);
    } catch(std::bad_cast &bc){
      streamlog_out(ERROR3) << "Unable to fill histogram: " << ssscatterxy.str().c_str() << endl;
    }
  }
  streamlog_out(DEBUG3) << "Made it to the end of kinkEstimate()" << endl;
}

void EUTelX0Processor::kinkGaussian(){
/*
  gStyle->SetOptFit(111);
  std::string _histoFile = "/scratch/hamnett/TestBeam/2013/data_X0MeasurementsOneWeek20012013/histo/000001-track-histo.root";
  TFile *file = new TFile(_histoFile.c_str(),"UPDATE");
  TDirectory *directory = (TDirectory*)file->GetDirectory("X0Scanner");
  directory->ls();
  TH1D *histogram = (TH1D*)directory->Get("KinkAngle");
  histogram->Draw();
  TF1 *fit = new TF1("fit","gaus",0,5);
  histogram->Fit("fit","R");
  file->Write();
*/
  
//  AIDA::IFitter * fit = AIDAProcessor::getIAnalysisFactory(this)->createFitFactory()->createFitter("Chi2")->fit(_histoThing["Kink Angle"],fitfunction);

}

void EUTelX0Processor::processEvent(LCEvent *evt)
{
  streamlog_out(DEBUG0) << "Running EUTelX0Processor::processEvent(LCEvent *evt) with evt = " << evt << std::endl;
  //Take track from input parameter
  //Work out kink angle from track
  //Put kink angle into histogram
  //Check for last event
    //Fit gaussian to histogram
    //Extract sigma value from histogram
    //Deduce radiation length from sigma value - See Nadler and Fruwurth paper

  try{
    if(_eventNumber == _maxRecords - 2){ //Not sure about the -2
      _finalEvent = true;
      kinkGaussian();
    } else{
      //_referenceHitVec = evt->getCollection(_referenceHitCollectionName);
      LCCollection* trackcollection = evt->getCollection(_trackCollectionName);
      int elementnumber = trackcollection->getNumberOfElements();
      for(int i = 0; i < elementnumber; ++i){
        Track* eventtrack = dynamic_cast< Track* >(trackcollection->getElementAt(i));
        streamlog_out(DEBUG0) << "Here is all the information about the track in run " << _runNumber << ", event " << _eventNumber << ", element " << i << std::endl << std::endl;
        printTrackParameters( eventtrack );
        kinkEstimate( eventtrack );
//        singlePointResolution( eventtrack );
//        threePointResolution( eventtrack );
      }
    }
  } catch(DataNotAvailableException &datanotavailable){
    streamlog_out(DEBUG4) << "Exception occured: " << datanotavailable.what() << std::endl
      << "Could not get collection '" << _trackCollectionName << "' from event " << _eventNumber << ", Skipping Event" << std::endl;
  } catch(...){
    streamlog_out(ERROR9) << "Unknown exception occured in EUTelX0Processor, here is some information which might help:" << std::endl
    << "Run Number: " << _runNumber << std::endl
    << "Event Number: " << _eventNumber << std::endl
    << "Throw occured somewhere within the try block which starts with the line: 'LCCollection* trackcollection = evt->getCollection(_trackCollectionName);'" << std::endl
    << "I hope that helps, if not then please submit a bug report or add more exception handling as appropriate, program will now exit" << std::endl;
    exit(1);
  }
  _eventNumber++;
}

int EUTelX0Processor::guessSensorID(const double * hit )
{
  int sensorID = -1;
  double minDistance = numeric_limits< double >::max() ;
  streamlog_out( DEBUG5 ) << "referencehit collection: " << _referenceHitCollectionName << " at "<< _referenceHitVec << std::endl;
  if( _referenceHitVec == 0)
  {
    streamlog_out( DEBUG5 ) << "_referenceHitVec is empty" << std::endl;
    return 0;
  }

  if( isFirstEvent() )
  {
    // message<DEBUG > ( log() << "number of elements : " << _referenceHitVec->getNumberOfElements() << std::endl );
  }

  for(int ii = 0 ; ii < _referenceHitVec->getNumberOfElements(); ii++)
  {
    EUTelReferenceHit* refhit = static_cast< EUTelReferenceHit*> ( _referenceHitVec->getElementAt(ii) ) ;
    if(refhit == 0 ) continue;
     TVector3 hit3d( hit[0], hit[1], hit[2] );
    TVector3 hitInPlane( refhit->getXOffset(), refhit->getYOffset(), refhit->getZOffset());
    TVector3 norm2Plane( refhit->getAlpha(), refhit->getBeta(), refhit->getGamma() );
    double distance = abs( norm2Plane.Dot(hit3d-hitInPlane) );
    streamlog_out( DEBUG5 ) << " hit " << hit[0] << " "<< hit[1] << " " << hit[2] << std::endl;
    streamlog_out( DEBUG5 ) << " " << refhit->getXOffset() << " " << refhit->getYOffset() << " " << refhit->getZOffset() << std::endl;
    streamlog_out( DEBUG5 ) << " " << refhit->getAlpha() << " " << refhit->getBeta() << " " << refhit->getGamma() << std::endl;
    streamlog_out( DEBUG5 ) << " distance " << distance << std::endl;
    if ( distance < minDistance )
    {
      minDistance = distance;
      sensorID = refhit->getSensorID();
    }
  }
  //Some useful debug printouts:
  bool debug = ( _debugCount>0 );

  if(debug)
  {
    for(int ii = 0 ; ii < _referenceHitVec->getNumberOfElements(); ii++)
    {
      EUTelReferenceHit* refhit = static_cast< EUTelReferenceHit*> ( _referenceHitVec->getElementAt(ii) ) ;
      if(refhit == 0 ) continue;
      // if( sensorID != refhit->getSensorID() ) continue;
      // streamlog_out( DEBUG5 ) << " _referenceHitVec " << _referenceHitVec << " " << _referenceHitCollectionName.c_str() << " " << refhit << " at "
      // << refhit->getXOffset() << " " << refhit->getYOffset() << " " << refhit->getZOffset() << " "
      //<< refhit->getAlpha() << " " << refhit->getBeta() << " " << refhit->getGamma() << std::endl ;
      //message<DEBUG5> ( log() << "iPlane " << refhit->getSensorID() << " hitPos: [" << hit[0] << " " << hit[1] << " " << hit[2] << "] distance: " << minDistance << std::endl );
    }
  }
  return sensorID;
}
 
void EUTelX0Processor::end()
{
  //Clean up memory
  //Set all values to zero or NULL
  //
  //calculateX0();
  double binsizex = (maxbinsx-minbinsx)/binsx;
  double binsizey = (maxbinsy-minbinsy)/binsy;

  for(double i = minbinsx; i <= maxbinsx; i += binsizex)
  {
    for(double j = minbinsy; j <= maxbinsy; j += binsizey)
    {
      
      //cout << "x = " << i->first.first << ", y = " << i->first.second << ", angle = " << i->second << endl;
    }
  }
  _hitInfo.clear();
  _projectedHits.clear();
  _histoThing.clear(); 
  _residual.clear(); 
  _residualAngle.clear();
  _residualProfile.clear();
  _inputHitCollectionVec = NULL;
  _inputTrackCollectionVec = NULL; 
  _referenceHitVec = NULL;  
}

double EUTelX0Processor::calculateX0()
{
  double X0(0);
  //Two methods for figuring out the material budget, X0, are used. The loss of energy and the scattering angle due to multiple scattering.
  //Potentially might use both a Forward-backward Kalmin filter and a Global Linear Estimator to work out the two above values, might give a more accurate meaurement of X0.
  double sigma_msi(0);  //Highland formula
  double sigma_measi(0);  //Sum of the variance of histogram alpha012 and alpha432
  double E = 5;  //The total particle energy (GeV) TODO: (Phillip Hamnett) is there some way to make this information (and the mass) into the steering file?
  double m = 0.000511;  //Electron Mass (GeV)
  double p = sqrt(E*E - m*m);  //The particle momentum (GeV)
  double theta012(0);  //The polar angle between the momentum vector and the z-axis
  double beta012(0);  //The difference between the azimuthal angle of the momentum vector and the azimuthal angle of a point on the track
  double theta432(0);
  double beta432(0);
  X0 = 100;  //An estimate of X0
  double thickness = 0.00005;  //Thickness of mimosa layer is 50 micrometers
  double X = X0/thickness;
  double BetterX = DBL_MAX;
  double currentlikelihood = DBL_MIN;
  double lastX0 = DBL_MIN;
  while(X0 != BetterX){
    sigma_measi = 5;  //HACK - This is actually the sum of the variances of theta in the forward and backward direction. 
    size_t size = _histoData["Theta012"].size();
    double loglikelihoodX0(0);
    for(size_t i = 0; i < size; ++i){
      theta012 = _histoData["Theta012"][i]*3.1415/180.0;
      theta432 = _histoData["Theta432"][i]*3.1415/180.0;
      beta012 = _histoData["Phi012"][i]*3.1415/180.0;
      beta432 = _histoData["Phi432"][i]*3.1415/180.0;
      double deltatheta = sqrt((theta432 - theta012)*(theta432 - theta012));  //Take the positive difference between the two angles in theta and beta
      double deltabeta = sqrt((beta432 - beta012)*(beta432 - beta012));
      sigma_msi = (0.015*E/(p*p))*sqrt(X/(sin(deltatheta)*cos(deltabeta)))*(1 + 0.038*std::log(X/(sin(deltatheta)*cos(deltabeta))));
      if(sigma_msi != sigma_msi){
        streamlog_out( ERROR5 ) << "Code failed with following values: " << std::endl << "E = " << E << std::endl << "p = " << p << std::endl << "X = " << X << std::endl << "theta012 = " << theta012 << std::endl << "beta012 = " << beta012 << std::endl << "theta432 = " << theta432 << std::endl << "beta432 = " << beta432 << std::endl << "deltatheta = " << deltatheta << std::endl << "deltabeta = " << deltabeta << std::endl << "sigma_msi = " << sigma_msi << std::endl << "sigma_measi = " << sigma_measi << std::endl << "loglikelihoodX0 = " << loglikelihoodX0 << std::endl << std::endl;
      }
      //Maximum log-likelihood:
      loglikelihoodX0 += -std::log(sigma_msi*sigma_msi + sigma_measi*sigma_measi) - (theta012*theta012)/(sigma_msi*sigma_msi + sigma_measi*sigma_measi); 
    }
    if(loglikelihoodX0 > currentlikelihood){
      BetterX = X0;
      currentlikelihood = loglikelihoodX0;
      if(lastX0 < X0){
        X++;
      }
      else{
        X--;
      }
    }
    else{
      if(lastX0 < X0){
        X--;
      }
      else{
        X++;
      }
    }
    X = X0/thickness;
    loglikelihoodX0 = 0;
  }
  return X0;
}

std::vector< TVector3* > EUTelX0Processor::getHitsFromTrack(Track *track){
 std::vector< TrackerHit* > trackhits = track->getTrackerHits();
 std::vector< TVector3* > hits;
 for(std::vector< TrackerHit* >::iterator it = trackhits.begin(); it != trackhits.end(); ++it){
    if((*it)->getType() < 32){  //Check if the hit type is appropriate
      Double_t x = (*it)->getPosition()[0];
      Double_t y = (*it)->getPosition()[1];
      Double_t z = (*it)->getPosition()[2];
      TVector3 *tempvec = new TVector3(x,y,z);
      hits.push_back(tempvec);
    } //End of if query for type check
  } //End of for loop running through trackhits
  return hits;
}

void EUTelX0Processor::singlePointResolution(Track *track){
//This function draws a line between two hits on adjacent planes and then projects in the forward and backward direction to the hit on the next plane.
//It then compares that hit with the actual hit on that next plane and the differencce is plotted
  streamlog_out(DEBUG1) << "Function EUTeoX0Processor::singlePointResolution(Track *" << &track << ") called" << std::endl;
  std::vector< TVector3* > hits = getHitsFromTrack(track);
  int i = 0;
  for(std::vector< TVector3* >::iterator it = hits.begin(); it != hits.end(); ++it){
    streamlog_out(DEBUG0) << "Entering for loop for Plane " << i << endl;
    stringstream histogramnamex,histogramnamey;
    histogramnamex << "SinglePointResidualXPlane" << i;
    histogramnamey << "SinglePointResidualYPlane" << i;
    double x0, y0, x1, y1, x2, y2, z0, z1, z2;
    x0 = (*it)->x();
    y0 = (*it)->y();
    z0 = (*it)->z();
    x1 = (*it + 1)->x();
    y1 = (*it + 1)->y();
    z1 = (*it + 1)->z();
    x2 = (*it + 2)->x();
    y2 = (*it + 2)->y();
    z2 = (*it + 2)->z();
    double predictx, predicty;
    if(it == hits.begin() || it == hits.begin() + 1){
      predictx = (x2-x1)*(z2-z0)/(z2-z1);
      predicty = (y2-y1)*(z2-z0)/(z2-z1);
      try{
        dynamic_cast< TH1D* > (_histoThing[histogramnamex.str().c_str()])->Fill(x0 - predictx);
      } catch(std::bad_cast &bc){
        streamlog_out(ERROR3) << "Unable to fill histogram: " << histogramnamex.str().c_str() << endl;
      }
      try{
        dynamic_cast< TH1D* > (_histoThing[histogramnamey.str().c_str()])->Fill(y0 - predicty);
      } catch(std::bad_cast &bc){
        streamlog_out(ERROR3) << "Unable to fill histogram: " << histogramnamey.str().c_str() << endl;
      }
    } else {
      predictx = (x1-x0)*(z2-z0)/(z1-z0);
      predicty = (y1-y0)*(z2-z0)/(z1-z0);
       try{
        dynamic_cast< TH1D* > (_histoThing[histogramnamex.str().c_str()])->Fill(x2 - predictx);
      } catch(std::bad_cast &bc){
        streamlog_out(ERROR3) << "Unable to fill histogram: " << histogramnamex.str().c_str() << endl;
      }
      try{
        dynamic_cast< TH1D* > (_histoThing[histogramnamey.str().c_str()])->Fill(y2 - predicty);
      } catch(std::bad_cast &bc){
        streamlog_out(ERROR3) << "Unable to fill histogram: " << histogramnamey.str().c_str() << endl;
      }     
    }
    ++i;
  }
}

void EUTelX0Processor::threePointResolution(Track *track){
//This function draws a line between the hits in plane i and i+2
//then it compares where the hit in plane i+1 is with the average of the other two planes
//this is then plotted as a residual for each plane.
  streamlog_out(DEBUG1) << "Function EUTelX0Processor::threePointResolution(Track *" << &track << ") called" << std::endl;

  std::vector< TrackerHit* > trackhits = track->getTrackerHits();
  std::vector< TVector3* > hits = getHitsFromTrack(track);
  int i = 1;
  for(std::vector< TVector3* >::iterator it = hits.begin(); it != hits.end() - 2; ++it){
    //This determines the guess of the position of the particle as it should hit the middle sensor
    streamlog_out(DEBUG1) << "In the for loop of the three point resolution function" << endl;
    double averagingfactor = ((*it + 1)->z() - (*it)->z())/((*it + 2)->z() - (*it)->z());
    if(averagingfactor == 0){
      streamlog_out(ERROR5) << "Averaging factor == 0)" << endl;
      continue;
    }
    double averagex = ((*it)->x() + (*it + 2)->x())/averagingfactor;
    double averagey = ((*it)->y() + (*it + 2)->x())/averagingfactor;
    double middlex = (*it + 1)->x();
    double middley = (*it + 1)->y();
    std::stringstream ResidualX, ResidualY;
    ResidualX << "ThreePointResidualXPlane" << i;
    ResidualX << "ThreePointResidualYPlane" << i;
    try{
      dynamic_cast< TH1D* > (_histoThing[ResidualX.str().c_str()])->Fill(averagex - middlex);
    } catch(std::bad_cast& bc){
      streamlog_out(ERROR5) << "Unable to fill histogram " << ResidualX.str() << " due to bad cast" << std::endl;
    }
    try{
      dynamic_cast< TH1D* > (_histoThing[ResidualY.str().c_str()])->Fill(averagey - middley);
    } catch(std::bad_cast& bc){
      streamlog_out(ERROR5) << "Unable to fill histogram " << ResidualY.str() << " due to bad cast" << std::endl;
    }
    streamlog_out(DEBUG1) << "Made it to the end of the for loop" << endl;
    i++;
  }
}

