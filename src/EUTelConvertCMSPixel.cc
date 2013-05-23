// Version: $Id: CMSPixelReader.cc 2557 2013-04-19 13:16:02Z spanns $

#ifdef USE_GEAR

// CMSPixel includes
#include "EUTelConvertCMSPixel.h"
#include "CMSPixelDecoder.h"

// EUTelescope includes
#include "EUTELESCOPE.h"
#include "EUTelSimpleSparsePixel.h"
#include "EUTelBaseSparsePixel.h"
#include "EUTelSparseDataImpl.h"
#include "EUTelExceptions.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"

// GEAR includes
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// Marlin includes
#include "marlin/Processor.h"
#include "marlin/DataSourceProcessor.h"
#include "marlin/ProcessorMgr.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
// Marlin AIDA include
#include "marlin/AIDAProcessor.h"
// AIDA includes
#include <AIDA/IBaseHistogram.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/ITree.h>
#endif

// LCIO includes
#include <IMPL/LCEventImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/LCTime.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <Exceptions.h>

// System includes
#include <vector>
#include <map>
#include <string.h>
#include <memory>

using namespace std;
using namespace marlin;
using namespace CMSPixel;
using namespace eutelescope;


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelConvertCMSPixel::_hitMapHistoName             	= "hitMap";
std::string EUTelConvertCMSPixel::_pulseHeightHistoName          	= "pulseHeight";
#endif


EUTelConvertCMSPixel::EUTelConvertCMSPixel ():DataSourceProcessor  ("EUTelConvertCMSPixel") {

  _description =
    "Reads PSI46 testboard data files and creates LCEvents (zero suppressed data as sparsePixel).\n"
    "Make sure to not specify any LCIOInputFiles in the steering in order to read CMSPixel files.";

  registerProcessorParameter ("FileName", "Input file",
                              _fileName, std::string ("mtb.bin"));

  registerProcessorParameter ("addressLevelsFile", "Address levels calibration file for the TBM and ROC address encoding levels.",
                              _levelsFile, std::string ("addressParameters.dat"));

  registerProcessorParameter ("runNumber", "RunNumber",
                              _srunNumber, std::string("000001"));

  registerProcessorParameter ("writeEmptyEvents", "Enable or disable the writing of events with no hit in all sensor planes.",
                              _writeEmptyEvents, static_cast < bool >(false));

  registerOutputCollection (LCIO::TRACKERDATA, "sparseDataCollectionName", "Name of the output sparsified data collection",
                            _sparseDataCollectionName, string("sparse"));

  registerProcessorParameter ("ROC_type", "Choose the ROC type. Can be: psi46v2, psi46xdb, psi46dig_trig, psi46dig, psi46digv2_b, psi46digv2.",
                              _ROC_type, static_cast < int >(0));

  registerProcessorParameter ("TB_type", "Choose the Testboard type. Can be: PSI_ATB, PSI_DTB, RAL",
                              _TB_type, static_cast < std::string >("RAL"));



  IntVec std_planes;
  for ( size_t iDetector = 0 ; iDetector < 16; ++iDetector ) {
    std_planes.push_back( iDetector );
  }
  
  registerOptionalParameter ("shufflePlanes", "int vector to hold the telescope plane IDs in the order in which they get the readout token.",
			     _shufflePlanes, std_planes);
                              
  //  registerOptionalParameter("eventSelection","Select the events to process: 0 - all, 1 - only with correct No. of ROC headers, 2 - only with corr. ROC headers and without bit errors in them.", _event_selection, static_cast< int > ( 0 ) );

registerOptionalParameter ("haveTBMheaders", "Switch TBM mode on and off. This gives the possibility to read data without real or emulated TBM headers as seen from soem testboard FPGAs. TRUE will look for TBM headers and trailers.",
			   _haveTBM, static_cast < bool >(true));

//gisterOptionalParameter ("useIPBus", "Switch between 4/12bit readout (PSI ATB/DTB) and IPBus 16bit readout (RAL board).", _useIPBus, static_cast < bool >(false));

 registerOptionalParameter("debugDecoder","Set decoder verbosity level: QUIET, SUMMARY, ERROR, WARNING, INFO, DEBUG, DEBUG1-4", _debugSwitch, static_cast< std::string > ( "SUMMARY" ) );
registerProcessorParameter("HistogramFilling","Switch on or off the histogram filling", _fillHistos, static_cast< bool > ( true ) );
                              
}


void EUTelConvertCMSPixel::init () {
 
  // Print the processor parameters:
  printParameters ();
    
  // Set processor back, initialize the runnumber:
  _runNumber = atoi(_srunNumber.c_str());
  _isFirstEvent = true;
  eventNumber = 0;
  iROC = 0;
    
  if(_haveTBM) flags += FLAG_HAVETBM;
  if(_useIPBus) flags += FLAG_16BITS_PER_WORD;

  // Set Decoder Logging level:
  Log::ReportingLevel() = Log::FromString(_debugSwitch);
    
}


void EUTelConvertCMSPixel::initializeGeometry() 
{
  streamlog_out( MESSAGE5 ) << "Initializing geometry" << endl;

  _noOfROC = 0;
  _noOfXPixel = 0;
  _noOfYPixel = 0;	

  _siPlanesParameters  = const_cast< gear::SiPlanesParameters*  > ( &(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast< gear::SiPlanesLayerLayout* > ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

  _layerIndexMap.clear();
  for ( int iLayer = 0; iLayer < _siPlanesLayerLayout->getNLayers(); ++iLayer ) {
    _layerIndexMap.insert( make_pair( _siPlanesLayerLayout->getID( iLayer ), iLayer ) );
  }

  _noOfROC = _siPlanesLayerLayout->getNLayers();
	
  // We only use identical telescope planes, so reading the parameters from the first should be fine:
  _noOfXPixel = _siPlanesLayerLayout->getSensitiveNpixelX( _layerIndexMap[0] );
  _noOfYPixel = _siPlanesLayerLayout->getSensitiveNpixelY( _layerIndexMap[0] );

  if ( _noOfROC == 0 || _noOfXPixel == 0 || _noOfYPixel == 0 ) {
    streamlog_out( WARNING ) << "Unable to initialize the geometry. Please check GEAR file." << endl;
    _isGeometryReady = false;
  } else {
    _isGeometryReady = true;
  }
  streamlog_out( MESSAGE5 ) << "Active SensorPlanes: " << _noOfROC << endl;
  streamlog_out( MESSAGE5 ) << "Pixels in X: " << _noOfXPixel << endl;
  streamlog_out( MESSAGE5 ) << "Pixels in Y: " << _noOfYPixel << endl;

}


void EUTelConvertCMSPixel::readDataSource (int Ntrig) 
{
   
  EUTelEventImpl *evt = NULL;
  CMSPixelFileDecoder *readout;
  // Initialize event vector:
  std::vector< event > event_data;
  int64_t timestamp;

  // Initialize geometry:
  initializeGeometry();
  if(!_isGeometryReady) throw InvalidGeometryException ("Wrong geometry file?");

  for(unsigned int i = 0; i < _noOfROC; i++) {
    streamlog_out( DEBUG5 ) << "ROC " << i << " -> TelescopePlane " << _shufflePlanes[i] << endl;
  }
       
  // Construct new decoder:
  if(strcmp(_TB_type.c_str(),"RAL") == 0) {
    streamlog_out( DEBUG5 ) << "Constructing RAL testboard decoder..." << endl;
    readout = new CMSPixelFileDecoderRAL(_fileName.c_str(),_noOfROC,flags,_ROC_type);
  }
  else if(strcmp(_TB_type.c_str(),"PSI_DTB") == 0) {
    streamlog_out( DEBUG5 ) << "Constructing PSI_DTB testboard decoder..." << endl;
    readout = new CMSPixelFileDecoderPSI_DTB(_fileName.c_str(),_noOfROC,flags,_ROC_type,_levelsFile.c_str());
  }
  else if(strcmp(_TB_type.c_str(),"PSI_ATB") == 0) {
    streamlog_out( DEBUG5 ) << "Constructing PSI_ATB testboard decoder..." << endl;
    readout = new CMSPixelFileDecoderPSI_ATB(_fileName.c_str(),_noOfROC,flags,_ROC_type,_levelsFile.c_str());
  }
  else {
    throw DataNotAvailableException("Could not determine correct testboard type. Check configuration.");
  }


  // Loop while we have input data
  while (true)
    {
      // Check if it's the first event:
      if(_isFirstEvent) {
	// We are in the first event, so type BORE. Write the run header.
	auto_ptr<IMPL::LCRunHeaderImpl> lcHeader  ( new IMPL::LCRunHeaderImpl );
	auto_ptr<EUTelRunHeaderImpl>    runHeader ( new EUTelRunHeaderImpl (lcHeader.get()) );
	runHeader->addProcessor( type() );
	runHeader->lcRunHeader()->setDescription(" Events read from CMSPixel input file: " + _fileName);
	runHeader->lcRunHeader()->setRunNumber (_runNumber);
	runHeader->setHeaderVersion (0.0001);
	runHeader->setDataType (EUTELESCOPE::CONVDATA);
	runHeader->setDateTime ();
	runHeader->addIntermediateFile (_fileName);
	runHeader->addProcessor (_processorName);
	runHeader->setNoOfDetector(_noOfROC);
	runHeader->setMinX(IntVec(_noOfROC, 0));
	runHeader->setMaxX(IntVec(_noOfROC, _noOfXPixel - 1));
	runHeader->setMinY(IntVec(_noOfROC, 0));
	runHeader->setMaxY(IntVec(_noOfROC, _noOfYPixel - 1));

	runHeader->lcRunHeader()->setDetectorName("CMSPixelTelescope");

	// Process the run header:
	ProcessorMgr::instance ()->processRunHeader ( static_cast<lcio::LCRunHeader*> ( lcHeader.release()) );
            
	// Book histogramms:
	if ( _fillHistos ) bookHistos();

	_isFirstEvent = false;
      }

      // Trigger counter:
      if(eventNumber >= Ntrig) {
	streamlog_out ( MESSAGE5 ) << " ---------------------------------------------------------" << endl;
	streamlog_out ( MESSAGE5 ) << "  End of processing: reached MaxRecordNumber (" << Ntrig << ")" << endl;
	streamlog_out ( MESSAGE5 ) << "  If you want to process more events check your steerfile." << endl;                                    
	break;
      }
      eventNumber++;

    
      // Clear vector and read next event from file, containing all ROCs / pixels for one trigger:
      event_data.clear();
        
      // Read next event from data source:
      status = readout->get_event(&event_data,timestamp);
        
      if(status <= DEC_ERROR_NO_MORE_DATA) {
	streamlog_out (ERROR) << "Decoder returned error " << status << std::endl;
	// We didn't write single event - it was just impossible to open/read the data file:
	if(eventNumber == 1) {
	  streamlog_out ( WARNING ) << "The data file contained no valid event." << endl;
	  throw DataNotAvailableException("Failed to read from data file.");
	}
	// Else: we just reached EOF.
	else break;
      }
      else if(status <= DEC_ERROR_NO_TBM_HEADER ) {
	streamlog_out (WARNING5) << "There was an exception while processing event " << eventNumber << ". Will continue with next event." << std::endl;
	continue;
      }
      else if(status <= DEC_ERROR_INVALID_ROC_HEADER) {
	streamlog_out (DEBUG5) << "Issue with ROC header or pixel address in event " << eventNumber << ". Event will be used anyway." << std::endl;
      }
      else if(_writeEmptyEvents && status == DEC_ERROR_EMPTY_EVENT) {
	streamlog_out (DEBUG5) << "Event " << eventNumber << " is empty. Continuing with next." << std::endl;
	continue;
      }


      LCCollectionVec * sparseDataCollection = new LCCollectionVec(LCIO::TRACKERDATA);

      TrackerDataImpl * sparse;
      EUTelSparseDataImpl<EUTelSimpleSparsePixel>  sparseData( sparse ) ;

      int16_t iROC = -1;
      for (unsigned int ipx = 0; ipx < event_data.size(); ipx++) {
	
	if(iROC != event_data[ipx].roc) {
	  iROC = event_data[ipx].roc;
	  sparse = new TrackerDataImpl();
	  CellIDEncoder<TrackerDataImpl> sparseDataEncoder(EUTELESCOPE::ZSDATADEFAULTENCODING, sparseDataCollection);
	  sparseDataEncoder["sensorID"]        = iROC;
	  sparseDataEncoder["sparsePixelType"] = static_cast<int> ( 1 );
	  sparseDataEncoder.setCellID(sparse);
	  
	}

	auto_ptr<EUTelSimpleSparsePixel> sparsePixel( new EUTelSimpleSparsePixel );
                    
	sparsePixel->setXCoord( static_cast<int> (event_data[ipx].col ));
	sparsePixel->setYCoord( static_cast<int> (event_data[ipx].row ));
	sparsePixel->setSignal( static_cast<short> (event_data[ipx].raw));
	streamlog_out ( DEBUG0 ) << (*sparsePixel.get()) << endl;
        
	//fill histogramms if necessary:
	if ( _fillHistos ) fillHistos ( event_data[ipx].col, event_data[ipx].row, event_data[ipx].raw, event_data[ipx].roc );
	sparseData.addSparsePixel( sparsePixel.get() );
	sparseDataCollection->push_back( sparse );
      }

      delete sparse;


      // Start constructing current event:
      evt = new EUTelEventImpl();
      evt->setDetectorName("CMSPixelTelescope");
      evt->setEventType(kDE);
      evt->setRunNumber (_runNumber);
      evt->setEventNumber (eventNumber);
      LCTime * now = new LCTime();
      evt->setTimeStamp(timestamp);
      delete now;
      // ...and write it out:
      evt->addCollection (sparseDataCollection, _sparseDataCollectionName);
      ProcessorMgr::instance ()->processEvent (static_cast<LCEventImpl*> (evt));
        
      delete evt;

    }; // end of while (true)


  // Write last event with type EORE
  eventNumber++;    
  evt = new EUTelEventImpl();
  evt->setDetectorName("CMSPixelTelescope");
  evt->setEventType(kEORE);
  evt->setRunNumber (_runNumber);
  evt->setEventNumber (eventNumber);
  LCTime * now = new LCTime();
  evt->setTimeStamp(now->timeStamp());
  delete now;
  ProcessorMgr::instance ()->processEvent (static_cast<LCEventImpl*> (evt));
  streamlog_out ( MESSAGE5 ) << " ---------------------------------------------------------" << endl;    
  streamlog_out ( MESSAGE5 ) << "  Write EORE as event " << evt->getEventNumber() << endl;
        
  // Print the readout statistics, invoked by the destructor:
  streamlog_out ( MESSAGE5 ) << " ---------------------------------------------------------" << endl;    
  readout->statistics.print();
  delete[] readout;
  streamlog_out ( MESSAGE5 ) << " ---------------------------------------------------------" << endl;    

  // Delete the EORE event:    
  delete evt;
}


void EUTelConvertCMSPixel::end () {
  message<MESSAGE5> ("Successfully finished") ;
}


#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
void EUTelConvertCMSPixel::fillHistos (int xCoord, int yCoord, int value, int sensorID) {

  string tempHistoName;
			
  tempHistoName = _hitMapHistoName + "_d" + to_string( sensorID );
  (dynamic_cast<AIDA::IHistogram2D*> (_aidaHistoMap[tempHistoName]))->fill(static_cast<double >(xCoord), static_cast<double >(yCoord), 1.);
			
  tempHistoName = _pulseHeightHistoName + "_d" + to_string( sensorID );
  (dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]))->fill(value);
}


void EUTelConvertCMSPixel::bookHistos() {
	
  streamlog_out ( MESSAGE5 )  << "Booking histograms " << endl;

  string tempHistoName;
  string basePath;
  for (unsigned int iDetector = 0; iDetector < _noOfROC; iDetector++) {

    basePath = "detector_" + to_string( iDetector );
    AIDAProcessor::tree(this)->mkdir(basePath.c_str());
    basePath.append("/");

    tempHistoName = _hitMapHistoName + "_d" + to_string( iDetector );
    AIDA::IHistogram2D * hitMapHisto =
      AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(), _noOfXPixel, 0, _noOfXPixel, _noOfYPixel, 0, _noOfYPixel);
    _aidaHistoMap.insert(make_pair(tempHistoName, hitMapHisto));
    hitMapHisto->setTitle("Hit map");

    string pulseHeightTitle = "pulse height ROC" + to_string( iDetector );
    tempHistoName = _pulseHeightHistoName + "_d" + to_string( iDetector );
    AIDA::IHistogram1D * pulseHeightHisto = AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), 525,-1050,1050);
    _aidaHistoMap.insert(make_pair(tempHistoName, pulseHeightHisto));
    pulseHeightHisto->setTitle(pulseHeightTitle.c_str());
				
  }
  streamlog_out ( MESSAGE5 )  << "end of Booking histograms " << endl;
}
#endif // USE_AIDA || MARLIN_USE_AIDA

#endif // USE_GEAR
