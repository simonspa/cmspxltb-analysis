// Version: $Id$
#ifndef EUTelProcessorFilteringHitFilter_h
#define EUTelProcessorFilteringHitFilter_h 1

// LCIO
#include "lcio.h"

#include "marlin/Processor.h"

#include "IMPL/TrackerHitImpl.h"

// C++
#include <string>


using namespace lcio;
using namespace marlin;

namespace eutelescope {

    class EUTelProcessorFilteringHitFilter : public Processor {
    public:

        virtual Processor* newProcessor() {
            return new EUTelProcessorFilteringHitFilter;
        }

        EUTelProcessorFilteringHitFilter();

        /** Called at the begin of the job before anything is read.
         * Use to initialize the processor, e.g. book histograms.
         */
        virtual void init();

        /** Called for every run.
         */
        virtual void processRunHeader(LCRunHeader* run);

        /** Called for every event - the working horse.
         */
        virtual void processEvent(LCEvent * evt);

        virtual void check(LCEvent * evt);

        /** Called after data processing for clean up.
         */
        virtual void end();


    protected:

        //! Input TrackerHit collection name
        std::string _hitInputCollectionName;
        
        //! Output TrackerHit collection name
        std::string _hitOutputCollectionName;

        IntVec _wantPlaneID;
        
        int _nProcessedRuns;
        int _nProcessedEvents;

    };

    //! A global instance of the processor
    EUTelProcessorFilteringHitFilter gEUTelProcessorFilteringHitFilter;
    
} // eutelescope

#endif // EUTelProcessorFilteringHitFilter_h
