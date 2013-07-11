/* 
 * File:   EUTelMagneticFieldFinder.h
 *
 * Created on July 2, 2013, 12:53 PM
 */

#ifndef EUTELMAGNETICFIELDFINDER_H
#define	EUTELMAGNETICFIELDFINDER_H

// system includes <>
#include <string>
#include <vector>

// EUTELESCOPE
#include "EUTelUtility.h"
#include "EUTelTrackFitter.h"

//LCIO
#include "lcio.h"
#include "IMPL/TrackStateImpl.h"
#include "IMPL/TrackerHitImpl.h"

#include "TLorentzVector.h"

class IMPL::TrackImpl;
class TrackerHit;

class MeasurementLayer;

namespace eutelescope {

    class EUTelKalmanFilter : public EUTelTrackFitter {
    private:
        DISALLOW_COPY_AND_ASSIGN(EUTelKalmanFilter) // prevent users from making (default) copies of processors

    public:
        EUTelKalmanFilter();

        explicit EUTelKalmanFilter(std::string name);

        virtual ~EUTelKalmanFilter();

        /** Fit supplied hits */
        virtual void FitTracks();

        /** Initialise Fitter */
        bool initialise();

        // Getters and Setters
    public:

        inline EVENT::TrackerHitVec getHits() const {
            return _allHits;
        };

        inline void setHits(EVENT::TrackerHitVec& hits) {
            this->_allHits = hits;
        }

        inline int getAllowedMissingHits() const {
            return _allowedMissingHits;
        };

        inline void setAllowedMissingHits(unsigned int allowedMissingHits) {
            this->_allowedMissingHits = allowedMissingHits;
        }

        inline int getMaxTrackCandidates() const {
            return _maxTrackCandidates;
        };

        inline void setMaxTrackCandidates(unsigned int maxTrackCandidates) {
            this->_maxTrackCandidates = maxTrackCandidates;
        }

        inline void setBeamMomentum(const TLorentzVector& beam) {
            this->_beamDir = beam;
        }

        inline TLorentzVector getBeamMomentum() const {
            return _beamDir;
        }

    private:
        /** Generate seed track candidates */
        void initialiseSeeds();

        /** Construct LCIO track object from internal track data */
        void prepareLCIOTrack();

        /** Sort hits according to particles propagation direction */
        void sortHitsByMeasurementLayers( const EVENT::TrackerHitVec& );
        

        // Kalman filter states and tracks
    private:
        /** Final set of tracks */
        std::vector< IMPL::TrackImpl* > _tracks;

        /** Kalman track states */
        std::vector< IMPL::TrackStateImpl* > _trackStates;

    private:
        /** Vector of hits to be processed */
        EVENT::TrackerHitVec _allHits;
        
        std::vector< MeasurementLayer* > _allMeasurements;

        // User supplied configuration of the fitter
    private:

        /** Validity of user input flag */
        bool _isReady;

        /** Maximum number of missing on a track candidate */
        int _allowedMissingHits;

        /** Maximum number of track candidates to be stored */
        int _maxTrackCandidates;

        /** Beam momentum vector */
        TLorentzVector _beamDir;
    };

    
    class MeasurementLayer {
    private:
        DISALLOW_COPY_AND_ASSIGN(MeasurementLayer)
        
    public:
        MeasurementLayer();
        
        explicit MeasurementLayer( int );
        
        virtual ~MeasurementLayer();

    public:
        /** Add hit */
        void addHit( EVENT::TrackerHit* );
        
    private:
        /** Measurement layer id */
        int _id;
        /** Set of hit belonging to the measurement layer */
        EVENT::TrackerHitVec _allHits;
    };
}

#endif	/* EUTELMAGNETICFIELDFINDER_H */

