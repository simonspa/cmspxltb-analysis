/* 
 * File:   EUTelTrackFinder.h
 * Contact: denys.lontkovskyi@desy.de
 *
 * Created on January 22, 2013, 2:13 PM
 */

#ifndef EUTELTRACKFINDER_H
#define	EUTELTRACKFINDER_H

// eutelescope includes ".h"
//#include "EUTelUtility.h"

// lcio includes <.h>
#include <IMPL/TrackerHitImpl.h>

// system includes <>
#include <map>
#include <string>
#include <vector>


namespace eutelescope {

    
    class EUTelTrackFinder {
        
    public:
        EUTelTrackFinder();
 //       EUTelTrackFinder( const EUTelTrackFinder& orig );
        EUTelTrackFinder( std::string name );
        
        virtual ~EUTelTrackFinder();

    public:
        enum State {
            kIsNotReady, kIsReady
        };

        enum SearchResult {
            kFailed, kSuccess
        };
        
    public:
        virtual void Reset();
        virtual EUTelTrackFinder::SearchResult SearchTracks();
        virtual std::map<int, EVENT::TrackerHitVec> GetTrackCandidates() const;
        
        void SetAllHits( const std::vector< EVENT::TrackerHitVec >& _allHits);
        std::vector< EVENT::TrackerHitVec > GetAllHits() const;

        inline void SetName( std::string name) { this->_name = name; }
        inline std::string GetName() const { return _name; }
        
    protected:
        virtual EUTelTrackFinder::SearchResult DoTrackSearch();
        
    protected:        
        State _isReady;
        SearchResult _searchResult;
                
        std::string _name;
        
        std::vector< EVENT::TrackerHitVec > _allHits;
        std::map< int, EVENT::TrackerHitVec > _trackCandidates;

    };

}
#endif	/* EUTELTRACKFINDER_H */

