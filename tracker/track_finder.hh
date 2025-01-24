#ifndef track_finder_HH
#define track_finder_HH

#include <iostream>
#include <memory>

#include "types.hh"

namespace Tracker
{
    class TrackFinder
    {
    public:
        TrackFinder();

    protected:
        HitList hits;
    };
} // namespace Tracker

#endif