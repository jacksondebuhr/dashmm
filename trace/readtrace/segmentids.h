#ifndef __TRACEUTILS_SEGMENT_IDS_H__
#define __TRACEUTILS_SEGMENT_IDS_H__


#include <string>


namespace traceutils {


namespace segment {


// This indicates that an event is not part of a segment
constexpr int kNoSegment = 0;

// The segment identifiers for network events
constexpr int kNetworkProgress = 1;

// The segment identifiers for DASHMM specific events
constexpr int kDASHMMStoT = 2;
constexpr int kDASHMMMtoT = 3;
constexpr int kDASHMMLtoT = 4;
constexpr int kDASHMMStoM = 5;
constexpr int kDASHMMStoL = 6;
constexpr int kDASHMMELCO = 7;
constexpr int kDASHMMMtoM = 8;
constexpr int kDASHMMMtoL = 9;
constexpr int kDASHMMLtoL = 10;
constexpr int kDASHMMMtoI = 11;
constexpr int kDASHMMItoI = 12;
constexpr int kDASHMMItoL = 13;


std::string name(int seg);


} // traceutils::segment


} // traceutils


#endif // __TRACEUTILS_SEGMENT_IDS_H__