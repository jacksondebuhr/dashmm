#include "segmentids.h"


namespace traceutils::segment {


std::string name(int seg) {
  switch (seg) {
    case kNoSegment:
      return std::string("NoSeg");
      break;
    case kNetworkProgress:
      return std::string("NetProg");
      break;
    case kDASHMMStoT:
      return std::string("StoT");
      break;
    case kDASHMMMtoT:
      return std::string("MtoT");
      break;
    case kDASHMMLtoT:
      return std::string("LtoT");
      break;
    case kDASHMMStoM:
      return std::string("StoM");
      break;
    case kDASHMMStoL:
      return std::string("StoL");
      break;
    case kDASHMMELCO:
      return std::string("ELCO");
      break;
    case kDASHMMMtoM:
      return std::string("MtoM");
      break;
    case kDASHMMMtoL:
      return std::string("MtoL");
      break;
    case kDASHMMLtoL:
      return std::string("LtoL");
      break;
    case kDASHMMMtoI:
      return std::string("MtoI");
      break;
    case kDASHMMItoI:
      return std::string("ItoI");
      break;
    case kDASHMMItoL:
      return std::string("ItoL");
      break;
    default:
      return std::string("Unknown");
      break;
  }
}


} // traceutils::segment
