#ifndef MEAS_MOSAIC_snapshot_h_INCLUDED
#define MEAS_MOSAIC_snapshot_h_INCLUDED

#include "lsst/meas/mosaic/mosaicfit.h"

namespace lsst { namespace meas { namespace mosaic {

void writeObsVec(std::string const & filename, ObsVec const & obsVec);

}}} // namespace lsst::meas::mosaic

#endif // !MEAS_MOSAIC_snapshot_h_INCLUDED
