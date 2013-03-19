%{
#include "lsst/meas/mosaic/mosaicfit.h"
%}

%include "std_vector.i"
%include "std_string.i"
%include "std_map.i"
%include "std_pair.i"

#if 1
%shared_ptr(lsst::afw::cameraGeom::Detector);
%shared_ptr(lsst::afw::cameraGeom::Ccd);
#else
%import "lsst/afw/cameraGeom/cameraGeom.i"
#endif

%shared_ptr(lsst::meas::mosaic::Source);
%shared_ptr(lsst::meas::mosaic::Coeff);
%shared_ptr(lsst::meas::mosaic::KDTree);
%shared_ptr(lsst::meas::mosaic::Obs);
%shared_ptr(lsst::meas::mosaic::FluxFitParams);

%include "lsst/meas/mosaic/mosaicfit.h"

%template(map_int_float) std::map<boost::int32_t, float>;
%template(map_int64_float) std::map<boost::int64_t, float>;

// These types need to match ChipType and ExpType in mosaicfit.h
// And have to be changed by hand if those types are changed
%pythoncode %{
    map_chiptype_float = map_int_float; 
%}
%pythoncode %{
    map_exptype_float = map_int64_float; 
%}

%template(SourceSet) std::vector<PTR(lsst::meas::mosaic::Source)>;
%template(SourceGroup) std::vector<std::vector<PTR(lsst::meas::mosaic::Source)> >;
%template(SourceMatch) std::pair<PTR(lsst::meas::mosaic::Source), PTR(lsst::meas::mosaic::Source)>;
%template(SourceMatchSet) std::vector<lsst::meas::mosaic::SourceMatch>;
%template(SourceMatchGroup) std::vector<std::vector<lsst::meas::mosaic::SourceMatch> >;

%template(WcsDic) std::map<lsst::meas::mosaic::ExpType, lsst::afw::image::Wcs::Ptr>;
%template(CcdSet) std::map<lsst::meas::mosaic::ChipType, lsst::afw::cameraGeom::Ccd::Ptr>;
%template(CoeffSet) std::map<lsst::meas::mosaic::ExpType, lsst::meas::mosaic::Coeff::Ptr>;
%template(ObsVec) std::vector<lsst::meas::mosaic::Obs::Ptr>;
