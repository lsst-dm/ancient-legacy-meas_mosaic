// -*- lsst-c++ -*-
%define fitLib_DOCSTRING
"
Python interface to lsst::meas::mosaic
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.meas.mosaic", docstring=fitLib_DOCSTRING) mosaicLib

%{
#include "lsst/afw/image.h"
#include "lsst/pex/logging.h"
#include "lsst/afw/detection.h"
#include "lsst/afw/math.h"


#define PY_ARRAY_UNIQUE_SYMBOL LSST_MEAS_MOSAIC_NUMPY_ARRAY_API
#include "numpy/arrayobject.h"
#include "ndarray/swig.h"
#include "ndarray/swig/eigen.h"
%}

%include "lsst/p_lsstSwig.i"
%include "ndarray.i"

%declareNumPyConverters(Eigen::MatrixXd);
%declareNumPyConverters(Eigen::VectorXd);

%import "lsst/afw/image/imageLib.i"
%import "lsst/afw/geom/geomLib.i"
%import "lsst/afw/detection/detectionLib.i"
%import "lsst/afw/math/mathLib.i"

%include "mosaic.i"
