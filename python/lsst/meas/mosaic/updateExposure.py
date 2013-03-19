from .mosaicLib import getFCorImg, FluxFitParams
import lsst.afw.image

__all__ = ("applyMosaicResults",)

def applyMosaicResults(dataRef, calexp=None):
    """Update an Exposure with the Wcs, Calib, and flux scaling from meas_mosaic.

    If None, the calexp will be loaded from the dataRef.
    """
    if calexp is None:
        calexp = dataRef.get("calexp", immediate=True)
    wcs_md = dataRef.get("wcs_md", immediate=True)
    wcs = lsst.afw.image.makeWcs(wcs_md)
    calib = lsst.afw.image.Calib(wcs_md)
    calexp.setWcs(wcs)
    calexp.setCalib(calib)
    ffp_md = dataRef.get("fcr_md", immediate=True)
    ffp = FluxFitParams(ffp_md)
    mi = calexp.getMaskedImage()
    fcor = getFCorImg(ffp, mi.getWidth(), mi.getHeight())
    mi *= fcor
    return calexp
