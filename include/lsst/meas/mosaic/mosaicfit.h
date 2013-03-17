// -*- lsst-c++ -*-
#if !defined(HSC_MEAS_MOSAIC_H)
#define HSC_MEAS_MOSAIC_H

#include <vector>
#include "lsst/afw/image.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/cameraGeom.h"
#include "lsst/afw/table.h"
#include "lsst/utils/ieee.h"
#include "boost/enable_shared_from_this.hpp"

namespace lsst {
    namespace meas {
	namespace mosaic {

            typedef boost::int64_t IdType;
            typedef boost::int32_t ChipType;
            typedef boost::int64_t ExpType;

            class Source {
            public:
                enum { UNSET = -1 };
                explicit Source(lsst::afw::table::SourceRecord const& record) :
                    _id(record.getId()), _chip(UNSET), _exp(UNSET), _sky(record.getRa(), record.getDec()),
                    _pixels(record.getX(), record.getY()), _flux(record.getApFlux()), _err(record.getApFluxErr()),
		    _xerr(sqrt(record.getCentroidErr()(0,0))), _yerr(sqrt(record.getCentroidErr()(1,1))),
                    _astromBad(record.getCentroidFlag() | record.getApFluxFlag()) {}
                Source(lsst::afw::table::SimpleRecord const& record, lsst::afw::image::Wcs const& wcs) :
                    _id(record.getId()), _chip(UNSET), _exp(UNSET), _sky(record.getRa(), record.getDec()),
                    _pixels(wcs.skyToPixel(_sky)),
                    _flux(record.get(record.getSchema().find<double>("flux").key)),
                    _err(record.get(record.getSchema().find<double>("flux.err").key)),
                    _astromBad(!lsst::utils::isfinite(_flux)) {}
                Source(lsst::afw::coord::Coord coord, double flux=std::numeric_limits<double>::quiet_NaN()) :
                    _id(-1), _chip(UNSET), _exp(UNSET), _sky(coord),
                    _pixels(lsst::afw::geom::Point2D(std::numeric_limits<double>::quiet_NaN(),
                                                     std::numeric_limits<double>::quiet_NaN())),
                    _flux(flux), _astromBad(false) {}

                IdType getId() const { return _id; }
                ChipType getChip() const { return _chip; }
                ExpType getExp() const { return _exp; }
                lsst::afw::coord::Coord getSky() const { return _sky; }
                lsst::afw::geom::Angle getRa() const { return getSky().getLongitude(); }
                lsst::afw::geom::Angle getDec() const { return getSky().getLatitude(); }
                lsst::afw::geom::Point2D getPixels() const { return _pixels; }
                double getX() const { return getPixels().getX(); }
                double getY() const { return getPixels().getY(); }
		double getXErr() const { return _xerr; }
		double getYErr() const { return _yerr; }
                double getFlux() const { return _flux; }
                double getFluxErr() const { return _err; }
                bool getAstromBad() const { return _astromBad; }

                void setChip(ChipType chip) { _chip = chip; }
                void setExp(ExpType exp) { _exp = exp; }
                
            private:
                IdType _id;                       // Identifier
                ChipType _chip;                   // Chip identifier
                ExpType _exp;                     // Exposure identifier
                lsst::afw::coord::Coord _sky;     // Sky coordinates
                lsst::afw::geom::Point2D _pixels; // Pixel coordinates
                double _flux;                     // Flux
                double _err;			  // Flux Err
		double _xerr;			  // x coordinate error
		double _yerr;			  // y coordinate error
                bool _astromBad;                  // Astrometry bad?
            };

            typedef std::pair<PTR(Source), PTR(Source)> SourceMatch;

            typedef std::vector<PTR(Source)> SourceSet;
	    typedef std::vector<std::vector<PTR(Source)> > SourceGroup;
            typedef std::vector<SourceMatch> SourceMatchSet;
	    typedef std::vector<std::vector<SourceMatch> > SourceMatchGroup;

	    typedef std::map<ExpType, lsst::afw::image::Wcs::Ptr> WcsDic;

	    class Poly {
	    public:
	        typedef boost::shared_ptr<Poly> Ptr;

		int order;
		int ncoeff;
		int *xorder;
		int *yorder;

		Poly(int order);
		~Poly(void);
		Poly(const Poly &p);
		int getIndex(int i, int j);
		int getXorder(int i) { return xorder[i]; }
		int getYorder(int i) { return yorder[i]; }
	    };

	    class Coeff {
	    public:
		typedef boost::shared_ptr<Coeff> Ptr;

		Poly::Ptr p;
		ExpType iexp;
		double *a;
		double *b;
		double *ap;
		double *bp;
		double A;
		double D;
		double x0;
		double y0;
		
		Coeff(int order);
		Coeff(Poly::Ptr const & p);
		~Coeff(void);
		Coeff(const Coeff &c);
		void show(void);
		void uvToXiEta(double u, double v, double *xi, double *eta);
		void xietaToUV(double xi, double eta, double *u, double *v);
		double get_a(int i) { return a[i]; }
		double get_b(int i) { return b[i]; }
		double get_ap(int i) { return ap[i]; }
		double get_bp(int i) { return bp[i]; }
		void set_a(int i, double v) { a[i] = v; }
		void set_b(int i, double v) { b[i] = v; }
		void set_ap(int i, double v) { ap[i] = v; }
		void set_bp(int i, double v) { bp[i] = v; }
		double xi(double u, double v);
		double eta(double u, double v);
		double dxidu(double u, double v);
		double dxidv(double u, double v);
		double detadu(double u, double v);
		double detadv(double u, double v);
		double detJ(double u, double v);
		int getNcoeff() { return p->ncoeff; }
		double pixelScale(void);

                // below lines for coeffFromTanWcs()
		void set_D(double v) { D = v; }
		void set_A(double v) { A = v; }
		void set_x0(double v) { x0 = v; }
		void set_y0(double v) { y0 = v; }
		void set_iexp(ExpType v) { iexp = v; }
		double get_D() { return D; }
		double get_A() { return A; }
		double get_x0() { return x0; }
		double get_y0() { return y0; }
		ExpType get_iexp() { return iexp; }
	    };

	    class Obs {
	    public:
		typedef boost::shared_ptr<Obs> Ptr;

		double ra, dec;
		double xi, eta;
		double xi_a, xi_d, eta_a, eta_d;
		double xi_A, xi_D, eta_A, eta_D;

		double x, y;
		double u, v;
		double u0, v0;

		double U, V;

		double xi_fit, eta_fit;
		double u_fit, v_fit;

		int id;
		int istar;
		int jstar;	/* index for fit */
		ExpType iexp;	/* exposure id or visit number */
		ChipType ichip;	/* ccdId */
		int jexp;	/* sequential index for exposure used for matrix*/
		int jchip;	/* sequential index for ccd used for matrix*/
		bool good;

		double xerr;
		double yerr;

		double mag;
		double mag0;
		double err;
		double mag_cat;
		double err_cat;

		Obs(int id, double ra, double dec, double x, double y, ChipType ichip, ExpType iexp);
		Obs(int id, double ra, double dec, ChipType ichip, ExpType iexp);
		void setUV(lsst::afw::cameraGeom::Ccd::Ptr const &ccd, double x0=0, double y0=0);
		void setXiEta(double ra_c, double dec_c);
		void setFitVal(Coeff::Ptr& c, Poly::Ptr p);
		void setFitVal2(Coeff::Ptr& c, Poly::Ptr p);
	    };

	    class KDTree
#if !defined(SWIG)
	      : public boost::enable_shared_from_this<KDTree>
#endif
	    {
	    public:
		typedef boost::shared_ptr<KDTree> Ptr;
		typedef boost::shared_ptr<const KDTree> ConstPtr;

		int depth;
		int axis;
                lsst::afw::geom::Angle location[2];
		lsst::afw::coord::Coord c;
		KDTree::Ptr left;
		KDTree::Ptr right;
                SourceSet set;

                KDTree(SourceSet& s, int depth) { _initializeSources(s, depth); }
		KDTree(PTR(Source) s, int depth);

                KDTree(SourceMatchSet m, int depth) { _initializeMatches(m, depth); }
                KDTree(SourceMatch const& m, int depth);

		~KDTree();
                ConstPtr search(lsst::afw::coord::Coord const& sky) const;
		ConstPtr findSource(Source const& s) const;
                void add(SourceMatch const& m);
		void add(PTR(Source) s,
                         lsst::afw::geom::Angle d_lim=lsst::afw::geom::Angle(0, lsst::afw::geom::degrees));
		int count(void);
		SourceGroup mergeMat() const;
		SourceGroup mergeSource();
		void printMat() const;
		void printSource() const;
		bool isLeaf(void) const { return left == NULL && right == NULL; }
		KDTree::Ptr findNearest(Source const& s);

		double distance(Source const& s) const {
                    return c.angularSeparation(lsst::afw::coord::Coord(s.getRa(), s.getDec())).asDegrees();
                }

            private:
                void _initializeSources(SourceSet& s, int depth);
                void _initializeMatches(SourceMatchSet& m, int depth);
            };

	    class FluxFitParams {
	    public:
		typedef boost::shared_ptr<FluxFitParams> Ptr;

		int order;
		bool chebyshev;
		int ncoeff;
		int *xorder;
		int *yorder;
		bool absolute;

		double *coeff;
		double u_max;
		double v_max;
		double x0;
		double y0;

		FluxFitParams(int order, bool absolute=false, bool chebyshev=false);
		FluxFitParams(lsst::daf::base::PropertySet::Ptr& metadata);
		~FluxFitParams();
		FluxFitParams(const FluxFitParams &p);
		double eval(double u, double v);
		int getXorder(int i) { return xorder[i]; }
		int getYorder(int i) { return yorder[i]; }
		int getCoeff(int i) { return coeff[i]; }
		int getIndex(int i, int j);
	    };

	    typedef std::map<ChipType, lsst::afw::cameraGeom::Ccd::Ptr> CcdSet;
	    typedef std::map<ExpType, Coeff::Ptr> CoeffSet;
	    typedef std::vector<Obs::Ptr> ObsVec;

	    KDTree::Ptr kdtreeMat(SourceMatchGroup &matchList);
	    KDTree::Ptr kdtreeSource(SourceGroup const &sourceSet,
				     KDTree::Ptr rootMat,
				     CcdSet &ccdSet,
				     lsst::afw::geom::Angle d_lim, unsigned int nbrightest);

	    ObsVec obsVecFromSourceGroup(SourceGroup const &all,
					 WcsDic &wcsDic,
					 CcdSet &ccdSet);

	    CoeffSet solveMosaic_CCD_shot(int order,
					  int nmatch,
					  ObsVec &matchVec,
					  WcsDic &wcsDic,
					  CcdSet &ccdSet,
					  FluxFitParams::Ptr &ffp,
					  std::map<ExpType, float> &fexp,
					  std::map<ChipType, float> &fchip,
					  bool solveCcd = true,
					  bool allowRotation = true,
					  bool verbose = false,
					  double catRMS = 0.0,
                                          bool writeSnapshots = false,
                                          std::string const & snapshotDir = ".");

	    CoeffSet solveMosaic_CCD(int order,
				     int nmatch,
				     int nsource,
				     ObsVec &matchVec,
				     ObsVec &sourceVec,
				     WcsDic &wcsDic,
				     CcdSet &ccdSet,
				     FluxFitParams::Ptr &ffp,
				     std::map<ExpType, float> &fexp,
				     std::map<ChipType, float> &fchip,
				     bool solveCcd = true,
				     bool allowRotation = true,
				     bool verbose = false,
				     double catRMS = 0.0,
                                     bool writeSnapshots = false,
                                     std::string const & snapshotDir = ".");

	    Coeff::Ptr convertCoeff(Coeff::Ptr& coeff,
				    lsst::afw::cameraGeom::Ccd::Ptr& ccd);

	    lsst::afw::image::TanWcs::Ptr wcsFromCoeff(Coeff::Ptr& coeff);

            Coeff::Ptr coeffFromTanWcs(lsst::afw::image::Wcs::Ptr& wcs);

	    FluxFitParams::Ptr
	      convertFluxFitParams(Coeff::Ptr& coeff,
				   lsst::afw::cameraGeom::Ccd::Ptr& ccd,
				   FluxFitParams::Ptr& ffp);

	    lsst::daf::base::PropertySet::Ptr
	      metadataFromFluxFitParams(FluxFitParams::Ptr& ffp);

	    lsst::afw::image::Image<float>::Ptr
	      getJImg(Coeff::Ptr& coeff,
		      lsst::afw::cameraGeom::Ccd::Ptr& ccd);

	    lsst::afw::image::Image<float>::Ptr
	      getJImg(lsst::afw::image::Wcs::Ptr& wcs,
		      int width, int height);

	    lsst::afw::image::Image<float>::Ptr
	      getJImg(lsst::afw::image::Wcs::Ptr& wcs,
		      lsst::afw::cameraGeom::Ccd::Ptr& ccd);

	    lsst::afw::image::Image<float>::Ptr
	      getFCorImg(FluxFitParams::Ptr& p,
			 lsst::afw::cameraGeom::Ccd::Ptr& ccd,
			 Coeff::Ptr& coeff);

	    lsst::afw::image::Image<float>::Ptr
	      getFCorImg(FluxFitParams::Ptr& p, int width, int height);

	    lsst::afw::image::Image<float>::Ptr
	      getFCorImg(FluxFitParams::Ptr& p,
			 lsst::afw::cameraGeom::Ccd::Ptr& ccd);

#include "chebyshev.h"
    }
  }
}

#endif
