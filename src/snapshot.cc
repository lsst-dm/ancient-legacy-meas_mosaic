#include "lsst/meas/mosaic/snapshot.h"
#include "lsst/afw/table.h"

namespace lsst { namespace meas { namespace mosaic {

struct ObsVecKeys : private boost::noncopyable {
    afw::table::Schema schema;
    afw::table::Key<int> id;
    afw::table::Key<double> ra, dec;
    afw::table::Key<double> xi, eta;
    afw::table::Key<double> x, y;
    afw::table::Key<double> u, v;
    afw::table::Key<double> u0, v0;
    afw::table::Key<int> iexp, ichip;
    afw::table::Key<int> jexp, jchip;

    static ObsVecKeys const & get() {
        static ObsVecKeys const instance;
        return instance;
    }

private:
    ObsVecKeys() :
        schema(),
        id(schema.addField<int>("id", "")),
        ra(schema.addField<double>("ra", "")), dec(schema.addField<double>("dec", "")),
        xi(schema.addField<double>("xi", "")), eta(schema.addField<double>("eta", "")),
        x(schema.addField<double>("x", "")), y(schema.addField<double>("y", "")),
        u(schema.addField<double>("u", "")), v(schema.addField<double>("v", "")),
        u0(schema.addField<double>("u0", "")), v0(schema.addField<double>("v0", "")),
        iexp(schema.addField<int>("iexp", "")), ichip(schema.addField<int>("ichip", "")),
        jexp(schema.addField<int>("jexp", "")), jchip(schema.addField<int>("jchip", ""))
        {}
};

void writeObsVec(std::string const & filename, ObsVec const & obsVec) {
    ObsVecKeys const & keys = ObsVecKeys::get();
    afw::table::BaseCatalog catalog(keys.schema);
    catalog.reserve(obsVec.size());
    for (ObsVec::const_iterator iter = obsVec.begin(); iter != obsVec.end(); ++iter) {
        Obs const & obs = **iter;
        PTR(afw::table::BaseRecord) record = catalog.addNew();
        record->set(keys.id, obs.id);
        record->set(keys.ra, obs.ra);
        record->set(keys.dec, obs.dec);
        record->set(keys.xi, obs.xi);
        record->set(keys.eta, obs.eta);
        record->set(keys.x, obs.x);
        record->set(keys.y, obs.y);
        record->set(keys.u, obs.u);
        record->set(keys.v, obs.v);
        record->set(keys.u0, obs.u0);
        record->set(keys.v0, obs.v0);
        record->set(keys.iexp, obs.iexp);
        record->set(keys.ichip, obs.ichip);
        record->set(keys.jexp, obs.jexp);
        record->set(keys.jchip, obs.jchip);
    }
    catalog.writeFits(filename);
}

}}} // namespace lsst::meas::mosaic
