#pragma once    

#include "conditionals.hpp"

#include <algorithm>
#include <vector>
#include "conditionals.hpp"
#include "datatype.hpp"
#include "qmatrix.hpp"
#include "partition.hpp"
#include "asrv.hpp"
#include "libhmsbeagle/beagle.h"
#include <boost/format.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <Eigen/Dense>

namespace lorad {
    
    class Likelihood;

    class Model { 
        
        friend class Likelihood;

        public:
            typedef std::vector<ASRV::SharedPtr>      asrv_vect_t;
            typedef std::vector<QMatrix::SharedPtr>   qmat_vect_t;
            typedef std::vector<unsigned>             subset_sizes_t;
            typedef std::vector<DataType>             subset_datatype_t;
            typedef std::vector<double>               subset_relrate_vect_t;
            typedef std::vector<QMatrix::SharedPtr>   state_freq_params_t;
            typedef std::vector<QMatrix::SharedPtr>   exchangeability_params_t;
            typedef std::vector<QMatrix::SharedPtr>   omega_params_t;
#if defined(HOLDER_ETAL_PRIOR)
            typedef std::vector<ASRV::SharedPtr>      shape_params_t;
#else
            typedef std::vector<ASRV::SharedPtr>      ratevar_params_t;
#endif
            typedef std::vector<ASRV::SharedPtr>      pinvar_params_t;
            typedef boost::shared_ptr<Model>          SharedPtr;
        
                                        Model();
                                        ~Model();

            void                        activate();
            void                        inactivate();
            
            std::string                 describeModel();
            
            void                        setModelToSampledPoint(unsigned i);
            void                        setSampledSubsetRelRates(unsigned i);
            void                        setSampledExchangeabilities(unsigned subset, unsigned i);
            void                        setSampledStateFreqs(unsigned subset, unsigned i);
            void                        setSampledOmega(unsigned subset, unsigned i);
            void                        setSampledPinvar(unsigned subset, unsigned i);
            void                        setSampledShape(unsigned subset, unsigned i);
            void                        setSampledRateVar(unsigned subset, unsigned i);
            

            void                        setSubsetDataTypes(const subset_datatype_t & datatype_vect);

#if defined(HOLDER_ETAL_PRIOR)
            void                        setSubsetShape(ASRV::shape_ptr_t shape, unsigned subset, bool fixed);
#else
            void                        setSubsetRateVar(ASRV::ratevar_ptr_t ratevar, unsigned subset, bool fixed);
#endif
            void                        setSubsetPinvar(ASRV::pinvar_ptr_t pinvar, unsigned subset, bool fixed);
            void                        setSubsetExchangeabilities(QMatrix::freq_xchg_ptr_t exchangeabilities, unsigned subset, bool fixed);
            void                        setSubsetStateFreqs(QMatrix::freq_xchg_ptr_t state_frequencies, unsigned subset, bool fixed);
            void                        setSubsetOmega(QMatrix::omega_ptr_t omega, unsigned subset, bool fixed);

            void                        setSubsetRelRates(subset_relrate_vect_t & relrates, bool fixed);
            subset_relrate_vect_t &     getSubsetRelRates();
            bool                        isFixedSubsetRelRates() const;
            double                      calcNormalizingConstantForSubsetRelRates() const;

            void                        setSubsetStateFreqRefDistParams(QMatrix::freq_xchg_ptr_t freq_refdist_params, unsigned subset);
            void                        setSubsetExchangeabilitiesRefDistParams(QMatrix::freq_xchg_ptr_t xchg_refdist_params, unsigned subset);
            void                        setSubsetPinvarRefDistParams(ASRV::pinvar_refdist_ptr_t pinvar_refdist_params, unsigned subset);
#if defined(HOLDER_ETAL_PRIOR)
            void                        setSubsetShapeRefDistParams(ASRV::shape_refdist_ptr_t shape_refdist_params, unsigned subset);
            void                        setEdgeLenRefDistParams(std::vector<double> & edgelen_refdist);
            std::vector<double>         getEdgeLenRefDistParamsVect() const;
#else
            void                        setSubsetRateVarRefDistParams(ASRV::ratevar_refdist_ptr_t ratevar_refdist_params, unsigned subset);
            void                        setTreeLengthRefDistParams(std::vector<double> & treelen_refdist);
            std::vector<double>         getTreeLengthRefDistParamsVect() const;
            void                        setEdgeProportionsRefDistParams(std::vector<double> & edgeprops_refdist);
            std::vector<double>         getEdgeProportionsRefDistParamsVect() const;
#endif
            void                        sampleParams();

            std::string                 saveReferenceDistributions(Partition::SharedPtr partition);
            std::string                 calcReferenceDistributions(Partition::SharedPtr partition, std::map<std::string, std::vector<double> > & refdist_map);
            
            void                        setSubsetRelRatesRefDistParams(std::vector<double> refdist_params);
            std::vector<double>         getSubsetRelRatesRefDistParamsVect();
            std::string                 calcBetaRefDist(std::string title, std::string subset_name, std::vector<double> & vect, std::vector<double> & params);
            std::string                 calcGammaRefDist(std::string title, std::string subset_name, std::vector<double> & vect, std::vector<double> & params);
            std::string                 calcDirichletRefDist(std::string title, std::string subset_name, std::vector< QMatrix::freq_xchg_t > & vect, std::vector<double> & params, bool relrates = false);

            void                        setTreeIndex(unsigned i, bool fixed);
            unsigned                    getTreeIndex() const;
            bool                        isFixedTree() const;
            
            void                        setTopologyPriorOptions(bool allow_polytomies, bool resclass, double C); 
            bool                        isResolutionClassTopologyPrior() const;
            double                      getTopologyPriorC() const;
            bool                        isAllowPolytomies() const; 

            unsigned                    getNumSubsets() const;
            unsigned                    getNumSites() const;
            unsigned                    getSubsetNumSites(unsigned subset) const;
            const QMatrix &             getQMatrix(unsigned subset) const;
            const ASRV &                getASRV(unsigned subset) const;

            void                        setSubsetIsInvarModel(bool is_invar, unsigned subset);
            bool                        getSubsetIsInvarModel(unsigned subset) const;

            void                        setSubsetSizes(const subset_sizes_t nsites_vect);
            subset_sizes_t &            getSubsetSizes();

            void                        setSubsetNumPatterns(const subset_sizes_t npatterns_vect);
            unsigned                    getSubsetNumPatterns(unsigned subset) const;

            void                        setSubsetNumCateg(unsigned ncateg, unsigned subset);
            unsigned                    getSubsetNumCateg(unsigned subset) const;

            state_freq_params_t &       getStateFreqParams();
//TODO: why is this commented out?
//            void                        setStateFreqRefDistParams(std::vector<double> refdist_params);
//            std::vector<double>         getStateFreqRefDistParams();

            exchangeability_params_t &  getExchangeabilityParams();
//            void                        setExchangeabilityRefDistParams(std::vector<double> refdist_params);
//            std::vector<double>         getExchangeabilityRefDistParams();

            omega_params_t &            getOmegaParams();
//            void                        setOmegaRefDistParams(std::vector<double> refdist_params);
//            std::vector<double>         getOmegaRefDistParams();

#if defined(HOLDER_ETAL_PRIOR)
            shape_params_t &            getShapeParams();
//            void                        setShapeRefDistParams(std::vector<double> refdist_params);
//            std::vector<double>         getShapeRefDistParams();
#else
            ratevar_params_t &          getRateVarParams();
//            void                        setRateVarRefDistParams(std::vector<double> refdist_params);
//            std::vector<double>         getRateVarRefDistParams();
#endif

            pinvar_params_t &           getPinvarParams();
//            void                        setPinvarRefDistParams(std::vector<double> refdist_params);
//            std::vector<double>         getPinvarRefDistParams();
        
            int                         setBeagleEigenDecomposition(int beagle_instance, unsigned subset, unsigned instance_subset);
            int                         setBeagleStateFrequencies(int beagle_instance, unsigned subset, unsigned instance_subset);
            int                         setBeagleAmongSiteRateVariationRates(int beagle_instance, unsigned subset, unsigned instance_subset);
            int                         setBeagleAmongSiteRateVariationProbs(int beagle_instance, unsigned subset, unsigned instance_subset);

            std::string                 paramNamesAsString(std::string sep, bool logscale) const;
            std::string                 paramValuesAsString(std::string sep, bool logscale, unsigned precision = 9) const;

            void                        saveParamNames(std::vector<std::string> & param_name_vect) const;
            double                      logRatioTransform(std::vector<double> & param_vect) const;
            double                      logRatioUntransform(std::vector<double> & param_vect) const;
            double                      logTransformParameters(std::vector<double> & param_vect) const;
            double                      setParametersFromLogTransformed(Eigen::VectorXd & param_vect, unsigned first, unsigned nparams);

        private:
        
            void                        clear();
        
            unsigned                    _num_subsets;
            unsigned                    _num_sites;
            subset_sizes_t              _subset_sizes;
            subset_sizes_t              _subset_npatterns;
            subset_datatype_t           _subset_datatypes;
            qmat_vect_t                 _qmatrix;
            asrv_vect_t                 _asrv;
        
            unsigned                    _tree_index;
            bool                        _tree_fixed;
            
            bool                        _allow_polytomies; 
            bool                        _resolution_class_prior;
            double                      _topo_prior_C; 

            bool                        _subset_relrates_fixed;
            subset_relrate_vect_t       _subset_relrates;
        
            state_freq_params_t         _state_freq_params;
            exchangeability_params_t    _exchangeability_params;
            omega_params_t              _omega_params;
#if defined(HOLDER_ETAL_PRIOR)
            shape_params_t              _shape_params;
#else
            ratevar_params_t            _ratevar_params;
#endif
            pinvar_params_t             _pinvar_params;
            
            std::vector<double>         _state_freq_refdist_params;
            std::vector<double>         _exchangeability_refdist_params;
            std::vector<double>         _omega_refdist_params;
            std::vector<double>         _pinvar_refdist_params;
            std::vector<double>         _subset_relrates_refdist_params;
            std::vector<double>         _treelen_refdist_params;
#if defined(HOLDER_ETAL_PRIOR)
            std::vector<double>         _shape_refdist_params;
            std::vector<double>         _edgelen_refdist_params;
#else
            std::vector<double>         _ratevar_refdist_params;
            std::vector<double>         _edgeprops_refdist_params;
#endif
            
            std::vector< QMatrix::freq_xchg_t>                      _sampled_subset_relrates;
            std::map<unsigned, std::vector<QMatrix::freq_xchg_t> >  _sampled_exchangeabilities;
            std::map<unsigned, std::vector<QMatrix::freq_xchg_t> >  _sampled_state_freqs;
            std::map<unsigned, std::vector<double> >                _sampled_omegas;
            std::map<unsigned, std::vector<double> >                _sampled_pinvars;
#if defined(HOLDER_ETAL_PRIOR)
            std::map<unsigned, std::vector<double> >                _sampled_shapes;
#else
            std::map<unsigned, std::vector<double> >                _sampled_ratevars;
#endif
        };
    
    
    inline Model::Model() {
        clear();
    }

    inline Model::~Model() {
    }

    inline void Model::clear() {    
        _state_freq_refdist_params.clear();
        _exchangeability_refdist_params.clear();
        _omega_refdist_params.clear();
        _pinvar_refdist_params.clear();
        _subset_relrates_refdist_params.clear();
        _treelen_refdist_params.clear();
#if defined(HOLDER_ETAL_PRIOR)
        _shape_refdist_params.clear();
        _edgelen_refdist_params.clear();
#else
        _ratevar_refdist_params.clear();
        _edgeprops_refdist_params.clear();
#endif
        _num_subsets = 0;
        _num_sites = 0;
        _tree_index = 0;
        _tree_fixed = false;
        _allow_polytomies = true;
        _resolution_class_prior = true; 
        _topo_prior_C = 1.0; 
        _subset_relrates_fixed = false;
        _subset_relrates.clear();
        _subset_sizes.clear();
        _subset_npatterns.clear();
        _subset_datatypes.clear();
        _qmatrix.clear();
        _asrv.clear();
    }   
    
    inline std::string Model::describeModel() {
        // Creates summary such as following and returns as a string:
        //
        // Partition information:
        //
        //          data subset           1           2           3
        //    -----------------------------------------------------
        //           num. sites          20          20          20
        //        num. patterns           7           5          17
        //          num. states           4           4           4
        //      rate categories           4           1           4
        //
        // Parameter linkage:
        //
        //          data subset           1           2           3
        //    -----------------------------------------------------
        //          state freqs           1           1           1
        //    exchangeabilities           1           1           2
        //        rate variance           1           2           3
        //               pinvar           1           2           -
        
        // Start with empty parameter vectors
        _state_freq_params.clear();
        _exchangeability_params.clear();
        _omega_params.clear();
#if defined(HOLDER_ETAL_PRIOR)
        _shape_params.clear();
#else
        _ratevar_params.clear();
#endif
        _pinvar_params.clear();

        // Sets used to determine which parameters are linked across subsets
        std::set<double *> freqset;
        std::set<double *> xchgset;
        std::set<double *> omegaset;
#if defined(HOLDER_ETAL_PRIOR)
        std::set<double *> shapeset;
#else
        std::set<double *> ratevarset;
#endif
        std::set<double *> pinvarset;
        std::set<double *> relrateset;

        // Vectors of pointers to distinct parameters
        std::vector<double *> unique_freq;
        std::vector<double *> unique_xchg;
        std::vector<double *> unique_omega;
#if defined(HOLDER_ETAL_PRIOR)
        std::vector<double *> unique_shape;
#else
        std::vector<double *> unique_ratevar;
#endif
        std::vector<double *> unique_pinvar;
        std::vector<double *> unique_relrate;

        // Map for storing strings that will contain the information for each row
        std::map<std::string, std::string> ss = {
            {"subset",    ""},
            {"dashes",    ""},
            {"freqs",     ""},
            {"xchg",      ""},
            {"omega",     ""},
            {"ratevar",   ""},
            {"pinvar",    ""},
            {"ncateg",    ""},
            {"nsites",    ""},
            {"npatterns", ""},
            {"nstates",   ""}
        };
        
        // Ensure that the subset relative rates are fixed if there is only one
        // subset; otherwise the subset relative rates will be added to the list
        // of free parameters that are updated, which makes no sense in this case
        if (_num_subsets == 1)
            _subset_relrates_fixed = true;

        // Loop through subsets, building up rows as we go
        for (unsigned i = 0; i < _num_subsets; i++) {
#if defined(HOLDER_ETAL_PRIOR)
            // Ensure that for subsets in which the number of rate categories is 1 that
            // the gamma shape is fixed; otherwise the gamma shape will
            // be added to the list of free parameters that are updated, which makes
            // no sense in this case
            if (_asrv[i]->getNumCateg() == 1) {
                _asrv[i]->fixShape(true);
            }
#else
            // Ensure that for subsets in which the number of rate categories is 1 that
            // the gamma rate variance is fixed; otherwise the gamma rate variance will
            // be added to the list of free parameters that are updated, which makes
            // no sense in this case
            if (_asrv[i]->getNumCateg() == 1) {
                _asrv[i]->fixRateVar(true);
            }
#endif
        
            unsigned index;
            ss["subset"] += boost::str(boost::format("%12d") % (i+1));
            ss["dashes"] += "------------";

            // Determine whether state freqs are unique for this subset
            QMatrix::freq_xchg_ptr_t pfreq = _qmatrix[i]->getStateFreqsSharedPtr();
            QMatrix::freq_xchg_t & freq = *pfreq;
            double * freq_addr = &freq[0];
            auto f = freqset.insert(freq_addr);
            if (f.second) {
                unique_freq.push_back(freq_addr);
                if (!_qmatrix[i]->isFixedStateFreqs())
                    _state_freq_params.push_back(_qmatrix[i]);
                index = (unsigned)unique_freq.size();
            }
            else {
                auto iter = std::find(unique_freq.begin(), unique_freq.end(), freq_addr);
                index = (unsigned)std::distance(unique_freq.begin(), iter) + 1;
            }
            ss["freqs"] += boost::str(boost::format("%12d") % index);

            // Determine whether exchangeabilities are unique for this subset
            if (_subset_datatypes[i].isNucleotide()) {
                QMatrix::freq_xchg_ptr_t pxchg = _qmatrix[i]->getExchangeabilitiesSharedPtr();
                QMatrix::freq_xchg_t & xchg = *pxchg;
                double * xchg_addr = &xchg[0];
                auto x = xchgset.insert(xchg_addr);
                if (x.second) {
                    unique_xchg.push_back(xchg_addr);
                    if (!_qmatrix[i]->isFixedExchangeabilities())
                        _exchangeability_params.push_back(_qmatrix[i]);
                    index = (unsigned)unique_xchg.size();
                }
                else {
                    auto iter = std::find(unique_xchg.begin(), unique_xchg.end(), xchg_addr);
                    index = (unsigned)std::distance(unique_xchg.begin(), iter) + 1;
                }
                ss["xchg"] += boost::str(boost::format("%12d") % index);
            }
            else {
                ss["xchg"] += boost::str(boost::format("%12s") % "-");
            }
            
            // Determine whether omega is unique for this subset
            if (_subset_datatypes[i].isCodon()) {
                QMatrix::omega_ptr_t pomega = _qmatrix[i]->getOmegaSharedPtr();
                QMatrix::omega_t omegavalue = *pomega;
                double * omega_addr = &omegavalue;
                auto o = omegaset.insert(omega_addr);
                if (o.second) {
                    unique_omega.push_back(omega_addr);
                    if (!_qmatrix[i]->isFixedOmega())
                        _omega_params.push_back(_qmatrix[i]);
                    index = (unsigned)unique_omega.size();
                }
                else {
                    auto iter = std::find(unique_omega.begin(), unique_omega.end(), omega_addr);
                    index = (unsigned)std::distance(unique_omega.begin(), iter) + 1;
                }
                ss["omega"] += boost::str(boost::format("%12d") % index);
            }
            else {
                ss["omega"] += boost::str(boost::format("%12s") % "-");
            }

#if defined(HOLDER_ETAL_PRIOR)
            // Determine whether gamma shape is unique for this subset
            ASRV::shape_ptr_t pshape = _asrv[i]->getShapeSharedPtr();
            double & shape = *pshape;
            double * shape_addr = &shape;
            auto r = shapeset.insert(shape_addr);
            if (r.second) {
                unique_shape.push_back(shape_addr);
                if (!_asrv[i]->isFixedShape())
                    _shape_params.push_back(_asrv[i]);
                index = (unsigned)unique_shape.size();
            }
            else {
                auto iter = std::find(unique_shape.begin(), unique_shape.end(), shape_addr);
                index = (unsigned)std::distance(unique_shape.begin(), iter) + 1;
            }
            ss["shape"] += boost::str(boost::format("%12d") % index);
#else
            // Determine whether rate variance is unique for this subset
            ASRV::ratevar_ptr_t pratevar = _asrv[i]->getRateVarSharedPtr();
            double & ratevar = *pratevar;
            double * ratevar_addr = &ratevar;
            auto r = ratevarset.insert(ratevar_addr);
            if (r.second) {
                unique_ratevar.push_back(ratevar_addr);
                if (!_asrv[i]->isFixedRateVar())
                    _ratevar_params.push_back(_asrv[i]);
                index = (unsigned)unique_ratevar.size();
            }
            else {
                auto iter = std::find(unique_ratevar.begin(), unique_ratevar.end(), ratevar_addr);
                index = (unsigned)std::distance(unique_ratevar.begin(), iter) + 1;
            }
            ss["ratevar"] += boost::str(boost::format("%12d") % index);
#endif
            
            // Determine whether pinvar is unique for this subset
            if (_asrv[i]->getIsInvarModel()) {
                ASRV::pinvar_ptr_t ppinvar = _asrv[i]->getPinvarSharedPtr();
                double & pinvar = *ppinvar;
                double * pinvar_addr = &pinvar;
                auto r = pinvarset.insert(pinvar_addr);
                if (r.second) {
                    unique_pinvar.push_back(pinvar_addr);
                    if (!_asrv[i]->isFixedPinvar())
                        _pinvar_params.push_back(_asrv[i]);
                    index = (unsigned)unique_pinvar.size();
                }
                else {
                    auto iter = std::find(unique_pinvar.begin(), unique_pinvar.end(), pinvar_addr);
                    index = (unsigned)std::distance(unique_pinvar.begin(), iter) + 1;
                }
                ss["pinvar"] += boost::str(boost::format("%12d") % index);
            }
            else {
                ss["pinvar"] += boost::str(boost::format("%12s") % "-");
            }
            
            // Save number of rate categories for this subset
            ss["ncateg"] += boost::str(boost::format("%12d") % _asrv[i]->getNumCateg());

            // Save number of sites for this subset
            ss["nsites"] += boost::str(boost::format("%12d") % _subset_sizes[i]);

            // Save number of patterns for this subset
            ss["npatterns"] += boost::str(boost::format("%12d") % _subset_npatterns[i]);

            // Save number of states for this subset
            if (_subset_datatypes.size() == _num_subsets)
                ss["nstates"] += boost::str(boost::format("%12d") % _subset_datatypes[i].getNumStates());
            else
                ss["nstates"] += boost::str(boost::format("%12s") % "?");

        }
        std::string s = "Partition information:\n\n";
        
        s += boost::str(boost::format("%20s%s\n") % "data subset" % ss["subset"]);
        s += boost::str(boost::format("%20s%s\n") % "-----------------" % ss["dashes"]);
        s += boost::str(boost::format("%20s%s\n") % "num. sites" % ss["nsites"]);
        s += boost::str(boost::format("%20s%s\n") % "num. patterns" % ss["npatterns"]);
        s += boost::str(boost::format("%20s%s\n") % "num. states" % ss["nstates"]);
        s += boost::str(boost::format("%20s%s\n") % "rate categories" % ss["ncateg"]);

        s += "\nParameter linkage:\n\n";
        
        s += boost::str(boost::format("%20s%s\n") % "data subset" % ss["subset"]);
        s += boost::str(boost::format("%20s%s\n") % "-----------------" % ss["dashes"]);
        s += boost::str(boost::format("%20s%s\n") % "state freqs" % ss["freqs"]);
        s += boost::str(boost::format("%20s%s\n") % "exchangeabilities" % ss["xchg"]);
        s += boost::str(boost::format("%20s%s\n") % "omega" % ss["omega"]);
        s += boost::str(boost::format("%20s%s\n") % "rate variance" % ss["ratevar"]);
        s += boost::str(boost::format("%20s%s\n") % "pinvar" % ss["pinvar"]);
        
        s += "\nParameter values for each subset:\n";

        s += "\n  relative rate:\n";
        for (unsigned i = 0; i < _num_subsets; i++) {
            s += boost::str(boost::format("  %12d: %g\n") % (i+1) % _subset_relrates[i]);
        }
        
        s += "\n  state freqs:\n";
        for (unsigned i = 0; i < _num_subsets; i++) {
            QMatrix::freq_xchg_t & freqs = *(_qmatrix[i]->getStateFreqsSharedPtr());
            std::vector<std::string> freqs_as_strings(freqs.size());
            std::transform(freqs.begin(), freqs.end(), freqs_as_strings.begin(), [](double freq) {return boost::str(boost::format("%g") % freq);});
            std::string tmp = boost::algorithm::join(freqs_as_strings, ",");
            s += boost::str(boost::format("  %12d: (%s)\n") % (i+1) % tmp);
        }

        s += "\n  exchangeabilities:\n";
        for (unsigned i = 0; i < _num_subsets; i++) {
            if (_subset_datatypes[i].isNucleotide()) {
                QMatrix::freq_xchg_t & xchg = *(_qmatrix[i]->getExchangeabilitiesSharedPtr());
                std::vector<std::string> xchg_as_strings(xchg.size());
                std::transform(xchg.begin(), xchg.end(), xchg_as_strings.begin(), [](double x) {return boost::str(boost::format("%g") % x);});
                std::string tmp = boost::algorithm::join(xchg_as_strings, ",");
                s += boost::str(boost::format("  %12d: (%s)\n") % (i+1) % tmp);
            }
            else {
                s += boost::str(boost::format("  %12d: -\n") % (i+1));
            }
        }

        s += "\n  omega:\n";
        for (unsigned i = 0; i < _num_subsets; i++) {
            if (_subset_datatypes[i].isCodon()) {
                double omega = *(_qmatrix[i]->getOmegaSharedPtr());
                s += boost::str(boost::format("  %12d: %g\n") % (i+1) % omega);
            }
            else {
                s += boost::str(boost::format("  %12d: -\n") % (i+1));
            }
        }

#if defined(HOLDER_ETAL_PRIOR)
        s += "\n  gamma shape:\n";
        for (unsigned i = 0; i < _num_subsets; i++) {
            if (_asrv[i]->getNumCateg() > 1) {
                double shape = *(_asrv[i]->getShapeSharedPtr());
                s += boost::str(boost::format("  %12d: %g\n") % (i+1) % shape);
            }
            else
                s += boost::str(boost::format("  %12d: -\n") % (i+1));
        }
#else
        s += "\n  rate variance:\n";
        for (unsigned i = 0; i < _num_subsets; i++) {
            if (_asrv[i]->getNumCateg() > 1) {
                double ratevar = *(_asrv[i]->getRateVarSharedPtr());
                s += boost::str(boost::format("  %12d: %g\n") % (i+1) % ratevar);
            }
            else
                s += boost::str(boost::format("  %12d: -\n") % (i+1));
        }
#endif

        s += "\n  pinvar:\n";
        for (unsigned i = 0; i < _num_subsets; i++) {
            double pinvar = *(_asrv[i]->getPinvarSharedPtr());
            bool is_invar_model = _asrv[i]->getIsInvarModel();
            if (is_invar_model)
                s += boost::str(boost::format("  %12d: %g\n") % (i+1) % pinvar);
            else
                s += boost::str(boost::format("  %12d: -\n") % (i+1));
        }

        return s;
    }
    
    inline unsigned Model::getSubsetNumPatterns(unsigned subset) const {
        assert(subset < _num_subsets);
        return _subset_npatterns[subset];
    }        
    
    inline unsigned Model::getSubsetNumSites(unsigned subset) const {
        assert(subset < _num_subsets);
        return _subset_sizes[subset];
    }
    
    inline unsigned Model::getNumSites() const {
        return _num_sites;
    }
    
    inline unsigned Model::getNumSubsets() const {
        return _num_subsets;
    }
    
    inline unsigned Model::getSubsetNumCateg(unsigned subset) const {
        assert(subset < _num_subsets);
        assert(_asrv.size() == _num_subsets);
        assert(_asrv[subset]);
        return _asrv[subset]->getNumCateg();
    }
    
    inline bool Model::getSubsetIsInvarModel(unsigned subset) const {
        assert(subset < _num_subsets);
        assert(_asrv.size() == _num_subsets);
        assert(_asrv[subset]);
        return _asrv[subset]->getIsInvarModel();
    }
    
    inline const QMatrix & Model::getQMatrix(unsigned subset) const {
        assert(subset < _num_subsets);
        return *(_qmatrix[subset]);
    }
    
    inline const ASRV & Model::getASRV(unsigned subset) const {
        assert(subset < _num_subsets);
        return *(_asrv[subset]);
    }
    
    inline Model::state_freq_params_t & Model::getStateFreqParams() {
        return _state_freq_params;
    }
    
    inline Model::exchangeability_params_t & Model::getExchangeabilityParams() {
        return _exchangeability_params;
    }

    inline Model::omega_params_t & Model::getOmegaParams() {
        return _omega_params;
    }
    
#if defined(HOLDER_ETAL_PRIOR)
    inline Model::shape_params_t & Model::getShapeParams() {
        return _shape_params;
    }
#else
    inline Model::ratevar_params_t & Model::getRateVarParams() {
        return _ratevar_params;
    }
#endif
    
    inline Model::pinvar_params_t & Model::getPinvarParams() {
        return _pinvar_params;
    }
    
    inline double Model::calcNormalizingConstantForSubsetRelRates() const {
        // normalize _relrates so that expected relative rate across subsets equals 1.0
        double normalizing_constant = 0.0;
        for (unsigned s = 0; s < _num_subsets; s++) {
            normalizing_constant += _subset_sizes[s]*_subset_relrates[s]/_num_sites;
        }
        return normalizing_constant;
    }

    inline Model::subset_sizes_t & Model::getSubsetSizes() {
        return _subset_sizes;
    }
    
    inline void Model::setSubsetSizes(const subset_sizes_t nsites_vect) {
        assert(nsites_vect.size() == _num_subsets);
        _subset_sizes.resize(_num_subsets);
        std::copy(nsites_vect.begin(), nsites_vect.end(), _subset_sizes.begin());
        _num_sites = std::accumulate(_subset_sizes.begin(), _subset_sizes.end(), 0);
    }

    inline void Model::setSubsetNumPatterns(const subset_sizes_t npatterns_vect) {
        assert(npatterns_vect.size() == _num_subsets);
        _subset_npatterns.resize(_num_subsets);
        std::copy(npatterns_vect.begin(), npatterns_vect.end(), _subset_npatterns.begin());
    }
    
    inline void Model::setSubsetDataTypes(const subset_datatype_t & datatype_vect) {
        _num_subsets = (unsigned)datatype_vect.size();

        _qmatrix.clear();
        _qmatrix.resize(_num_subsets);

        _asrv.clear();
        _asrv.resize(_num_subsets);

        _subset_datatypes.resize(_num_subsets);
        std::copy(datatype_vect.begin(), datatype_vect.end(), _subset_datatypes.begin());
        
		_subset_relrates.assign(_num_subsets, 1.0);
        
        for (unsigned s = 0; s < _num_subsets; s++) {
            _asrv[s].reset(new ASRV());
            if (_subset_datatypes[s].isNucleotide())
                _qmatrix[s].reset(new QMatrixNucleotide());
            else if (_subset_datatypes[s].isCodon()) {
                GeneticCode::SharedPtr gcptr = _subset_datatypes[s].getGeneticCode();
                _qmatrix[s].reset(new QMatrixCodon(gcptr));
                }
            else
                throw XLorad(boost::format("Only nucleotide or codon data allowed in this version, you specified data type \"%s\" for subset %d") % _subset_datatypes[s].getDataTypeAsString() % (s+1));
        }
    }

    inline void Model::setSubsetNumCateg(unsigned ncateg, unsigned subset) {
        assert(subset < _num_subsets);
        if (ncateg < 1) {
            throw XLorad(boost::str(boost::format("number of categories used for among-site rate variation must be greater than zero but the value %d was supplied") % ncateg));
        }
        _asrv[subset]->setNumCateg(ncateg);
    }
    
#if defined(HOLDER_ETAL_PRIOR)
    inline void Model::setSubsetShape(ASRV::shape_ptr_t shape, unsigned subset, bool fixed) {
        assert(subset < _num_subsets);
        assert(shape);
        if (*shape < 0.0)
            throw XLorad(boost::str(boost::format("gamma shape must be greater than or equal to zero but the value %.5f was supplied") % *shape));
        _asrv[subset]->setShapeSharedPtr(shape);
        _asrv[subset]->fixShape(fixed);
    }
#else
    inline void Model::setSubsetRateVar(ASRV::ratevar_ptr_t ratevar, unsigned subset, bool fixed) {
        assert(subset < _num_subsets);
        assert(ratevar);
        if (*ratevar < 0.0)
            throw XLorad(boost::str(boost::format("rate variance must be greater than or equal to zero but the value %.5f was supplied") % *ratevar));
        _asrv[subset]->setRateVarSharedPtr(ratevar);
        _asrv[subset]->fixRateVar(fixed);
    }
#endif
    
    inline void Model::setSubsetPinvar(ASRV::pinvar_ptr_t pinvar, unsigned subset, bool fixed) {
        assert(subset < _num_subsets);
        assert(pinvar);
        if (*pinvar < 0.0)
            throw XLorad(boost::str(boost::format("proportion of invariable sites must be greater than or equal to zero but the value %.5f was supplied") % *pinvar));
        if (*pinvar >= 1.0)
            throw XLorad(boost::str(boost::format("proportion of invariable sites must be less than one but the value %.5f was supplied") % *pinvar));
        _asrv[subset]->setPinvarSharedPtr(pinvar);
        _asrv[subset]->fixPinvar(fixed);
    }
    
    inline void Model::setSubsetIsInvarModel(bool is_invar, unsigned subset) {
        assert(subset < _num_subsets);
        _asrv[subset]->setIsInvarModel(is_invar);
    }
    
    inline void Model::setSubsetExchangeabilities(QMatrix::freq_xchg_ptr_t exchangeabilities, unsigned subset, bool fixed) {
        assert(subset < _num_subsets);
        if (!_subset_datatypes[subset].isCodon()) {
            double first_xchg = (*exchangeabilities)[0];
            if (first_xchg == -1)
                _qmatrix[subset]->setEqualExchangeabilities(exchangeabilities);
            else
                _qmatrix[subset]->setExchangeabilitiesSharedPtr(exchangeabilities);
            _qmatrix[subset]->fixExchangeabilities(fixed);
        }
        //TODO: there doesn't seem to be code supporting a codon model here
    }
    
    inline void Model::setSubsetStateFreqs(QMatrix::freq_xchg_ptr_t state_frequencies, unsigned subset, bool fixed) {
        assert(subset < _num_subsets);
        double first_freq = (*state_frequencies)[0];
        if (first_freq == -1)
            _qmatrix[subset]->setEqualStateFreqs(state_frequencies);
        else
            _qmatrix[subset]->setStateFreqsSharedPtr(state_frequencies);
        _qmatrix[subset]->fixStateFreqs(fixed);
    }
    
    inline void Model::setSubsetStateFreqRefDistParams(QMatrix::freq_xchg_ptr_t freq_refdist_params, unsigned subset) {
        assert(subset < _num_subsets);
        _qmatrix[subset]->setStateFreqRefDistParamsSharedPtr(freq_refdist_params);
    }
    
    inline void Model::setSubsetExchangeabilitiesRefDistParams(QMatrix::freq_xchg_ptr_t xchg_refdist_params, unsigned subset) {
        assert(subset < _num_subsets);
        _qmatrix[subset]->setExchangeabilityRefDistParamsSharedPtr(xchg_refdist_params);
    }
#if defined(HOLDER_ETAL_PRIOR)
    inline void Model::setSubsetShapeRefDistParams(ASRV::shape_refdist_ptr_t shape_refdist_params, unsigned subset) {
        assert(subset < _num_subsets);
        _asrv[subset]->setShapeRefDistParamsSharedPtr(shape_refdist_params);
    }
#else
    inline void Model::setSubsetRateVarRefDistParams(ASRV::ratevar_refdist_ptr_t ratevar_refdist_params, unsigned subset) {
        assert(subset < _num_subsets);
        _asrv[subset]->setRateVarRefDistParamsSharedPtr(ratevar_refdist_params);
    }
#endif
    inline void Model::setSubsetPinvarRefDistParams(ASRV::pinvar_refdist_ptr_t pinvar_refdist_params, unsigned subset) {
        assert(subset < _num_subsets);
        _asrv[subset]->setPinvarRefDistParamsSharedPtr(pinvar_refdist_params);
    }
    
    inline void Model::setSubsetOmega(QMatrix::omega_ptr_t omega, unsigned subset, bool fixed) {
        assert(subset < _num_subsets);
        assert(*omega > 0.0);
        if (_subset_datatypes[subset].isCodon()) {
            _qmatrix[subset]->setOmegaSharedPtr(omega);
            _qmatrix[subset]->fixOmega(fixed);
        }
    }
    
    inline void Model::activate() {
        for (auto q : _qmatrix)
            q->setActive(true);
    }

    inline void Model::inactivate() {
        for (auto q : _qmatrix)
            q->setActive(false);
    }

    inline int Model::setBeagleEigenDecomposition(int beagle_instance, unsigned subset, unsigned instance_subset) {
        assert(subset < _qmatrix.size());
        const double * pevec = _qmatrix[subset]->getEigenvectors();
        const double * pivec = _qmatrix[subset]->getInverseEigenvectors();
        const double * pival = _qmatrix[subset]->getEigenvalues();
        int code = beagleSetEigenDecomposition(
            beagle_instance,    // Instance number (input)
            instance_subset,    // Index of eigen-decomposition buffer (input)
            pevec,              // Flattened matrix (stateCount x stateCount) of eigen-vectors (input)
            pivec,              // Flattened matrix (stateCount x stateCount) of inverse-eigen- vectors (input)
            pival);             // Vector of eigenvalues

        return code;
    }
    
    inline int Model::setBeagleStateFrequencies(int beagle_instance, unsigned subset, unsigned instance_subset) {
        assert(subset < _qmatrix.size());
        const double * pfreq = _qmatrix[subset]->getStateFreqs();
        int code = beagleSetStateFrequencies(
             beagle_instance,   // Instance number (input)
             instance_subset,   // Index of state frequencies buffer (input)
             pfreq);            // State frequencies array (stateCount) (input)

        return code;
    }
    
    inline int Model::setBeagleAmongSiteRateVariationRates(int beagle_instance, unsigned subset, unsigned instance_subset) {
        assert(subset < _asrv.size());
        const double * prates = _asrv[subset]->getRates();
        int code = beagleSetCategoryRatesWithIndex(
            beagle_instance,    // Instance number (input)
            instance_subset,    // Index of category rates buffer (input)
            prates);            // Array containing categoryCount rate scalers (input)

        return code;
    }
    
    inline int Model::setBeagleAmongSiteRateVariationProbs(int beagle_instance, unsigned subset, unsigned instance_subset) {
        assert(subset < _asrv.size());
        const double * pprobs = _asrv[subset]->getProbs();
        int code = beagleSetCategoryWeights(
            beagle_instance,    // Instance number (input)
            instance_subset,    // Index of category weights buffer (input)
            pprobs);            // Category weights array (categoryCount) (input)

        return code;
    }
        
    inline void Model::saveParamNames(std::vector<std::string> & param_name_vect) const {
        // Almost the same as Model::paramNamesAsString, but Model::paramNamesAsString is for output
        // to parameter log files where redundancy does not matter (for example, TL can be output
        // in addition to all edge lengths). Model::saveParamNames, on the other hand, is for output
        // to a file used in marginal likelihood estimation, and thus does not have
        // columns for redundant parameters (e.g. no name for 6th exchangeability or 4th frequency)
        if (_num_subsets > 1) {
            for (unsigned i = 1; i <= _num_subsets - 1; ++i)
                param_name_vect.push_back(boost::str(boost::format("subsetrate-%d") % i));
        }
        for (unsigned k = 1; k <= _num_subsets; k++) {
            if (_subset_datatypes[k-1].isNucleotide()) {
                if (!_qmatrix[k-1]->isFixedExchangeabilities()) {
                    param_name_vect.push_back(boost::str(boost::format("xchg-%d-1") % k));
                    param_name_vect.push_back(boost::str(boost::format("xchg-%d-2") % k));
                    param_name_vect.push_back(boost::str(boost::format("xchg-%d-3") % k));
                    param_name_vect.push_back(boost::str(boost::format("xchg-%d-4") % k));
                    param_name_vect.push_back(boost::str(boost::format("xchg-%d-5") % k));
                }
                
                if (!_qmatrix[k-1]->isFixedStateFreqs()) {
                    param_name_vect.push_back(boost::str(boost::format("freq-%d-1") % k));
                    param_name_vect.push_back(boost::str(boost::format("freq-%d-2") % k));
                    param_name_vect.push_back(boost::str(boost::format("freq-%d-3") % k));
                }
            }
            else if (_subset_datatypes[k-1].isCodon()) {
                if (!_qmatrix[k-1]->isFixedOmega()) {
                    param_name_vect.push_back(boost::str(boost::format("omega-%d") % k));
                }
                
                if (!_qmatrix[k-1]->isFixedStateFreqs()) {
                    for (unsigned i = 1; i <= 60; ++i)
                        param_name_vect.push_back(boost::str(boost::format("freq-%d-%d") % k % i));
                }
            }
            if (_asrv[k-1]->getIsInvarModel()) {
                if (!_asrv[k-1]->isFixedPinvar()) {
                    param_name_vect.push_back(boost::str(boost::format("pinvar-%d") % k));
                }
            }
#if defined(HOLDER_ETAL_PRIOR)
            if (_asrv[k-1]->getNumCateg() > 1) {
                if (!_asrv[k-1]->isFixedShape()) {
                    param_name_vect.push_back(boost::str(boost::format("shape-%d") % k));
                }
            }
#else
            if (_asrv[k-1]->getNumCateg() > 1) {
                if (!_asrv[k-1]->isFixedRateVar()) {
                    param_name_vect.push_back(boost::str(boost::format("ratevar-%d") % k));
                }
            }
#endif
        }
    }
    
    inline std::string Model::paramNamesAsString(std::string sep, bool logscale) const {
        unsigned k;
        std::string s = "";
        if (_num_subsets > 1) {
            for (k = (logscale ? 1 : 0); k < _num_subsets; k++) {
                s += boost::str(boost::format("relrate-%d%s") % (k+1) % sep);
            }
        }
        for (k = 0; k < _num_subsets; k++) {
            if (_subset_datatypes[k].isNucleotide()) {
                if (!_qmatrix[k]->isFixedExchangeabilities()) {
                    if (logscale)
                        s += boost::str(boost::format("rAG-%d%srAT-%d%srCG-%d%srCT-%d%srGT-%d%s") % k % sep % k % sep % k % sep % k % sep % k % sep);
                    else
                        s += boost::str(boost::format("rAC-%d%srAG-%d%srAT-%d%srCG-%d%srCT-%d%srGT-%d%s") % k % sep % k % sep % k % sep % k % sep % k % sep % k % sep);
                }
                if (!_qmatrix[k]->isFixedStateFreqs()) {
                    if (logscale)
                        s += boost::str(boost::format("piC-%d%spiG-%d%spiT-%d%s") % k % sep % k % sep % k % sep);
                    else
                        s += boost::str(boost::format("piA-%d%spiC-%d%spiG-%d%spiT-%d%s") % k % sep % k % sep % k % sep % k % sep);
                }
            }
            else if (_subset_datatypes[k].isCodon()) {
                if (!_qmatrix[k]->isFixedOmega()) {
                    s += boost::str(boost::format("omega-%d%s") % k % sep);
                }
                if (!_qmatrix[k]->isFixedStateFreqs()) {
                    bool first = true;
                    for (std::string codon : _subset_datatypes[0].getGeneticCode()->_codons) {
                        if (!logscale || !first)
                            s += boost::str(boost::format("pi%s-%d%s") % codon % k % sep);
                        first = false;
                    }
                }
            }
            if (_asrv[k]->getIsInvarModel()) {
                if (!_asrv[k]->isFixedPinvar()) {
                    s += boost::str(boost::format("pinvar-%d%s") % k % sep);
                }
            }
#if defined(HOLDER_ETAL_PRIOR)
            if (_asrv[k]->getNumCateg() > 1) {
                if (!_asrv[k]->isFixedShape()) {
                    s += boost::str(boost::format("shape-%d%s") % k % sep);
                }
            }
#else
            if (_asrv[k]->getNumCateg() > 1) {
                if (!_asrv[k]->isFixedRateVar()) {
                    s += boost::str(boost::format("ratevar-%d%s") % k % sep);
                }
            }
#endif
        }
        return s;
    }

    inline std::string Model::paramValuesAsString(std::string sep, bool logscale, unsigned precision) const {
        std::string fmtonestr  = boost::str(boost::format("%%.%df%s") % precision % sep);
        boost::format fmtone   = boost::format(fmtonestr);

        std::string fmtthreestr = boost::str(boost::format("%%.%df%s%%.%df%s%%.%df%s") % precision % sep % precision % sep % precision % sep);
        boost::format fmtthree  = boost::format(fmtthreestr);

        std::string fmtfourstr = boost::str(boost::format("%%.%df%s%%.%df%s%%.%df%s%%.%df%s") % precision % sep % precision % sep % precision % sep % precision % sep);
        boost::format fmtfour  = boost::format(fmtfourstr);

        std::string fmtfivestr  = boost::str(boost::format("%%.%df%s%%.%df%s%%.%df%s%%.%df%s%%.%df%s") % precision % sep % precision % sep % precision % sep % precision % sep % precision % sep);
        boost::format fmtfive  = boost::format(fmtfivestr);
        
        std::string fmtsixstr  = boost::str(boost::format("%%.%df%s%%.%df%s%%.%df%s%%.%df%s%%.%df%s%%.%df%s") % precision % sep % precision % sep % precision % sep % precision % sep % precision % sep % precision % sep);
        boost::format fmtsix  = boost::format(fmtsixstr);

        unsigned k;
        std::string s = "";
        if (_num_subsets > 1) {
            for (k = (logscale ? 1 : 0); k < _num_subsets; k++) {
                s += boost::str(fmtone % _subset_relrates[k]);
            }
        }
        for (k = 0; k < _num_subsets; k++) {
            if (_subset_datatypes[k].isNucleotide()) {
                if (!_qmatrix[k]->isFixedExchangeabilities()) {
                    QMatrix::freq_xchg_t x = *_qmatrix[k]->getExchangeabilitiesSharedPtr();
                    if (logscale) {
                        assert(x[0] > 0.0);
                        assert(x[1] > 0.0);
                        assert(x[2] > 0.0);
                        assert(x[3] > 0.0);
                        assert(x[4] > 0.0);
                        assert(x[5] > 0.0);
                        double logAC = log(x[0]);
                        double logAG = log(x[1]) - logAC;
                        double logAT = log(x[2]) - logAC;
                        double logCG = log(x[3]) - logAC;
                        double logCT = log(x[4]) - logAC;
                        double logGT = log(x[5]) - logAC;
                        s += boost::str(fmtfive % logAG % logAT % logCG % logCT % logGT);
                    }
                    else
                        s += boost::str(fmtsix % x[0] % x[1] % x[2] % x[3] % x[4] % x[5]);
                }
                if (!_qmatrix[k]->isFixedStateFreqs()) {
                    QMatrix::freq_xchg_t f = *_qmatrix[k]->getStateFreqsSharedPtr();
                    if (logscale) {
                        assert(f[0] > 0.0);
                        assert(f[1] > 0.0);
                        assert(f[2] > 0.0);
                        assert(f[3] > 0.0);
                        double logA = log(f[0]);
                        double logC = log(f[1]) - logA;
                        double logG = log(f[2]) - logA;
                        double logT = log(f[3]) - logA;
                        s += boost::str(fmtthree % logC % logG % logT);
                    }
                    else
                        s += boost::str(fmtfour % f[0] % f[1] % f[2] % f[3]);
                }
            }
            else if (_subset_datatypes[k].isCodon()) {
                if (!_qmatrix[k]->isFixedOmega()) {
                    double omega = _qmatrix[k]->getOmega();
                    if (logscale) {
                        assert(omega > 0.0);
                        double logomega = log(omega);
                        s += boost::str(fmtone % logomega);
                    }
                    else
                        s += boost::str(fmtone % omega);
                }
                if (!_qmatrix[k]->isFixedStateFreqs()) {
                    QMatrix::freq_xchg_t f = *_qmatrix[k]->getStateFreqsSharedPtr();
                    double log_first = 0.0;
                    for (unsigned m = 0; m < _subset_datatypes[0].getNumStates(); m++) {
                        double state_freq = f[m];
                        if (logscale) {
                            assert(f[m] > 0.0);
                            if (m == 0)
                                log_first = log(f[m]);
                            else {
                                double logfreq = log(f[m]) - log_first;
                                s += boost::str(fmtone % logfreq);
                            }
                        }
                        else
                            s += boost::str(fmtone % state_freq);
                    }
                }
            }
            if (_asrv[k]->getIsInvarModel()) {
                if (!_asrv[k]->isFixedPinvar()) {
                    double pinv = _asrv[k]->getPinvar();
                    if (logscale) {
                        assert(pinv > 0.0);
                        double logpinv = log(pinv);
                        s += boost::str(fmtone % logpinv);
                    }
                    else
                        s += boost::str(fmtone % pinv);
                }
            }
#if defined(HOLDER_ETAL_PRIOR)
            if (_asrv[k]->getNumCateg() > 1) {
                if (!_asrv[k]->isFixedShape()) {
                    double shape = _asrv[k]->getShape();
                    if (logscale) {
                        assert(shape > 0);
                        double logshape = log(shape);
                        s += boost::str(fmtone % logshape);
                    }
                    else
                        s += boost::str(fmtone % shape);
                }
            }
#else
            if (_asrv[k]->getNumCateg() > 1) {
                if (!_asrv[k]->isFixedRateVar()) {
                    double ratevar = _asrv[k]->getRateVar();
                    if (logscale) {
                        assert(ratevar > 0);
                        double logratevar = log(ratevar);
                        s += boost::str(fmtone % logratevar);
                    }
                    else
                        s += boost::str(fmtone % ratevar);
                }
            }
#endif
        }
        return s;
    }
    
    inline void Model::setSubsetRelRates(subset_relrate_vect_t & relrates, bool fixed) {
        assert(_num_subsets > 0);
        assert(relrates.size() > 0);
        if (relrates[0] == -1)
            _subset_relrates.assign(_num_subsets, 1.0);
        else
            _subset_relrates.assign(relrates.begin(), relrates.end());
        _subset_relrates_fixed = fixed;
    }

    inline Model::subset_relrate_vect_t & Model::getSubsetRelRates() {
        return _subset_relrates;
    }
    
    inline bool Model::isFixedSubsetRelRates() const {
        return _subset_relrates_fixed;
    }
    
    inline void Model::setTreeIndex(unsigned i, bool fixed) {
        _tree_index = i;
        _tree_fixed = fixed;
    }

    inline unsigned Model::getTreeIndex() const {
        return _tree_index;
    }
    
    inline bool Model::isFixedTree() const {
        return _tree_fixed;
    }

    inline void Model::setTopologyPriorOptions(bool allow_polytomies, bool resclass, double C) {  
        _allow_polytomies       = allow_polytomies;
        _resolution_class_prior = resclass;
        _topo_prior_C           = C;
    }   

    inline bool Model::isAllowPolytomies() const {  
        return _allow_polytomies;
    }   

    inline bool Model::isResolutionClassTopologyPrior() const {  
        return _resolution_class_prior;
    }   

    inline double Model::getTopologyPriorC() const {  
        return _topo_prior_C;
    }   

    // Suppose param_vect = {a, b, c, d} and the sum of elements = 1.
    // Replaces param_vect with {log(b/a), log(c/a), log(d/a)}.
    // phi = b/a + c/a + d/a = (1-a)/a.
    inline double Model::logRatioTransform(std::vector<double> & param_vect) const {
        unsigned sz = (unsigned)param_vect.size();
        assert(sz > 0);
        std::vector<double> result_vect(sz - 1);
        double first = param_vect[0];
        double log_first = log(first);
        double log_jacobian = log_first;
        for (unsigned i = 1; i < sz; ++i) {
            double element = param_vect[i];
            double log_element = log(element);
            log_jacobian += log_element;
            result_vect[i-1] = log_element - log_first;
        }
        
        // result_vect = {log(b/a), log(c/a), log(d/a)}
        // log_jacobian = log(a) + log(b) + log(c) + log(d)
        param_vect.clear();
        param_vect.resize(sz - 1);
        std::copy(result_vect.begin(), result_vect.end(), param_vect.begin());
        return log_jacobian;
    }

    // Suppose param_vect = {log(b/a), log(c/a), log(d/a)}.
    // If phi = b/a + c/a + d/a = (1-a)/a, then a = 1/(1+phi).
    // Replaces param_vect with {a, b, c, d}.
    inline double Model::logRatioUntransform(std::vector<double> & param_vect) const {
        unsigned sz = (unsigned)param_vect.size();
        assert(sz > 0);
        std::vector<double> result_vect(sz + 1);
        double phi = 0.0;
        result_vect[0] = 1.0;
        for (unsigned i = 1; i < sz + 1; ++i) {
            double r = exp(param_vect[i-1]);
            result_vect[i] = r;
            phi += r;
        }
        
        // result_vect = {1, b/a, c/a, d/a}
        // phi = b/a + c/a + d/a = (b+c+d)/a = (1-a)/a
        double log_jacobian = 0.0;
        for (unsigned i = 0; i < sz + 1; ++i) {
            result_vect[i] /= (1.0 + phi);
            log_jacobian -= log(result_vect[i]);  //POLMOD1 2022-01-05 (+=), 2022-01-30 (-=)
        }

        // result_vect = {a, b, c, d}
        // log_jacobian = -log(a) - log(b) - log(c) - log(d)
        param_vect.clear();
        param_vect.resize(sz + 1);
        std::copy(result_vect.begin(), result_vect.end(), param_vect.begin());
        return log_jacobian;
    }

    inline double Model::logTransformParameters(std::vector<double> & param_vect) const {
        unsigned k;
        double log_jacobian = 0.0;
        if (_num_subsets > 1 && !isFixedSubsetRelRates()) {
            // Suppose the subset relative rates are r1, r2, and r3
            // and the probabilities associated with these relative rates
            // are p1, p2, and p3. For example, p1 = p2 = p3 = 1/3 if
            // partitioning by codon position. In order to log-ratio-transform
            // (r1, r2, r3), we first must transform to a Dirichlet-distributed
            // random variable (y1,y2,y3) = (p1*r1, p2*r2, p3*r3)
            // (log jacobian = -log(p1)-log(p2)) and then perform
            // the log ratio transformation (log jacobian = log(p1*r1) + log(p2*r2) + log(p3*r3)). The total log-jacobian
            // is log(p1) + log(r1) + log(p2) + log(r2) + log(p3) + log(r3) - log(p1) - log(p2)
            // = log(r1) + log(r2) + log(r3) + log(p3)
            std::vector<double> tmp(_subset_relrates.begin(), _subset_relrates.end());
            // for transformation of relrate vector to dirichlet-distributed random variable
            double log_jacobian_correction = 0.0;               //POLMOD3
            assert(_num_subsets == tmp.size());                 //POLMOD3
            for (unsigned i = 0; i < _num_subsets; i++) {       //POLMOD3
                double p = 1.0*_subset_sizes[i]/_num_sites;     //POLMOD3
                tmp[i] *= p;                                    //POLMOD3
                if (i < _num_subsets - 1) {                     //POLMOD3: changed, was if (i > 0)...
                    log_jacobian_correction -= std::log(p);     //POLMOD3
                }                                               //POLMOD3
            }                                                   //POLMOD3
            // reduces length of tmp by 1                       //POLMOD3
            log_jacobian += logRatioTransform(tmp);             //POLMOD3 2022-01-05 (no correction)
            log_jacobian += log_jacobian_correction;            //POLMOD3 2022-01-30 (correction)
            param_vect.insert(param_vect.end(), tmp.begin(), tmp.end());
        }
        for (k = 0; k < _num_subsets; k++) {
            if (_subset_datatypes[k].isNucleotide()) {
                if (!_qmatrix[k]->isFixedExchangeabilities()) {
                    QMatrix::freq_xchg_t x = *_qmatrix[k]->getExchangeabilitiesSharedPtr();
                    log_jacobian += logRatioTransform(x);
                    param_vect.insert(param_vect.end(), std::begin(x), std::end(x));
                }
                
                if (!_qmatrix[k]->isFixedStateFreqs()) {
                    QMatrix::freq_xchg_t f = *_qmatrix[k]->getStateFreqsSharedPtr();
                    log_jacobian += logRatioTransform(f);
                    param_vect.insert(param_vect.end(), std::begin(f), std::end(f));
                }
            }
            else if (_subset_datatypes[k].isCodon()) {
                if (!_qmatrix[k]->isFixedOmega()) {
                    double log_omega = log(_qmatrix[k]->getOmega());
                    log_jacobian += log_omega;
                    param_vect.push_back(log_omega);
                }
                
                if (!_qmatrix[k]->isFixedStateFreqs()) {
                    QMatrix::freq_xchg_t f = *_qmatrix[k]->getStateFreqsSharedPtr();
                    log_jacobian += logRatioTransform(f);
                    param_vect.insert(param_vect.end(), std::begin(f), std::end(f));
                }
            }
            if (_asrv[k]->getIsInvarModel()) {
                if (!_asrv[k]->isFixedPinvar()) {
                    double pinvar = _asrv[k]->getPinvar();
                    assert(pinvar > 0.0 && pinvar < 1.0);
                    double log_pinvar = log(pinvar);
                    double log_one_minus_pinvar = log(1.0 - pinvar);
                    log_jacobian += log_pinvar;
                    log_jacobian += log_one_minus_pinvar;
                    param_vect.push_back(log_pinvar - log_one_minus_pinvar);
                }
            }
#if defined(HOLDER_ETAL_PRIOR)
            if (_asrv[k]->getNumCateg() > 1) {
                if (!_asrv[k]->isFixedShape()) {
                    double log_shape = log(_asrv[k]->getShape());
                    log_jacobian += log_shape;
                    param_vect.push_back(log_shape);
                }
            }
#else
            if (_asrv[k]->getNumCateg() > 1) {
                if (!_asrv[k]->isFixedRateVar()) {
                    double log_ratevar = log(_asrv[k]->getRateVar());
                    log_jacobian += log_ratevar;
                    param_vect.push_back(log_ratevar);
                }
            }
#endif
        }
        return log_jacobian;
    }

    inline double Model::setParametersFromLogTransformed(Eigen::VectorXd & param_vect, unsigned first, unsigned nparams) {
        unsigned k;
        unsigned cursor = first;
        double log_jacobian = 0.0;
        if (_num_subsets > 1) {
            assert(param_vect.rows() >= cursor + _num_subsets - 1);

            // Copy log-ratio-transformed subset relative rates to temporary vector
            std::vector<double> tmp(_num_subsets-1);
            for (unsigned i = 0; i < _num_subsets - 1; ++i) {
                tmp[i] = param_vect(cursor + i);
            }
            log_jacobian += logRatioUntransform(tmp);
            assert(tmp.size() == _num_subsets); // tmp should have increased in size by 1
            
            // Recover relative rates                         //POLMOD3
            double log_jacobian_correction = 0.0;             //POLMOD3 2022-01-05 (no correction)
            for (unsigned i = 0; i < _num_subsets; ++i) {     //POLMOD3 2022-01-30 (correction)
                double p = 1.0*_subset_sizes[i]/_num_sites;   //POLMOD3
                tmp[i] /= p;                                  //POLMOD3
                if (i < _num_subsets - 1) {                   //POLMOD3: changed, was if (i > 0)...
                    log_jacobian_correction += std::log(p);   //POLMOD3
                }                                             //POLMOD3
            }                                                 //POLMOD3
            
            // log_jacobian_correction accounts for transformation from (p1*r1, p2*r2, p3*r3) --> (r1,r2,r3)
            log_jacobian += log_jacobian_correction;
            
            // Copy detransformed subset relative rates to model
            std::copy(tmp.begin(), tmp.end(), _subset_relrates.begin());
            cursor += _num_subsets - 1;
        }
        for (k = 0; k < _num_subsets; k++) {
            if (_subset_datatypes[k].isNucleotide()) {
                if (!_qmatrix[k]->isFixedExchangeabilities()) {
                    assert(param_vect.rows() >= cursor + 5);
                    QMatrix::freq_xchg_t x(5);
                    x[0] = param_vect(cursor++);
                    x[1] = param_vect(cursor++);
                    x[2] = param_vect(cursor++);
                    x[3] = param_vect(cursor++);
                    x[4] = param_vect(cursor++);
                    log_jacobian += logRatioUntransform(x);
                    assert(x.size() == 6); // x should have increased in size by 1
                    _qmatrix[k]->setExchangeabilities(x);
                }

                if (!_qmatrix[k]->isFixedStateFreqs()) {
                    assert(param_vect.rows() >= cursor + 3);
                    QMatrix::freq_xchg_t f(3);
                    f[0] = param_vect(cursor++);
                    f[1] = param_vect(cursor++);
                    f[2] = param_vect(cursor++);
                    log_jacobian += logRatioUntransform(f);
                    assert(f.size() == 4); // f should have increased in size by 1
                    _qmatrix[k]->setStateFreqs(f);
                }
            }
            else if (_subset_datatypes[k].isCodon()) {
                if (!_qmatrix[k]->isFixedOmega()) {
                    assert(param_vect.rows() >= cursor + 1);
                    double log_omega = param_vect(cursor++);
                    log_jacobian += log_omega;
                    double omega = exp(log_omega);
                    _qmatrix[k]->setOmega(omega);
                }

                if (!_qmatrix[k]->isFixedStateFreqs()) {
                    assert(param_vect.rows() >= cursor + 60);
                    QMatrix::freq_xchg_t f(60);
                    for (unsigned i = 0; i < 60; ++i)
                        f[i] = param_vect(cursor++);
                    log_jacobian += logRatioUntransform(f);
                    assert(f.size() == 61); // f should have increased in size by 1
                    _qmatrix[k]->setStateFreqs(f);
                }
            }
            if (_asrv[k]->getIsInvarModel()) {
                if (!_asrv[k]->isFixedPinvar()) {
                    assert(param_vect.rows() >= cursor + 1);
                    double logit_pinvar = param_vect(cursor++);
                    double log_pinvar = logit_pinvar - log(1.0 + exp(logit_pinvar));
                    log_jacobian += log_pinvar;
                    double pinvar = exp(log_pinvar);
                    log_jacobian += log(1.0 - pinvar);
                    _asrv[k]->setPinvar(pinvar);
                }
            }
#if defined(HOLDER_ETAL_PRIOR)
            if (_asrv[k]->getNumCateg() > 1) {
                if (!_asrv[k]->isFixedShape()) {
                    assert(param_vect.rows() >= cursor + 1);
                    double log_shape = param_vect(cursor++);
                    log_jacobian += log_shape;
                    double shape = exp(log_shape);
                    _asrv[k]->setShape(shape);
                }
            }
#else
            if (_asrv[k]->getNumCateg() > 1) {
                if (!_asrv[k]->isFixedRateVar()) {
                    assert(param_vect.rows() >= cursor + 1);
                    double log_ratevar = param_vect(cursor++);
                    log_jacobian += log_ratevar;
                    double ratevar = exp(log_ratevar);
                    _asrv[k]->setRateVar(ratevar);
                }
            }
#endif
        }
        return log_jacobian;
    }

#if defined(HOLDER_ETAL_PRIOR)
    inline void Model::setEdgeLenRefDistParams(std::vector<double> & edgelen_refdist_params) {
        _edgelen_refdist_params = edgelen_refdist_params;
    }
    
    inline std::vector<double> Model::getEdgeLenRefDistParamsVect() const {
        return _edgelen_refdist_params;
    }
#else
    inline void Model::setEdgeProportionsRefDistParams(std::vector<double> & edgeprops_refdist_params) {
        _edgeprops_refdist_params = edgeprops_refdist_params;
    }
    
    inline std::vector<double> Model::getEdgeProportionsRefDistParamsVect() const {
        return _edgeprops_refdist_params;
    }

    inline void Model::setTreeLengthRefDistParams(std::vector<double> & treelen_refdist_params) {
        _treelen_refdist_params = treelen_refdist_params;
    }
    
    inline std::vector<double> Model::getTreeLengthRefDistParamsVect() const {
        return _treelen_refdist_params;
    }
#endif

    inline void Model::setSampledSubsetRelRates(unsigned i) {
        assert(_num_subsets > 0);
        assert(_sampled_subset_relrates.size() > i);
        assert(_sampled_subset_relrates[i].size() > 0);
        _subset_relrates.assign(_sampled_subset_relrates[i].begin(), _sampled_subset_relrates[i].end());
    }

    inline void Model::setSampledExchangeabilities(unsigned subset, unsigned i) {
        assert(_num_subsets > subset);
        assert(_sampled_exchangeabilities[subset].size() > i);
        assert(_sampled_exchangeabilities[subset][i].size() > 0);
        _qmatrix[subset]->setExchangeabilities(_sampled_exchangeabilities[subset][i]);
    }
    
    inline void Model::setSampledStateFreqs(unsigned subset, unsigned i) {
        assert(_num_subsets > subset);
        assert(_sampled_state_freqs[subset].size() > i);
        assert(_sampled_state_freqs[subset][i].size() > 0);
        _qmatrix[subset]->setStateFreqs(_sampled_state_freqs[subset][i]);
    }
    
    inline void Model::setSampledOmega(unsigned subset, unsigned i) {
        assert(_num_subsets > subset);
        assert(_sampled_omegas[subset].size() > i);
        _qmatrix[subset]->setOmega(_sampled_omegas[subset][i]);
    }
    
    inline void Model::setSampledPinvar(unsigned subset, unsigned i) {
        assert(_num_subsets > subset);
        assert(_sampled_pinvars[subset].size() > i);
        _asrv[subset]->setPinvar(_sampled_pinvars[subset][i]);
    }

#if defined(HOLDER_ETAL_PRIOR)
    inline void Model::setSampledShape(unsigned subset, unsigned i) {
        assert(_num_subsets > subset);
        assert(_sampled_shapes[subset].size() > i);
        assert(_sampled_shapes[subset][i].size() > 0);
        _asrv[subset]->setShape(_sampled_shapes[subset][i]);
    }
#else
    inline void Model::setSampledRateVar(unsigned subset, unsigned i) {
        assert(_num_subsets > subset);
        assert(_sampled_ratevars[subset].size() > i);
        _asrv[subset]->setRateVar(_sampled_ratevars[subset][i]);
    }
#endif

    inline void Model::setModelToSampledPoint(unsigned i) {
        if (_num_subsets > 1) {
            setSampledSubsetRelRates(i);
        }
        for (unsigned k = 0; k < _num_subsets; k++) {
            if (_subset_datatypes[k].isNucleotide()) {
                setSampledExchangeabilities(k, i);
                setSampledStateFreqs(k, i);
            }
            else if (_subset_datatypes[k].isCodon()) {
                setSampledOmega(k, i);
                setSampledStateFreqs(k, i);
            }
            if (_asrv[k]->getIsInvarModel()) {
                setSampledPinvar(k, i);
            }
#if defined(HOLDER_ETAL_PRIOR)
            if (_asrv[k]->getNumCateg() > 1) {
                setSampledShape(k, i);
            }
#else
            if (_asrv[k]->getNumCateg() > 1) {
                setSampledRateVar(k, i);
            }
#endif
        }
    }

    inline void Model::sampleParams() {
        unsigned k;
        if (_num_subsets > 1) {
            _sampled_subset_relrates.push_back(_subset_relrates);
        }
        for (k = 0; k < _num_subsets; k++) {
            if (_subset_datatypes[k].isNucleotide()) {
                QMatrix::freq_xchg_t & x = *_qmatrix[k]->getExchangeabilitiesSharedPtr();
                _sampled_exchangeabilities[k].push_back(x);
                
                QMatrix::freq_xchg_t & f = *_qmatrix[k]->getStateFreqsSharedPtr();
                _sampled_state_freqs[k].push_back(f);
            }
            else if (_subset_datatypes[k].isCodon()) {
                _sampled_omegas[k].push_back(_qmatrix[k]->getOmega());

                QMatrix::freq_xchg_t & f = *_qmatrix[k]->getStateFreqsSharedPtr();
                _sampled_state_freqs[k].push_back(f);
            }
            if (_asrv[k]->getIsInvarModel()) {
                _sampled_pinvars[k].push_back(_asrv[k]->getPinvar());
            }
#if defined(HOLDER_ETAL_PRIOR)
            if (_asrv[k]->getNumCateg() > 1) {
                _sampled_shapes[k].push_back(_asrv[k]->getShape());
            }
#else
            if (_asrv[k]->getNumCateg() > 1) {
                _sampled_ratevars[k].push_back(_asrv[k]->getRateVar());
            }
#endif
        }
    }

    inline std::string Model::calcGammaRefDist(std::string title, std::string subset_name, std::vector<double> & vect, std::vector<double> & params) {
        //TODO: nearly identical to TreeManip::calcGammaRefDist - make one version that can be used by both Model and TreeManip
        // Compute sums and sums-of-squares for each component
        unsigned n = (unsigned)vect.size();
        double sumv   = 0.0;
        double sumsqv = 0.0;
        for (unsigned i = 0; i < n; i++) {
            double v = vect[i];
            sumv += v;
            sumsqv += v*v;
        }
        
        // Compute mean and variance
        double mu;
        double s;
        mu = sumv/n;
        s = (sumsqv - mu*mu*n)/(n-1);
        
        // Compute parameters of reference distribution and save each
        // as an element of the string vector svect
        // s = shape*scale^2
        // mu = shape*scale
        double scale = s/mu;
        double shape = mu/scale;
        params.resize(2);
        params[0] = shape;
        params[1] = scale;
        std::string refdiststr = boost::str(boost::format("%s = %s:%.3f, %.3f\n") % title % subset_name % shape % scale);
        
        return refdiststr;
    }
    
    inline std::string Model::calcBetaRefDist(std::string title, std::string subset_name, std::vector<double> & vect, std::vector<double> & params) {
        // mean = a/(a+b)  1-mean = b/(a+b)  var = ab/[(a+b)^2 (a+b+1)]
        // phi = a+b = mean*(1-mean)/var   a = mean*phi   b = (1-mean)*phi
        
        // Calculate mean and var
        double mean = 0.0;
        double var = 0.0;
        for (auto p : vect) {
            mean += p;
            var += p*p;
        }
        double n = (double)vect.size();
        assert(n > 1.0);
        mean /= n;
        var = (var - n*mean*mean)/(n-1);
        
        // Calculate a and b
        double phi = mean*(1.0 - mean)/var;
        double a = mean*phi;
        double b = (1.0-mean)*phi;
        params.resize(2);
        params[0] = a;
        params[1] = b;
        std::string refdiststr = boost::str(boost::format("%s = %s:%.3f, %.3f\n") % title % subset_name % a % b);
        return refdiststr;
    }
    
    inline std::string Model::calcDirichletRefDist(std::string title, std::string subset_name, std::vector< QMatrix::freq_xchg_t > & vect, std::vector<double> & params, bool relrates) {
        // Sanity check: also calculate means and variances using Boost accumulator
        // see https://www.nu42.com/2016/12/descriptive-stats-with-cpp-boost.html
        //boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::variance> > acc;
        
        //TODO: nearly identical to TreeManip::calcGammaRefDist - make one version that can be used by both Model and TreeManip
        // Ming-Hui Chen method of matching component variances
        // mu_i = phi_i/phi is mean of component i (estimate using sample mean)
        // s_i^2 is sample variance of component i
        //
        //       sum_i mu_i^2 (1 - mu_i)^2
        // phi = --------------------------- - 1
        //       sum_i s_i^2 mu_i (1 - mu_i)
        //
        // phi_i = phi mu_i
        unsigned n = (unsigned)vect.size();
        assert(n > 0);
        unsigned k = (unsigned)vect[0].size();
        
#if defined(POLTMPPARAMS)
        std::ofstream outf(boost::str(boost::format("doof-%s.txt") % title));
        if (title == "exchangerefdist") {
            outf << "rAC-0	rAG-0	rAT-0	rCG-0	rCT-0	rGT-0\n";
        }
        else if (title == "statefreqrefdist") {
            outf << "piA-0	piC-0	piG-0	piT-0\n";
        }
#endif
        // Compute sums and sums-of-squares for each component
        std::vector<double> sums(k, 0.0);
        std::vector<double> sumsq(k, 0.0);
        for (unsigned i = 0; i < n; i++) {
            QMatrix::freq_xchg_t & dir = vect[i];
            //acc(dir[0]);
            for (unsigned j = 0; j < k; j++) {
                double v = dir[j];
#if defined(POLTMPPARAMS)
                outf << boost::format("%.5f\t") % v;
#endif
                sums[j] += v;
                sumsq[j] += v*v;
            }
#if defined(POLTMPPARAMS)
            outf << "\n";
#endif
        }
#if defined(POLTMPPARAMS)
        outf.close();
#endif
        
        // Compute means and variances for each component
        std::vector<double> mu(k, 0.0);
        std::vector<double> s2(k, 0.0);
        double numer = 0.0;
        double denom = 0.0;
        for (unsigned j = 0; j < k; j++) {
            mu[j] = sums[j]/n;
            s2[j] = (sumsq[j] - mu[j]*mu[j]*n)/(n-1);
            if (relrates) {
                double pj_inv = 1.0*_num_sites/_subset_sizes[j];
                numer += mu[j]*mu[j]*(pj_inv - mu[j])*(pj_inv - mu[j]);
                denom += s2[j]*mu[j]*(pj_inv - mu[j]);
            }
            else {
                numer += mu[j]*mu[j]*(1.0 - mu[j])*(1.0 - mu[j]);
                denom += s2[j]*mu[j]*(1.0 - mu[j]);
            }
        }
        
        // Compute phi
        double phi = numer/denom - 1.0;
#if defined(POLTMPPARAMS)
        ::om.outputConsole(boost::format("%s:\n") % title;
        ::om.outputConsole("  phi = %g\n") % phi);
        for (unsigned j = 0; j < k; j++) {
            ::om.outputConsole(boost::format("  mu[%d] = %g\n") % j % mu[j]);
        }
        ::om.outputConsole();
#endif

        // Compute parameters of reference distribution and save each
        // as an element of the string vector svect
        params.clear();
        std::vector<std::string> svect;
        for (unsigned j = 0; j < k; j++) {
            double c = phi*mu[j];
            if (relrates) {
                double pj = 1.0*_subset_sizes[j]/_num_sites;
                c *= pj;
            }
            params.push_back(c);
            std::string stmp = boost::str(boost::format("%.3f") % c);
            svect.push_back(stmp);
        }
        std::string refdiststr = boost::str(boost::format("%s = %s:%s\n") % title % subset_name % boost::algorithm::join(svect, ","));
        
        return refdiststr;
    }
    
    inline std::string Model::saveReferenceDistributions(Partition::SharedPtr partition) {
        std::map<std::string, std::vector<double> > dummy_refdist_map;
        std::string s = calcReferenceDistributions(partition, dummy_refdist_map);
        return s;
    }
    
    inline std::string Model::calcReferenceDistributions(Partition::SharedPtr partition, std::map<std::string, std::vector<double> > & refdist_map) {

        // Calculate and save reference distribution parameters in a conf file that can be used
        // in a subsequent generalized steppingstone analysis
        unsigned k;
        std::string s;
        if (_num_subsets > 1) {
            std::vector<double> & v = refdist_map["Subset Relative Rates"];
            s += calcDirichletRefDist("relratesrefdist", "default", _sampled_subset_relrates, v, true /* relative rate distribution*/);
        }
        for (k = 0; k < _num_subsets; k++) {
            if (_subset_datatypes[k].isNucleotide()) {
                if (!_qmatrix[k]->isFixedExchangeabilities()) {
                    std::vector<double> & v = refdist_map["Exchangeabilities"];
                    s += calcDirichletRefDist("exchangerefdist", partition->getSubsetName(k), _sampled_exchangeabilities[k], v);
                }
                if (!_qmatrix[k]->isFixedStateFreqs()) {
                    std::vector<double> & v = refdist_map["State Frequencies"];
                    s += calcDirichletRefDist("statefreqrefdist", partition->getSubsetName(k), _sampled_state_freqs[k], v);
                }
            }
            else if (_subset_datatypes[k].isCodon()) {
                if (!_qmatrix[k]->isFixedOmega()) {
                    std::vector<double> & v = refdist_map["Omega"];
                    s += calcGammaRefDist("omegarefdist", partition->getSubsetName(k), _sampled_omegas[k], v);
                }
                if (!_qmatrix[k]->isFixedStateFreqs()) {
                    std::vector<double> & v = refdist_map["State Frequencies"];
                    s += calcDirichletRefDist("statefreqrefdist", partition->getSubsetName(k), _sampled_state_freqs[k], v);
                }
            }
            if (_asrv[k]->getIsInvarModel()) {
                if (!_asrv[k]->isFixedPinvar()) {
                    std::vector<double> & v = refdist_map["Proportion of Invariable Sites"];
                    s += calcBetaRefDist("pinvarrefdist", partition->getSubsetName(k), _sampled_pinvars[k], v);
                }
            }
#if defined(HOLDER_ETAL_PRIOR)
            if (_asrv[k]->getNumCateg() > 1) {
                if (!_asrv[k]->isFixedShape()) {
                    std::vector<double> & v = refdist_map["Gamma Shape"];
                    s += calcGammaRefDist("shaperefdist", partition->getSubsetName(k), _sampled_shapes[k], v);
                }
            }
#else
            if (_asrv[k]->getNumCateg() > 1) {
                if (!_asrv[k]->isFixedRateVar()) {
                    std::vector<double> & v = refdist_map["Gamma Rate Variance"];
                    s += calcGammaRefDist("ratevarrefdist", partition->getSubsetName(k), _sampled_ratevars[k], v);
                }
            }
#endif
        }

        return s;
    }
        
    inline void Model::setSubsetRelRatesRefDistParams(std::vector<double> refdist_params) {
        _subset_relrates_refdist_params.resize(refdist_params.size());
        std::copy(refdist_params.begin(), refdist_params.end(), _subset_relrates_refdist_params.begin());
    }

    inline std::vector<double> Model::getSubsetRelRatesRefDistParamsVect() {
        return _subset_relrates_refdist_params;
    }
}
