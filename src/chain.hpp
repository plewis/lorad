#pragma once

//#define POLTMPPRIOR

#include "conditionals.hpp"

#include <memory>
#include <boost/format.hpp>
#include "lot.hpp"
#include "data.hpp"
#include "tree.hpp"
#include "model.hpp"
#include "likelihood.hpp"
#include "tree_manip.hpp"
#include "updater.hpp"
#include "gamma_ratevar_updater.hpp"
#include "omega_updater.hpp"
#include "pinvar_updater.hpp"
#include "statefreq_updater.hpp"
#include "exchangeability_updater.hpp"
#include "subset_relrate_updater.hpp"
#include "tree_updater.hpp"
#include "polytomy_updater.hpp"
#include "tree_length_updater.hpp"
#include "edge_proportion_updater.hpp"

namespace strom {

    class Chain {
    
        friend class Likelihood;

        public:
        
            typedef std::vector<Updater::SharedPtr> updater_vect_t;
            typedef std::shared_ptr<Chain>          SharedPtr;

                                                    Chain();
                                                    ~Chain();

            void                                    clear();

            void                                    startTuning();
            void                                    stopTuning();

            void                                    setTreeFromNewick(std::string & newick);
            unsigned                                createUpdaters(Model::SharedPtr model, Lot::SharedPtr lot, Likelihood::SharedPtr likelihood);

            TreeManip::SharedPtr                    getTreeManip();
            Model::SharedPtr                        getModel();
            double                                  getLogLikelihood() const;


            void                                    setHeatingPower(double p);
            double                                  getHeatingPower() const;

            void                                    setNextHeatingPower(double p);
            void                                    storeLogLikelihood();
            double                                  calcLogSteppingstoneRatio() const;

            void                                    setChainIndex(unsigned idx);
            double                                  getChainIndex() const;

            Updater::SharedPtr                      findUpdaterByName(std::string name);
            std::vector<std::string>                getUpdaterNames() const;
            std::vector<double>                     getAcceptPercentages() const;
            std::vector<unsigned>                   getNumUpdates() const;
            std::vector<double>                     getLambdas() const;
            void                                    setLambdas(std::vector<double> & v);

            double                                  calcLogLikelihood() const;
#if defined(ALWAYS_UPDATE_EDGE_PROPORTIONS)
            double                                  calcLogJointPrior(bool verbose = false) const;
#else
            double                                  calcLogJointPrior() const;
#endif
#if defined(POLGSS)
            double                                  calcLogReferenceDensity() const;
            void                                    setSteppingstoneMode(unsigned mode);
#endif

            void                                    start();
            void                                    stop();
            void                                    nextStep(int iteration);

        private:

            Model::SharedPtr                        _model;
            Lot::SharedPtr                          _lot;
            TreeManip::SharedPtr                    _tree_manipulator;

            updater_vect_t                          _updaters;
            
#if defined(ALWAYS_UPDATE_EDGE_PROPORTIONS)
            updater_vect_t                          _prior_calculators;
#endif
            unsigned                                _chain_index;
            double                                  _heating_power;

            bool                                    _heat_likelihood_only;
            double                                  _next_heating_power;
            std::vector<double>                     _ss_loglikes;

#if defined(POLGSS)
            std::vector<double>                     _ss_logpriors;
            std::vector<double>                     _ss_logrefdists;
            unsigned                                _ss_mode;
#endif
            double                                  _log_likelihood;
    };
    
    inline Chain::Chain() {
        //std::cout << "Chain being created" << std::endl;
        clear();
    }

    inline Chain::~Chain() {
        //std::cout << "Chain being destroyed" << std::endl;
    } 

    inline void Chain::clear() {
        _log_likelihood = 0.0;
        _updaters.clear();
        _chain_index = 0;
        setHeatingPower(1.0);
        _heat_likelihood_only = false;
        _next_heating_power = 1.0;
        _ss_loglikes.clear();
#if defined(POLGSS)
        _ss_logpriors.clear();
        _ss_logrefdists.clear();
        _ss_mode = 0;
#endif
        startTuning();
    }

    inline void Chain::startTuning() {
        for (auto u : _updaters)
            u->setTuning(true);
    }

    inline void Chain::stopTuning() {
        for (auto u : _updaters)
            u->setTuning(false);
    }

    inline void Chain::setTreeFromNewick(std::string & newick) {
        assert(_updaters.size() > 0);
        if (!_tree_manipulator)
            _tree_manipulator.reset(new TreeManip);
        _tree_manipulator->buildFromNewick(newick, false, true);
        for (auto u : _updaters)
            u->setTreeManip(_tree_manipulator);
    }

    inline unsigned Chain::createUpdaters(Model::SharedPtr model, Lot::SharedPtr lot, Likelihood::SharedPtr likelihood) {
        _model = model;
        _lot = lot;
        _updaters.clear();
#if defined(ALWAYS_UPDATE_EDGE_PROPORTIONS)
        _prior_calculators.clear();
#endif

        double wstd             = 1.0;
        double wtreelength      = 1.0;
#if defined(ALWAYS_UPDATE_EDGE_PROPORTIONS)
        double wtreetopology    = 10.0;
        double wedgelengths     = 10.0;
#else
        double wtreetopology    = 19.0;
        double wedgelengths     = 1.0;
#endif
        double wpolytomy        = 0.0;
        double sum_weights      = 0.0;
        
        if (_model->isAllowPolytomies()) {
            wstd             = 1.0;
            wtreelength      = 2.0;
            wtreetopology    = 9.0;
            wpolytomy        = 9.0;
        }
        
        // Add state frequency parameter updaters to _updaters
        Model::state_freq_params_t & statefreq_shptr_vect = _model->getStateFreqParams();
        for (auto statefreq_shptr : statefreq_shptr_vect) {
            Updater::SharedPtr u = StateFreqUpdater::SharedPtr(new StateFreqUpdater(statefreq_shptr));
            u->setLikelihood(likelihood);
            u->setLot(lot);
            u->setLambda(0.0001);
            u->setTargetAcceptanceRate(0.3);
            u->setPriorParameters(std::vector<double>(statefreq_shptr->getStateFreqsSharedPtr()->size(), 1.0));
#if defined(POLGSS)
            u->setRefDistParameters(statefreq_shptr->getStateFreqRefDistParamsVect());
#endif
            u->setWeight(wstd); sum_weights += wstd;
            _updaters.push_back(u);
#if defined(ALWAYS_UPDATE_EDGE_PROPORTIONS)
            _prior_calculators.push_back(u);
#endif
        }

        // Add exchangeability parameter updaters to _updaters 
        Model::exchangeability_params_t & exchangeability_shptr_vect = _model->getExchangeabilityParams();
        for (auto exchangeability_shptr : exchangeability_shptr_vect) {
            Updater::SharedPtr u = ExchangeabilityUpdater::SharedPtr(new ExchangeabilityUpdater(exchangeability_shptr));
            u->setLikelihood(likelihood);
            u->setLot(lot);
            u->setLambda(0.0001);
            u->setTargetAcceptanceRate(0.3);
            u->setPriorParameters({1.0, 1.0, 1.0, 1.0, 1.0, 1.0});
#if defined(POLGSS)
            u->setRefDistParameters(exchangeability_shptr->getExchangeabilityRefDistParamsVect());
#endif
            u->setWeight(wstd); sum_weights += wstd;
            _updaters.push_back(u);
#if defined(ALWAYS_UPDATE_EDGE_PROPORTIONS)
            _prior_calculators.push_back(u);
#endif
        }

        // Add rate variance parameter updaters to _updaters
        Model::ratevar_params_t & ratevar_shptr_vect = _model->getRateVarParams();
        for (auto ratevar_shptr : ratevar_shptr_vect) {
            Updater::SharedPtr u = GammaRateVarUpdater::SharedPtr(new GammaRateVarUpdater(ratevar_shptr));
            u->setLikelihood(likelihood);
            u->setLot(lot);
            u->setLambda(1.0);
            u->setTargetAcceptanceRate(0.3);
            u->setPriorParameters({1.0, 1.0});
#if defined(POLGSS)
            u->setRefDistParameters(ratevar_shptr->getRateVarRefDistParamsVect());
#endif
            u->setWeight(wstd); sum_weights += wstd;
            _updaters.push_back(u);
#if defined(ALWAYS_UPDATE_EDGE_PROPORTIONS)
            _prior_calculators.push_back(u);
#endif
        }
        
        // Add pinvar parameter updaters to _updaters
        Model::pinvar_params_t & pinvar_shptr_vect = _model->getPinvarParams();
        for (auto pinvar_shptr : pinvar_shptr_vect) {
            Updater::SharedPtr u = PinvarUpdater::SharedPtr(new PinvarUpdater(pinvar_shptr));
            u->setLikelihood(likelihood);
            u->setLot(lot);
            u->setLambda(0.5);
            u->setTargetAcceptanceRate(0.3);
            u->setPriorParameters({1.0, 1.0});
//#if defined(POLGSS)
//            u->setRefDistParameters(pinvar_shptr->getPinvarRefDistParamsVect());
//#endif
            u->setWeight(wstd); sum_weights += wstd;
            _updaters.push_back(u);
#if defined(ALWAYS_UPDATE_EDGE_PROPORTIONS)
            _prior_calculators.push_back(u);
#endif
        }
        
        // Add omega parameter updaters to _updaters
        Model::omega_params_t & omega_shptr_vect = _model->getOmegaParams();
        for (auto omega_shptr : omega_shptr_vect) {
            Updater::SharedPtr u = OmegaUpdater::SharedPtr(new OmegaUpdater(omega_shptr));
            u->setLikelihood(likelihood);
            u->setLot(lot);
            u->setLambda(1.0);
            u->setTargetAcceptanceRate(0.3);
            u->setPriorParameters({1.0, 1.0});
//#if defined(POLGSS)
//            u->setRefDistParameters(omega_shptr->getOmegaRefDistParamsVect());
//#endif
            u->setWeight(wstd); sum_weights += wstd;
            _updaters.push_back(u);
#if defined(ALWAYS_UPDATE_EDGE_PROPORTIONS)
            _prior_calculators.push_back(u);
#endif
        }
        
        // Add subset relative rate parameter updater to _updaters
        if (!_model->isFixedSubsetRelRates()) {
            Updater::SharedPtr u = SubsetRelRateUpdater::SharedPtr(new SubsetRelRateUpdater(_model));
            u->setLikelihood(likelihood);
            u->setLot(lot);
            u->setLambda(1.0);
            u->setTargetAcceptanceRate(0.3);
            u->setPriorParameters(std::vector<double>(_model->getNumSubsets(), 1.0));
//#if defined(POLGSS)
//            u->setRefDistParameters(_model->getSubsetRelRatesRefDistParamsVect());
//#endif
            u->setWeight(wstd); sum_weights += wstd;
            _updaters.push_back(u);
#if defined(ALWAYS_UPDATE_EDGE_PROPORTIONS)
            _prior_calculators.push_back(u);
#endif
        }
        
        // Add tree updater and tree length updater to _updaters
        double tree_length_shape = 1.0;
        double tree_length_scale = 10.0;
        double dirichlet_param   = 1.0;
                    
#if defined(ALWAYS_UPDATE_EDGE_PROPORTIONS)
        Updater::SharedPtr uu = EdgeProportionUpdater::SharedPtr(new EdgeProportionUpdater());
        uu->setLikelihood(likelihood);
        uu->setLot(lot);
        uu->setLambda(0.0001);
        uu->setTargetAcceptanceRate(0.3);
        uu->setPriorParameters({tree_length_shape, tree_length_scale, dirichlet_param});
#if defined(POLGSS)
        uu->setRefDistParameters(_model->getEdgeProportionsRefDistParamsVect());
#endif
        uu->setWeight(wedgelengths); sum_weights += wedgelengths;
        _updaters.push_back(uu);
        
        if (!_model->isFixedTree()) {
            Updater::SharedPtr u = TreeUpdater::SharedPtr(new TreeUpdater());
            u->setLikelihood(likelihood);
            u->setLot(lot);
            u->setLambda(0.5);
            u->setTargetAcceptanceRate(0.3);
            u->setPriorParameters({tree_length_shape, tree_length_scale, dirichlet_param});
            u->setTopologyPriorOptions(_model->isResolutionClassTopologyPrior(), _model->getTopologyPriorC());
            u->setWeight(wtreetopology); sum_weights += wtreetopology;
            _updaters.push_back(u);

            if (_model->isAllowPolytomies()) {
                Updater::SharedPtr u = PolytomyUpdater::SharedPtr(new PolytomyUpdater());
                u->setLikelihood(likelihood);
                u->setLot(lot);
                u->setLambda(0.05);
                u->setTargetAcceptanceRate(0.5);
                u->setPriorParameters({tree_length_shape, tree_length_scale, dirichlet_param});
                u->setTopologyPriorOptions(_model->isResolutionClassTopologyPrior(), _model->getTopologyPriorC());
                u->setWeight(wpolytomy); sum_weights += wpolytomy;
                _updaters.push_back(u);
            }
        }
#else
        if (_model->isFixedTree()) {
            Updater::SharedPtr u = EdgeProportionUpdater::SharedPtr(new EdgeProportionUpdater());
            u->setLikelihood(likelihood);
            u->setLot(lot);
            u->setLambda(0.2);
            u->setTargetAcceptanceRate(0.3);
            u->setPriorParameters({tree_length_shape, tree_length_scale, dirichlet_param});
#if defined(POLGSS)
            u->setRefDistParameters(_model->getEdgeProportionsRefDistParamsVect());
#endif
            u->setWeight(wedgelengths); sum_weights += wedgelengths;
            _updaters.push_back(u);
        }
        else {
            Updater::SharedPtr u = TreeUpdater::SharedPtr(new TreeUpdater());
            u->setLikelihood(likelihood);
            u->setLot(lot);
            u->setLambda(0.5);
            u->setTargetAcceptanceRate(0.3);
            u->setPriorParameters({tree_length_shape, tree_length_scale, dirichlet_param});
            u->setTopologyPriorOptions(_model->isResolutionClassTopologyPrior(), _model->getTopologyPriorC());
            u->setWeight(wtreetopology); sum_weights += wtreetopology;
            _updaters.push_back(u);

            if (_model->isAllowPolytomies()) {
                Updater::SharedPtr u = PolytomyUpdater::SharedPtr(new PolytomyUpdater());
                u->setLikelihood(likelihood);
                u->setLot(lot);
                u->setLambda(0.05);
                u->setTargetAcceptanceRate(0.5);
                u->setPriorParameters({tree_length_shape, tree_length_scale, dirichlet_param});
                u->setTopologyPriorOptions(_model->isResolutionClassTopologyPrior(), _model->getTopologyPriorC());
                u->setWeight(wpolytomy); sum_weights += wpolytomy;
                _updaters.push_back(u);
            }
        }
#endif
        
        Updater::SharedPtr u = TreeLengthUpdater::SharedPtr(new TreeLengthUpdater());
        u->setLikelihood(likelihood);
        u->setLot(lot);
        u->setLambda(0.2);
        u->setTargetAcceptanceRate(0.3);
        u->setPriorParameters({tree_length_shape, tree_length_scale, dirichlet_param});
        u->setTopologyPriorOptions(_model->isResolutionClassTopologyPrior(), _model->getTopologyPriorC());
#if defined(POLGSS)
        u->setRefDistParameters(_model->getTreeLengthRefDistParamsVect());
#endif
        u->setWeight(wtreelength); sum_weights += wtreelength;
        _updaters.push_back(u);
#if defined(ALWAYS_UPDATE_EDGE_PROPORTIONS)
        _prior_calculators.push_back(u);
#endif
        
        for (auto u : _updaters) {
            u->calcProb(sum_weights);
        }
        
        return (unsigned)_updaters.size();
    }

    inline TreeManip::SharedPtr Chain::getTreeManip() {
        return _tree_manipulator;
    }

    inline Model::SharedPtr Chain::getModel() {
        return _model;
    }

    inline double Chain::getHeatingPower() const {
        return _heating_power;
    }

    inline void Chain::setHeatingPower(double p) {
        _heating_power = p;
        for (auto u : _updaters)
            u->setHeatingPower(p);
    }

    inline void Chain::setNextHeatingPower(double p) {
        _heat_likelihood_only = true; // next heating power only set if doing steppingstone
#if defined(POLGSS)
        // Steppingstone mode:
        //   0: no steppingstone
        //   1: steppingstone (Xie et al. 2011)
        //   2: generalized steppingstone (Fan et al. 2011)
        for (auto u : _updaters)
            u->setSteppingstoneMode(_ss_mode);
#else
        for (auto u : _updaters)
            u->setHeatLikelihoodOnly(true);
#endif
        _next_heating_power = p;
    }
    
    inline void Chain::storeLogLikelihood() {
        double logLike = getLogLikelihood();
        _ss_loglikes.push_back(logLike);
#if defined(POLGSS)
        if (_ss_mode == 2) {
            double logPrior = calcLogJointPrior();
            _ss_logpriors.push_back(logPrior);
            double logRefDist = calcLogReferenceDensity();
            _ss_logrefdists.push_back(logRefDist);
        }
#endif
    }
    
    inline double Chain::calcLogSteppingstoneRatio() const {
        // Find the maximum log likelihood sampled by this chain
        unsigned sample_size = (unsigned)_ss_loglikes.size();
        assert(sample_size > 0);

        double log_ratio = 0.0;

#if defined(POLGSS)
        if (_ss_mode == 2) {
            assert(_ss_logpriors.size() == sample_size);
            assert(_ss_logrefdists.size() == sample_size);
            double beta_diff = _next_heating_power - _heating_power;
            std::vector<double> log_ratios(sample_size, 0.0);
            double max_log_ratio = 0.0;
            for (unsigned i = 0; i < sample_size; ++i) {
                double log_ratio = _ss_loglikes[i] + _ss_logpriors[i] - _ss_logrefdists[i];
                if (i == 0 || log_ratio > max_log_ratio)
                    max_log_ratio = log_ratio;
                log_ratios[i] = log_ratio;
            }
            
            // Compute sum, factoring out maxLnL
            double sum_of_terms = 0.0;
            for (auto logR : log_ratios) {
                sum_of_terms += exp(beta_diff*(logR - max_log_ratio));
            }
            
            // Compute the log of the steppingstone ratio
            assert(sum_of_terms > 0.0);
            log_ratio = beta_diff*max_log_ratio + log(sum_of_terms) - log(sample_size);
            }
        else {
#endif
            double maxLogL = *(std::max_element(_ss_loglikes.begin(), _ss_loglikes.end()));
            
            // Compute sum, factoring out maxLnL
            double sum_of_terms = 0.0;
            for (auto logL : _ss_loglikes) {
                sum_of_terms += exp((_next_heating_power - _heating_power)*(logL - maxLogL));
            }
            
            // Compute the log of the steppingstone ratio
            assert(sum_of_terms > 0.0);
            log_ratio = (_next_heating_power - _heating_power)*maxLogL + log(sum_of_terms) - log(sample_size);
#if defined(POLGSS)
        }
#endif
        return log_ratio;
    }

    inline double Chain::getChainIndex() const {
        return _chain_index;
    }

    inline void Chain::setChainIndex(unsigned idx) {
        _chain_index = idx;
    }
        
    inline Updater::SharedPtr Chain::findUpdaterByName(std::string name) {
        Updater::SharedPtr retval = nullptr;
        for (auto u : _updaters) {
            if (u->getUpdaterName() == name) {
                retval = u;
                break;
            }
        }
        assert(retval != nullptr);
        return retval;
    } 

    inline std::vector<std::string> Chain::getUpdaterNames() const {
        std::vector<std::string> v;
        for (auto u : _updaters)
            v.push_back(u->getUpdaterName());
        return v;
    }

    inline std::vector<double> Chain::getAcceptPercentages() const {
        std::vector<double> v;
        for (auto u : _updaters)
            v.push_back(u->getAcceptPct());
        return v;
    }

    inline std::vector<unsigned> Chain::getNumUpdates() const {
        std::vector<unsigned> v;
        for (auto u : _updaters)
            v.push_back(u->getNumUpdates());
        return v;
    }

    inline std::vector<double> Chain::getLambdas() const {
        std::vector<double> v;
        for (auto u : _updaters)
            v.push_back(u->getLambda());
        return v;
    }

    inline void Chain::setLambdas(std::vector<double> & v) {
        assert(v.size() == _updaters.size());
        unsigned index = 0;
        for (auto u : _updaters) {
            u->setLambda(v[index++]);
        }
    }
    
    inline double Chain::calcLogLikelihood() const {
        return _updaters[0]->calcLogLikelihood();
    }

#if defined(ALWAYS_UPDATE_EDGE_PROPORTIONS)
    inline double Chain::calcLogJointPrior(bool verbose) const {
        double lnP = 0.0;
        for (auto u : _prior_calculators) {
            std::string this_name = u->getUpdaterName();
            if (this_name == "Tree Length") {
                auto edgelen_prior = u->calcLogEdgeLengthPrior();
                lnP += edgelen_prior.first;
                lnP += edgelen_prior.second;
                if (verbose) {
                    std::cerr << boost::format("%12.5f <-- Tree Length\n") % edgelen_prior.first;
                    std::cerr << boost::format("%12.5f <-- Edge Length Proportions\n") % edgelen_prior.second;
                }
                if (!_model->isFixedTree()) {
                    double topology_prior = u->calcLogTopologyPrior();
                    if (verbose)
                        std::cerr << boost::format("%12.5f <-- Tree Topology\n") % topology_prior;
                    lnP += topology_prior;
                }
            }
            else {
                double this_log_prior = u->calcLogPrior();
                if (verbose)
                    std::cerr << boost::format("%12.5f <-- %s\n") % this_log_prior % this_name;
                lnP += this_log_prior;
            }
        }
        if (verbose)
            std::cerr << boost::format("%12.5f <-- Log Joint Prior\n") % lnP;
        return lnP;
    }
#else
    inline double Chain::calcLogJointPrior() const {
        double lnP = 0.0;
#if defined(POLTMPPRIOR)
        std::cerr << "\nChain::calcLogJointPrior():\n";
        for (auto u : _updaters) {
            if (u->_name != "Tree Length" && u->_name != "Polytomies" ) {
                std::string this_name = u->getUpdaterName();
                double this_log_prior = u->calcLogPrior();
                std::cerr << boost::format("%12.5f <-- %s\n") % this_log_prior % this_name;
                lnP += this_log_prior;
            }
        }
        std::cerr << boost::format("%12.5f <-- joint log prior\n") % lnP;
#else
        for (auto u : _updaters) {
            if (u->_name != "Tree Length" && u->_name != "Polytomies" )
                lnP += u->calcLogPrior();
        }
#endif
        return lnP;
    }
#endif

#if defined(POLGSS)
    inline double Chain::calcLogReferenceDensity() const {
        double lnP = 0.0;
#if defined(POLTMPPRIOR)
        std::cerr << "\nChain::calcLogReferenceDensity():\n";
        for (auto u : _updaters) {
            assert(u->_name != "Polytomies");
            double log_reference_density = u->calcLogRefDist();
            std::cerr << boost::format("%12.5f <-- %s\n") % log_reference_density % u->getUpdaterName();
            lnP += log_reference_density;
        }
        std::cerr << boost::format("%12.5f <-- joint log reference density\n") % lnP;
#else
        for (auto u : _updaters) {
            assert(u->_name != "Polytomies");
            double log_reference_density = u->calcLogRefDist();
            lnP += log_reference_density;
        }
#endif
        return lnP;
    }
    
    inline void Chain::setSteppingstoneMode(unsigned mode) {
        _ss_mode = mode;
    }
#endif

    inline void Chain::start() {
        _tree_manipulator->selectAllPartials();
        _tree_manipulator->selectAllTMatrices();
        _log_likelihood = calcLogLikelihood();
    }

    inline void Chain::stop() { 
    } 

    inline void Chain::nextStep(int iteration) {
        assert(_lot);
        double u = _lot->uniform();
        double cumprob = 0.0;
        unsigned i = 0;
        for (auto updater : _updaters) {
            cumprob += updater->_prob;
            if (u <= cumprob)
                break;
            i++;
        }
        assert(i < _updaters.size());
        _log_likelihood = _updaters[i]->update(_log_likelihood);
    } 

    inline double Chain::getLogLikelihood() const {
        return _log_likelihood;
    }

}
