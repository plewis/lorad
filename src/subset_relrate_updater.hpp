#pragma once

#include "dirichlet_updater.hpp"

namespace lorad {
    
    class SubsetRelRateUpdater : public DirichletUpdater {
        public:
            typedef std::shared_ptr< SubsetRelRateUpdater > SharedPtr;

                                            SubsetRelRateUpdater(Model::SharedPtr model);
                                            ~SubsetRelRateUpdater();
        
            virtual void                    debugPriorCalculation();
            virtual double                  calcLogPrior();
            virtual double                  calcLogRefDist();

        private:
        
            void                            pullFromModel();
            void                            pushToModel();

            Model::SharedPtr                _model;
        };

    inline SubsetRelRateUpdater::SubsetRelRateUpdater(Model::SharedPtr model) {
        DirichletUpdater::clear();
        _name = "Subset Relative Rates";
        _model = model;
    }

    inline SubsetRelRateUpdater::~SubsetRelRateUpdater() {
    }

    inline void SubsetRelRateUpdater::debugPriorCalculation() {
        ::om.outputConsole(boost::format("\n%24s %12s %12s %12s %12s\n") % _name % "p_i" % "value" % "parameter" % "log(term)");
        double log_prior = calcLogPrior();
        assert(_curr_point.size() == _prior_parameters.size());

        unsigned num_subsets = _model->getNumSubsets();
        assert(_curr_point.size() == num_subsets);

#if defined(RELRATE_DIRICHLET_PRIOR)
        double sum_log_terms = 0.0;
        double sum_params = 0.0;
        double sum_lgamma_params = 0.0;
        for (unsigned i = 0; i < num_subsets; ++i) {
            double logyi = log(_curr_point[i]);
            double yi = exp(logyi);
            double logterm  = (_prior_parameters[i] - 1.0)*logyi;
            sum_log_terms  += logterm;
            sum_params     += _prior_parameters[i];
            sum_lgamma_params += lgamma(_prior_parameters[i]);
            std::string s   = boost::str(boost::format("m%d") % i);
            ::om.outputConsole(boost::format("%24s %12.5f %12.5f\n") % s % _prior_parameters[i] % logterm);
        }

        double logC = lgamma(sum_params) - sum_lgamma_params;
        sum_log_terms += logC;
        ::om.outputConsole(boost::format("%24s %51.5f\n") % "log-constant" % logC);
        ::om.outputConsole(boost::format("%24s %12.5f\n") % "sum" % sum_log_terms);
        ::om.outputConsole(boost::format("%24s %51.5f\n") % "log-prior" % log_prior);
#else
        Model::subset_sizes_t & subset_sizes = _model->getSubsetSizes();
        double log_num_sites = std::log(_model->getNumSites());

        double sum_log_terms = 0.0;
        double sum_params = 0.0;
        double sum_lgamma_params = 0.0;
        double sum_log_pi = 0.0;
        double sum_pi = 0.0;
        for (unsigned i = 0; i < num_subsets; ++i) {
            double logyi = log(_curr_point[i]);
            double yi = exp(logyi);
            double logpi = log(subset_sizes[i]) - log_num_sites;
            double pi = exp(logpi);
            sum_pi += pi;
            if (i < num_subsets - 1)
                sum_log_pi += logpi;
            double logterm  = (_prior_parameters[i] - 1.0)*logyi;
            sum_log_terms  += logterm;
            sum_params     += _prior_parameters[i];
            sum_lgamma_params += lgamma(_prior_parameters[i]);
            std::string s   = boost::str(boost::format("m%d") % i);
            ::om.outputConsole(boost::format("%24s %12.5f %12.5f %12.5f %12.5f\n") % s % pi % (yi/pi) % _prior_parameters[i] % logterm);
        }

        double logC = sum_log_pi + lgamma(sum_params) - sum_lgamma_params;
        sum_log_terms += logC;
        ::om.outputConsole(boost::format("%24s %51.5f\n") % "log-constant" % logC);
        ::om.outputConsole(boost::format("%24s %12.5f %38.5f\n") % "sum" % sum_pi % sum_log_terms);
        ::om.outputConsole(boost::format("%24s %51.5f\n") % "log-prior" % log_prior);
#endif
    }

    inline double SubsetRelRateUpdater::calcLogPrior() {
        double log_prior = DirichletUpdater::calcLogPrior();
#if defined(RELRATE_DIRICHLET_PRIOR)
        // nothing else to be done
#else
        Model::subset_sizes_t & subset_sizes = _model->getSubsetSizes();
        double log_num_sites = std::log(_model->getNumSites());
        unsigned num_subsets = _model->getNumSubsets();
        double log_jacobian = 0.0;
        for (unsigned i = 0; i < num_subsets - 1; i++) {    //POLMOD2  2022-01-05 (i = 0), 2022-01-30 (i = 1)
            log_jacobian += std::log(subset_sizes[i]);
        }
        log_jacobian -= log_num_sites*(num_subsets - 1);
        log_prior += log_jacobian;
#endif
        return log_prior;
    }

    inline double SubsetRelRateUpdater::calcLogRefDist() {
        double log_refdist = DirichletUpdater::calcLogRefDist();
#if defined(RELRATE_DIRICHLET_PRIOR)
        // nothing else to be done
#else
        Model::subset_sizes_t & subset_sizes = _model->getSubsetSizes();
        unsigned num_subsets = (unsigned)subset_sizes.size();
        double log_num_sites = std::log(_model->getNumSites());
        for (unsigned i = 0; i < num_subsets - 1; i++) {  //POLMOD2 2022-01-05 (i = 0), 2022-01-30 (i = 1)
            log_refdist += std::log(subset_sizes[i]) - log_num_sites;
        }
#endif
        return log_refdist;
    }

#if defined(RELRATE_DIRICHLET_PRIOR)
    inline void SubsetRelRateUpdater::pullFromModel() {
        // Subset relative rate parameters are already Dirichlet distributed
        // Just copy them directly from the model to _curr_point so that the
        // base class DirichletUpdater::update function can update them
        Model::subset_relrate_vect_t & relative_rates = _model->getSubsetRelRates();
        unsigned num_subsets = _model->getNumSubsets();
        _curr_point.resize(num_subsets);
        for (unsigned i = 0; i < num_subsets; i++) {
            _curr_point[i] = relative_rates[i];
        }
    }
    
    inline void SubsetRelRateUpdater::pushToModel() {
        // Subset relative rate parameters are already Dirichlet distributed
        // Just copy them directly from _curr_point (which the base class
        // DirichletUpdater::update function has updated them) back to the model
        _model->setSubsetRelRates(_curr_point, /*fixed*/false);
    }
#else
    inline void SubsetRelRateUpdater::pullFromModel() {
        // Subset relative rate parameters have mean 1, so transform them to a simplex
        // before letting the base class DirichletUpdater::update function operate on them
        Model::subset_relrate_vect_t & relative_rates = _model->getSubsetRelRates();
        Model::subset_sizes_t &        subset_sizes   = _model->getSubsetSizes();
        unsigned num_sites   = _model->getNumSites();
        unsigned num_subsets = _model->getNumSubsets();
        assert(subset_sizes.size() == num_subsets);
        assert(relative_rates.size() == num_subsets);
        _curr_point.resize(num_subsets);
        for (unsigned i = 0; i < num_subsets; i++) {
            _curr_point[i] = relative_rates[i]*subset_sizes[i]/num_sites;
        }
    }
    
    inline void SubsetRelRateUpdater::pushToModel() {
        // Subset relative rate parameters have mean 1, but the base class
        // DirichletUpdater::update function has been operating on a simplex,
        // so translate them back to relative rate parameter space before
        // copying back to the model
        Model::subset_sizes_t & subset_sizes = _model->getSubsetSizes();
        unsigned num_sites   = _model->getNumSites();
        unsigned num_subsets = _model->getNumSubsets();
        point_t tmp(num_subsets);
        for (unsigned i = 0; i < num_subsets; i++) {
            tmp[i] = _curr_point[i]*num_sites/subset_sizes[i];
        }
        _model->setSubsetRelRates(tmp, /*fixed*/false);
    }
#endif
    
}
