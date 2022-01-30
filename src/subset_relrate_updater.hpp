#pragma once

#include "dirichlet_updater.hpp"

namespace lorad {
    
    class SubsetRelRateUpdater : public DirichletUpdater {
        public:
            typedef std::shared_ptr< SubsetRelRateUpdater > SharedPtr;

                                            SubsetRelRateUpdater(Model::SharedPtr model);
                                            ~SubsetRelRateUpdater();
        
            virtual double                  calcLogPrior();
            
#if defined(POLGSS)
            virtual double                  calcLogRefDist();
#endif

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

    inline double SubsetRelRateUpdater::calcLogPrior() {
        Model::subset_sizes_t & subset_sizes = _model->getSubsetSizes();
        double log_num_sites = std::log(_model->getNumSites());
        unsigned num_subsets = _model->getNumSubsets();
        double log_prior = DirichletUpdater::calcLogPrior();
        for (unsigned i = 1; i < num_subsets; i++) {
            log_prior += std::log(subset_sizes[i]) - log_num_sites;
        }
        return log_prior;
    }

#if defined(POLGSS)
    inline double SubsetRelRateUpdater::calcLogRefDist() {
        Model::subset_sizes_t & subset_sizes = _model->getSubsetSizes();
        unsigned num_subsets = (unsigned)subset_sizes.size();
        double log_num_sites = std::log(_model->getNumSites());
        double log_refdist = DirichletUpdater::calcLogRefDist();
        for (unsigned i = 1; i < num_subsets; i++) {
            log_refdist += std::log(subset_sizes[i]) - log_num_sites;
        }
        return log_refdist;
    }
#endif

    inline void SubsetRelRateUpdater::pullFromModel() {
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
        Model::subset_sizes_t & subset_sizes = _model->getSubsetSizes();
        unsigned num_sites   = _model->getNumSites();
        unsigned num_subsets = _model->getNumSubsets();
        point_t tmp(num_subsets);
        for (unsigned i = 0; i < num_subsets; i++) {
            tmp[i] = _curr_point[i]*num_sites/subset_sizes[i];
        }
        _model->setSubsetRelRates(tmp, false);
    }
    
}
