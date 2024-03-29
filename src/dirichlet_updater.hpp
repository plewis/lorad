#pragma once

#include "updater.hpp"

namespace lorad {
    
    class Chain;
    
    class DirichletUpdater : public Updater {

        friend class Chain;
        
        public:
            typedef std::vector<double>                 point_t;
            typedef std::shared_ptr< DirichletUpdater > SharedPtr;

                                                DirichletUpdater();
            virtual                             ~DirichletUpdater();
        
            void                                clear();
            virtual void                        debugPriorCalculation();
            virtual double                      calcLogPrior();
            double                              calcLogRefDist();
        
        protected:
        
            virtual void                        pullFromModel() = 0;
            virtual void                        pushToModel() = 0;

            void                                proposeNewState();
            void                                revert();
        
            point_t                             _curr_point;
            point_t                             _prev_point;
    };
    
    inline DirichletUpdater::DirichletUpdater() {
        clear();
    }

    inline DirichletUpdater::~DirichletUpdater() {
    }
    
    inline void DirichletUpdater::clear() {
        Updater::clear();
        _prev_point.clear();
    }
    
    inline double DirichletUpdater::calcLogRefDist() {
        pullFromModel();
        assert(_curr_point.size() > 0);
        assert(_curr_point.size() == _refdist_parameters.size());
        double log_refdist = 0.0;
        double refdist_param_sum = 0.0;
        for (unsigned i = 0; i < _curr_point.size(); ++i) {
            log_refdist += (_refdist_parameters[i] - 1.0)*std::log(_curr_point[i]);
            log_refdist -= std::lgamma(_refdist_parameters[i]);
            refdist_param_sum += _refdist_parameters[i];
        }
        log_refdist += std::lgamma(refdist_param_sum);
        return log_refdist;
    }

    inline void DirichletUpdater::debugPriorCalculation() {
        double log_prior = calcLogPrior();
        ::om.outputConsole(boost::format("%s:\n") % _name);
        ::om.outputConsole("  value:\n");
        for (unsigned i = 0; i < _curr_point.size(); ++i) {
            ::om.outputConsole(boost::format("  %12.5f\n") % _curr_point[i]);
        }
        ::om.outputConsole("  parameters:\n");
        for (unsigned i = 0; i < _prior_parameters.size(); ++i) {
            ::om.outputConsole(boost::format("  %12.5f\n") % _prior_parameters[i]);
        }
        ::om.outputConsole("  log-prior:\n");
        ::om.outputConsole(boost::format("  %12.5f\n") % log_prior);
    }

#define POL_2022_09_02

    inline double DirichletUpdater::calcLogPrior() {
        pullFromModel();
        assert(_curr_point.size() > 0);
        assert(_curr_point.size() == _prior_parameters.size());
#if defined(POL_2022_09_02)
        // playing it safe and just calculating full prior each time
#else
        bool flat_prior = true;
#endif
        bool bad_point = false;
        double log_prior = 0.0;
        double prior_param_sum = 0.0;
        for (unsigned i = 0; i < _curr_point.size(); ++i) {
#if defined(POL_2022_09_02)
        // playing it safe and just calculating full prior each time
#else
            if (_prior_parameters[i] != 1.0)
                flat_prior = false;
#endif
            if (_curr_point[i] == 0.0)
                bad_point = true;
            log_prior += (_prior_parameters[i] - 1.0)*std::log(_curr_point[i]);
            log_prior -= std::lgamma(_prior_parameters[i]);
            prior_param_sum += _prior_parameters[i];
        }
#if defined(POL_2022_09_02)
        // playing it safe and just calculating full prior each time
        if (bad_point)
            return Updater::_log_zero;
        else
            log_prior += std::lgamma(prior_param_sum);
#else
        if (flat_prior)
            return std::lgamma(prior_param_sum);
        else if (bad_point)
            return Updater::_log_zero;
        else
            log_prior += std::lgamma(prior_param_sum);
#endif
        return log_prior;
    }

    inline void DirichletUpdater::proposeNewState() {
        // Save length of _curr_point.
        pullFromModel();
        unsigned dim = (unsigned)_curr_point.size();
        
        // Save copy of _curr_point in case revert is necessary.
        _prev_point.assign(_curr_point.begin(), _curr_point.end());
        
        // Determine parameters of Dirichlet forward proposal distribution and, at the same time,
        // draw gamma deviates that will be used to form the proposed point.
        std::vector<double> forward_params(dim, 0.0);
        for (unsigned i = 0; i < dim; ++i) {
            // Calculate ith forward parameter
            double alpha_i = 1.0 + _prev_point[i]/_lambda;
            if (alpha_i < 1.e-12)
                alpha_i = 1.e-12;
            forward_params[i] = alpha_i;
            
            // Draw ith gamma deviate
            _curr_point[i] = 0.0;
            if (alpha_i > 0.0)
                _curr_point[i] = _lot->gamma(alpha_i, 1.0);
        }
        
        double sum_gamma_deviates     = std::accumulate(_curr_point.begin(), _curr_point.end(), 0.0);
        double sum_forward_parameters = std::accumulate(forward_params.begin(), forward_params.end(), 0.0);

        // Choose new state by sampling from forward proposal distribution.
        // We've already stored gamma deviates in _curr_point, now just need to normalize them.
        for (unsigned i = 0; i < dim; ++i) {
            _curr_point[i] /= sum_gamma_deviates;
        }
        
        // Determine probability density of the forward proposal
        double log_forward_density = 0.0;
        for (unsigned i = 0; i < dim; ++i) {
            log_forward_density += (forward_params[i] - 1.0)*std::log(_prev_point[i]);
            log_forward_density -= std::lgamma(forward_params[i]);
        }
        log_forward_density += std::lgamma(sum_forward_parameters);
        
        // Determine parameters of Dirichlet reverse proposal distribution
        std::vector<double> reverse_params(dim, 0.0);
        for (unsigned i = 0; i < dim; ++i) {
            reverse_params[i] = 1.0 + _curr_point[i]/_lambda;
        }
        
        double sum_reverse_parameters = std::accumulate(reverse_params.begin(), reverse_params.end(), 0.0);

        // determine probability density of the reverse proposal
        double log_reverse_density = 0.0;
        for (unsigned i = 0; i < dim; ++i) {
            log_reverse_density += (reverse_params[i] - 1.0)*std::log(_curr_point[i]);
            log_reverse_density -= std::lgamma(reverse_params[i]);
        }
        log_reverse_density += std::lgamma(sum_reverse_parameters);
        
        // calculate the logarithm of the Hastings ratio
        _log_hastings_ratio = log_reverse_density - log_forward_density;
        
        pushToModel();

        // This proposal invalidates all transition matrices and partials
        _tree_manipulator->selectAllPartials();
        _tree_manipulator->selectAllTMatrices();
    }
    
    inline void DirichletUpdater::revert() {
        std::copy(_prev_point.begin(), _prev_point.end(), _curr_point.begin());
        pushToModel();
    }

}
