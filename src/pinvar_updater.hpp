#pragma once

#include "model.hpp"
#include "updater.hpp"
#include "asrv.hpp"

namespace lorad {
    
    class PinvarUpdater : public Updater {
    
        public:
            typedef std::shared_ptr<PinvarUpdater> SharedPtr;

                                        PinvarUpdater(ASRV::SharedPtr asrv);
                                        ~PinvarUpdater();
        
            virtual void                clear();
            double                      getCurrentPoint() const;

            // mandatory overrides of pure virtual functions
            virtual double              calcLogPrior();
            double                      calcLogRefDist();
            virtual void                revert();
            virtual void                proposeNewState();
        
        private:
        
            double                      _prev_point;
            ASRV::SharedPtr             _asrv;
    };

    inline PinvarUpdater::PinvarUpdater(ASRV::SharedPtr asrv) {
        clear();
        _name = "Proportion of Invariable Sites";
        assert(asrv);
        _asrv = asrv;
    }

    inline PinvarUpdater::~PinvarUpdater() {
        _asrv.reset();
    }

    inline void PinvarUpdater::clear() {
        Updater::clear();
        _prev_point = 0.0;
        _asrv = nullptr;
        reset();
    }

    inline double PinvarUpdater::getCurrentPoint() const {
        return *(_asrv->getPinvarSharedPtr());
    }
    
    inline double PinvarUpdater::calcLogRefDist() {
        // Assumes Beta(a,b) reference distribution
        assert(_refdist_parameters.size() == 2);
        double refdist_a = _refdist_parameters[0];
        double refdist_b = _refdist_parameters[1];
        
        double log_refdist = 0.0;
        double curr_point = getCurrentPoint();
        if (curr_point > 0.0 && curr_point < 1.0) {
            log_refdist += (refdist_a - 1.0)*std::log(curr_point);
            log_refdist += (refdist_b - 1.0)*std::log(1.0 - curr_point);
            log_refdist += std::lgamma(refdist_a + refdist_b);
            log_refdist -= std::lgamma(refdist_a);
            log_refdist -= std::lgamma(refdist_b);
        }
        else
            log_refdist = Updater::_log_zero;
        return log_refdist;
    }

    inline double PinvarUpdater::calcLogPrior() {
        // Assumes Beta(a,b) prior with mean a/(a+b) and variance a*b/((a + b + 1)*(a + b)^2)
        assert(_prior_parameters.size() == 2);
        double prior_a = _prior_parameters[0];
        double prior_b = _prior_parameters[1];
        
        double log_prior = 0.0;
        double curr_point = getCurrentPoint();
        if (curr_point > 0.0 && curr_point < 1.0) {
            log_prior += (prior_a - 1.0)*std::log(curr_point);
            log_prior += (prior_b - 1.0)*std::log(1.0 - curr_point);
            log_prior += std::lgamma(prior_a + prior_b);
            log_prior -= std::lgamma(prior_a);
            log_prior -= std::lgamma(prior_b);
        }
        else
            log_prior = Updater::_log_zero;
        return log_prior;
    }

    inline void PinvarUpdater::revert() {
        _asrv->setPinvar(_prev_point);
    }

    inline void PinvarUpdater::proposeNewState() {
        if (_lambda > 1.0)
            _lambda = 1.0;
        
        // Save copy of _curr_point in case revert is necessary.
        _prev_point = getCurrentPoint();
        
        // Propose new value using uniform window of width _lambda centered over _prev_point
        double p = (_prev_point - _lambda/2.0) + _lambda*_lot->uniform();
        if (p < 0.0)
            p = -p;
        else if (p > 1.0)
            p = 1.0 - (p - 1.0);
        _asrv->setPinvar(p);
        
        _log_hastings_ratio = 0.0;  //symmetric proposal
        
        // This proposal invalidates all transition matrices and partials
        _tree_manipulator->selectAllPartials();
        _tree_manipulator->selectAllTMatrices();
    }

}
