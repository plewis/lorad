#pragma once

#include "conditionals.hpp"

#include "model.hpp"
#include "updater.hpp"
#include "asrv.hpp"

#if defined(HOLDER_ETAL_PRIOR)
namespace lorad {
    
    class GammaShapeUpdater : public Updater {
    
        public:
            typedef std::shared_ptr<GammaShapeUpdater> SharedPtr;

                                        GammaShapeUpdater(ASRV::SharedPtr asrv);
                                        ~GammaShapeUpdater();
        
            virtual void                clear();
            double                      getCurrentPoint() const;

            // mandatory overrides of pure virtual functions
            virtual double              calcLogPrior();
            virtual void                revert();
            virtual void                proposeNewState();

            double                      calcLogRefDist();

        private:
        
            double                      _prev_point;
            ASRV::SharedPtr             _asrv;
    };

    inline GammaShapeUpdater::GammaShapeUpdater(ASRV::SharedPtr asrv) {
        clear();
        _name = "Gamma Shape";
        assert(asrv);
        _asrv = asrv;
    }

    inline GammaShapeUpdater::~GammaShapeUpdater() {
        _asrv.reset();
    }

    inline void GammaShapeUpdater::clear() {
        Updater::clear();
        _prev_point = 0.0;
        _asrv = nullptr;
        reset();
    }

    inline double GammaShapeUpdater::getCurrentPoint() const {
        return *(_asrv->getShapeSharedPtr());
    }
    
    inline double GammaShapeUpdater::calcLogRefDist() {
        // Assumes Gamma(a,b) reference distribution with mean a*b and variance a*b^2
        assert(_refdist_parameters.size() == 2);
        double refdist_a = _refdist_parameters[0];
        double refdist_b = _refdist_parameters[1];
        
        double log_refdist = 0.0;
        double curr_point = getCurrentPoint();
        if (curr_point > 0.0) {
            log_refdist += (refdist_a - 1.0)*std::log(curr_point);
            log_refdist -= curr_point/refdist_b;
            log_refdist -= refdist_a*std::log(refdist_b);
            log_refdist -= std::lgamma(refdist_a);
        }
        else if (curr_point == 0.0) {
            if (refdist_a == 1.0) {
                assert(refdist_b > 0.0);
                return -std::log(refdist_b);
            }
            else if (refdist_a > 1.0) {
                log_refdist = Updater::_log_zero;
            }
            else {
                // refdist_a < 1.0
                log_refdist = -Updater::_log_zero;
            }
        }
        else
            log_refdist = Updater::_log_zero;
        return log_refdist;
    }

    inline double GammaShapeUpdater::calcLogPrior() {
        // Assumes Gamma(a,b) prior with mean a*b and variance a*b^2
        assert(_prior_parameters.size() == 2);
        double prior_a = _prior_parameters[0];
        double prior_b = _prior_parameters[1];
        
        double log_prior = 0.0;
        double curr_point = getCurrentPoint();
        if (curr_point > 0.0) {
            log_prior += (prior_a - 1.0)*std::log(curr_point);
            log_prior -= curr_point/prior_b;
            log_prior -= prior_a*std::log(prior_b);
            log_prior -= std::lgamma(prior_a);
        }
        else if (curr_point == 0.0) {
            if (prior_a == 1.0) {
                assert(prior_b > 0.0);
                return -std::log(prior_b);
            }
            else if (prior_a > 1.0) {
                log_prior = Updater::_log_zero;
            }
            else {
                // prior_a < 1.0
                log_prior = -Updater::_log_zero;
            }
        }
        else
            log_prior = Updater::_log_zero;
        return log_prior;
    }

    inline void GammaShapeUpdater::revert() {
        _asrv->setShape(_prev_point);
    }

    inline void GammaShapeUpdater::proposeNewState() {
        // Save copy of _curr_point in case revert is necessary.
        _prev_point = getCurrentPoint();

        // Propose new value using window of width _lambda centered over _prev_point
        double u = _lot->uniform();
        double new_point = (_prev_point - _lambda/2.0) + u*_lambda;
        assert(new_point != 0.0);
        if (new_point < 0.0)
            new_point *= -1.0;
        _asrv->setShape(new_point);
        
        // Calculate log of Hastings ratio
        _log_hastings_ratio = 0.0;  // symmetric proposal
        
        // This proposal invalidates all transition matrices and partials
        _tree_manipulator->selectAllPartials();
        _tree_manipulator->selectAllTMatrices();
    }

}
#endif

