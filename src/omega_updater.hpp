#pragma once

#include "model.hpp"
#include "updater.hpp"
#include "qmatrix.hpp"

namespace strom {
    
    class OmegaUpdater : public Updater {
    
        public:
            typedef std::shared_ptr<OmegaUpdater> SharedPtr;

                                        OmegaUpdater(QMatrix::SharedPtr q);
                                        ~OmegaUpdater();
        
            virtual void                clear();
            double                      getCurrentPoint() const;

            // mandatory overrides of pure virtual functions
            virtual double              calcLogPrior();
#if defined(POLGSS)
            double                      calcLogRefDist();
#endif
            virtual void                revert();
            virtual void                proposeNewState();
        
        private:
        
            double                      _prev_point;
            QMatrix::SharedPtr          _q;
    };

    // member function bodies go here    
    inline OmegaUpdater::OmegaUpdater(QMatrix::SharedPtr q) {
        //std::cout << "OmegaUpdater being created" << std::endl;
        clear();
        _name = "Omega";
        assert(q);
        _q = q;
    }

    inline OmegaUpdater::~OmegaUpdater() {
        //std::cout << "OmegaUpdater being destroyed" << std::endl;
        _q.reset();
    }

    inline void OmegaUpdater::clear() {
        Updater::clear();
        _prev_point = 0.0;
        _q = nullptr;
        reset();
    }

    inline double OmegaUpdater::getCurrentPoint() const {
        return *(_q->getOmegaSharedPtr());
    }
    
#if defined(POLGSS)
    inline double OmegaUpdater::calcLogRefDist() {
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
#endif

    inline double OmegaUpdater::calcLogPrior() {
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
        else
            log_prior = Updater::_log_zero;
        return log_prior;
    }

    inline void OmegaUpdater::revert() {
        _q->setOmega(_prev_point);
    }

    inline void OmegaUpdater::proposeNewState() {
        // Save copy of _curr_point in case revert is necessary.
        _prev_point = getCurrentPoint();
        
        // Propose new value using multiplier with boldness _lambda
        double m = exp(_lambda*(_lot->uniform() - 0.5));
        _q->setOmega(m*_prev_point);
        
        // Calculate log of Hastings ratio
        _log_hastings_ratio = log(m);
        
        // This proposal invalidates all transition matrices and partials
        _tree_manipulator->selectAllPartials();
        _tree_manipulator->selectAllTMatrices();
    }

}
