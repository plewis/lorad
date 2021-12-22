#pragma once

#include "conditionals.hpp"
#include "updater.hpp"

namespace lorad {

    class TreeLengthUpdater : public Updater {
    
        public:
        
            typedef std::shared_ptr< TreeLengthUpdater > SharedPtr;

                                        TreeLengthUpdater();
                                        ~TreeLengthUpdater();

            virtual void                clear();
            virtual void                proposeNewState();
            virtual void                revert();

            virtual double              calcLogPrior();
#if defined(POLGSS)
            double                      calcLogRefDist();
#endif

            void                        pullFromModel();
            void                        pushToModel() const;
        
        private:

            double                      _prev_point;
            double                      _curr_point;
    };

    inline TreeLengthUpdater::TreeLengthUpdater() {
        clear();
        _name = "Tree Length";
    }

    inline TreeLengthUpdater::~TreeLengthUpdater() {
    }

    inline void TreeLengthUpdater::clear() {
        Updater::clear();
        _prev_point     = 0.0;
        _curr_point     = 0.0;
        reset();
    }
    
    inline void TreeLengthUpdater::proposeNewState() {
        // Save copy of _curr_point in case revert is necessary.
        pullFromModel();
        _prev_point = _curr_point;

        // Let _curr_point be proposed value
        double m = exp(_lambda*(_lot->uniform() - 0.5));
        _curr_point = m*_prev_point;
        pushToModel();

#if defined(HOLDER_ETAL_PRIOR)
        // calculate log of Hastings ratio under independent exponential priors
        double num_edges = _tree_manipulator->countEdges();
        _log_hastings_ratio = num_edges*log(m);
#else
        // calculate log of Hastings ratio under GammaDir parameterization
        _log_hastings_ratio = log(m);
#endif

        // This proposal invalidates all transition matrices and partials
        _tree_manipulator->selectAllPartials();
        _tree_manipulator->selectAllTMatrices();
    }

    inline void TreeLengthUpdater::revert() {
        // swap _curr_point and _prev_point so that edge length scaler
        // in pushCurrentStateToModel will be correctly calculated
        double tmp = _curr_point;
        _curr_point = _prev_point;
        _prev_point = tmp;
        pushToModel();
    }

    inline double TreeLengthUpdater::calcLogPrior() {
#if defined(HOLDER_ETAL_PRIOR)
        return Updater::calcLogEdgeLengthPrior();
#else
        return Updater::calcLogEdgeLengthPrior().first;
#endif
    }

#if defined(POLGSS)
    inline double TreeLengthUpdater::calcLogRefDist() {
        Tree::SharedPtr tree = _tree_manipulator->getTree();
        assert(tree);

#if defined(HOLDER_ETAL_PRIOR)
        assert(_refdist_parameters.size() == 1);
        double refdist_b = _refdist_parameters[0];  // rate parameter for independent edge length exponential refdist

        double num_edges = _tree_manipulator->countEdges();
        double curr_point = _tree_manipulator->calcTreeLength();
        double log_refdist = num_edges*std::log(refdist_b) - curr_point*refdist_b;
#else
        assert(_refdist_parameters.size() == 2);
        double refdist_a = _refdist_parameters[0];  // shape parameter for tree length exponential refdist
        double refdist_b = _refdist_parameters[1];  // scale parameter for tree length exponential refdist
        
        double log_refdist = 0.0;
        double curr_point = _tree_manipulator->calcTreeLength();
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
#endif
        return log_refdist;
    }
#endif

    inline void TreeLengthUpdater::pullFromModel() {
        _curr_point = _tree_manipulator->calcTreeLength();
    }

    inline void TreeLengthUpdater::pushToModel() const {
        double scaler = _curr_point/_prev_point;
        _tree_manipulator->scaleAllEdgeLengths(scaler);
    }

}
