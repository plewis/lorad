#pragma once

#include "conditionals.hpp"
#include "updater.hpp"

#if defined(HOLDER_ETAL_PRIOR)
namespace lorad {

    class EdgeLengthUpdater : public Updater {
    
        public:
        
            typedef std::shared_ptr< EdgeLengthUpdater > SharedPtr;

                                        EdgeLengthUpdater();
                                        ~EdgeLengthUpdater();

            virtual double              calcLogPrior();
            virtual void                proposeNewState();
            virtual void                revert();
            virtual void                reset();

#if defined(POLGSS)
            double                      calcLogRefDist();
#endif

            void                        pullFromModel();
            void                        pushToModel();
            
        private:
        
            Node *                      _focal_node;
            double                      _curr_point;
            double                      _prev_point;
        
    };
    
    inline EdgeLengthUpdater::EdgeLengthUpdater() {
        _focal_node = 0;
        _curr_point = 0.0;
        _prev_point = 0.0;
        _name = "Edge Length";
    }

    inline EdgeLengthUpdater::~EdgeLengthUpdater() {
    }

    inline void EdgeLengthUpdater::proposeNewState() {
        // Choose a random edge to modify
        _focal_node = _tree_manipulator->randomEdge(_lot);
        
        // Save copy of _curr_point in case revert is necessary.
        pullFromModel();
        _prev_point = _curr_point;

        // Let _curr_point be proposed value
        double m = exp(_lambda*(_lot->uniform() - 0.5));
        _curr_point = m*_prev_point;
        pushToModel();

        _log_hastings_ratio = log(m);

        // This proposal invalidates all transition matrices and partials
        _tree_manipulator->selectAllPartials();
        _tree_manipulator->selectAllTMatrices();
        
        //TODO: try these versions once working
        //_tree_manipulator->selectPartialsHereToRoot(_focal_node);
        //_focal_node->selectTMatrix();
    }

    inline void EdgeLengthUpdater::revert() {
        _curr_point = _prev_point;
        pushToModel();
    }

    inline void EdgeLengthUpdater::reset() {
        Updater::reset();
        _curr_point = 0.0;
        _prev_point = 0.0;
        _focal_node = 0;
    }

    inline void EdgeLengthUpdater::pullFromModel() {
        assert(_focal_node);
        _curr_point = _focal_node->getEdgeLength();
    }

    inline void EdgeLengthUpdater::pushToModel() {
        _focal_node->setEdgeLength(_curr_point);
        //TODO: invalidate from this node back to root
    }

    inline double EdgeLengthUpdater::calcLogPrior() {
        return Updater::calcLogEdgeLengthPrior();
        //TODO: just compute prior for this one edge length?
    }

#if defined(POLGSS)
    inline double EdgeLengthUpdater::calcLogRefDist() {
#if defined(HOLDER_ETAL_PRIOR)
        // Assumes Exp(r) reference distribution with rate r
        assert(_refdist_parameters.size() == 1);
        assert(_refdist_parameters[0] > 0.0);
        double refdist_a = 1.0;
        double refdist_b = 1.0/_refdist_parameters[0];
#else
        // Assumes Gamma(a,b) reference distribution with mean a*b and variance a*b^2
        assert(_refdist_parameters.size() == 2);
        double refdist_a = _refdist_parameters[0];
        double refdist_b = _refdist_parameters[1];
#endif
        
        double log_refdist = 0.0;
        if (_curr_point > 0.0) {
            log_refdist += (refdist_a - 1.0)*std::log(_curr_point);
            log_refdist -= _curr_point/refdist_b;
            log_refdist -= refdist_a*std::log(refdist_b);
            log_refdist -= std::lgamma(refdist_a);
        }
        else if (_curr_point == 0.0) {
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
}
#endif

