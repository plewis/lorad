#pragma once

#define POLNEW

#include "updater.hpp"

namespace strom {

#if defined(POLNEW)
    class EdgeLengthUpdater : public Updater {
    
        public:
        
            typedef std::shared_ptr< EdgeLengthUpdater > SharedPtr;

                                        EdgeLengthUpdater();
                                        ~EdgeLengthUpdater();

            virtual void                clear();
            virtual void                proposeNewState();
            virtual void                revert();

            virtual double              calcLogPrior();

            void                        pullFromModel();
            void                        pushToModel() const;
        
        private:

            double                      _prev_edgelen;
            double                      _curr_edgelen;
            Node *                      _focal_node;
    };

    // Member function bodies go here
    
    inline EdgeLengthUpdater::EdgeLengthUpdater() {
        // std::cout << "Creating an EdgeLengthUpdater..." << std::endl;
        clear();
        _name = "Edge Length";
    }

    inline EdgeLengthUpdater::~EdgeLengthUpdater() {
        // std::cout << "Destroying an EdgeLengthUpdater..." << std::endl;
    }

    inline void EdgeLengthUpdater::clear() {
        Updater::clear();
        _focal_node       = 0;
        _prev_edgelen     = 0.0;
        _curr_edgelen     = 0.0;
        reset();
    }
    
    inline void EdgeLengthUpdater::pullFromModel() {
        _focal_node = _tree_manipulator->randomEdge(_lot);
        _curr_edgelen = _focal_node->getEdgeLength();
    }

    inline void EdgeLengthUpdater::pushToModel() const {
        _focal_node->setEdgeLength(_curr_edgelen);
    }

    inline void EdgeLengthUpdater::proposeNewState() {
        // Save copy of _curr_edgelen in case revert is necessary.
        pullFromModel();
        _prev_edgelen = _curr_edgelen;

        // Let _curr_edgelen be proposed value
        double m = exp(_lambda*(_lot->uniform() - 0.5));
        _curr_edgelen = m*_prev_edgelen;
        pushToModel();

        // calculate log of Hastings ratio under GammaDir parameterization
        _log_hastings_ratio = log(m);

        // This proposal invalidates partials from _focal_node down to the root of the tree
        _tree_manipulator->selectPartialsHereToRoot(_focal_node);
        
        // This proposal invalidates the transition matrix of the _focal_node only
        _focal_node->selectTMatrix();
    }

    inline void EdgeLengthUpdater::revert() {
        _curr_edgelen = _prev_edgelen;
        pushToModel();
    }

    inline double EdgeLengthUpdater::calcLogPrior() {
        return Updater::calcLogEdgeLengthPrior();
    }

#endif

}
