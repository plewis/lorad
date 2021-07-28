#pragma once

#define POLNEW

#include "dirichlet_updater.hpp"

namespace strom {

#if defined(POLNEW)
    class EdgeProportionUpdater : public DirichletUpdater {
    
        public:
        
            typedef std::shared_ptr< EdgeProportionUpdater > SharedPtr;

                                        EdgeProportionUpdater();
                                        ~EdgeProportionUpdater();

            virtual double              calcLogPrior();

            void                        pullFromModel();
            void                        pushToModel();
        
        private:

            double                      _tree_length;
    };

    // Member function bodies go here
    
    inline EdgeProportionUpdater::EdgeProportionUpdater() {
        // std::cout << "Creating an EdgeProportionUpdater..." << std::endl;
        _tree_length = 0.0;
        _name = "Edge Proportions";
    }

    inline EdgeProportionUpdater::~EdgeProportionUpdater() {
        // std::cout << "Destroying an EdgeProportionUpdater..." << std::endl;
    }

    inline void EdgeProportionUpdater::pullFromModel() {
        _tree_length = _tree_manipulator->copyEdgeProportionsTo(_curr_point);
    }

    inline void EdgeProportionUpdater::pushToModel() {
        _tree_manipulator->copyEdgeProportionsFrom(_tree_length, _curr_point);
    }

    inline double EdgeProportionUpdater::calcLogPrior() {
        return Updater::calcLogEdgeLengthPrior();
    }

#endif

}
