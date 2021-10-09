#pragma once

#include "dirichlet_updater.hpp"

namespace strom {

    class EdgeProportionUpdater : public DirichletUpdater {
    
        public:
        
            typedef std::shared_ptr< EdgeProportionUpdater > SharedPtr;

                                        EdgeProportionUpdater();
                                        ~EdgeProportionUpdater();

            virtual double              calcLogPrior();

#if defined(POLGSS)
            double                      calcLogRefDist();
#endif

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
        return Updater::calcLogEdgeLengthPrior().second;
    }

#if defined(POLGSS)
    inline double EdgeProportionUpdater::calcLogRefDist() {
        //double log_refdist = 0.0;
        Tree::SharedPtr tree = _tree_manipulator->getTree();
        assert(tree);

        double TL = _tree_manipulator->calcTreeLength();
        double num_edges = _tree_manipulator->countEdges();
        assert(_refdist_parameters.size() == num_edges);

        // Calculate Dirichlet refdist on edge length proportions
        //
        // Note that, for n edges, the Dirichlet refdist density is
        //
        // p1^{c-1} p2^{c-1} ... pn^{c-1}
        // ------------------------------
        //    Gamma(c)^n / Gamma(n*c)
        //
        // where n = num_edges, pk = edge length k / TL and Gamma is the Gamma function.
        // If c == 1, then both numerator and denominator equal 1, so it is pointless
        // do loop over edge lengths.
        double log_edge_length_proportions_refdist = 0.0;
        unsigned i = 0;
        double csum = 0.0;
        for (auto nd : tree->_preorder) {
            double c = _refdist_parameters[i++];
            csum += c;
            if (c != 1.0) {
                double edge_length_proportion = nd->_edge_length/TL;
                log_edge_length_proportions_refdist += (c - 1.0)*log(edge_length_proportion);
                log_edge_length_proportions_refdist -= std::lgamma(c);
            }
        }
        log_edge_length_proportions_refdist += std::lgamma(csum);
        return log_edge_length_proportions_refdist;
    }
#endif
}
