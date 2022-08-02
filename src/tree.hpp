#pragma once    

#include <stack>    
#include <memory>
#include <iostream>
#include "node.hpp"

namespace lorad {

    class TreeManip;
    class Likelihood;
    class Updater;
    class TreeUpdater;
    class PolytomyUpdater;  
    class EdgeProportionUpdater;

    class Tree {

            friend class TreeManip;
            friend class Likelihood;
            friend class Updater;
            friend class TreeUpdater;
            friend class PolytomyUpdater;   
            friend class EdgeProportionUpdater;

        public:

                                        Tree();
                                        ~Tree();

            bool                        isRooted() const;
            unsigned                    numLeaves() const;
            unsigned                    numInternals() const;
            unsigned                    numNodes() const;
        
        private:

            void                        clear();

            bool                        _is_rooted;
            Node *                      _root;
            unsigned                    _nleaves;
            unsigned                    _ninternals;
            Node::PtrVector             _preorder;
            Node::PtrVector             _levelorder;
            Node::Vector                _nodes;
            std::vector<Node *>         _unused_nodes;  

        public:

            typedef std::shared_ptr<Tree> SharedPtr;
    };
    

    inline Tree::Tree() {
        clear();
    }

    inline Tree::~Tree() {
    }

    inline void Tree::clear() {
        _is_rooted = false;
        _root = 0;
        _nodes.clear();
        _preorder.clear();
        _levelorder.clear();
    }

    inline bool Tree::isRooted() const {
        return _is_rooted;
    }

    inline unsigned Tree::numLeaves() const {
        return _nleaves;
    }

    inline unsigned Tree::numInternals() const {
        return _ninternals;
    }

    inline unsigned Tree::numNodes() const {
        return (unsigned)_nodes.size();
    }

}
