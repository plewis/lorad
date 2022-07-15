#pragma once    

#include "conditionals.hpp"

#include <cassert>
#include <memory>
#include <stack>
#include <queue>
#include <set>
#if defined(USE_BOOST_REGEX)
#   include <boost/regex.hpp>
#else
#   include <regex>
#endif
#include <boost/range/adaptor/reversed.hpp>
#include <boost/format.hpp>
#include <Eigen/Dense>
#include "tree.hpp"
#include "lot.hpp"
#include "conditional_clade_store.hpp"
#include "xlorad.hpp"

namespace lorad {

    class TreeManip {
            
        public:
        
                                        TreeManip();
                                        TreeManip(Tree::SharedPtr t);
                                        ~TreeManip();

            void                        setTree(Tree::SharedPtr t);
            Tree::SharedPtr             getTree();

            double                      calcTreeLength() const;
            unsigned                    calcResolutionClass() const;    
            unsigned                    countEdges() const;
            unsigned                    countInternals() const; 
            void                        scaleAllEdgeLengths(double scaler);
            
            void                        createTestTree();
            std::string                 makeNewick(unsigned precision, bool use_names = false) const;

            void                        buildFromNewick(const std::string newick, bool rooted, bool allow_polytomies);
            void                        storeSplits(std::set<Split> & splitset, bool include_trivial_splits = false);
            void                        storeClades(ConditionalCladeStore::SharedPtr ccs);
            double                      calcEmpiricalCladeProb(ConditionalCladeStore::SharedPtr ccs);
            double                      calcLogReferenceCladeProb(ConditionalCladeStore::SharedPtr ccs);
            void                        rerootAtNodeNumber(int node_number);
        
            Node *                      randomEdge(Lot::SharedPtr lot);
            Node *                      randomInternalEdge(Lot::SharedPtr lot);
            Node *                      randomChild(Lot::SharedPtr lot, Node * x, Node * avoid, bool parent_included);
            void                        LargetSimonSwap(Node * a, Node * b);
            
            bool                        isPolytomy(Node * nd) const;   
            void                        nniNodeSwap(Node * a, Node * b); 
            unsigned                    countChildren(Node * nd) const;
            Node *                      findLeftSib(Node * nd);
            Node *                      findRightmostChild(Node * nd);
            Node *                      findLastPreorderInClade(Node * start);
            void                        insertSubtreeOnLeft(Node * s, Node * u);
            void                        insertSubtreeOnRight(Node * s, Node * u);
            void                        detachSubtree(Node * s);
            void                        rectifyNumInternals(int incr);
            void                        refreshNavigationPointers();
            Node *                      getUnusedNode(Node * sought = 0);
            void                        putUnusedNode(Node * nd);      
            
            void                        selectAll();
            void                        deselectAll();
            void                        selectAllPartials();
            void                        deselectAllPartials();
            void                        selectAllTMatrices();
            void                        deselectAllTMatrices();
            void                        selectPartialsHereToRoot(Node * a);
            void                        flipPartialsAndTMatrices();

            void                        clear();
            
            void                        edgeLengthsAsString(std::string & receptacle, unsigned precision = 5, char separator = '\t') const;
            double                      copyEdgeLengthsTo(std::vector<double> & receptacle) const;
            void                        copyEdgeLengthsFrom(const std::vector<double> & new_edgelens);
            
            void                        edgeProportionsAsString(std::string & receptacle, unsigned precision = 5, char separator = '\t') const;
            double                      copyEdgeProportionsTo(std::vector<double> & receptacle) const;
            void                        copyEdgeProportionsFrom(double TL, const std::vector<double> & new_props);
            
            void                        saveParamNames(std::vector<std::string> & param_name_vect) const;
            double                      logTransformEdgeLengths(std::vector<double> & param_vect) const;
            double                      setEdgeLengthsFromLogTransformed(Eigen::VectorXd & param_vect, double TL, unsigned first, unsigned nedges);

#if defined(POLGSS)
            void                        sampleTree();
            std::string                 calcExpRefDist(std::string title, std::vector<double> & vect, std::vector<double> & params);
            std::string                 calcGammaRefDist(std::string title, std::vector<double> & vect, std::vector<double> & params);
            std::string                 calcDirichletRefDist(std::string title, std::vector< std::vector<double> > & vect, std::vector<double> & params);
            std::string                 saveReferenceDistributions();
            std::string                 calcReferenceDistributions(std::map<std::string, std::vector<double> > & refdist_map);
#endif

#if defined(POLGHM)
            void                        setModelToSampledPoint(unsigned i);
#endif

        private:

#if defined(POLGSS) || defined(POLGHM)
#if defined(HOLDER_ETAL_PRIOR)
            std::vector<double>                 _sampled_edge_lengths;
#else
            std::vector< std::vector<double> >  _sampled_edge_proportions;
            std::vector<double>                 _sampled_tree_lengths;
#endif
#endif

            Node *                      findNextPreorder(Node * nd);
            void                        refreshPreorder();
            void                        refreshLevelorder();
            void                        renumberInternals();
            void                        rerootAtNode(Node * prospective_root);
            void                        extractNodeNumberFromName(Node * nd, std::set<unsigned> & used);
            void                        extractEdgeLen(Node * nd, std::string edge_length_string);
            unsigned                    countNewickLeaves(const std::string newick);
            void                        stripOutNexusComments(std::string & newick);
            bool                        canHaveSibling(Node * nd, bool rooted, bool allow_polytomies);

            Tree::SharedPtr             _tree;

        public:

            typedef std::shared_ptr< TreeManip > SharedPtr;
    };

    inline TreeManip::TreeManip() {
        clear();
    }

    inline TreeManip::TreeManip(Tree::SharedPtr t) {
        clear();
        setTree(t);
    }

    inline TreeManip::~TreeManip() {
    }

    inline void TreeManip::clear() {
#if defined(POLGSS)
#if defined(HOLDER_ETAL_PRIOR)
        _sampled_edge_lengths.clear();
#else
        _sampled_edge_proportions.clear();
        _sampled_tree_lengths.clear();
#endif
#endif
        _tree.reset();
    }

    inline void TreeManip::setTree(Tree::SharedPtr t) {
        assert(t);
        _tree = t;
    }

    inline Tree::SharedPtr TreeManip::getTree() {
        return _tree;
    }

    inline double TreeManip::calcTreeLength() const {
        double TL = 0.0;
        for (auto nd : _tree->_preorder) {
            TL += nd->_edge_length;
        }
        return TL;
    }

    inline unsigned TreeManip::calcResolutionClass() const {    
        return _tree->_ninternals;
    }   

    inline unsigned TreeManip::countEdges() const {
        return (unsigned)_tree->_preorder.size();
    }

    inline unsigned TreeManip::countInternals() const { 
        unsigned m = 0;
        for (auto nd : _tree->_preorder) {
            if (nd->_left_child)
                m++;
        }
        return m;
    }   

    inline unsigned TreeManip::countChildren(Node * nd) const { 
        assert(nd);
        unsigned nchildren = 0;
        Node * child = nd->getLeftChild();
        while (child) {
            nchildren++;
            child = child->getRightSib();
        }
        return nchildren;
    }   

    inline void TreeManip::scaleAllEdgeLengths(double scaler) {
        for (auto nd : _tree->_preorder) {
            nd->setEdgeLength(scaler*nd->_edge_length);
        }
    }

    inline void TreeManip::createTestTree() {
        clear();
        _tree = Tree::SharedPtr(new Tree());
        _tree->_nodes.resize(6);

        Node * root_node       = &_tree->_nodes[0];
        Node * first_internal  = &_tree->_nodes[1];
        Node * second_internal = &_tree->_nodes[2];
        Node * first_leaf      = &_tree->_nodes[3];
        Node * second_leaf     = &_tree->_nodes[4];
        Node * third_leaf      = &_tree->_nodes[5];

        // Here is the structure of the tree (numbers in
        // parentheses are node numbers, other numbers
        // are edge lengths):
        //
        // first_leaf (0)   second_leaf (1)   third_leaf (2)
        //      \              /                  /
        //       \ 0.1        / 0.1              /
        //        \          /                  /
        //     second_internal (3)             / 0.2
        //             \                      /
        //              \ 0.1                /
        //               \                  /
        //                first_internal (4)
        //                        |
        //                        | 0.1
        //                        |
        //                    root_node (5)
        //
        root_node->_parent = 0;
        root_node->_left_child = first_internal;
        root_node->_right_sib = 0;
        root_node->_number = 5;
        root_node->_name = "root_node";
        root_node->_edge_length = 0.0;

        first_internal->_parent = root_node;
        first_internal->_left_child = second_internal;
        first_internal->_right_sib = 0;
        first_internal->_number = 4;
        first_internal->_name = "first_internal_node";
        first_internal->_edge_length = 0.1;

        second_internal->_parent = first_internal;
        second_internal->_left_child = first_leaf;
        second_internal->_right_sib = third_leaf;
        second_internal->_number = 3;
        second_internal->_name = "second_internal_node";
        second_internal->_edge_length = 0.1;

        first_leaf->_parent = second_internal;
        first_leaf->_left_child = 0;
        first_leaf->_right_sib = second_leaf;
        first_leaf->_number = 0;
        first_leaf->_name = "first_leaf";
        first_leaf->_edge_length = 0.1;

        second_leaf->_parent = second_internal;
        second_leaf->_left_child = 0;
        second_leaf->_right_sib = 0;
        second_leaf->_number = 1;
        second_leaf->_name = "second_leaf";
        second_leaf->_edge_length = 0.1;

        third_leaf->_parent = first_internal;
        third_leaf->_left_child = 0;
        third_leaf->_right_sib = 0;
        third_leaf->_number = 2;
        third_leaf->_name = "third_leaf";
        third_leaf->_edge_length = 0.2;

        _tree->_is_rooted = true;
        _tree->_root = root_node;
        _tree->_nleaves = 3;

        // Note that root node is not included in _preorder
        _tree->_preorder.push_back(first_internal);
        _tree->_preorder.push_back(second_internal);
        _tree->_preorder.push_back(first_leaf);
        _tree->_preorder.push_back(second_leaf);
        _tree->_preorder.push_back(third_leaf);

        _tree->_levelorder.push_back(first_internal);
        _tree->_levelorder.push_back(second_internal);
        _tree->_levelorder.push_back(third_leaf);
        _tree->_levelorder.push_back(first_leaf);
        _tree->_levelorder.push_back(second_leaf);
    }
        
    inline std::string TreeManip::makeNewick(unsigned precision, bool use_names) const {
        std::string newick;
        const boost::format tip_node_name_format( boost::str(boost::format("%%s:%%.%df") % precision) );
        const boost::format tip_node_number_format( boost::str(boost::format("%%d:%%.%df") % precision) );
        const boost::format internal_node_format( boost::str(boost::format("):%%.%df") % precision) );
        std::stack<Node *> node_stack;

        Node * root_tip = (_tree->_is_rooted ? 0 : _tree->_root);
        for (auto nd : _tree->_preorder) {
            if (nd->_left_child) {
                newick += "(";
                node_stack.push(nd);
                if (root_tip) {
                    if (use_names) {
                        newick += boost::str(boost::format(tip_node_name_format)
                            % root_tip->_name
                            % nd->_edge_length);
                    } else {
                        newick += boost::str(boost::format(tip_node_number_format)
                            % (root_tip->_number + 1)
                            % nd->_edge_length);
                    }
                    newick += ",";
                    root_tip = 0;
                }
            }
            else {
                if (use_names) {
                    newick += boost::str(boost::format(tip_node_name_format)
                        % nd->_name
                        % nd->_edge_length);
                } else {
                    newick += boost::str(boost::format(tip_node_number_format)
                        % (nd->_number + 1)
                        % nd->_edge_length);
                }
                if (nd->_right_sib)
                    newick += ",";
                else {
                    Node * popped = (node_stack.empty() ? 0 : node_stack.top());
                    while (popped && !popped->_right_sib) {
                        node_stack.pop();
                        if (node_stack.empty()) {
                            newick += ")";
                            popped = 0;
                        }
                        else {
                            newick += boost::str(boost::format(internal_node_format) % popped->_edge_length);
                            popped = node_stack.top();
                        }
                    }
                    if (popped && popped->_right_sib) {
                        node_stack.pop();
                        newick += boost::str(boost::format(internal_node_format) % popped->_edge_length);
                        newick += ",";
                    }
                }
            }
        }

        return newick;
    }

    inline void TreeManip::extractNodeNumberFromName(Node * nd, std::set<unsigned> & used) {
        assert(nd);
        bool success = true;
        unsigned x = 0;
        try {
            x = std::stoi(nd->_name);
        }
        catch(std::invalid_argument &) {
            // node name could not be converted to an integer value
            success = false;
        }

        if (success) {
            // conversion succeeded
            // attempt to insert x into the set of node numbers already used
            std::pair<std::set<unsigned>::iterator, bool> insert_result = used.insert(x);
            if (insert_result.second) {
                // insertion was made, so x has NOT already been used
                nd->_number = x - 1;
            }
            else {
                // insertion was not made, so set already contained x
                throw XLorad(boost::str(boost::format("leaf number %d used more than once") % x));
            }
        }
        else
            throw XLorad(boost::str(boost::format("node name (%s) not interpretable as a positive integer") % nd->_name));
    }

    inline void TreeManip::extractEdgeLen(Node * nd, std::string edge_length_string) {
        assert(nd);
        bool success = true;
        double d = 0.0;
        try {
            d = std::stod(edge_length_string);
        }
        catch(std::invalid_argument &) {
            // edge_length_string could not be converted to a double value
            success = false;
        }

        if (success) {
            // conversion succeeded
            nd->setEdgeLength(d);
        }
        else
            throw XLorad(boost::str(boost::format("%s is not interpretable as an edge length") % edge_length_string));

    }

    inline unsigned TreeManip::countNewickLeaves(const std::string newick) {
#if defined(USE_BOOST_REGEX)
        boost::regex taxonexpr("[(,]\\s*(\\d+|\\S+?|['].+?['])\\s*(?=[,):])");
        boost::sregex_iterator m1(newick.begin(), newick.end(), taxonexpr);
        boost::sregex_iterator m2;
#else
        std::regex taxonexpr("[(,]\\s*(\\d+|\\S+?|['].+?['])\\s*(?=[,):])");
        std::sregex_iterator m1(newick.begin(), newick.end(), taxonexpr);
        std::sregex_iterator m2;
#endif
        return (unsigned)std::distance(m1, m2);
    }

    inline void TreeManip::stripOutNexusComments(std::string & newick) {
#if defined(USE_BOOST_REGEX)
        boost::regex commentexpr("\\[.*?\\]");
        newick = boost::regex_replace(newick, commentexpr, std::string(""));
#else
        std::regex commentexpr("\\[.*?\\]");
        newick = std::regex_replace(newick, commentexpr, std::string(""));
#endif
    }
    
    inline void TreeManip::refreshPreorder() {  
        // Create vector of node pointers in preorder sequence
        _tree->_preorder.clear();
        _tree->_preorder.reserve(_tree->_nodes.size() - 1); // _preorder does not include root node

        if (!_tree->_root)
            return;

        Node * first_preorder = _tree->_root->_left_child;

        // sanity check: first preorder node should be the only child of the root node
        assert(first_preorder->_right_sib == 0);

        Node * nd = first_preorder;
        _tree->_preorder.push_back(nd);

        while (true) {
            nd = findNextPreorder(nd);  
            if (nd)
                _tree->_preorder.push_back(nd);
            else
                break;  
        }   // end while loop
    }   

    //                            1. start by adding only descendant of root node to buffer queue
    //                               queue = [1], stack = []
    //                            2. move node at front of buffer queue to back of stack vector
    //                               queue = [], stack = [1]
    //                 8    9     3. add this node's immediate children to back of buffer queue
    //                  \  /         queue = [2,3], stack = [1]
    //                   \/       4. repeat 2 and 3 until all nodes are processed
    //      4  5    6    7           (2) queue = [3], stack = [1,2]
    //       \ |     \  /            (3) queue = [3,4,5], stack = [1,2]
    //        \|      \/             (2) queue = [4,5], stack = [1,2,3]
    //         2      3              (3) queue = [4,5,6,7], stack = [1,2,3]
    //          \    /               (2) queue = [5,6,7], stack = [1,2,3,4]
    //           \  /                (3) no-op: 4 has no children
    //            \/                 (2) queue = [6,7], stack = [1,2,3,4,5]
    //            1                  (3) no-op: 5 has no children
    //            |                  (2) queue = [7], stack = [1,2,3,4,5,6]
    //            0                  (3) no-op: 6 has no children
    //                               (2) queue = [], stack = [1,2,3,4,5,6,7]
    //                               (3) queue = [8,9], stack = [1,2,3,4,5,6,7]
    //                               (2) queue = [9], stack = [1,2,3,4,5,6,7,8]
    //                               (3) no-op: 8 has no children
    //                               (2) queue = [], stack = [1,2,3,4,5,6,7,8,9]
    //                               (3) no-op: 9 has no children
    //                            5. stack vector is now in level order
    inline void TreeManip::refreshLevelorder() {
        if (!_tree->_root)
            return;

        // q is the buffer queue
        std::queue<Node *> q;

        // _tree->_levelorder is the stack vector
        _tree->_levelorder.clear();
        _tree->_levelorder.reserve(_tree->_nodes.size() - 1);

        Node * nd = _tree->_root->_left_child;

        // sanity check: first node should be the only child of the root node
        assert(nd->_right_sib == 0);

        // Push nd onto back of queue
        q.push(nd);

        while (!q.empty()) {
            // pop nd off front of queue
            nd = q.front(); q.pop();

            // and push it onto the stack
            _tree->_levelorder.push_back(nd);

            // add all children of nd to back of queue
            Node * child = nd->_left_child;
            if (child) {
                q.push(child);
                child = child->_right_sib;
                while (child) {
                    q.push(child);
                    child = child->_right_sib;
                }
            }
        }   // end while loop
    }

    inline void TreeManip::renumberInternals() {    
        assert(_tree->_preorder.size() > 0);
        
        // Renumber internal nodes in postorder sequence
        unsigned curr_internal = _tree->_nleaves;
        for (auto nd : boost::adaptors::reverse(_tree->_preorder)) {
            if (nd->_left_child) {
                // nd is an internal node
                nd->_number = curr_internal++;
            }
        }
        
        // Root node is not included in _tree->_preorder, so if the root node
        // is an internal node we need to number it here
        if (_tree->_is_rooted)
            _tree->_root->_number = curr_internal++;
            
        _tree->_ninternals = curr_internal - _tree->_nleaves;
        
        _tree->_unused_nodes.clear();   
        for (; curr_internal < (unsigned)_tree->_nodes.size(); curr_internal++) {
            Node * curr = &_tree->_nodes[curr_internal];
            putUnusedNode(curr);
            assert(curr->_number == -1);
            curr->_number = curr_internal;
        }   
        
    }   
    
    inline bool TreeManip::canHaveSibling(Node * nd, bool rooted, bool allow_polytomies) {
        assert(nd);
        if (!nd->_parent) {
            // trying to give root node a sibling
            return false;
        }

        if (allow_polytomies)
            return true;

        bool nd_can_have_sibling = true;
        if (nd != nd->_parent->_left_child) {
            if (nd->_parent->_parent) {
                // trying to give a sibling to a sibling of nd, and nd's parent is not the root
                nd_can_have_sibling = false;
            }
            else {
                if (rooted) {
                    // root node has exactly 2 children in rooted trees
                    nd_can_have_sibling = false;
                }
                else if (nd != nd->_parent->_left_child->_right_sib) {
                    // trying to give root node more than 3 children
                    nd_can_have_sibling = false;
                }
            }
        }

        return nd_can_have_sibling;
    }

    inline void TreeManip::rerootAtNodeNumber(int node_number) {
        // Locate node having _number equal to node_number
        Node * nd = 0;
        for (auto & curr : _tree->_nodes) {
            if (curr._number == node_number) {
                nd = &curr;
                break;
            }
        }

        if (!nd)
            throw XLorad(boost::str(boost::format("no node found with number equal to %d") % node_number));

        if (nd != _tree->_root) {
            if (nd->_left_child)
                throw XLorad(boost::str(boost::format("cannot currently root trees at internal nodes (e.g. node %d)") % nd->_number));
            rerootAtNode(nd);
        }
    }

    inline void TreeManip::rerootAtNode(Node * prospective_root) {
        Node * a = prospective_root;
        Node * b = prospective_root->_parent;
        Node * c = 0;
        Node * d = 0;
        Node * p = 0;
        a->_parent = 0;
        double tmp_edgelen  = 0.0;
        double prev_edgelen = a->getEdgeLength();

        while (b) {
            // Prune node a from b
            if (a == b->_left_child) {
                if (a->_right_sib) {
                    b->_left_child = a->_right_sib;
                    a->_right_sib = 0;
                }
                else {
                    b->_left_child = 0;
                }
            }
            else {
                c = b->_left_child;
                while (c->_right_sib != a)
                    c = c->_right_sib;
                d = a->_right_sib;
                c->_right_sib = d;
                a->_right_sib = 0;
            }

            // Graft node b onto node a (but don't unhook node b from its parent just yet)
            if (a->_left_child) {
                c = a->_left_child;
                while (c->_right_sib)
                    c = c->_right_sib;
                c->_right_sib = b;
            }
            else {
                a->_left_child = b;
            }

            // Rotate
            p = a;
            a = b;
            b = b->_parent;
            a->_parent = p;

            // Swap nd's edge length with its new parent's edge length
            tmp_edgelen = a->getEdgeLength();
            a->setEdgeLength(prev_edgelen);
            prev_edgelen = tmp_edgelen;
        }
        prospective_root->setEdgeLength(0.0);
        _tree->_root = prospective_root;
        refreshPreorder();
        refreshLevelorder();
    }

    inline void TreeManip::buildFromNewick(const std::string newick, bool rooted, bool allow_polytomies) {
        _tree.reset(new Tree());
        _tree->_is_rooted = rooted;

        std::set<unsigned> used; // used to ensure that no two leaf nodes have the same number
        unsigned curr_leaf = 0;
        unsigned num_edge_lengths = 0;
        unsigned curr_node_index = 0;

        // Remove comments from the supplied newick string
        std::string commentless_newick = newick;
        stripOutNexusComments(commentless_newick);

        // Resize the _nodes vector
        _tree->_nleaves = countNewickLeaves(commentless_newick);
        if (_tree->_nleaves < 4)
            throw XLorad("Expecting newick tree description to have at least 4 leaves");
        unsigned max_nodes = 2*_tree->_nleaves - (rooted ? 0 : 2);
        _tree->_nodes.resize(max_nodes);
        for (auto & nd : _tree->_nodes )
            nd._number = -1;

        try {
            // Root node
            Node * nd = &_tree->_nodes[curr_node_index];
            _tree->_root = nd;

            if (_tree->_is_rooted) {
                nd = &_tree->_nodes[++curr_node_index];
                nd->_parent = &_tree->_nodes[curr_node_index - 1];
                nd->_parent->_left_child = nd;
            }

            // Some flags to keep track of what we did last
            enum {
                Prev_Tok_LParen		= 0x01,	// previous token was a left parenthesis ('(')
                Prev_Tok_RParen		= 0x02,	// previous token was a right parenthesis (')')
                Prev_Tok_Colon		= 0x04,	// previous token was a colon (':')
                Prev_Tok_Comma		= 0x08,	// previous token was a comma (',')
                Prev_Tok_Name		= 0x10,	// previous token was a node name (e.g. '2', 'P._articulata')
                Prev_Tok_EdgeLen	= 0x20	// previous token was an edge length (e.g. '0.1', '1.7e-3')
            };
            unsigned previous = Prev_Tok_LParen;

            // Some useful flag combinations
            unsigned LParen_Valid = (Prev_Tok_LParen | Prev_Tok_Comma);
            unsigned RParen_Valid = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
            unsigned Comma_Valid  = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
            unsigned Colon_Valid  = (Prev_Tok_RParen | Prev_Tok_Name);
            unsigned Name_Valid   = (Prev_Tok_RParen | Prev_Tok_LParen | Prev_Tok_Comma);

            // Set to true while reading an edge length
            bool inside_edge_length = false;
            std::string edge_length_str;
            unsigned edge_length_position = 0;

            // Set to true while reading a node name surrounded by (single) quotes
            bool inside_quoted_name = false;

            // Set to true while reading a node name not surrounded by (single) quotes
            bool inside_unquoted_name = false;

            // Set to start of each node name and used in case of error
            unsigned node_name_position = 0;

            // loop through the characters in newick, building up tree as we go
            unsigned position_in_string = 0;
            for (auto ch : commentless_newick) {
                position_in_string++;

                if (inside_quoted_name) {
                    if (ch == '\'') {
                        inside_quoted_name = false;
                        node_name_position = 0;
                        if (!nd->_left_child) {
                            extractNodeNumberFromName(nd, used);
                            curr_leaf++;
                        }
                        previous = Prev_Tok_Name;
                    }
                    else if (iswspace(ch))
                        nd->_name += ' ';
                    else
                        nd->_name += ch;

                    continue;
                }
                else if (inside_unquoted_name) {
                    if (ch == '(')
                        throw XLorad(boost::str(boost::format("Unexpected left parenthesis inside node name at position %d in tree description") % node_name_position));

                    if (iswspace(ch) || ch == ':' || ch == ',' || ch == ')') {
                        inside_unquoted_name = false;

                        // Expect node name only after a left paren (child's name), a comma (sib's name) or a right paren (parent's name)
                        if (!(previous & Name_Valid))
                            throw XLorad(boost::str(boost::format("Unexpected node name (%s) at position %d in tree description") % nd->_name % node_name_position));

                        if (!nd->_left_child) {
                            extractNodeNumberFromName(nd, used);
                            curr_leaf++;
                        }

                        previous = Prev_Tok_Name;
                    }
                    else {
                        nd->_name += ch;
                        continue;
                    }
                }
                else if (inside_edge_length) {
                    if (ch == ',' || ch == ')' || iswspace(ch)) {
                        inside_edge_length = false;
                        edge_length_position = 0;
                        extractEdgeLen(nd, edge_length_str);
                        ++num_edge_lengths;
                        previous = Prev_Tok_EdgeLen;
                    }
                    else {
                        bool valid = (ch =='e' || ch == 'E' || ch =='.' || ch == '-' || ch == '+' || isdigit(ch));
                        if (!valid)
                            throw XLorad(boost::str(boost::format("Invalid branch length character (%c) at position %d in tree description") % ch % position_in_string));
                        edge_length_str += ch;
                        continue;
                    }
                }

                if (iswspace(ch))
                    continue;

                switch(ch) {
                    case ';':
                        break;

                    case ')':
                        // If nd is bottommost node, expecting left paren or semicolon, but not right paren
                        if (!nd->_parent)
                            throw XLorad(boost::str(boost::format("Too many right parentheses at position %d in tree description") % position_in_string));

                        // Expect right paren only after an edge length, a node name, or another right paren
                        if (!(previous & RParen_Valid))
                            throw XLorad(boost::str(boost::format("Unexpected right parenthesisat position %d in tree description") % position_in_string));

                        // Go down a level
                        nd = nd->_parent;
                        if (!nd->_left_child->_right_sib)
                            throw XLorad(boost::str(boost::format("Internal node has only one child at position %d in tree description") % position_in_string));
                        previous = Prev_Tok_RParen;
                        break;

                    case ':':
                        // Expect colon only after a node name or another right paren
                        if (!(previous & Colon_Valid))
                            throw XLorad(boost::str(boost::format("Unexpected colon at position %d in tree description") % position_in_string));
                        previous = Prev_Tok_Colon;
                        break;

                    case ',':
                        // Expect comma only after an edge length, a node name, or a right paren
                        if (!nd->_parent || !(previous & Comma_Valid))
                            throw XLorad(boost::str(boost::format("Unexpected comma at position %d in tree description") % position_in_string));

                        // Check for polytomies
                        if (!canHaveSibling(nd, rooted, allow_polytomies)) {
                            throw XLorad(boost::str(boost::format("Polytomy found in the following tree description but polytomies prohibited:\n%s") % newick));
                        }

                        // Create the sibling
                        curr_node_index++;
                        if (curr_node_index == _tree->_nodes.size())
                            throw XLorad(boost::str(boost::format("Too many nodes specified by tree description (%d nodes allocated for %d leaves)") % _tree->_nodes.size() % _tree->_nleaves));
                        nd->_right_sib = &_tree->_nodes[curr_node_index];
                        nd->_right_sib->_parent = nd->_parent;
                        nd = nd->_right_sib;
                        previous = Prev_Tok_Comma;
                        break;

                    case '(':
                        // Expect left paren only after a comma or another left paren
                        if (!(previous & LParen_Valid))
                            throw XLorad(boost::str(boost::format("Not expecting left parenthesis at position %d in tree description") % position_in_string));

                        // Create new node above and to the left of the current node
                        assert(!nd->_left_child);
                        curr_node_index++;
                        if (curr_node_index == _tree->_nodes.size())
                            throw XLorad(boost::str(boost::format("malformed tree description (more than %d nodes specified)") % _tree->_nodes.size()));
                        nd->_left_child = &_tree->_nodes[curr_node_index];
                        nd->_left_child->_parent = nd;
                        nd = nd->_left_child;
                        previous = Prev_Tok_LParen;
                        break;

                    case '\'':
                        // Encountered an apostrophe, which always indicates the start of a
                        // node name (but note that node names do not have to be quoted)

                        // Expect node name only after a left paren (child's name), a comma (sib's name)
                        // or a right paren (parent's name)
                        if (!(previous & Name_Valid))
                            throw XLorad(boost::str(boost::format("Not expecting node name at position %d in tree description") % position_in_string));

                        // Get the rest of the name
                        nd->_name.clear();

                        inside_quoted_name = true;
                        node_name_position = position_in_string;

                        break;

                    default:
                        // Get here if ch is not one of ();:,'

                        // Expecting either an edge length or an unquoted node name
                        if (previous == Prev_Tok_Colon) {
                            // Edge length expected (e.g. "235", "0.12345", "1.7e-3")
                            inside_edge_length = true;
                            edge_length_position = position_in_string;
                            edge_length_str = ch;
                        }
                        else {
                            // Get the node name
                            nd->_name = ch;

                            inside_unquoted_name = true;
                            node_name_position = position_in_string;
                        }
                }   // end of switch statement
            }   // loop over characters in newick string

            if (inside_unquoted_name)
                throw XLorad(boost::str(boost::format("Tree description ended before end of node name starting at position %d was found") % node_name_position));
            if (inside_edge_length)
                throw XLorad(boost::str(boost::format("Tree description ended before end of edge length starting at position %d was found") % edge_length_position));
            if (inside_quoted_name)
                throw XLorad(boost::str(boost::format("Expecting single quote to mark the end of node name at position %d in tree description") % node_name_position));

            if (_tree->_is_rooted) {
                refreshPreorder();
                refreshLevelorder();
            }
            else {
                // Root at leaf whose _number = 0
                // refreshPreorder() and refreshLevelorder() called after rerooting
                rerootAtNodeNumber(0);
            }
            renumberInternals();
        }
        catch(XLorad & x) {
            clear();
            throw x;
        }
    }

    inline void TreeManip::storeSplits(std::set<Split> & splitset, bool include_trivial_splits) {
        // Start by clearing and resizing all splits
        for (auto & nd : _tree->_nodes) {
            nd._split.resize(_tree->_nleaves);
        }

        // Now do a postorder traversal and add the bit corresponding
        // to the current node in its parent node's split
        for (auto nd : boost::adaptors::reverse(_tree->_preorder)) {
            if (nd->_left_child) {
                // add this internal node's split to splitset
                splitset.insert(nd->_split);
            }
            else {
                // set bit corresponding to this leaf node's number
                nd->_split.setBitAt(nd->_number);
                
                if (include_trivial_splits) {
                    splitset.insert(nd->_split);
                }
            }

            if (nd->_parent) {
                // parent's bits are the union of the bits set in all its children
                nd->_parent->_split.addSplit(nd->_split);
            }
        }
    }
    
    inline double TreeManip::calcEmpiricalCladeProb(ConditionalCladeStore::SharedPtr ccs) {
        // Performs preorder traversal, mutiplying together all non-trivial conditional clade
        // probabilities
        double prob = 1.0;
        
        std::set<Split> splitset;
        storeSplits(splitset);
        
        for (auto nd : _tree->_preorder) {
            Node * lchild = nd->_left_child;
            if (lchild) {
                // nd is internal
                Node * rchild = lchild->_right_sib;
                
                // we assume that there is a right child too
                assert(rchild);
                
                // assume no polytomies
                assert(!rchild->_right_sib);
                
                // find largest of the two clades
                Split & lsplit = lchild->_split;
                Split & rsplit = rchild->_split;
                unsigned lcount = lsplit.countBitsSet();
                unsigned rcount = rsplit.countBitsSet();
                if (lcount > 1 || rcount > 1) {
                    bool left_larger = (lcount > rcount) || ((lcount == rcount) && (rsplit < lsplit));
                    if (left_larger) {
                        prob *= ccs->getEmpiricalProb(nd->_split, lchild->_split);
                    }
                    else {
                        prob *= ccs->getEmpiricalProb(nd->_split, rchild->_split);
                    }
                }
            }
        }
        
        return prob;
    }
    
    inline double TreeManip::calcLogReferenceCladeProb(ConditionalCladeStore::SharedPtr ccs) {
        // Performs preorder traversal, mutiplying together all non-trivial conditional clade
        // probabilities
        double log_prob = 0.0;
        
        std::set<Split> splitset;
        storeSplits(splitset);
        
        for (auto nd : _tree->_preorder) {
            Node * lchild = nd->_left_child;
            if (lchild) {
                // nd is internal
                Node * rchild = lchild->_right_sib;
                
                // we assume that there is a right child too
                assert(rchild);
                
                // assume no polytomies
                assert(!rchild->_right_sib);
                
                // find largest of the two clades
                Split & lsplit = lchild->_split;
                Split & rsplit = rchild->_split;
                unsigned lcount = lsplit.countBitsSet();
                unsigned rcount = rsplit.countBitsSet();
                if (lcount > 1 || rcount > 1) {
                    bool left_larger = (lcount > rcount) || ((lcount == rcount) && (rsplit < lsplit));
                    double reference_prob = 0.0;
                    if (left_larger) {
                        reference_prob = ccs->getReferenceProb(nd->_split, lchild->_split);
                    }
                    else {
                        reference_prob = ccs->getReferenceProb(nd->_split, rchild->_split);
                    }
                    assert(reference_prob > 0.0);
                    log_prob += log(reference_prob);
                }
            }
        }
        
        return log_prob;
    }
    
    inline void TreeManip::storeClades(ConditionalCladeStore::SharedPtr ccs) {
        // Performs a preorder traversal to add conditional clades to ccs.
        // Assumes storeSplits has already been called.
        // Assumes tree is binary and rooted at tip 0.
        //
        //    01000  00100  00010  00001
        //        2     3     4     5
        //         \   /     /     /
        //          \ /     /     /
        //    01100  8     /     /
        //            \   /     /
        //             \ /     /
        //       01110  7     /           ccs->_parent_map[01110] = 01100
        //               \   /
        //                \ /
        //          01111  6              ccs->_parent_map[01111] = 01110
        //                 |
        //                 |
        //          11111  1
        
        // Preorder traversal to add internal nodes with at least one internal child to ccs
        for (auto nd : _tree->_preorder) {
            Node * lchild = nd->_left_child;
            if (lchild) {
                // nd is internal, so make an entry for the largest of the two child clades
                // or, if the child clades are the same size, the one with the greater split
                Node * rchild = lchild->_right_sib;
                assert(rchild);
                assert(!rchild->_right_sib);
                Split & lsplit = lchild->_split;
                Split & rsplit = rchild->_split;
                unsigned lcount = lsplit.countBitsSet();
                unsigned rcount = rsplit.countBitsSet();
                if (lcount > 1 || rcount > 1) {
                    // only bother storing conditional clades if there is a possibility of variation
                    bool add_left = (lcount > rcount) || ((lcount == rcount) && (rsplit < lsplit));
                    if (add_left)
                        ccs->addParentChildSplit(nd->_split, lchild->_split);
                    else
                        ccs->addParentChildSplit(nd->_split, rchild->_split);
                }
            }
        }
        
        // Add entry for root node
        //Node * root = _tree->_root;
        //assert(root);
        //Node * subroot = _tree->_root->_left_child;
        //assert(subroot);
        //ccs->addParentChildSplit(root->_split, subroot->_split);
    }

    inline Node * TreeManip::randomEdge(Lot::SharedPtr lot) {
        // Unrooted case:                        Rooted case:
        //
        // 2     3     4     5                   1     2     3     4
        //  \   /     /     /                     \   /     /     /
        //   \ /     /     /                       \ /     /     /
        //    8     /     /                         7     /     /
        //     \   /     /                           \   /     /
        //      \ /     /                             \ /     /
        //       7     /                               6     /
        //        \   /                                 \   /
        //         \ /                                   \ /
        //          6   nleaves = 5                       5    nleaves = 4
        //          |   preorder length = 7               |    preorder length = 7
        //          |   num_edges = 7 - 0 = 7             |    num_edges = 7 - 1 = 6
        //          1   choices: 2,3,4,5,6,7,8           root  choices: 1,2,3,4,6,7
        //
        // _preorder = [6, 7, 8, 2, 3, 4, 5]     _preorder = [5, 6, 7, 1, 2, 3, 4]
        //
        // Note: _preorder is actually a vector of Node *, but is shown here as a
        // vector of integers solely to illustrate the algorithm below.

        // If tree is rooted, add one to skip first node in _preorder vector (and to subtract
        // one from number of valid edges) because, for rooted tree case, this first noe is
        // an internal node whose edge is not a valid choice.
        int rooted_offset = (_tree->_is_rooted ? 1 : 0);
        
        int num_edges = (unsigned)_tree->_preorder.size() - rooted_offset;
        double uniform_deviate = lot->uniform();
        unsigned index_of_chosen = rooted_offset + (unsigned)std::floor(uniform_deviate*num_edges);

        unsigned nodes_visited = 0;
        Node * chosen_node = 0;
        for (auto nd : _tree->_preorder) {
            if (nodes_visited == index_of_chosen) {
                chosen_node = nd;
                break;
            }
            else
                ++nodes_visited;
        }
        assert(chosen_node);
        return chosen_node;
    }

    inline Node * TreeManip::randomInternalEdge(Lot::SharedPtr lot) {   
        // Unrooted case:                        Rooted case:
        //
        // 2     3     4     5                   1     2     3     4
        //  \   /     /     /                     \   /     /     /
        //   \ /     /     /                       \ /     /     /
        //    8     /     /                         7     /     /
        //     \   /     /                           \   /     /
        //      \ /     /                             \ /     /
        //       7     /                               6     /
        //        \   /                                 \   /
        //         \ /                                   \ /
        //          6   nleaves = 5                       5    nleaves = 4
        //          |   preorder length = 7               |    preorder length = 7
        //          |   num_internal_edges = 7 - 5 = 2    |    num_internal_edges = 7 - 4 - 1 = 2
        //          1   choose node 7 or node 8          root  choose node 6 or node 7
        //
        // _preorder = [6, 7, 8, 2, 3, 4, 5]     _preorder = [5, 6, 7, 1, 2, 3, 4]
        //
        // Note: _preorder is actually a vector of T *, but is shown here as a
        // vector of integers solely to illustrate the algorithm below
        
        int num_internal_edges = (unsigned)_tree->_preorder.size() - _tree->_nleaves - (_tree->_is_rooted ? 1 : 0);
        if (num_internal_edges == 0) {  
            // Star tree: return hub node, which is the first node in the preorder sequence
            return _tree->_preorder[0];
        }   

        // Add one to skip first node in _preorder vector, which is an internal node whose edge
        // is either a terminal edge (if tree is unrooted) or invalid (if tree is rooted)
        double uniform_deviate = lot->uniform();
        unsigned index_of_chosen = 1 + (unsigned)std::floor(uniform_deviate*num_internal_edges);

        unsigned internal_nodes_visited = 0;
        Node * chosen_node = 0;
        for (auto nd : _tree->_preorder) {
            if (nd->_left_child) {
                if (internal_nodes_visited == index_of_chosen) {
                    chosen_node = nd;
                    break;
                }
                else
                    ++internal_nodes_visited;
            }
        }
        assert(chosen_node);
        return chosen_node;
    }   

    inline Node * TreeManip::randomChild(Lot::SharedPtr lot, Node * x, Node * avoid, bool parent_included) {
        // Count number of children of x
        unsigned n = 0;
        Node * child = x->getLeftChild();
        while (child) {
            if (child != avoid)
                n++;
            child = child->getRightSib();
    }

        // Choose random child index
        unsigned upper = n + (parent_included ? 1 : 0);
        unsigned chosen = lot->randint(0,upper - 1);
        
        // If chosen < n, then find and return that particular child
        if (chosen < n) {
            child = x->getLeftChild();
            unsigned i = 0;
            while (child) {
                if (child != avoid && i == chosen)
                    return child;
                else if (child != avoid)
                    i++;
                child = child->getRightSib();
            }
        }

        // If chosen equals n, then the parent was chosen, indicated by returning NULL
        return NULL;
    }

    inline void TreeManip::LargetSimonSwap(Node * a, Node * b) {
        // a and b are the ends of the selected 3-edge path in a Larget-Simon move
        // The 3-edge path is indicated by parentheses around the nodes involved.
        // x is always the parent of a
        // y can be the parent of b (case 1) or the child of b (case 2)
        
        Node * x = a->_parent;
        assert(x);
        
        Node * y = x->_parent;
        assert(y);
        
        if (y == b->_parent) {
            // Case 1: y is the parent of b
            //
            //    (a) d  e             (b) d  e
            //      \ | /                \ | /
            //       \|/                  \|/
            //       (x) f (b)            (x) f (a)    Swap a and b, leaving everything
            //         \ | /                \ | /      else as is
            //          \|/     ==>          \|/
            //          (y)                  (y)
            //           |                    |
            //           |                    |
            //           c                    c
            //

            // Detach a from tree
            if (a == x->_left_child) {
                x->_left_child = a->_right_sib;
            } else {
                Node * child = x->_left_child;
                while (child->_right_sib != a)
                    child = child->_right_sib;
                child->_right_sib = a->_right_sib;
            }
            a->_parent = 0;
            a->_right_sib = 0;
            
            // Detach b from tree
            if (b == y->_left_child) {
                y->_left_child = b->_right_sib;
            } else {
                Node * child = y->_left_child;
                while (child->_right_sib != b)
                    child = child->_right_sib;
                child->_right_sib = b->_right_sib;
            }
            b->_parent = 0;
            b->_right_sib = 0;

            // Reattach a to y
            a->_right_sib = y->_left_child;
            y->_left_child = a;
            a->_parent = y;
            
            // Reattach b to x
            b->_right_sib = x->_left_child;
            x->_left_child = b;
            b->_parent = x;
        }
        else {
            // Case 2: y is the child of b
            //
            //    (a) d  e             (a) f  c
            //      \ | /                \ | /
            //       \|/                  \|/
            //       (x) f  c            (x) d  e    swap x's children (except a)
            //         \ | /               \ | /     with y's children (except x)
            //          \|/     ==>         \|/
            //          (y)                 (y)
            //           |                   |
            //           |                   |
            //          (b)                 (b)
            assert(b == y->_parent);
            
            // Remove x's children from tree and store in xchildren stack
            std::stack<Node *> xchildren;
            Node * child = x->_left_child;
            Node * prevchild = 0;
            while (child) {
                if (child == a) {
                    prevchild = child;
                    child = child->_right_sib;
                } else {
                    if (child == x->_left_child) {
                        x->_left_child = child->_right_sib;
                        child->_right_sib = 0;
                        child->_parent = 0;
                        xchildren.push(child);
                        child = x->_left_child;
                    } else if (child->_right_sib) {
                        prevchild->_right_sib = child->_right_sib;
                        child->_right_sib = 0;
                        child->_parent = 0;
                        xchildren.push(child);
                        child = prevchild->_right_sib;
                    } else {
                        assert(prevchild == a);
                        a->_right_sib = 0;
                        child->_parent = 0;
                        xchildren.push(child);
                        child = 0;
                        prevchild = 0;
                    }
                }
            }
            
            // Remove y's children from tree and store in ychildren stack
            std::stack<Node *> ychildren;
            child = y->_left_child;
            prevchild = 0;
            while (child) {
                if (child == x) {
                    prevchild = child;
                    child = child->_right_sib;
                } else {
                    if (child == y->_left_child) {
                        y->_left_child = child->_right_sib;
                        child->_right_sib = 0;
                        child->_parent = 0;
                        ychildren.push(child);
                        child = y->_left_child;
                    } else if (child->_right_sib) {
                        prevchild->_right_sib = child->_right_sib;
                        child->_right_sib = 0;
                        child->_parent = 0;
                        ychildren.push(child);
                        child = prevchild->_right_sib;
                    } else {
                        assert(prevchild == x);
                        x->_right_sib = 0;
                        child->_parent = 0;
                        ychildren.push(child);
                        child = 0;
                        prevchild = 0;
                    }
                }
            }
            
            // Reattach xchildren to y
            while (!xchildren.empty()) {
                Node * popped = xchildren.top();
                xchildren.pop();
                popped->_right_sib = y->_left_child;
                y->_left_child = popped;
                popped->_parent = y;
            }

            // Reattach ychildren to x
            while (!ychildren.empty()) {
                Node * popped = ychildren.top();
                ychildren.pop();
                popped->_right_sib = x->_left_child;
                x->_left_child = popped;
                popped->_parent = x;
            }
        }
        
        refreshPreorder();
        refreshLevelorder();
    }
    
    inline void TreeManip::selectAll() {
        for (auto & nd : _tree->_nodes) {
            nd.select();
        }
    }

    inline void TreeManip::deselectAll() {
        for (auto & nd : _tree->_nodes) {
            nd.deselect();
        }
    }

    inline void TreeManip::selectAllPartials() {
        for (auto & nd : _tree->_nodes)
            nd.selectPartial();
    }

    inline void TreeManip::deselectAllPartials() {
        for (auto & nd : _tree->_nodes) {
            nd.deselectPartial();
        }
    }

    inline void TreeManip::selectAllTMatrices() {
        for (auto & nd : _tree->_nodes)
            nd.selectTMatrix();
    }

    inline void TreeManip::deselectAllTMatrices() {
        for (auto & nd : _tree->_nodes) {
            nd.deselectTMatrix();
        }
    }

    inline void TreeManip::selectPartialsHereToRoot(Node * a) {
        a->selectPartial();
        while (a->_parent) {
            a = a->_parent;
            a->selectPartial();
        }
    }

    inline void TreeManip::flipPartialsAndTMatrices() {
        for (auto & nd : _tree->_nodes) {
            if (nd.isSelPartial())
                nd.flipPartial();
            
            if (nd.isSelTMatrix())
                nd.flipTMatrix();
        }
    }
    
    inline Node * TreeManip::findNextPreorder(Node * nd) {  
        assert(nd);
        Node * next = 0;
        if (!nd->_left_child && !nd->_right_sib) {
            // nd has no children and no siblings, so next preorder is the right sibling of
            // the first ancestral node that has a right sibling.
            Node * anc = nd->_parent;
            while (anc && !anc->_right_sib)
                anc = anc->_parent;
            if (anc) {
                // We found an ancestor with a right sibling
                next = anc->_right_sib;
            }
            else {
                // nd is last preorder node in the tree
                next = 0;
            }
        }
        else if (nd->_right_sib && !nd->_left_child) {
            // nd has no children (it is a tip), but does have a sibling on its right
            next = nd->_right_sib;
        }
        else if (nd->_left_child && !nd->_right_sib) {
            // nd has children (it is an internal node) but no siblings on its right
            next = nd->_left_child;
        }
        else {
            // nd has both children and siblings on its right
            next = nd->_left_child;
        }
        return next;
    }   
    
    inline Node * TreeManip::findLeftSib(Node * nd) {   
        assert(nd);
        assert(nd->_parent);
        Node * child = nd->_parent->_left_child;
        while (child && child->_right_sib != nd)
            child = child->_right_sib;
        return child;
    }   
    
    inline Node * TreeManip::findRightmostChild(Node * nd) {    
        assert(nd);
        Node * child = nd->getLeftChild();
        while (child->getRightSib())
            child = child->getRightSib();
        return child;
    }   
    
    inline Node * TreeManip::findLastPreorderInClade(Node * start) {    
        assert(start);
        Node * curr = start;
        Node * rchild = findRightmostChild(curr);
        while (rchild) {
            curr = rchild;
            rchild = findRightmostChild(curr);
        }
        return curr;
    }   
    
    inline void TreeManip::insertSubtreeOnLeft(Node * s, Node * u) {    
        assert(u);
        assert(s);
        s->_right_sib  = u->_left_child;
        s->_parent     = u;
        u->_left_child = s;
    }   

    inline void TreeManip::insertSubtreeOnRight(Node * s, Node * u) {   
        assert(u);
        assert(s);

        s->_right_sib = 0;
        s->_parent    = u;
        if (u->_left_child) {
            Node * u_rchild = findRightmostChild(u);
            u_rchild->_right_sib = s;
        }
        else
            u->_left_child = s;
    }   
    
    inline void TreeManip::detachSubtree(Node * s) {    
        assert(s);
        assert(s->_parent);
        
        // Save pointers to relevant nodes
        Node * s_leftsib  = findLeftSib(s);
        Node * s_rightsib = s->_right_sib;
        Node * s_parent   = s->_parent;

        // Completely detach s and seal up the wound
        s->_parent = 0;
        s->_right_sib = 0;
        if (s_leftsib)
            s_leftsib->_right_sib = s_rightsib;
        else
            s_parent->_left_child = s_rightsib;
    }   
    
    inline void TreeManip::rectifyNumInternals(int incr) {  
        assert(_tree->_nodes.size() == _tree->_unused_nodes.size() + _tree->_nleaves + _tree->_ninternals + incr);
        _tree->_ninternals += incr;
    }   
    
    inline void TreeManip::refreshNavigationPointers() {    
        refreshPreorder();
        refreshLevelorder();
    }   
    
    inline Node * TreeManip::getUnusedNode(Node * sought) {  
        assert(!_tree->_unused_nodes.empty());
        Node * nd = 0;
        if (sought) {
            unsigned i = 0;
            for (Node * und : _tree->_unused_nodes) {
                if (und == sought) {
                    nd = und;
                    _tree->_unused_nodes.erase(_tree->_unused_nodes.begin()+i);
                    break;
                }
                ++i;
            }
            assert(nd);
        }
        else {
            nd = _tree->_unused_nodes.back();
            _tree->_unused_nodes.pop_back();
        }
        nd->clearPointers();
        return nd;
    }   
    
    inline void TreeManip::putUnusedNode(Node * nd) {   
        nd->clearPointers();
        _tree->_unused_nodes.push_back(nd);
    }   
    
    inline bool TreeManip::isPolytomy(Node * nd) const {   
        Node * lchild = nd->_left_child;
        assert(lchild);    // should only call this function for internal nodes
        
        Node * rchild = lchild->_right_sib;
        if (rchild && rchild->_right_sib)
            return true;
        return false;
    }   
    
    inline void TreeManip::copyEdgeLengthsFrom(const std::vector<double> & new_edgelens) {
        assert(new_edgelens.size() == _tree->_preorder.size());
        unsigned i = 0;
        for (auto nd : _tree->_preorder) {
            nd->setEdgeLength(new_edgelens[i++]);
        }
    }
    
    inline void TreeManip::edgeLengthsAsString(std::string & receptacle, unsigned precision, char separator) const {
        receptacle = "";
        std::string fmtstr = boost::str(boost::format("%%.%df%s") % precision % separator);
        for (auto nd : _tree->_preorder) {
            assert(nd->_edge_length > 0.0);
            receptacle += boost::str(boost::format(fmtstr) % nd->_edge_length);
        }
    }

    inline double TreeManip::copyEdgeLengthsTo(std::vector<double> & receptacle) const {
        receptacle.resize(_tree->_preorder.size());
        unsigned i = 0;
        double TL = 0.0;
        for (auto nd : _tree->_preorder) {
            assert(nd->_edge_length > 0.0);
            receptacle[i++] = nd->_edge_length;
            TL += nd->_edge_length;
        }

        return TL;
    }

    inline void TreeManip::copyEdgeProportionsFrom(double TL, const std::vector<double> & new_props) {
        assert(new_props.size() == _tree->_preorder.size());
        unsigned i = 0;
        for (auto nd : _tree->_preorder) {
            nd->setEdgeLength(TL*new_props[i++]);
        }
    }
    
    inline void TreeManip::edgeProportionsAsString(std::string & receptacle, unsigned precision, char separator) const {
        // Compute and store edge proportions in temporary vector tmp
        std::vector<double> tmp(_tree->_preorder.size());
        unsigned i = 0;
        double TL = 0.0;
        for (auto nd : _tree->_preorder) {
            assert(nd->_edge_length > 0.0);
            double v = nd->_edge_length;
            tmp[i++] = v;
            TL += v;
        }
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), [TL](double v){return v/TL;});

        // Concatenate edge proportions onto receptacle
        std::string fmtstr = boost::str(boost::format("%%.%df%s") % precision % separator);
        receptacle = "";
        for (double v : tmp) {
            receptacle += boost::str(boost::format(fmtstr) % v);
        }
    }

    inline double TreeManip::copyEdgeProportionsTo(std::vector<double> & receptacle) const {
        receptacle.resize(_tree->_preorder.size());
        unsigned i = 0;
        double TL = 0.0;
        for (auto nd : _tree->_preorder) {
            assert(nd->_edge_length > 0.0);
            receptacle[i++] = nd->_edge_length;
            TL += nd->_edge_length;
        }

        // Convert edge lengths into proportions by dividing each by TL
        std::transform(receptacle.begin(), receptacle.end(), receptacle.begin(), [TL](double v) {return v/TL;});
        
        return TL;
    }

    inline void TreeManip::saveParamNames(std::vector<std::string> & param_name_vect) const {
#if defined(HOLDER_ETAL_PRIOR)
        unsigned num_edgelens = (unsigned)(_tree->_preorder.size());
        for (unsigned i = 0; i < num_edgelens; ++i) {
            param_name_vect.push_back(boost::str(boost::format("v-%d") % (i+1)));
        }
#else
        param_name_vect.push_back("TL");
        
        // the number of transformed edge length proportions is one less than the number of edges
        unsigned num_edgelen_proportions = (unsigned)(_tree->_preorder.size() - 1);
        for (unsigned i = 0; i < num_edgelen_proportions; ++i) {
            param_name_vect.push_back(boost::str(boost::format("edgeprop-%d") % (i+1)));
        }
#endif
    }
    
    inline double TreeManip::logTransformEdgeLengths(std::vector<double> & param_vect) const {
#if defined(HOLDER_ETAL_PRIOR)
#   if defined(LORAD_VARIABLE_TOPOLOGY)
        // Record each split and log of associated edge length in a vector
        std::vector< std::pair<Split, double> > split_logedgelen_vect;
        double log_jacobian = 0.0;
        for (auto nd : _tree->_preorder) {
            assert(nd->_edge_length > 0.0);
            double logv = log(nd->_edge_length);
            log_jacobian += logv;
            split_logedgelen_vect.push_back(std::make_pair(nd->getSplit(),logv));
        }
        
        // Sort split_logedgelen_vect by split (this works because split occurs first in each pair)
        std::sort(split_logedgelen_vect.begin(), split_logedgelen_vect.end());
        
        // Save log edge lengths in param_vect in order of sorted splits so that identical tree topologies
        // will have edge length parameters always saved in the same order
        param_vect.resize(split_logedgelen_vect.size());
        unsigned i = 0;
        for (auto & split_logv : split_logedgelen_vect) {
            param_vect[i++] = split_logv.second;
        }
#   else
        // Record each edge length in param_vect and compute the log of the Jacobian
        double log_jacobian = 0.0;
        for (auto nd : _tree->_preorder) {
            assert(nd->_edge_length > 0.0);
            double logv = log(nd->_edge_length);
            log_jacobian += logv;
            param_vect.push_back(logv);
        }
#   endif
#else
#   if defined(LORAD_VARIABLE_TOPOLOGY)
        // Record tree length and each edge proportion in param_vect and compute the log of the Jacobian
        double TL = 0.0;
        std::vector< std::pair<Split, double> > split_edgeprop_vect;
        //unsigned i = 0;
        for (auto nd : _tree->_preorder) {
            double v = nd->_edge_length;
            assert(v > 0.0);
            split_edgeprop_vect.push_back(std::make_pair(nd->getSplit(), v));
            TL += v;
        }
        
        // Sort split_edgeprop_vect by split (since split occurs first in each pair)
        std::sort(split_edgeprop_vect.begin(), split_edgeprop_vect.end());
        
        // Convert edge lengths into proportions by dividing each by TL
        for (auto & p : split_edgeprop_vect) {
            p.second /= TL;
        }
        
        // Record everything in param_vect and compute the log of the Jacobian
        double log_TL = log(TL);
        double log_jacobian = log_TL;
        param_vect.push_back(log_TL);
        
        // Suppose there are 5 edge length proportions: p1, p2, p3, p4, p5
        // These are stored as 4 values in param_vect:
        //   param_vect[0] = log(p2) - log(p1)
        //   param_vect[1] = log(p3) - log(p1)
        //   param_vect[2] = log(p4) - log(p1)
        //   param_vect[3] = log(p5) - log(p1)
        // To later retrieve the original 5 proportions from param_vect:
        //   phi = (p1/p1) + (p2/p1) + (p3/p1) + (p4/p1) + (p5/p1) = (p1 + p2 + p3 + p4 + p5)/p1 = 1/p1
        //     = 1.0 + exp(param_vect[0]) + exp(param_vect[1]) + exp(param_vect[2]) + exp(param_vect[3])
        //   p1 = 1.0/phi
        //   p2 = exp(param_vect[0])/phi
        //   p3 = exp(param_vect[1])/phi
        //   p4 = exp(param_vect[2])/phi
        //   p5 = exp(param_vect[3])/phi
        auto & first_pair = split_edgeprop_vect[0];
        double logprop_first = log(first_pair.second);
        log_jacobian += logprop_first;
        for (unsigned i = 1; i < split_edgeprop_vect.size(); ++i) {
            auto & p = split_edgeprop_vect[i];
            double logprop = log(p.second);
            double transformed = logprop - logprop_first;
            log_jacobian += logprop;
            param_vect.push_back(transformed);
        }
#   else
        // Record tree length and each edge proportion in param_vect and compute the log of the Jacobian
        double TL = 0.0;
        std::vector<double> edge_length_proportions(_tree->_preorder.size());
        unsigned i = 0;
        for (auto nd : _tree->_preorder) {
            assert(nd->_edge_length > 0.0);
            edge_length_proportions[i++] = nd->_edge_length;
            TL += nd->_edge_length;
        }
        
        // Convert edge lengths into proportions by dividing each by TL
        std::transform(edge_length_proportions.begin(), edge_length_proportions.end(), edge_length_proportions.begin(), [TL](double v) {return v/TL;});
        
        // Record everything in param_vect and compute the log of the Jacobian
        double log_TL = log(TL);
        double log_jacobian = log_TL;
        param_vect.push_back(log_TL);
        
        double first = edge_length_proportions[0];
        double log_first = log(first);
        log_jacobian += log_first;
        for (unsigned i = 1; i < edge_length_proportions.size(); ++i) {
            double p = edge_length_proportions[i];
            double logp = log(p);
            double transformed = logp - log_first;
            log_jacobian += logp;
            param_vect.push_back(transformed);
        }
#   endif
#endif

        return log_jacobian;
    }

    inline double TreeManip::setEdgeLengthsFromLogTransformed(Eigen::VectorXd & param_vect, double TL, unsigned start_at, unsigned nedges) {
        double log_jacobian = 0.0;
#if defined(HOLDER_ETAL_PRIOR)
        // In this case, edge length parameters are saved in param_vect:
        //      TL       = 0.0
        //      start_at = 0
        //      nedges   = total number of edges
#   if defined(LORAD_VARIABLE_TOPOLOGY)
        // Edge lengths stored in standardized order (splits sorted), so first need to assign edge lengths
        // from param_vect to splits in a map, then walk through tree and set edge lengths using the map
        
        // Store splits in a set, which is conveniently kept in sorted order
        std::set<Split> splitset;
        storeSplits(splitset, /*include_trivial_splits*/true);
        assert(splitset.size() == nedges);
                
        // Create map assigning edge lengths (values) to splits (keys)
        std::map<Split, double> split_map;
        unsigned i = 0;
        for (Split s : splitset) {
            double logv = param_vect[start_at + i];
            log_jacobian += logv;
            double v = exp(logv);
            split_map[s] = v;
            ++i;
        }

        // Assign edge lengths according to split
        for (auto nd : _tree->_preorder) {
            Split s = nd->getSplit();
            double v = split_map[s];
            nd->setEdgeLength(v);
        }
#   else
        for (unsigned i = 0; i < nedges; ++i) {
            double logv = param_vect[start_at + i];
            double v = exp(logv);
            Node * nd = _tree->_preorder[i];
            nd->setEdgeLength(v);
            log_jacobian += logv;
        }
#   endif
#else   // !defined(HOLDER_ETAL_PRIOR)
        // In this case, TL and edge length proportions are saved in param_vect:
        //      TL       = tree length
        //      start_at = 1
        //      nedges   = total number of edges minus 1
#   if defined(LORAD_VARIABLE_TOPOLOGY)
        // Edge length proportions stored in standardized order (splits sorted),
        // so first need to assign edge length proportions from param_vect to
        // splits in a map, then walk through tree and set edge length
        // proportions using the map
        
        // Store splits in a set, which is conveniently kept in sorted order
        std::set<Split> splitset;
        storeSplits(splitset, /*include_trivial_splits*/true);
        assert(splitset.size() == nedges + 1);
        
        // Suppose there are 5 edge length proportions: p1, p2, p3, p4, p5
        // These were originally stored as 4 values in param_vect:
        //   param_vect[0] = log(p2) - log(p1)
        //   param_vect[1] = log(p3) - log(p1)
        //   param_vect[2] = log(p4) - log(p1)
        //   param_vect[3] = log(p5) - log(p1)
        // To retrieve the original 5 proportions from param_vect:
        //   phi = (p1/p1) + (p2/p1) + (p3/p1) + (p4/p1) + (p5/p1) = (p1 + p2 + p3 + p4 + p5)/p1 = 1/p1
        //     = 1.0 + exp(param_vect[0]) + exp(param_vect[1]) + exp(param_vect[2]) + exp(param_vect[3])
        //   p1 = 1.0/phi
        //   p2 = exp(param_vect[0])/phi
        //   p3 = exp(param_vect[1])/phi
        //   p4 = exp(param_vect[2])/phi
        //   p5 = exp(param_vect[3])/phi
        double phi = 1.0;
        for (unsigned i = 0; i < nedges; ++i) {
            double logp = param_vect[start_at + i];
            phi += exp(logp);
        }

        // Create map assigning edge lengths (values) to splits (keys)
        std::map<Split, double> split_map;
        auto split_iter = splitset.begin();
        Split s = *split_iter;
        double proportion = 1.0/phi;
        double edge_length = TL*proportion;
        split_map[s] = edge_length;
        log_jacobian = log(TL);
        unsigned i = 0;
        for (split_iter++; split_iter != splitset.end(); split_iter++) {
            Split s = *split_iter;
            double logr = param_vect[start_at + i];
            proportion = exp(logr)/phi;
            log_jacobian += log(proportion);
            edge_length = TL*proportion;
            split_map[s] = edge_length;
            ++i;
        }

        // Assign edge lengths according to split
        for (auto nd : _tree->_preorder) {
            Split s = nd->getSplit();
            double v = split_map[s];
            nd->setEdgeLength(v);
        }
#   else
        double phi = 1.0;
        Node * nd = _tree->_preorder[0];
        nd->setEdgeLength(1.0);
        for (unsigned i = 0; i < nedges; ++i) {
            double t = param_vect[start_at + i];
            double r = exp(t);  // r is ratio of this edge proportion to first edge proportion
            nd = _tree->_preorder[i + 1];
            nd->setEdgeLength(r);
            phi += r;
        }
        double first = 1.0/phi;
        log_jacobian = log(TL);
        for (auto nd : _tree->_preorder) {
            double proportion = first*nd->getEdgeLength(); // edge length currently equal to ratio of this edge proportion to first
            double edgelen = TL*proportion; // convert edge proportion to edge length
            nd->setEdgeLength(edgelen);
            log_jacobian += log(proportion); //TODO: I think first edge proportion should be skipped here
        }
#   endif
#endif
        return log_jacobian;
    }

#if defined(POLGSS)
    inline void TreeManip::sampleTree() {
        std::vector<double> tmp;
#if defined(HOLDER_ETAL_PRIOR)
        copyEdgeLengthsTo(tmp);
        for (auto v : tmp) {
            _sampled_edge_lengths.push_back(v);
        }
#else
        double TL = copyEdgeProportionsTo(tmp);
        _sampled_edge_proportions.push_back(tmp);
        _sampled_tree_lengths.push_back(TL);
#endif
    }
    
    inline std::string TreeManip::calcExpRefDist(std::string title, std::vector<double> & vect, std::vector<double> & params) {
        // Compute sums and sums-of-squares for each component
        unsigned n = (unsigned)vect.size();
        double sumv   = 0.0;
        double sumsqv = 0.0;
        for (unsigned i = 0; i < n; i++) {
            double v = vect[i];
            sumv += v;
            sumsqv += v*v;
        }
        
        // Compute mean and variance
        double mu;
        //double s;
        mu = sumv/n;
        //s = (sumsqv - mu*mu*n)/(n-1);
        
        // Compute parameters of reference distribution and save each
        // as an element of the string vector svect
        double rate = 1.0/mu;
        params.resize(1);
        params[0] = rate;
        std::string refdiststr = boost::str(boost::format("%s = default:%.3f\n") % title % rate);
        
        return refdiststr;
    }
    
    inline std::string TreeManip::calcGammaRefDist(std::string title, std::vector<double> & vect, std::vector<double> & params) {
        //TODO: nearly identical to Model::calcGammaRefDist - make one version that can be used by both Model and TreeManip
        // Compute sums and sums-of-squares for each component
        unsigned n = (unsigned)vect.size();
        double sumv   = 0.0;
        double sumsqv = 0.0;
        for (unsigned i = 0; i < n; i++) {
            double v = vect[i];
            sumv += v;
            sumsqv += v*v;
        }
        
        // Compute mean and variance
        double mu;
        double s;
        mu = sumv/n;
        s = (sumsqv - mu*mu*n)/(n-1);
        
        // Compute parameters of reference distribution and save each
        // as an element of the string vector svect
        // s = shape*scale^2
        // mu = shape*scale
        double scale = s/mu;
        double shape = mu/scale;
        params.resize(2);
        params[0] = shape;
        params[1] = scale;
        std::string refdiststr = boost::str(boost::format("%s = default:%.3f, %.3f\n") % title % shape % scale);
        
        return refdiststr;
    }
    
    inline std::string TreeManip::calcDirichletRefDist(std::string title, std::vector< std::vector<double> > & vect, std::vector<double> & params) {
        //TODO: identical to Model::calcGammaRefDist - make one version that can be used by both Model and TreeManip
        // Ming-Hui Chen method of matching component variances
        // mu_i = phi_i/phi is mean of component i (estimate using sample mean)
        // s_i^2 is sample variance of component i
        //
        //       sum_i mu_i^2 (1 - mu_i)^2
        // phi = --------------------------- - 1
        //       sum_i s_i^2 mu_i (1 - mu_i)
        //
        // phi_i = phi mu_i
        unsigned n = (unsigned)vect.size();
        assert(n > 0);
        unsigned k = (unsigned)vect[0].size();
        
        // Compute sums and sums-of-squares for each component
        std::vector<double> sums(k, 0.0);
        std::vector<double> sumsq(k, 0.0);
        for (unsigned i = 0; i < n; i++) {
            std::vector<double> & dir = vect[i];
            for (unsigned j = 0; j < k; j++) {
                double v = dir[j];
                sums[j] += v;
                sumsq[j] += v*v;
            }
        }
        
        // Compute means and variances for each component
        std::vector<double> mu(k, 0.0);
        std::vector<double> s(k, 0.0);
        double numer = 0.0;
        double denom = 0.0;
        for (unsigned j = 0; j < k; j++) {
            mu[j] = sums[j]/n;
            numer += mu[j]*mu[j]*(1.0 - mu[j])*(1.0 - mu[j]);
            s[j] = (sumsq[j] - mu[j]*mu[j]*n)/(n-1);
            denom += s[j]*mu[j]*(1.0 - mu[j]);
        }
        
        // Compute phi
        double phi = numer/denom - 1.0;

        // Compute parameters of reference distribution and save each
        // as an element of the string vector svect
        std::vector<std::string> svect;
        params.clear();
        for (unsigned j = 0; j < k; j++) {
            double c = phi*mu[j];
            params.push_back(c);
            std::string stmp = boost::str(boost::format("%.3f") % c);
            svect.push_back(stmp);
        }
        std::string refdiststr = boost::str(boost::format("%s = default:%s\n") % title % boost::algorithm::join(svect, ","));
        
        return refdiststr;
    }

    inline std::string TreeManip::saveReferenceDistributions() {
        std::map<std::string, std::vector<double> > dummy_refdist_map;
        return calcReferenceDistributions(dummy_refdist_map);
    }
    
    inline std::string TreeManip::calcReferenceDistributions(std::map<std::string, std::vector<double> > & refdist_map) {
        // Calculate and save reference distribution parameters in a conf file that can be used
        // in a subsequent generalized steppingstone analysis
#if defined(HOLDER_ETAL_PRIOR)
        std::vector<double> & v = refdist_map["Edge Length"];
        std::string s = calcExpRefDist("edgelenrefdist", _sampled_edge_lengths, v);
#else
        std::vector<double> & v = refdist_map["Edge Proportions"];
        std::string s = calcDirichletRefDist("edgeproprefdist", _sampled_edge_proportions, v);
        std::vector<double> & vv = refdist_map["Tree Length"];
        s += calcGammaRefDist("treelenrefdist", _sampled_tree_lengths, vv);
#endif
        
        return s;
    }
#endif

#if defined(POLGHM)
    inline void TreeManip::setModelToSampledPoint(unsigned i) {
#if defined(HOLDER_ETAL_PRIOR)
        assert(_sampled_edge_lengths.size() > i);
        copyEdgeLengthsFrom(_sampled_edge_lengths[i]);
#else
        assert(_sampled_edge_proportions.size() > i);
        assert(_sampled_tree_lengths.size() > i);
        double TL = _sampled_tree_lengths[i];
        copyEdgeProportionsFrom(TL, _sampled_edge_proportions[i]);
#endif
        
    }
#endif

}
