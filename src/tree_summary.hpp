#pragma once

#include <set>
#include <map>
#include <vector>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <boost/format.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include "split.hpp"
#include "tree_manip.hpp"
#include "xlorad.hpp"

#include "ncl/nxsmultiformat.h"

#include "output_manager.hpp"
extern lorad::OutputManager om;

namespace lorad {

    class TreeSummary {
        public:
                                        TreeSummary();
                                        ~TreeSummary();

            void                        readTreefile(const std::string filename, unsigned skip);
            void                        showSummary() const;
            typename Tree::SharedPtr    getTree(unsigned index);
            std::string                 getNewick(unsigned index);
            void                        clear();
            
            
            void                        showResClassSummary() const;

        private:

            Split::treemap_t            _treeIDs;
            std::vector<std::string>    _newicks;

        public:

            typedef std::shared_ptr< TreeSummary > SharedPtr;
    };

    // insert member function bodies here

    inline TreeSummary::TreeSummary() {
    }

    inline TreeSummary::~TreeSummary() {
    }

    inline Tree::SharedPtr TreeSummary::getTree(unsigned index) {
        if (index >= _newicks.size())
            throw XLorad("getTree called with index >= number of stored trees");

        TreeManip tm;

        // build the tree
        tm.buildFromNewick(_newicks[index], false, true);

        return tm.getTree();
    }

    inline std::string TreeSummary::getNewick(unsigned index) {
        if (index >= _newicks.size())
            throw XLorad("getNewick called with index >= number of stored trees");

        return _newicks[index];
    }

    inline void TreeSummary::clear() {
        _treeIDs.clear();
        _newicks.clear();
    }

    inline void TreeSummary::readTreefile(const std::string filename, unsigned skip) {
        TreeManip tm;
        Split::treeid_t splitset;

        // See http://phylo.bio.ku.edu/ncldocs/v2.1/funcdocs/index.html for NCL documentation

        MultiFormatReader nexusReader(-1, NxsReader::WARNINGS_TO_STDERR);
        try {
            nexusReader.ReadFilepath(filename.c_str(), MultiFormatReader::NEXUS_FORMAT);
        }
        catch(...) {
            nexusReader.DeleteBlocksFromFactories();
            throw;
        }

        int numTaxaBlocks = nexusReader.GetNumTaxaBlocks();
        for (int i = 0; i < numTaxaBlocks; ++i) {
            clear();
            NxsTaxaBlock * taxaBlock = nexusReader.GetTaxaBlock(i);
            std::string taxaBlockTitle = taxaBlock->GetTitle();

            const unsigned nTreesBlocks = nexusReader.GetNumTreesBlocks(taxaBlock);
            for (unsigned j = 0; j < nTreesBlocks; ++j) {
                const NxsTreesBlock * treesBlock = nexusReader.GetTreesBlock(taxaBlock, j);
                unsigned ntrees = treesBlock->GetNumTrees();
                if (skip < ntrees) {
                    for (unsigned t = skip; t < ntrees; ++t) {
                        const NxsFullTreeDescription & d = treesBlock->GetFullTreeDescription(t);

                        // store the newick tree description
                        std::string newick = d.GetNewick();
                        _newicks.push_back(newick);
                        unsigned tree_index = (unsigned)_newicks.size() - 1;

                        // build the tree
                        tm.buildFromNewick(newick, false, true);

                        // store set of splits
                        splitset.clear();
                        tm.storeSplits(splitset);

                        // iterator iter will point to the value corresponding to key splitset
                        // or to end (if splitset is not already a key in the map)
                        Split::treemap_t::iterator iter = _treeIDs.lower_bound(splitset);

                        if (iter == _treeIDs.end() || iter->first != splitset) {
                            // splitset key not found in map, need to create an entry
                            std::vector<unsigned> v(1, tree_index);  // vector of length 1 with only element set to tree_index
                            _treeIDs.insert(iter, Split::treemap_t::value_type(splitset, v));
                        }
                        else {
                            // splitset key was found in map, need to add this tree's index to vector
                            iter->second.push_back(tree_index);
                        }
                    } // trees loop
                } // if skip < ntrees
            } // TREES block loop
        } // TAXA block loop

        // No longer any need to store raw data from nexus file
        nexusReader.DeleteBlocksFromFactories();
    }

    inline void TreeSummary::showSummary() const {
        // Produce some output to show that it works
        ::om.outputConsole(boost::format("\nRead %d trees from file\n") % _newicks.size());

        // Show all unique topologies with a list of the trees that have that topology
        // Also create a map that can be used to sort topologies by their sample frequency
        typedef std::pair<unsigned, unsigned> sorted_pair_t;
        std::vector< sorted_pair_t > sorted;
        int t = 0;
        for (auto & key_value_pair : _treeIDs) {
            unsigned topology = ++t;
            unsigned ntrees = (unsigned)key_value_pair.second.size();
            sorted.push_back(std::pair<unsigned, unsigned>(ntrees,topology));
            ::om.outputConsole(boost::format("Topology %d seen in these %d trees:\n  ") % topology % ntrees);
            std::vector<std::string> tree_indices(key_value_pair.second.size());
            std::transform(key_value_pair.second.begin(), key_value_pair.second.end(), tree_indices.begin(), [](unsigned i){return std::to_string(i);});
            ::om.outputConsole(boost::algorithm::join(tree_indices, ","));
            ::om.outputConsole();
        }

        // Show sorted histogram data
        std::sort(sorted.begin(), sorted.end());
        //unsigned npairs = (unsigned)sorted.size();
        ::om.outputConsole("\nTopologies sorted by sample frequency:\n");
        ::om.outputConsole(boost::format("%20s %20s\n") % "topology" % "frequency");
        for (auto & ntrees_topol_pair : boost::adaptors::reverse(sorted)) {
            unsigned n = ntrees_topol_pair.first;
            unsigned t = ntrees_topol_pair.second;
            ::om.outputConsole(boost::format("%20d %20d\n") % t % n);
        }
    }

    inline void TreeSummary::showResClassSummary() const {
        std::set<unsigned> mset;
        std::map<unsigned, std::vector<unsigned> > mmap;
        for (auto & key_value_pair : _treeIDs) {
            Split::treeid_t splitset = key_value_pair.first;
            unsigned m = (unsigned)splitset.size();
            mset.insert(m);
            unsigned ntrees = (unsigned)key_value_pair.second.size();
            mmap[m].push_back(ntrees);
        }
        
        for (auto m : mset) {
            std::vector<unsigned> & v = mmap[m];
            unsigned msum = std::accumulate(v.begin(), v.end(), 0);
            std::vector<std::string> vstr(v.size());
            std::transform(v.begin(), v.end(), vstr.begin(), [](unsigned i){return std::to_string(i);});
            ::om.outputConsole(boost::format("m = %d (%d): barplot(c(") % m % msum);
            ::om.outputConsole(boost::algorithm::join(vstr, ","));
            ::om.outputConsole("))\n");
        }
    }

}
