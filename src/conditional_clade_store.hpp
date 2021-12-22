#pragma once

#include <map>
#include <set>
#include "split.hpp"

namespace lorad {

    struct CondCladeInfo {
        CondCladeInfo() : _count(0), _empirical_prob(0.0), _reference_prob(0.0) {}
        
        unsigned _count;            // number of times this clade was seen
        double   _empirical_prob;   // count divided by total count for all clades with same parent
        double   _reference_prob;   // probability used for reference distribution (non-zero even if count==0)
    };
    
    struct ParentMapValue {
        ParentMapValue() : _unseen_reference_prob(0.0) {}
        double _unseen_reference_prob;
        std::map< Split, CondCladeInfo> _child_map;
    };
    
    class ConditionalCladeStore {

        public:
        
                     ConditionalCladeStore();
                     ~ConditionalCladeStore();
                                                            
            void     addParentChildSplit(Split & parent, Split & child);
            unsigned getCount(Split & parent, Split & child);
            double   getEmpiricalProb(Split & parent, Split & child);
            double   getReferenceProb(Split & parent, Split & child);
            void     finalize(double unseen_fraction = 0.1);
            void     summarize();
            
            typedef  std::shared_ptr<ConditionalCladeStore>  SharedPtr;

        private:

            //double calcLogNRootedTopologies(unsigned n) const;
            double calcTotalSplits(unsigned parent_clade_size) const;

            std::map< Split, ParentMapValue > _parent_map;
    };
        
    inline ConditionalCladeStore::ConditionalCladeStore() {
    }
        
    inline ConditionalCladeStore::~ConditionalCladeStore() {
    }
    
    inline double ConditionalCladeStore::calcTotalSplits(unsigned parent_clade_size) const {
        // Example for parent_clade_size = 3
        // half = 1, remainder = 1, start with 2
        // Total number of splits is
        //
        // 000  1 = {3 choose 0} <-- omit (zero bits set)
        //
        // 001
        // 010  3 = {3 choose 1} <-- omit (fewer than half of bits set)
        // 100
        //
        // 011
        // 101  3 = {3 choose 2} <-- keep all
        // 110
        //
        // 111  1 = {3 choose 3} <-- omit (all bits set)
        
        // Example for parent_clade_size = 4
        // half = 2, remainder = 0, start with 2
        // Total number of splits is 6/2 + 4 = 7
        //
        // 0000  1 = {4 choose 0} <-- omit (zero bits set)
        //
        // 0001
        // 0010  4 = {4 choose 1} <-- omit (fewer than half of bits set)
        // 0100
        // 1000
        //
        // 0011
        // 0101
        // 1001 6 = {4 choose 2} <-- half of bits set, keep half of these
        // 0110
        // 1010
        // 1100
        //
        // 0111
        // 1011  4 = {4 choose 3} <-- keep all
        // 1101
        // 1110
        //
        // 1111  1 = {4 choose 4} <-- omit (all bits set)
        
        unsigned half      = parent_clade_size / 2;
        unsigned remainder = parent_clade_size % 2;
        double nsplits = 0;
        if (remainder == 0) {
            // deal with half
            unsigned k = half;
            double nchoosek = exp(lgamma(parent_clade_size + 1) - lgamma(k + 1) - lgamma(parent_clade_size - k + 1));
            nsplits += nchoosek/2.;
            
            // deal with rest
            for (k = half + 1; k <= parent_clade_size - 1; k++) {
                double nchoosek = exp(lgamma(parent_clade_size + 1) - lgamma(k + 1) - lgamma(parent_clade_size - k + 1));
                nsplits += nchoosek;
            }
        }
        else {
            for (unsigned k = half + remainder; k <= parent_clade_size - 1; k++) {
                // add parent_clade_size choose k
                double nchoosek = exp(lgamma(parent_clade_size + 1) - lgamma(k + 1) - lgamma(parent_clade_size - k + 1));
                nsplits += nchoosek;
            }
        }
        return nsplits;
    }
    
    inline void ConditionalCladeStore::finalize(double unseen_fraction) {
        // compute empirical and reference probabilities from counts
        for (auto & pm : _parent_map) {
            auto & parent_map_value = pm.second;
            
            // loop through all child splits, totalling up counts
            unsigned total_count = 0;
            for (auto & clade : parent_map_value._child_map) {
                CondCladeInfo & c = clade.second;
                total_count += c._count;
            }
            
            // loop through child splits again, computing probabilities
            for (auto & clade : parent_map_value._child_map) {
                CondCladeInfo & c = clade.second;
                c._empirical_prob = (double)c._count/total_count;
                c._reference_prob = c._empirical_prob*(1.0 - unseen_fraction);
            }
            
            // set reference probability that will be used if child split lookup fails
            unsigned parent_clade_size = pm.first.countBitsSet();
            double total_splits = calcTotalSplits(parent_clade_size);
            double observed_splits = (double)parent_map_value._child_map.size();
            double unobserved_splits = total_splits - observed_splits;
            parent_map_value._unseen_reference_prob = unseen_fraction/unobserved_splits;
        }
    }
    
    inline unsigned ConditionalCladeStore::getCount(Split & parent, Split & child) {
        unsigned count = 0;
        auto parent_map_iter = _parent_map.find(parent);
        if (parent_map_iter != _parent_map.end()) {
            ParentMapValue & v = parent_map_iter->second;
            auto child_map_iter = v._child_map.find(child);
            if (child_map_iter != v._child_map.end()) {
                count = child_map_iter->second._count;
            }
        }
        return count;
    }
        
    inline double ConditionalCladeStore::getEmpiricalProb(Split & parent, Split & child) {
        double empirical_prob = 0.0;
        auto parent_map_iter = _parent_map.find(parent);
        if (parent_map_iter != _parent_map.end()) {
            ParentMapValue & v = parent_map_iter->second;
            auto cond_clade_iter = v._child_map.find(child);
            if (cond_clade_iter != v._child_map.end()) {
                empirical_prob = cond_clade_iter->second._empirical_prob;
            }
        }
        return empirical_prob;
    }
    
//    inline double ConditionalCladeStore::calcLogNRootedTopologies(unsigned n) const {
//        //  n   unrooted
//        //  4     1*3
//        //  5     1*3*5
//        //  6     1*3*5*7
//        //  n     1*3*5*7*...*(2n-5)
//        //
//        //       1*2*3*4*5*6*7       (2n-5)!         7!
//        // n=6  --------------- = -------------- = ------
//        //      (2*1)(2*2)(2*3)   2^(n-3) (n-3)!   2^3 3!
//        //
//        //        1*2*3*4*5*6*7*8*9       (2n-5)!         9!
//        // n=7  -------------------- = -------------- = ------
//        //      (2*1)(2*2)(2*3)(2*4)   2^(n-3) (n-3)!   2^4 4!
//
//        unsigned nplus1 = n + 1;
//        double logntopol = lgamma(2*nplus1 - 5 + 1);
//        logntopol -= lgamma(nplus1 - 3 + 1);
//        logntopol -= (double)(nplus1 - 3)*log(2);
//        return logntopol;
//    }
        
    inline double ConditionalCladeStore::getReferenceProb(Split & parent, Split & child) {
        double reference_prob = 0.0;
        auto parent_map_iter = _parent_map.find(parent);
        if (parent_map_iter != _parent_map.end()) {
            ParentMapValue & v = parent_map_iter->second;
            auto cond_clade_iter = v._child_map.find(child);
            if (cond_clade_iter != v._child_map.end()) {
                reference_prob = cond_clade_iter->second._reference_prob;
            }
            else {
                reference_prob = v._unseen_reference_prob;
            }
        }
        else {
            double total_splits = calcTotalSplits(parent.countBitsSet());
            return 1.0/total_splits;
        }
        return reference_prob;
    }
        
    inline void ConditionalCladeStore::addParentChildSplit(Split & parent, Split & child) {
        auto & v = _parent_map[parent];
        CondCladeInfo & c = v._child_map[child];
        c._count++;
    }
        
    inline void ConditionalCladeStore::summarize() {
        std::cout << "\nConditional Clade Store:" << std::endl;
        for (auto & pm : _parent_map) {
            auto & parent_split = pm.first;
            std::cout << boost::format("+------- %s\n") % parent_split.createPatternRepresentation();
            auto & v = pm.second;
            for (auto & clade : v._child_map) {
                auto & child_split = clade.first;
                CondCladeInfo & c = clade.second;
                std::cout << boost::format("| %6d %s\n") % c._count % child_split.createPatternRepresentation();
            }
        }
    }
        
}
