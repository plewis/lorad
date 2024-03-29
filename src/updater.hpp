#pragma once

#include "conditionals.hpp"

#include "tree.hpp"
#include "tree_manip.hpp"
#include "lot.hpp"
#include "xlorad.hpp"
#include "likelihood.hpp"
#include "topo_prior_calculator.hpp"
#include "conditional_clade_store.hpp"

namespace lorad {
    class Chain;

    class Updater { 
    
        friend class Chain;

        public:
            typedef std::shared_ptr<Updater>        SharedPtr;
        
            TreeManip::SharedPtr                    getTreeManip() const;

                                                    Updater();
            virtual                                 ~Updater();

            void                                    setLikelihood(Likelihood::SharedPtr likelihood);
            void                                    setTreeManip(TreeManip::SharedPtr treemanip);
            void                                    setLot(Lot::SharedPtr lot);
            void                                    setLambda(double lambda);
            void                                    setHeatingPower(double p);
            void                                    setSteppingstoneMode(unsigned mode);
            void                                    setTuning(bool on);
            void                                    setTargetAcceptanceRate(double target);
            void                                    setPriorParameters(const std::vector<double> & c);
            void                                    setConditionalCladeStore(ConditionalCladeStore::SharedPtr ccs);
            void                                    setRefDistParameters(const std::vector<double> & c);
            void                                    setTopologyPriorOptions(bool resclass, double C);
            void                                    setWeight(double w);
            void                                    calcProb(double wsum);

            double                                  getLambda() const;
            double                                  getWeight() const;
            double                                  getProb() const;
            double                                  getAcceptPct() const;
            double                                  getNumUpdates() const;
            std::string                             getUpdaterName() const;

            virtual void                            clear();

            virtual void                            debugPriorCalculation();
            virtual double                          calcLogPrior() = 0;
            double                                  calcLogTopologyPrior() const;
#if defined(HOLDER_ETAL_PRIOR)
            double                                  calcLogEdgeLengthPrior() const;
#else
            std::pair<double,double>                calcLogEdgeLengthPrior() const;
#endif
            //double                                  calcLogEdgeLengthRefDist() const;
            virtual double                          calcLogRefDist() = 0;
            double                                  calcLogLikelihood() const;
            virtual double                          update(double prev_lnL);

            static double                           getLogZero();
            
        protected:

            virtual void                            reset();
            virtual void                            tune(bool accepted);

            virtual void                            revert() = 0;
            virtual void                            proposeNewState() = 0;

            Lot::SharedPtr                          _lot;
            Likelihood::SharedPtr                   _likelihood;
            TreeManip::SharedPtr                    _tree_manipulator;
            std::string                             _name;
            double                                  _weight;
            double                                  _prob;
            double                                  _lambda;
            double                                  _log_hastings_ratio;
            double                                  _log_jacobian;
            double                                  _target_acceptance;
            unsigned                                _naccepts;
            unsigned                                _nattempts;
            bool                                    _tuning;
            std::vector<double>                     _prior_parameters;
            ConditionalCladeStore::SharedPtr        _conditional_clade_store;
            std::vector<double>                     _refdist_parameters;
            unsigned                                _ss_mode;
            double                                  _heating_power;
            mutable PolytomyTopoPriorCalculator     _topo_prior_calculator;
            
            static const double                     _log_zero;
    }; 
 
    inline Updater::Updater() {
        clear();
    } 

    inline Updater::~Updater() {
    }

    inline void Updater::clear() { 
        _name                   = "updater";
        _tuning                 = true;
        _lambda                 = 0.0001;
        _weight                 = 1.0;
        _prob                   = 0.0;
        _target_acceptance      = 0.3;
        _naccepts               = 0;
        _nattempts              = 0;
        _heating_power          = 1.0;
        _prior_parameters.clear();
        _refdist_parameters.clear();
        _ss_mode                = 0;    // no steppingstone
        reset();
    } 

    inline void Updater::reset() { 
        _log_hastings_ratio = 0.0;
        _log_jacobian       = 0.0;
    } 

    inline void Updater::setLikelihood(Likelihood::SharedPtr likelihood) { 
        _likelihood = likelihood;
    } 
    
    inline void Updater::setTreeManip(TreeManip::SharedPtr treemanip) { 
        _tree_manipulator = treemanip;
    } 

    inline TreeManip::SharedPtr Updater::getTreeManip() const { 
        return _tree_manipulator;
    } 

    inline void Updater::setLot(Lot::SharedPtr lot) { 
        _lot = lot;
    } 
    
    inline void Updater::setHeatingPower(double p) { 
        _heating_power = p;
    } 

    inline void Updater::setSteppingstoneMode(unsigned mode) {
        // Steppingstone mode:
        //   0: no steppingstone
        //   1: steppingstone (Xie et al. 2011)
        //   2: generalized steppingstone (Fan et al. 2011)
        _ss_mode = mode;
    }

    inline void Updater::setLambda(double lambda) {
        _lambda = lambda;
    } 

    void Updater::setTuning(bool do_tune) { 
        _tuning = do_tune;
        _naccepts = 0;
        _nattempts = 0;
    } 

    inline void Updater::tune(bool accepted) { 
        _nattempts++;
        if (_tuning) {
            double gamma_n = 10.0/(100.0 + (double)_nattempts);
            if (accepted)
                _lambda *= 1.0 + gamma_n*(1.0 - _target_acceptance)/(2.0*_target_acceptance);
            else
                _lambda *= 1.0 - gamma_n*0.5;

            // Prevent run-away increases in boldness for low-information marginal densities
            if (_lambda > 1000.0)
                _lambda = 1000.0;
        }
    } 

    inline void Updater::setTargetAcceptanceRate(double target) { 
        _target_acceptance = target;
    } 

    inline void Updater::setPriorParameters(const std::vector<double> & c) { 
        _prior_parameters.clear();
        _prior_parameters.assign(c.begin(), c.end());
    } 
    
    inline void Updater::setConditionalCladeStore(ConditionalCladeStore::SharedPtr ccs) {
        _conditional_clade_store = ccs;
    }
    
    inline void Updater::setRefDistParameters(const std::vector<double> & c) {
        _refdist_parameters.clear();
        _refdist_parameters.assign(c.begin(), c.end());
    }
    
    inline void Updater::setWeight(double w) {
        _weight = w;
    } 
    
    inline void Updater::calcProb(double wsum) { 
        assert(wsum > 0.0);
        _prob = _weight/wsum;
    } 

    inline double Updater::getLambda() const { 
        return _lambda;
    } 
    
    inline double Updater::getWeight() const { 
        return _weight;
    } 
    
    inline double Updater::getProb() const { 
        return _prob;
    } 
    
    inline double Updater::getAcceptPct() const { 
        return (_nattempts == 0 ? 0.0 : (100.0*_naccepts/_nattempts));
    } 

    inline double Updater::getNumUpdates() const { 
        return _nattempts;
    } 

    inline std::string Updater::getUpdaterName() const { 
        return _name;
    } 

    inline double Updater::calcLogLikelihood() const { 
        return _likelihood->calcLogLikelihood(_tree_manipulator->getTree());
    } 

    inline double Updater::update(double prev_lnL) { 
        double prev_log_prior = calcLogPrior();
        double prev_log_refdist = 0.0;
        if (_ss_mode == 2) {
            // Steppingstone mode:
            //   0: no steppingstone
            //   1: steppingstone (Xie et al. 2011)
            //   2: generalized steppingstone (Fan et al. 2011)
            prev_log_refdist = calcLogRefDist();
        }
        
        // Clear any nodes previously selected so that we can detect those nodes
        // whose partials and/or transition probabilities need to be recalculated
        _tree_manipulator->deselectAllPartials();
        _tree_manipulator->deselectAllTMatrices();
        
        // Set model to proposed state and calculate _log_hastings_ratio
        proposeNewState();
        
        // Use alternative partials and transition probability buffer for any selected nodes
        // This allows us to easily revert to the previous values if the move is rejected
        _tree_manipulator->flipPartialsAndTMatrices();

        // Calculate the log-likelihood and log-prior for the proposed state
        double log_likelihood = calcLogLikelihood();
        double log_prior = calcLogPrior();
        
        // Decide whether to accept or reject the proposed state
        bool accept = true;
        if (log_prior > _log_zero) {
            double log_R = 0.0;
            if (_ss_mode == 1) {
                // Xie et al. 2011 steppingstone
                log_R += _heating_power*(log_likelihood - prev_lnL);
                log_R += (log_prior - prev_log_prior);
            }
            else if (_ss_mode == 2) {
                // Fan et al. 2011 generalized steppingstone
                double log_refdist = calcLogRefDist();
                log_R += _heating_power*(log_likelihood - prev_lnL);
                log_R += _heating_power*(log_prior - prev_log_prior);
                log_R += (1.0 - _heating_power)*(log_refdist - prev_log_refdist);
            }
            else {
                // normal heated chain
                assert(_ss_mode == 0);
                log_R += _heating_power*(log_likelihood - prev_lnL);
                log_R += _heating_power*(log_prior - prev_log_prior);
            }
            log_R += _log_hastings_ratio;
            log_R += _log_jacobian;

            double logu = _lot->logUniform();
            if (logu > log_R)
                accept = false;
        }
        else
            accept = false;

        if (accept) {
            _naccepts++;
        }
        else {
            revert();
            _tree_manipulator->flipPartialsAndTMatrices();
            log_likelihood = prev_lnL;
        }

        tune(accept);
        reset();

        return log_likelihood;
    } 
    
    inline void Updater::setTopologyPriorOptions(bool resclass, double C) {
        _topo_prior_calculator.setC(C);
        if (resclass)
            _topo_prior_calculator.chooseResolutionClassPrior();
        else
            _topo_prior_calculator.choosePolytomyPrior();
    }   
    
    inline double Updater::calcLogTopologyPrior() const {
        Tree::SharedPtr tree = _tree_manipulator->getTree();
        assert(tree);
        if (tree->isRooted())
            _topo_prior_calculator.chooseRooted();
        else
            _topo_prior_calculator.chooseUnrooted();

        // Avoid recalculating if there has been no change in the number of leaves
        unsigned nleaves = tree->numLeaves();
        unsigned ntax = _topo_prior_calculator.getNTax();
        if (ntax != nleaves)
            _topo_prior_calculator.setNTax(nleaves);
        
        double log_topology_prior = 0.0;
        if (_likelihood->getModel()->isAllowPolytomies()) {
            unsigned m = tree->numInternals();
            log_topology_prior = _topo_prior_calculator.getLogNormalizedTopologyPrior(m);
        }
        else {
            log_topology_prior = -1.0*_topo_prior_calculator.getLogSaturatedCount(nleaves);
        }
        return log_topology_prior;
    }

    inline void Updater::debugPriorCalculation() {
        if (_name == "Tree Length") {
            Tree::SharedPtr tree = _tree_manipulator->getTree();
            assert(tree);

            double TL = _tree_manipulator->calcTreeLength();
            double num_edges = _tree_manipulator->countEdges();

#if defined(HOLDER_ETAL_PRIOR)
            double log_prior = calcLogEdgeLengthPrior();
            
            // Assume prior on each edge length is Exp(r) independently
            double exponential_rate = _prior_parameters[0];
            double check_log_prior = num_edges*log(exponential_rate) - exponential_rate*TL;
            ::om.outputConsole(boost::format("\n%24s: log prior = %.5f") % "Edge Lengths" % log_prior);
            ::om.outputConsole(boost::format("\n%37s %d*log(%.5f) - %.5f*%.5f") % "=" % num_edges % exponential_rate % exponential_rate % TL);
            ::om.outputConsole(boost::format("\n%37s %.5f\n") % "=" % check_log_prior);
#else
            std::pair<double, double> log_prior = calcLogEdgeLengthPrior();
            double a = _prior_parameters[0];    // shape of Gamma prior on TL
            double b = _prior_parameters[1];    // scale of Gamma prior on TL
            double c = _prior_parameters[2];    // parameter of Dirichlet prior on edge length proportions

            // Calculate Gamma prior on tree length (TL)
            double check_log_TL_prior = (a - 1.0)*log(TL) - TL/b - a*log(b) - std::lgamma(a);
            ::om.outputConsole(boost::format("\n%24s: Gamma shape = %.5f, scale = %.5f") % "Tree Length" % a % b);
            ::om.outputConsole(boost::format("\n%24s  value = %.5f") % " " % TL);
            ::om.outputConsole(boost::format("\n%24s  log prior = %.5f") % " " % log_prior.first);
            ::om.outputConsole(boost::format("\n%37s (%.5f - 1)*log(%.5f) - %.5f/%.5f - %.5f*log(%.5f) - lgamma(%.5f)") % "=" % a % TL % TL % b % a % b % a);
            ::om.outputConsole(boost::format("\n%37s %.5f\n") % "=" % check_log_TL_prior);
            
            ::om.outputConsole(boost::format("\n%24s %12s %12s %12s\n") % "Edge Proportions" % "value" % "parameter" %  "log(term)");
            
            // Calculate Dirichlet prior on edge length proportions
            //
            // Note that, for n edges, the Dirichlet prior density is
            //
            // p1^{c-1} p2^{c-1} ... pn^{c-1}
            // ------------------------------
            //    Gamma(c)^n / Gamma(n*c)
            //
            // where n = num_edges, pk = edge length k / TL and Gamma is the Gamma function.
            // If c == 1, then both numerator and denominator equal 1, so it is pointless
            // do loop over edge lengths.
            double sum_log_terms = 0.0;
            unsigned i = 0;
            for (auto nd : tree->_preorder) {
                double edge_length_proportion = nd->_edge_length/TL;
                double log_term = (c - 1.0)*log(edge_length_proportion);
                sum_log_terms += log_term;
                std::string s = boost::str(boost::format("%d") % i);
                ::om.outputConsole(boost::format("%24s %12.5f %12.5f %12.5f\n") % s % edge_length_proportion % c % log_term);
                ++i;
            }
            assert(i == num_edges);
            double logC = std::lgamma(num_edges*c) - std::lgamma(c)*num_edges;
            ::om.outputConsole(boost::format("%24s %38.5f\n") % "log-constant" % logC);
            sum_log_terms += logC;
            ::om.outputConsole(boost::format("%24s %38.5f\n") % "sum" % sum_log_terms);
            ::om.outputConsole(boost::format("%24s %38.5f\n") % "log-prior" % log_prior.second);
#endif
        }
        else {
            ::om.outputConsole(boost::format("%24s: debugPriorCalculation not yet implemented\n") % _name);
        }
    }

#if defined(HOLDER_ETAL_PRIOR)
    inline double Updater::calcLogEdgeLengthPrior() const {
#else
    inline std::pair<double,double> Updater::calcLogEdgeLengthPrior() const {
#endif
        Tree::SharedPtr tree = _tree_manipulator->getTree();
        assert(tree);

        double TL = _tree_manipulator->calcTreeLength();
        double num_edges = _tree_manipulator->countEdges();

#if defined(HOLDER_ETAL_PRIOR)
        // Assume prior on each edge length is Exp(r) independently
        double exponential_rate = _prior_parameters[0];
        double log_exp_prior = num_edges*log(exponential_rate) - exponential_rate*TL;
        return log_exp_prior;
#else
        double a = _prior_parameters[0];    // shape of Gamma prior on TL
        double b = _prior_parameters[1];    // scale of Gamma prior on TL
        double c = _prior_parameters[2];    // parameter of Dirichlet prior on edge length proportions

        // Calculate Gamma prior on tree length (TL)
        double log_gamma_prior_on_TL = (a - 1.0)*log(TL) - TL/b - a*log(b) - std::lgamma(a);
        
        // Calculate Dirichlet prior on edge length proportions
        //
        // Note that, for n edges, the Dirichlet prior density is
        //
        // p1^{c-1} p2^{c-1} ... pn^{c-1}
        // ------------------------------
        //    Gamma(c)^n / Gamma(n*c)
        //
        // where n = num_edges, pk = edge length k / TL and Gamma is the Gamma function.
        // If c == 1, then both numerator and denominator equal 1, so it is pointless
        // do loop over edge lengths.
        double log_edge_length_proportions_prior = std::lgamma(num_edges*c);
        if (c != 1.0) {
            for (auto nd : tree->_preorder) {
                double edge_length_proportion = nd->_edge_length/TL;
                log_edge_length_proportions_prior += (c - 1.0)*log(edge_length_proportion);
            }
            log_edge_length_proportions_prior -= std::lgamma(c)*num_edges;
        }

        return std::make_pair(log_gamma_prior_on_TL, log_edge_length_proportions_prior);
#endif
    }

    inline double Updater::getLogZero() {
        return _log_zero;
    }
    
}
