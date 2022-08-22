#pragma once

#include "dirichlet_updater.hpp"

namespace lorad {
    
    class StateFreqUpdater : public DirichletUpdater {
        public:
            typedef std::shared_ptr< StateFreqUpdater > SharedPtr;

                                            StateFreqUpdater(QMatrix::SharedPtr qmatrix);
                                            ~StateFreqUpdater();
                                            
            virtual void                    debugPriorCalculation();

        private:
        
            void                            pullFromModel();
            void                            pushToModel();

            QMatrix::SharedPtr              _qmatrix;
        };

    inline StateFreqUpdater::StateFreqUpdater(QMatrix::SharedPtr qmatrix) {
        DirichletUpdater::clear();
        _name = "State Frequencies";
        assert(qmatrix);
        _qmatrix = qmatrix;
    }

    inline StateFreqUpdater::~StateFreqUpdater() {
    }
    
    inline void StateFreqUpdater::debugPriorCalculation() {
        ::om.outputConsole(boost::format("\n%24s %12s %12d %12s\n") % _name % "value" % "parameter" % "log(term)");
        double log_prior = calcLogPrior();
        
        assert(_curr_point.size() == _prior_parameters.size());
        assert(_curr_point.size() == 4);    // doesn't yet handle codon models

        double sum_log_terms = 0.0;
        
        double logtermA = (_prior_parameters[0] - 1.0)*log(_curr_point[0]);
        sum_log_terms += logtermA;
        ::om.outputConsole(boost::format("%24s %12.5f %12.5f %12.5f\n") % "piA" % _curr_point[0] % _prior_parameters[0] % logtermA);

        double logtermC = (_prior_parameters[1] - 1.0)*log(_curr_point[1]);
        sum_log_terms += logtermC;
        ::om.outputConsole(boost::format("%24s %12.5f %12.5f %12.5f\n") % "piC" % _curr_point[1] % _prior_parameters[1] % logtermC);

        double logtermG = (_prior_parameters[2] - 1.0)*log(_curr_point[2]);
        sum_log_terms += logtermG;
        ::om.outputConsole(boost::format("%24s %12.5f %12.5f %12.5f\n") % "piG" % _curr_point[2] % _prior_parameters[2] % logtermG);

        double logtermT = (_prior_parameters[3] - 1.0)*log(_curr_point[3]);
        sum_log_terms += logtermT;
        ::om.outputConsole(boost::format("%24s %12.5f %12.5f %12.5f\n") % "piT" % _curr_point[3] % _prior_parameters[3] % logtermT);

        double logC = lgamma(_prior_parameters[0] + _prior_parameters[1] + _prior_parameters[2] + _prior_parameters[3]) - lgamma(_prior_parameters[0]) - lgamma(_prior_parameters[1]) - lgamma(_prior_parameters[2]) - lgamma(_prior_parameters[3]);
        sum_log_terms += logC;
        ::om.outputConsole(boost::format("%24s %38.5f\n") % "log-constant" % logC);
        ::om.outputConsole(boost::format("%24s %38.5f\n") % "sum" % sum_log_terms);
        ::om.outputConsole(boost::format("%24s %38.5f\n") % "log-prior" % log_prior);
    }

    inline void StateFreqUpdater::pullFromModel() {
        QMatrix::freq_xchg_ptr_t freqs = _qmatrix->getStateFreqsSharedPtr();
        _curr_point.assign(freqs->begin(), freqs->end());
    }
    
    inline void StateFreqUpdater::pushToModel() {
        _qmatrix->setStateFreqs(_curr_point);
    }
    
}
