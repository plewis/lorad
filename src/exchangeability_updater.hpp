#pragma once

#include "dirichlet_updater.hpp"

namespace lorad {
    
    class ExchangeabilityUpdater : public DirichletUpdater {
        public:
            typedef std::shared_ptr< ExchangeabilityUpdater > SharedPtr;

                                            ExchangeabilityUpdater(QMatrix::SharedPtr qmatrix);
                                            ~ExchangeabilityUpdater();
                                        
            void                            debugPriorCalculation();
            
        private:
        
            void                            pullFromModel();
            void                            pushToModel();

            QMatrix::SharedPtr              _qmatrix;
        };

    inline ExchangeabilityUpdater::ExchangeabilityUpdater(QMatrix::SharedPtr qmatrix) { 
        DirichletUpdater::clear();
        _name = "Exchangeabilities";
        assert(qmatrix);
        _qmatrix = qmatrix;
    }

    inline ExchangeabilityUpdater::~ExchangeabilityUpdater() {
    }
    
    inline void ExchangeabilityUpdater::debugPriorCalculation() {
        ::om.outputConsole(boost::format("\n%24s %12s %12d %12s\n") % _name % "value" % "parameter" % "log(term)");
        double log_prior = calcLogPrior();
        
        assert(_curr_point.size() == _prior_parameters.size());
        assert(_curr_point.size() == 6);    // doesn't yet handle codon models

        double sum_log_terms = 0.0;
        double sum_params = 0.0;
        double sum_lgamma_params = 0.0;
        
        double logtermAC = (_prior_parameters[0] - 1.0)*log(_curr_point[0]);
        sum_log_terms += logtermAC;
        sum_params += _prior_parameters[0];
        sum_lgamma_params += lgamma(_prior_parameters[0]);
        ::om.outputConsole(boost::format("%24s %12.5f %12.5f %12.5f\n") % "rAC" % _curr_point[0] % _prior_parameters[0] % logtermAC);

        double logtermAG = (_prior_parameters[1] - 1.0)*log(_curr_point[1]);
        sum_log_terms += logtermAG;
        sum_params += _prior_parameters[1];
        sum_lgamma_params += lgamma(_prior_parameters[1]);
        ::om.outputConsole(boost::format("%24s %12.5f %12.5f %12.5f\n") % "rAG" % _curr_point[1] % _prior_parameters[1] % logtermAG);

        double logtermAT = (_prior_parameters[2] - 1.0)*log(_curr_point[2]);
        sum_log_terms += logtermAT;
        sum_params += _prior_parameters[2];
        sum_lgamma_params += lgamma(_prior_parameters[2]);
        ::om.outputConsole(boost::format("%24s %12.5f %12.5f %12.5f\n") % "rAT" % _curr_point[2] % _prior_parameters[2] % logtermAT);

        double logtermCG = (_prior_parameters[3] - 1.0)*log(_curr_point[3]);
        sum_log_terms += logtermCG;
        sum_params += _prior_parameters[3];
        sum_lgamma_params += lgamma(_prior_parameters[3]);
        ::om.outputConsole(boost::format("%24s %12.5f %12.5f %12.5f\n") % "rCG" % _curr_point[3] % _prior_parameters[3] % logtermCG);

        double logtermCT = (_prior_parameters[4] - 1.0)*log(_curr_point[4]);
        sum_log_terms += logtermCT;
        sum_params += _prior_parameters[4];
        sum_lgamma_params += lgamma(_prior_parameters[4]);
        ::om.outputConsole(boost::format("%24s %12.5f %12.5f %12.5f\n") % "rCT" % _curr_point[4] % _prior_parameters[4] % logtermCT);

        double logtermGT = (_prior_parameters[5] - 1.0)*log(_curr_point[5]);
        sum_log_terms += logtermGT;
        sum_params += _prior_parameters[5];
        sum_lgamma_params += lgamma(_prior_parameters[5]);
        ::om.outputConsole(boost::format("%24s %12.5f %12.5f %12.5f\n") % "rGT" % _curr_point[5] % _prior_parameters[5] % logtermGT);

        double logC = lgamma(sum_params) - sum_lgamma_params;
        sum_log_terms += logC;
        ::om.outputConsole(boost::format("%24s %38.5f\n") % "log-constant" % logC);
        ::om.outputConsole(boost::format("%24s %38.5f\n") % "sum" % sum_log_terms);
        ::om.outputConsole(boost::format("%24s %38.5f\n") % "log-prior" % log_prior);
    }

    inline void ExchangeabilityUpdater::pullFromModel() {
        QMatrix::freq_xchg_ptr_t xchg = _qmatrix->getExchangeabilitiesSharedPtr();
        _curr_point.assign(xchg->begin(), xchg->end());
    }
    
    inline void ExchangeabilityUpdater::pushToModel() {
        _qmatrix->setExchangeabilities(_curr_point);
    }

}
