#pragma once

#include "dirichlet_updater.hpp"

namespace lorad {
    
    class ExchangeabilityUpdater : public DirichletUpdater {
        public:
            typedef std::shared_ptr< ExchangeabilityUpdater > SharedPtr;

                                            ExchangeabilityUpdater(QMatrix::SharedPtr qmatrix);
                                            ~ExchangeabilityUpdater();
                                        
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
    
    inline void ExchangeabilityUpdater::pullFromModel() {
        QMatrix::freq_xchg_ptr_t xchg = _qmatrix->getExchangeabilitiesSharedPtr();
        _curr_point.assign(xchg->begin(), xchg->end());
    }
    
    inline void ExchangeabilityUpdater::pushToModel() {
        _qmatrix->setExchangeabilities(_curr_point);
    }

}
