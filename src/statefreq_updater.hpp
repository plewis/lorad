#pragma once

#include "dirichlet_updater.hpp"

namespace lorad {
    
    class StateFreqUpdater : public DirichletUpdater {
        public:
            typedef std::shared_ptr< StateFreqUpdater > SharedPtr;

                                            StateFreqUpdater(QMatrix::SharedPtr qmatrix);
                                            ~StateFreqUpdater();
                                            
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
    
    inline void StateFreqUpdater::pullFromModel() {
        QMatrix::freq_xchg_ptr_t freqs = _qmatrix->getStateFreqsSharedPtr();
        _curr_point.assign(freqs->begin(), freqs->end());
    }
    
    inline void StateFreqUpdater::pushToModel() {
        _qmatrix->setStateFreqs(_curr_point);
    }
    
}
