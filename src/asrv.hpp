#pragma once

#include "conditionals.hpp"

#include <vector>
#include <boost/math/distributions/gamma.hpp>

namespace lorad {

    class ASRV {

        public:
#if defined(POLGSS)
#if defined(HOLDER_ETAL_PRIOR)
            typedef std::vector<double>                 shape_refdist_t;
            typedef std::shared_ptr<shape_refdist_t>    shape_refdist_ptr_t;
#else
            typedef std::vector<double>                 ratevar_refdist_t;
            typedef std::shared_ptr<ratevar_refdist_t>  ratevar_refdist_ptr_t;
#endif
            typedef std::vector<double>                 pinvar_refdist_t;
            typedef std::shared_ptr<pinvar_refdist_t>   pinvar_refdist_ptr_t;
#endif

            typedef std::vector<double>                 rate_prob_t;
            typedef std::shared_ptr<double>             relrate_ptr_t;
#if defined(HOLDER_ETAL_PRIOR)
            typedef std::shared_ptr<double>             shape_ptr_t;
#else
            typedef std::shared_ptr<double>             ratevar_ptr_t;
#endif
            typedef std::shared_ptr<double>             pinvar_ptr_t;
            typedef std::shared_ptr<ASRV>               SharedPtr;
        
                                                ASRV();
            virtual                             ~ASRV();
        
            void                                clear();
        
            void                                setNumCateg(unsigned ncateg);
            unsigned                            getNumCateg() const;

#if defined(HOLDER_ETAL_PRIOR)
            void                                setShapeSharedPtr(shape_ptr_t shape);
            void                                setShape(double v);
            const shape_ptr_t                   getShapeSharedPtr() const;
            double                              getShape() const;
            void                                fixShape(bool is_fixed);
            bool                                isFixedShape() const;
#else
            void                                setRateVarSharedPtr(ratevar_ptr_t ratevar);
            void                                setRateVar(double v);
            const ratevar_ptr_t                 getRateVarSharedPtr() const;
            double                              getRateVar() const;
            void                                fixRateVar(bool is_fixed);
            bool                                isFixedRateVar() const;
#endif

            void                                setPinvarSharedPtr(pinvar_ptr_t pinvar);
            void                                setPinvar(double p);
            const pinvar_ptr_t                  getPinvarSharedPtr() const;
            double                              getPinvar() const;
            void                                fixPinvar(bool is_fixed);
            bool                                isFixedPinvar() const;

#if defined(POLGSS)
#if defined(HOLDER_ETAL_PRIOR)
            void                                setShapeRefDistParamsSharedPtr(shape_refdist_ptr_t shape_refdist_ptr);
            std::vector<double>                 getShapeRefDistParamsVect() const;
#else
            void                                setRateVarRefDistParamsSharedPtr(ratevar_refdist_ptr_t ratevar_refdist_ptr);
            std::vector<double>                 getRateVarRefDistParamsVect() const;
#endif
            void                                setPinvarRefDistParamsSharedPtr(pinvar_refdist_ptr_t pinvar_refdist_ptr);
            std::vector<double>                 getPinvarRefDistParamsVect() const;
#endif

            void                                setIsInvarModel(bool is_invar_model);
            bool                                getIsInvarModel() const;

            const double *                      getRates() const;
            const double *                      getProbs() const;

        private:
        
            virtual void                        recalcASRV();

            unsigned                            _num_categ;
            bool                                _invar_model;
        
#if defined(HOLDER_ETAL_PRIOR)
            bool                                _shape_fixed;
            shape_ptr_t                         _shape;
#else
            bool                                _ratevar_fixed;
            ratevar_ptr_t                       _ratevar;
#endif
            pinvar_ptr_t                        _pinvar;
        
            bool                                _pinvar_fixed;
        
            rate_prob_t                         _rates;
            rate_prob_t                         _probs;
#if defined(POLGSS)
#if defined(HOLDER_ETAL_PRIOR)
            shape_refdist_ptr_t                 _shape_refdist;
#else
            ratevar_refdist_ptr_t               _ratevar_refdist;
#endif
            pinvar_refdist_ptr_t                _pinvar_refdist;
#endif
    };
    
    inline ASRV::ASRV() {
        clear();
    }

    inline ASRV::~ASRV() {
    }

    inline void ASRV::clear() {
        // Rate homogeneity is the default
        _invar_model = false;
#if defined(HOLDER_ETAL_PRIOR)
        _shape_fixed = false;
        _shape = std::make_shared<double>(1.0);
#else
        _ratevar_fixed = false;
        _ratevar = std::make_shared<double>(1.0);
#endif
        _pinvar_fixed = false;
        _pinvar = std::make_shared<double>(0.0);
        _num_categ = 1;
#if defined(POLGSS)
#if defined(HOLDER_ETAL_PRIOR)
        ASRV::shape_refdist_t tmp = {1.0, 1.0};
        _shape_refdist = std::make_shared<ASRV::shape_refdist_t>(tmp);
#else
        ASRV::ratevar_refdist_t tmp = {1.0, 1.0};
        _ratevar_refdist = std::make_shared<ASRV::ratevar_refdist_t>(tmp);
#endif
        ASRV::pinvar_refdist_t tmp2 = {1.0, 1.0};
        _pinvar_refdist = std::make_shared<ASRV::pinvar_refdist_t>(tmp2);
#endif
        recalcASRV();
    }

#if defined(HOLDER_ETAL_PRIOR)
    inline const ASRV::shape_ptr_t ASRV::getShapeSharedPtr() const {
        return _shape;
    }

    inline double ASRV::getShape() const {
        assert(_shape);
        return *_shape;
    }
    
    inline void ASRV::setShapeSharedPtr(shape_ptr_t shape) {
        _shape = shape;
        recalcASRV();
    }
    
    inline void ASRV::setShape(double v) {
        *_shape = v;
        recalcASRV();
    }
    
    inline void ASRV::fixShape(bool is_fixed) {
        _shape_fixed = is_fixed;
    }

    inline bool ASRV::isFixedShape() const {
        return _shape_fixed;
    }
#else
    inline const ASRV::ratevar_ptr_t ASRV::getRateVarSharedPtr() const {
        return _ratevar;
    }

    inline double ASRV::getRateVar() const {
        assert(_ratevar);
        return *_ratevar;
    }
    
    inline void ASRV::setRateVarSharedPtr(ratevar_ptr_t ratevar) {
        _ratevar = ratevar;
        recalcASRV();
    }
    
    inline void ASRV::setRateVar(double v) {
        *_ratevar = v;
        recalcASRV();
    }
    
    inline void ASRV::fixRateVar(bool is_fixed) {
        _ratevar_fixed = is_fixed;
    }

    inline bool ASRV::isFixedRateVar() const {
        return _ratevar_fixed;
    }
#endif

    inline const ASRV::pinvar_ptr_t ASRV::getPinvarSharedPtr() const {
        return _pinvar;
    }

    inline double ASRV::getPinvar() const {
        assert(_pinvar);
        return *_pinvar;
    }

    inline const double * ASRV::getRates() const {
        return &_rates[0];
    }

    inline const double * ASRV::getProbs() const {
        return &_probs[0];
    }

    inline bool ASRV::getIsInvarModel() const {
        return _invar_model;
    }

    inline unsigned ASRV::getNumCateg() const {
        return _num_categ;
    }
    
    inline void ASRV::setNumCateg(unsigned ncateg) {
        _num_categ = ncateg;
        recalcASRV();
    }
    
    inline void ASRV::setPinvarSharedPtr(pinvar_ptr_t pinvar) {
        _pinvar = pinvar;
        recalcASRV();
    }
    
    inline void ASRV::setPinvar(double p) {
        *_pinvar = p;
        recalcASRV();
    }
    
    inline void ASRV::setIsInvarModel(bool is_invar_model) {
        _invar_model = is_invar_model;
        recalcASRV();
    }

    inline void ASRV::fixPinvar(bool is_fixed) {
        _pinvar_fixed = is_fixed;
    }

    inline bool ASRV::isFixedPinvar() const {
        return _pinvar_fixed;
    }
    
    inline void ASRV::recalcASRV() {
        // This implementation assumes discrete gamma among-site rate heterogeneity
        // using a _num_categ category discrete gamma distribution with equal category
        // probabilities and Gamma density with mean 1.0 and variance _rate_var.
        // If _invar_model is true, then rate probs will sum to 1 - _pinvar rather than 1
        // and the mean rate will be 1/(1 - _pinvar) rather than 1; the rest of the invariable
        // sites component of the model is handled outside the ASRV class.
        
        // _num_categ, _rate_var, and _pinvar must all have been assigned in order to compute rates and probs
#if defined(HOLDER_ETAL_PRIOR)
        if ( (!_shape) || (!_num_categ) || (!_pinvar) )
            return;
#else
        if ( (!_ratevar) || (!_num_categ) || (!_pinvar) )
            return;
#endif
        
        double pinvar = *_pinvar;
        assert(pinvar >= 0.0);
        assert(pinvar <  1.0);

        assert(_num_categ > 0);

        double equal_prob = 1.0/_num_categ;
        double mean_rate_variable_sites = 1.0;
        if (_invar_model)
            mean_rate_variable_sites /= (1.0 - pinvar);
        
        _rates.assign(_num_categ, mean_rate_variable_sites);
        _probs.assign(_num_categ, equal_prob);

#if defined(HOLDER_ETAL_PRIOR)
        double gamma_shape = *_shape;
        assert(gamma_shape >= 0.0);
        
        if (_num_categ == 1 || gamma_shape > 1000.0)
            return;

        double alpha = gamma_shape;
        double beta = 1.0/gamma_shape;
#else
        double rate_variance = *_ratevar;
        assert(rate_variance >= 0.0);
        
        if (_num_categ == 1 || rate_variance == 0.0)
            return;

        double alpha = 1.0/rate_variance;
        double beta = rate_variance;
#endif
    
        boost::math::gamma_distribution<> my_gamma(alpha, beta);
        boost::math::gamma_distribution<> my_gamma_plus(alpha + 1.0, beta);

        double cum_upper        = 0.0;
        double cum_upper_plus   = 0.0;
        double upper            = 0.0;
        double cum_prob         = 0.0;
        for (unsigned i = 1; i <= _num_categ; ++i) {
            double cum_lower_plus       = cum_upper_plus;
            double cum_lower            = cum_upper;
            cum_prob                    += equal_prob;

            if (i < _num_categ) {
                upper                   = boost::math::quantile(my_gamma, cum_prob);
                cum_upper_plus          = boost::math::cdf(my_gamma_plus, upper);
                cum_upper               = boost::math::cdf(my_gamma, upper);
            }
            else {
                cum_upper_plus          = 1.0;
                cum_upper               = 1.0;
            }

            double numer                = cum_upper_plus - cum_lower_plus;
            double denom                = cum_upper - cum_lower;
            double r_mean               = (denom > 0.0 ? (alpha*beta*numer/denom) : 0.0);
            _rates[i-1]        = r_mean*mean_rate_variable_sites;
        }
    }

#if defined(POLGSS)
#if defined(HOLDER_ETAL_PRIOR)
    inline void ASRV::setShapeRefDistParamsSharedPtr(ASRV::shape_refdist_ptr_t shape_refdist_params_ptr) {
        if (shape_refdist_params_ptr->size() != 2)
            throw XLorad(boost::format("Expecting 2 shape reference distribution parameters and got %d") % shape_refdist_params_ptr->size());
        _shape_refdist = shape_refdist_params_ptr;
    }
    
    inline std::vector<double> ASRV::getShapeRefDistParamsVect() const {
        return std::vector<double>(_shape_refdist->begin(), _shape_refdist->end());
    }
#else
    inline void ASRV::setRateVarRefDistParamsSharedPtr(ASRV::ratevar_refdist_ptr_t ratevar_refdist_params_ptr) {
        if (ratevar_refdist_params_ptr->size() != 2)
            throw XLorad(boost::format("Expecting 2 rate variance reference distribution parameters and got %d") % ratevar_refdist_params_ptr->size());
        _ratevar_refdist = ratevar_refdist_params_ptr;
    }
    
    inline std::vector<double> ASRV::getRateVarRefDistParamsVect() const {
        return std::vector<double>(_ratevar_refdist->begin(), _ratevar_refdist->end());
    }
#endif

    inline void ASRV::setPinvarRefDistParamsSharedPtr(ASRV::pinvar_refdist_ptr_t pinvar_refdist_params_ptr) {
        if (pinvar_refdist_params_ptr->size() != 2)
            throw XLorad(boost::format("Expecting 2 pinvar reference distribution parameters and got %d") % pinvar_refdist_params_ptr->size());
        _pinvar_refdist = pinvar_refdist_params_ptr;
    }
    
    inline std::vector<double> ASRV::getPinvarRefDistParamsVect() const {
        return std::vector<double>(_pinvar_refdist->begin(), _pinvar_refdist->end());
    }
#endif
}
