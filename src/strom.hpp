#pragma once

//#define USE_BOOST_REGEX
#define HPD

#include <iostream>
#include "data.hpp"
#include "likelihood.hpp"
#include "tree_summary.hpp"
#include "partition.hpp"
#include "lot.hpp"
#include "chain.hpp"
#include "output_manager.hpp"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

namespace strom {

#if defined(HPD)
    struct Kernel {
        std::string _title;
        double      _log_likelihood;
        double      _log_prior;
        double      _log_jacobian_log_transformation;
        double      _log_jacobian_standardization;
        
        Kernel() : _log_likelihood(0.0), _log_prior(0.0), _log_jacobian_log_transformation(0.0), _log_jacobian_standardization(0.0) {}
        
        Kernel(double lnL, double lnP, double lnJlog, double lnJstd)
          : _log_likelihood(lnL),
            _log_prior(lnP),
            _log_jacobian_log_transformation(lnJlog),
            _log_jacobian_standardization(lnJstd) {
        }
            
        double logKernel() const {
            return _log_likelihood + _log_prior + _log_jacobian_log_transformation + _log_jacobian_standardization;
        }
        
        double logJacobian() const {
            return _log_jacobian_log_transformation + _log_jacobian_standardization;
        }
    };
    
    std::ostream & operator<<(std::ostream & os, const Kernel & k) {
        if (k._title.size() > 0)
            os << "\n" << k._title << ":\n";
        else
            os << "\n(untitled Kernel object):\n";
        os << boost::format("  _log_likelihood = %.5f\n") % k._log_likelihood;
        os << boost::format("  _log_prior      = %.5f\n") % k._log_prior;
        os << boost::format("  _log_jacobian_log_transformation = %.5f\n") % k._log_jacobian_log_transformation;
        os << boost::format("  _log_jacobian_standardization    = %.5f\n") % k._log_jacobian_standardization;
        os << boost::format("  logKernel() = %.5f\n") % k.logKernel();
        os << std::endl;
        return os;
    }
    
    struct ParameterSample {
        unsigned         _iteration;
        Kernel           _kernel;
        Eigen::VectorXd  _param_vect;
        
        // Define greater-than operator so that a vector of ParameterSample objects can be sorted
        bool operator>(const ParameterSample & other) const {
            return _kernel.logKernel() > other._kernel.logKernel();
        }
    };

    struct WorkingSpaceSubset {
        double              _lower_logq;
        double              _upper_logq;
        double              _norm_min;
        double              _norm_max;
        double              _representative;
        double              _log_base_volume;
        unsigned            _begin;
        unsigned            _end;
        unsigned            _ndarts_thrown;
        unsigned            _ndarts_inside;
        std::vector<double> _norms;
        
        WorkingSpaceSubset() : _lower_logq(0.0), _upper_logq(0.0), _norm_min(0.0), _norm_max(0.0), _representative(0.0), _log_base_volume(0.0), _begin(0), _end(0), _ndarts_thrown(0), _ndarts_inside(0) {}
    };

#endif

    class Strom {
        public:
                                                    Strom();
                                                    ~Strom();

            void                                    clear();
            void                                    processCommandLineOptions(int argc, const char * argv[]);
            void                                    run();
                    
        private:
            bool                                    processAssignmentString(Model::SharedPtr m, const std::string & which, const std::string & definition);
            void                                    handleAssignmentStrings(Model::SharedPtr m, const boost::program_options::variables_map & vm, std::string label, const std::vector<std::string> & definitions, std::string default_definition); 
            bool                                    splitAssignmentString(const std::string & definition, std::vector<std::string> & vector_of_subset_names, std::vector<double>  & vector_of_values);
            void                                    sample(unsigned iter, Chain & chain);

            void                                    readData();
            void                                    readTrees();
            void                                    showPartitionInfo();
            void                                    showBeagleInfo();
            void                                    showMCMCInfo();
            void                                    calcHeatingPowers();
            void                                    calcMarginalLikelihood();
            void                                    initChains();
            void                                    startTuningChains();
            void                                    stopTuningChains();
            void                                    stepChains(unsigned iteration, bool sampling);
            void                                    swapChains();
            void                                    stopChains();
            void                                    swapSummary() const;
            void                                    showChainTuningInfo() const;

#if defined(HPD)
            void                                    saveLogTransformedParameters(unsigned iteration, double logLike, double logPrior, Model::SharedPtr model, TreeManip::SharedPtr tm);
            void                                    saveStandardizedSamples() const;
            void                                    inputStandardizedSamples();
            void                                    standardizeParameters();
            void                                    defineShells();
            void                                    throwDarts();
            Kernel                                  calcLogTransformedKernel(Eigen::VectorXd & x);
            double                                  calcLogSum(const std::vector<double> & logx_vect);
            double                                  histogramMethod();
#endif

            double                                  _expected_log_likelihood;
            
            double                                  _topo_prior_C;
            bool                                    _allow_polytomies;
            bool                                    _resolution_class_prior;

            std::string                             _data_file_name;
            std::string                             _tree_file_name;
            Partition::SharedPtr                    _partition;

            Data::SharedPtr                         _data;
            std::vector<Likelihood::SharedPtr>      _likelihoods;
            TreeSummary::SharedPtr                  _tree_summary;
            Lot::SharedPtr                          _lot;

            unsigned                                _random_seed;
            unsigned                                _num_iter;
            unsigned                                _print_freq;
            unsigned                                _sample_freq;

            unsigned                                _num_burnin_iter; 
            bool                                    _using_stored_data;
            bool                                    _use_gpu;
            unsigned                                _nstones;
            double                                  _ss_alpha;

            bool                                    _ambig_missing;
            unsigned                                _nchains;
            double                                  _heating_lambda;
            std::vector<Chain>                      _chains;
            std::vector<double>                     _heating_powers;
            std::vector<unsigned>                   _swaps;

            bool                                    _use_underflow_scaling;

            static std::string                      _program_name;
            static unsigned                         _major_version;
            static unsigned                         _minor_version;
            
            OutputManager::SharedPtr                _output_manager;
            
#if defined(HPD)
            bool                                    _skipMCMC;
            double                                  _coverage;
            unsigned                                _nshells;
            unsigned                                _ndarts;
            unsigned                                _nparams;
            unsigned                                _nsamples;
            std::string                             _param_file_name;

            typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> eigenMatrixXd_t;
            eigenMatrixXd_t                         _S; // var-cov matrix
            eigenMatrixXd_t                         _sqrtS;
            eigenMatrixXd_t                         _invSqrtS;
            double                                  _logDetSqrtS;
            
            typedef Eigen::VectorXd eigenVectorXd_t;
            eigenVectorXd_t                         _mean_transformed;
            eigenVectorXd_t                         _mode_transformed;
            
            std::vector< ParameterSample >          _log_transformed_parameters;
            std::vector< ParameterSample >          _standardized_parameters;
            std::vector<WorkingSpaceSubset>         _A;
#endif

    };
    
    inline Strom::Strom() {
        //std::cout << "Constructing a Strom" << std::endl;
        clear();
    }

    inline Strom::~Strom() {
        //std::cout << "Destroying a Strom" << std::endl;
    }

    inline void Strom::clear() {
        _data_file_name             = "";
        _tree_file_name             = "";
        _tree_summary               = nullptr;
        _partition.reset(new Partition());
        _use_gpu                    = true;
        _nstones                    = 0;
        _ss_alpha                   = 0.25;
        _ambig_missing              = true;
        _expected_log_likelihood    = 0.0;
        _data                       = nullptr;
        _use_underflow_scaling      = false;
        _lot                        = nullptr;
        _random_seed                = 1;
        _num_iter                   = 1000;
        _print_freq                 = 1;
        _sample_freq                = 1;
        _output_manager             = nullptr;
        
        _topo_prior_C               = 1.0;
        _allow_polytomies           = true;
        _resolution_class_prior     = true;

        _using_stored_data          = true;
        _likelihoods.clear();
        _num_burnin_iter            = 1000;
        _heating_lambda             = 0.5;
        _nchains                    = 1;
        _chains.resize(0);
        _heating_powers.resize(0);
        _swaps.resize(0);

#if defined(HPD)
        _skipMCMC                   = false;
        _coverage                   = 0.1;
        _nshells                    = 0;
        _ndarts                     = 10000;
        _nparams                    = 0;
        _nsamples                   = 0;
        _param_file_name            = "standardized_params.txt";
#endif
        
    }

    inline void Strom::processCommandLineOptions(int argc, const char * argv[]) {
        std::vector<std::string> partition_statefreq;
        std::vector<std::string> partition_rmatrix;
        std::vector<std::string> partition_omega;
        std::vector<std::string> partition_ratevar;
        std::vector<std::string> partition_pinvar;
        std::vector<std::string> partition_ncateg;
        std::vector<std::string> partition_subsets;
        std::vector<std::string> partition_relrates;
        std::vector<std::string> partition_tree;
        boost::program_options::variables_map vm;
        boost::program_options::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("version,v", "show program version")
            ("seed,z",        boost::program_options::value(&_random_seed)->default_value(1),   "pseudorandom number seed")
            ("niter,n",       boost::program_options::value(&_num_iter)->default_value(1000),   "number of MCMC iterations")
            ("printfreq",  boost::program_options::value(&_print_freq)->default_value(1),   "skip this many iterations before reporting progress")
            ("samplefreq",  boost::program_options::value(&_sample_freq)->default_value(1),   "skip this many iterations before sampling next")
            ("datafile,d",  boost::program_options::value(&_data_file_name)->required(), "name of a data file in NEXUS format")
            ("treefile,t",  boost::program_options::value(&_tree_file_name)->required(), "name of a tree file in NEXUS format")
            ("subset",  boost::program_options::value(&partition_subsets), "a string defining a partition subset, e.g. 'first:1-1234\3' or 'default[codon:standard]:1-3702'")
            ("ncateg,c", boost::program_options::value(&partition_ncateg), "number of categories in the discrete Gamma rate heterogeneity model")
            ("statefreq", boost::program_options::value(&partition_statefreq), "a string defining state frequencies for one or more data subsets, e.g. 'first,second:0.1,0.2,0.3,0.4'")
            ("omega", boost::program_options::value(&partition_omega), "a string defining the nonsynonymous/synonymous rate ratio omega for one or more data subsets, e.g. 'first,second:0.1'")
            ("rmatrix", boost::program_options::value(&partition_rmatrix), "a string defining the rmatrix for one or more data subsets, e.g. 'first,second:1,2,1,1,2,1'")
            ("ratevar", boost::program_options::value(&partition_ratevar), "a string defining the among-site rate variance for one or more data subsets, e.g. 'first,second:2.5'")
            ("pinvar", boost::program_options::value(&partition_pinvar), "a string defining the proportion of invariable sites for one or more data subsets, e.g. 'first,second:0.2'")
            ("relrate", boost::program_options::value(&partition_relrates), "a string defining the (unnormalized) relative rates for all data subsets (e.g. 'default:3,1,6').")
            ("tree", boost::program_options::value(&partition_tree), "the index of the tree in the tree file (first tree has index = 1)")
            ("topopriorC", boost::program_options::value(&_topo_prior_C)->default_value(1.0), "topology prior C: tree (or resolution class) with m internal nodes has probability C time greater than tree (or resolution class) with m+1 internal nodes.")
            ("allowpolytomies", boost::program_options::value(&_allow_polytomies)->default_value(true), "yes or no; if yes, then topopriorC and polytomyprior are used, otherwise topopriorC and polytomyprior are ignored")
            ("resclassprior", boost::program_options::value(&_resolution_class_prior)->default_value(true), "if yes, topologypriorC will apply to resolution classes; if no, topologypriorC will apply to individual tree topologies")
            ("expectedLnL", boost::program_options::value(&_expected_log_likelihood)->default_value(0.0), "log likelihood expected")
            ("nchains",       boost::program_options::value(&_nchains)->default_value(1),                "number of chains")
            ("heatfactor",    boost::program_options::value(&_heating_lambda)->default_value(0.5),          "determines how hot the heated chains are")
            ("burnin",        boost::program_options::value(&_num_burnin_iter)->default_value(100),         "number of iterations used to burn in chains")
            ("usedata",       boost::program_options::value(&_using_stored_data)->default_value(true),      "use the stored data in calculating likelihoods (specify no to explore the prior)")
            ("gpu",           boost::program_options::value(&_use_gpu)->default_value(true),                "use GPU if available")
            ("ambigmissing",  boost::program_options::value(&_ambig_missing)->default_value(true),          "treat all ambiguities as missing data")
            ("underflowscaling",  boost::program_options::value(&_use_underflow_scaling)->default_value(true),          "scale site-likelihoods to prevent underflow (slower but safer)")
            ("nstones", boost::program_options::value(&_nstones)->default_value(false),                "use heated chains to compute marginal likelihood with the steppingstone method using nstones steppingstone ratios")
            ("ssalpha", boost::program_options::value(&_ss_alpha)->default_value(0.25),                "determines how bunched steppingstone chain powers are toward the prior: chain k of K total chains has power (k/K)^{1/ssalpha}")
#if defined(HPD)
            ("nshells", boost::program_options::value(&_nshells)->default_value(5), "the number of subsets of the working parameter space")
            ("coverage", boost::program_options::value(&_coverage)->default_value(0.95), "the fraction of samples used to construct the working parameter space")
            ("ndarts", boost::program_options::value(&_ndarts)->default_value(1000), "the number of \"darts\" to throw at each shell to determine what fraction of that shell's volume is inside the working parameter space subset")
            ("skipmcmc", boost::program_options::value(&_skipMCMC)->default_value(false),                "estimate marginal likelihood using the HPD histogram method from parameter vectors previously saved in paramfile (only used if marglike is yes)")
#endif
        ;
        boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
        try {
            const boost::program_options::parsed_options & parsed = boost::program_options::parse_config_file< char >("hpdml.conf", desc, false);
            boost::program_options::store(parsed, vm);
        }
        catch(boost::program_options::reading_file & x) {
            std::cout << "Note: configuration file (hpdml.conf) not found" << std::endl;
        }
        boost::program_options::notify(vm);

        // If user specified --help on command line, output usage summary and quit
        if (vm.count("help") > 0) {
            std::cout << desc << "\n";
            std::exit(1);
        }

        // If user specified --version on command line, output version and quit
        if (vm.count("version") > 0) {
            std::cout << boost::str(boost::format("This is %s version %d.%d") % _program_name % _major_version % _minor_version) << std::endl;
            std::exit(1);
        }
    
        // If user specified --subset on command line, break specified partition subset 
        // definition into name and character set string and add to _partition
        if (vm.count("subset") > 0) {
            _partition.reset(new Partition());
            for (auto s : partition_subsets) {
                _partition->parseSubsetDefinition(s);
            }
        }

        // Be sure number of chains is greater than or equal to 1
        if (_nchains < 1)
            throw XStrom("nchains must be a positive integer greater than 0");

        // Be sure number of stones is greater than or equal to 0
        if (_nstones < 0)
            throw XStrom("nstones must be a positive integer greater than or equal to 0");

#if defined(HPD)
        // Can't set nstones > 0 and nshells > 0 at the same time
        if (_nstones > 0 && _nshells > 0)
            throw XStrom("Cannot specify the steppingstone marginal likelihood method (nstones > 0) and the HPD marginal likelihood method (nshells > 0) at the same time; please choose one or the other");
#endif

        // If number of stones is greater than 0, then set _nchains to that value
        if (_nstones > 0) {
            std::cout << (boost::format("\nNumber of chains was set to the specified number of stones (%d)\n") % _nstones) << std::endl;
            _nchains = (unsigned)_nstones;
        }

        // Be sure heatfactor is between 0 and 1
        if (_heating_lambda <= 0.0 || _heating_lambda > 1.0)
            throw XStrom("heatfactor must be a real number in the interval (0.0,1.0]");
        
        if (!_using_stored_data)
            std::cout << "\n*** Not using stored data (posterior = prior) ***\n" << std::endl;
            
        // Allocate a separate model for each chain
        for (unsigned c = 0; c < _nchains; c++) {
            Likelihood::SharedPtr likelihood = Likelihood::SharedPtr(new Likelihood());
            likelihood->setPreferGPU(_use_gpu);
            likelihood->setAmbiguityEqualsMissing(_ambig_missing);
            Model::SharedPtr m = likelihood->getModel();
            m->setSubsetDataTypes(_partition->getSubsetDataTypes());
            handleAssignmentStrings(m, vm, "statefreq", partition_statefreq, "default:equal");
            handleAssignmentStrings(m, vm, "rmatrix",   partition_rmatrix,   "default:equal");
            handleAssignmentStrings(m, vm, "omega",     partition_omega,     "default:0.1"  );
            handleAssignmentStrings(m, vm, "ncateg",    partition_ncateg,    "default:1"    );
            handleAssignmentStrings(m, vm, "ratevar",   partition_ratevar,   "default:1.0"  );
            handleAssignmentStrings(m, vm, "pinvar",    partition_pinvar,    "default:0.0"  );
            handleAssignmentStrings(m, vm, "relrate",   partition_relrates,  "default:equal");
            handleAssignmentStrings(m, vm, "tree",      partition_tree,      "default:1");
            _likelihoods.push_back(likelihood);
        }
    }
    
    inline void Strom::handleAssignmentStrings(Model::SharedPtr m, const boost::program_options::variables_map & vm, std::string label, const std::vector<std::string> & definitions, std::string default_definition) {
        if (vm.count(label) > 0) {
            bool first = true;
            for (auto s : definitions) {
                bool is_default = processAssignmentString(m, label, s);
                if (is_default && !first)
                    throw XStrom(boost::format("default specification must be first %s encountered") % label);
                first = false;
            }
        }
        else {
            processAssignmentString(m, label, default_definition);
        }
    }
    
    inline bool Strom::processAssignmentString(Model::SharedPtr m, const std::string & which, const std::string & definition) {
        unsigned num_subsets_defined = _partition->getNumSubsets();
        std::vector<std::string> vector_of_subset_names;
        std::vector<double> vector_of_values;
        bool fixed = splitAssignmentString(definition, vector_of_subset_names, vector_of_values);
        
        if (vector_of_values.size() == 1 && vector_of_values[0] == -1 && !(which == "statefreq" || which == "rmatrix" || which == "relrate"))
            throw XStrom("Keyword equal is only allowed for statefreq, rmatrix, and relrate");

        // Assign values to subsets in model
        bool default_found = false;
        if (which == "statefreq") {
            QMatrix::freq_xchg_ptr_t freqs = std::make_shared<QMatrix::freq_xchg_t>(vector_of_values);
            if (vector_of_subset_names[0] == "default") {
                default_found = true;
                for (unsigned i = 0; i < num_subsets_defined; i++)
                    m->setSubsetStateFreqs(freqs, i, fixed);
            }
            else {
                for (auto s : vector_of_subset_names) {
                    m->setSubsetStateFreqs(freqs, _partition->findSubsetByName(s), fixed);
                }
            }
        }
        else if (which == "rmatrix") {
            QMatrix::freq_xchg_ptr_t xchg = std::make_shared<QMatrix::freq_xchg_t>(vector_of_values);
            if (vector_of_subset_names[0] == "default") {
                default_found = true;
                for (unsigned i = 0; i < num_subsets_defined; i++)
                    m->setSubsetExchangeabilities(xchg, i, fixed);
            }
            else {
                for (auto s : vector_of_subset_names) {
                    m->setSubsetExchangeabilities(xchg, _partition->findSubsetByName(s), fixed);
                }
            }
        }
        else if (which == "omega") {
            if (vector_of_values.size() > 1)
                throw XStrom(boost::format("expecting 1 value for omega, found %d values") % vector_of_values.size());
            QMatrix::omega_ptr_t omega = std::make_shared<QMatrix::omega_t>(vector_of_values[0]);
            if (vector_of_subset_names[0] == "default") {
                default_found = true;
                for (unsigned i = 0; i < num_subsets_defined; i++)
                    m->setSubsetOmega(omega, i, fixed);
            }
            else {
                for (auto s : vector_of_subset_names) {
                    m->setSubsetOmega(omega, _partition->findSubsetByName(s), fixed);
                }
            }
        }
        else if (which == "pinvar") {
            if (vector_of_values.size() > 1)
                throw XStrom(boost::format("expecting 1 value for pinvar, found %d values") % vector_of_values.size());
            ASRV::pinvar_ptr_t p = std::make_shared<double>(vector_of_values[0]);
            bool invar_model = (*p > 0);
            if (vector_of_subset_names[0] == "default") {
                default_found = true;
                for (unsigned i = 0; i < num_subsets_defined; i++) {
                    m->setSubsetIsInvarModel(invar_model, i);
                    m->setSubsetPinvar(p, i, fixed);
                }
            }
            else {
                for (auto s : vector_of_subset_names) {
                    unsigned i = _partition->findSubsetByName(s);
                    m->setSubsetIsInvarModel(invar_model, i);
                    m->setSubsetPinvar(p, i, fixed);
                }
            }
        }
        else if (which == "ratevar") {
            if (vector_of_values.size() > 1)
                throw XStrom(boost::format("expecting 1 value for ratevar, found %d values") % vector_of_values.size());
            ASRV::ratevar_ptr_t rv = std::make_shared<double>(vector_of_values[0]);
            if (vector_of_subset_names[0] == "default") {
                default_found = true;
                for (unsigned i = 0; i < num_subsets_defined; i++)
                    m->setSubsetRateVar(rv, i, fixed);
            }
            else {
                for (auto s : vector_of_subset_names) {
                    m->setSubsetRateVar(rv, _partition->findSubsetByName(s), fixed);
                }
            }
        }
        else if (which == "ncateg") {
            if (vector_of_values.size() > 1)
                throw XStrom(boost::format("expecting 1 value for ncateg, found %d values") % vector_of_values.size());
            unsigned ncat = vector_of_values[0];
            if (vector_of_subset_names[0] == "default") {
                default_found = true;
                for (unsigned i = 0; i < num_subsets_defined; i++)
                    m->setSubsetNumCateg(ncat, i);
            }
            else {
                for (auto s : vector_of_subset_names) {
                    m->setSubsetNumCateg(ncat, _partition->findSubsetByName(s));
                }
            }
        }
        else if (which == "tree") {
            if (vector_of_values.size() > 1)
                throw XStrom(boost::format("expecting 1 value for tree, found %d values") % vector_of_values.size());
            unsigned tree_index = vector_of_values[0];
            assert(tree_index > 0);
            m->setTreeIndex(tree_index - 1, fixed);
            if (vector_of_subset_names[0] != "default")
                throw XStrom("tree must be assigned to default only");
        }
        else {
            assert(which == "relrate");
            if (vector_of_subset_names[0] != "default")
                throw XStrom("relrate must be assigned to default only");
            m->setSubsetRelRates(vector_of_values, fixed);
        }

        return default_found;
    }

    inline bool Strom::splitAssignmentString(const std::string & definition, std::vector<std::string> & vector_of_subset_names, std::vector<double>  & vector_of_values) {
        // Split subset names from definition
        std::vector<std::string> twoparts;
        boost::split(twoparts, definition, boost::is_any_of(":"));
        if (twoparts.size() != 2)
            throw XStrom("Expecting exactly one colon in assignment");
        std::string comma_delimited_subset_names_string = twoparts[0];
        std::string comma_delimited_value_string = twoparts[1];
        boost::to_lower(comma_delimited_value_string);
        
        // Check for brackets indicating that parameter should be fixed
        // now see if before_colon contains a data type specification in square brackets
        bool fixed = false;
        const char * pattern_string = R"(\s*\[(.+?)\]\s*)";
#if defined(USE_BOOST_REGEX)
        boost::regex re(pattern_string);
        boost::smatch match_obj;
        bool matched = boost::regex_match(comma_delimited_value_string, match_obj, re);
#else
        std::regex re(pattern_string);
        std::smatch match_obj;
        bool matched = std::regex_match(comma_delimited_value_string, match_obj, re);
#endif
        if (matched) {
            comma_delimited_value_string = match_obj[1];
            fixed = true;
        }
        
        if (comma_delimited_value_string == "equal") {
            vector_of_values.resize(1);
            vector_of_values[0] = -1;
        }
        else {
            // Convert comma_delimited_value_string to vector_of_strings
            std::vector<std::string> vector_of_strings;
            boost::split(vector_of_strings, comma_delimited_value_string, boost::is_any_of(","));

            // Convert vector_of_strings to vector_of_values (single values in case of ratevar, ncateg, and pinvar)
            vector_of_values.resize(vector_of_strings.size());
            std::transform(
                vector_of_strings.begin(),      // start of source vector
                vector_of_strings.end(),        // end of source vector
                vector_of_values.begin(),       // start of destination vector
                [](const std::string & vstr) {  // anonymous function that translates
                    return std::stod(vstr);     // each string element to a double
                }
            );
            assert(vector_of_values.size() > 0);
        }
        
        // Split comma_delimited_subset_names_string into vector_of_subset_names
        boost::split(vector_of_subset_names, comma_delimited_subset_names_string, boost::is_any_of(","));
        
        // Sanity check: at least one subset name must be provided
        if (vector_of_subset_names.size() == 0) {
            XStrom("At least 1 subset must be provided in assignments (use \"default\" if not partitioning)");
        }
        
        // Sanity check: if no partition was defined, then values should be assigned to "default" subset
        // and if "default" is in the list of subset names, it should be the only thing in that list
        unsigned num_subsets_defined = _partition->getNumSubsets();
        std::vector<std::string>::iterator default_iter = std::find(vector_of_subset_names.begin(), vector_of_subset_names.end(), std::string("default"));
        bool default_found = (default_iter != vector_of_subset_names.end());
        if (default_found) {
            if (vector_of_subset_names.size() > 1)
                throw XStrom("The \"default\" specification cannot be mixed with other subset-specific parameter specifications");
        }
        else if (num_subsets_defined == 0) {
            throw XStrom("Must specify partition before assigning values to particular subsets (or assign to subset \"default\")");
        }
        return fixed;
    }

    inline void Strom::sample(unsigned iteration, Chain & chain) {
        if (_nstones > 0) {
            bool time_to_sample = (bool)(iteration % _sample_freq == 0);
            if (time_to_sample && iteration > 0) {
                chain.storeLogLikelihood();
            }

            assert(_heating_powers.size() > 0);
            double largest_power = *(_heating_powers.rbegin());
            if (chain.getHeatingPower() != largest_power)
                return;
                
            bool time_to_report = (bool)(iteration % _print_freq == 0);
            if (time_to_report) {
                double logLike = chain.getLogLikelihood();
                double logPrior = chain.calcLogJointPrior();
                double TL = chain.getTreeManip()->calcTreeLength();
                unsigned m = chain.getTreeManip()->calcResolutionClass();
                if (time_to_report) {
                    if (logPrior == Updater::getLogZero())
                        _output_manager->outputConsole(boost::str(boost::format("%12d %12d %12.5f %12s %12.5f") % iteration % m % logLike % "-infinity" % TL));
                    else
                        _output_manager->outputConsole(boost::str(boost::format("%12d %12d %12.5f %12.5f %12.5f") % iteration % m % logLike % logPrior % TL));
                }
            }
        }
        else {
            if (chain.getHeatingPower() < 1.0)
                return;
                
            bool time_to_sample = (bool)(iteration % _sample_freq == 0);
            bool time_to_report = (bool)(iteration % _print_freq == 0);
            if (time_to_sample || time_to_report) {
                double logLike = chain.getLogLikelihood();
                double logPrior = chain.calcLogJointPrior();
                double TL = chain.getTreeManip()->calcTreeLength();
                unsigned m = chain.getTreeManip()->calcResolutionClass();
                if (time_to_report) {
                    if (logPrior == Updater::getLogZero())
                        _output_manager->outputConsole(boost::str(boost::format("%12d %12d %12.5f %12s %12.5f") % iteration % m % logLike % "-infinity" % TL));
                    else
                        _output_manager->outputConsole(boost::str(boost::format("%12d %12d %12.5f %12.5f %12.5f") % iteration % m % logLike % logPrior % TL));
                }
                if (time_to_sample) {
                    _output_manager->outputTree(iteration, chain.getTreeManip());
                    _output_manager->outputParameters(iteration, logLike, logPrior, TL, m, chain.getModel());
#if defined(HPD)
                    if (_nshells > 0 && iteration > 0) {
                        // Save parameters for marginal likelihood estimation
                        saveLogTransformedParameters(iteration, logLike, logPrior, chain.getModel(), chain.getTreeManip());
                    }
#endif
                }
            }
        }
    }

    inline void Strom::calcHeatingPowers() {
        if (_nstones > 0) {
            // Specify chain heating power for steppingstone
            // For K = 5 chains and alpha = 0.25 (1/alpha = 4):
            //   k   chain power
            // ---------------------
            //   0   (0/5)^4 = 0.0000 <-- prior
            //   1   (1/5)^4 = 0.0016
            //   2   (2/5)^4 = 0.0256
            //   3   (3/5)^4 = 0.1296
            //   4   (4/5)^4 = 0.4096
            //   5   (5/5)^4 = 1.0000 <-- posterior not used
            double inv_alpha = 1.0/_ss_alpha;
            double k = 0.0;
            double K = (double)_heating_powers.size();
            for (auto & h : _heating_powers) {
                h = pow(k++/K, inv_alpha);
            }
        }
        else {
            // Specify chain heating power (e.g. _heating_lambda = 0.2)
            // chain_index  power
            //      0       1.000 = 1/(1 + 0.2*0)
            //      1       0.833 = 1/(1 + 0.2*1)
            //      2       0.714 = 1/(1 + 0.2*2)
            //      3       0.625 = 1/(1 + 0.2*3)
            unsigned i = 0;
            for (auto & h : _heating_powers) {
                h = 1.0/(1.0 + _heating_lambda*i++);
            }
        }
    }

    inline void Strom::showChainTuningInfo() const {
        for (unsigned idx = 0; idx < _nchains; ++idx) {
            for (auto & c : _chains) {
                if (c.getChainIndex() == idx) {
                    _output_manager->outputConsole(boost::str(boost::format("\nChain %d (power %.5f)") % idx % c.getHeatingPower()));
                    std::vector<std::string> names = c.getUpdaterNames();
                    std::vector<double> lambdas    = c.getLambdas();
                    std::vector<double> accepts    = c.getAcceptPercentages();
                    std::vector<unsigned> nupdates = c.getNumUpdates();
                    unsigned n = (unsigned)names.size();
                    _output_manager->outputConsole(boost::str(boost::format("%35s %15s %15s %15s") % "Updater" % "Tuning Param." % "Accept %" % "No. Updates"));
                    for (unsigned i = 0; i < n; ++i) {
                        _output_manager->outputConsole(boost::str(boost::format("%35s %15.3f %15.1f %15d") % names[i] % lambdas[i] % accepts[i] % nupdates[i]));
                    }
                }
            }
        }
    }

    inline void Strom::calcMarginalLikelihood() {
        if (_nstones > 0) {
            // Calculate the log ratio for each steppingstone
            std::vector<std::pair<double, double> > log_ratio;
            for (auto & c : _chains) {
                log_ratio.push_back(std::make_pair(c.getHeatingPower(), c.calcLogSteppingstoneRatio()));
            }
            
            // Sort log_ratio vector from lowest to highest power
            std::sort(log_ratio.begin(), log_ratio.end());
            
            _output_manager->outputConsole("\nSteppingstone results:");
            _output_manager->outputConsole(boost::str(boost::format("%20s %20s %20s") % "beta" % "log(ratio)" % "cumulative"));
            double log_marginal_likelihood = 0.0;
            for (auto p : log_ratio) {
                double beta = p.first;
                double logratio = p.second;
                log_marginal_likelihood += logratio;
                _output_manager->outputConsole(boost::str(boost::format("%20.5f %20.5f %20.5f") % beta % logratio % log_marginal_likelihood));
            }
            _output_manager->outputConsole(boost::str(boost::format("\nlog(marginal likelihood) = %.5f") % log_marginal_likelihood));
        }
#if defined(HPD)
        else if (_nshells > 0) {
            std::cout << "\nEstimating marginal likelihood using HPD Histogram method:" << std::endl;
            if (_skipMCMC) {
                inputStandardizedSamples();
            }
            if (!_skipMCMC) {
                standardizeParameters();
                saveStandardizedSamples();
            }
            defineShells();
            throwDarts();
            histogramMethod();
        }
#endif
    }
    
    inline void Strom::startTuningChains() {
        _swaps.assign(_nchains*_nchains, 0);
        for (auto & c : _chains) {
            c.startTuning();
        }
    } 
    
    inline void Strom::stopTuningChains() {
        _swaps.assign(_nchains*_nchains, 0);
        for (auto & c : _chains) {
            c.stopTuning();
        }
    }
    
    inline void Strom::stepChains(unsigned iteration, bool sampling) {
        for (auto & c : _chains) {
             c.nextStep(iteration);
            if (sampling)
                sample(iteration, c);
        }
    }

    inline void Strom::swapChains() {
        if (_nchains == 1 || _nstones > 0)
            return;
            
        // Select two chains at random to swap
        // If _nchains = 3...
        //  i  j  = (i + 1 + randint(0,1)) % _nchains
        // ---------------------------------------------
        //  0  1  = (0 + 1 +      0      ) %     3
        //     2  = (0 + 1 +      1      ) %     3
        // ---------------------------------------------
        //  1  2  = (1 + 1 +      0      ) %     3
        //     0  = (1 + 1 +      1      ) %     3
        // ---------------------------------------------
        //  2  0  = (2 + 1 +      0      ) %     3
        //     1  = (2 + 1 +      1      ) %     3
        // ---------------------------------------------
        unsigned i = (unsigned)_lot->randint(0, _nchains-1);
        unsigned j = i + 1 + (unsigned)_lot->randint(0, _nchains-2);
        j %= _nchains;

        assert(i != j && i >=0 && i < _nchains && j >= 0 && j < _nchains);

        // Determine upper and lower triangle cells in _swaps vector
        unsigned smaller = _nchains;
        unsigned larger  = _nchains;
        double index_i   = _chains[i].getChainIndex();
        double index_j   = _chains[j].getChainIndex();
        if (index_i < index_j) {
            smaller = index_i;
            larger  = index_j;
        }
        else {
            smaller = index_j;
            larger  = index_i;
        }
        unsigned upper = smaller*_nchains + larger;
        unsigned lower = larger*_nchains  + smaller;
        _swaps[upper]++;

        // Propose swap of chains i and j
        // Proposed state swap will be successful if a uniform random deviate is less than R, where
        //    R = Ri * Rj = (Pi(j) / Pi(i)) * (Pj(i) / Pj(j))
        // Chain i: power = a, kernel = pi
        // Chain j: power = b, kernel = pj
        //      pj^a         pi^b
        // Ri = ----    Rj = ----
        //      pi^a         pj^b
        // log R = (a-b) [log(pj) - log(pi)]

        double heat_i       = _chains[i].getHeatingPower();
        double log_kernel_i = _chains[i].calcLogLikelihood() + _chains[i].calcLogJointPrior();

        double heat_j       = _chains[j].getHeatingPower();
        double log_kernel_j = _chains[j].calcLogLikelihood() + _chains[j].calcLogJointPrior();

        double logR = (heat_i - heat_j)*(log_kernel_j - log_kernel_i);

        double logu = _lot->logUniform();
        if (logu < logR) {
            // accept swap
            _swaps[lower]++;
            _chains[j].setHeatingPower(heat_i);
            _chains[i].setHeatingPower(heat_j);
            _chains[j].setChainIndex(index_i);
            _chains[i].setChainIndex(index_j);
            std::vector<double> lambdas_i = _chains[i].getLambdas();
            std::vector<double> lambdas_j = _chains[j].getLambdas();
            _chains[i].setLambdas(lambdas_j);
            _chains[j].setLambdas(lambdas_i);
        }
    }

    inline void Strom::stopChains() {
        for (auto & c : _chains)
            c.stop();
    }

    inline void Strom::swapSummary() const {
        if (_nchains > 1 && _nstones == 0) {
            unsigned i, j;
            std::cout << "\nSwap summary (upper triangle = no. attempted swaps; lower triangle = no. successful swaps):" << std::endl;

            // column headers
            std::cout << boost::str(boost::format("%12s") % " ");
            for (i = 0; i < _nchains; ++i)
                std::cout << boost::str(boost::format(" %12d") % i);
            std::cout << std::endl;

            // top line
            std::cout << boost::str(boost::format("%12s") % "------------");
            for (i = 0; i < _nchains; ++i)
                std::cout << boost::str(boost::format("-%12s") % "------------");
            std::cout << std::endl;

            // table proper
            for (i = 0; i < _nchains; ++i) {
                std::cout << boost::str(boost::format("%12d") % i);
                for (j = 0; j < _nchains; ++j) {
                    if (i == j)
                        std::cout << boost::str(boost::format(" %12s") % "---");
                    else
                        std::cout << boost::str(boost::format(" %12.5f") % _swaps[i*_nchains + j]);
                }
                std::cout << std::endl;
            }

            // bottom line
            std::cout << boost::str(boost::format("%12s") % "------------");
            for (i = 0; i < _nchains; ++i)
                std::cout << boost::str(boost::format("-%12s") % "------------");
            std::cout << std::endl;
        }
    }

    inline void Strom::initChains() {
        // Create _nchains chains
        _chains.resize(_nchains);
        
        // Create _nchains by _nchains swap matrix
        _swaps.assign(_nchains*_nchains, 0);

        // Create heating power vector
        _heating_powers.assign(_nchains, 1.0);
        calcHeatingPowers();
        
        // Initialize chains
        for (unsigned chain_index = 0; chain_index < _nchains; ++chain_index) {
            auto & c        = _chains[chain_index];
            auto likelihood = _likelihoods[chain_index];
            auto m          = likelihood->getModel();
            
            // Finish setting up models
            m->setTopologyPriorOptions(_allow_polytomies, _resolution_class_prior, _topo_prior_C);
            m->setSubsetNumPatterns(_data->calcNumPatternsVect());
            m->setSubsetSizes(_partition->calcSubsetSizes());
            m->activate();
            if (chain_index == 0)
                std::cout << "\n" << m->describeModel() << std::endl;
            else
                m->describeModel();
                
            // Finish setting up likelihoods
            likelihood->setData(_data);
            likelihood->useUnderflowScaling(_use_underflow_scaling);
            likelihood->initBeagleLib();
            likelihood->useStoredData(_using_stored_data);
            
            // Build list of updaters, one for each free parameter in the model
            unsigned num_free_parameters = c.createUpdaters(m, _lot, likelihood);
            if (num_free_parameters == 0)
                throw XStrom("MCMC skipped because there are no free parameters in the model");

            // Tell the chain that it should adapt its updators (at least initially)
            c.startTuning();

            // Set heating power to precalculated value
            c.setChainIndex(chain_index);
            c.setHeatingPower(_heating_powers[chain_index]);
            if (_nstones > 0) {
                if (chain_index == _nchains - 1)
                    c.setNextHeatingPower(1.0);
                else
                    c.setNextHeatingPower(_heating_powers[chain_index + 1]);
            }
                        
            // Give the chain a starting tree
            std::string newick = _tree_summary->getNewick(m->getTreeIndex());
            c.setTreeFromNewick(newick);

            // Print headers in output files and make sure each updator has its starting value
            c.start();
        }
    }

    inline void Strom::readData() {
        std::cout << "\n*** Reading and storing the data in the file " << _data_file_name << std::endl;
        _data = Data::SharedPtr(new Data());
        _data->setPartition(_partition);
        _data->getDataFromFile(_data_file_name);
    }
    
    inline void Strom::readTrees() {
        assert(_data);
        assert(_likelihoods.size() > 0 && _likelihoods[0]);
        auto m = _likelihoods[0]->getModel();
        unsigned tree_index = m->getTreeIndex();
        std::cout << "\n*** Reading and storing tree number " << (tree_index + 1) << " in the file " << _tree_file_name << std::endl;
        _tree_summary = TreeSummary::SharedPtr(new TreeSummary());
        _tree_summary->readTreefile(_tree_file_name, 0);

        Tree::SharedPtr tree = _tree_summary->getTree(tree_index);
        if (tree->numLeaves() != _data->getNumTaxa())
            throw XStrom(boost::format("Number of taxa in tree (%d) does not equal the number of taxa in the data matrix (%d)") % tree->numLeaves() % _data->getNumTaxa());
    }

    inline void Strom::showPartitionInfo() {
        // Report information about data partition subsets
        unsigned nsubsets = _data->getNumSubsets();
        std::cout << "\nNumber of taxa: " << _data->getNumTaxa() << std::endl;
        std::cout << "Number of partition subsets: " << nsubsets << std::endl;
        for (unsigned subset = 0; subset < nsubsets; subset++) {
            DataType dt = _partition->getDataTypeForSubset(subset);
            std::cout << "  Subset " << (subset+1) << " (" << _data->getSubsetName(subset) << ")" << std::endl;
            std::cout << "    data type: " << dt.getDataTypeAsString() << std::endl;
            std::cout << "    sites:     " << _data->calcSeqLenInSubset(subset) << std::endl;
            std::cout << "    patterns:  " << _data->getNumPatternsInSubset(subset) << std::endl;
            std::cout << "    ambiguity: " << (_ambig_missing || dt.isCodon() ? "treated as missing data (faster)" : "handled appropriately (slower)") << std::endl;
        }
    }

    inline void Strom::showBeagleInfo() {
        assert(_likelihoods.size() > 0 && _likelihoods[0]);
        std::cout << "\n*** BeagleLib " << _likelihoods[0]->beagleLibVersion() << " resources:\n";
        std::cout << "Preferred resource: " << (_use_gpu ? "GPU" : "CPU") << std::endl;
        std::cout << "Available resources:" << std::endl;
        std::cout << _likelihoods[0]->availableResources() << std::endl;
        std::cout << "Resources used:" << std::endl;
        std::cout << _likelihoods[0]->usedResources() << std::endl;
    }
    
    inline void Strom::showMCMCInfo() {
        assert(_likelihoods.size() > 0 && _likelihoods[0]);
        std::cout << "\n*** MCMC analysis beginning..." << std::endl;
        if (_likelihoods[0]->usingStoredData()) {
            unsigned tree_index = _likelihoods[0]->getModel()->getTreeIndex();
            Tree::SharedPtr tree = _tree_summary->getTree(tree_index);
            double lnL = _chains[0].getLogLikelihood();
            std::cout << boost::str(boost::format("Starting log likelihood = %.5f") % lnL) << std::endl;
        }
        else
            std::cout << "Exploring prior" << std::endl;
    
        if (_expected_log_likelihood != 0.0)
            std::cout << boost::str(boost::format("      (expecting %.3f)") % _expected_log_likelihood) << std::endl;
    
        std::cout << "Number of chains is " << _nchains << std::endl;
        std::cout << "Burning in for " << _num_burnin_iter << " iterations." << std::endl;
        std::cout << "Running after burn-in for " << _num_iter << " iterations." << std::endl;
        std::cout << "Sampling every " << _sample_freq << " iterations." << std::endl;
        std::cout << "Sample size is " << (int)(_num_iter/_sample_freq) << std::endl;
                
    }
    
    inline void Strom::run() {
        std::cout << "Starting..." << std::endl;
        std::cout << "Current working directory: " << boost::filesystem::current_path() << std::endl;
        
        try {
            readData();
            readTrees();
            showPartitionInfo();

            // Create a Lot object that generates (pseudo)random numbers
            _lot = Lot::SharedPtr(new Lot);
            _lot->setSeed(_random_seed);

            // Create  Chain objects
            initChains();
            
            showBeagleInfo();
            showMCMCInfo();

            // Create an output manager and open output files
            _output_manager.reset(new OutputManager);
            _output_manager->outputConsole(boost::str(boost::format("\n%12s %12s %12s %12s %12s") % "iteration" % "m" % "logLike" % "logPrior" % "TL"));

#if defined(HPD)
            if (_skipMCMC) {
                calcMarginalLikelihood();
            }
            else {
#endif
                if (_nstones == 0) {
                    _output_manager->openTreeFile("trees.tre", _data);
                    _output_manager->openParameterFile("params.txt", _chains[0].getModel());
                }
                sample(0, _chains[0]);
                
                // Burn-in the chains
                startTuningChains();
                for (unsigned iteration = 1; iteration <= _num_burnin_iter; ++iteration) {
                    stepChains(iteration, false);
                    swapChains();
                }
                stopTuningChains();

#if defined(HPD)
                _log_transformed_parameters.clear();
#endif

                // Sample the chains
                for (unsigned iteration = 1; iteration <= _num_iter; ++iteration) {
                    stepChains(iteration, true);
                    swapChains();
                }
                showChainTuningInfo();
                stopChains();
                
                // Create swap summary
                swapSummary();
                
                // Estimate the marginal likelihood if doing steppingstone
                calcMarginalLikelihood();
                
                // Close output files
                if (_nstones == 0) {
                    _output_manager->closeTreeFile();
                    _output_manager->closeParameterFile();
                }
#if defined(HPD)
            }
#endif
        }
        catch (XStrom & x) {
            std::cerr << "Strom encountered a problem:\n  " << x.what() << std::endl;
        }

        std::cout << "\nFinished!" << std::endl;
    }
    
#if defined(HPD)
    inline void Strom::saveLogTransformedParameters(unsigned iteration, double logLike, double logPrior, Model::SharedPtr model, TreeManip::SharedPtr tm) {
        std::vector<double> params;
        
        // Record log-transformed tree length and log-ratio-transformed edge length proportions
        double log_jacobian = tm->logTransformEdgeLengths(params);
        
        // Record log-transformed parameters
        log_jacobian += model->logTransformParameters(params);
        
        if (_nparams == 0)
            _nparams = (unsigned)params.size();
        assert(_nparams == (unsigned)params.size());
        
        ParameterSample v;
        v._iteration = iteration;
        v._kernel = Kernel(logLike, logPrior, log_jacobian, 0.0);
        v._param_vect = Eigen::Map<Eigen::VectorXd>(params.data(),_nparams);
        _log_transformed_parameters.push_back(v);
    }
    
    // Input standardized parameter samples from file _param_file_name so that marginal
    // likelihood can be recomputed without having to resample
    inline void Strom::inputStandardizedSamples() {
        std::string line;
        double param_value;
        std::ifstream inf(_param_file_name);
        
        // Input number of parameters and number of samples
        inf >> _nparams >> _nsamples;
        
        unsigned expected_nsamples = (unsigned)floor(_num_iter/_sample_freq);
        if (expected_nsamples != _nsamples) {
            throw XStrom(boost::format("Expecting to find %d samples in file \"%s\" but instead found %d.\nDid you modify niter or samplefreq settings since the \"%s\" file was created?") % expected_nsamples % _param_file_name % _nsamples % _param_file_name);
        }

        unsigned i;
        unsigned sz = _nparams*_nparams;
        
        // Input _logDetSqrtS
        inf >> _logDetSqrtS;
        
        // Input variance-covariance matrix _S
        _S.resize(_nparams, _nparams);
        for (i = 0; i < sz; ++i)
            inf >> _S(i);
        
        // Input square root matrix _sqrtS
        _sqrtS.resize(_nparams, _nparams);
        for (i = 0; i < sz; ++i)
            inf >> _sqrtS(i);
        
        // Input inverse square root matrix _invSqrtS
        _invSqrtS.resize(_nparams, _nparams);
        for (i = 0; i < sz; ++i)
            inf >> _invSqrtS(i);
        
        // Input mean vector _mean_transformed
        _mean_transformed.resize(_nparams);
        for (i = 0; i < _nparams; ++i)
            inf >> _mean_transformed(i);
        
        // Input mode vector _mode_transformed
        _mode_transformed.resize(_nparams);
        for (i = 0; i < _nparams; ++i)
            inf >> _mode_transformed(i);
        
        // Input sample vectors
        _standardized_parameters.clear();
        _standardized_parameters.resize(_nsamples);
        i = 0;
        std::getline(inf, line); // swallow up newline after _mode_transformed values input
        while (std::getline(inf, line)) {
            ParameterSample sample;
            std::istringstream iss(line);
            iss >> _standardized_parameters[i]._iteration;
            iss >> _standardized_parameters[i]._kernel._log_likelihood;
            iss >> _standardized_parameters[i]._kernel._log_prior;
            iss >> _standardized_parameters[i]._kernel._log_jacobian_log_transformation;
            iss >> _standardized_parameters[i]._kernel._log_jacobian_standardization;
            assert(i < _nsamples);
            _standardized_parameters[i]._param_vect.resize(_nparams);
            for (unsigned j = 0; j < _nparams; ++j) {
                iss >> param_value;
                _standardized_parameters[i]._param_vect(j) = param_value;
            }
            ++i;
            if (i == _nsamples)
                break;
        }
        inf.close();
    }
    
    inline void Strom::standardizeParameters() {
        std::cout << "  Standardizing parameters..." << std::endl;
        
        // Start off by zeroing mean vector (_mean_transformed), mode vector (_mode_transformed), and variance-covariance matrix (_S)
        assert(_nparams > 0);
        _mean_transformed = eigenVectorXd_t::Zero(_nparams);
        _mode_transformed = eigenVectorXd_t::Zero(_nparams);
        _S.resize(_nparams, _nparams);
        _S = eigenMatrixXd_t::Zero(_nparams, _nparams);
        
        // Calculate mean vector _mean_transformed and mode vector _mode_transformed
        assert(_nsamples == 0);
        _nsamples = (unsigned)_log_transformed_parameters.size();
        assert(_nsamples > 1);
        double best_log_kernel = _log_transformed_parameters[0]._kernel.logKernel();
        _mode_transformed = _log_transformed_parameters[0]._param_vect;
        for (auto & v : _log_transformed_parameters) {
            _mean_transformed += v._param_vect;    // adds v._param_vect elementwise to _mean_transformed
            if (v._kernel.logKernel() > best_log_kernel) {
                best_log_kernel = v._kernel.logKernel();
                _mode_transformed = v._param_vect;
            }
        }
        _mean_transformed /= _nsamples;

        // Calculate variance-covariance matrix _S
        for (auto & v : _log_transformed_parameters) {
            eigenVectorXd_t  x = v._param_vect - _mean_transformed;
            _S += x*x.transpose();
        }
        _S /= _nsamples - 1;
        
        // Can use efficient eigensystem solver because S is positive definite symmetric
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(_S);
        if (solver.info() != Eigen::Success) {
            throw XStrom("Error in the calculation of eigenvectors and eigenvalues of the variance-covariance matrix");
        }
        Eigen::VectorXd L = solver.eigenvalues();
        L = L.array().abs();    // in case some values are "negative zero"
        Eigen::MatrixXd V = solver.eigenvectors();
        
        L = L.array().sqrt();   // component-wise square root
        _sqrtS = V*L.asDiagonal()*V.transpose();
        _invSqrtS = _sqrtS.inverse();
        _logDetSqrtS = log(_sqrtS.determinant());
        
        std::cout << boost::format("\n_logDetSqrtS = %.5f\n") % _logDetSqrtS;
        
        _standardized_parameters.clear();
        for (auto & v : _log_transformed_parameters) {
            ParameterSample s;
            eigenVectorXd_t  x = v._param_vect - _mean_transformed;
            s._iteration = v._iteration;
            s._param_vect = _invSqrtS*x;
            s._kernel = v._kernel;
            s._kernel._log_jacobian_standardization = _logDetSqrtS;
            _standardized_parameters.push_back(s);
        }
        
        // Sort log-transformed and standardized parameter vectors from highest to lowest posterior kernel
        std::sort(_standardized_parameters.begin(), _standardized_parameters.end(), std::greater<ParameterSample>());
    }
    
    // Output standardized parameter samples to a file _param_file_name so that marginal
    // likelihood can be recomputed without having to resample
    inline void Strom::saveStandardizedSamples() const {
        Eigen::IOFormat fmt(Eigen::FullPrecision, Eigen::DontAlignCols,"\t", "\t", "", "", "", "");
        std::ofstream outf(_param_file_name);
        outf << boost::format("%d\t%d\n") % _nparams % _nsamples;
        outf << boost::format("%.9f\n") % _logDetSqrtS;
        outf << _S.format(fmt) << "\n";
        outf << _sqrtS.format(fmt) << "\n";
        outf << _invSqrtS.format(fmt) << "\n";
        outf << _mean_transformed.format(fmt) << "\n";
        outf << _mode_transformed.format(fmt) << "\n";
        for (auto & s : _standardized_parameters) {
            outf << boost::format("%d\t%.9f\t%.9f\t%.9f\t%.9f\t")
                % s._iteration
                % s._kernel._log_likelihood
                % s._kernel._log_prior
                % s._kernel._log_jacobian_log_transformation
                % s._kernel._log_jacobian_standardization;
            outf << s._param_vect.format(fmt) << "\n";
        }
        outf.close();
    }

    inline void Strom::defineShells() {
        std::cout << "  Partitioning working parameter space..." << std::endl;
        
        // Determine how many sample vectors to use for working parameter space
        unsigned nretained = (unsigned)floor(_coverage*_nsamples);
        assert(nretained > 1);
        std::cout << boost::format("    fraction of samples used: %.3f\n") % _coverage;
        std::cout << boost::format("    retaining %d of %d total samples\n") % nretained % _nsamples;

        // Partition the working parameter space into _nshells equal subsets
        // A small easily-understood example:
        //
        // _nsamples = 20  (total sample size)
        // _coverage = 0.8 (i.e. 16 samples will be retained)
        // _nshells  = 4   (working partition comprises 4 subsets)
        // subset_size = 20/4 = 4
        // sample    log-kernel
        //      0    -547.41060 logq[0] A1 k=0 <-- largest kernel value
        //      1    -547.54520         A1
        //      2    -547.54643         A1
        //      3    -547.59144         A1
        //      4    -547.61429 logq[1] A2 k=1
        //      5    -547.61473         A2
        //      6    -547.73307         A2
        //      7    -547.74205         A2
        //      8    -547.85183 logq[2] A3 k=2
        //      9    -547.87130         A3
        //     10    -547.97558         A3
        //     11    -548.00282         A3
        //     12    -548.11777 logq[3] A4 k=3
        //     13    -548.66624         A4
        //     14    -549.01043         A4
        //     15    -549.07118         A4
        //     16    -549.25179 logq[4]    k=4 <-- nretained = 16 (16 % 4 = 0)
        //     17    -549.28753
        //     18    -549.95452
        //     19    -553.81266 <-------------- smallest kernel value
        //
        // What if nretained % _nshells != 0?
        // _nshells  = 4
        // _nsamples = 1111
        // _coverage = 0.95
        // nretained = 1055
        // nretained / _nshells = 1055/4 = 263.75 (263*4 = 1052, which
        // nretained % _nshells = 1055%4 = 3      is 3 less than 1055)
        // solution --> subtract 3 from nretained:
        //     nretained / _nshells = (1055-3)/4 = 263
        //     nretained % _nshells = (1055-3)%4 = 263
        unsigned remainder = nretained % _nshells;
        if (remainder > 0) {
            nretained -= remainder;
            std::cout << boost::format("    reduced number of samples retained by %d in order to divide samples evenly into %d subsets\n") % remainder % _nshells;
        }
        unsigned subset_size = (unsigned)(nretained/_nshells);
        std::cout << boost::format("    working space comprises %d subsets each containing %d samples\n") % _nshells % subset_size;
        
        std::vector<std::string> qvr_norms;
        std::vector<std::string> qvr_logkernels;
        
        // Calculate properties of each of the _nshells partition subsets and store in vector _A
        // _logq[k] holds the log kernel value at the upper (inclusive) boundary of shell with index k
        // -logq[k+1] holds the log kernel value at the lower (non-inclusive) boundary of the shell with index k
        // all samples included have log kernel values greater than _logq[_nshells]
        std::vector<double> logq;
        _A.clear();
        std::cout << boost::format("%12s %12s %12s %12s %12s %12s %12s %12s %12s\n") % "k" % "begin" % "end" % "n" % "qlower" % "qupper" % "qrep" % "min(norm)" % "max(norm)";

        //temporary!
        unsigned ninside = 0;
        unsigned ntotal = (unsigned)_standardized_parameters.size();
        double norm_upper_bound = 0.0;
        
        for (unsigned k = 0; k <= _nshells; ++k) {
            unsigned gk = k*subset_size;

            // Skip the first point in _standardized_parameters, which is the mode
            if (gk == 0)
                gk = 1;

            logq.push_back(_standardized_parameters[gk]._kernel.logKernel());
            if (k > 0) {
                WorkingSpaceSubset A;
                A._lower_logq     = logq[k];
                A._upper_logq     = logq[k-1];
                A._begin          = gk - subset_size;
                A._end            = gk;
                A._representative = logq[k-1] + log(1.0 + exp(logq[k] - logq[k-1])) - log(2.0);
                
                // Skip the first point in _standardized_parameters, which is the mode
                if (A._begin == 0)
                    A._begin = 1;
                    
                // Calculate norms of all points in subset k-1
                A._norms.resize(subset_size);
                unsigned j = 0;
                A._norm_min = -1.0;
                A._norm_max = 0.0;
                for (unsigned i = A._begin; i < A._end; ++i) {
                    eigenVectorXd_t centered = _standardized_parameters[i]._param_vect;
                    double norm = centered.norm();
                    
                    //temporary!
                    if (norm <= norm_upper_bound)
                        ninside++;
                    
                    qvr_norms.push_back(boost::str(boost::format("%g") % norm));
                    qvr_logkernels.push_back(boost::str(boost::format("%g") % _standardized_parameters[i]._kernel.logKernel()));

                    A._norms[j++] = norm;
                    if (A._norm_min < 0.0 || norm < A._norm_min)
                        A._norm_min = norm;
                    if (norm > A._norm_max)
                        A._norm_max = norm;
                }
                
                //temporary!
                if (k == 1) {
                    norm_upper_bound = A._norm_max;
                    ninside = subset_size;
                }
                
                // Calculate volume of ring (base area of shell)
                A._log_base_volume = log(pow(A._norm_max,_nparams) - pow(A._norm_min,_nparams)) + (_nparams/2.0)*log(M_PI) - std::lgamma(_nparams/2.0 + 1.0);
                
                _A.push_back(A);
                
                std::cout << boost::format("%12d %12d %12d %12d %12.5f %12.5f %12.5f %12.5f %12.5f\n") % k % (A._begin + 1) % A._end % (A._end - A._begin) % A._lower_logq % A._upper_logq % A._representative % A._norm_min % A._norm_max;
            }
        }

        std::ofstream plotf("qvr.R");
        plotf << boost::format("# %d = number of points with norm <= %.5f\n") % ninside % norm_upper_bound;
        plotf << boost::format("# %d = number of points with norm <= %.5f and log kernel >= %.5f\n") % subset_size % norm_upper_bound % _A[0]._lower_logq;
        plotf << "r <- c(" << boost::algorithm::join(qvr_norms, ",") << ")\n";
        plotf << "logq <- c(" << boost::algorithm::join(qvr_logkernels, ",") << ")\n";
        plotf << "plot(r, logq, type=\"p\", xlab=\"r\", ylab=\"log kernel\")\n";
        plotf << boost::format("abline(h=%.5f, lty=\"dotted\")\n") % _A[0]._lower_logq;
        plotf << boost::format("abline(v=%.5f, lty=\"dotted\")\n") % _A[0]._norm_max;
        plotf.close();
    }

    inline double Strom::calcLogSum(const std::vector<double> & logx_vect) {
        double max_logx = *(std::max_element(logx_vect.begin(), logx_vect.end()));
        double sum_terms = 0.0;
        for (auto logx : logx_vect) {
            sum_terms += exp(logx - max_logx);
        }
        double logsum = max_logx + log(sum_terms);
        return logsum;
    }

    inline void Strom::throwDarts() {
        std::cout << "  Throwing darts to determine base volume for each shell..." << std::endl;
        
        std::vector<std::string> dart_norms;
        std::vector<std::string> dart_logkernels;
        
        double nparams_inverse = 1.0/_nparams;
        for (unsigned k = 0; k < _nshells; ++k) {
            WorkingSpaceSubset & A = _A[k];
            A._ndarts_thrown = 0;
            A._ndarts_inside = 0;
            double tmp_min_dart_norm = A._norm_max;
            double tmp_max_dart_norm = 0.0;
            
            std::ofstream tmpf(boost::str(boost::format("A%d-darts.txt") % k));
            
            while (A._ndarts_thrown < _ndarts) {
                // Draw standard multivariatenormal deviate z with dimension _nparams
                Eigen::VectorXd z(_nparams);
                for (unsigned i = 0; i < _nparams; ++i) {
                    z(i) = _lot->normal();
                }
                
                // Calculate norm of z
                double znorm = z.norm();
                
                // Draw uniform random deviates based on minimum and maximum radius for this shell
                double rmax = A._norm_max;
                double rmin = A._norm_min;
                //double loc = pow(rmin/rmax, _nparams);
                double loc = 0.0;
                double scale = 1.0 - loc;
                double u = loc + _lot->uniform()*scale;
                
                // Construct dart
                Eigen::VectorXd dart(_nparams);
                dart = z*(rmax*pow(u,nparams_inverse)/znorm);
                
                double dart_norm = dart.norm();
                if (dart_norm < rmin)
                    continue;

                if (dart_norm < tmp_min_dart_norm)
                    tmp_min_dart_norm = dart_norm;
                if (dart_norm > tmp_max_dart_norm)
                    tmp_max_dart_norm = dart_norm;

                // Test if dart falls inside HPD shell
                Kernel kernel = calcLogTransformedKernel(dart);
                double qdart = kernel.logKernel();
                tmpf << boost::format("%12.5f %12.5f") % dart_norm % kernel.logKernel();
                A._ndarts_thrown++;
                if (qdart > A._lower_logq && qdart <= A._upper_logq) {
                    tmpf << "1\n";
                    A._ndarts_inside++;
                }
                else
                    tmpf << "0\n";
                    
                if (k == 0) {
                    dart_norms.push_back(boost::str(boost::format("%g") % dart_norm));
                    dart_logkernels.push_back(boost::str(boost::format("%g") % kernel.logKernel()));
                }
            }
            tmpf.close();
            
            if (k == 0) {
                std::ofstream plotf("darts.R");
                plotf << "r <- c(" << boost::algorithm::join(dart_norms, ",") << ")\n";
                plotf << "logq <- c(" << boost::algorithm::join(dart_logkernels, ",") << ")\n";
                plotf << "plot(r, logq, type=\"p\", xlab=\"r\", ylab=\"log kernel\", xlim=c(0,15), ylim=c(-6845,-6813))\n";
                plotf << boost::format("abline(h=%.5f, lty=\"dotted\")\n") % _A[0]._lower_logq;
                plotf << boost::format("abline(v=%.5f, lty=\"dotted\")\n") % _A[0]._norm_max;
                plotf.close();
            }
        }
    }

    inline double Strom::histogramMethod() {
        std::vector<double> log_numer;
        std::vector<double> log_denom;
        std::cout << "    estimated volumes of working parameter space subsets:\n";
        std::cout << boost::format("%12s %12s %12s %12s %12s\n") % "k" % "log(volume)" % "log(ring vol.)" % "log(fraction)" % "fraction";
        for (unsigned k = 0; k < _nshells; ++k) {
            assert(_A[k]._ndarts_thrown > 0);
            if (_A[k]._ndarts_inside == 0) {
                throw XStrom(boost::format("None of the %d darts were inside shell %d") % _A[k]._ndarts_thrown % k);
            }
            double log_ring_volume = _A[k]._log_base_volume;
            double log_fraction_inside = log(_A[k]._ndarts_inside) - log(_A[k]._ndarts_thrown);
            double log_Vk = log_ring_volume + log_fraction_inside;
            std::cout << boost::format("%12d %12.5f %12.5f %12.5f %12.5f\n") % k % log_Vk % log_ring_volume % log_fraction_inside % exp(log_fraction_inside);
            log_denom.push_back(_A[k]._representative + log_Vk);
            for (unsigned i = _A[k]._begin; i < _A[k]._end; ++i) {
                double log_qratio = _A[k]._representative - _standardized_parameters[i]._kernel.logKernel();
                log_numer.push_back(log_qratio);
            }
        }
        double log_Delta = calcLogSum(log_denom);
        std::cout << boost::format("log(Delta) = %.5f\n") % log_Delta;
        double log_marginal_likelihood = -(calcLogSum(log_numer) - log(_nsamples) - log_Delta);
        _output_manager->outputConsole(boost::str(boost::format("\nlog(marginal likelihood) = %.5f") % log_marginal_likelihood));
        return log_marginal_likelihood;
    }

    inline Kernel Strom::calcLogTransformedKernel(Eigen::VectorXd & standardized_logtransformed) {
        // Get reference to cold chain
        Chain & chain = _chains[0];
        
        // Destandardize parameter vector
        Eigen::VectorXd destandardized = _sqrtS*standardized_logtransformed + _mean_transformed;

        // Set edge lengths
        TreeManip::SharedPtr tm = chain.getTreeManip();
        double TL = exp(destandardized[0]);
        unsigned nedges = tm->countEdges();
        double log_jacobian = tm->setEdgeLengthsFromLogTransformed(destandardized, TL, 1, nedges-1);
        
        // Parameterize model
        Model::SharedPtr model = chain.getModel();
        unsigned nparams = (unsigned)(destandardized.rows() - nedges);
        log_jacobian += model->setParametersFromLogTransformed(destandardized, nedges, nparams);

        tm->selectAllPartials();
        tm->selectAllTMatrices();
        double log_likelihood = chain.calcLogLikelihood();
        double log_prior = chain.calcLogJointPrior();
        return Kernel(log_likelihood, log_prior, log_jacobian, _logDetSqrtS);
    }

#endif  //HPD


}
