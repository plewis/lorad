#pragma once

#include "conditionals.hpp"

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
#include <boost/math/special_functions/gamma.hpp>

namespace lorad {

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
        double           _norm;
        Eigen::VectorXd  _param_vect;
        Split::treeid_t  _treeID;
        static bool      _sort_by_topology;

        // Define greater-than operator so that a vector of ParameterSample objects can be sorted
        bool operator>(const ParameterSample & other) const {
            if (_sort_by_topology)
                return _treeID > other._treeID;
            else
                return _norm > other._norm;
                //return _kernel.logKernel() > other._kernel.logKernel();
        }
    };
    
    class LoRaD {
        public:
                                                    LoRaD();
                                                    ~LoRaD();

            void                                    clear();
            void                                    processCommandLineOptions(int argc, const char * argv[]);
            void                                    run();
                    
        private:
            bool                                    processAssignmentString(Model::SharedPtr m, const std::string & which, const std::string & definition);
            void                                    handleAssignmentStrings(Model::SharedPtr m, const boost::program_options::variables_map & vm, std::string label, const std::vector<std::string> & definitions, std::string default_definition); 
#if defined(POLGSS)
            bool                                    processReferenceDistribution(Model::SharedPtr m, const std::string & which, const std::string & definition);
            void                                    handleReferenceDistributions(Model::SharedPtr m, const boost::program_options::variables_map & vm, std::string label, const std::vector<std::string> & definitions);
#endif
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

            void                                    saveParameterNames(Model::SharedPtr model, TreeManip::SharedPtr tm);
            void                                    saveLogTransformedParameters(unsigned iteration, double logLike, double logPrior, Model::SharedPtr model, TreeManip::SharedPtr tm);
            void                                    saveStandardizedSamples();
            void                                    inputStandardizedSamples();
            void                                    standardizeParameters();
            void                                    kernelNormPlot();
            Kernel                                  calcLogTransformedKernel(Eigen::VectorXd & x);
            double                                  calcLogSum(const std::vector<double> & logx_vect);
            double                                  loradMethod(double coverage);
            void                                    createNormLogratioPlot(std::string fnprefix, std::vector< std::pair<double, double> > &  norm_logratios) const;

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
            
            std::string                             _fnprefix;

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

#if defined(POLGSS)
            bool                                    _use_gss;
            bool                                    _fixed_tree_topology;
#endif
            
            bool                                    _lorad;
            bool                                    _skipMCMC;
            std::vector<double>                     _coverages;

            unsigned                                _nparams;
            unsigned                                _nsamples;
            std::string                             _param_file_name;
            std::string                             _trimmed_param_file_name;

            typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> eigenMatrixXd_t;
            eigenMatrixXd_t                         _S; // var-cov matrix
            eigenMatrixXd_t                         _sqrtS;
            eigenMatrixXd_t                         _invSqrtS;
            double                                  _logDetSqrtS;
            
            typedef Eigen::VectorXd eigenVectorXd_t;
            eigenVectorXd_t                         _mean_transformed;
            eigenVectorXd_t                         _mode_transformed;
            
            std::vector<std::string>                _param_names;

#if defined(LORAD_VARIABLE_TOPOLOGY)
            unsigned                                _nsamples_total;
            unsigned                                _focal_topol_count;
            std::string                             _focal_newick;
            std::map< Split::treeid_t, unsigned >   _topology_count;
            std::map< Split::treeid_t, unsigned >   _topology_identity;
            std::map< Split::treeid_t, std::string >   _topology_newick;
            unsigned                                _ntopologies;
            std::deque< ParameterSample >           _log_transformed_parameters;
            std::set<Split::treeid_t>               _treeIDset;
#else
            std::vector< ParameterSample >          _log_transformed_parameters;
#endif
            std::vector< ParameterSample >          _standardized_parameters;
    };
    
    inline LoRaD::LoRaD() {
        clear();
    }

    inline LoRaD::~LoRaD() {
    }

    inline void LoRaD::clear() {
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
        _fnprefix                   = "";
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

#if defined(POLGSS)
        _use_gss                    = false;
        _fixed_tree_topology        = false;
#endif

        _skipMCMC                   = false;
        _lorad                      = false;
        //_coverage                   = 0.1;
        _nparams                    = 0;
        _nsamples                   = 0;
        _param_file_name            = "standardized_params.txt";
        _trimmed_param_file_name    = "standardized_params_trimmed.txt";

#if defined(LORAD_VARIABLE_TOPOLOGY)
        _topology_count.clear();
        _topology_identity.clear();
        _topology_newick.clear();
        _treeIDset.clear();
        _ntopologies = 0;
        _nsamples_total = 0;
        _focal_topol_count = 0;
        _focal_newick = "";
        _param_names.clear();
#endif
    }

    inline void LoRaD::processCommandLineOptions(int argc, const char * argv[]) {
        std::vector<std::string> partition_statefreq;
        std::vector<std::string> partition_rmatrix;
        std::vector<std::string> partition_omega;
        std::vector<std::string> partition_ratevar;
        std::vector<std::string> partition_pinvar;
        std::vector<std::string> partition_ncateg;
        std::vector<std::string> partition_subsets;
        std::vector<std::string> partition_relrates;
        std::vector<std::string> partition_tree;
        std::vector<std::string> coverage_values;
#if defined(POLGSS)
        std::vector<std::string> refdist_statefreq;
        std::vector<std::string> refdist_rmatrix;
        std::vector<std::string> refdist_ratevar;
        std::vector<std::string> refdist_edgeprop;
        std::vector<std::string> refdist_treelen;
        std::vector<std::string> refdist_subsetrelrates;
#endif
        boost::program_options::variables_map vm;
        boost::program_options::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("version,v", "show program version")
            ("fnprefix",        boost::program_options::value(&_fnprefix)->default_value(""),   "prefix for output files (empty by default)")
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
            ("nstones", boost::program_options::value(&_nstones)->default_value(0),                "use heated chains to compute marginal likelihood with the steppingstone method using nstones steppingstone ratios")
            ("ssalpha", boost::program_options::value(&_ss_alpha)->default_value(0.25),                "determines how bunched steppingstone chain powers are toward the prior: chain k of K total chains has power (k/K)^{1/ssalpha}")
#if defined(POLGSS)
            ("usegss", boost::program_options::value(&_use_gss)->default_value(false),                   "use generalized steppingstone (note: ignored unless tree topology is fixed)")
            ("statefreqrefdist", boost::program_options::value(&refdist_statefreq), "a string defining parameters for the state frequency Dirichlet reference distribution for one or more data subsets, e.g. 'first,second:492.0,364.3,347.1,525.1'")
            ("exchangerefdist", boost::program_options::value(&refdist_rmatrix), "a string defining parameters for the rmatrix Dirichlet reference distribution for one or more data subsets, e.g. 'first,second:288.0,129.6,310.3,296.8,223.6,224.8'")
            ("ratevarrefdist", boost::program_options::value(&refdist_ratevar), "a string defining parameters for the exchangeability Gamma reference distribution for one or more data subsets, e.g. 'first,second:213.543,0.018'")
            ("edgeproprefdist", boost::program_options::value(&refdist_edgeprop), "a string defining parameters for the edge length proportions Dirichlet reference distribution, e.g. '509.4,569.4,...,184.7' (note: ellipses used to simplify presentation)")
            ("treelenrefdist", boost::program_options::value(&refdist_treelen), "a string defining parameters for the tree length Gamma reference distribution, e.g. '163.900, 0.011'")
            ("relratesrefdist", boost::program_options::value(&refdist_subsetrelrates), "a string defining parameters for the subset relative rates reference distribution, e.g. '0.37,0.13,2.5'")
#endif
            ("lorad", boost::program_options::value(&_lorad)->default_value(false),                   "use LoRaD marginal likelihood method")
            ("coverage",  boost::program_options::value(&coverage_values), "the fraction of samples used to construct the working parameter space (can specify this option more than once to evaluate several coverage values)")
            ("skipmcmc", boost::program_options::value(&_skipMCMC)->default_value(false),                "estimate marginal likelihood using the LoRaD method from parameter vectors previously saved in paramfile (only used if marglike is yes)")
        ;
        boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);

        // Options in config file one directory up take precedence because this file
        // is read first (if it exists)
        try {
            const boost::program_options::parsed_options & parsed = boost::program_options::parse_config_file< char >("../lorad.conf", desc, false);
            boost::program_options::store(parsed, vm);
        }
        catch(boost::program_options::reading_file & x) {
            ::om.outputConsole("Note: higher-level configuration file (../lorad.conf) not found\n");
        }

#if defined(POLGSS)
        // Read in reference distributions (if the file exists)
        try {
            const boost::program_options::parsed_options & parsed = boost::program_options::parse_config_file< char >("refdist.conf", desc, false);
            boost::program_options::store(parsed, vm);
        }
        catch(boost::program_options::reading_file & x) {
            ::om.outputConsole("Note: no reference distribution configuration file (refdist.conf) was found\n");
        }
#endif
        try {
            const boost::program_options::parsed_options & parsed = boost::program_options::parse_config_file< char >("lorad.conf", desc, false);
            boost::program_options::store(parsed, vm);
        }
        catch(boost::program_options::reading_file & x) {
            ::om.outputConsole("Note: configuration file (lorad.conf) not found\n");
        }
        boost::program_options::notify(vm);

        // If user specified --help on command line, output usage summary and quit
        if (vm.count("help") > 0) {
            ::om.outputConsole(desc);
            ::om.outputConsole();
            std::exit(1);
        }

        // If user specified --version on command line, output version and quit
        if (vm.count("version") > 0) {
            ::om.outputConsole(boost::format("This is %s version %d.%d\n") % _program_name % _major_version % _minor_version);
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
            throw XLorad("nchains must be a positive integer greater than 0");

        // Be sure number of stones is greater than or equal to 0
        if (_nstones < 0)
            throw XLorad("nstones must be a positive integer greater than or equal to 0");

        // Can't set nstones > 0 and _lorad at the same time
        if (_nstones > 0 && _lorad) {
            throw XLorad("Cannot specify the steppingstone marginal likelihood method (nstones > 0) and the LoRaD marginal likelihood method at the same time; please choose one or the other");
        }
            
        // Can't specify skipmcmc unless using LoRaD method
        if (_skipMCMC && _lorad) {
            throw XLorad("Cannot specify skipmcmc unless the LoRaD marginal likelihood method is also specified");
        }

        // If number of stones is greater than 0, then set _nchains to that value
        if (_nstones > 0) {
            ::om.outputConsole(boost::format("\nNumber of chains was set to the specified number of stones (%d)\n\n") % _nstones);
            _nchains = (unsigned)_nstones;
        }

        // If user specified --coverage on command line, save coverage value specified in vector _coverages
        if (_lorad) {
            double c = 0.5;
            if (vm.count("coverage") > 0) {
                for (auto s : coverage_values) {
                    try {
                        c = std::stod(s);
                    }
                    catch (const std::invalid_argument & ia) {
                        throw XLorad(boost::format("specified coverage (%s) was not able to be converted to a floating point number)") % s);
                    }
                    if (c == 0.0)
                        throw XLorad("coverage values must be greater than zero");
                    _coverages.push_back(c);
                }
            }
            else {
                ::om.outputConsole(boost::format("\n*** No coverage specified: using %.3f by default ***\n\n") % c);
                _coverages.push_back(c);
            }
        }

        // Be sure heatfactor is between 0 and 1
        if (_heating_lambda <= 0.0 || _heating_lambda > 1.0)
            throw XLorad("heatfactor must be a real number in the interval (0.0,1.0]");
        
        if (!_using_stored_data)
            ::om.outputConsole("\n*** Not using stored data (posterior = prior) ***\n\n");
            
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
#if defined(POLGSS)
            if (_use_gss && _nstones > 0) {
                handleReferenceDistributions(m, vm, "statefreqrefdist", refdist_statefreq);
                handleReferenceDistributions(m, vm, "exchangerefdist",  refdist_rmatrix);
                handleReferenceDistributions(m, vm, "ratevarrefdist",   refdist_ratevar);
                handleReferenceDistributions(m, vm, "edgeproprefdist",  refdist_edgeprop);
                handleReferenceDistributions(m, vm, "treelenrefdist",   refdist_treelen);
                handleReferenceDistributions(m, vm, "relratesrefdist",  refdist_subsetrelrates);
            }
#endif
            _likelihoods.push_back(likelihood);
        }

#if defined(LORAD_VARIABLE_TOPOLOGY)
        // This version allows tree to be variable, so no sanity check required
#else
        // This version requires tree topology to be fixed if carrying out LoRaD marg. like. estim.
        assert(_likelihoods.size() > 0);
        if (_lorad && !_likelihoods[0]->getModel()->isFixedTree()) {
            throw XLorad("Tree topology must be fixed for LoRaD marginal likelihood method");
        }
#endif
    }
    
#if defined(POLGSS)
    inline void LoRaD::handleReferenceDistributions(Model::SharedPtr m, const boost::program_options::variables_map & vm, std::string label, const std::vector<std::string> & definitions) {
        if (vm.count(label) > 0) {
            for (auto s : definitions) {
                bool ok = processReferenceDistribution(m, label, s);
                if (!ok) {
                    throw XLorad(boost::format("Problem processing reference distribution for %s") % label);
                }
            }
        }
    }
    
    inline bool LoRaD::processReferenceDistribution(Model::SharedPtr m, const std::string & which, const std::string & definition) {
        unsigned num_subsets_defined = _partition->getNumSubsets();
        std::vector<std::string> vector_of_subset_names;
        std::vector<double> vector_of_values;
        bool fixed = splitAssignmentString(definition, vector_of_subset_names, vector_of_values);
        if (fixed) {
            throw XLorad("Square brackets found in %s declaration, but square brackets have no meaning in reference distribution specification");
        }
        
        // Assign values to subsets in model
        bool ok = true;
        if (which == "statefreqrefdist") {
            QMatrix::freq_xchg_ptr_t freq_params = std::make_shared<QMatrix::freq_xchg_t>(vector_of_values);
            if (vector_of_subset_names[0] == "default") {
                for (unsigned i = 0; i < num_subsets_defined; i++)
                    m->setSubsetStateFreqRefDistParams(freq_params, i);
            }
            else {
                for (auto s : vector_of_subset_names) {
                    m->setSubsetStateFreqRefDistParams(freq_params, _partition->findSubsetByName(s));
                }
            }
        }
        else if (which == "exchangerefdist") {
            QMatrix::freq_xchg_ptr_t xchg = std::make_shared<QMatrix::freq_xchg_t>(vector_of_values);
            if (vector_of_subset_names[0] == "default") {
                for (unsigned i = 0; i < num_subsets_defined; i++)
                    m->setSubsetExchangeabilitiesRefDistParams(xchg, i);
            }
            else {
                for (auto s : vector_of_subset_names) {
                    m->setSubsetExchangeabilitiesRefDistParams(xchg, _partition->findSubsetByName(s));
                }
            }
        }
//        else if (which == "omegarefdist") {
//            if (vector_of_values.size() > 1)
//                throw XLorad(boost::format("expecting 1 value for omega, found %d values") % vector_of_values.size());
//            QMatrix::omega_ptr_t omega = std::make_shared<QMatrix::omega_t>(vector_of_values[0]);
//            if (vector_of_subset_names[0] == "default") {
//                for (unsigned i = 0; i < num_subsets_defined; i++)
//                    m->setSubsetOmega(omega, i, fixed);
//            }
//            else {
//                for (auto s : vector_of_subset_names) {
//                    m->setSubsetOmega(omega, _partition->findSubsetByName(s), fixed);
//                }
//            }
//        }
//        else if (which == "pinvarrefdist") {
//            if (vector_of_values.size() > 1)
//                throw XLorad(boost::format("expecting 1 value for pinvar, found %d values") % vector_of_values.size());
//            ASRV::pinvar_ptr_t p = std::make_shared<double>(vector_of_values[0]);
//            bool invar_model = (*p > 0);
//            if (vector_of_subset_names[0] == "default") {
//                for (unsigned i = 0; i < num_subsets_defined; i++) {
//                    m->setSubsetIsInvarModel(invar_model, i);
//                    m->setSubsetPinvar(p, i, fixed);
//                }
//            }
//            else {
//                for (auto s : vector_of_subset_names) {
//                    unsigned i = _partition->findSubsetByName(s);
//                    m->setSubsetIsInvarModel(invar_model, i);
//                    m->setSubsetPinvar(p, i, fixed);
//                }
//            }
//        }
        else if (which == "ratevarrefdist") {
            if (vector_of_values.size() != 2)
                throw XLorad(boost::format("expecting 2 parameter values for the rate variance reference distribution, found %d values") % vector_of_values.size());
            ASRV::ratevar_refdist_ptr_t rv = std::make_shared<ASRV::ratevar_refdist_t>(vector_of_values);
            if (vector_of_subset_names[0] == "default") {
                for (unsigned i = 0; i < num_subsets_defined; i++)
                    m->setSubsetRateVarRefDistParams(rv, i);
            }
            else {
                for (auto s : vector_of_subset_names) {
                    m->setSubsetRateVarRefDistParams(rv, _partition->findSubsetByName(s));
                }
            }
        }
        else if (which == "treelenrefdist") {
            if (vector_of_values.size() != 2)
                throw XLorad(boost::format("expecting 2 parameter values for the tree length reference distribution, found %d values") % vector_of_values.size());
            if (vector_of_subset_names[0] == "default") {
                m->setTreeLengthRefDistParams(vector_of_values);
            }
            else {
                throw XLorad("treelenrefdist must be assigned to the default subset");
            }
        }
        else if (which == "edgeproprefdist") {
            // Must put off checking size (i.e. no. edges) until we have a tree
            if (vector_of_subset_names[0] == "default") {
                m->setEdgeProportionsRefDistParams(vector_of_values);
            }
            else {
                throw XLorad("edgeproprefdist must be assigned to the default subset");
            }
        }
        else if (which == "relratesrefdist") {
            if (vector_of_subset_names[0] == "default") {
                m->setSubsetRelRatesRefDistParams(vector_of_values);
            }
            else {
                throw XLorad("relratesrefdist must be assigned to the default subset");
            }
        }
        else {
            ok = false;
        }

        return ok;
    }
#endif
    
    inline void LoRaD::handleAssignmentStrings(Model::SharedPtr m, const boost::program_options::variables_map & vm, std::string label, const std::vector<std::string> & definitions, std::string default_definition) {
        if (vm.count(label) > 0) {
            bool first = true;
            for (auto s : definitions) {
                bool is_default = processAssignmentString(m, label, s);
                if (is_default && !first)
                    throw XLorad(boost::format("default specification must be first %s encountered") % label);
                first = false;
            }
        }
        else {
            processAssignmentString(m, label, default_definition);
        }
    }
    
    inline bool LoRaD::processAssignmentString(Model::SharedPtr m, const std::string & which, const std::string & definition) {
        unsigned num_subsets_defined = _partition->getNumSubsets();
        std::vector<std::string> vector_of_subset_names;
        std::vector<double> vector_of_values;
        bool fixed = splitAssignmentString(definition, vector_of_subset_names, vector_of_values);
        
        if (vector_of_values.size() == 1 && vector_of_values[0] == -1 && !(which == "statefreq" || which == "rmatrix" || which == "relrate"))
            throw XLorad("Keyword equal is only allowed for statefreq, rmatrix, and relrate");

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
                throw XLorad(boost::format("expecting 1 value for omega, found %d values") % vector_of_values.size());
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
                throw XLorad(boost::format("expecting 1 value for pinvar, found %d values") % vector_of_values.size());
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
                throw XLorad(boost::format("expecting 1 value for ratevar, found %d values") % vector_of_values.size());
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
                throw XLorad(boost::format("expecting 1 value for ncateg, found %d values") % vector_of_values.size());
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
                throw XLorad(boost::format("expecting 1 value for tree, found %d values") % vector_of_values.size());
            unsigned tree_index = vector_of_values[0];
            assert(tree_index > 0);
            m->setTreeIndex(tree_index - 1, fixed);
#if defined(POLGSS)
            _fixed_tree_topology = fixed;
#endif
            if (vector_of_subset_names[0] != "default")
                throw XLorad("tree must be assigned to default only");
        }
        else {
            assert(which == "relrate");
            if (vector_of_subset_names[0] != "default")
                throw XLorad("relrate must be assigned to default only");
            m->setSubsetRelRates(vector_of_values, fixed);
        }

        return default_found;
    }

    inline bool LoRaD::splitAssignmentString(const std::string & definition, std::vector<std::string> & vector_of_subset_names, std::vector<double>  & vector_of_values) {
        // Split subset names from definition
        std::vector<std::string> twoparts;
        boost::split(twoparts, definition, boost::is_any_of(":"));
        if (twoparts.size() != 2)
            throw XLorad("Expecting exactly one colon in assignment");
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
            XLorad("At least 1 subset must be provided in assignments (use \"default\" if not partitioning)");
        }
        
        // Sanity check: if no partition was defined, then values should be assigned to "default" subset
        // and if "default" is in the list of subset names, it should be the only thing in that list
        unsigned num_subsets_defined = _partition->getNumSubsets();
        std::vector<std::string>::iterator default_iter = std::find(vector_of_subset_names.begin(), vector_of_subset_names.end(), std::string("default"));
        bool default_found = (default_iter != vector_of_subset_names.end());
        if (default_found) {
            if (vector_of_subset_names.size() > 1)
                throw XLorad("The \"default\" specification cannot be mixed with other subset-specific parameter specifications");
        }
        else if (num_subsets_defined == 0) {
            throw XLorad("Must specify partition before assigning values to particular subsets (or assign to subset \"default\")");
        }
        return fixed;
    }

    inline void LoRaD::sample(unsigned iteration, Chain & chain) {
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
                    std::string newick = chain.getTreeManip()->makeNewick(5);
                    _output_manager->outputTree(iteration, newick);
                    std::string param_values = chain.getModel()->paramValuesAsString("\t");
                    _output_manager->outputParameters(iteration, logLike, logPrior, TL, m, param_values);
#if defined(POLGSS)
                    if (_fixed_tree_topology && _nstones == 0 && _use_gss && iteration > 0) {
                        // Save parameters and edge proportions/TL so that reference distributions
                        // can be computed at the end of a posterior sampling run
                        chain.getModel()->sampleParams();
                        chain.getTreeManip()->sampleTree();
                    }
#endif
                    if (_lorad) {
                        if (iteration == 0)
                            saveParameterNames(chain.getModel(), chain.getTreeManip());
                        else {
                            // Save parameters for marginal likelihood estimation
                            saveLogTransformedParameters(iteration, logLike, logPrior, chain.getModel(), chain.getTreeManip());
                        }
                    }
                }
            }
        }
    }

    inline void LoRaD::calcHeatingPowers() {
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

    inline void LoRaD::showChainTuningInfo() const {
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

    inline void LoRaD::calcMarginalLikelihood() {
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
        else if (_lorad) {
            _output_manager->outputConsole("\nEstimating marginal likelihood using LoRaD method:");
            if (_skipMCMC) {
                inputStandardizedSamples();
            }
            if (!_skipMCMC) {
                standardizeParameters();
                saveStandardizedSamples();
            }
            kernelNormPlot();
            for (auto coverage : _coverages)
                loradMethod(coverage);
        }
    }
    
    inline void LoRaD::startTuningChains() {
        _swaps.assign(_nchains*_nchains, 0);
        for (auto & c : _chains) {
            c.startTuning();
        }
    } 
    
    inline void LoRaD::stopTuningChains() {
        _swaps.assign(_nchains*_nchains, 0);
        for (auto & c : _chains) {
            c.stopTuning();
        }
    }
    
    inline void LoRaD::stepChains(unsigned iteration, bool sampling) {
        for (auto & c : _chains) {
             c.nextStep(iteration);
            if (sampling)
                sample(iteration, c);
        }
    }

    inline void LoRaD::swapChains() {
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

    inline void LoRaD::stopChains() {
        for (auto & c : _chains)
            c.stop();
    }

    inline void LoRaD::swapSummary() const {
        if (_nchains > 1 && _nstones == 0) {
            unsigned i, j;
            ::om.outputConsole("\nSwap summary (upper triangle = no. attempted swaps; lower triangle = no. successful swaps):\n");

            // column headers
            ::om.outputConsole(boost::format("%12s") % " ");
            for (i = 0; i < _nchains; ++i)
                ::om.outputConsole(boost::format(" %12d") % i);
            ::om.outputConsole();

            // top line
            ::om.outputConsole(boost::format("%12s") % "------------");
            for (i = 0; i < _nchains; ++i)
                ::om.outputConsole(boost::format("-%12s") % "------------");
            ::om.outputConsole();

            // table proper
            for (i = 0; i < _nchains; ++i) {
                ::om.outputConsole(boost::format("%12d") % i);
                for (j = 0; j < _nchains; ++j) {
                    if (i == j)
                        ::om.outputConsole(boost::format(" %12s") % "---");
                    else
                        ::om.outputConsole(boost::format(" %12.5f") % _swaps[i*_nchains + j]);
                }
                ::om.outputConsole();
            }

            // bottom line
            ::om.outputConsole(boost::format("%12s") % "------------");
            for (i = 0; i < _nchains; ++i)
                ::om.outputConsole(boost::format("-%12s") % "------------");
            ::om.outputConsole();
        }
    }

    inline void LoRaD::initChains() {
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
                ::om.outputConsole(boost::format("\n%s\n") % m->describeModel());
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
                throw XLorad("MCMC skipped because there are no free parameters in the model");

            // Tell the chain that it should adapt its updators (at least initially)
            c.startTuning();

#if defined(POLGSS)
            // Set steppingstone status:
            //   0: no steppingstone
            //   1: steppingstone (Xie et al. 2011)
            //   2: generalized steppingstone (Fan et al. 2011)
            c.setSteppingstoneMode(_nstones == 0 ? 0 : (_use_gss ? 2 : 1) );
#endif

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

    inline void LoRaD::readData() {
        ::om.outputConsole(boost::format("\n*** Reading and storing the data in the file %s\n") % _data_file_name);
        _data = Data::SharedPtr(new Data());
        _data->setPartition(_partition);
        _data->getDataFromFile(_data_file_name);
    }
    
    inline void LoRaD::readTrees() {
        assert(_data);
        assert(_likelihoods.size() > 0 && _likelihoods[0]);
        auto m = _likelihoods[0]->getModel();
        unsigned tree_index = m->getTreeIndex();
        ::om.outputConsole(boost::format("\n*** Reading and storing tree number %d in the file %s\n") % (tree_index + 1) % _tree_file_name);
        _tree_summary = TreeSummary::SharedPtr(new TreeSummary());
        _tree_summary->readTreefile(_tree_file_name, 0);

        Tree::SharedPtr tree = _tree_summary->getTree(tree_index);
        if (tree->numLeaves() != _data->getNumTaxa())
            throw XLorad(boost::format("Number of taxa in tree (%d) does not equal the number of taxa in the data matrix (%d)") % tree->numLeaves() % _data->getNumTaxa());
    }

    inline void LoRaD::showPartitionInfo() {
        // Report information about data partition subsets
        unsigned nsubsets = _data->getNumSubsets();
        ::om.outputConsole(boost::format("\nNumber of taxa: %d\n") % _data->getNumTaxa());
        ::om.outputConsole(boost::format("Number of partition subsets: %d\n") % nsubsets);
        for (unsigned subset = 0; subset < nsubsets; subset++) {
            DataType dt = _partition->getDataTypeForSubset(subset);
            ::om.outputConsole(boost::format("  Subset %d (%s)\n") % (subset+1) % _data->getSubsetName(subset));
            ::om.outputConsole(boost::format("    data type: %s\n") % dt.getDataTypeAsString());
            ::om.outputConsole(boost::format("    sites:     %d\n") % _data->calcSeqLenInSubset(subset));
            ::om.outputConsole(boost::format("    patterns:  %d\n") % _data->getNumPatternsInSubset(subset));
            ::om.outputConsole(boost::format("    ambiguity: %s\n") % (_ambig_missing || dt.isCodon() ? "treated as missing data (faster)" : "handled appropriately (slower)"));
        }
    }

    inline void LoRaD::showBeagleInfo() {
        assert(_likelihoods.size() > 0 && _likelihoods[0]);
        ::om.outputConsole(boost::format("\n*** BeagleLib %s resources:\n") % _likelihoods[0]->beagleLibVersion());
        ::om.outputConsole(boost::format("Preferred resource: %s\n") % (_use_gpu ? "GPU" : "CPU"));
        ::om.outputConsole("Available resources:\n");
        ::om.outputConsole(boost::format("%s\n") % _likelihoods[0]->availableResources());
        ::om.outputConsole("Resources used:\n");
        ::om.outputConsole(boost::format("%s\n") % _likelihoods[0]->usedResources());
    }
    
    inline void LoRaD::showMCMCInfo() {
        assert(_likelihoods.size() > 0 && _likelihoods[0]);
        ::om.outputConsole("\n*** MCMC analysis beginning...\n");
        if (_likelihoods[0]->usingStoredData()) {
            unsigned tree_index = _likelihoods[0]->getModel()->getTreeIndex();
            Tree::SharedPtr tree = _tree_summary->getTree(tree_index);
            double lnL = _chains[0].getLogLikelihood();
            ::om.outputConsole(boost::format("Starting log likelihood = %.5f\n") % lnL);
            ::om.outputConsole("Starting log joint prior:\n");
            _chains[0].calcLogJointPrior(true);
        }
        else
            ::om.outputConsole("Exploring prior\n");
    
        if (_expected_log_likelihood != 0.0)
            ::om.outputConsole(boost::format("      (expecting %.3f)\n") % _expected_log_likelihood);
    
        ::om.outputConsole(boost::format("Number of chains is %d\n") % _nchains);
        ::om.outputConsole(boost::format("Burning in for %d iterations.\n") % _num_burnin_iter);
        ::om.outputConsole(boost::format("Running after burn-in for %d iterations.\n") % _num_iter);
        ::om.outputConsole(boost::format("Sampling every %d iterations.\n") % _sample_freq);
        ::om.outputConsole(boost::format("Sample size is %d\n") % (int)(_num_iter/_sample_freq));
                
    }
    
    inline void LoRaD::run() {
        ::om.outputConsole("Starting...\n");
        ::om.outputConsole(boost::format("Pseudorandom number seed: %d\n") % _random_seed);
        ::om.outputConsole(boost::format("Current working directory: %s\n") % boost::filesystem::current_path());
        
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

            if (_skipMCMC) {
                calcMarginalLikelihood();
            }
            else {
                _output_manager->outputConsole(boost::str(boost::format("\n%12s %12s %12s %12s %12s") % "iteration" % "m" % "logLike" % "logPrior" % "TL"));
                if (_nstones == 0) {
                    std::string taxa_block = _data->createTaxaBlock();
                    std::string translate_statement = _data->createTranslateStatement();
                    _output_manager->openTreeFile(boost::str(boost::format("%strees.tre") % _fnprefix), taxa_block, translate_statement);
                    std::string param_names = _chains[0].getModel()->paramNamesAsString("\t");
                    _output_manager->openParameterFile(boost::str(boost::format("%sparams.txt") % _fnprefix), param_names);
                }
                sample(0, _chains[0]);
                
                // Burn-in the chains
                startTuningChains();
                for (unsigned iteration = 1; iteration <= _num_burnin_iter; ++iteration) {
                    stepChains(iteration, false);
                    swapChains();
                }
                stopTuningChains();

                _log_transformed_parameters.clear();

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
#if defined(POLGSS)
                    if (_fixed_tree_topology && _use_gss) {
                        std::string s;
                        s += _chains[0].getModel()->saveReferenceDistributions(_partition);
                        s += _chains[0].getTreeManip()->saveReferenceDistributions();
                                                
                        // Append new reference distribution commands in lorad.conf
                        std::ofstream outf("refdist.conf");
                        outf << s;
                        outf.close();
                    }
#endif
                }
            }
        }
        catch (XLorad & x) {
            std::cerr << "LoRaD encountered a problem:\n  " << x.what() << std::endl;
        }

        ::om.outputConsole("\nFinished!\n");
    }
    
    inline void LoRaD::saveParameterNames(Model::SharedPtr model, TreeManip::SharedPtr tm) {
        _param_names.clear();
        tm->saveParamNames(_param_names);
        model->saveParamNames(_param_names);
    }
    
    inline void LoRaD::saveLogTransformedParameters(unsigned iteration, double logLike, double logPrior, Model::SharedPtr model, TreeManip::SharedPtr tm) {
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
#if defined(LORAD_VARIABLE_TOPOLOGY)
        // Record tree topology in the form of a set of Split objects (i.e. the tree ID)
        tm->storeSplits(v._treeID);
        
        // Increment count associated with this tree topology
        _topology_count[v._treeID]++;

        // Assign next number if topology is distinct
        auto tmp = _treeIDset.insert(v._treeID);
        if (tmp.second) {
            // insertion was successful because v._treeID not found in set
            _ntopologies++;
            _topology_identity[v._treeID] = _ntopologies;
            _topology_newick[v._treeID] = tm->makeNewick(5);
        }
        
#endif
        _log_transformed_parameters.push_back(v);
    }
    
    // Input standardized parameter samples from file _param_file_name so that marginal
    // likelihood can be recomputed without having to resample
    inline void LoRaD::inputStandardizedSamples() {
        std::string line;
        double param_value;
        std::string fn = boost::str(boost::format("%s%s") % _fnprefix % _param_file_name);
        std::ifstream inf(fn.c_str());
        
        // Input number of parameters and number of samples
#if defined(LORAD_VARIABLE_TOPOLOGY)
        inf >> _nparams >> _nsamples_total >> _focal_topol_count;
        inf >> _focal_newick;
        _nsamples = _focal_topol_count;
#else
        inf >> _nparams >> _nsamples;
        unsigned expected_nsamples = (unsigned)floor(_num_iter/_sample_freq);
        if (expected_nsamples != _nsamples) {
            throw XLorad(boost::format("Expecting to find %d samples in file \"%s\" but instead found %d.\nDid you modify niter or samplefreq settings since the \"%s\" file was created?") % expected_nsamples % _param_file_name % _nsamples % _param_file_name);
        }
#endif

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
        
        // Input parameter names
        _param_names.clear();
        std::getline(inf, line);
        boost::split(_param_names, line, boost::is_any_of("\t"));

        // Input sample vectors
        _standardized_parameters.clear();
        _standardized_parameters.resize(_nsamples);
        std::getline(inf, line); // swallow up newline after _mode_transformed values input
        double prev_logkernel = 0.0;
        double curr_logkernel = 0.0;
        while (std::getline(inf, line)) {
            std::istringstream iss(line);
            iss >> i;
            assert(i > 0);
            assert(i <= _nsamples);
            iss >> _standardized_parameters[i-1]._iteration;
            iss >> _standardized_parameters[i-1]._kernel._log_likelihood;
            iss >> _standardized_parameters[i-1]._kernel._log_prior;
            iss >> _standardized_parameters[i-1]._kernel._log_jacobian_log_transformation;
            iss >> _standardized_parameters[i-1]._kernel._log_jacobian_standardization;
            iss >> _standardized_parameters[i-1]._norm;
            curr_logkernel = _standardized_parameters[i-1]._kernel.logKernel();
            assert(i == 1 || curr_logkernel <= prev_logkernel);
            prev_logkernel = curr_logkernel;
            _standardized_parameters[i-1]._param_vect.resize(_nparams);
            for (unsigned j = 0; j < _nparams; ++j) {
                iss >> param_value;
                _standardized_parameters[i-1]._param_vect(j) = param_value;
            }
        }
        assert(i == _nsamples);
        inf.close();
    }
    
#if defined(LORAD_VARIABLE_TOPOLOGY)
    bool topolCountCompare(std::pair<Split::treeid_t,unsigned> elem1, std::pair<Split::treeid_t,unsigned> elem2) {
        return elem1.second < elem2.second;
    }
    
    bool parameterSampleCompare(const ParameterSample & p1, const ParameterSample & p2) {
        return p1._treeID < p2._treeID;
    }
#endif

    inline void LoRaD::standardizeParameters() {
        ::om.outputConsole("  Standardizing parameters...\n");
        
#if defined(LORAD_VARIABLE_TOPOLOGY)
        // Sort _log_transformed_parameters by topology
        ParameterSample::_sort_by_topology = true;
        std::sort(_log_transformed_parameters.begin(), _log_transformed_parameters.end(), std::greater<ParameterSample>());
        
        ::om.outputConsole(boost::format("  Sample size is %d\n") % _log_transformed_parameters.size());
        
        // Identify the dominant tree topology
        auto dominant_iter = std::max_element(_topology_count.begin(), _topology_count.end(), topolCountCompare);

        ::om.outputConsole(boost::format("  Most frequently sampled topology was %d and it occurred %d times\n") % _topology_identity[dominant_iter->first] % dominant_iter->second);
        
        _focal_topol_count = dominant_iter->second;
        _focal_newick = _topology_newick[dominant_iter->first];
        
        // Create a dummy ParameterSample having the correct tree ID (the other fields do not matter)
        ParameterSample dummy;
        dummy._treeID = dominant_iter->first;

        // Obtain pair of iterators to first and last element of _log_transformed_parameters
        // having the same tree ID as dummy
        auto iter_pair = std::equal_range(_log_transformed_parameters.begin(), _log_transformed_parameters.end(), dummy, std::greater<ParameterSample>());
        
        ::om.outputConsole(boost::format("  Distance to lower bound = %d\n") % std::distance(_log_transformed_parameters.begin(), iter_pair.first));
        ::om.outputConsole(boost::format("  Distance from lower to upper bound = %d\n") % std::distance(iter_pair.first,iter_pair.second));
        ::om.outputConsole(boost::format("  Distance from upper bound to end = %d\n") % std::distance(iter_pair.second,_log_transformed_parameters.end()));
        
        // Remove all elements before and after the range of elements corresponding to the most frequently
        // sampled topology
        _nsamples_total = (unsigned)_log_transformed_parameters.size();
        _log_transformed_parameters.erase(_log_transformed_parameters.begin(), iter_pair.first);
        _log_transformed_parameters.erase(iter_pair.second, _log_transformed_parameters.end());
        ::om.outputConsole(boost::format("  Length of _log_transformed_parameters after filtering by topology = %d\n") % _log_transformed_parameters.size());
#endif
        
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
            throw XLorad("Error in the calculation of eigenvectors and eigenvalues of the variance-covariance matrix");
        }
        Eigen::VectorXd L = solver.eigenvalues();
        L = L.array().abs();    // in case some values are "negative zero"
        Eigen::MatrixXd V = solver.eigenvectors();
        
        L = L.array().sqrt();   // component-wise square root
        _sqrtS = V*L.asDiagonal()*V.transpose();
        _invSqrtS = _sqrtS.inverse();
        _logDetSqrtS = log(_sqrtS.determinant());
        
        ::om.outputConsole(boost::format("  _logDetSqrtS = %.5f\n") % _logDetSqrtS);
        
        //_parameter_map.clear();
        _standardized_parameters.clear();
        for (auto & v : _log_transformed_parameters) {
            ParameterSample s;
#if defined(MODE_CENTER)
            eigenVectorXd_t  x = v._param_vect - _mode_transformed;
#else
            eigenVectorXd_t  x = v._param_vect - _mean_transformed;
#endif
            s._iteration = v._iteration;
            s._param_vect = _invSqrtS*x;
            s._norm = s._param_vect.norm();
            s._kernel = v._kernel;
            s._kernel._log_jacobian_standardization = _logDetSqrtS;
            _standardized_parameters.push_back(s);
        }
        
        // Sort log-transformed and standardized parameter vectors from highest to lowest norm
#if defined(LORAD_VARIABLE_TOPOLOGY)
        ParameterSample::_sort_by_topology = false;
#endif
        std::sort(_standardized_parameters.begin(), _standardized_parameters.end(), std::greater<ParameterSample>());
    }
    
    // Output standardized parameter samples to a file _param_file_name so that marginal
    // likelihood can be recomputed without having to resample
    inline void LoRaD::saveStandardizedSamples() {
    
        Eigen::IOFormat fmt(Eigen::FullPrecision, Eigen::DontAlignCols,"\t", "\t", "", "", "", "");
        std::string fn = boost::str(boost::format("%s%s") % _fnprefix % _param_file_name);
        std::ofstream outf(fn.c_str());
#if defined(LORAD_VARIABLE_TOPOLOGY)
        outf << boost::format("%d\t%d\t%d\n") % _nparams % _nsamples_total % _focal_topol_count;
        outf << boost::format("%s\n") % _focal_newick;
#else
        outf << boost::format("%d\t%d\n") % _nparams % _nsamples;
#endif
        outf << boost::format("%.9f\n") % _logDetSqrtS;
        outf << _S.format(fmt) << "\n";
        outf << _sqrtS.format(fmt) << "\n";
        outf << _invSqrtS.format(fmt) << "\n";
        outf << _mean_transformed.format(fmt) << "\n";
        outf << _mode_transformed.format(fmt) << "\n";
        unsigned i = 0;
        outf << boost::format("%s\t%s\t%s\t%s\t%s\t%s") % "row" % "iter" % "lnL" % "lnP" % "lnJtrans" % "lnJstd";
        for (auto & s : _param_names)
            outf << boost::format("\t%s") % s;
        outf << "\n";
        for (auto & s : _standardized_parameters) {
            outf << boost::format("%d\t%d\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t")
                % (i+1)
                % s._iteration
                % s._kernel._log_likelihood
                % s._kernel._log_prior
                % s._kernel._log_jacobian_log_transformation
                % s._kernel._log_jacobian_standardization
                % s._norm;
            outf << s._param_vect.format(fmt) << "\n";
            ++i;
        }
        outf.close();
    }

    inline void LoRaD::kernelNormPlot() {
        std::vector<std::string> qvr_all_norms;
        std::vector<std::string> qvr_all_logkernels;
        
        for (auto & p : _standardized_parameters) {
            eigenVectorXd_t centered = p._param_vect;
            qvr_all_norms.push_back(boost::str(boost::format("%g") % centered.norm()));
            qvr_all_logkernels.push_back(boost::str(boost::format("%g") % p._kernel.logKernel()));
        }
        
        std::ofstream plotf(boost::str(boost::format("%sqvr-all.R") % _fnprefix));

        plotf << "cwd = system('cd \"$( dirname \"$0\" )\" && pwd', intern = TRUE)\n";
        plotf << "setwd(cwd)\n";
        plotf << "pdf(\"qvrall.pdf\")\n";

        plotf << "r <- c(" << boost::algorithm::join(qvr_all_norms, ",") << ")\n";
        plotf << "logq <- c(" << boost::algorithm::join(qvr_all_logkernels, ",") << ")\n";
        plotf << "plot(r, logq, type=\"p\", xlab=\"r\", ylab=\"log kernel\")\n";

        plotf << "dev.off()\n";
        plotf.close();
    }
    
    inline double LoRaD::calcLogSum(const std::vector<double> & logx_vect) {
        double max_logx = *(std::max_element(logx_vect.begin(), logx_vect.end()));
        double sum_terms = 0.0;
        for (auto logx : logx_vect) {
            sum_terms += exp(logx - max_logx);
        }
        double logsum = max_logx + log(sum_terms);
        return logsum;
    }
    
    inline void LoRaD::createNormLogratioPlot(std::string fnprefix, std::vector< std::pair<double, double> > & norm_logratios) const {
        std::string fn = boost::str(boost::format("%s.R") % fnprefix);
        std::ofstream outf(fn.c_str());
        outf << "cwd = system('cd \"$( dirname \"$0\" )\" && pwd', intern = TRUE)\n";
        outf << "setwd(cwd)\n";
        outf << "pdf(\"" << fnprefix << ".pdf\")\n";
        unsigned n = (unsigned)norm_logratios.size();
        outf << "x <- c(";
        outf << norm_logratios[0].first;
        for (unsigned i = 1; i < n; i++) {
            outf << "," << norm_logratios[i].first;
        }
        outf << ")\n";
        outf << "y <- c(";
        outf << norm_logratios[0].second;
        for (unsigned i = 1; i < n; i++) {
            outf << "," << norm_logratios[i].second;
        }
        outf << ")\n";
        outf << "plot(x, y, type=\"p\", bty=\"L\", xlab=\"norm\", ylab=\"log-ratio\")\n";
        outf << "dev.off()\n";
        outf.close();
    }

    inline double LoRaD::loradMethod(double coverage) {
        // Determine how many sample vectors to use for working parameter space
        unsigned nretained = (unsigned)floor(coverage*_nsamples);
        assert(nretained > 1);
        
//        // Calculate norms of all points in the retained sample
//        std::vector<double> norms;
//        norms.resize(nretained);
//        double rmin = -1.0;
//        double rmax = 0.0;
//        unsigned j = 0;
//        for (unsigned i = 0; i < nretained; ++i) {
//            eigenVectorXd_t centered = _standardized_parameters[i]._param_vect;
//            double norm = centered.norm();
//            norms[j++] = norm;
//            if (rmin < 0.0 || norm < rmin)
//                rmin = norm;
//            if (norm > rmax)
//                rmax = norm;
//        }
        double norm_max = _standardized_parameters[nretained-1]._norm;

        // Determine Delta, the integral of the multivariate standard normal distribution
        // from 0.0 to norm_max. This makes use of the formula for the cumulative
        // distribution of radial error at the very bottom of p. 11 in Edmundson, HP. 1961.
        // The distribution of radial error and its statistical application in war gaming.
        // Operations Research 9(1):8-21.
        // The incomplete gamma function used by Edmundson is:
        //   f_t(s)/Gamma(s) = [ int_0^t u^{s-1} e^{-u} du ]/Gamma(s)
        // The Boost gamma_p function is defined (using the same notation):
        //      gamma_p(s,t) = [ int_0^t u^{s-1} e^{-u} du ]/Gamma(s)
        // Edmundson used t = r^2/(2*sigma^2) and s = p/2, and the cumulative distribution
        // function is 2*f_t(s)/Gamma(s). Hence, we need 2*gamma_p(p/2, r^2/(2*sigma^2))
        // That said, Edmundson seems to have been incorrect to include that factor of 2,
        // so we use gamma_p(p/2, r^2/(2*sigma^2)) instead
        double s = _nparams/2.0;
        double t = norm_max*norm_max/2.0; // sigma = 1 for standard normal
        double log_Delta = /* log(2.0) + */ log(boost::math::gamma_p(s, t));
        
        // Calculate the sum of ratios in the PWK method, using the multivariate standard normal
        // density as the reference (rather than a constant representative height as in PWK)
        double log_mvnorm_constant = 0.5*log(2.*M_PI)*_nparams;
        std::vector<double> log_ratios(nretained, 0.0);

        // For plotting norm (x-axis) vs. log_ratio (y-axis)
        std::vector< std::pair<double, double> > norm_logratios(nretained);

        for (unsigned i = 0; i < nretained; ++i) {
            double log_kernel = _standardized_parameters[i]._kernel.logKernel();
            double norm = _standardized_parameters[i]._norm;
            double log_reference = -0.5*pow(norm,2.0) - log_mvnorm_constant;
            log_ratios[i] = log_reference - log_kernel;
            norm_logratios[i] = std::make_pair(norm, log_reference - log_kernel);
        }
        
        // Create plot file showing norm on x-axis and logratio on y-axis
        unsigned coverage_percent = (unsigned)(100.0*coverage);
        std::string fnprefix = boost::str(boost::format("%snorm-logratio-%d") % _fnprefix % coverage_percent);
        createNormLogratioPlot(fnprefix, norm_logratios);

        double log_marginal_likelihood = log_Delta - (calcLogSum(log_ratios) - log(_nsamples));

        _output_manager->outputConsole(boost::str(boost::format("\n  Determining working parameter space for coverage = %.3f...") % coverage));
        _output_manager->outputConsole(boost::str(boost::format("    fraction of samples used: %.3f") % coverage));
        _output_manager->outputConsole(boost::str(boost::format("    retaining %d of %d total samples") % nretained % _nsamples));
        _output_manager->outputConsole(boost::str(boost::format("    number of parameters = %d") % _nparams));
        _output_manager->outputConsole(boost::str(boost::format("    log(marginal likelihood) = %.5f") % log_marginal_likelihood));

        return log_marginal_likelihood;
    }

    inline Kernel LoRaD::calcLogTransformedKernel(Eigen::VectorXd & standardized_logtransformed) {
        // Get reference to cold chain
        Chain & chain = _chains[0];
        
        // Destandardize parameter vector
#if defined(MODE_CENTER)
        Eigen::VectorXd destandardized = _sqrtS*standardized_logtransformed + _mode_transformed;
#else
        Eigen::VectorXd destandardized = _sqrtS*standardized_logtransformed + _mean_transformed;
#endif

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
}
