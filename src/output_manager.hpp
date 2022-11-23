#pragma once    

//#include "data.hpp"
//#include "tree_manip.hpp"
//#include "model.hpp"
#include "xlorad.hpp"
#include <fstream>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace lorad {

    class OutputManager {
        public:
            typedef std::shared_ptr< OutputManager >            SharedPtr;

                                                                OutputManager();
                                                                ~OutputManager();

            void                                                openTreeFile(const std::string & filename, const std::string & taxa_block, const std::string & translate_statement);
            void                                                openDistinctTopologiesFile(const std::string & filename, const std::string & taxa_block, const std::string & translate_statement);
            void                                                openParameterFile(const std::string & filename, const std::string & parameter_names, unsigned nedges, bool incl_refdist);
            void                                                openLogtransformedParameterFile(const std::string & filename, const std::string & parameter_names, unsigned nedges);
            
            void                                                closeTreeFile();
            void                                                closeDistinctTopologiesFile();
            void                                                closeParameterFile();
            void                                                closeLogtransformedParameterFile();

            void                                                outputConsole() const;
            void                                                outputConsole(const std::string & s) const;
            void                                                outputConsole(const boost::format & fmt) const;
            void                                                outputConsole(const boost::program_options::options_description & description) const;
            void                                                outputTree(unsigned iter, const std::string & newick);
            void                                                outputDistinctTopology(unsigned iter, unsigned topol, const std::string & newick);
            void                                                outputParameters(unsigned iter, double logL, double logP, double TL, const std::string & parameter_values, std::string & edgelen_values);
            void                                                outputParametersAlt(unsigned iter, double logL, double logP, double logR, double TL, const std::string & parameter_values, std::string & edgelen_values);
            void                                                outputLogtransformedParameters(unsigned iter, double logL, double logP, double logJ, unsigned topol, double logTL, const std::string & parameter_values, std::string & edgelen_values);

        private:

            std::string                                         _standard_tree_file_name;
            std::ofstream                                       _standard_tree_file;

            std::string                                         _standard_param_file_name;
            std::ofstream                                       _standard_param_file;

            std::string                                         _distinct_topol_file_name;
            std::ofstream                                       _distinct_topol_file;

            std::string                                         _logtransformed_param_file_name;
            std::ofstream                                       _logtransformed_param_file;
    };
    
    
    inline OutputManager::OutputManager() {
        _standard_tree_file_name = "trees.t";
        _standard_param_file_name = "params.p";
    }

    inline OutputManager::~OutputManager() {
    }

    inline void OutputManager::openTreeFile(const std::string & filename, const std::string & taxa_block, const std::string & translate_statement) {
        assert(!_standard_tree_file.is_open());
        
        // Create any directories in path that do not already exist
        boost::filesystem::path p(filename);
        boost::filesystem::path pp = p.parent_path();
        if (!pp.empty() && !boost::filesystem::exists(pp)) {
            bool ok = boost::filesystem::create_directories(pp);
            assert(ok);
            outputConsole(boost::format("Created directories that did not exist in \"%s\"\n") % filename);
        }
        
        _standard_tree_file_name = filename;
        _standard_tree_file.open(_standard_tree_file_name.c_str());
        if (!_standard_tree_file.is_open())
            throw XLorad(boost::str(boost::format("Could not open tree file \"%s\"") % _standard_tree_file_name));

        _standard_tree_file << "#nexus\n\n";
        _standard_tree_file << taxa_block << std::endl;
       
        _standard_tree_file << "begin trees;\n";
        _standard_tree_file << translate_statement << std::endl;
    }

    inline void OutputManager::closeTreeFile() {
        assert(_standard_tree_file.is_open());
        _standard_tree_file << "end;\n";
        _standard_tree_file.close();
    }

    inline void OutputManager::openDistinctTopologiesFile(const std::string & filename, const std::string & taxa_block, const std::string & translate_statement) {
        assert(!_distinct_topol_file.is_open());
        
        // Create any directories in path that do not already exist
        boost::filesystem::path p(filename);
        boost::filesystem::path pp = p.parent_path();
        if (!pp.empty() && !boost::filesystem::exists(pp)) {
            bool ok = boost::filesystem::create_directories(pp);
            assert(ok);
            outputConsole(boost::format("Created directories that did not exist in \"%s\"\n") % filename);
        }
        
        _distinct_topol_file_name = filename;
        _distinct_topol_file.open(_distinct_topol_file_name.c_str());
        if (!_distinct_topol_file.is_open())
            throw XLorad(boost::str(boost::format("Could not open distinct topologies tree file \"%s\"") % _distinct_topol_file_name));

        _distinct_topol_file << "#nexus\n\n";
        _distinct_topol_file << taxa_block << std::endl;
       
        _distinct_topol_file << "begin trees;\n";
        _distinct_topol_file << translate_statement << std::endl;
    }

    inline void OutputManager::closeDistinctTopologiesFile() {
        assert(_distinct_topol_file.is_open());
        _distinct_topol_file << "end;\n";
        _distinct_topol_file.close();
    }

    inline void OutputManager::openParameterFile(const std::string & filename, const std::string & parameter_names, unsigned nedges, bool incl_refdist) {
        assert(!_standard_param_file.is_open());

        // Create any directories in path that do not already exist
        boost::filesystem::path p(filename);
        boost::filesystem::path pp = p.parent_path();
        if (!pp.empty() && !boost::filesystem::exists(pp)) {
            bool ok = boost::filesystem::create_directories(pp);
            assert(ok);
            outputConsole(boost::format("Created directories that did not exist in \"%s\"\n") % filename);
        }
        
        _standard_param_file_name = filename;
        _standard_param_file.open(_standard_param_file_name.c_str());
        if (!_standard_param_file.is_open())
            throw XLorad(boost::str(boost::format("Could not open parameter file \"%s\"") % _standard_param_file_name));
        if (incl_refdist) {
            _standard_param_file << boost::str(boost::format("%s\t%s\t%s\t%s\t%s\t") % "iter" % "logL" % "logP" % "logR" % "TL");
        }
        else {
            _standard_param_file << boost::str(boost::format("%s\t%s\t%s\t%s\t") % "iter" % "logL" % "logP" % "TL");
        }
        if (nedges > 0) {
            for (unsigned v = 1; v <= nedges; v++)
                _standard_param_file << boost::str(boost::format("edgeProp_%d\t") % v);
        }
        _standard_param_file << parameter_names << std::endl;
    }

    inline void OutputManager::closeParameterFile() {
        if (_standard_param_file.is_open())
            _standard_param_file.close();
    }

    inline void OutputManager::openLogtransformedParameterFile(const std::string & filename, const std::string & parameter_names, unsigned nedges) {
        assert(!_logtransformed_param_file.is_open());

        // Create any directories in path that do not already exist
        boost::filesystem::path p(filename);
        boost::filesystem::path pp = p.parent_path();
        if (!pp.empty() && !boost::filesystem::exists(pp)) {
            bool ok = boost::filesystem::create_directories(pp);
            assert(ok);
            outputConsole(boost::format("Created directories that did not exist in \"%s\"\n") % filename);
        }
        
        _logtransformed_param_file_name = filename;
        _logtransformed_param_file.open(_logtransformed_param_file_name.c_str());
        if (!_logtransformed_param_file.is_open())
            throw XLorad(boost::str(boost::format("Could not open parameter file \"%s\"") % _logtransformed_param_file_name));
        _logtransformed_param_file << boost::str(boost::format("%s\t%s\t%s\t%s\t%s\t%s\t") % "iter" % "logL" % "logP" % "logJ" % "topology" % "logTL");
#if defined(HOLDER_ETAL_PRIOR)
        if (nedges > 0) {
            for (unsigned v = 1; v <= nedges; v++)
                _logtransformed_param_file << boost::str(boost::format("logEdgeLength_%d\t") % v);
        }
#else
        if (nedges > 0) {
            for (unsigned v = 2; v <= nedges; v++)
                _logtransformed_param_file << boost::str(boost::format("logEdgeLenProp_%d\t") % v);
        }
#endif
        _logtransformed_param_file << parameter_names << std::endl;
    }

    inline void OutputManager::closeLogtransformedParameterFile() {
        if (_logtransformed_param_file.is_open())
            _logtransformed_param_file.close();
    }

    inline void OutputManager::outputConsole() const {
        std::cout << std::endl;
    }
    
    inline void OutputManager::outputConsole(const std::string & s) const {
        std::cout << s;
    }
    
    inline void OutputManager::outputConsole(const boost::format & fmt) const {
        std::cout << boost::str(fmt);
    }
    
    inline void OutputManager::outputConsole(const boost::program_options::options_description & description) const {
        std::cout << description << std::endl;
    }
    
    inline void OutputManager::outputTree(unsigned iter, const std::string & newick) {
        assert(_standard_tree_file.is_open());
        _standard_tree_file << boost::str(boost::format("  tree iter_%d = [&U] %s;") % iter % newick) << std::endl;
    }
    
    inline void OutputManager::outputDistinctTopology(unsigned iter, unsigned topol, const std::string & newick) {
        assert(_distinct_topol_file.is_open());
        _distinct_topol_file << boost::str(boost::format("  tree iter_%d_topol_%d = [&U] %s;") % iter % topol % newick) << std::endl;
    }
    
    inline void OutputManager::outputParameters(unsigned iter, double logL, double logP, double TL, const std::string & parameter_values, std::string & edgelen_values) {
        assert(_standard_param_file.is_open());
        if (edgelen_values.length() > 0) {
            _standard_param_file << boost::str(boost::format("%d\t%.5f\t%.5f\t%.5f\t") % iter % logL % logP % TL);
            _standard_param_file << edgelen_values;
            _standard_param_file << parameter_values << std::endl;
        } else {
            _standard_param_file << boost::str(boost::format("%d\t%.5f\t%.5f\t%.5f\t%s") % iter % logL % logP % TL % parameter_values) << std::endl;
        }
    }

    inline void OutputManager::outputParametersAlt(unsigned iter, double logL, double logP, double logR, double TL, const std::string & parameter_values, std::string & edgelen_values) {
        assert(_standard_param_file.is_open());
        if (edgelen_values.length() > 0) {
            _standard_param_file << boost::str(boost::format("%d\t%.5f\t%.5f\t%.5f\t%.5f\t") % iter % logL % logP % logR % TL);
            _standard_param_file << edgelen_values;
            _standard_param_file << parameter_values << std::endl;
        } else {
            _standard_param_file << boost::str(boost::format("%d\t%.5f\t%.5f\t%.5f\t%.5f\t%s") % iter % logL % logP % logR % TL % parameter_values) << std::endl;
        }
    }

    inline void OutputManager::outputLogtransformedParameters(unsigned iter, double logL, double logP, double logJ, unsigned topol, double logTL, const std::string & parameter_values, std::string & edgelen_values) {
        // First save the parameters
        assert(_logtransformed_param_file.is_open());
        if (edgelen_values.length() > 0) {
            _logtransformed_param_file << boost::str(boost::format("%d\t%.5f\t%.5f\t%.5f\t%d\t%.5f\t") % iter % logL % logP % logJ % topol % logTL);
            _logtransformed_param_file << edgelen_values;
            _logtransformed_param_file << parameter_values << std::endl;
        } else {
            _logtransformed_param_file << boost::str(boost::format("%d\t%.5f\t%.5f\t%.5f\t%d\t%.5f\t%s") % iter % logL % logP % logJ % topol % logTL % parameter_values) << std::endl;
        }
    }

}
