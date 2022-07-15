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
            void                                                openParameterFile(const std::string & filename, const std::string & parameter_names, unsigned nedges);
            
            void                                                closeTreeFile();
            void                                                closeParameterFile();

            void                                                outputConsole() const;
            void                                                outputConsole(const std::string & s) const;
            void                                                outputConsole(const boost::format & fmt) const;
            void                                                outputConsole(const boost::program_options::options_description & description) const;
            void                                                outputTree(unsigned iter, const std::string & newick);
//#if defined(POLGHM)
//            void                                                outputParameters(unsigned iter, double //lnL, double lnP, double lnR, double TL, unsigned m, const std::string & //parameter_values, std::string & edgelen_values);
//#else
            void                                                outputParameters(unsigned iter, double lnL, double lnP, double TL, unsigned m, const std::string & parameter_values, std::string & edgelen_values);
//#endif

        private:

            //TreeManip::SharedPtr                                _tree_manip;
            //Model::SharedPtr                                    _model;
            std::ofstream                                       _treefile;
            std::ofstream                                       _parameterfile;
            std::string                                         _tree_file_name;
            std::string                                         _param_file_name;
    };
    
    
    inline OutputManager::OutputManager() {
        _tree_file_name = "trees.t";
        _param_file_name = "params.p";
    }

    inline OutputManager::~OutputManager() {
    }

    inline void OutputManager::openTreeFile(const std::string & filename, const std::string & taxa_block, const std::string & translate_statement) {
        assert(!_treefile.is_open());
        
        // Create any directories in path that do not already exist
        boost::filesystem::path p(filename);
        boost::filesystem::path pp = p.parent_path();
        if (!pp.empty() && !boost::filesystem::exists(pp)) {
            bool ok = boost::filesystem::create_directories(pp);
            assert(ok);
            outputConsole(boost::format("Created directories that did not exist in \"%s\"\n") % filename);
        }
        
        _tree_file_name = filename;
        _treefile.open(_tree_file_name.c_str());
        if (!_treefile.is_open())
            throw XLorad(boost::str(boost::format("Could not open tree file \"%s\"") % _tree_file_name));

        _treefile << "#nexus\n\n";
        _treefile << taxa_block << std::endl;
       
        _treefile << "begin trees;\n";
        _treefile << translate_statement << std::endl;
    }

    inline void OutputManager::closeTreeFile() {
        assert(_treefile.is_open());
        _treefile << "end;\n";
        _treefile.close();
    }

    inline void OutputManager::openParameterFile(const std::string & filename, const std::string & parameter_names, unsigned nedges) {
        assert(!_parameterfile.is_open());

        // Create any directories in path that do not already exist
        boost::filesystem::path p(filename);
        boost::filesystem::path pp = p.parent_path();
        if (!pp.empty() && !boost::filesystem::exists(pp)) {
            bool ok = boost::filesystem::create_directories(pp);
            assert(ok);
            outputConsole(boost::format("Created directories that did not exist in \"%s\"\n") % filename);
        }
        
        _param_file_name = filename;
        _parameterfile.open(_param_file_name.c_str());
        if (!_parameterfile.is_open())
            throw XLorad(boost::str(boost::format("Could not open parameter file \"%s\"") % _param_file_name));
#if defined(POLGHM)
        _parameterfile << boost::str(boost::format("%s\t%s\t%s\t%s\t%s\t%s\t") % "iter" % "lnL" % "lnPr" % "lnRef" % "TL" % "m");
#else
        _parameterfile << boost::str(boost::format("%s\t%s\t%s\t%s\t%s\t") % "iter" % "lnL" % "lnPr" % "TL" % "m");
#endif
        if (nedges > 0) {
            for (unsigned v = 1; v <= nedges; v++)
                _parameterfile << boost::str(boost::format("v_%d\t") % v);
        }
        _parameterfile << parameter_names << std::endl;
    }

    inline void OutputManager::closeParameterFile() {
#if defined(POLGHM)
        if (_parameterfile.is_open())
            _parameterfile.close();
#else
        assert(_parameterfile.is_open());
        _parameterfile.close();
#endif
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
        assert(_treefile.is_open());
        _treefile << boost::str(boost::format("  tree iter_%d = [&U] %s;") % iter % newick) << std::endl;
    }
    
//#if defined(POLGHM)
//    // Save log of the joint reference distribution along with the log likelihood and log prior so that
//    // GHME can be computed
//    inline void OutputManager::outputParameters(unsigned iter, double lnL, double lnP, double //lnR, double TL, unsigned m, const std::string & parameter_values, std::string & //edgelen_values) {
//        assert(_parameterfile.is_open());
//        if (edgelen_values.length() > 0) {
//            _parameterfile << boost::str(boost::format("%d\t%.5f\t%.5f\t%.5f\t%.5f\t%d\t") % iter % lnL //% lnP % lnR % TL % m);
//            _parameterfile << edgelen_values;
//            _parameterfile << parameter_values << std::endl;
//        } else {
//            _parameterfile << boost::str(boost::format("%d\t%.5f\t%.5f\t%.5f\t%.5f\t%d\t%s") % iter % //lnL % lnP % lnR % TL % m % parameter_values) << std::endl;
//        }
//    }
//#else
    inline void OutputManager::outputParameters(unsigned iter, double lnL, double lnP, double TL, unsigned m, const std::string & parameter_values, std::string & edgelen_values) {
        assert(_parameterfile.is_open());
        if (edgelen_values.length() > 0) {
            _parameterfile << boost::str(boost::format("%d\t%.5f\t%.5f\t%.5f\t%d\t") % iter % lnL % lnP % TL % m);
            _parameterfile << edgelen_values;
            _parameterfile << parameter_values << std::endl;
        } else {
            _parameterfile << boost::str(boost::format("%d\t%.5f\t%.5f\t%.5f\t%d\t%s") % iter % lnL % lnP % TL % m % parameter_values) << std::endl;
        }
    }
//#endif

}
