#pragma once    

//#include "data.hpp"
//#include "tree_manip.hpp"
//#include "model.hpp"
#include "xlorad.hpp"
#include <fstream>
#include <boost/program_options.hpp>

namespace lorad {

    class OutputManager {
        public:
            typedef std::shared_ptr< OutputManager >            SharedPtr;

                                                                OutputManager();
                                                                ~OutputManager();
            
            //void                                                setModel(Model::SharedPtr model) {_model = model;}
            //void                                                setTreeManip(TreeManip::SharedPtr tm) {_tree_manip = tm;}
            
            void                                                openTreeFile(const std::string & filename, const std::string & taxa_block, const std::string & translate_statement);
            //void                                                openTreeFile(std::string filename, Data::SharedPtr data);
            void                                                openParameterFile(const std::string & filename, const std::string & parameter_names);
            //void                                                openParameterFile(std::string filename, Model::SharedPtr model);
            
            void                                                closeTreeFile();
            void                                                closeParameterFile();

            void                                                outputConsole() const;
            void                                                outputConsole(const std::string & s) const;
            void                                                outputConsole(const boost::format & fmt) const;
            void                                                outputConsole(const boost::program_options::options_description & description) const;
            //void                                                outputTree(unsigned iter, TreeManip::SharedPtr tm);
            void                                                outputTree(unsigned iter, const std::string & newick);
            //void                                                outputParameters(unsigned iter, double lnL, double lnP, double TL, unsigned m, Model::SharedPtr model);
            void                                                outputParameters(unsigned iter, double lnL, double lnP, double TL, unsigned m, const std::string & parameter_values);

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
        _tree_file_name = filename;
        _treefile.open(_tree_file_name.c_str());
        if (!_treefile.is_open())
            throw XLorad(boost::str(boost::format("Could not open tree file \"%s\"") % _tree_file_name));

        _treefile << "#nexus\n\n";
        //_treefile << data->createTaxaBlock() << std::endl;
        _treefile << taxa_block << std::endl;
       
        _treefile << "begin trees;\n";
        //_treefile << data->createTranslateStatement() << std::endl;
        _treefile << translate_statement << std::endl;
    }

    inline void OutputManager::closeTreeFile() {
        assert(_treefile.is_open());
        _treefile << "end;\n";
        _treefile.close();
    }

    //inline void OutputManager::openParameterFile(std::string filename, Model::SharedPtr model) {
    inline void OutputManager::openParameterFile(const std::string & filename, const std::string & parameter_names) {
        //assert(model);
        assert(!_parameterfile.is_open());
        _param_file_name = filename;
        _parameterfile.open(_param_file_name.c_str());
        if (!_parameterfile.is_open())
            throw XLorad(boost::str(boost::format("Could not open parameter file \"%s\"") % _param_file_name));
        //_parameterfile << boost::str(boost::format("%s\t%s\t%s\t%s\t%s\t%s") % "iter" % "lnL" % "lnPr" % "TL" % "m" % model->paramNamesAsString("\t")) << std::endl;
        _parameterfile << boost::str(boost::format("%s\t%s\t%s\t%s\t%s\t%s") % "iter" % "lnL" % "lnPr" % "TL" % "m" % parameter_names) << std::endl;
    }

    inline void OutputManager::closeParameterFile() {
        assert(_parameterfile.is_open());
        _parameterfile.close();
    }

    inline void OutputManager::outputConsole() const {
        std::cout << std::endl;
    }
    
    inline void OutputManager::outputConsole(const std::string & s) const {
        std::cout << s << std::endl;
    }
    
    inline void OutputManager::outputConsole(const boost::format & fmt) const {
        std::cout << boost::str(fmt) << std::endl;
    }
    
    inline void OutputManager::outputConsole(const boost::program_options::options_description & description) const {
        std::cout << description << std::endl;
    }
    
    //inline void OutputManager::outputTree(unsigned iter, TreeManip::SharedPtr tm) {
    inline void OutputManager::outputTree(unsigned iter, const std::string & newick) {
        assert(_treefile.is_open());
        //assert(tm);
        //_treefile << boost::str(boost::format("  tree iter_%d = %s;") % iter % tm->makeNewick(5)) << std::endl;
        _treefile << boost::str(boost::format("  tree iter_%d = %s;") % iter % newick) << std::endl;
    }
    
    //inline void OutputManager::outputParameters(unsigned iter, double lnL, double lnP, double TL, unsigned m, Model::SharedPtr model) {
    inline void OutputManager::outputParameters(unsigned iter, double lnL, double lnP, double TL, unsigned m, const std::string & parameter_values) {
        //assert(model);
        assert(_parameterfile.is_open());
        //_parameterfile << boost::str(boost::format("%d\t%.5f\t%.5f\t%.5f\t%d\t%s") % iter % lnL % lnP % TL % m % model->paramValuesAsString("\t")) << std::endl;
        _parameterfile << boost::str(boost::format("%d\t%.5f\t%.5f\t%.5f\t%d\t%s") % iter % lnL % lnP % TL % m % parameter_values) << std::endl;
    }

}
