// Setting up the lorad.conf file for different types of marginal likelihood estimation
//
// GSS (Generalized SteppingStone)
//   1. run standard MCMC on posterior first to compute reference distributions
//        saverefdists = yes (to save reference distributions)
//        nstones      = 0   (don't try to calculate GSS estimate during this run)
//   2. run MCMC on series of power posteriors (GSS will be computed at the end of this run)
//        saverefdists = no  (reference distributions already computed)
//        usegss       = yes (calculate GSS estimate at end)
//        nstones      = 10  (or some other number)
//
// LoRaD (Lowest Radial Distance)
//   1. Run standard MCMC on posterior
//        usegss           = no
//        nstones          = 0
//        lorad            = yes
//        useregression    = no  (set to yes if you want to use regression correction)
//        linearregression = no  (can set to yes to use simpler regression correction)
//        coverage         = 0.1 (can specify several coverage values by repeating)
//        coverage         = 0.5
//        coverage         = 0.9
//        skipmcmc         = no  (set to yes if you've already done MCMC once to just read in sampled values)
//
// GHME (Generalized Harmonic Mean Estimator)
//   1. run standard MCMC on posterior first to compute reference distributions
//        saverefdists = yes (to save reference distributions)
//        nstones      = 0   (don't try to calculate GSS estimate during this run)
//   2. run standard MCMC on posterior again to compute reference density for each point
//        saverefdists = no  (reference distributions already computed)
//        nstones      = 0   (don't try to calculate GSS estimate during this run)
//      note: this is indeed wasteful to require two MCMC runs for GHME
//   3. compute GHME using saved parameter file
//        saverefdists = no  (reference distributions already computed)
//        ghme         = yes
//
// Version history:
// 1.0 used for initial submission to Systematic Biology
// 1.1 (7-July-2022) used for revision (added ability to computer GHME)

#include "conditionals.hpp"

#include <limits>
#include <iostream>
#include "lorad.hpp"

using namespace lorad;

// static data member initializations
bool         ParameterSample::_sort_by_topology = false;
std::string  LoRaD::_program_name        = "lorad";
unsigned     LoRaD::_major_version       = 1;
unsigned     LoRaD::_minor_version       = 1;
const double Node::_smallest_edge_length = 1.0e-12;
const double Updater::_log_zero          = -std::numeric_limits<double>::max();
GeneticCode::genetic_code_definitions_t GeneticCode::_definitions = { // codon order is alphabetical: i.e. AAA, AAC, AAG, AAT, ACA, ..., TTT
    {"standard",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"},
    {"vertmito",             "KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"yeastmito",            "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"moldmito",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"invertmito",           "KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"ciliate",              "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLF"},
    {"echinomito",           "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"euplotid",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLF"},
    {"plantplastid",         "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"},
    {"altyeast",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"},
    {"ascidianmito",         "KNKNTTTTGSGSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"altflatwormmito",      "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLF"},
    {"blepharismamacro",     "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YQYSSSS*CWCLFLF"},
    {"chlorophyceanmito",    "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLYSSSS*CWCLFLF"},
    {"trematodemito",        "NNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"scenedesmusmito",      "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLY*SSS*CWCLFLF"},
    {"thraustochytriummito", "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWC*FLF"}
};

OutputManager om;

int main(int argc, const char * argv[]) {

    LoRaD lorad;
    try {
        lorad.processCommandLineOptions(argc, argv);
        lorad.run();
    }
    catch(std::exception & x) {
        std::cerr << "Exception: " << x.what() << std::endl;
        std::cerr << "Aborted." << std::endl;
    }
    catch(...) {
        std::cerr << "Exception of unknown type!\n";
    }

    return 0;
}
