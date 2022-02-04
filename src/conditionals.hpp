#pragma once

// Uncomment the following define to enable LoRaD for variable topologies as well as fixed topologies.
//#define LORAD_VARIABLE_TOPOLOGY

// Uncomment the following define to implement the generalized steppingstone method described in
// Y Fan, R Wu, MH Chen, L Kuo, and PO Lewis. 2011. Choosing among Partition Models in Bayesian
// Phylogenetics. Molecular Biology and Evolution 28(1):523-532.
#define POLGSS

// Uncomment the following to use edge length and rate heterogeneity priors described in this chapter:
// MT Holder, PO Lewis, DL Swofford, and D Bryant. 2014. Variable tree topology stepping-stone marginal
// likelihood estimation. Chapter 5, pp. 95-111, in: MH Chen, L Kuo, and PO Lewis (eds.), Bayesian
// phylogenetics: methods, algorithms, and applications. Chapman and Hall/CRC Mathematical and Computational
// Biology Series, 365 pages. ISBN-13: 978-1-4665-00790-2
//#define HOLDER_ETAL_PRIOR

//
// Conditional defines below here will probably be eliminated soon because they are no longer useful
//

// Uncomment the following define to track the lengths of runs of the focal topology and exclude the
// first focalburnin iterations from each run. The focal topology is always the topology supplied
// via the treefile option (i.e. start the run with the topology you are interested in tracking.
// The focalburnin option does not seem to have much of an effect.
//#define TRACK_RUNS_OF_FOCAL_TREE

// Uncomment the following define to standardize the log-transformed parameter vector using the mode
// rather than the mean. The mode used is the highest point on the log-transformed posterior kernel surface.
//#define MODE_CENTER

// Uncomment the following define only if std::regex is not available
//#define USE_BOOST_REGEX

