#pragma once

// Uncomment the following define to fix bug in topology prior calculation (polytomy prior or
// resolution class prior being used even if allowpolytomies = no. This bug fix will be made
// permanent (and this define removed) after the next commit/push.
#define POL_TOPOPRIOR_BUGFIX

// Uncomment the following define to enable options (nshells, ndarts, coverage, and skipmcmc)
// that allow the use of the HPD variant of the PWK method for marginal likelihood estimation.
// This version of HPD PWK only works on a fixed tree topology.
#define HPD_PWK_METHOD

// Uncomment the following define to enable HPD PWK for variable topologies as well as fixed topologies.
#define HPD_VARIABLE_TOPOLOGY

// Uncomment the following define to include EdgeProportionUpdater updates for both fixed tree topology
// and variable tree topology analyses. The code turned on by this define will be made permanent (and
// this define removed) after the next commit/push.
#define ALWAYS_UPDATE_EDGE_PROPORTIONS

// Uncomment the following define to implement the generalized steppingstone method described in
// Y Fan, R Wu, MH Chen, L Kuo, and PO Lewis. 2011. Choosing among Partition Models in Bayesian
// Phylogenetics. Molecular Biology and Evolution 28(1):523-532.
#define POLGSS

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

