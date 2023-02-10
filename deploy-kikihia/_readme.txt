The deploy.py script allow analyses reported in Table 1 and Figure 3
in the following paper to be recreated:

Wang, YB, A Milkey, A Li, MH Chen, L Kuo, and PO Lewis. 
LoRaD: marginal likelihood estimation with haste (but no waste).
Systematic Biology (in press).

The software at https://github.com/plewis/lorad is required to carry
out the MCMC analyses reported in the paper. Separate software at
https://github.com/plewis/loradML is needed for estimating the
marginal likelihood using LoRaD (and GHM if desired).

This analysis requires the conditional compilation macro
HOLDER_ETAL_PRIOR be undefined in conditionals.hpp when 
the lorad (MCMC) software is compiled. Example meson build scripts are
provided in the src directory. The conditional compilation macro
GHM must be defined in conditionals.hpp when the loradML software
is compiled if you want it to estimate the marginal likelihood 
using GHM.

S1679.nex was downloaded from treebase.org (study ID 1679)

Set include_lorad, include_gss, and include_ghm prior to running deploy.py
to control whether lorad, gss, or ghm analyses are performed. You should not
set include_gss or include_ghm to True without also setting include_lorad to
True because GSS and GHM analyses depend on LoRaD having been run first.

You can also set include_rev to True to create Rev scripts for MCMC or 
steppingstone if you prefer to use RevBayes (https://revbayes.github.io) 
to do the analyses (but note that RevBayes analyses were not reported in
the above paper).

Set include_unpart, include_bycodon, include_bygene, and include_byboth, 
as desired, to include analyses of unpartitioned, partitioned by codon position, 
partitioned by gene, and partitioned by both gene and codon position.

You will also need to specify different usernames and email
addresses in deploy.py as well as modify the slurm-*.txt 
templates here to specify the correct location of the 
lorad and loradML executables if they were not placed in a
directory (e.g. /usr/local/bin) that is included in the PATH
environmental variable.

To generate a directory g1 containing files needed for one replicate: 

python3 deploy.py

Issue this command to start LoRaD analyses using slurm:

cd g1
. submit-lorad.sh

The script epideploy.sh is provided to make it easy to start multiple
replicates (20 replicates were used in the paper). The file go.sh is
copied by epideploy.sh into the destination directory to make it easier
to start all runs on a cluster that uses slurm. Parts of go.sh can be
uncommented or (perhaps safer) copied into a separate bash script to
summarize the LoRaD or GSS results using loradML or egrep, respectively.

The file coverage-series.sh can be used to start a slurm array job that
runs loradML for 11 values of the coverage fraction (as was done in the
paper). You will need to change "unpart" to "bycodon", "bygene", or "byboth"
in order to start this analyses for partition schemes other than 
bycodon.




