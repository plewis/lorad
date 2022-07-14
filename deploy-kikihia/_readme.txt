The deploy.py script allow analyses performed in the 
Wang et al. paper entitled  "LoRaD: accurate marginal likelihood 
with haste (but no waste)" submitted to Systematic 
Biology to be recreated.

This analysis requires the conditional compilation macro
HOLDER_ETAL_PRIOR be undefined in conditionals.hpp when 
lorad is compiled.

S1679.nex was downloaded from treebase.org (study ID 1679)

Set include_lorad, include_gss, and include_ghme prior to running deploy.py
to control whether lorad, gss, or ghme analyses are performed.

You will also need to specify different usernames and email
addresses in deploy.py as well as modify the slurm-*.txt 
templates here to specify the correct location of the 
lorad executable.

To generate a directory g1 containing files needed: 

python3 deploy.py

Issue this command to start all analyses using slurm:

cd g1
. submit-all.sh

Once all runs are finished, the results can be summarized
as follows:

python3 summary.py

Below is summary output from an estimation of the normalizing 
constant when data is present using the default settings.

From the paper:

 Partition Scheme  parameters            SS                      LoRaD
    Unpartitioned          70  -10334.84064 (0.04541)   -10335.89042 (0.05034)
          By gene         100  -10368.37473 (0.07783)   -10367.14373 (0.07633)
         By codon          90   -9826.67881 (0.06662)    -9826.84305 (0.22573)
By gene and codon         180   -9882.77806 (9.10763)    -9884.25645 (0.30798)

partition   seed       secs           gss   seed       secs | cov1        lorad1 | cov2        lorad2 | cov3        lorad3
   unpart  12345  16707.887  -10075.33020  12345  16053.237 |  0.1  -10335.88990 |  0.5  -10335.94243 |  0.9  -10335.89228
   bygene  12345  17965.562  -10108.35254  12345  17541.901 |  0.1  -10369.80493 |  0.5  -10369.69795 |  0.9  -10369.65036
  bycodon  12345  17477.906   -9564.99931  12345  17007.578 |  0.1   -9827.90786 |  0.5   -9827.92508 |  0.9   -9827.94793
   byboth  12345  23404.099   -9622.75667  12345  23909.455 |  0.1   -9887.84049 |  0.5   -9887.74592 |  0.9   -9887.71839

