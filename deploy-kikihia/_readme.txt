The deploy.py script allow analyses performed in the 
Wang et al. paper entitled  "LoRaD: marginal likelihood 
from a single posterior sample" submitted to Systematic 
Biology to be recreated.

This analysis requires the conditional compilation macro
HOLDER_ETAL_PRIOR be undefined in conditionals.hpp when 
lorad is compiled.

S1679.nex was downloaded from treebase.org (study ID 1679)

Set include_lorad and include_gss prior to running deploy.py
to control whether lorad or gss analyses are performed.

You will also need to specify different usernames and email
addresses in deploy.py as well as modify the slurm-*.txt 
templates here to specify the correct location of the 
lorad executable.

To generate a directory g1 containing files needed to 

python deploy.py

Issue this command to start all analyses using slurm:

cd g1
. submit-all.sh

Once all runs are finished, the results can be summarized
as follows:

python summary.py

Below is summary output from an estimation of the normalizing 
constant when data is present using the default settings.

Without regression (i.e. as described in Syst. Biol. paper):

From paper:
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

Regression used to improve fit of reference function (not discussed in the Syst. Biol. paper):

partition   seed       secs           gss   seed       secs  cov1       beta01   beta11    beta21        lorad1   cov2       beta02    beta12    beta22        lorad2  cov3       beta03    beta13     beta23        lorad3
--------9 -----6 --------10 -----------13 -----6 --------10 ----5 ----------12 -------8 --------9 -----------13 -----6 -----------12 -------8 --------9 -----------13 ----5 -----------12 -------8 --------10 -----------13
   unpart  12345  16707.887  -10075.33020  12345  16053.237   0.1  10328.48735  2.01060  -0.15443  -10335.89005    0.5  10323.75450   3.36350  -0.25084  -10335.93656   0.9  10312.06102   6.50177   -0.46043  -10335.85604
   bygene  12345  17965.562  -10108.35254  12345  17541.901   0.1  10355.97732  3.24006  -0.20781  -10369.79967    0.5  10346.73197   5.33085  -0.32550  -10369.67971   0.9  10333.84899   8.13281   -0.47729  -10369.63807
  bycodon  12345  17477.906   -9564.99931  12345  17007.578   0.1   9817.56046  2.48359  -0.16911   -9827.90540    0.5   9809.39975   4.52199  -0.29604   -9827.93038   0.9   9790.85406   8.81400   -0.54338   -9827.87750
   byboth  12345  23404.099   -9622.75667  12345  23909.455   0.1   9835.78248  9.01950  -0.41425   -9887.82835    0.5   9827.31964  10.42972  -0.47291   -9887.73106   0.9   9786.92155  16.88018   -0.72987   -9887.41042
