The deploy.py script allow analyses performed in the 
Wang et al. paper entitled  "LoRaD: marginal likelihood 
from a single posterior sample" submitted to Systematic 
Biology to be recreated.

This analysis requires the conditional compilation macros
LORAD_VARIABLE_TOPOLOGY and HOLDER_ETAL_PRIOR be defined
in conditionals.hpp when lorad is compiled.

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
constant when data is absent (expecting log marginal likelihood
to be 0.0 in this case. The MCMC run conditions were set to: 
    burnin     = 100000
    niter      = 3000000000
    samplefreq = 1000
The longest run required 3h18m50s on an Intel SkyLake 2.70Ghz processor.

model     seed  |  cov1    lorad1 | cov2     lorad2 |  cov3    lorad3
JC      3954141 | 0.900   0.00198 | 0.500  -0.01406 | 0.100  -0.02943
JCI     3954141 | 0.900  -0.00682 | 0.500  -0.01888 | 0.100  -0.01924
JCG     3954141 | 0.900  -0.01408 | 0.500  -0.01729 | 0.100  -0.03375
JCIG    3954141 | 0.900  -0.38547 | 0.500  -0.01466 | 0.100  -0.02539
GTR     3954141 | 0.900   0.02256 | 0.500  -0.00497 | 0.100  -0.01363
GTRI    3954141 | 0.900   0.07994 | 0.500   0.01471 | 0.100  -0.05501
GTRG    3954141 | 0.900   0.07354 | 0.500   0.02798 | 0.100  -0.01552
GTRIG   3954141 | 0.900   0.03680 | 0.500  -0.01647 | 0.100  -0.04395
3JC     3954141 | 0.900   0.01498 | 0.500  -0.02169 | 0.100  -0.04638
3JCI    3954141 | 0.900   0.00505 | 0.500   0.00040 | 0.100  -0.03416
3JCG    3954141 | 0.900   0.06432 | 0.500  -0.02027 | 0.100   0.00747
3JCIG   3954141 | 0.900   0.06363 | 0.500  -0.00781 | 0.100  -0.03908
3GTR    3954141 | 0.900   0.02418 | 0.500  -0.04076 | 0.100  -0.04527
3GTRI   3954141 | 0.900   0.09459 | 0.500   0.03068 | 0.100  -0.07046
3GTRG   3954141 | 0.900   0.14493 | 0.500   0.06548 | 0.100  -0.12992
3GTRIG  3954141 | 0.900   0.02845 | 0.500  -0.00274 | 0.100   0.04882

Below is summary output from an estimation of the normalizing 
constant when data is present using the default settings:

model      seed       gss |  cov1       lorad1     diff1 |  cov2       lorad2     diff2 |  cov3       lorad3     diff3
JC      3954141  -2776.52 | 0.900  -2776.51041   0.00959 | 0.500  -2776.53583  -0.01583 | 0.100  -2776.51250   0.00750
JCI     3954141  -2744.59 | 0.900  -2744.54688   0.04312 | 0.500  -2744.59773  -0.00773 | 0.100  -2744.57384   0.01616
JCG     3954141  -2747.44 | 0.900  -2747.43269   0.00731 | 0.500  -2747.44287  -0.00287 | 0.100  -2747.41476   0.02524
JCIG    3954141  -2743.56 | 0.900  -2743.50855   0.05145 | 0.500  -2743.57589  -0.01589 | 0.100  -2743.58137  -0.02137
GTR     3954141  -2714.20 | 0.900  -2714.28060  -0.08060 | 0.500  -2714.32212  -0.12212 | 0.100  -2714.29485  -0.09485
GTRI    3954141  -2681.00 | 0.900  -2681.15386  -0.15386 | 0.500  -2681.22450  -0.22450 | 0.100  -2681.21624  -0.21624
GTRG    3954141  -2682.73 | 0.900  -2682.92394  -0.19394 | 0.500  -2682.95596  -0.22596 | 0.100  -2682.93688  -0.20688
GTRIG   3954141  -2680.29 | 0.900  -2680.35500  -0.06500 | 0.500  -2680.38685  -0.09685 | 0.100  -2680.38629  -0.09629
3JC     3954141  -2681.79 | 0.900  -2681.77103   0.01897 | 0.500  -2681.80590  -0.01590 | 0.100  -2681.83674  -0.04674
3JCI    3954141  -2668.38 | 0.900  -2668.33014   0.04986 | 0.500  -2668.36467   0.01533 | 0.100  -2668.36910   0.01090
3JCG    3954141  -2668.99 | 0.900  -2668.93644   0.05356 | 0.500  -2668.96676   0.02324 | 0.100  -2669.02708  -0.03708
3JCIG   3954141  -2667.19 | 0.900  -2667.00391   0.18609 | 0.500  -2667.07251   0.11749 | 0.100  -2667.15399   0.03601
3GTR    3954141  -2551.10 | 0.900  -2551.43303  -0.33303 | 0.500  -2551.44330  -0.34330 | 0.100  -2551.40931  -0.30931
3GTRI   3954141  -2535.57 | 0.900  -2536.15774  -0.58774 | 0.500  -2536.30243  -0.73243 | 0.100  -2536.18738  -0.61738
3GTRG   3954141  -2536.75 | 0.900  -2537.06411  -0.31411 | 0.500  -2537.12071  -0.37071 | 0.100  -2537.21766  -0.46766
3GTRIG  3954141  -2534.66 | 0.900  -2534.97682  -0.31682 | 0.500  -2535.18014  -0.52014 | 0.100  -2535.27264  -0.61264
