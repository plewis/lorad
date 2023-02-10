The deploy.py script allow analyses reported in Table 2
in the following paper to be recreated:

Wang, YB, A Milkey, A Li, MH Chen, L Kuo, and PO Lewis. 
LoRaD: marginal likelihood estimation with haste (but no waste).
Systematic Biology (in press)

This analysis requires the software from the repository 
https://github.com/plewis/lorad. The conditional 
compilation macro HOLDER_ETAL_PRIOR must be defined in 
conditionals.hpp when the software is compiled.

You will also need to specify different usernames and email
addresses in deploy.py as well as modify the slurm-*.txt 
templates here to specify the correct location of the 
lorad executable.

To generate a directory g1 containing files needed to perform
one replicate:

python3 deploy.py

Issue this command to start analyses of all 16 models using slurm:

cd g1
. submit-all.sh

Once all runs are finished, the results can be summarized
as follows:

python3 summary.py
