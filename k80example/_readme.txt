The k80lorad.py script recreates the illustrative example.

The variable do_jc determines whether JC or K80 model is used.
The variable do_steppingstone determines whether Generalized Steppingstone is used.
The variable do_plots determines whether the plots used in Figures 1 and 2 of the paper are generated.

Run the program as follows:

python3 k80lorad.py

The files mcmc-sample.txt, mcmc-info.txt, and (if do_steppingstone is True) ss-results.txt
will be generated. If these files are present when the program is run again, it will 
not perform simulation or MCMC or steppingstone, and instead simply reuse the data and
samples already generated. This is useful if you are playing around with LoRaD settings
such as trainingfrac or coverage, but be aware that you need to remove these files if
you wish to run everything from scratch.



