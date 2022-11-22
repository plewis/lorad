import plotly.express as px
import plotly.figure_factory as ff
import plotly.graph_objs as go
import pandas as pd
import random, sys, os, re
import math
import numpy
import scipy

# If True, Jukes-Cantor (1969) model will be evaluated rather than the Kimura (1980) model
do_jc = False

# True estimates marginal likelihood using generalized steppingstone method; False skips steppingstone analysis
do_steppingstone = False

# True creates figures (slow), False skips figures (fast)
do_plots = False

# LoRaD marginal likelihood estimation method
trainingfrac = 0.5 # fraction of sample to be used for training (i.e. working parameter space determination)
coverage     = 0.1 # fraction of training sample used to define limits of working parameter space

# The variable script_dir holds the directory in which this python script is located
script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))

# True creates files out of figures, False shows figures in web browser
plots_saved_to_file = True

# The height of all figures in pixels assuming plotting_scale = 1.
plotting_height = 1000

# Scaling factor that multiplies the nominal height
plotting_scale = 1

# Determines how fine the grid will be for surface plotting
# (the higher this value, the smoother the surface will look).
plotting_increments = 500

# If True, a file named mcmc-samples.txt will be created in which the MCMC sample is saved
# If it exists, the mcmc-samples.txt file will be read in and the samples therein used.
save_samples_to_file = True

# If True, the Edmundson (1961) cumulative radial error distribution formula will be tested
# using test_sample_size draws from a standard multivariate normal distribution
test_sample_size = 10000
test_edmundson = False

# Pseudorandom number seed
rnseed = 13579

# Levels of verbosity
do_show_sequences                 = False  # show the simulated sequences
do_show_transformed_mean_vector   = False  # show mean vector of log-transformed sample (used for standardization transformation)
do_show_transformed_varcov_matrix = False  # show variance-covariance matrix of log-transformed sample (used for standardization transformation)
do_show_sorted_transformed        = False  # caution: output will have number of lines equal to the number of MCMC samples)

# Parameters of data simulation (data comprises two sequences separated by an edge)
true_edgelen =  0.2  # expected number of substitutions per site in K80 model
true_kappa   =  5.0  # transition/transversion rate ratio in K80 model
seqlen       =  200  # length of simulated sequences

# Axis limits for plots on linear scale
maxv_linear_plot = 0.4
maxk_linear_plot = 10.0

# Data (these will be set by function simulateData)
nsame = 0  # number of sites that are identical
ntrs  = 0  # number of sites that show a transition type substitution (A <-> G, C <-> T)
ntrv  = 0  # number of sites that show a transversion type substitution (A <-> C, A <-> T, C <-> G, C <-> T, G <-> C, G <-> C, T <-> A, T <-> G)

# Parameters
edgelen = true_edgelen
kappa   = true_kappa

if do_jc:
    assumed_kappa = 1.0     # Jukes-Cantor model assumes kappa = 1
    kappa = assumed_kappa

# Both edge length and kappa have Gamma priors with mean equal to prior_shape*prior_scale
prior_shape = 1.0
prior_scale = 50.0

# MCMC settings
burnin            =          1000   # number of burnin iterations
niters            =       1000000   # number of MCMC iterations
thinby            =           100   # save sample every thinby iterations
target_acceptance =           0.3   # deltak and deltav will be tuned to make acceptance rate close to this
deltak            =          50.0   # width of window centered over kappa; used for proposing new value
deltav            =           2.0   # width of window centered over edgelen; used for proposing new value
vaccepts          =             0   # number of times a proposed new edge length was accepted
vupdates          =             0   # number of times edge length was updated
kaccepts          =             0   # number of times a proposed new kappa was accepted
kupdates          =             0   # number of times kappa was updated

# Steppingstone (SS) settings
ssalpha          =         1.0  # shape of the Beta(ssalpha,1) distribution used to choose beta values
nstones          =           5  # number of stones (i.e. ratios) used in steppingstone analysis (used only if do_steppingstone is True)
niters_per_stone =      100000  # number of iterations to use for estimating each ratio
thin_per_stone   =         100  # sample saved every thin_per_stone iterations
burnin_per_stone =        1000  # number of burn-in iterations to use for each stone (used to modify tuning parameters for updaters)
ss_v_alpha       = prior_shape  # shape parameter for the edge length reference distribution
ss_v_beta        = prior_scale  # scale parameter for the edge length reference distribution
ss_kappa_alpha   = prior_shape  # shape parameter for the kappa reference distribution
ss_kappa_beta    = prior_scale  # scale parameter for the kappa reference distribution

# Settings used for both MCMC and SS
starting_edgelen  =  true_edgelen
starting_kappa    =    true_kappa
if do_jc:
    starting_kappa = assumed_kappa
    
# LoRaD related    
lop_off      = 1.0 - coverage   # fraction of sample discarded
mode_center  = False            # if True, use mode vector for standardization; otherwise, use mean vector

# Computes log of the prior for supplied edge length
def logPriorV(x):
    return (prior_shape - 1.0)*math.log(x) - x/prior_scale - math.log(prior_scale)*prior_shape - math.lgamma(prior_shape)

# Computes log of the prior for supplied kappa value
def logPriorK(x):
    return (prior_shape - 1.0)*math.log(x) - x/prior_scale - math.log(prior_scale)*prior_shape - math.lgamma(prior_shape)

# Computes log of the reference distribution for edge lengths used in generalized steppingstone method
def logRefDistV(v):
    return (ss_v_alpha - 1.0)*math.log(v) - v/ss_v_beta - math.log(ss_v_beta)*ss_v_alpha - math.lgamma(ss_v_alpha);

# Computes log of the reference distribution for edge lengths used in generalized steppingstone method
def logRefDistK(k):
    return (ss_kappa_alpha - 1.0)*math.log(k) - k/ss_kappa_beta - math.log(ss_kappa_beta)*ss_kappa_alpha - math.lgamma(ss_kappa_alpha);

# Computes the transition probabilities representing no change (same), 
# a transition-type substitution (trs), or a transversion-type substitution (trv)
# These are used both for simulating data as well as computing the likelihood
def transitionProbs(v, k):
    # exp1 = math.exp(-v/(k+2.0))
    # exp2 = math.exp(-(v/2.0)*(k+1.0)/(k+2.0))
    # trs = 0.25 + 0.25*exp1 - 0.5*exp2
    # trv = 0.25 - 0.25*exp1
    # same = 0.25 + 0.25*exp1 + 0.5*exp2
    trv  = 0.25 - 0.25*math.exp(-4.*v/3)
    trs  = 0.25 + 0.25*math.exp(-4.*v/3) - 0.5*math.exp((-4.*v/3)*0.5*(k+1.))
    same = 0.25 + 0.25*math.exp(-4.*v/3) + 0.5*math.exp((-4.*v/3)*0.5*(k+1.))
    return (same, trs, trv)
    
# Computes the log likelihood of the parameters v (edgelen) and k (kappa) given the data (nsame, ntrs, ntrv) 
def logLikelihood(nsame, ntrs, ntrv, v, k):
    same, trs, trv = transitionProbs(v, k)
    try:
        loglike = math.log(0.25)*(nsame + ntrs + ntrv) + math.log(same)*nsame + math.log(trs)*ntrs + math.log(trv)*ntrv
    except ValueError:
        sys.exit('math domain error in logLikelihood(%d, %d, %d, %.5f, %.5f)' % (nsame, ntrs, ntrv, v, k))
    return loglike
    
# Used to choose a state ('A', 'C', 'G', 'T') for the "root" node given the K80 model, 
# which assumes equal state frequencies at the root. The root in this case is one of the 
# two nodes (chosen arbitrarily) at either end of the single edge connecting the two sequences
def chooseRandomState():
    u = random.random()
    s = None
    if u < 0.25:
        s = 'A'
    elif u < 0.5:
        s = 'C'
    elif u < 0.75:
        s = 'G'
    else:
        s = 'T'
    return s

# Uses transition probabilities psame, ptrs, and ptrv to choose state at opposite end
# of the edge given the starting state s0 (possible states are 'A', 'C', 'G', 'T')
# Returns end state s1 as well as nsame, ntrs, and ntrv, each of which is either 0 or 1
# depending on the type of substitution (if any) that occurred
def chooseStateConditional(s0, psame, ptrs, ptrv):
    s1 = None
    nsame = 0
    ntrs  = 0
    ntrv  = 0
    u = random.random()
    if u < psame:
        nsame += 1
        s1 = s0
    elif u < psame + ptrs:
        ntrs = 1
        if s0 == 'A':
            s1 = 'G'
        elif s0 == 'C':
            s1 = 'T'
        elif s0 == 'G':
            s1 = 'A'
        else:
            s1 = 'C'
    else:
        ntrv = 1
        u = random.random()
        if s0 == 'A' or s0 == 'G':
            if u < 0.5:
                s1 = 'C'
            else:
                s1 = 'T'
        else:
            if u < 0.5:
                s1 = 'A'
            else:
                s1 = 'G'
    return s1,nsame,ntrs,ntrv

# Simulates two sequences of length n nucleotides, one at each end of the single edge
# of length v0, given the K80 model with transition/transversion rate ratio k0
def simulateData(v0, k0, n):
    sequence0 = ''
    sequence1 = ''
    psame, ptrs, ptrv = transitionProbs(v0, k0)
    nsame = 0
    ntrs = 0
    ntrv = 0
    for i in range(n):
        s0 = chooseRandomState()
        s1,is_same,is_trs,is_trv = chooseStateConditional(s0, psame, ptrs, ptrv)
        sequence0 += s0
        sequence1 += s1
        nsame += is_same
        ntrs  += is_trs
        ntrv  += is_trv
    return (sequence0, sequence1, nsame, ntrs, ntrv)
    
# Method for tuning deltav based on this paper:
# V Prokaj. 2009. Proposal selection for MCMC simulation. pp. 61-65 in: 
# Applied Stochastic Models and Data Analysis. XIII international conference 
# on applied stochastic models and data analysis. Vilnius, Lithuania.
def tunev(accepted):
    global deltav, vupdates, vaccepts
    vupdates += 1
    gamma_n = 10.0/(100.0 + vupdates)
    if accepted:
        vaccepts += 1
        deltav *= 1.0 + gamma_n*(1.0 - target_acceptance)/(2.0*target_acceptance)
    else:
        deltav *= 1.0 - gamma_n*0.5

    # Prevent run-away increases in boldness for low-information marginal densities
    if deltav > 1000.0:
        deltav = 1000.0

# Method for tuning deltav based on this paper:
# V Prokaj. 2009. Proposal selection for MCMC simulation. pp. 61-65 in: 
# Applied Stochastic Models and Data Analysis. XIII international conference 
# on applied stochastic models and data analysis. Vilnius, Lithuania.
def tunek(accepted):
    global deltak, kupdates, kaccepts
    kupdates += 1
    gamma_n = 10.0/(100.0 + kupdates)
    if accepted:
        kaccepts += 1
        deltak *= 1.0 + gamma_n*(1.0 - target_acceptance)/(2.0*target_acceptance)
    else:
        deltak *= 1.0 - gamma_n*0.5

    # Prevent run-away increases in boldness for low-information marginal densities
    if deltak > 1000.0:
        deltak = 1000.0

# Updates the kappa parameter using a proposal window of width deltak centered over
# the current value of kappa. Globals lnL0 and lnPk0 hold the log likelihood and 
# log prior of the current state, and both are updated if a proposed new value of 
# kappa is accepted. The parameter beta is the power of the power posterior distribution
# being sampled and is used only during steppingstone sampling (for MCMC it equals 1.0)
def updatek(beta, tune):
    global kappa, lnL0, lnPk0, kaccepts, kupdates, deltak

    k0 = kappa
    k = (k0 - deltak/2.0) + random.random()*deltak
    if k < 0.0:
        k = -k
    lnL1 = logLikelihood(nsame, ntrs, ntrv, edgelen, k)
    lnPk1 = logPriorK(k)
    
    logratio = beta*((lnL1 + lnPk1) - (lnL0 + lnPk0))
    if beta < 1.0:
        lnR0 = logRefDistK(k0)
        lnR1 = logRefDistK(k)
        logratio += (1.-beta)*(lnR1 - lnR0)

    logu = math.log(random.random())
    accepted = False
    if logu < logratio:
        # accept
        accepted = True
        kappa = k
        lnL0  = lnL1
        lnPk0 = lnPk1
    else:
        # reject
        kappa = k0

    if tune:        
        tunek(accepted)        
    
# Updates the edgelen parameter using a proposal window of width deltav centered over
# the current value of edgelen. Globals lnL0 and lnPk0 hold the log likelihood and 
# log prior of the current state, and both are updated if a proposed new value of 
# kappa is accepted. The parameter beta is the power of the power posterior distribution
# being sampled and is used only during steppingstone sampling (for MCMC it equals 1.0)
def updatev(beta, tune):
    global edgelen, lnL0, lnPv0, vaccepts, vupdates, deltav

    v0 = edgelen
    v = (v0 - deltav/2.0) + random.random()*deltav
    if v < 0.0:
        v = -v
    lnL1 = logLikelihood(nsame, ntrs, ntrv, v, kappa)
    lnPv1 = logPriorV(v)

    logratio = beta*((lnL1 + lnPv1) - (lnL0 + lnPv0))
    if beta < 1.0:
        lnR0 = logRefDistV(v0)
        lnR1 = logRefDistV(v)
        logratio += (1.-beta)*(lnR1 - lnR0)

    logu = math.log(random.random())
    accepted = False
    if logu < logratio:
        # accept
        accepted = True
        edgelen = v
        lnL0  = lnL1
        lnPv0 = lnPv1
    else:
        # reject
        edgelen = v0

    if tune:
        tunev(accepted)        

# Conducts an MCMC simulation that runs for niters iterations, updating all parameters each iteration.
# A sample is taken of all parameters every thinby iterations.
def MCMC():
    global sample, lnL0, lnPv0, lnPk0, edgelen, kappa
    print('  Using %d iterations (thinning by %d)' % (niters, thinby))
    edgelen = starting_edgelen
    kappa = starting_kappa
    vsample = []
    ksample = []
    lnL0 = logLikelihood(nsame, ntrs, ntrv, edgelen, kappa)
    lnPv0 = logPriorV(edgelen)
    if not do_jc:
        lnPk0 = logPriorK(kappa)
    for i in range(burnin):
        updatev(1.0, True)
        if not do_jc:
            updatek(1.0, True)
    for i in range(niters):
        updatev(1.0, False)
        if not do_jc:
            updatek(1.0, False)
        if (i+1) % thinby == 0:
            lnP0 = lnPv0
            if not do_jc:
                lnP0 += lnPk0
            sample.append((lnL0 + lnP0, edgelen, kappa))
            
def createReferenceDistributionPlotFile(log_marglike):
    vmax = 1.0
    kmax = 50.0
    nincr = 100
    
    if not do_jc:
        # 2-parameter version not yet implemented
        return
    
    outf = open('edgelen-refdist.R', 'w')
    vvect = []
    pvect = []
    rvect = []
    for i in range(1, nincr+1):
        v = vmax*i/nincr
        vvect.append(v)
        lnL = logLikelihood(nsame, ntrs, ntrv, v, assumed_kappa)
        lnP = logPriorV(v)
        lnK = lnL + lnP - log_marglike
        posterior_prob = math.exp(lnK)
        pvect.append(posterior_prob)
        lnR = logRefDistV(v)
        refdist_prob = math.exp(lnR)
        rvect.append(refdist_prob)
    ymax = 1.05*max([max(pvect),max(rvect)])
        
    outf.write('cwd = system(\'cd "$( dirname "$0" )" && pwd\', intern = TRUE)\n')
    outf.write('setwd(cwd)\n')
    outf.write('pdf("edgelength.pdf")\n')
    outf.write('v <- c(%s)\n' % ','.join(['%g' % v for v in vvect]))
    outf.write('p <- c(%s)\n' % ','.join(['%g' % p for p in pvect])) # posterior probability density normalized with specified log_marglike
    outf.write('r <- c(%s)\n' % ','.join(['%g' % r for r in rvect])) # reference probability density
    outf.write('plot(v, p, type="l", lwd=2, lty="solid", col="black", xlim=c(0,%g), ylim=c(0,%g), xlab="v", ylab="prob. density")\n' % (vmax, ymax))
    outf.write('lines(v, r, lwd=2, lty="dotted", col="red")\n')
    outf.write('dev.off()\n')
    outf.close()
                    
# Conducts a generalized steppingstone analysis starting with edgelen equal to starting_edgelen and 
# kappa equal to starting_kappa. nstones MCMC analyses will be performed, each exploring a different
# power posterior distribution, and each running for niters_per_stone iterations after burning in for
# burnin_per_stone iterations. Sampling occurs every thin_per_stone iterations. 
# Assumes an MCMC analysis has already been performed for the posterior
# distribution and the global variable sample holds the samples from that analysis.
def SS():
    global ss_kappa_alpha, ss_kappa_beta, ss_v_alpha, ss_v_beta, vaccepts, vupdates, kaccepts, kupdates, edgelen, kappa

    inv_alpha = 1.0/ssalpha
    beta_k = [math.pow(float(k)/nstones, inv_alpha) for k in range(nstones+1)]
    
    print('  Using %d stone%s, %d iterations/stone (thinning by %d)' % (nstones, nstones == 1 and '' or 's', niters_per_stone, thin_per_stone))
    print('  Power for power posterior from Beta(%g,1) distribution' % ssalpha)
    print('  Tuning parameters used for steppingstone:')
    print('    deltav = %.1f' % deltav)
    if not do_jc:
        print('    deltak = %.1f' % deltak)

    # Compute reference distribution for edge lengths
    edgelens = [v for (lnp,v,k) in sample]
    d = scipy.stats.describe(edgelens)
    ss_v_alpha = d.mean*d.mean/d.variance
    ss_v_beta = d.variance/d.mean
    print('  Reference distribution for edgelen is Gamma(%.5f, %.5f)' % (ss_v_alpha, ss_v_beta))
    
    # Compute reference distribution for kappas
    if not do_jc:
        kappas = [k for (lnp,v,k) in sample]
        d = scipy.stats.describe(kappas)
        ss_kappa_alpha = d.mean*d.mean/d.variance
        ss_kappa_beta = d.variance/d.mean
        print('  Reference distribution for kappa is Gamma(%.5f, %.5f)' % (ss_kappa_alpha, ss_kappa_beta))
    
    edgelen = starting_edgelen
    kappa = starting_kappa
    #ssincr = 1.0/nstones
    log_marglike = 0.0
    print('  Starting steppingstone analysis')
    for i in range(nstones):
        #ssbeta = float(i)/nstones
        ssbeta = beta_k[i]
        ssincr = beta_k[i+1] - beta_k[i]
        
        # Burn-in and adapt updaters
        vaccepts = 0 
        vupdates = 0
        kaccepts = 0
        kupdates = 0
        for j in range(burnin_per_stone):
            updatev(ssbeta, True)
            if not do_jc:
                updatek(ssbeta, True)
                            
        # Sample from power posterior
        tmp = []
        for j in range(niters_per_stone):
            updatev(ssbeta, False)
            if not do_jc:
                updatek(ssbeta, False)

            #input('.. v = %.5f ..' % (edgelen,))    

            if (j+1) % thin_per_stone == 0:
                lnL = logLikelihood(nsame, ntrs, ntrv, edgelen, kappa)
                lnPv = logPriorV(edgelen)
                lnRv = logRefDistV(edgelen)
                logr = lnL + lnPv - lnRv
                if not do_jc:
                    lnPk = logPriorK(kappa)
                    lnRk = logRefDistK(kappa)
                    logr += lnPk - lnRk
                tmp.append(ssincr*logr)

        log_sum_terms = logSum(tmp)
        logn = math.log(len(tmp))
        logrk = log_sum_terms - logn            
        log_marglike += logrk
        print('    Step %d of %d: n = %d, beta = %.5f, diff = %.5f, logrk = %.5f, logmarglike = %.5f' % (i+1,nstones,len(tmp),ssbeta,ssincr,logrk,log_marglike))

    createReferenceDistributionPlotFile(log_marglike)
        
    return log_marglike

# Assumes the global variable sample contains a sample from the posterior distribution and
# is a list of 3-tuples (log posterior kernel, sampled edgelen, sampled kappa). This function
# log-transforms the sampled values and then standardizes them by subtracting the mean vector 
# (meanvect) and premultiplying by the inverse of the standard deviation matrix (sqrtS). If
# training is True, the meanvect and sqrtS (as well as the log-determininant of sqrtS) are 
# calculated before being used. If training is False, the meanvect and sqrtS are simply 
# used (it is assumed they have been previously calculated). 
#
# The log transformation involves creating random variable Y by taking the log of a 
# random variable X:
#
#    Y = log(X)   X = exp(Y)   |dX/dY| = exp(Y)
#
# The Jacobian for this transformation is thus exp(y), which equals just y on log scale.
#
# The standardization transformation i:
#
#    Y = S^{-0.5} (X - meanvect)  X = S^{0.5} Y + meanvect   |dX/dY| = |S^{0.5}|
#
# The Jacobian for this transformation is thus the determinant of S^{0.5}, which equals 
# log det S^{0.5} on log scale. This value is stored as the global logdetsqrtS.
def transformSample(sample, training):
    global meanvect, S, sqrtS, invsqrtS, logdetsqrtS
    
    n = len(sample)
    
    if do_jc:
        # Using JC69 model

        if training:
            # calculate mean vector
            logvmean = 0.0
            for i in range(n):
                logv = math.log(sample[i][1])
                logvmean += logv
            logvmean /= n
            meanvect = logvmean
        
            # log-transform the edge length parameter and compute variance
            S = 0.0
            for logkernel,v,k in sample:
                a = math.log(v) - meanvect
                S += a*a
            S /= n-1;
            sqrtS = math.sqrt(S)
            invsqrtS = 1.0/sqrtS
            detsqrtS = math.fabs(sqrtS)
            logdetsqrtS = math.log(detsqrtS)

        if training and do_show_transformed_mean_vector:
            print('Mean vector:')
            print(meanvect)
    
        if training and do_show_transformed_varcov_matrix:
            print('Variance-covariance matrix:')
            print(S)
    
        transformed = []
        edge_lengths = []
        kappas = []
        for i in range(n):
            kernel = sample[i][0]
            v = sample[i][1]
            logv = math.log(v)
            vect = logv
            stdvect = (vect - meanvect)/sqrtS
            edge_lengths.append(stdvect)
            kappas.append(1.0)
            newkernel = kernel + logv + logdetsqrtS
            transformed.append((newkernel, stdvect, 0.0))
    else:
        # Using K80 model
        
        if training:
            # calculate mean vector
            logvmean = 0.0
            logkmean = 0.0
            n = len(sample)
            for i in range(n):
                logv = math.log(sample[i][1])
                logk = math.log(sample[i][2])
                logvmean += logv
                logkmean += logk
            logvmean /= n
            logkmean /= n
            meanvect = numpy.array([[logvmean], [logkmean]]) # note shape of numpy array is (2,1) (i.e. vector shape)
    
            # log-transform both parameters and compute sample variance-covariance matrix
            S = numpy.zeros((2,2))
            for logkernel,v,k in sample:
                a = numpy.array([[math.log(v)],[math.log(k)]]) - meanvect
                aT = a.transpose()
                S += a.dot(aT)
            S /= n-1;
            sqrtS = scipy.linalg.sqrtm(S)
            invsqrtS = scipy.linalg.inv(sqrtS)
            detsqrtS = scipy.linalg.det(sqrtS)
            logdetsqrtS = math.log(detsqrtS)
    
        if training and do_show_transformed_mean_vector:
            print('Mean vector:')
            print(meanvect)
    
        if training and do_show_transformed_varcov_matrix:
            print('Variance-covariance matrix:')
            print(S)
    
        transformed = []
        edge_lengths = []
        kappas = []
        for i in range(n):
            kernel = sample[i][0]
            v = sample[i][1]
            k = sample[i][2]
            logv = math.log(v)
            logk = math.log(k)
            vect = numpy.array([[logv],[logk]])
            stdvect = invsqrtS.dot(vect - meanvect)
            edge_lengths.append(stdvect[0][0])
            kappas.append(stdvect[1][0])
            newkernel = kernel + logv + logk + logdetsqrtS
            transformed.append((newkernel, stdvect[0][0], stdvect[1][0]))

    dv = scipy.stats.describe(edge_lengths, 0, 1)
    dk = scipy.stats.describe(kappas, 0, 1)
    return (dv.mean, dk.mean, transformed)

# Reverses the standardization and log-transformation for supplied transformed edgelen
# (stdlogv) and transformed kappa (stdlogk).
#
# if is_standardized = False, then it is assumed that the supplied points have been
# log transformed but not standardized.
# 
# If X = (v,k) is a vector of untransformed values, logX = (log(v), log(k)) the vector
# of log-transformed values, the standardization transformation creates vector Y
# such that
#   Y = S^{-0.5} (logX - meanLogX) = sqrtS^{-1} (logX - meanLogX)
# where meanLogX is the mean vector of the logX vectors and S is the variance-covariance
# matrix of the log-transformed parameter vectors. To reverse both transformations, first
# destandardize:
#   logX = sqrtS Y + meanLogX
# Then exponentiate both elements of logX to create the X vector.
def detransform(stdlogv, stdlogk, is_standardized = True):
    if do_jc:
        if is_standardized:
            logx = sqrtS.dot(stdlogv) + meanvect
        else:
            logx = stdlogv
        v = math.exp(logx)
        return (v,assumed_kappa)
    else:
        # destandardize point p
        if is_standardized:
            y = numpy.array([[stdlogv],[stdlogk]])
            logx = sqrtS.dot(y) + meanvect
        else:
            logx = [stdlogv, stdlogk]
    
        # de-log-transform point
        v = math.exp(logx[0])
        k = math.exp(logx[1])
        return (v,k)

# Computes log posterior kernel for untransformed point (v, k).
def calcLogUntransformedKernel(v, k):
    logkernel = logLikelihood(nsame, ntrs, ntrv, v, k)
    logkernel += logPriorV(v)
    if not do_jc:
        logkernel += logPriorK(k)
    return logkernel

# Computes log posterior kernel for transformed point (stdlogv, stdlogk).
# The values logv and logk represent the Jacobian for the log transformation on log scale.
# the value logdetsqrtS is the Jacobian for the standardization transformation on log scale.
def calcLogTransformedStandardizedKernel(stdlogv, stdlogk):
    v,k = detransform(stdlogv, stdlogk)
    logkernel = logLikelihood(nsame, ntrs, ntrv, v, k)
    logkernel += logPriorV(v)
    if not do_jc:
        logkernel += logPriorK(k)
    logv = math.log(v)
    transformed_logkernel = logkernel + logv + logdetsqrtS
    if not do_jc:
        logk = math.log(k)
        transformed_logkernel += logk
    return transformed_logkernel

# Computes log posterior kernel for the transformed point (logv, logk).
# The values logv and logk represent the Jacobian for the log transformation on log scale.
def calcLogTransformedKernel(logv, logk):
    v,k = detransform(logv, logk, False)  # False argument says logv,logk are not standardized
    logkernel = logLikelihood(nsame, ntrs, ntrv, v, k)
    logkernel += logPriorV(v)
    if not do_jc:
        logkernel += logPriorK(k)
    #logv = math.log(v)
    transformed_logkernel = logkernel + logv
    if not do_jc:
        #logk = math.log(k)
        transformed_logkernel += logk
    return transformed_logkernel

# Computes the log of the sum of a vector of values that are on log scale. The largest 
# values in the vector is factored out to avoid underflow issues. 
# 
# For example, let a + b + c be a sum of interest and c is the largest element. 
# Factor out c, yielding
#   c (a/c + b/c + 1) = c (e^{log a - log c} + e^{log b - log c} + 1)
# This function returns the sum on log scale, i.e.:
#   log(c) + log(e^{log a - log c} + e^{log b - log c} + 1)
# There is no underflow because the sum is guaranteed to be greater than 1.
def logSum(logx_vect):
    max_logx = max(logx_vect)
    sum_terms = 0.0
    for logx in logx_vect:
        sum_terms += math.exp(logx - max_logx)
    logsum = max_logx + math.log(sum_terms)
    return logsum
    
# Computes mean, minimum, and maximum of one parameter in sample_vect, with the parameter
# determined by the argument supplied to the "which" parameter. If cutoff is supplied, the
# first cutoff samples are ignored.
def calcModeMeanMinMaxVar(sample_vect, which, cutoff = 0):
    vect = []
    mode = None
    qmode = None
    for s in sample_vect[cutoff:]:
        vect.append(s[which])
        if qmode is None or s[0] > qmode:
            qmode = s[0]
            mode = s[which]
    d = scipy.stats.describe(vect)
    return (mode, d.mean, d.minmax[0], d.minmax[1], d.variance)
    
# Computes the effective sample size (ESS) of one parameter in the sample vector 
# (determined by argument supplied to which parameter) 
# Assumes AR1 autoregression
# ESS = N (1/(1 + 2 sum_{k=1}^{\infty} \rho_k))
#     = N (1/(1 + 2 (1/(1-\rho) - 1)))
#     = N (1-\rho)/(1+\rho)
def effectiveSampleSizeSimple(which):
    earlier = [s[which] for s in sample[:-1]]
    later   = [s[which] for s in sample[1:]]
    r1, pvalue = scipy.stats.pearsonr(earlier, later)
    N = float(len(sample))
    ESS = N*(1.0 - r1)/(1.0 + r1)
    tau = (1.0 + r1)/(1.0 - r1)
    return ESS,tau

# Comparison of ESS methods (simple method indicated by 0):
#    niters  thinby          ESSv maxlag       ESSv0         ESSk  maxlag        ESSk0
# 100000000   10000   10375.54282      5  9761.20970   9955.55812       6   9772.50023
#  10000000    1000   10180.65908      5 10112.46689  10514.80579       5  10308.73645
#   1000000     100   10212.94018      5 10112.27825  10684.40472       5  10225.68492
#    100000      10    5167.36978     10  5661.20041   5244.16875      10   5394.01535
#     10000       1     529.70452     95   966.35427    569.69365      88    799.41454

# Computes the effective sample size (ESS) of one parameter in the sample vector 
# (determined by argument supplied to which parameter) 
# see https://dfm.io/posts/autocorr/
def effectiveSampleSize(which):
    f = [s[which] for s in sample]
    N = len(f)
    mu_f = sum(f)/N
    cf0 = sum([math.pow(fn - mu_f,2.) for fn in f])/N
    rhoftau = []
    M = 1
    C = 5.0
    taufM = None
    for tau in range(1,N+1):
        cftau = 0.0
        for n in range(N-tau):
            cftau += (f[n] - mu_f)*(f[n+tau] - mu_f)
        cftau /= (N-tau)
        rhoftau.append(cftau/cf0)
        taufM = 1.0 + 2.0*sum(rhoftau)
        if M >= C*taufM:
            break
        M += 1
    ESS = float(N)/taufM
    return ESS, M, taufM

# Used by surface plotting functions, this function creates a square 2-dimensional list of 
# density heights.
# nincr is the number of rows and also the number of columns
# linear_scale: if True, use density; if False, use log-density
# vmin, vmax, kmin, kmax are the axis limits in the space appropriate for plot_type
# plot_type determines how density height is defined, and should be one of these strings:
#   'transformed-standardized-posterior' defines height to be posterior of transformed and standardized parameters (used in Figures 1b and 2 in paper)
#   'transformed-unstandardized-posterior' defines height to be posterior of the transformed (but not standardized) parameters
#   'untransformed-posterior' defines height to be the posterior of the untransformed parameters (used in Figure 1a in the paper)
#   'untransformed-refdist' defines height to be the reference distribution density
#   'mvnorm' defines height to be the multivariate normal density
#   'mvstdnorm' defines height to be the multivariate standard normal density (used in Figure 2 in the paper)
# fn is a file name that, if not None, cuases the 2-dimensional list to be stored in csv format
def createZArray(nincr, linear_scale, log_norm_const, vmin, vmax, kmin, kmax, plot_type, fn = None):
    dv = (vmax - vmin)/nincr
    dk = (kmax - kmin)/nincr

    xx = []
    yy = []
    zz = []
    zzmin = None
    zzmax = None
    
    if fn is not None:
        zfile = open(fn, 'w')
    xx = [(vmin + dv*j) for j in range(1,nincr-1)]
    yy = [(kmin + dk*j) for j in range(1,nincr-1)]
    for i in range(1,nincr-1): 
        kappa_value = kmin + dk*i
        row = []
        for j in range(1,nincr-1):
            edgelen_value = vmin + dv*j
            log_scale = not linear_scale
            if plot_type == 'transformed-standardized-posterior' and log_scale:
                log_kernel = calcLogTransformedStandardizedKernel(edgelen_value, kappa_value)
                h = log_kernel - log_norm_const
            elif plot_type == 'transformed-standardized-posterior' and linear_scale:
                log_kernel = calcLogTransformedStandardizedKernel(edgelen_value, kappa_value)
                h = math.exp(log_kernel - log_norm_const)
            elif plot_type == 'transformed-unstandardized-posterior' and log_scale:
                log_kernel = calcLogTransformedKernel(edgelen_value, kappa_value)
                h = log_kernel - log_norm_const
            elif plot_type == 'transformed-unstandardized-posterior' and linear_scale:
                log_kernel = calcLogTransformedKernel(edgelen_value, kappa_value)
                h = math.exp(log_kernel - log_norm_const)
            elif plot_type == 'untransformed-posterior' and linear_scale:
                log_kernel = calcLogUntransformedKernel(edgelen_value, kappa_value)
                h = math.exp(log_kernel - log_norm_const)
            elif plot_type == 'untransformed-refdist' and linear_scale:
                log_kernel  = logRefDistV(edgelen_value)
                log_kernel += logRefDistK(kappa_value)
                h = math.exp(log_kernel - log_norm_const)
            elif plot_type == 'mvnorm' and log_scale:
                meanv = meanvect[0]
                meank = meanvect[1]
                sigmav = sqrtS[0][0]
                sigmak = sqrtS[1][1]
                rho = S[0][1]/(sigmav*sigmak)
                rho_term = 1. - pow(rho, 2.)
                kappa_term = (kappa_value-meank)/sigmak
                edgelen_term = (edgelen_value-meanv)/sigmav
                rsquared  = math.pow(kappa_term,2.)
                rsquared += math.pow(edgelen_term,2.)
                rsquared -= 2.*rho*kappa_term*edgelen_term
                h = -0.5*rsquared/rho_term - math.log(2.*math.pi*sigmav*sigmak) - 0.5*math.log(rho_term) - log_norm_const
            elif plot_type == 'mvnorm' and linear_scale:
                meanv = meanvect[0]
                meank = meanvect[1]
                sigmav = sqrtS[0][0]
                sigmak = sqrtS[1][1]
                rho = S[0][1]/(sigmav*sigmak)
                rho_term = 1. - pow(rho, 2.)
                kappa_term = (kappa_value-meank)/sigmak
                edgelen_term = (edgelen_value-meanv)/sigmav
                rsquared  = math.pow(kappa_term,2.)
                rsquared += math.pow(edgelen_term,2.)
                rsquared -= 2.*rho*kappa_term*edgelen_term
                h = math.exp( -0.5*rsquared/rho_term - math.log(2.*math.pi*sigmav*sigmak) - 0.5*math.log(rho_term) - log_norm_const )
            elif plot_type == 'mvstdnorm' and log_scale:
                rsquared = math.pow(kappa_value,2.) + math.pow(edgelen_value,2.)
                h = -0.5*rsquared - math.log(2.*math.pi) - log_norm_const
            elif plot_type == 'mvstdnorm' and linear_scale:
                rsquared = math.pow(kappa_value,2.) + math.pow(edgelen_value,2.)
                h = math.exp( -0.5*rsquared - math.log(2.*math.pi) - log_norm_const )
            else:
                sys.exit('combination of plot_type "%s" and %s scale is not recognized in createZArray' % (plot_type, linear_scale and "linear" or "log"))

            if zzmin is None or h < zzmin:
                zzmin = h

            if zzmax is None or h > zzmax:
                zzmax = h

            row.append(h)
        zz.append(row)
            
        if fn is not None:
            zfile.write('%s\n' % (','.join(['%.5f' % z for z in row])))
            
    if fn is not None:
        zfile.close()
    
    return (xx, yy, zz, zzmin, zzmax)

# Uses plotly.graph_objs to create a plot showing the posterior or reference distribution
# density surface. See https://plotly.com/python/3d-surface-plots/
# fn is the file name to use for storing the plot if plots_saved_to_file is True
# linear_scale: if True, show density; if False, show log-density
# vmin, vmax, kmin, kmax are the axis limits (these should be chosen to correspond to
#   the selected plot_type (see below)
# zsep is the amount to subtract from the calculated zmin in order to separate 
#   contours from the bottom of the plotted surface (zsep = 0.0 means no separation)
# xeye, yeye, zeye specify the position of the eye viewing the plot (1.8, 1.8, 1.0 are good)
# plot_type determines how density height is defined, and should be one of these strings:
#   'transformed-standardized-posterior' defines height to be posterior of transformed and standardized parameters
#   'transformed-unstandardized-posterior' defines height to be posterior of the transformed (but not standardized) parameters
#   'untransformed-posterior' defines height to be the posterior of the untransformed parameters
#   'untransformed-refdist' defines height to be the reference distribution density
#   'mvnorm' defines height to be the multivariate normal density
#   'mvstdnorm' defines height to be the multivariate standard normal density
def plotSurfaces(fn, 
        linear_scale, 
        indep_color_scales, 
        vmin, vmax, kmin, kmax, zsep,
        vticks, kticks,
        xeye, yeye, zeye,
        plot_types, 
        color_scales, 
        log_normalization_constants, 
        opacities):
        
    assert len(plot_types) > 0
    if do_jc:
        return
        
    nincr = plotting_increments
    surfaces = []
    zaxismin = None
    zaxismax = None
    
    for plot_type,color_scale,log_norm_const,op in zip(plot_types, color_scales, log_normalization_constants, opacities):   
        if linear_scale:
            norm_const = math.exp(log_norm_const)
         
        # Create surface
        xx, yy, zz, zzmin, zzmax = createZArray(nincr, linear_scale, log_norm_const, vmin, vmax, kmin, kmax, plot_type)
        #print('    zzmax for plot_type="%s" is' % plot_type, zzmax)

        if zaxismin is None:
            zaxismin = zzmin
        else:
            zaxismin = min([zaxismin, zzmin])
        
        if zaxismax is None:
            zaxismax = zzmax
        else:
            zaxismax = max([zaxismax, zzmax])
        
        if indep_color_scales:
            surf = go.Surface(x=xx, y=yy, z=zz,
                colorscale = color_scale,
                cmin       = zzmin, 
                cmax       = zzmax,
                opacity    = op
            )
        else:
            surf = go.Surface(x=xx, y=yy, z=zz,
                colorscale = color_scale,
                opacity    = op
            )
        surfaces.append(surf)
        
    if not indep_color_scales:
        # Make all plots use same color scale
        for surf in surfaces:
            surf.cmin = zaxismin
            surf.cmax = zaxismax
  
    # Add surface plots to the figure
    fig = go.Figure(data=surfaces)
    
    fig.update_traces(
        showscale  = False,
        contours_z = dict(
            show           = True, 
            usecolormap    = True, 
            highlightcolor = "limegreen", 
            project_z      = True)
    )

    # Construct tick labels for the x-axis (edge length parameter v)
    if vticks == 'auto':
        vtickvals = [(vmin + (vmax-vmin)*i/5.) for i in range(6)]
    else:
        vtickvals = vticks
    vticktext = ['%.1f' % vtick for vtick in vtickvals]

    # Construct tick labels for the y-axis (transition-transversion rate ratio parameter k)
    if kticks == 'auto':
        ktickvals = [(kmin + (kmax-kmin)*i/5.) for i in range(6)]
    else:
        ktickvals = kticks
    kticktext = ['%.1f' % ktick for ktick in ktickvals]
    
    fig.update_layout(
        #title                        = main_title,
        #height=1000,  # must make height = width in addition to making range the same
        #width=1000,  # for both x and y axes in order to get aspect ratio 1
        #autosize=False

        scene_camera_up               = dict(x=0, y=0, z=1),
        scene_camera_center           = dict(x=0, y=0, z=0),
        scene_camera_eye              = dict(x=xeye, y=yeye, z=zeye),
        scene_camera_projection_type  = 'perspective',
        #scene_camera_projection_type  = 'orthographic',

        scene_xaxis_tickvals          = vtickvals,
        scene_xaxis_ticktext          = vticktext,
        scene_xaxis_tickmode          = 'array',
        scene_xaxis_title_text        = 'edge length',
        #scene_xaxis_title_text        = '&nu', # supposed to work (https://plotly.com/chart-studio-help/adding-HTML-and-links-to-charts/), but doesn't
        #scene_xaxis_title_font_family = 'Helvetica',
        scene_xaxis_title_font_size   = 12,
        #scene_xaxis_title_font_color  = 'navy',

        yaxis_scaleanchor             = "x",
        yaxis_scaleratio              = 1,
        scene_yaxis_tickvals          = ktickvals,
        scene_yaxis_ticktext          = kticktext,
        scene_yaxis_tickmode          = 'array',
        #scene_yaxis_title_text        = '&kappa', # supposed to work (https://plotly.com/chart-studio-help/adding-HTML-and-links-to-charts/), but doesn't
        scene_yaxis_title_text        = 'trs/trv rate ratio',
        #scene_yaxis_title_font_family = 'Helvetica',
        scene_yaxis_title_font_size   = 12,
        #scene_yaxis_title_font_color  = 'navy',

        scene_zaxis_range             = [zaxismin - zsep, zaxismax],
        scene_zaxis_title             = 'posterior kernel',
        scene_zaxis_visible           = False,
        #scene_zaxis_title_font_family = 'Times',
        #scene_zaxis_title_font_size   = 12,
        #scene_zaxis_title_font_color  = 'navy',

        autosize                      = False,
        paper_bgcolor                 = "white",
        #margin_l                      = 100,
        #margin_r                      = 100,
        #margin_t                      = 100,
        #margin_b                       = 100,
        #margin_pad                     = 100
    )
    
    #fig.write_json('%s.json' % fn, pretty=True)
    if plots_saved_to_file:
        fig.write_image(fn, height=plotting_height,scale=plotting_scale) # this creates a static file  
    else:
        fig.show(renderer='browser') # this opens your browser to show you the plot now
        
    return fig, zaxismax

# Creates 2-D plot showing points sampled from a bivariate standard normal distribution
# Those points whose norms are less than rmax are shown with arrows.
# A circle shows the limit defined by rmax.
# No-op if do_jc is True.
# See https://plotly.com/python/reference/scattergl/ for go.Scatter options
def plotNorms(fn, xmin, xmax, ymin, ymax, rmax, points):
    if do_jc:
        return
        
    xvect = [p[0] for p in points]
    yvect = [p[1] for p in points]
    
    cxvect = []
    cyvect = []
    nincr = 100
    for i in range(nincr+1):
        theta = 2.0*math.pi*i/nincr
        cxvect.append(rmax*math.cos(theta))
        cyvect.append(rmax*math.sin(theta))
        
    scatter = go.Scatter(x=xvect, y=yvect, mode='markers', marker=dict(size=5, opacity=1.0, color="black"), showlegend=False)
    circle  = go.Scatter(x=cxvect, y=cyvect, mode='lines', line=dict(color="navy",dash="dot"), showlegend=False)
    
    fig = go.Figure(data=[scatter,circle])
    for p in points:
        # annotate with arrow if norm of p < rmax
        x = p[0]
        y = p[1]
        norm = math.sqrt(x*x + y*y)
        if norm <= rmax:
            print('adding arrow from 0,0 to %g,%g' % (x,y))
            fig.add_shape(
                    type="line",
                    x0=0.0,
                    y0=0.0,
                    x1=x,
                    y1=y,
                    xref="x",
                    yref="y",
                    line=dict(color="DarkOrange",width=2),
                    )

    fig.update_xaxes(range=[xmin, xmax])
    fig.update_yaxes(range=[ymin, ymax])
    fig.update_layout(
        height=1000,  # must make height = width in addition to making range the same
        width=1000,  # for both x and y axes in order to get aspect ratio 1
        autosize=False
    )
    #fig.show()
    fig.write_image(fn)
    
def doSimulationMCMC():
    global sample, rnseed, nsamples, T
    global nsame, ntrs, ntrv
    global lnL0, lnPv0, lnPk0
    global edgelen, vaccepts, vupdates, deltav
    global kappa, kaccepts, kupdates, deltak
    
    if do_jc:
        print('Using Jukes-Cantor (1969) model (estimating edge length only)')
    else:
        print('Using Kimura (1980) model (estimating edge length and kappa)')

    # Set the random number seed          
    random.seed(rnseed)
    print('Using random number seed: %d' % rnseed)

    # Attempt to read previously-generated results
    mcmc_samples_exist = False
    simulation_settings_modified = False
    sample = []
    nsamples = 0
    in_rnseed = 0
    in_true_edgelen = 0.0
    in_true_kappa = 0.0
    in_seqlen = 0
    in_sequence0 = ''
    in_sequence1 = ''
    in_nsame = 0
    in_ntrs = 0
    in_ntrv = 0
    in_nsamples = 0
    lnL0 = 0.0
    lnPv0 = 0.0
    lnPk0 = 0.0
    last_sampled_edgelen = 0.0
    last_sampled_kappa = 0.0
    if os.path.exists('mcmc-info.txt'):
        mcmc_samples_exist = True
        mcmc_lines = open('mcmc-info.txt').readlines()
        for line in mcmc_lines:
            m = re.match('rnseed:\s+(\d+)', line)
            if m is not None:
                in_rnseed = int(m.group(1))            
                #print('Found rnseed (%d) in mcmc-samples.txt' % in_rnseed)

            m = re.match('true edge length:\s+([.0-9]+)', line)
            if m is not None:
                in_true_edgelen = float(m.group(1))            
                #print('Found true_edgelen (%g) in mcmc-samples.txt' % in_true_edgelen)

            m = re.match('true kappa:\s+([.0-9]+)', line)
            if m is not None:
                in_true_kappa = float(m.group(1))            
                #print('Found true_kappa (%g) in mcmc-samples.txt' % in_true_kappa)

            m = re.match('sequence length:\s+([0-9]+)', line)
            if m is not None:
                in_seqlen = float(m.group(1))            
                #print('Found seqlen (%d) in mcmc-samples.txt' % in_seqlen)

            m = re.match('sequence0:\s+([ACGT]+)', line)
            if m is not None:
                #print('Found sequence0 in mcmc-samples.txt')
                in_sequence0 = m.group(1)

            m = re.match('sequence1:\s+([ACGT]+)', line)
            if m is not None:
                #print('Found sequence1 in mcmc-samples.txt')
                in_sequence1 = m.group(1)

            m = re.match('nsame:\s+(\d+)', line)
            if m is not None:
                in_nsame = int(m.group(1))
                #print('Found nsame (%d) in mcmc-samples.txt' % in_nsame)

            m = re.match('ntrs:\s+(\d+)', line)
            if m is not None:
                in_ntrs = int(m.group(1))
                #print('Found ntrs (%d) in mcmc-samples.txt' % in_ntrs)

            m = re.match('ntrv:\s+(\d+)', line)
            if m is not None:
                in_ntrv = int(m.group(1))
                #print('Found ntrv (%d) in mcmc-samples.txt' % in_ntrv)

            m = re.match('nsamples:\s+(\d+)', line)
            if m is not None:
                in_nsamples = int(m.group(1))
                #print('Found nsamples (%d) in mcmc-samples.txt' % in_nsamples)

            m = re.match('MCMC samples:', line)
            if m is not None:
                #print('Found "MCMC samples:" in mcmc-samples.txt')
                pass
                            
    if os.path.exists('mcmc-samples.txt'):
        mcmc_lines = open('mcmc-samples.txt').readlines()
        for line in mcmc_lines:
            m = re.match('\s+log-kernel\s+edgelen\s+kappa', line)
            if m is not None:
                #print('Found column headers in mcmc-samples.txt')
                pass
                            
            m = re.match('\s*([-.0-9]+)\s+([.0-9]+)\s+([.0-9]+)', line)
            if m is not None:
                in_f = float(m.group(1))
                last_sampled_edgelen = float(m.group(2))
                last_sampled_kappa = float(m.group(3))
                sample.append((in_f, last_sampled_edgelen, last_sampled_kappa))
    else:
        mcmc_samples_exist = False

    if mcmc_samples_exist:
        print('wwwww')
        lnL0 = logLikelihood(nsame, ntrs, ntrv, last_sampled_edgelen, last_sampled_kappa)
        lnPv0 = logPriorV(last_sampled_edgelen)
        lnPk0 = logPriorK(last_sampled_kappa)
        
        nsamples = niters//thinby
        assert nsamples == in_nsamples, 'expecting nsamples to be equal to %d in mcmc-samples.txt, but instead found %d' % (nsamples,in_nsamples)
        assert nsamples == len(sample), 'expecting nsamples (%d) to be equal to the length of the samples vector (%d)' % (nsamples, len(samples))
        if in_rnseed != rnseed or in_true_edgelen != true_edgelen or in_true_kappa != true_kappa or in_seqlen != seqlen: 
            simulation_settings_modified = True
        else:
            nsame = in_nsame
            ntrv  = in_ntrv
            ntrs  = in_ntrs
            
    # If necessary, simulate data and perform MCMC
    if simulation_settings_modified or not mcmc_samples_exist:
        sample = []
        
        # Simulate the data
        sequence0, sequence1, nsame, ntrs, ntrv = simulateData(true_edgelen, true_kappa, seqlen)
        print('\nSimulated data:')
        print('  %d sites are the same' % nsame)
        print('  %d sites show a transition' % ntrs)
        print('  %d sites show a transversion' % ntrv)
        if do_show_sequences:
            print('sequence0: %s' % sequence0)
            print('sequence1: %s' % sequence1)

        vaccepts = 0
        vupdates = 0
        kaccepts = 0
        kupdates = 0
        
        # Perform MCMC analysis starting with true parameter values to avoid need for burn-in
        MCMC()
        
        nsamples = len(sample)

        if save_samples_to_file:
            # Save samples to file mcmc-samples.txt
            outf = open('mcmc-info.txt', 'w')
            outf.write('rnseed:           %d\n' % rnseed)
            outf.write('true edge length: %.5f\n' % true_edgelen)
            outf.write('true kappa:       %.5f\n' % true_kappa)
            outf.write('sequence length:  %d\n' % seqlen)
            outf.write('sequence0:        %s\n' % sequence0)
            outf.write('sequence1:        %s\n' % sequence1)
            outf.write('nsame:            %d\n' % nsame)
            outf.write('ntrs:             %d\n' % ntrs)
            outf.write('ntrv:             %d\n' % ntrv)
            outf.write('nsamples:         %d\n' % nsamples)
            outf.close()
            
            outf = open('mcmc-samples.txt', 'w')
            outf.write('%s\t%s\t%s\n' % ('log-kernel','edgelen','kappa'))
            if do_jc:
                for f,v,dummyk in sample:
                    outf.write('%.9f\t%.5f\t%.5f\n' % (f,v,dummyk))
            else:
                for f,v,k in sample:
                    outf.write('%.9f\t%.5f\t%.5f\n' % (f,v,k))
            outf.close()

        if burnin > 0:
            print('  Tuning parameters:')
            print('    deltav = %.1f (accept %% = %.1f)' % (deltav, 100.*vaccepts/vupdates))
            if not do_jc:
                print('    deltak = %.1f (accept %% = %.1f)' % (deltak, 100.*kaccepts/kupdates))

    # Show MCMC results
    T = float(len(sample))
    #v_accept_pct = 100.0*vaccepts/vupdates
    modev, meanv, minv, maxv, varv = calcModeMeanMinMaxVar(sample, 1)
    ESSv,maxlagv,tauv = effectiveSampleSize(1)
    ESSv0,tauv0 = effectiveSampleSizeSimple(1)
    if not do_jc:
        #k_accept_pct = 100.0*kaccepts/kupdates
        modek, meank, mink, maxk, vark = calcModeMeanMinMaxVar(sample, 2)
        ESSk,maxlagk,tauk = effectiveSampleSize(2)
        ESSk0,tauk0 = effectiveSampleSizeSimple(2)
        
    print('\nMCMC analysis:')
    if mcmc_samples_exist and not simulation_settings_modified:
        print('  Using previous MCMC samples found in file mcmc-samples.txt')
        
    print('  Sample size = %.0f' % T)
    print('  Edge length:')
    print('    mean = %.5f' % meanv)
    print('    mode = %.5f' % modev)
    print('    var  = %.5f' % varv)
    print('    ESS  = %.5f (max. lag. = %d, autocorr. time = %.5f)' % (ESSv,maxlagv,tauv))
    print('    ESS  = %.5f (simple formula, autocorr. time = %.5f)' % (ESSv0,tauv0))
    if not do_jc:
        print('  Kappa:')
        print('    mean = %.5f' % meank)
        print('    mode = %.5f' % modek)
        print('    var  = %.5f' % vark)
        print('    ESS  = %.5f (max. lag. = %d, autocorr. time = %.5f)' % (ESSk,maxlagk,tauk))
        print('    ESS  = %.5f (simple formula, autocorr. time = %.5f)' % (ESSk0,tauk0))
        
def doSS():
    if not do_steppingstone:
        return
        
    global ss_kappa_alpha, ss_kappa_beta, kaccepts, kupdates, kappa
    global ss_v_alpha, ss_v_beta, vaccepts, vupdates, edgelen
    global log_marglike_ss
        
    # Attempt to read previously-generated results
    ss_results_exist = False
    ss_settings_modified = False
    in_ssalpha = 0.0
    in_nstones = 0
    in_niters_per_stone = 0
    in_thin_per_stone = 0
    in_burnin_per_stone = 0
    in_log_marglike_ss = 0.0
    if os.path.exists('ss-results.txt'):
        ss_results_exist = True
        ss_lines = open('ss-results.txt').readlines()
        for line in ss_lines:
            m = re.match('ssalpha:\s+([.0-9]+)', line)
            if m is not None:
                in_ssalpha = float(m.group(1))            
                #print('Found ssalpha (%g) in ss-results.txt' % in_ssalpha)

            m = re.match('nstones:\s+(\d+)', line)
            if m is not None:
                in_nstones = int(m.group(1))            
                #print('Found nstones (%d) in ss-results.txt' % in_nstones)

            m = re.match('niters_per_stone:\s+(\d+)', line)
            if m is not None:
                in_niters_per_stone = int(m.group(1))            
                #print('Found niters_per_stone (%d) in ss-results.txt' % in_niters_per_stone)

            m = re.match('thin_per_stone:\s+(\d+)', line)
            if m is not None:
                in_thin_per_stone = int(m.group(1))            
                #print('Found thin_per_stone (%d) in ss-results.txt' % in_thin_per_stone)

            m = re.match('burnin_per_stone:\s+(\d+)', line)
            if m is not None:
                in_burnin_per_stone = int(m.group(1))            
                #print('Found burnin_per_stone (%d) in ss-results.txt' % in_burnin_per_stone)

            m = re.match('log_marglike_ss:\s+([-.0-9]+)', line)
            if m is not None:
                in_log_marglike_ss = float(m.group(1))            
                #print('Found log_marglike_ss (%g) in ss-results.txt' % in_log_marglike_ss)

        if in_ssalpha != ssalpha or in_nstones != nstones or in_niters_per_stone != niters_per_stone or in_thin_per_stone != thin_per_stone or in_burnin_per_stone != burnin_per_stone: 
            ss_settings_modified = True

    # If requested, perform steppingstone to obtain estimate of the marginal likelihood
    print('\nSteppingstone analysis:')
    if ss_settings_modified or not ss_results_exist:
        log_marglike_ss = SS()
        print('  Log marginal likelihood (steppingstone) = %.5f' % log_marglike_ss)

        # Save steppingstone settings to ss-results.txt        
        outf = open('ss-results.txt', 'w')
        outf.write('ssalpha:                 %.5f\n' % ssalpha)
        outf.write('nstones:                 %d\n' % nstones)
        outf.write('niters_per_stone:        %d\n' % niters_per_stone)
        outf.write('thin_per_stone:          %d\n' % thin_per_stone)
        outf.write('burnin_per_stone:        %d\n' % burnin_per_stone)
        outf.write('log_marglike_ss:         %.5f\n' % log_marglike_ss)
        outf.close()
    else:
        print('  Using previous SS results found in file ss-results.txt')
        log_marglike_ss = in_log_marglike_ss
        print('  Log marginal likelihood (steppingstone) = %.5f' % log_marglike_ss)
    
def doPartitionSample():
    global training_sample, estimation_sample, ntraining, T0, nestimation, T1

    ntraining = int(trainingfrac*nsamples)
    T0 = float(ntraining)
    nestimation = nsamples - ntraining
    T1 = float(nestimation)
    training_sample = sample[:ntraining]
    estimation_sample = sample[ntraining:]

    print('  Partitioning sample:')
    print('    Training sample:   %d sites from 1 to %d' % (ntraining, ntraining))
    print('    Estimation sample: %d sites from %d to %d' % (nestimation, ntraining+1, nsamples))

def doTransformSample():
    global transformed_training, transformed_estimation
    
    print('  Transforming samples')

    # Transform the training sample
    transformed_training = []
    meanv_training, meank_training, transformed_training = transformSample(training_sample, True)

    print('    Training sample:')
    print('      Mean edge length: %.5f' % meanv_training)
    print('      Mean kappa: %.5f' % meank_training)
    
    transformed_training.sort()
    cutoff = int(math.floor(lop_off*ntraining))

    # Transform the estimation sample
    transformed_estimation = []
    meanv_estim, meank_estim, transformed_estimation = transformSample(estimation_sample, False)

    print('    Estimation sample:')
    print('      Mean edge length: %.5f' % meanv_estim)
    print('      Mean kappa: %.5f' % meank_estim)
    
def calculateNorms():
    global rmax, lorad_index
    
    # Determine working parameter space and its partition
    # A small easily-understood example:
    #
    # ntraining = 20  (training sample size)
    # coverage = 0.8 (i.e. 16 samples will be retained)
    #
    #  index     norm    edgelen
    #      0  0.02605    0.33865 <--+ <-- smallest norm
    #      1  0.05481    0.25779    |
    #      2  0.15008    0.46267    |
    #      3  0.17257    0.14002    |
    #      4  0.26212    0.57472    |
    #      5  0.35802    0.67061    |
    #      6  0.38802    0.70062    |
    #      7  0.44639    0.75899    | The 16 sampled points (80% of T=20) having the
    #      8  0.52640    0.83899    | lowest norm will be retained
    #      9  0.53953   -0.22694    |
    #     10  0.63665    0.94924    |
    #     11  0.72527   -0.41268    |
    #     12  0.86608   -0.55348    |
    #     13  1.00532   -0.69273    |
    #     14  1.15997   -0.84737    |
    #     15  1.28285   -0.97026    |
    #     16  1.38335   -1.07075 <--+ <-- lorad_index = 16
    #     17  1.48334   -1.17075
    #     18  1.87464   -1.56204
    #     19  2.21317   -1.90058 <-- largest norm
    
    lorad_index = int(coverage*ntraining) 

    # Calculate norms of training sample
    allnorms = []
    for i in range(ntraining):
        if do_jc:
            norm = math.fabs(transformed_training[i][1])
        else:
            norm = math.sqrt(math.pow(transformed_training[i][1],2.) + math.pow(transformed_training[i][2],2.))
        allnorms.append(norm)
    allnorms.sort()
    rmax = allnorms[lorad_index] 
    print('  rmax = %g' % rmax)
    
def testEdmundson():
    if not test_edmundson:
        return 
        
    # generate test_sample_size standard normal deviates and compare the fraction having
    # norm less than or equal to rmax to delta
    ndarts = test_sample_size
    if do_jc:
        mvnorm_sample = scipy.stats.multivariate_normal.rvs(mean=0, cov=1, size=ndarts)
    else:
        mvnorm_sample = scipy.stats.multivariate_normal.rvs(mean=[0.,0.], cov=[[1.,0.],[0.,1.]], size=ndarts)
    xmean = 0.0
    ymean = 0.0
    xmin = None
    ymin = None
    xmax = None
    ymax = None
    for s in mvnorm_sample:
        x = s[0]
        xmean += x
        if xmin is None or x < xmin:
            xmin = x
        if xmax is None or x > xmax:
            xmax = x
        if not do_jc:
            y = s[1]
            ymean += y
            if ymin is None or y < ymin:
                ymin = y
            if ymax is None or y > ymax:
                ymax = y
    xmean /= ndarts
    if not do_jc:
        ymean /= ndarts
    ninside = 0
    for s in mvnorm_sample:
        if do_jc:
            norm = math.fabs(s[0] - xmean)
        else:
            norm = math.sqrt(math.pow(s[0] - xmean,2.) + math.pow(s[1] - ymean,2.))
        if norm < rmax:
            ninside += 1
    print('Test of Edmundson formula:')
    print('  sample size = %d' % ndarts)
    print('  xmean       = %.5f' % xmean)
    if not do_jc:
        print('  ymean       = %.5f' % ymean)
    print('  test  = %.5f' % (float(ninside)/ndarts,))

    # Create two-dimensional plot showing first 50 sampled points, with arrows 
    # indicated values inside the rmax circle
    if not do_jc:
        minval = xmin < ymin and xmin or ymin
        maxval = xmax > ymax and xmax or ymax
        plotNorms("norm-plot.pdf", minval, maxval, minval, maxval, rmax, mvnorm_sample[:50])
    sys.exit('debug abort.')
    
def calcDelta():
    global delta, p
    
    # Find the cumulative probability delta for multivariate standard normal radial error 
    # from 0 to rmax using the formula at the very bottom of p. 11 in:
    # Edmundson, HP. 1961. The distribution of radial error and its statistical 
    # application in war gaming. Operations Research 9(1):8-21.
    #
    # scipy.special.gammainc(a,x) returns (int_0^x t^{a-1} e^-1 dt)/Gamma(a)
    #
    # Note: Edmundson's formula is incorrect: the factor 2 should not be there
    # and has been eliminated in the code below.
    p = 2.0
    if do_jc:
        p = 1.0
    delta = scipy.special.gammainc(p/2., pow(rmax,2.)/2.)
    print('  Delta = %g' % delta)

def doLoRaD():
    global log_marglike_lorad
    
    print('\nLoRaD analysis:')
    print('  training fraction (psi): %.5f' % trainingfrac)
    print('  coverage (phi): %.5f' % coverage)
    
    doPartitionSample()
    doTransformSample()
    calculateNorms()
    calcDelta()
    testEdmundson()
    
    # Compute sum of ratios used in the LoRaD method
    log_normalizing_constant = math.log(2.*math.pi)*p/2.
    log_ratio_vect = []
    for i in range(nestimation):
        if do_jc:
            norm = math.fabs(transformed_estimation[i][1])
        else:
            norm = math.sqrt(math.pow(transformed_estimation[i][1],2.) + math.pow(transformed_estimation[i][2],2.))
        if norm < rmax:
            # logh is multivariate standard normal density
            logh = -pow(norm,2.)/2. - log_normalizing_constant
    
            # logq is log-kernel
            logq = transformed_estimation[i][0]
    
            log_ratio_vect.append(logh - logq)

    max_log_ratio = max(log_ratio_vect)
    sum_terms = sum([math.exp(log_ratio - max_log_ratio) for log_ratio in log_ratio_vect])
    log_numerator = max_log_ratio + math.log(sum_terms) - math.log(T1)
    log_marglike_lorad = math.log(delta) - log_numerator
    print('  log marginal likelihood = %.5f' % log_marglike_lorad)
    
def doPlots():
    
    if do_jc:
        axis_min = -rmax
        axis_max =   rmax
    else:
        axis_min = min([-rmax, -rmax])
        axis_max = max([ rmax,  rmax])
        
    test_color = [[0, "rgb(166,206,227,.5)"], [1, "rgb(166,206,227,.5)"]]
    cylinder_color = [[0, "rgb(166,206,227)"], [1, "rgb(166,206,227)"]]
    aerobie_color  = [[0, "rgb(31,120,180)"], [1, "rgb(31,120,180)"]]
    zmax = -math.log(2.*math.pi)

    ################## log scale below here ######################        

    if False and do_plots:
        # Transformed but NOT standardized posterior on log scale (rainbow)
        # Multivariate standard normal on log scale (monochrome)
        linear_scale = False
        indep_color_scales = False
        plotSurfaces('linear-untransformed-comparison.png', linear_scale, indep_color_scales,
            axis_min, axis_max, axis_min, axis_max, 
            ['transformed-unstandardized-posterior','mvstdnorm'], 
            ['portland', cylinder_color],   # color schemes
            [log_marglike_lorad, 0.0],      # normalizing constants (specify on log scale)
            [1.,1.])                        # opacity (1 = opaque, 0 = transparent)

    if True and do_plots: 
        # Figure 1b in the paper
        # Transformed and standardized posterior on log scale (rainbow)
        linear_scale = False
        indep_color_scales = False
        vmin = -2.6
        vmax = 2.6
        vticks = [-2, -1, 0, 1, 2]
        kmin = -2.6
        kmax = 2.6
        kticks = [-2, -1, 0, 1, 2]
        zsep = 2.0
        fig,zmax = plotSurfaces('transformed-standardized.png', linear_scale, indep_color_scales,
            vmin, vmax, kmin, kmax, zsep,
            vticks, kticks,
            -1.8, 1.8, 0.5,
            ['transformed-standardized-posterior'], 
            ['portland'],          # color scheme
            [log_marglike_lorad],  # normalizing constant (specify on log scale)
            [1.])                  # opacity (1 = opaque, 0 = transparent)
        
    if True and do_plots:
        # Figure 2 in the paper
        # Transformed and standardized posterior on log scale (rainbow)
        # Multivariate standard normal on log scale (monochrome)
        linear_scale = False
        indep_color_scales = False
        vmin = -2.6
        vmax = 2.6
        vticks = [-2, -1, 0, 1, 2]
        kmin = -2.6
        kmax = 2.6
        kticks = [-2, -1, 0, 1, 2]
        zsep = 2.0
        fig,zmax = plotSurfaces('transformed-standardized-mvstdnorm.png', linear_scale, indep_color_scales,
            vmin, vmax, kmin, kmax, zsep,
            vticks, kticks,
            -1.8, 1.8, 0.5,
            ['transformed-standardized-posterior','mvstdnorm'], 
            ['portland', cylinder_color],   # color schemes
            [log_marglike_lorad, 0.0],      # normalizing constants (specify on log scale)
            [1.,.5])                        # opacity (1 = opaque, 0 = transparent)
        
    ################## linear scale below here ######################        

    if False and do_plots:
        # Transformed but NOT standardized posterior on linear scale (rainbow)
        # Multivariate normal (fit to posterior) on linear scale (rainbow)
        vaxis_min =  0.0
        vaxis_max = -3.0
        kaxis_min =  0.0
        kaxis_max =  3.0
        linear_scale = True
        indep_color_scales = True
        h2 = plotSurfaces('transformed-only-mvnorm.png', linear_scale, indep_color_scales, 
            vaxis_min, vaxis_max, kaxis_min, kaxis_max, 
            ['transformed-unstandardized-posterior','mvnorm'], 
            ['portland',test_color],    # color schemes
            [log_marglike, 0.0],        # normalizing constants (specify on log scale)
            [1.,1.])                    # opacity (1 = opaque, 0 = transparent)

    if False and do_plots:
        # Transformed but NOT standardized posterior on linear scale (rainbow)
        # Multivariate standard normal on linear scale (rainbow)
        #
        # Illustrates that doing log-transformation along is not enough to make the
        # posterior density similar to a multivariate standard normal (MVSN) distribution.
        # the variances are much smaller than the MVSN and the means are not at zero.
        axis_min = -3.0
        axis_max = 3.0
        linear_scale = True
        indep_color_scales = True
        h2 = plotSurfaces('transformed-only-mvnorm.png', linear_scale, indep_color_scales, 
            axis_min, axis_max, axis_min, axis_max, 
            ['transformed-unstandardized-posterior','mvstdnorm'], 
            ['portland', 'portland'],   # color schemes
            [log_marglike, 0.0],        # normalizing constants (specify on log scale)
            [1.,1.])                    # opacity (1 = opaque, 0 = transparent)

    if True and do_plots:
        # Figure 1a in the paper
        plotfn = 'untransformed-posterior.png'
        print('  File "%s" plots untransformed posterior surface on linear scale (rainbow):' % plotfn)
        vmin = 0.0
        vmax = 1.0
        vticks = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
        kmin = 0.0
        kmax = 8.0
        kticks = [1,2,3,4,5,6,7]
        zsep = 2.0
        linear_scale = True
        indep_color_scales = False
        fig,zmax = plotSurfaces(plotfn, linear_scale, indep_color_scales,
            vmin, vmax, kmin, kmax, zsep,
            vticks, kticks,
            -1.8, 1.8, 0.5,
            ['untransformed-posterior'], # plot types (see createZArray)
            ['portland'],                # color schemes
            [log_marglike_lorad],        # normalizing constants (specify on log scale)
            [1.0])                       # opacity (1 = opaque, 0 = transparent)
        #fig.write_json('junk.json', pretty=True)

    if False and do_plots:
        plotfn = 'transformed-and-standardized.png'
        print('  File "%s" plots two surfaces:' % plotfn)
        print('    (1) Transformed and standardized posterior on linear scale (rainbow)')
        print('    (2) Multivariate standard normal on linear scale (monochrome)')
        #
        # Standardizing as well as log-transforming causes the posterior density to be
        # nearly equal to the MVSN (multivariate standard normal) distribution. The
        # MVSN is shown in monochrome and semi-transparent to allow it to be 
        # distinguished from the posterior density
        axis_min = -3.0
        axis_max = 3.0
        zsep = 0.0
        linear_scale = True
        indep_color_scales = False
        plotSurfaces(plotfn, linear_scale, indep_color_scales,
            axis_min, axis_max, axis_min, axis_max, zsep,
            'auto','auto',
            -1.8, 1.8, 1.0,
            ['transformed-standardized-posterior','mvstdnorm'], # plot types (see createZArray)
            [cylinder_color,                       'portland'], # color schemes
            [log_marglike_lorad,                          0.0], # normalizing constants (specify on log scale)
            [1.0,                                         0.5]) # opacity (1 = opaque, 0 = transparent)

##########################################################################################
############################## main program starts here ##################################
##########################################################################################
doSimulationMCMC()
doSS()
doLoRaD()
doPlots()
