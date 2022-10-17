from scipy.special import gammainc
from scipy.stats import describe,multivariate_normal,multivariate_t
from scipy.linalg import fractional_matrix_power
from scipy.optimize import minimize_scalar
import scipy.integrate as integrate
from math import sqrt,log,exp,lgamma,pow,pi,fabs
import subprocess as sub
import sys, os.path, numpy as np

debugging = False

# Note: regression currently should remain False because some parts of this script
#       currently do not take account of regression (e.g. KL divergence)
regression    = False   # adjust log ratios using polynomial regression

# if True, use multivariate t distribution as reference; if False, use multivariate normal
mvtreference = True  

# Set multivariate t degrees of freedom mvtv to None if you want it to be estimated
mvtv = 2.1

# Note: modecentering should remain False (not used in manuscript and not well tested)
modecentering = False

# If True, creates plot of logz - logq (y-axis) as a function of norm (x-axis)
logratioplot = False

# The following variable specifies how many columns are present before
# the first parameter column: e.g. skip 0:iter, 1:logL, 2:logP, 3:logJ, and 4:topology columns
first_param_index = 5       

# The following variable specifies which columns should be summed to yield the log posterior kernel
# e.g. 1=log-likelihood, 2=log-prior, 3=log-jacobian
posterior_indices = [1,2,3] 

def calcLogMultivariateNormal(Y, M, V):
    # Computes log of the multivariate normal density function
    # Y is the nparameters X nsamples data matrix
    # M is the nparameters mean vector
    # V is the variance-covariance matrix
    
    # Get dimensions
    p = len(M)
    if len(Y.shape) == 1:
        n = 1
        assert p == len(Y)
    else:
        n = Y.shape[1]
        assert p == Y.shape[0]
        
    if debugging:
        print('n=',n)
        print('p=',p)
    
    # Compute the log normalizing constant
    #Vscaled = V*pow(2.0*pi, p)
    #(sign, logdetVscaled) = np.linalg.slogdet(Vscaled)
    (sign, logdetV) = np.linalg.slogdet(V)
    assert sign == 1
    #log_normalizing_constant = -0.5*n*logdetVscaled
    log_normalizing_constant = -0.5*n*(float(p)*log(2.0*pi) + logdetV)
    
    # Center the Ys by subtracting the mean vector M from each column    
    Ycentered = Y - M  # (p x n) - (p x 1) = p x n
    
    # Compute V^{-1} 
    Vinv = fractional_matrix_power(V, -1.0)  # p x p    
    
    # Sum log kernel over all sampled vectors
    log_density = 0.0
    if n == 1:
        X = Ycentered
    
        # Multiply (1 x p) X^T by (p x p) V^{-1} to yield (1 x p) transposed vector
        YtSinv = np.dot(np.transpose(X),Vinv)
    
        # Multiply above (1 x p) transposed vector by ((p x 1) X to yield (1 x 1) scalar
        YtSinvY = np.dot(YtSinv,X)
        assert YtSinvY.shape == (), 'Expecting scalar, but got shape = (%d,%d)' % YtSinvY.shape
    
        log_density -= 0.5*YtSinvY
    else:
        for i in range(n):
            # Extract ith sampled vector from Ycentered
            X = Ycentered[:,i]
        
            # Multiply (1 x p) X^T by (p x p) V^{-1} to yield (1 x p) transposed vector
            YtSinv = np.dot(np.transpose(X),Vinv)
        
            # Multiply above (1 x p) transposed vector by ((p x 1) X to yield (1 x 1) scalar
            YtSinvY = np.dot(YtSinv,X)
            assert YtSinvY.shape == (), 'Expecting scalar, but got shape = (%d,%d)' % YtSinvY.shape
        
            log_density -= 0.5*YtSinvY

    return log_density + log_normalizing_constant

def calcLogUnnormalizedMultivariateT(Y, M, v, Sinv):
    # Computes log of the multivariate t density function
    # Y is the nparameters X nsamples data matrix
    # M is the nparameters mean vector
    # v is the degrees of freedom
    # S is the scale matrix
    
    # Get dimensions
    if len(Y.shape) == 1:
        n = 1
    else:
        n = Y.shape[1]
    p = len(M)
    assert p == Y.shape[0]
        
    # Center the Ys by subtracting the mean vector M from each column    
    Ycentered = Y - M  # (p x n) - (p x 1) = p x n
        
    # Sum log kernel over all sampled vectors
    log_density = 0.0
    for i in range(n):
        # Extract ith sampled vector from Ycentered
        X = Ycentered[:,i]
        
        # Multiply (1 x p) X^T by (p x p) S^{-1} to yield (1 x p) transposed vector
        YtSinv = np.dot(np.transpose(X),Sinv)
        
        # Multiply above (1 x p) transposed vector by ((p x 1) X to yield (1 x 1) scalar
        YtSinvY = np.dot(YtSinv,X)
        assert YtSinvY.shape == (), 'Expecting scalar, but got shape = (%d,%d)' % YtSinvY.shape
        
        log_density -= 0.5*(v + p)*log(1.0 + YtSinvY/v)
    return log_density

def calcMultivariateTConstants(n, p, v, S):
    # Computes and returns log_normalizing_constant and Sinv matrix
    # p is the number of parameters
    # v is the degrees of freedom
    # S is the scale matrix
    
    # Compute the log normalizing constant if necessary
    lognumer1 = lgamma(v + p) - log(2.0)
    logdenom1 = lgamma(0.5*v)
    logdenom2 = 0.5*p*log(v)
    logdenom3 = 0.5*p*log(pi)
    (sign, logdetS) = np.linalg.slogdet(S)
    assert sign == 1
    logdenom4 = 0.5*logdetS
    log_normalizing_constant = float(n)*(lognumer1 - logdenom1 - logdenom2 - logdenom3 - logdenom4)

    if debugging:
        print('\nInside calcMultivariateTConstants:')
        print('%12.5f n' % n)
        print('%12.5f v' % v)
        print('%12.5f p' % p)
        print('%12.5f lognumer1' % lognumer1)
        print('%12.5f logdenom1' % logdenom1)
        print('%12.5f logdenom2' % logdenom2)
        print('%12.5f logdenom3' % logdenom3)
        print('%12.5f logdenom4' % logdenom4)
        print('%12.5f log_normalizing_constant' % log_normalizing_constant)
        print('\n')
    
    # Compute S^{-1} 
    Sinv = fractional_matrix_power(S, -1.0)  # p x p    
    
    return (log_normalizing_constant, Sinv)

def calcLogMultivariateT(Y, M, v, S):
    # Computes log of the multivariate t density function
    # Y is the nparameters X nsamples data matrix
    # M is the nparameters mean vector
    # v is the degrees of freedom
    # S is the scale matrix
    
    # Get dimensions
    if len(Y.shape) == 1:
        n = 1
    else:
        n = Y.shape[1]
    p = len(M)
    assert p == Y.shape[0]
    
    if debugging:
        n = 1
    
    # Compute the log normalizing constant if necessary
    lognumer1 = lgamma(v + p) - log(2.0)
    logdenom1 = lgamma(0.5*v)
    logdenom2 = 0.5*p*log(v)
    logdenom3 = 0.5*p*log(pi)
    (sign, logdetS) = np.linalg.slogdet(S)
    assert sign == 1
    logdenom4 = 0.5*logdetS
    log_normalizing_constant = float(n)*(lognumer1 - logdenom1 - logdenom2 - logdenom3 - logdenom4)
    
    # Center the Ys by subtracting the mean vector M from each column    
    Ycentered = Y - M  # (p x n) - (p x 1) = p x n
    
    # Compute S^{-1} 
    Sinv = fractional_matrix_power(S, -1)  # p x p    
    
    # Sum log kernel over all sampled vectors
    log_density = 0.0
    for i in range(n):
        # Extract ith sampled vector from Ycentered
        X = Ycentered[:,i]
        
        # Multiply (1 x p) X^T by (p x p) S^{-1} to yield (1 x p) transposed vector
        XtSinv = np.dot(np.transpose(X),Sinv)
                
        # Multiply above (1 x p) transposed vector by ((p x 1) X to yield (1 x 1) scalar
        XtSinvX = np.dot(XtSinv,X)
        assert XtSinvX.shape == (), 'Expecting scalar, but got shape = (%d,%d)' % XtSinvX.shape
        
        log_density -= 0.5*(v + p)*log(1.0 + XtSinvX/v)
        
        if debugging:
            print('log_density using my function: %.5f' % log_density)
            print('  n             = %d' % n)
            print('  S trace       = %.5f' % (np.trace(S),))
            print('  Sinv trace    = %.5f' % (np.trace(Sinv),))
            print('  X.shape       =',X.shape)
            print('  Sinv.shape    =', Sinv.shape)
            print('  XtSinv.shape  =',XtSinv.shape)
            print('  XtSinvX.shape =',XtSinvX.shape)
            print('  log_normalizing_constant = %.5f' % log_normalizing_constant)
            break

    if debugging:
        print('\nInside calcLogMultivariateT:')
        print('%12.5f n' % n)
        print('%12.5f v' % v)
        print('%12.5f p' % p)
        print('%12.5f lognumer1' % lognumer1)
        print('%12.5f logdenom1' % logdenom1)
        print('%12.5f logdenom2' % logdenom2)
        print('%12.5f logdenom3' % logdenom3)
        print('%12.5f logdenom4' % logdenom4)
        print('%12.5f log_normalizing_constant' % log_normalizing_constant)
        print('%12.5f log_density' % log_density)
        print('%12.5f total' % (log_density + log_normalizing_constant,))
        print('\n')
    
    return log_density + log_normalizing_constant

def multivariate_t_numer_func(r, v, p):
    # r is radius (norm)
    # v is multivariate t degrees of freedom
    # p is number of parameters
    if r > 0.0:
        log_normalzing_constant = log(2.0) + lgamma(0.5*(v + p)) - lgamma(0.5*v) - lgamma(0.5*p) - 0.5*p*log(v-2.0)
        log_integrand = -0.5*(v+p)*log(1.0 + r*r/(v - 2.0)) + (p - 1.0)*log(r)
        return exp(log_normalzing_constant + log_integrand)    
    else:
        return 0.0

def multivariate_t_minimize_func(mvtv):
    m = (mvtv-2.0)/mvtv
    scale_matrix = m*np.identity(nparameters)
    mean_vector = np.zeros((nparameters,1))
    
    log_density0 = calcLogMultivariateT(Ystd, mean_vector, mvtv, scale_matrix)
    
    if debugging:
        print('mvtv: %.5f' % (mvtv,))
        print('scale_matrix trace: %.5f' % (np.trace(scale_matrix),))
        print('mean_vector sum: %.5f' % (np.sum(mean_vector),))

    mvt = multivariate_t(shape=scale_matrix, df=mvtv)
    n = Ystd.shape[1]
    log_density = 0.0
    for i in range(n):
        # Extract ith sampled vector from Ycentered
        X = Ystd[:,i]
        log_density += mvt.logpdf(X)
        if debugging:
            print('log_density using scipy: %.5f' % log_density)
            break
        
    return -log_density

def calcKL(Dtrim, nsamples, phi, Delta, mvn, logc):
    # Dtrim first row is ||theta||, second row is q(theta), and the remaining rows 
    #   represent the p-dimensional theta vector
    # Delta is total probability (less than 1 because Dtrim is trimmed)
    # mvn is a standard multivariate normal of same dimension (p) as theta
    KL0 = log(Delta) - log(phi)
    
    nrows,ncols = np.shape(Dtrim)
    KL1 = 0.0
    for i in range(ncols):
        v = Dtrim[2:,i]
        logz = mvn.logpdf(v)
        logq = Dtrim[1,i] - logc
        KL1 += logq - logz
    KL1 /= nsamples
    KL1 /= phi
    
    KL = KL0 + KL1

    return (KL0,KL1,KL)

def calcLogSum(logv):
    #  Want log(a + b + c), where c is (arbitrarily) the largest value
    #  log(c (a/c + b/c + 1))
    #  logc + log(a/c + b/c + 1)
    #  logc + log(exp{loga - logc} + exp{logb - logc} + 1)
    
    # First find max element
    maxlogv = max(logv)
    
    # Now sum ratios v/maxv
    sum_of_exponentials = 0.0
    for logx in logv:
        sum_of_exponentials += exp(logx - maxlogv)
        
    # Return sum of elements v on log scale
    return maxlogv + log(sum_of_exponentials)
    
def regressLogRatiosOntoNorms(log_ratios, norms):
    # log_ratios is a vector of log(reffunc(theta)) - log(kernel(theta)) for sampled theta values
    #   number of sampled theta values included is phi*nsamples
    # norms is a vector of norms for sampled theta values
    # This function returns beta0 (intercept), beta1 (norm), and beta2 (norm^2) regression
    # coefficients allowing log_ratio to be predicted given a norm
    n = len(log_ratios)
    assert n == len(norms), 'Expecting length of log_ratios vector (%d) to be same as length of norms vector (%d) in regressLogRatiosOntoNorms' % (len(log_ratios), len(norms))
    Xmatrix = []  # matrix of dimension n X 3
    Yvector = []  # vector of dimension n
    for i in range(n):
        Xmatrix.append([1.0, norms[i], pow(norms[i],2.0)])
        Yvector.append(log_ratios[i])
    X = np.array(Xmatrix).reshape((n, 3))
    Y = np.array(Yvector)
    XtX = np.dot(np.transpose(X),X) # XtX is (3 X n)(n X 3) = (3 X 3)
    XX = np.linalg.inv(XtX)         # XX  is                  (3 X 3)
    XY = np.dot(np.transpose(X),Y)  # XY  is (3 X n)(n X 1) = (3 X 1)
    B = np.dot(XX,XY)               # B   is (3 X 3)(3 X 1) = (3 X 1)
    return tuple(B)

def logIntegralOld(rstart, rstop, tol):
    p = float(nparameters)
    result = integrate.quad(lambda r: pow(r,p-1.0)*exp(-r*r/2.0), rstart, rstop)
    #print('logIntegralOld:')
    #print('  %15.9f %15.9f' % result)
    assert result[0] > 0.0
    return log(result[0])

def logIntegralNew(beta1, beta2, rstart, rstop, tol):
    p = float(nparameters)
    result = integrate.quad(lambda r: pow(r,p-1.0)*exp(-beta1*r - beta2*r*r - r*r/2.0), rstart, rstop)
    #print('logIntegralNew:')
    #print('  %15.9f %15.9f' % result)
    assert result[0] > 0.0
    return log(result[0])
    
def kernelDifferencePlot(plotfnprefix, D, nsamples, phi):
        # odd fractions will be indicated with a red vertical dotted line
        odd  = [int(0.1*nsamples),int(0.3*nsamples),int(0.5*nsamples),int(0.7*nsamples),int(0.9*nsamples)]
        # even fractions will be indicated with a gray vertical dotted line
        even = [int(0.2*nsamples),int(0.4*nsamples),int(0.6*nsamples),int(0.8*nsamples)]
        # phi will be indicated with a thicker navy vertical solid line
        curr  = [int(phi*nsamples)]
        # the 99th percentile will be indicated with a black vertical dotted line
        xtra  = [int(0.99*nsamples)]
        # the estimated log-marginal-likelihood will be indicated with a horizontal navy dotted line
        oddvect = []
        evenvect = []
        xtravect = []
        currvect = []
        nrmvect  = []
        lnzqvect = []
        mvn = multivariate_normal(mean=np.zeros((nparameters,)), cov=np.identity(nparameters))
        for i in range(nsamples):
            nrm = D[0,i]
            nrmvect.append(nrm)
            if i in odd:
                oddvect.append(nrm)
            if i in even:
                evenvect.append(nrm)
            if i in xtra:
                xtravect.append(nrm)
            if i in curr:
                currvect.append(nrm)
            v = D[2:,i]
            logz = mvn.logpdf(v)
            logq = D[1,i]
            lnzqvect.append(logz - logq)
        maxnorm = max(nrmvect)
        rfile = open('%s.R' % plotfnprefix, 'w')
        rfile.write("cwd = system('cd \"$( dirname \"$0\" )\" && pwd', intern = TRUE)\n")
        rfile.write("setwd(cwd)\n")
        rfile.write('pdf("%s.pdf")\n' % plotfnprefix)
        rfile.write('n <- c(%s)\n' % ','.join(['%g' % nrm for nrm in nrmvect]))
        rfile.write('logzq <- c(%s)\n' % ','.join(['%g' % lnzq for lnzq in lnzqvect]))
        rfile.write('plot(n, logzq, type="p", bty="L", xlab="norm", ylab="log(z/q)", xlim=c(0,%g))\n' % maxnorm)
        for m in oddvect:
            rfile.write('abline(v=%g,lty="dotted", col="red", lwd=1)\n' % m)
        for m in evenvect:
            rfile.write('abline(v=%g,lty="dotted", col="gray", lwd=1)\n' % m)
        for m in xtravect:
            rfile.write('abline(v=%g,lty="dotted", col="black", lwd=1)\n' % m)
        for m in currvect:
            rfile.write('abline(v=%g,lty="solid", col="navy", lwd=2)\n' % m)
        rfile.write('abline(h=%g,lty="solid", col="navy", lwd=1)\n' % (-1.0*lorad_c,))
        if regression:
            rfile.write('x <- seq(0.0, %g, 0.1)\n' % maxnorm)
            rfile.write('y <- (%g) + (%g)*x + (%g)*x*x\n' % (beta0, beta1, beta2))
            rfile.write('lines(x,y,lwd=2,col="yellow")\n')
        rfile.write("dev.off()\n")

#################################################################
#################################################################
##################### LoRaD function ############################
#################################################################
#################################################################

def LoRaD(npz, fnprefix, phi, verbose):
    global Ystd, M, V, S, nparameters, Dtrim, imax, mvtv, mvtscale
    if npz:
        # Load results from previous run (coverage different, but everything else is the same)
        fn = fnprefix + '.npz'
        try:
            npzar = np.load(fn)
            P = npzar['P']
            Y = npzar['Y']
            M = npzar['M']
            W = npzar['W']
            V = npzar['V']
            (nparameters,nsamples) = np.shape(Y)
            if verbose:
                print('File "%s":' % fn)
                print('  Shape of log posterior density vector P is ',np.shape(P))
                print('  Shape of data matrix Y is ',np.shape(Y))
                print('  Shape of mean vector M is ',np.shape(M))
                print('  Shape of mode vector W is ',np.shape(W))
                print('  Shape of variance-covariance matrix V is ',np.shape(V))
        except:
            sys.exit('Could not load %s' % fn)
    else:
        # Starting from scratch
        # Open the supplied log-transformed parameter sample file
        fn = fnprefix + '.txt'
        try:
            lines = open(fn, 'r').readlines()
        except:
            sys.exit('Could not load %s' % fn)

        # Here's what the beginning of the file should look like for 32 taxa (2*32-3 = 61 edge lengths):
        # 
        # iter           logL        logP         logJ   topology     logTL   logEdgeLenProp_2 ... logEdgeLenProp_61         rAG-0          rAT-0          rCG-0         rCT-0          rGT-0          piC-0          piG-0         piT-0     ratevar-0	
        #    1   -10205.37424   186.13124   -302.26081          1   0.07456        1.093373278 ...      -0.474942723   2.065022859   -1.024028312   -0.708630807   2.239913666   -1.697753626   -1.333158885   -1.364931693   0.144475022   1.896998740
        nlines = len(lines)
        nsamples = nlines - 1
        headers = lines[0].strip().split('\t')
        nheaders = len(headers)
        parameter_columns = headers[first_param_index:]
        parameters = parameter_columns[:]
        nparameters = len(parameters)
        if verbose:
            print('File "%s":' % fn)
            print('%12d lines' % nlines)
            print('%12d columns' % nheaders)
            print('%12d columns represent log-transformed parameters' % nparameters)
            print('%12d other columns are "iter", "logL", "logP", "logJ", and "topology"' % 5)
    
        # Create data structure: each row is a vector of sampled values for one parameter
        data = []
        for p in parameters:
            data.append([])
    
        # Record data
        log_posterior = []
        max_log_posterior = None
        peakj = None
        for j,line in enumerate(lines[1:]):
            # Split line into nheaders parts
            parts = line.strip().split('\t')
            assert len(parts) == len(headers)
        
            # Compute the log posterior by summing the log-likelihood (logL), log-prior (logP), and log-Jacobian (logJ)
            jth_log_posterior = sum([float(parts[k]) for k in posterior_indices])
            log_posterior.append(jth_log_posterior)

            # See if this log_posterior is the largest seen so far        
            if max_log_posterior is None or jth_log_posterior > max_log_posterior:
                max_log_posterior = jth_log_posterior
                peakj = j

            for i,p in enumerate(parameter_columns):
                # Grab ith log-transformed parameter value
                k = first_param_index + i
                logv = float(parts[k])
            
                # Append ith log-transformed parameter value to ith parameter vector
                data[i].append(logv)
            
        # Save point closest to the posterior mode
        mode_vector = np.array([data[i][peakj] for i,p in enumerate(parameters)])
        
        # Create map with parameter names as keys and scipy.stats.describe objects as values
        d = {}
        for i,p in enumerate(parameters):
            d[p] = describe(data[i], 0, 1) # 0 is axis, 1 is df for variance calculation

        # Output summary statistics for each parameter 
        if (verbose):
            print('\nParameter statistics (note: log scale):')
            print('%12s %20s %12s %12s' % ('Row','Description','Mean','Std. Dev.'))
            for i,p in enumerate(parameters):
                mean = d[p].mean
                sd = sqrt(d[p].variance)
                print('%12s %20s %12.5f %12.5f' % (i+1, p, mean, sd))
            
        # Create log posterior density vector
        P = np.array(log_posterior).reshape((1,nsamples))
        if verbose:
            print('\nShape of log posterior density vector P is ',np.shape(P))

        # Create data matrix Y
        # Example:
        # Y = [ 1  2  3  4   5 ]  shape = (2,5)
        #     [ 6  7  8  9  10 ]  2 parameters, 5 sampled values for each
        Y = np.array(data).reshape((nparameters, nsamples))
        if verbose:
            print('\nShape of data matrix Y is ',np.shape(Y))

        # Create mean vector M
        # Example:
        # M = [3] shape = (2,1)
        #     [8] axis=0 yields means of columns, axis=1 yields means of rows
        M = Y.mean(axis=1).reshape((nparameters, 1))
        if verbose:
            print('\nShape of mean vector M is ',np.shape(M))
            print('Mean vector M:')
            print(M)

        # Create mode vector W
        # Example:
        W = mode_vector.reshape((nparameters, 1))
        if verbose:
            print('\nShape of mode vector W is ',np.shape(W))
            print('Mode vector W:')
            print(W)

        # Estimate variance-covariance matrix V
        # Example:
        # V = [2.5  2.5] shape = (2,2)
        #     [2.5  2.5]
        V = np.cov(Y).reshape((nparameters, nparameters))
        if verbose:
            print('\nShape of variance-covariance matrix V is ',np.shape(V))
            print('Variance-covariance matrix V:')
            print(V)
    
        # Save everything in a numpy archive named <prefix>.npz, where <prefix> is the 
        # file name prefix of the original lorad parameter file
        np.savez('%s.npz' % fnprefix, Y=Y, P=P, M=M, W=W, V=V)

    ######## could part below here also go into npz file? ########
    
    # Compute square root of V and call it S
    S = fractional_matrix_power(V, 0.5)
    if verbose:
        print('\nShape of standard deviation matrix S is ',np.shape(S))

    # Compute inverse square root of V and call it Sinv
    Sinv = fractional_matrix_power(V, -0.5)
    if verbose:
        print('Shape of inverse standard deviation matrix S is ',np.shape(Sinv))

    # The Jacobian for the standardization transformation is det(S)
    (sign, logdetS) = np.linalg.slogdet(S)
    np.set_printoptions(precision=3, suppress=True) # suppress scientific notation
    #if sign != 1:
    #    print('sign = %g' % sign)
    #assert sign == 1
    if verbose:
        print('log(det(S)) = ',logdetS)

    # Standardize Y by subtracting M (or W) and premultiplying by Sinv
    if modecentering:
        Ycentered = Y - W
    else:
        Ycentered = Y - M
    Ystd = np.dot(Sinv, Ycentered)
    if verbose:
        print('\nShape of standardized data matrix Ystd is ',np.shape(Ystd))
        
    if mvtreference:
        if mvtv is None:
            # determine degrees of freedom mvtv for the multivariate t distribution
            # scale matrix equals V with each element multiplied by (mvtv-2)/mvtv
            #logfmvt = calcLogMultivariateNormal(Ystd, np.zeros_like(M), np.identity(nparameters))
            #print('log-density at infinity = %12.5f' % logfmvt)
            result = minimize_scalar(multivariate_t_minimize_func, bounds=[2.1,10000000.0], method="bounded")
            print('log-density at optimum  = %.3f' % (-result.fun,))
            print('optimal mvtv            = %.3f' % result.x)
            mvtv = result.x
        mvtscale = (mvtv-2.0)*np.identity(nparameters)/mvtv

    # Construct matrix D in which first row is norm, second row is P + logdetS, and the remaining rows come from Ystd
    Nstd = np.linalg.norm(Ystd, axis=0).reshape((1,nsamples))
    #print('shape of Nstd is ', np.shape(Nstd))
    Pstd = P + logdetS
    D = np.concatenate((Nstd,Pstd,Ystd)).reshape(nparameters+2,nsamples)
    #print('shape of D before sorting is ', np.shape(D))

    if verbose:
        print('\nSaving file "D.txt" (matrix in which 1st row is norm, 2nd row is log-posterior, and remaining rows are from Ystd)')
        np.savetxt('D.txt', D, fmt='%.6f', delimiter='\t', newline='\n')

    # Sort D so that columns with smallest norm are first
    sort_indices = D[0,].argsort() # new column indices are based on sort_indices
    D = D[:, sort_indices]         # calculated for the first (0th) row
    if verbose:
        print('\nMatrix D sorted so that columns with smallest norm are first')
    #print(D)
    #print('shape of D after sorting is ', np.shape(D))

    ######## could part above here also go into npz file? ########

    # Determine rmax: the norm of the (phi)th quantile and trim D
    imax = int(phi*nsamples)
    Dtrim = D[:,:imax]
    rmax = D[0,imax-1]
    if verbose:
        #print('imax = %d' % imax)
        #print('Dtrim:')
        #print(Dtrim)
        #print('shape of Dtrim is ', np.shape(Dtrim))
        print('\nNorm (radius) to %dth sampled point is %g' % (imax,rmax))
        print('The value %d equals the number of samples (%d) times the specified coverage fraction (%g)' % (imax,nsamples, phi))
        print('The value %g defines the limits of the working parameter space' % rmax)

    # Determine Delta, the volume under the reference density
    # within a distance rmax of the mean
    if mvtreference:
        Delta_numer = integrate.quad(multivariate_t_numer_func, 0.0, rmax, args=(mvtv, nparameters))
        Delta_denom = integrate.quad(multivariate_t_numer_func, 0.0, np.inf, args=(mvtv, nparameters))
        logDelta = log(Delta_numer[0]) - log(Delta_denom[0])
        Delta = exp(logDelta)
        if verbose:
            print('\nThe volume (Delta) under the multivariate standard t')
            print('reference function within the working parameter space is:')
            print('  mvt d.f.    =', mvtv)
            #print('  nparameters =', nparameters)
            #print('  Delta_numer =', Delta_numer)
            #print('  Delta_denom =', Delta_denom)
            #print('  logDelta    =', logDelta)
            print('  Delta       =', Delta)
    else:
        Delta = gammainc(nparameters/2.0, rmax*rmax/2.0)
        logDelta = log(Delta)
        if verbose:
            print('\nThe volume (Delta) under the multivariate standard normal')
            print('reference function within the working parameter space is:')
            print('  Delta = %g' % Delta)

    # Compute log ratios
    log_ratios = []
    norms = []
    if mvtreference:
        # temp
        mvt = multivariate_t(shape=mvtscale, df=mvtv) 
        mvn = multivariate_normal(mean=np.zeros((nparameters,)), cov=np.identity(nparameters)) 
        fitmvnorm = []
        fitmvt = []
        fitmvt2 = []
        
        log_normalizing_const, Sinv = calcMultivariateTConstants(1, nparameters, mvtv, mvtscale)
        for i in range(imax):
            #print('Computing ratio %d of %d...' % (i,imax))
            ith_norm = Dtrim[0,i]
            norms.append(ith_norm)
            v = Dtrim[2:,i]
            #logz = log_normalizing_const + calcLogUnnormalizedMultivariateT(v, np.zeros_like(M), mvtv, Sinv)
            logz = calcLogMultivariateT(v, np.zeros_like(M), mvtv, mvtscale)

            # temp            
            logz2 = mvt.logpdf(v)
            logzmvnorm = mvn.logpdf(v)

            logq = Dtrim[1,i]
            log_ratios.append(logz2 - logq)
            
            # temp
            fitmvnorm.append(fabs(logzmvnorm - logq))
            fitmvt.append(fabs(logz - logq))
            fitmvt2.append(fabs(logz2 - logq))
        
        # temp
        descnorm = describe(fitmvnorm, 0, 1)  # 2nd arg 0 is axis, 3rd arg 1 is df
        print('##### mean difference for multivariate normal: %.5f' % descnorm.mean)
        desct = describe(fitmvt, 0, 1)  # 2nd arg 0 is axis, 3rd arg 1 is df
        print('##### mean difference for multivariate t: %.5f' % desct.mean)
        desct2 = describe(fitmvt2, 0, 1)  # 2nd arg 0 is axis, 3rd arg 1 is df
        print('##### mean difference for multivariate t (using scipy): %.5f' % desct2.mean)
        
    else:
        mvn = multivariate_normal(mean=np.zeros((nparameters,)), cov=np.identity(nparameters))
        for i in range(imax):
            ith_norm = Dtrim[0,i]
            norms.append(ith_norm)
            v = Dtrim[2:,i]
            logz = mvn.logpdf(v)
            logq = Dtrim[1,i]
            log_ratios.append(logz - logq)
    
    if regression:
        # Estimate regression coefficients beta0, beta1, and beta2
        beta0, beta1, beta2 = regressLogRatiosOntoNorms(log_ratios, norms)
        if verbose:
            print('\nRegression of log-ratios onto norms:')
            print('  beta0 = %g' % beta0)
            print('  beta1 = %g' % beta1)
            print('  beta2 = %g' % beta2)
    
        # Modify the entries of the log_ratios vector by subtracting predicted log ratio
        orig_log_ratios = log_ratios[:]
        n = len(log_ratios)
        for i in range(n):
            r = norms[i]
            log_ratios[i] -= beta0 + beta1*r + beta2*r*r
    
        # Modify log(Delta) by swapping integrals
        # Main idea is to replace log(q(theta_i)) with log(q*(theta_i)) = log(q(theta-i)) - (beta0 + beta1*r_i + beta2*r_i*r_i)
        # This straightens out the log(z/q) vs. norm plot and creates a reference function that fits
        # the standardized posterior kernel as closely as possible. This change in the reference function means
        # that the volume Delta must be adjusted. This is done by dividing by the old volume and multiplying
        # by the new volume.
        # see /Users/plewis/Documents/manuscripts/wang/meetings/2021-11-11/proposed_density_method.pdf for details
        log_old_integral = logIntegralOld(0.0, rmax, 0.000001)
        log_new_integral = -beta0 + logIntegralNew(beta1, beta2, 0.0, rmax, 0.000001)
        logDelta += log_new_integral
        logDelta -= log_old_integral
    
        denom_log_sum = calcLogSum(log_ratios)
    else:
        denom_log_sum = calcLogSum(log_ratios)

    lorad_c = logDelta - denom_log_sum + log(nsamples)

    # Compute KL divergence between truncated standard normal and Dtrim
    # TODO: currently ignores regression
    KL0 = 0.0
    KL1 = 0.0
    KL = 0.0
    #KL0,KL1,KL = calcKL(Dtrim, nsamples, phi, Delta, mvn, lorad_c)

    # Create plot of logz - logq (y-axis) as a function of norm (x-axis)
    # The hope is that the difference between the kernel and reference will be close to constant
    #if logratioplot and not npz:
    #    kernelDifferencePlot('kernel-difference-plot', D, nsamples, phi)
        
    return {'phi':phi, 'n':imax, 'fit':KL1, 'KL':KL, 'logc':lorad_c}

#################################################################
#################################################################
####################### main loop ###############################
#################################################################
#################################################################

if __name__ == '__main__':
    if False:
        debugging = False
        print('Compare normal and t distributions as function of r:')
        r = 0.0
        mvtv0 = 436.455
        mvtv1 = 2.1
        mvtv2 = 5.0
        nparameters = 180
        
        rvect   = []
        ft0vect = []
        ft1vect = []
        ft2vect = []
        fnvect  = []
        fnastvect  = []
        
        nuptor  = 0.0
        tuptor0 = 0.0
        tuptor1 = 0.0
        tuptor2 = 0.0
        
        rincr = 0.01
        ymin = None
        ymax = None
        
        print('%12s %12s %12s %12s %12s' % ('r', 'fn', 'ft0', 'ft1', 'ft2'))
        while r < 20.0:
            nuptor_prev = nuptor
            nuptor = gammainc(1.0*nparameters/2.0, r*r/2.0)
            fn = (nuptor - nuptor_prev)/rincr
            
            tuptor_prev = tuptor0
            tuptor0 = integrate.quad(multivariate_t_numer_func, 0.0, r, args=(mvtv0, float(nparameters)))[0]
            ft0 = (tuptor0 - tuptor_prev)/rincr
            
            tuptor_prev = tuptor1
            tuptor1 = integrate.quad(multivariate_t_numer_func, 0.0, r, args=(mvtv1, float(nparameters)))[0]
            ft1 = (tuptor1 - tuptor_prev)/rincr
            
            tuptor_prev = tuptor2
            tuptor2 = integrate.quad(multivariate_t_numer_func, 0.0, r, args=(mvtv2, float(nparameters)))[0]
            ft2 = (tuptor2 - tuptor_prev)/rincr
            
            invrincr = int(1.0/rincr)
            if int(r*invrincr) % invrincr == 0:
                print('%12.1f %12.5f %12.5f %12.5f %12.5f' % (r, fn, ft0, ft1, ft2))
            
            if ymin is None or fn < ymin:
                ymin = fn
            if ft0 < ymin:
                ymin = ft0
            if ft1 < ymin:
                ymin = ft1
            if ft2 < ymin:
                ymin = ft2

            if ymax is None or fn > ymax:
                ymax = fn
            if ft0 > ymax:
                ymax = ft0
            if ft1 > ymax:
                ymax = ft1
            if ft2 > ymax:
                ymax = ft2
                
            rvect.append(r)
            fnvect.append(fn)
            ft0vect.append(ft0)
            ft1vect.append(ft1)
            ft2vect.append(ft2)
            r += rincr
            
        rfile = open('tmp.R', 'w')
        rfile.write('cwd = system(\'cd "$( dirname "$0" )" && pwd\', intern = TRUE)\n')
        rfile.write('pdf("reference-comparison.pdf")\n')
        rfile.write('r <- c(%s)\n' % ','.join(['%g' % r for r in rvect]))
        rfile.write('fn <- c(%s)\n' % ','.join(['%g' % fn for fn in fnvect]))
        rfile.write('ft0 <- c(%s)\n' % ','.join(['%g' % ft for ft in ft0vect]))
        rfile.write('ft1 <- c(%s)\n' % ','.join(['%g' % ft for ft in ft1vect]))
        rfile.write('ft2 <- c(%s)\n' % ','.join(['%g' % ft for ft in ft2vect]))
        rfile.write('plot( r, fn, type="l", lwd=2, lty="solid",  col="black",  xlab="r", ylab="PDF", ylim=c(%g,%g))\n' % (ymin, ymax))    
        rfile.write('lines(r, ft0,          lwd=2, lty="solid",  col="red")\n')    
        rfile.write('lines(r, ft1,          lwd=2, lty="dotted", col="red")\n')    
        rfile.write('lines(r, ft2,          lwd=2, lty="dotted", col="red")\n')    
        rfile.write('text(3.0, 0.32, labels="v=2.1", col="red")\n')    
        rfile.write('text(10.0, 0.15, "v=3.0", col="red")\n')    
        rfile.write('text(16.0, 0.4, "v=436.455", col="red")\n')    
        rfile.write('dev.off()\n')
        rfile.close()
        sub.Popen(['rscript', 'tmp.R']).communicate()
        sys.exit()

    assert not (regression and mvtreference), 'Cannot specify both regression and mvtreference to be True'

    nargs = len(sys.argv)
    if not nargs in [3,4]:
        print('Expecting number of arguments to be 3 or 4 (you specified %d)' % nargs)
        print('Usage: python3 lorad.py <filename> <fraction> [verbose]')
        print(' <filename> is the name of the file containing log-transformed parameter values')
        print(' <fraction> is a float between 0.0 and 1.0 specifying the fraction of the sample to retain')
        print(' verbose is an optional literal string that, if specified, will result in copious output')
        sys.exit()

    # Grab command line arguments    
    fnpath = sys.argv[1]
    phi = float(sys.argv[2])
    assert phi > 0.0 and phi < 1.0
    verbose = False
    if nargs == 4:
        assert sys.argv[3] == 'verbose', 'Expecting 4th command line argument to be equal to "verbose"'
        verbose = True

    # Split fnpath into a directory path and file name
    dir_file = os.path.split(fnpath)
    assert len(dir_file) == 2
    fdir = dir_file[0]
    fn = dir_file[1]

    # Now split fn into a root file name and an extension
    root_ext = os.path.splitext(fn)
    assert len(root_ext) == 2
    fnprefix = root_ext[0]
    fnext = root_ext[1]

    npz = (fnext == '.npz') and True or False
    
    results = LoRaD(npz, fnprefix, phi, verbose)
    print('%20.5f = fraction of sample used' % results['phi'])
    print('%20d   = sample points used' % results['n'])
    print('%20.5f = fit to reference (closest to 1.0 best)' % results['fit'])
    print('%20.5f = KL divergence' % results['KL'])
    print('%20.5f = log(marginal likelihood)' % results['logc'])

    if False:
        #debugging = True
        print('Plotting sample as rug on top of density plot...')

        log_marg_like = -9884.0
        mvtv = 436.455
        nparameters = 180

        step = imax//100
        scale_matrix = ((mvtv-2.0)/mvtv)*np.identity(nparameters)
        log_normalizing_constant, scale_matrix_inv = calcMultivariateTConstants(1, nparameters, mvtv, scale_matrix)
        mean_vector = np.zeros_like(M)

        rvect = []
        mvtvect = []
        mvt2vect = []
        mvnvect = []
        normvect = []
        postvect = []
        
        ymin = None
        ymax = None

        for i in range(0, imax, step):
            norm = Dtrim[0,i]
            post = Dtrim[1,i] - log_marg_like
            Y = Dtrim[2:,i]
            tpdf = calcLogUnnormalizedMultivariateT(Y, mean_vector, mvtv, scale_matrix_inv)
            tpdf += log_normalizing_constant
            tpdf2 = calcLogMultivariateT(Y, mean_vector, mvtv, scale_matrix)
            npdf = calcLogMultivariateNormal(Y, mean_vector, np.identity(nparameters))
            
            
            if ymax is None or tpdf > ymax:
                ymax = tpdf
            if tpdf2 > ymax:
                ymax = tpdf2
            if npdf > ymax:
                ymax = npdf
                
            if ymin is None or tpdf < ymin:
                ymin = tpdf
            if tpdf2 < ymin:
                ymin = tpdf2
            if npdf < ymin:
                ymin = npdf

            normvect.append(norm)
            postvect.append(post)
            mvtvect.append(tpdf)
            mvt2vect.append(tpdf2)
            mvnvect.append(npdf)
            
        rfile = open('tmp.R', 'w')
        rfile.write('cwd = system(\'cd "$( dirname "$0" )" && pwd\', intern = TRUE)\n')
        rfile.write('pdf("rugplot.pdf")\n')
        rfile.write('mvtpdf <- c(%s)\n' % ','.join(['%g' % mvt for mvt in mvtvect]))
        rfile.write('mvtpdf2 <- c(%s)\n' % ','.join(['%g' % mvt for mvt in mvt2vect]))
        rfile.write('mvnpdf <- c(%s)\n' % ','.join(['%g' % mvn for mvn in mvnvect]))
        rfile.write('norms <- c(%s)\n' % ','.join(['%g' % n for n in normvect]))
        #rfile.write('posts <- c(%s)\n' % ','.join(['%g' % p for p in postvect]))
        rfile.write('plot(norms, posts, type="p", col="gray", xlab="r", ylab="logpdf", ylim=c(%g,%g))\n' % (ymin,ymax))
        rfile.write('points(norms, mvtpdf, col="red")\n')    
        rfile.write('points(norms, mvtpdf2, col="orange")\n')    
        rfile.write('points(norms, mvnpdf, col="blue")\n')    
        rfile.write('dev.off()\n')
        rfile.close()
        sub.Popen(['rscript', 'tmp.R']).communicate()
        sys.exit()

    if False:
        debugging = False
        print('Plot samples as rug beneath candidate reference marginal densities:')
        print('  Multivariate-t degrees of freedom: %g' % mvtv)
        r = 0.0
        mvtv0 = mvtv
        mvtv1 = 2.1
        mvtv2 = 5.0
        nparameters = 180
        
        rvect   = []
        ft0vect = []
        ft1vect = []
        ft2vect = []
        fnvect  = []
        fnastvect  = []
        
        nuptor  = 0.0
        tuptor0 = 0.0
        tuptor1 = 0.0
        tuptor2 = 0.0
        
        rincr = 0.01
        ymin = None
        ymax = None
        
        print('%12s %12s %12s %12s %12s' % ('r', 'fn', 'ft0', 'ft1', 'ft2'))
        while r < 20.0:
            nuptor_prev = nuptor
            nuptor = gammainc(1.0*nparameters/2.0, r*r/2.0)
            fn = (nuptor - nuptor_prev)/rincr
            
            tuptor_prev = tuptor0
            tuptor0 = integrate.quad(multivariate_t_numer_func, 0.0, r, args=(mvtv0, float(nparameters)))[0]
            ft0 = (tuptor0 - tuptor_prev)/rincr
            
            tuptor_prev = tuptor1
            tuptor1 = integrate.quad(multivariate_t_numer_func, 0.0, r, args=(mvtv1, float(nparameters)))[0]
            ft1 = (tuptor1 - tuptor_prev)/rincr
            
            tuptor_prev = tuptor2
            tuptor2 = integrate.quad(multivariate_t_numer_func, 0.0, r, args=(mvtv2, float(nparameters)))[0]
            ft2 = (tuptor2 - tuptor_prev)/rincr
            
            invrincr = int(1.0/rincr)
            if int(r*invrincr) % invrincr == 0:
                print('%12.1f %12.5f %12.5f %12.5f %12.5f' % (r, fn, ft0, ft1, ft2))
            
            if ymin is None or fn < ymin:
                ymin = fn
            if ft0 < ymin:
                ymin = ft0
            if ft1 < ymin:
                ymin = ft1
            if ft2 < ymin:
                ymin = ft2

            if ymax is None or fn > ymax:
                ymax = fn
            if ft0 > ymax:
                ymax = ft0
            if ft1 > ymax:
                ymax = ft1
            if ft2 > ymax:
                ymax = ft2
                
            rvect.append(r)
            fnvect.append(fn)
            ft0vect.append(ft0)
            ft1vect.append(ft1)
            ft2vect.append(ft2)
            r += rincr
        
        # Collect data for rug plot
        step = imax//100
        normvect = []
        for i in range(0, imax, step):
            norm = Dtrim[0,i]
            normvect.append(norm)
            
        rfile = open('tmp.R', 'w')
        rfile.write('cwd = system(\'cd "$( dirname "$0" )" && pwd\', intern = TRUE)\n')
        rfile.write('pdf("reference-comparison-rug.pdf")\n')
        rfile.write('r <- c(%s)\n' % ','.join(['%g' % r for r in rvect]))
        rfile.write('fn <- c(%s)\n' % ','.join(['%g' % fn for fn in fnvect]))
        rfile.write('ft0 <- c(%s)\n' % ','.join(['%g' % ft for ft in ft0vect]))
        rfile.write('ft1 <- c(%s)\n' % ','.join(['%g' % ft for ft in ft1vect]))
        rfile.write('ft2 <- c(%s)\n' % ','.join(['%g' % ft for ft in ft2vect]))
        rfile.write('norms <- c(%s)\n' % ','.join(['%g' % n for n in normvect]))
        rfile.write('plot( r, fn, type="l", lwd=2, lty="solid",  col="black",  xlab="r", ylab="PDF", ylim=c(%g,%g))\n' % (ymin, ymax))    
        rfile.write('lines(r, ft0,          lwd=2, lty="solid",  col="red")\n')    
        rfile.write('lines(r, ft1,          lwd=2, lty="dotted", col="red")\n')    
        rfile.write('lines(r, ft2,          lwd=2, lty="dotted", col="red")\n')    
        rfile.write('text(3.0, 0.32, labels="v=2.1", col="red")\n')    
        rfile.write('text(10.0, 0.15, "v=3.0", col="red")\n')    
        rfile.write('text(16.0, 0.4, "v=436.455", col="red")\n') 
        rfile.write('rug(norms, ticksize=0.03, side=1, lwd=1, col="blue")\n')   
        rfile.write('dev.off()\n')
        rfile.close()
        sub.Popen(['rscript', 'tmp.R']).communicate()
        sys.exit()
    
