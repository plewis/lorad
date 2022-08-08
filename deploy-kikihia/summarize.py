import re, os, glob, math, sys, subprocess, scipy.stats

def roundUpToNearest(x, p):
    return p*math.floor(float(x)/p + 1.0)

def roundDownToNearest(x, p):
    return p*math.floor(float(x)/p)

def phiPlot(objective, xphi, yobj):
    # Sanity check
    assert objective in ['fit', 'KL', 'diff']
        
    for scheme in ['unpart','bycodon','bygene','byboth']:
        ymax = None
        ymin = None
        
        fnprefix = '%s-%s' % (scheme, objective)
        fn = '%s.R' % fnprefix
        rfile = open(fn, 'w')
        rfile.write('cwd = system(\'cd "$( dirname "$0" )" && pwd\', intern = TRUE)\n')
        rfile.write('setwd(cwd)\n')
        rfile.write('pdf("%s.pdf")\n' % fnprefix)

        xstr = ['%g' % x for x in xphi]
        rfile.write('x <- c(%s)\n' % ','.join(xstr))

        maxy = max(yobj[scheme][0])
        if ymax is None or maxy > ymax:
            ymax = maxy
        miny = min(yobj[scheme][0])
        if ymin is None or miny > ymin:
            ymin = miny
        ystr = ['%g' % y for y in yobj[scheme][0]]
        rfile.write('y1 <- c(%s)\n' % ','.join(ystr))

        n = len(yobj[scheme])
        k = 2
        for v in yobj[scheme][1:]:
            if v is None:
                sys.exit('v is None for k = %d, scheme = "%s", objective = "%s"' % (k,scheme, objective))
            maxy = max(v)
            if ymax is None or maxy > ymax:
                ymax = maxy
            miny = min(v)
            if ymin is None or miny > ymin:
                ymin = miny
            s = ','.join(['%g' % y for y in v])
            rfile.write('y%d <- c(%s)\n' % (k,s))
            k += 1
            
        ymax = roundUpToNearest(ymax, 0.1)
        if objective == 'diff':
            ymin = roundDownToNearest(ymin, 0.1)
        else:
            ymin = 0.0

        rfile.write('plot(x, y1, type="l", lwd=1, main="%s", xlab="", ylab="", ylim=c(%g,%g))\n' % (fnprefix,ymin,ymax))
        for k in range(2,n+1):
            rfile.write('lines(x, y%d, type="l", lwd=1)\n' % k)
        
        rfile.write('dev.off()\n')
        rfile.close()
    
        subprocess.Popen(['Rscript', fn]).communicate()

def extractXYFromRPlot(rfile):
    stuff = open(rfile, 'r').read()

    regex = 'x = c\(([-.,0-9e]+?)\)'
    m = re.search(regex, stuff, re.S | re.M)
    assert m is not None, 'could not match "%s" in "%s"' % (regex, rfile)
    xvect = None
    xvect = [float(s) for s in m.group(1).split(',')]

    regex = 'y = c\(([-.,0-9e]+?)\)'
    m = re.search(regex, stuff, re.S | re.M)
    assert m is not None, 'could not match "%s" in "%s"' % (regex, rfile)
    yvect = None
    yvect = [float(s) for s in m.group(1).split(',')]

    return xvect, yvect
    
def checkPhi(x,scheme):
    global x_phi
    if x_phi is None:
        x_phi = x
    else:
        n = len(x)
        nphi = len(x_phi)
        if nphi != n:
            print('x_phi vector has length %d which conflicts with phi values from scheme "%s" (%d values):' % (nphi,scheme,n))

def processDir(dirname, scheme):
    global x_phi, y_KL, y_fit, y_diff
    
    # Extract x and y values from fit plot file (x is phi, y is fit)
    subdir = os.path.join(dirname, scheme, 'lorad')

    x,yfit = extractXYFromRPlot(os.path.join(subdir, 'fit-plot.R'))
    checkPhi(x,scheme)
    y_fit[scheme].append(yfit)

    x,yKL = extractXYFromRPlot(os.path.join(subdir, 'KL-plot.R'))
    checkPhi(x,scheme)
    y_KL[scheme].append(yKL)
    
    x,ydiff = extractXYFromRPlot(os.path.join(subdir, 'Delta-over-phi-plot.R'))
    checkPhi(x,scheme)
    y_diff[scheme].append(ydiff)
    
    # Choose best fitting phi
    nphi = len(x)
    best_phi = None
    for i in range(nphi):
        phi = x[i]
        if use_fit:
            obj = math.fabs(yfit[i] - 1.0) # choose phi whose fit is closest to 1.0
        else:
            obj = math.fabs(yKL[i]) # choose phi whose KL divergence is closest to 0.0
        if best_phi is None or obj < best_phi[1]:
            best_phi = (phi, obj)
    
    # Extract logc for best_phi from output file
    fn = glob.glob(os.path.join(subdir, '%s-lorad-*.out' % scheme))
    stuff = open(fn[0], 'r').read()
    m = re.search('Summarizing\.\.\.\s+phi\s+fit\s+logc\s+(.+)', stuff, re.S | re.M)
    assert m is not None
    summary_lines = m.group(1).split('\n')
    best_phi_str = '%.5f' % best_phi[0]
    best_logc = None
    for line in summary_lines:
        stripped = line.strip()
        if len(stripped) > 0:
            parts = line.strip().split()
            phi = parts[0]
            if phi == best_phi_str:
                best_logc = float(parts[2])
                
    # Extract GHM estimate from output file
    m = re.search('Estimating marginal likelihood using the GHME method:\s+log Pr\(data\|focal topol\.\) = ([-.0-9]+) \(GHM estimator\)', stuff, re.S | re.M)
    assert m is not None    
    ghm = float(m.group(1))
    
    return (best_phi,best_logc,ghm)
    
if __name__ == '__main__':
    # Use fit component of KL divergence to select phi (if False, use full KL divergence)
    use_fit = False
    if len(sys.argv) > 1 and sys.argv[1] == 'fit':
        use_fit = True
        
    if use_fit:
        print('Using fit component of KL divergence to choose optimal phi')
    else:
        print('Using KL divergence to choose optimal phi')
    
    # Create list of all directory names (e.g. g1, g2, ..., g100)
    dirnames = glob.glob('g*')
    dirnames.remove('go.sh')
    dirnames.sort()
    ndirnames = len(dirnames)
    print('found %d directories to process' % ndirnames)

    unpart_logc  = []
    bycodo_logc = []
    bygene_logc  = []
    byboth_logc  = []
    unpart_ghm  = []
    bycodo_ghm = []
    bygene_ghm  = []
    byboth_ghm  = []
    x_phi = None
    y_KL   = {'unpart':[], 'bycodon':[], 'bygene':[], 'byboth':[]}
    y_fit  = {'unpart':[], 'bycodon':[], 'bygene':[], 'byboth':[]}
    y_diff = {'unpart':[], 'bycodon':[], 'bygene':[], 'byboth':[]}
    for dirname in dirnames:
        unpart_phi_logc_ghm  = processDir(dirname, 'unpart')
        unpart_logc.append(unpart_phi_logc_ghm[1])
        unpart_ghm.append(unpart_phi_logc_ghm[2])

        bycodo_phi_logc_ghm = processDir(dirname, 'bycodon')
        bycodo_logc.append(bycodo_phi_logc_ghm[1])
        bycodo_ghm.append(bycodo_phi_logc_ghm[2])

        bygene_phi_logc_ghm  = processDir(dirname, 'bygene')
        bygene_logc.append(bygene_phi_logc_ghm[1])
        bygene_ghm.append(bygene_phi_logc_ghm[2])

        byboth_phi_logc_ghm  = processDir(dirname, 'byboth')
        byboth_logc.append(byboth_phi_logc_ghm[1])
        byboth_ghm.append(byboth_phi_logc_ghm[2])
        
    phiPlot('fit', x_phi, y_fit)
    phiPlot('KL', x_phi, y_KL)
    phiPlot('diff', x_phi, y_diff)
    
    lorad_unpart = scipy.stats.describe(unpart_logc)
    lorad_bycodo = scipy.stats.describe(bycodo_logc)
    lorad_bygene = scipy.stats.describe(bygene_logc)
    lorad_byboth = scipy.stats.describe(byboth_logc)
    
    print('\nLoRaD summary:')
    print('      scheme         mean       stderr')
    print('      unpart %12.5f %12.5f' % (lorad_unpart.mean, math.sqrt(lorad_unpart.variance)))
    print('     bycodon %12.5f %12.5f' % (lorad_bycodo.mean, math.sqrt(lorad_bycodo.variance)))
    print('      bygene %12.5f %12.5f' % (lorad_bygene.mean, math.sqrt(lorad_bygene.variance)))
    print('      byboth %12.5f %12.5f' % (lorad_byboth.mean, math.sqrt(lorad_byboth.variance)))

    ghm_unpart = scipy.stats.describe(unpart_ghm)
    ghm_bycodo = scipy.stats.describe(bycodo_ghm)
    ghm_bygene = scipy.stats.describe(bygene_ghm)
    ghm_byboth = scipy.stats.describe(byboth_ghm)
    
    print('\nGHM summary:')
    print('      scheme         mean       stderr')
    print('      unpart %12.5f %12.5f' % (ghm_unpart.mean, math.sqrt(ghm_unpart.variance)))
    print('     bycodon %12.5f %12.5f' % (ghm_bycodo.mean, math.sqrt(ghm_bycodo.variance)))
    print('      bygene %12.5f %12.5f' % (ghm_bygene.mean, math.sqrt(ghm_bygene.variance)))
    print('      byboth %12.5f %12.5f' % (ghm_byboth.mean, math.sqrt(ghm_byboth.variance)))
    
    