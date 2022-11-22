import re, os, glob, math, sys, subprocess, scipy.stats

#exclude_dirs = ['g8','g9','g10','g11','g12','g13','g14','g15','g16','g17','g18','g19']
exclude_dirs = []

include_unpart  = True
include_bycodon = True
include_bygene  = True
include_byboth  = True

include_gss     = False
include_lorad   = True
include_ghm     = False
include_rev     = False

gss_no_output = []
gss_not_finished = [] 

def roundUpToNearest(x, p):
    return p*math.floor(float(x)/p + 1.0)

def roundDownToNearest(x, p):
    return p*math.floor(float(x)/p)

def useFit():
    # See if user has supplied the keyword "fit" on the command line
    # If so, use fit component of KL divergence to select phi 
    # Otherwise, use full KL divergence to select phi
    use_fit = False
    if len(sys.argv) > 1 and sys.argv[1] == 'fit':
        use_fit = True
        
    if use_fit:
        print('Using fit component of KL divergence to choose optimal phi')
    else:
        print('Using KL divergence to choose optimal phi')
        
    return use_fit
    
def phiPlot(objective, xphi, yobj):
    # Sanity check
    assert objective in ['fit', 'KL', 'diff']
        
    for scheme in ['unpart','bycodon','bygene','byboth']:
        if scheme == 'unpart' and not include_unpart:
            continue
        elif scheme == 'bycodon' and not include_unpart:
            continue
        elif scheme == 'bygene' and not include_bygene:
            continue
        elif scheme == 'byboth' and not include_byboth:
            continue
        
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

def processGSSDir(dirname, scheme): 
    global gss_no_output, gss_not_finished   
    # Extract x and y values from fit plot file (x is phi, y is fit)
    subdir = os.path.join(dirname, scheme, 'gss')
    outfname = glob.glob(os.path.join(subdir, '*.out'))
    if len(outfname) == 0:
        gss_no_output.append(dirname)
        return (0.0,False)
    assert len(outfname) == 1, 'There are %d output file names that fit glob pattern "%s", expecting just 1' % (len(outfname),os.path.join(subdir, '*.out'))
    stuff = open(outfname[0], 'r').read()
    m = re.search('log\(marginal likelihood\) = ([-0-9.]+)', stuff, re.S | re.M)
    if m is None:
        gss_not_finished.append(dirname)
        return (0.0,False)
    logL = float(m.group(1))
    return (logL,True)
    
def processLoRaDDir(dirname, scheme):
    global x_phi, y_KL, y_fit, y_diff
    
    # Extract x and y values from fit plot file (x is phi, y is fit)
    subdir = os.path.join(dirname, scheme, 'lorad')
    
    # Here are the files that need to exist before we can proceed
    fit_plot_R            = os.path.join(subdir, 'fit-plot.R')
    KL_plot_R             = os.path.join(subdir, 'KL-plot.R')
    delta_over_phi_plot_R = os.path.join(subdir, 'Delta-over-phi-plot.R')
    output_file_glob      = glob.glob(os.path.join(subdir, '%s-lorad-*.out' % scheme))
    output_file           = None
    if len(output_file_glob) == 1:
        output_file = output_file_glob[0]
    ok = (output_file is not None) and os.path.exists(output_file)
    ok = ok and os.path.exists(fit_plot_R)
    ok = ok and os.path.exists(KL_plot_R)
    ok = ok and os.path.exists(delta_over_phi_plot_R)
    if not ok:
        return {'phi':0.0, 'logc':0.0, 'ghm':0.0, 'ok':False}

    x,yfit = extractXYFromRPlot(fit_plot_R)
    checkPhi(x,scheme)
    y_fit[scheme].append(yfit)

    x,yKL = extractXYFromRPlot(KL_plot_R)
    checkPhi(x,scheme)
    y_KL[scheme].append(yKL)
    
    x,ydiff = extractXYFromRPlot(delta_over_phi_plot_R)
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
    stuff = open(output_file, 'r').read()
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
    
    return {'phi':best_phi, 'logc':best_logc, 'ghm':ghm, 'ok':True}
    
def summarizeLoRaD():
    global x_phi, y_KL, y_fit, y_diff, use_fit
    use_fit = useFit()

    unpart_logc  = []
    bycodon_logc = []
    bygene_logc  = []
    byboth_logc  = []
    unpart_ghm  = []
    bycodon_ghm = []
    bygene_ghm  = []
    byboth_ghm  = []
    x_phi = None
    y_KL   = {'unpart':[], 'bycodon':[], 'bygene':[], 'byboth':[]}
    y_fit  = {'unpart':[], 'bycodon':[], 'bygene':[], 'byboth':[]}
    y_diff = {'unpart':[], 'bycodon':[], 'bygene':[], 'byboth':[]}
    for dirname in dirnames:
        if include_unpart:
            unpart_phi_logc_ghm  = processLoRaDDir(dirname, 'unpart')
            if unpart_phi_logc_ghm['ok']:
                unpart_logc.append(unpart_phi_logc_ghm['logc'])
                unpart_ghm.append(unpart_phi_logc_ghm['ghm'])

        if include_bycodon:
            bycodon_phi_logc_ghm = processLoRaDDir(dirname, 'bycodon')
            if bycodon_phi_logc_ghm['ok']:
                bycodon_logc.append(bycodon_phi_logc_ghm['logc'])
                bycodon_ghm.append(bycodon_phi_logc_ghm['ghm'])

        if include_bygene:
            bygene_phi_logc_ghm  = processLoRaDDir(dirname, 'bygene')
            if bygene_phi_logc_ghm['ok']:
                bygene_logc.append(bygene_phi_logc_ghm['logc'])
                bygene_ghm.append(bygene_phi_logc_ghm['ghm'])

        if include_byboth:
            byboth_phi_logc_ghm  = processLoRaDDir(dirname, 'byboth')
            if byboth_phi_logc_ghm['ok']:
                byboth_logc.append(byboth_phi_logc_ghm['logc'])
                byboth_ghm.append(byboth_phi_logc_ghm['ghm'])
        
    phiPlot('fit', x_phi, y_fit)
    phiPlot('KL', x_phi, y_KL)
    phiPlot('diff', x_phi, y_diff)
    
    if include_lorad:
        if include_unpart:
            lorad_unpart = scipy.stats.describe(unpart_logc)
        if include_bycodon:
            lorad_bycodon = scipy.stats.describe(bycodon_logc)
        if include_bygene:
            lorad_bygene = scipy.stats.describe(bygene_logc)
        if include_byboth:
            lorad_byboth = scipy.stats.describe(byboth_logc)
    
        print('\nLoRaD summary:')
        print('      scheme         mean       stderr         nobs')
        if include_unpart:
            print('      unpart %12.5f %12.5f %12d' % (lorad_unpart.mean, math.sqrt(lorad_unpart.variance), lorad_unpart.nobs))
        if include_bycodon:
            print('     bycodon %12.5f %12.5f %12d' % (lorad_bycodon.mean, math.sqrt(lorad_bycodon.variance), lorad_bycodon.nobs))
        if include_bygene:
            print('      bygene %12.5f %12.5f %12d' % (lorad_bygene.mean, math.sqrt(lorad_bygene.variance), lorad_bygene.nobs))
        if include_byboth:
            print('      byboth %12.5f %12.5f %12d' % (lorad_byboth.mean, math.sqrt(lorad_byboth.variance), lorad_byboth.nobs))

    if include_ghm:
        if include_unpart:
            ghm_unpart = scipy.stats.describe(unpart_ghm)
        if include_bycodon:
            ghm_bycodon = scipy.stats.describe(bycodon_ghm)
        if include_bygene:
            ghm_bygene = scipy.stats.describe(bygene_ghm)
        if include_byboth:
            ghm_byboth = scipy.stats.describe(byboth_ghm)
    
        print('\nGHM summary:')
        print('      scheme         mean       stderr         nobs')
        if include_unpart:
            print('      unpart %12.5f %12.5f %12d' % (ghm_unpart.mean, math.sqrt(ghm_unpart.variance), ghm_unpart.nobs))
        if include_bycodon:
            print('     bycodon %12.5f %12.5f %12d' % (ghm_bycodon.mean, math.sqrt(ghm_bycodon.variance), ghm_bycodon.nobs))
        if include_bygene:
            print('      bygene %12.5f %12.5f %12d' % (ghm_bygene.mean, math.sqrt(ghm_bygene.variance), ghm_bygene.nobs))
        if include_byboth:
            print('      byboth %12.5f %12.5f %12d' % (ghm_byboth.mean, math.sqrt(ghm_byboth.variance), ghm_byboth.nobs))

def summarizeGSS():
    unpart_gss_logc  = []
    bycodon_gss_logc = []
    bygene_gss_logc  = []
    byboth_gss_logc  = []
    for dirname in dirnames:
        if include_unpart:
            logc,ok = processGSSDir(dirname, 'unpart')
            if ok:
                unpart_gss_logc.append(logc)

        if include_bycodon:
            logc,ok = processGSSDir(dirname, 'bycodon')
            if ok:
                bycodon_gss_logc.append(logc)

        if include_bygene:
            logc,ok = processGSSDir(dirname, 'bygene')
            if ok:
                bycodon_gss_logc.append(logc)

        if include_byboth:
            logc,ok = processGSSDir(dirname, 'byboth')
            if ok:
                byboth_gss_logc.append(logc)

    if include_unpart:
        gss_unpart = scipy.stats.describe(unpart_gss_logc)
    if include_bycodon:
        gss_bycodon = scipy.stats.describe(bycodon_gss_logc)
    if include_bygene:
        gss_bygene = scipy.stats.describe(bygene_gss_logc)
    if include_byboth:
        gss_byboth = scipy.stats.describe(byboth_gss_logc)
    
    print('\nGSS summary:')
    print('      scheme         mean       stderr         nobs')
    if include_unpart:
        print('      unpart %12.5f %12.5f %12d' % (gss_unpart.mean, math.sqrt(gss_unpart.variance), gss_unpart.nobs))
    if include_bycodon:
        print('     bycodon %12.5f %12.5f %12d' % (gss_bycodon.mean, math.sqrt(gss_bycodon.variance), gss_bycodon.nobs))
    if include_bygene:
        print('      bygene %12.5f %12.5f %12d' % (gss_bygene.mean, math.sqrt(gss_bygene.variance), gss_bygene.nobs))
    if include_byboth:
        print('      byboth %12.5f %12.5f %12d' % (gss_byboth.mean, math.sqrt(gss_byboth.variance), gss_byboth.nobs))
    
if __name__ == '__main__':
    # Create list of all directory names (e.g. g1, g2, ..., g100)
    dirnames = glob.glob('g*')
    dirnames.remove('go.sh')
    for x in exclude_dirs:
        dirnames.remove(x)
    dirnames.sort()
    ndirnames = len(dirnames)
    print('found %d directories to process' % ndirnames)
    
    if includdeploy-kikihia/summarize.pye_lorad or include_ghm:
        summarizeLoRaD()
    if include_gss:
        summarizeGSS()
        
    if len(gss_no_output) > 0:
        print('\nNo output file found for these GSS analyses:')
        print(' '.join(['%s' % d[1:] for d in gss_no_output]))
            
    if len(gss_not_finished) > 0:
        print('\nThese GSS analyses did not complete:')
        print(' '.join(['%s' % d[1:] for d in gss_not_finished]))
