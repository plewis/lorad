import sys,shutil,os,glob,re

use_comma_separated_values = False

def summarizeLoRaD(model):
    lc_model = model.lower()
    
    # set default values in case there are no results yet for this analysis
    rnseed   = 0
    cov1     = 0.0
    logL1    = 0.0
    cov2     = 0.0
    logL2    = 0.0
    cov3     = 0.0
    logL3    = 0.0
    
    # read output files
    outfilenames = glob.glob('%s/frt-%s.out' % (model,lc_model))
    errfilenames = glob.glob('%s/frt-%s.err' % (model,lc_model))
    if len(outfilenames) == 1 and len(errfilenames) == 1:
        print('%s...' % lc_model)
        outfn = outfilenames[0]
        errfn = errfilenames[0]

        # read output and error files
        outstuff = open(outfn, 'r').read()
        errstuff = open(errfn, 'r').read()
        stuff = outstuff + errstuff

        # get seed
        m = re.search('Pseudorandom number seed: (\d+)', stuff, re.M | re.S)
        if m is not None:
            rnseed = int(m.group(1))
    
        # grab times
        #m = re.search('user-seconds\s+([.0-9]+)', stuff, re.M | re.S)
        #if m is not None:
        #    secs = float(m.group(1))

        # grab marginal likelihood estimate for each of the three coverage values
        results = re.findall(' Determining working parameter space for coverage = ([.0-9]+?)[.][.][.].+?log Pr\(data\)\s+=\s+([-.0-9]+)', stuff, re.M | re.S)
        nresults = len(results)               
        if nresults == 3:
            cov1     = float(results[0][0])
            cov2     = float(results[1][0])
            cov3     = float(results[2][0])
        
            logL1    = float(results[0][1])
            logL2    = float(results[1][1])
            logL3    = float(results[2][1])
        else:
            print('  nresults was %d (expecting 3) so did not process' % nresults)
    else:
        print('%s: Did not process because there were %d outfilenames and %d errfilenames' % (lc_model,len(outfilenames),len(errfilenames)))
        
    return {
        'rnseed':rnseed, 
        'cov1':cov1, 
        'logL1':logL1, 
        'cov2':cov2, 
        'logL2':logL2, 
        'cov3':cov3, 
        'logL3':logL3
    }

models = ['JC', 'JCI', 'JCG', 'JCIG', 'GTR', 'GTRI', 'GTRG', 'GTRIG', '3JC', '3JCI', '3JCG', '3JCIG', '3GTR', '3GTRI', '3GTRG', '3GTRIG']

lorad = {}
for m in models:
    lorad[m] = summarizeLoRaD(m)

gss = {}
gss['JC']     =  -2776.52 
gss['JCI']    =  -2744.59 
gss['JCG']    =  -2747.44 
gss['JCIG']   =  -2743.56 
gss['GTR']    =  -2714.20 
gss['GTRI']   =  -2681.00 
gss['GTRG']   =  -2682.73 
gss['GTRIG']  =  -2680.29 
gss['3JC']    =  -2681.79 
gss['3JCI']   =  -2668.38 
gss['3JCG']   =  -2668.99 
gss['3JCIG']  =  -2667.19 
gss['3GTR']   =  -2551.10 
gss['3GTRI']  =  -2535.57 
gss['3GTRG']  =  -2536.75 
gss['3GTRIG'] =  -2534.66 

outf = open('output-summary.txt','w')
if use_comma_separated_values:
    outf.write('model,seed,gss,cov1,lorad1,diff1,cov2,lorad2,diff2,cov3,lorad3,diff3\n')
    for m in models:
        outf.write('%s,%d,%.2f,%.3f,%.5f,%.5f,%.3f,%.5f,%.5f,%.3f,%.5f,%.5f\n' % (
            m, 
            lorad[m]['rnseed'],
            gss[m],
            lorad[m]['cov1'],
            lorad[m]['logL1'],
            lorad[m]['logL1'] - gss[m],
            lorad[m]['cov2'],
            lorad[m]['logL2'],
            lorad[m]['logL2'] - gss[m],
            lorad[m]['cov3'],
            lorad[m]['logL3'],
            lorad[m]['logL3'] - gss[m]
        ))
else:
    outf.write('model\tseed\tgss\tcov1\tlorad1\tdiff1\tcov2\tlorad2\tdiff2\tcov3\tlorad3\tdiff3\n')
    for m in models:
        outf.write('%s\t%d\t%.2f\t%.3f\t%.5f\t%.5f\t%.3f\t%.5f\t%.5f\t%.3f\t%.5f\t%.5f\n' % (
            m, 
            lorad[m]['rnseed'],
            gss[m],
            lorad[m]['cov1'],
            lorad[m]['logL1'],
            lorad[m]['logL1'] - gss[m],
            lorad[m]['cov2'],
            lorad[m]['logL2'],
            lorad[m]['logL2'] - gss[m],
            lorad[m]['cov3'],
            lorad[m]['logL3'],
            lorad[m]['logL3'] - gss[m]
        ))
outf.close()
