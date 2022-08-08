import subprocess
from lorad import LoRaD
import sys

assert len(sys.argv) == 2, 'Specify fnprefix on command line'

ntasks = 31
fnprefix = sys.argv[1]

def Rplot(fnprefix, xvect, yvect):
    fn = '%s.R' % fnprefix
    rfile = open(fn, 'w')
    rfile.write('cwd = system(\'cd "$( dirname "$0" )" && pwd\', intern = TRUE)\n')
    rfile.write('setwd(cwd)\n')
    rfile.write('pdf("%s.pdf")\n' % fnprefix)

    xstr = ['%g' % x for x in xvect]
    rfile.write('x = c(%s)\n' % ','.join(xstr))

    ystr = ['%g' % y for y in yvect]
    rfile.write('y = c(%s)\n' % ','.join(ystr))

    rfile.write('plot(x, y, type="l", lwd=2, main="%s", xlab="", ylab="", ylim=c(0,max(y)))\n' % fnprefix)
    rfile.write('dev.off()\n')
    rfile.close()
    
    subprocess.Popen(['Rscript', fn]).communicate()

verbose = False
phivect = [1.0*(i+1)/ntasks for i in range(ntasks)]
fitvect = []
deltaphivect = []
KLvect = []

suffix = ntasks == 1 and '' or 's'
print('Evaluating %d task%s (phi value%s)' % (ntasks, suffix, suffix))
print('First phi value (%.5f) will be handled separately to create the npz file' % phivect[0])
    
results = {}

# First call to LoRaD creates npz file
phi = phivect[0]
phikey = '%.5f' % phi
print('working on phi = %s...' % phikey)
r = LoRaD(False, fnprefix, phi, verbose)    
results[phikey] = r

for i,phi in enumerate(phivect[1:-1]):
    task = i + 1
    phikey = '%.5f' % phi
    print('working on phi = %s...' % phikey)
    r = LoRaD(True, fnprefix, phi, verbose)
    results[phikey] = r

print('Summarizing...')
print('%12s %12s %12s' % ('phi','fit','logc'))
for i,phi in enumerate(phivect[:-1]):
    phikey = '%.5f' % phi
    print('%12s %12.5f %12.5f' % (phikey, results[phikey]['fit'],results[phikey]['logc']))
    fitvect.append(results[phikey]['fit'])
    KLvect.append(results[phikey]['KL'])
    deltaphivect.append(results[phikey]['KL'] - results[phikey]['fit'])
Rplot('fit-plot', phivect[:-1], fitvect)
Rplot('KL-plot', phivect[:-1], KLvect)
Rplot('Delta-over-phi-plot', phivect[:-1], deltaphivect)

