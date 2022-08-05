import subprocess
import multiprocessing
from lorad import LoRaD
import sys

assert len(sys.argv) == 2, 'Specify number of cores on command line'

ntasks = 30
ncores = int(sys.argv[1])
tasks_per_core = ntasks//ncores

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

    rfile.write('plot(x, y, type="l", lwd=2, main="%s", xlab="", ylab="")\n' % fnprefix)
    rfile.write('dev.off()\n')
    rfile.close()
    
    subprocess.Popen(['Rscript', fn]).communicate()

fnprefix = 'unpart-lorad-logtransformed-params'
verbose = False
phivect = [1.0*(i+1)/ntasks for i in range(ntasks)]
fitvect = []
deltaphivect = []
KLvect = []
use_npz = False

corebill = [list() for i in range(ncores)]

for i,phi in enumerate(phivect):
    task = i + 1
    core = i//tasks_per_core
    corebill[core].append((i,phi))
    
print('\n')
print('%12s  %s' % ('core','phi values handled'))
for i,iphi in enumerate(corebill):
    print('%12d  [%s]' % (i,', '.join(['%.3f' % phi for (i,phi) in iphi])))
    
def doLoRaD(iphi, lock):
    for i,phi in iphi:
        r = LoRaD(False, fnprefix, phi, verbose)
        lock.acquire()
        outf = open('doof.txt', 'a')
        outf.write('%12.3f %12d %12.5f %12.5f %12.5f\n' % (r['phi'],r['n'],r['fit'],r['KL'],r['logc']))
        outf.close()
        lock.release()

open('tmp.txt', 'w').close()

if ncores == 1:
    use_npz = False
    outf = open('tmp.txt', 'w')
    for i in range(ntasks):
        r = LoRaD(use_npz, fnprefix, phi, verbose)
        outf.write('%12.3f %12d %12.5f %12.5f %12.5f\n' % (r['phi'],r['n'],r['fit'],r['KL'],r['logc']))
        use_npz = True
    outf.close()
else:
    jobs = []
    lock = multiprocessing.Lock()
    for i in range(ncores):
        process = multiprocessing.Process(
            target=doLoRaD, 
            args=(corebill[i],lock)
        )
        jobs.append(process)
    
    print('starting jobs...')
    for j in jobs:
        j.start()
    
    print('joining jobs...')
    for j in jobs:
        j.join()
    
lines = open('tmp.txt', 'r').readlines()
results = {}
for line in lines:
    stripped = line.strip()
    parts = stripped.split()
    assert len(parts) == 5
    phi  = float(parts[0])
    n    = float(parts[1])
    fit  = float(parts[2])
    KL   = float(parts[3])
    logc = float(parts[4])
    results[phi] = (n, fit, KL, logc)

print('finishing...')
print('%12s %12s %12s' % ('phi','fit','logc'))
for i,phi in enumerate(phivect):
    print('%12.3f %12.5f %12.5f' % (phi, results[phi]['fit'],results[phi]['logc']))
    fitvect.append(results[phi]['fit'])
    KLvect.append(results[phi]['KL'])
    deltaphivect.append(results[phi]['KL'] - results[phi]['fit'])
Rplot('fit-plot', phivect, fitvect)
Rplot('KL-plot', phivect, KLvect)
Rplot('Delta-over-phi-plot', phivect, deltaphivect)

