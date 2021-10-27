import sys,shutil,os,glob,re

def summarizeGSS(partition_scheme):
    # find output file
    filenames = glob.glob('%s/gss/slurm-*.out' % partition_scheme)
    assert len(filenames) == 1
    fn = filenames[0]

    # read output file 
    stuff = open(fn, 'r').read()

    # get seed
    m = re.search('Pseudorandom number seed: (\d+)', stuff, re.M | re.S)
    assert m is not None, 'could not extract seed from file "%s"' % fn
    rnseed = int(m.group(1))

    # grab marginal likelihood estimate
    m = re.search('log\(marginal likelihood\) = ([-.0-9]+)', stuff, re.M | re.S)
    assert m is not None, 'could not find marginal likelihood in file "%s"' % fn
    logL = float(m.group(1))

    # grab times
    m = re.search('user-seconds\s+(.+)', stuff, re.M | re.S)
    assert m is not None, 'could not find user-seconds in file "%s"' % fn
    secs = float(m.group(1))
    
    return {'rnseed':rnseed, 'secs':secs, 'logL':logL}

def summarizeLoRaD(partition_scheme):
    # find output file
    filenames = glob.glob('%s/lorad/slurm-*.out' % partition_scheme)
    assert len(filenames) == 1
    fn = filenames[0]

    # read output file 
    stuff = open(fn, 'r').read()

    # get seed
    m = re.search('Pseudorandom number seed: (\d+)', stuff, re.M | re.S)
    assert m is not None, 'could not extract seed from file "%s"' % fn
    rnseed = int(m.group(1))

    # grab marginal likelihood estimate for coverage value 1
    findall = re.findall(' Determining working parameter space for coverage = ([.0-9]+?)[.][.][.].+?log\(marginal likelihood\) = ([-.0-9]+)', stuff, re.M | re.S)
    assert len(findall) == 3
    cov1 = float(findall[0][0])
    cov2 = float(findall[1][0])
    cov3 = float(findall[2][0])
    logL1 = float(findall[0][1])
    logL2 = float(findall[1][1])
    logL3 = float(findall[2][1])

    # grab times
    m = re.search('user-seconds\s+(.+)', stuff, re.M | re.S)
    assert m is not None, 'could not find user-seconds in file "%s"' % fn
    secs = float(m.group(1))
    
    return {'rnseed':rnseed, 'secs':secs, 'cov1':cov1, 'cov2':cov2, 'cov3':cov3, 'logL1':logL1, 'logL2':logL2, 'logL3':logL3}
    
lorad = {}
lorad['unpart']  = summarizeLoRaD('unpart')
lorad['bygene']  = summarizeLoRaD('bygene')
lorad['bycodon'] = summarizeLoRaD('bycodon')
lorad['byboth']  = summarizeLoRaD('byboth')

gss = {}
gss['unpart']  = summarizeGSS('unpart')
gss['bygene']  = summarizeGSS('bygene')
gss['bycodon'] = summarizeGSS('bycodon')
gss['byboth']  = summarizeGSS('byboth')

assert gss['unpart']['rnseed']  == lorad['unpart']['rnseed']
assert gss['bygene']['rnseed']  == lorad['bygene']['rnseed']
assert gss['bycodon']['rnseed'] == lorad['bycodon']['rnseed']
assert gss['byboth']['rnseed']  == lorad['byboth']['rnseed']

outf = open('summary.txt','w')
outf.write('partition\tseed\tsecs\tcov1\tcov2\tcov3\tlorad1\tlorad2\tlorad3\tsecs\tgss\n')
outf.write('unpart\t%d\t%.3f\t%.1f\t%.1f\t%.1f\t%.5f\t%.5f\t%.5f\t%.3f\t%.5f\n'  % (lorad['unpart']['rnseed'],lorad['unpart']['secs'],lorad['unpart']['cov1'],lorad['unpart']['cov2'],lorad['unpart']['cov3'],lorad['unpart']['logL1'],lorad['unpart']['logL2'],lorad['unpart']['logL3'],gss['unpart']['secs'],gss['unpart']['logL']))
outf.write('bygene\t%d\t%.3f\t%.1f\t%.1f\t%.1f\t%.5f\t%.5f\t%.5f\t%.3f\t%.5f\n'  % (lorad['bygene']['rnseed'],lorad['bygene']['secs'],lorad['bygene']['cov1'],lorad['bygene']['cov2'],lorad['bygene']['cov3'],lorad['bygene']['logL1'],lorad['bygene']['logL2'],lorad['bygene']['logL3'],gss['bygene']['secs'],gss['bygene']['logL']))
outf.write('bycodon\t%d\t%.3f\t%.1f\t%.1f\t%.1f\t%.5f\t%.5f\t%.5f\t%.3f\t%.5f\n' % (lorad['bycodon']['rnseed'],lorad['bycodon']['secs'],lorad['bycodon']['cov1'],lorad['bycodon']['cov2'],lorad['bycodon']['cov3'],lorad['bycodon']['logL1'],lorad['bycodon']['logL2'],lorad['bycodon']['logL3'],gss['bycodon']['secs'],gss['bycodon']['logL']))
outf.write('byboth\t%d\t%.3f\t%.1f\t%.1f\t%.1f\t%.5f\t%.5f\t%.5f\t%.3f\t%.5f\n'  % (lorad['byboth']['rnseed'],lorad['byboth']['secs'],lorad['byboth']['cov1'],lorad['byboth']['cov2'],lorad['byboth']['cov3'],lorad['byboth']['logL1'],lorad['byboth']['logL2'],lorad['byboth']['logL3'],gss['byboth']['secs'],gss['byboth']['logL']))
outf.close()
