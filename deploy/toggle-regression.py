import os,re,glob

unpart  = {'name':'unpart',  'contents':open(os.path.join('unpart','lorad','lorad.conf'),'r').read()}
bygene  = {'name':'bygene',  'contents':open(os.path.join('bygene','lorad','lorad.conf'),'r').read()}
bycodon = {'name':'bycodon', 'contents':open(os.path.join('bycodon','lorad','lorad.conf'),'r').read()}
byboth  = {'name':'byboth',  'contents':open(os.path.join('byboth','lorad','lorad.conf'),'r').read()}
for stuff in [unpart, bygene, bycodon, byboth]:
    print('working on "%s"' % stuff['name'])
    m = re.search('skipmcmc\s+=\s+(no|yes)', stuff['contents'], re.S | re.M)
    assert m is not None, 'could not parse skipmcmc'
    skipmcmc = m.group(1)
    if skipmcmc == 'no':
        stuff['contents'] = re.sub('skipmcmc\s+=\s+no','skipmcmc = yes', stuff['contents'])
    
    m = re.search('useregression\s+=\s+(no|yes)', stuff['contents'], re.S | re.M)
    assert m is not None, 'could not parse useregression'
    useregression = m.group(1)
    if useregression == 'yes':
        print('  toggling useregression from "yes" to "no"')
        stuff['contents'] = re.sub('useregression\s+=\s+yes','useregression = no', stuff['contents'])
    else:
        print('  toggling useregression from "no" to "yes"')
        stuff['contents'] = re.sub('useregression\s+=\s+no','useregression = yes', stuff['contents'])

    f = open(os.path.join(stuff['name'],'lorad','lorad.conf'),'w')
    f.write(stuff['contents'])
    f.close()

    filenames = glob.glob(os.path.join(stuff['name'], 'lorad', '*.out'))
    for fn in filenames:
        os.remove(fn)
        
    filenames = glob.glob(os.path.join(stuff['name'], 'lorad', '*.err'))
    for fn in filenames:
        os.remove(fn)
