import os,re,glob

stuff = open('lorad.conf','r').read()
m = re.search('skipmcmc\s+=\s+(no|yes)', stuff, re.S | re.M)
assert m is not None, 'could not parse skipmcmc'
skipmcmc = m.group(1)
if skipmcmc == 'no':
    stuff = re.sub('skipmcmc\s+=\s+no','skipmcmc = yes', stuff)

m = re.search('useregression\s+=\s+(no|yes)', stuff, re.S | re.M)
assert m is not None, 'could not parse useregression'
useregression = m.group(1)
if useregression == 'yes':
    print('  toggling useregression from "yes" to "no"')
    stuff = re.sub('useregression\s+=\s+yes','useregression = no', stuff)
else:
    print('  toggling useregression from "no" to "yes"')
    stuff = re.sub('useregression\s+=\s+no','useregression = yes', stuff)

f = open('lorad.conf','w')
f.write(stuff)
f.close()
