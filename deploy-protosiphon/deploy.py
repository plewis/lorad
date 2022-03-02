import sys,os,re,shutil

email = {}
email['pol02003'] = 'paul.o.lewis@gmail.com'
email['aam21005'] = 'analisa.milkey@uconn.edu'

userid                 = 'pol02003'
dest_dir_prefix        = 'g'                       # prefix of name of directory to be created 
dest_dir_index         = 1                         # appended to dest_dir_prefix (e.g. 'g1' if dest_dir_prefix='g' and dest_dir_index=1)
rnseed                 = '12345'                   # the pseudorandom number seed to use for all analyses

nargs = len(sys.argv)
if nargs == 3:
    dest_dir_index = int(sys.argv[1])
    rnseed         = sys.argv[2]
    print('Setting dest_dir_index to %d and rnseed to %s from your command line input' % (dest_dir_index, rnseed))
elif nargs == 1:
    print('Using default values of dest_dir_index (%d) and rnseed (%s)' % (dest_dir_index, rnseed))
else:
    print('You must specify either 0 or 2 extra command line arguments, you specified %d' % nargs)
    print('Usage: python %d [<index> <rnseed>]')
    print('where <index> is integer index appended to destination directory name (e.g. 1 yields "g1"')
    print('and <rnseed> is the integer pseudorandom number seed to use for the entire set of analyses.')
    sys.exit('Aborting. Please try again with correct number of command line arguments.')

# Usage:
#   python3 deploy.py
#
# Creates the following directory structure:
# <dest_dir_prefix><dest_dir_index>
# 3GTR
#     go.sh
#     lorad.conf
#     s.sh
# 3GTRG
#     go.sh
#     lorad.conf
#     s.sh
# 3GTRI
#     go.sh
#     lorad.conf
#     s.sh
# 3GTRIG
#     go.sh
#     lorad.conf
#     s.sh
# 3JC
#     go.sh
#     lorad.conf
#     s.sh
# 3JCG
#     go.sh
#     lorad.conf
#     s.sh
# 3JCI
#     go.sh
#     lorad.conf
#     s.sh
# 3JCIG
#     go.sh
#     lorad.conf
#     s.sh
# common
#     FRT2000rbcL.nex
#     frtmle.tre
# GTR
#     go.sh
#     lorad.conf
#     s.sh
# GTRG
#     go.sh
#     lorad.conf
#     s.sh
# GTRI
#     go.sh
#     lorad.conf
#     s.sh
# GTRIG
#     go.sh
#     lorad.conf
#     s.sh
# JC
#     go.sh
#     lorad.conf
#     s.sh
# JCG
#     go.sh
#     lorad.conf
#     s.sh
# JCI
#     go.sh
#     lorad.conf
#     s.sh
# JCIG
#     go.sh
#     lorad.conf
#     s.sh
# submit-all.sh
# summary.py

models = [
    '3GTR',
    '3GTRG',
    '3GTRI',
    '3GTRIG',
    '3JC',
    '3JCG',
    '3JCI',
    '3JCIG',
    'GTR',
    'GTRG',
    'GTRI',
    'GTRIG',
    'JC',
    'JCG',
    'JCI',
    'JCIG'
]

# Create destination directory
dest_dir = '%s%d' % (dest_dir_prefix, dest_dir_index)
if os.path.exists(dest_dir):
    sys.exit('destination directory (%s) exists; please rename, delete, or move it and try again' % dest_dir)
os.mkdir(dest_dir)

for m in models:
    mlc = m.lower()
    
    # Create directory for model m
    m_dir = os.path.join(dest_dir, m)
    os.mkdir(m_dir)
    
    # Create go.sh file
    go_contents = open('go-template.sh', 'r').read()
    f = open(os.path.join(m_dir,'go.sh'),'w')
    go_contents = re.sub('__PREFIX__', m.lower(), go_contents)
    f.write(go_contents)
    f.close()

    # Create s.sh file
    s_contents = open('slurm-template.sh', 'r').read()
    f = open(os.path.join(m_dir,'s.sh'),'w')
    s_contents = re.sub('__PREFIX__', mlc, s_contents)
    f.write(s_contents)
    f.close()
    
    # Create lorad.conf file
    shutil.copyfile('lorad-%s-template.conf' % mlc, os.path.join(m_dir, 'lorad.conf'))    

# Create directory for data and tree files
data_dir = os.path.join(dest_dir, 'common')
os.mkdir(data_dir)

# Copy data file to data_dir
shutil.copyfile('FRT2000rbcL.nex', os.path.join(dest_dir, 'common', 'FRT2000rbcL.nex'))

# Copy tree file to data_dir
shutil.copyfile('frtmle.tre', os.path.join(dest_dir, 'common', 'frtmle.tre'))

# Copy main lorad.conf file to dest_dir
shutil.copyfile('lorad-template.conf', os.path.join(dest_dir, 'lorad.conf'))

# Copy summary.py file to dest_dir
shutil.copyfile('summary.py', os.path.join(dest_dir, 'summary.py'))

# Copy submit-all.sh file to dest_dir
shutil.copyfile('submit-all-template.sh', os.path.join(dest_dir, 'submit-all.sh'))

