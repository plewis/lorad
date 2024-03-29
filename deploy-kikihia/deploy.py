import sys,os,re

include_lorad = True
include_gss   = True
include_ghm   = True
include_rev   = False

include_unpart  = True
include_bycodon = True
include_bygene  = True
include_byboth  = True

notify = False
email = {}
email['pol02003'] = 'paul.o.lewis@gmail.com'
email['aam21005'] = 'analisa.milkey@uconn.edu'

userid                 = 'pol02003'
dest_dir_prefix        = 'g'                 # prefix of name of directory to be created 
dest_dir_index         = 1                         # appended to dest_dir_prefix (e.g. 'g1' if dest_dir_prefix='g' and dest_dir_index=1)
rnseed                 = '14219'                   # the pseudorandom number seed to use for all analyses

double_data = False  # experimental, do not set to True

nargs = len(sys.argv)
if nargs == 3:
    # python3 deploy.py 2  13579
    #                   |  seed
    #                   index (i.e. will create directory g2)
    dest_dir_index = int(sys.argv[1])
    rnseed         = sys.argv[2]
    print('Setting dest_dir_index to %d and rnseed to %s from your command line input' % (dest_dir_index, rnseed))
elif nargs == 1:
    # python3 deploy.py
    print('Using default values of dest_dir_index (%d) and rnseed (%s)' % (dest_dir_index, rnseed))
else:
    print('You must specify either 0 or 2 extra command line arguments, you specified %d' % nargs)
    print('Usage: python %d [<index> <rnseed>]')
    print('where <index> is integer index appended to destination directory name (e.g. 1 yields "g1"')
    print('and <rnseed> is the integer pseudorandom number seed to use for the entire set of analyses.')
    sys.exit('Aborting. Please try again with correct number of command line arguments.')

# LoRaD/GHM method settings
lorad_burnin           = '100000'                  # the burnin used by all LoRaD analyses  (was 1110000)
lorad_niter            = '10000000'                # the number of iterations used by all LoRaD analyses
lorad_samplefreq       = '100'                     # the sampling frequency used by all LoRaD analyses
lorad_printfreq        = '100000'                  # the print frequency used by all LoRaD analyses
lorad_coverage10       = '0.1'                     # coverage probability used by all LoRaD analyses
lorad_coverage20       = '0.2'                     # coverage probability used by all LoRaD analyses
lorad_coverage30       = '0.3'                     # coverage probability used by all LoRaD analyses
lorad_coverage40       = '0.4'                     # coverage probability used by all LoRaD analyses
lorad_coverage50       = '0.5'                     # coverage probability used by all LoRaD analyses
lorad_coverage60       = '0.6'                     # coverage probability used by all LoRaD analyses
lorad_coverage70       = '0.7'                     # coverage probability used by all LoRaD analyses
lorad_coverage80       = '0.8'                     # coverage probability used by all LoRaD analyses
lorad_coverage90       = '0.9'                     # coverage probability used by all LoRaD analyses
lorad_coverage95       = '0.95'                    # coverage probability used by all LoRaD analyses
lorad_coverage99       = '0.99'                    # coverage probability used by all LoRaD analyses
lorad_simple           = 'no'                      # used when LoRaD will not use regression to tune the reference function
lorad_regression       = 'yes'                     # whether LoRaD will use regression to tune the reference function
lorad_linearregression = 'no'                      # whether LoRaD will use linear regression (alternative is polytomial)

# Generalized Steppingstone settings
gss_burnin             = '10000'                   # the burnin used by all GSS analyses 
gss_niter              = '1000000'                 # the number of itertions used by all GSS analyses
gss_samplefreq         = '100'                     # the sampling frequency used by all GSS analyses
gss_printfreq          = '10000'                   # the print frequency used by all GSS analyses
gss_nstones            = '20'                      # the number of steppingstone ratios used by all GSS analyses
gss_alpha              = '1.0'                     # the alpha value used by all GSS analyses to choose power posterior powers

# RevBayes Steppingstone settings (note: not used in Wang et al. 2021 paper)
rev_burnin             = '1000'                    # the burnin used by all SS analyses 
rev_niter              = '10000'                   # the number of iterations used by all SS analyses
rev_samplefreq         = '100'                     # the sampling frequency used by all SS analyses
rev_printfreq          = '1000'                    # the print frequency used by all SS analyses
rev_tuning_interval    = '50'                      # the tuning interval used during burnin for all SS analyses
rev_nstones            = '50'                      # the number of steppingstone ratios used by all SS analyses
rev_alpha              = '0.25'                    # the alpha value used by all SS analyses to choose power posterior powers

# True produces data used in the Fan et al. 2011 paper
# This should be False for compatibility with the Wang et al. 2022 LoRaD paper
fan_etal_2011 = False

# True uses all 32 taxa as in Fan et al. 2011
# False excludes Kikihia muta east, whick lacks data for ATPase6 and ATPase8
# This should be True for compatibility with the Wang et al. 2022 LoRaD paper
thirtytwo = True

# Usage:
#   python3 deploy.py
#
# Reads S1679.nex and creates from it 23 data files, distributed as follows within
# the specified destination directory (using dest_dir_prefix and dest_dir_index variables).
# Note that fewer directories will be created if any of the include_lorad, include_gss, 
# include_ghm, and include_ref variables are False.
# 
# Directory structure:
# <dest_dir_prefix><dest_dir_index>
#    unpart 
#       data
#           unpart.nex
#       lorad 
#           lorad.conf
#           s.sh
#       gss  
#           lorad.conf
#           s.sh
#       rev
#           ss.Rev
#           s.sh
#    bycodon
#       data
#           bycodon.nex
#           codon1st.nex
#           codon2nd.nex
#           codon3rd.nex
#       lorad 
#           lorad.conf
#           s.sh
#       gss  
#           lorad.conf
#           s.sh
#       rev
#           ss.Rev
#           s.sh
#    bygene 
#       data
#           bygene.nex
#           COI.nex
#           COII.nex
#           ATPase6.nex
#           ATPase8.nex
#       lorad 
#           lorad.conf
#           s.sh
#       gss  
#           lorad.conf
#           s.sh
#       rev
#           ss.Rev
#           s.sh
#    byboth 
#       data
#           byboth.nex
#           COI-1st.nex
#           COI-2nd.nex
#           COI-3rd.nex, 
#           COII-1st.nex
#           COII-2nd.nex
#           COII-3rd.nex, 
#           ATPase6-1st.nex
#           ATPase6-2nd.nex
#           ATPase6-3rd.nex
#           ATPase8-1st.nex
#           ATPase8-2nd.nex
#           ATPase8-3rd.nex
#       lorad 
#           lorad.conf
#           s.sh
#       gss
#           lorad.conf
#           s.sh
#       rev
#           ss.Rev
#           s.sh
#    gtrg-31taxa.tre
#    gtrg-32taxa.tre
#    submit-lorad.sh  
#    submit-gss.sh  
#    submit-rev.sh  

excluded_taxa = ['Kikihia muta east']
if thirtytwo:
    excluded_taxa = []
    
tree_file_name = 'gtrg-31taxa.tre'
if thirtytwo:
    tree_file_name = 'gtrg-32taxa.tre'

# This taxon ordering is from S1679.nex
taxa_in_order = [
  'Kikihia_acoustica',
  'Kikihia_angusta',
  'Kikihia_aotea_east',
  'Kikihia_aotea_west',
  'Kikihia_astragali',
  'Kikihia_balaena',
  'Kikihia_cauta',
  'Kikihia_convicta',
  'Kikihia_cutora_cumberi',
  'Kikihia_cutora_cutora',
  'Kikihia_cutora_exulis',
  'Kikihia_dugdalei',
  'Kikihia_flemingi',
  'Kikihia_horologium',
  'Kikihia_laneorum',
  'Kikihia_longula',
  'Kikihia_murihikua',
  'Kikihia_muta',
  'Kikihia_muta_east',
  'Kikihia_nelsonensis',
  'Kikihia_ochrina',
  'Kikihia_paxillulae',
  'Kikihia_peninsularis',
  'Kikihia_rosea',
  'Kikihia_scutellaris',
  'Kikihia_subalpina',
  'Kikihia_tasmani',
  'Kikihia_tuta',
  'Kikihia_westlandica_north',
  'Kikihia_westlandica_south',
  'Maoricicada_cassiope',
  'Rhodopsalta_microdora'
]

genes_used = ['COI','COII','ATPase8','ATPase6']
gene_boundaries = {}
gene_boundaries['COI']     = (1,774)     #  774 -    0 = 774 (258 codons)
gene_boundaries['COII']    = (775,1476)  # 1476 -  774 = 702 (234 codons)
gene_boundaries['tRNA']    = (1477,1538) # 1538 - 1476 =  62
if fan_etal_2011:
    gene_boundaries['ATPase8'] = (1539,1689) # (1477,1627) 1627 - 1476 = 151 (50.333 codons)  <-- includes AT from beginning of ATPase6 at end of ATPase8
    gene_boundaries['ATPase6'] = (1690,2152) # (1628,2090) 2090 - 1627 = 463 (154.333 codons) <-- ATPase6 starts with 3rd position A
else:
    gene_boundaries['ATPase8'] = (1539,1687) # (1477,1625) 1625 - 1476 = 149 (49.667 codons)
    gene_boundaries['ATPase6'] = (1691,2152) # (1629,2090) 2090 - 1628 = 462 (154 codons)     <-- note AT from ATPase6 and G from ATPase6 deleted (and ATPase6 lacks final 3rd position)

def readNexusFile(fn):
    '''
    Reads nexus file whose name is specified by fn and returns ntax, nchar, taxa, and a
    sequences dictionary with taxon names as keys. The values ntax and nchar are integers,
    while taxa is a list of taxon names in the order they were found in the taxa block or
    data block. Any underscores in taxon names are converted to spaces before being saved
    in the taxa list or as a key in the sequences dictionary. Also all nexus comments
    (text in square brackets) will be ignored.
    '''
    stuff = open(fn, 'r').read()
    mask = None

    # determine if taxa block exists
    taxa_block = None
    m = re.search('(?:BEGIN|Begin|begin)\s+(?:TAXA|Taxa|taxa)\s*;(.+?)(?:END|End|end)\s*;', stuff, re.M | re.S)
    if m is not None:
        taxa_block = m.group(1).strip()

    # determine if characters block exists
    characters_block = None
    m = re.search('(?:BEGIN|Begin|begin)\s+(?:CHARACTERS|Characters|characters)\s*;(.+?)(?:END|End|end)\s*;', stuff, re.M | re.S)
    if m is not None:
        characters_block = m.group(1).strip()        

    # determine if data block exists
    data_block = None
    m = re.search('(?:BEGIN|Begin|begin)\s+(?:DATA|Data|data)\s*;(.+?)(?:END|End|end)\s*;', stuff, re.M | re.S)
    if m is not None:
        data_block = m.group(1).strip()

    if data_block is not None:
        # get ntax and nchar
        m = re.search('(?:DIMENSIONS|dimensions|Dimensions)\s+(?:NTAX|ntax|Ntax|NTax)\s*=\s*(\d+)\s+(?:NCHAR|nchar|Nchar|NChar)\s*=\s*(\d+)\s*;', data_block, re.M | re.S)
        assert m, 'Could not decipher dimensions statement in data block'
        ntax = int(m.group(1))
        nchar = int(m.group(2))

        # get matrix
        m = re.search('(?:MATRIX|matrix|Matrix)\s+(.+?)\s*;', data_block, re.M | re.S)
        assert m, 'Could not decipher matrix statement in data block'
        lines = m.group(1).strip().split('\n')
        taxa = []
        sequences = {}
        for line in lines:
            m = re.match('\[([-*]+)\]', line.strip())
            if m is not None:
                mask = m.group(1)
            else:
                stripped_line = re.sub('\[.+?\]', '', line).strip()
                if len(stripped_line) > 0:
                    parts = line.split()
                    assert len(parts) == 2, 'Found more than 2 parts to this line:\n%s' % line
                    taxon_name = re.sub('_', ' ', parts[0]).strip()
                    taxa.append(taxon_name)
                    sequences[taxon_name] = parts[1]
    else:
        assert characters_block is not None and taxa_block is not None, 'Assuming nexus file contains either a data block or a taxa block and characters block'

        # strip nexus comments from characters block
        characters_block = re.sub('\[.+?\]', '', characters_block, re.S | re.M)        

        # strip TITLE and LINK commands from characters block if there is one
        characters_block = re.sub('(?:TITLE|Title|title)\s.+?;', '', characters_block, re.S | re.M)        
        characters_block = re.sub('(?:LINK|Link|link)\s.+?;', '', characters_block, re.S | re.M)
        
        # strip blank lines from characters block
        characters_block = re.sub('\n\n+', '\n', characters_block, re.S | re.M)        
        
        # strip leading and trailing whitespace
        characters_block = re.sub('^\s+', '', characters_block, re.S | re.M)
        characters_block = re.sub('\s+$', '', characters_block, re.S | re.M)

        # replace sites with multiple states with ambiguity codes
        characters_block = re.sub('{AG}', 'R', characters_block, re.S | re.M)
        characters_block = re.sub('{CT}', 'Y', characters_block, re.S | re.M)
        characters_block = re.sub('{AC}', 'M', characters_block, re.S | re.M)
        characters_block = re.sub('{GT}', 'K', characters_block, re.S | re.M)
        characters_block = re.sub('{CG}', 'S', characters_block, re.S | re.M)
        characters_block = re.sub('{AT}', 'W', characters_block, re.S | re.M)
        characters_block = re.sub('{ACT}', 'H', characters_block, re.S | re.M)
        characters_block = re.sub('{CGT}', 'B', characters_block, re.S | re.M)
        characters_block = re.sub('{ACG}', 'V', characters_block, re.S | re.M)
        characters_block = re.sub('{AGT}', 'D', characters_block, re.S | re.M)
        characters_block = re.sub('{ACGT}', 'N', characters_block, re.S | re.M)        

        # print('*************** characters block *****************')
        # print(characters_block)
        # doof = open('doof.txt','w')
        # doof.write(characters_block)
        # doof.write('\n')
        # doof.close()
        # sys.exit('debug stop')

        # get ntax from taxa block
        m = re.search('(?:DIMENSIONS|dimensions|Dimensions)\s+(?:NTAX|ntax|Ntax|NTax)\s*=\s*(\d+)\s*;', taxa_block, re.M | re.S)
        assert m, 'Could not decipher dimensions statement in taxa block'
        ntax = int(m.group(1))

        # get nchar from characters block
        m = re.search('(?:DIMENSIONS|dimensions|Dimensions)\s+(?:NCHAR|nchar|Nchar|NChar)\s*=\s*(\d+)\s*;', characters_block, re.M | re.S)
        assert m, 'Could not decipher dimensions statement in characters block'
        nchar = int(m.group(1))

        # get matrix from characters block
        m = re.search('(?:MATRIX|matrix|Matrix)\s+(.+?)\s*;', characters_block, re.M | re.S)
        assert m, 'Could not decipher matrix statement in characters block'
        lines = m.group(1).strip().split('\n')
        taxa = []
        sequences = {}
        for line in lines:
            m = re.match('\[([-*]+)\]', line.strip())
            if m is not None:
                mask = m.group(1)
            else:
                stripped_line = re.sub('\[.+?\]', '', line).strip()
                if len(stripped_line) > 0:
                    parts = stripped_line.split()
                    assert len(parts) == 2, 'Found more than 2 parts to this line:\n%s' % line
                    taxon_name = re.sub('_', ' ', parts[0]).strip()
                    taxa.append(taxon_name)
                    sequences[taxon_name] = parts[1]

    return (ntax, nchar, mask, taxa, sequences)

def writeNexusFile(fn, taxa, mask_vect, nchar_vect, sequences_vect):
    # - fn is filename to use for the nexus file saved
    # - taxa is a list of taxon names serving as indices into the maps stored in the sequences vector
    # - masks_vect is a vector holding masks for each element of sequences (a mask is a string containing
    #   '-' for every site that should be included and '*' for every site that should be excluded.
    # - nchar_vect is a vector holding the number of sites for each element of sequences
    # - sequences_vect is a vector of maps, each with taxa as keys and sequences as values
    # The file saved will be a concatenation of all elements of the sequences vector, and
    # a vector will be returned containing the index of the first site beyond the end of each subset:
    # For example: [300, 600, 900] would be returned if sequences contained 1st, 2nd, and 3rd
    # codon positions and the total gene length was 900.
    if os.path.exists(fn):
        os.rename(fn, '%s.bak' % fn)
    nsubsets = len(sequences_vect)
    assert mask_vect is None or len(mask_vect) == nsubsets
    assert len(nchar_vect) == nsubsets
    
    # compute boundaries
    boundaries = []
    ntotal = 0
    for n in nchar_vect:
        ntotal += n
        if double_data:
            ntotal += n
        boundaries.append(ntotal)
    
    ntax = len(taxa) - len(excluded_taxa)
    nchar = sum(nchar_vect)
    if double_data:
        nchar *= 2
    mask = None
    if mask_vect is not None:
        mask = ''
        for m in mask_vect:
            mask += m
    longest = max([len(t) for t in taxa])
    taxonfmt = '  %%%ds' % longest
    f = open(fn, 'w')
    f.write('#nexus\n\n')
    f.write('begin data;\n')
    f.write('  dimensions ntax=%d nchar=%d;\n' % (ntax, nchar))
    f.write('  format datatype=dna gap=- missing=?;\n')
    f.write('  matrix\n')
    if mask is not None:
        f.write(taxonfmt % ' ')
        f.write('[%s]\n' % mask)
    for t in taxa:
        if not t in excluded_taxa:
            taxon_name = re.sub('\s+', '_', t)
            f.write(taxonfmt % taxon_name)
            f.write(' ')
            for s in sequences_vect:
                f.write('%s' % s[t])
                if double_data:
                    f.write('%s' % s[t])
            f.write('\n')
    f.write('  ;\n')
    f.write('end;\n')
    f.close()
    return boundaries

ntax, nchar, mask, taxa, sequences = readNexusFile('S1679.nex')

########################           
# Create sequence maps #
########################           

unpartseqs   = {}
codon1seqs   = {}
codon2seqs   = {}
codon3seqs   = {}
COIseqs      = {}
COIIseqs     = {}
ATPase6seqs  = {}
ATPase8seqs  = {}
COIseqs1     = {}
COIseqs2     = {}
COIseqs3     = {}
COIIseqs1    = {}
COIIseqs2    = {}
COIIseqs3    = {}
ATPase6seqs1 = {}
ATPase6seqs2 = {}
ATPase6seqs3 = {}
ATPase8seqs1 = {}
ATPase8seqs2 = {}
ATPase8seqs3 = {}
nchar0 = 0
nchar1 = 0
nchar2 = 0
nchar3 = 0
ncharCOI     = 0
ncharCOII    = 0
ncharATPase6 = 0
ncharATPase8 = 0
ncharCOI1       = 0
ncharCOI2       = 0
ncharCOI3       = 0
ncharCOII1      = 0
ncharCOII2      = 0
ncharCOII3      = 0
ncharATPase61   = 0
ncharATPase62   = 0
ncharATPase63   = 0
ncharATPase81   = 0
ncharATPase82   = 0
ncharATPase83   = 0

for t in taxa:
    unpartseqs[t] = ''
    nchar0 = 0
    
    codon1seqs[t] = ''
    codon2seqs[t] = ''
    codon3seqs[t] = ''
    nchar1 = 0
    nchar2 = 0
    nchar3 = 0
    
    COIseqs[t]     = ''
    COIIseqs[t]    = ''
    ATPase6seqs[t] = ''
    ATPase8seqs[t] = ''
    ncharCOI     = 0
    ncharCOII    = 0
    ncharATPase6 = 0
    ncharATPase8 = 0

    COIseqs1[t]     = ''
    COIseqs2[t]     = ''
    COIseqs3[t]     = ''
    COIIseqs1[t]    = ''
    COIIseqs2[t]    = ''
    COIIseqs3[t]    = ''
    ATPase6seqs1[t] = ''
    ATPase6seqs2[t] = ''
    ATPase6seqs3[t] = ''
    ATPase8seqs1[t] = ''
    ATPase8seqs2[t] = ''
    ATPase8seqs3[t] = ''
    ncharCOI1       = 0
    ncharCOI2       = 0
    ncharCOI3       = 0
    ncharCOII1      = 0
    ncharCOII2      = 0
    ncharCOII3      = 0
    ncharATPase61   = 0
    ncharATPase62   = 0
    ncharATPase63   = 0
    ncharATPase81   = 0
    ncharATPase82   = 0
    ncharATPase83   = 0

    for g in genes_used:
        first,last = gene_boundaries[g]

        # j is index relative to start of gene
        j = 0
        
        # fan
        starting_codon = 1
        if fan_etal_2011 and g == 'ATPase6':
            starting_codon = 3
            
        # i is index relative to start of matrix
        for i in range(first-1,last):
            # unpart
            unpartseqs[t] += sequences[t][i]
            nchar0 += 1
        
            # by gene
            if g == 'COI':
                COIseqs[t] += sequences[t][i]
                ncharCOI += 1
            elif g == 'COII':
                COIIseqs[t] += sequences[t][i]
                ncharCOII += 1
            elif g == 'ATPase6':
                ATPase6seqs[t] += sequences[t][i]
                ncharATPase6 += 1
            else:
                assert g == 'ATPase8'
                ATPase8seqs[t] += sequences[t][i]
                ncharATPase8 += 1

            # by codon
            site = starting_codon + j
            if site % 3 == 1:
                codon1seqs[t] += sequences[t][i]
                nchar1 += 1
            elif site % 3 == 2:
                codon2seqs[t] += sequences[t][i]
                nchar2 += 1
            elif site % 3 == 0:
                codon3seqs[t] += sequences[t][i]
                nchar3 += 1
                
            # by both gene and codon
            if g == 'COI':
                if site % 3 == 1:
                    COIseqs1[t] += sequences[t][i]
                    ncharCOI1 += 1
                elif site % 3 == 2:
                    COIseqs2[t] += sequences[t][i]
                    ncharCOI2 += 1
                elif site % 3 == 0:
                    COIseqs3[t] += sequences[t][i]
                    ncharCOI3 += 1
            elif g == 'COII':
                if site % 3 == 1:
                    COIIseqs1[t] += sequences[t][i]
                    ncharCOII1 += 1
                elif site % 3 == 2:
                    COIIseqs2[t] += sequences[t][i]
                    ncharCOII2 += 1
                elif site % 3 == 0:
                    COIIseqs3[t] += sequences[t][i]
                    ncharCOII3 += 1
            elif g == 'ATPase6':
                if site % 3 == 1:
                    ATPase6seqs1[t] += sequences[t][i]
                    ncharATPase61 += 1
                elif site % 3 == 2:
                    ATPase6seqs2[t] += sequences[t][i]
                    ncharATPase62 += 1
                elif site % 3 == 0:
                    ATPase6seqs3[t] += sequences[t][i]
                    ncharATPase63 += 1
            else:
                assert g == 'ATPase8'
                if site % 3 == 1:
                    ATPase8seqs1[t] += sequences[t][i]
                    ncharATPase81 += 1
                elif site % 3 == 2:
                    ATPase8seqs2[t] += sequences[t][i]
                    ncharATPase82 += 1
                elif site % 3 == 0:
                    ATPase8seqs3[t] += sequences[t][i]
                    ncharATPase83 += 1
                
            j += 1
 
##########################
# Create directory paths #
##########################

dest_dir = '%s%d' % (dest_dir_prefix, dest_dir_index)

if include_unpart:
    unpart_dir       = os.path.join(dest_dir, 'unpart')
    unpart_data_dir  = os.path.join(unpart_dir, 'data')
    if not fan_etal_2011:
        if include_lorad or include_ghm:
            unpart_lorad_dir        = os.path.join(unpart_dir, 'lorad')
            unpart_lorad_dir_short  = os.path.join('unpart', 'lorad')
        if include_gss:
            unpart_gss_dir       = os.path.join(unpart_dir, 'gss')
            unpart_gss_dir_short = os.path.join('unpart', 'gss')
    if include_rev:
        unpart_rev_dir       = os.path.join(unpart_dir, 'rev')
        unpart_rev_dir_short = os.path.join('unpart', 'rev')

if include_bycodon:
    bycodon_dir       = os.path.join(dest_dir, 'bycodon')
    bycodon_data_dir  = os.path.join(bycodon_dir, 'data')
    if not fan_etal_2011:
        if include_lorad or include_ghm:
            bycodon_lorad_dir        = os.path.join(bycodon_dir, 'lorad')
            bycodon_lorad_dir_short  = os.path.join('bycodon', 'lorad')
        if include_gss:
            bycodon_gss_dir       = os.path.join(bycodon_dir, 'gss')
            bycodon_gss_dir_short = os.path.join('bycodon', 'gss')
    if include_rev:
        bycodon_rev_dir       = os.path.join(bycodon_dir, 'rev')
        bycodon_rev_dir_short = os.path.join('bycodon', 'rev')
       
if include_bygene:
    bygene_dir       = os.path.join(dest_dir, 'bygene')
    bygene_data_dir  = os.path.join(bygene_dir, 'data')
    if not fan_etal_2011:
        if include_lorad or include_ghm:
            bygene_lorad_dir        = os.path.join(bygene_dir, 'lorad')
            bygene_lorad_dir_short  = os.path.join('bygene', 'lorad')
        if include_gss:
            bygene_gss_dir       = os.path.join(bygene_dir, 'gss')
            bygene_gss_dir_short = os.path.join('bygene', 'gss')
    if include_rev:
        bygene_rev_dir       = os.path.join(bygene_dir, 'rev')
        bygene_rev_dir_short = os.path.join('bygene', 'rev')
       
if include_byboth:
    byboth_dir       = os.path.join(dest_dir, 'byboth')
    byboth_data_dir  = os.path.join(byboth_dir, 'data')
    if not fan_etal_2011:
        if include_lorad or include_ghm:
            byboth_lorad_dir       = os.path.join(byboth_dir, 'lorad')
            byboth_lorad_dir_short = os.path.join('byboth', 'lorad')
        if include_gss:
            byboth_gss_dir       = os.path.join(byboth_dir, 'gss')
            byboth_gss_dir_short = os.path.join('byboth', 'gss')
    if include_rev:
        byboth_rev_dir       = os.path.join(byboth_dir, 'rev')
        byboth_rev_dir_short = os.path.join('byboth', 'rev')

######################
# Create directories #
######################

if os.path.exists(dest_dir):
    sys.exit('destination directory (%s) exists; please rename, delete, or move it and try again' % dest_dir)
os.mkdir(dest_dir)

if include_unpart:
    os.mkdir(unpart_dir     )
    os.mkdir(unpart_data_dir)
    if not fan_etal_2011:
        if include_lorad or include_ghm:
            os.mkdir(unpart_lorad_dir )
        if include_gss:
            os.mkdir(unpart_gss_dir  )
    if include_rev:
        os.mkdir(unpart_rev_dir  )

if include_bycodon:
    os.mkdir(bycodon_dir     )
    os.mkdir(bycodon_data_dir)
    if not fan_etal_2011:
        if include_lorad or include_ghm:
            os.mkdir(bycodon_lorad_dir )
        if include_gss:
            os.mkdir(bycodon_gss_dir  )
    if include_rev:
        os.mkdir(bycodon_rev_dir  )

if include_bygene:
    os.mkdir(bygene_dir     )
    os.mkdir(bygene_data_dir)
    if not fan_etal_2011:
        if include_lorad or include_ghm:
            os.mkdir(bygene_lorad_dir )
        if include_gss:
            os.mkdir(bygene_gss_dir  )
    if include_rev:
        os.mkdir(bygene_rev_dir  )

if include_byboth:
    os.mkdir(byboth_dir     )
    os.mkdir(byboth_data_dir)
    if not fan_etal_2011:
        if include_lorad or include_ghm:
            os.mkdir(byboth_lorad_dir )
        if include_gss:
            os.mkdir(byboth_gss_dir  )
    if include_rev:
        os.mkdir(byboth_rev_dir  )
       
####################
# Write data files #
####################

if include_unpart:
    unpart_boundaries  = writeNexusFile(os.path.join(unpart_data_dir, 'unpart.nex'),  taxa, None, [nchar0], [unpartseqs])
    print('unpart:')
    print('  nchar0 = %d' % nchar0)

if include_bycodon:
    bycodon_boundaries  = writeNexusFile(os.path.join(bycodon_data_dir,'bycodon.nex'), taxa, None, [nchar1,nchar2,nchar3], [codon1seqs,codon2seqs,codon3seqs])
    codon1st_boundaries = writeNexusFile(os.path.join(bycodon_data_dir,'codon1st.nex'), taxa, None, [nchar1], [codon1seqs])
    codon2nd_boundaries = writeNexusFile(os.path.join(bycodon_data_dir,'codon2nd.nex'), taxa, None, [nchar2], [codon2seqs])
    codon3rd_boundaries = writeNexusFile(os.path.join(bycodon_data_dir,'codon3rd.nex'), taxa, None, [nchar3], [codon3seqs])
    print('bycodon:')
    print('  nchar1 = %d' % nchar1)
    print('  nchar2 = %d' % nchar2)
    print('  nchar3 = %d' % nchar3)

if include_bygene:
    bygene_boundaries  = writeNexusFile(os.path.join(bygene_data_dir, 'bygene.nex'),  taxa, None, [ncharCOI,ncharCOII,ncharATPase6,ncharATPase8], [COIseqs,COIIseqs,ATPase6seqs,ATPase8seqs])
    COI_boundaries     = writeNexusFile(os.path.join(bygene_data_dir, 'COI.nex'),     taxa, None, [ncharCOI],     [COIseqs])
    COII_boundaries    = writeNexusFile(os.path.join(bygene_data_dir, 'COII.nex'),    taxa, None, [ncharCOII],    [COIIseqs])
    ATPase6_boundaries = writeNexusFile(os.path.join(bygene_data_dir, 'ATPase6.nex'), taxa, None, [ncharATPase6], [ATPase6seqs])
    ATPase8_boundaries = writeNexusFile(os.path.join(bygene_data_dir, 'ATPase8.nex'), taxa, None, [ncharATPase8], [ATPase8seqs])
    print('bygene:')
    print('  ncharCOI     = %d' % ncharCOI)
    print('  ncharCOII    = %d' % ncharCOII)
    print('  ncharATPase6 = %d' % ncharATPase6)
    print('  ncharATPase8 = %d' % ncharATPase8)

if include_byboth:
    byboth_boundaries   = writeNexusFile(os.path.join(byboth_data_dir, 'byboth.nex'),   taxa, None, [ncharCOI1, ncharCOI2, ncharCOI3, ncharCOII1, ncharCOII2, ncharCOII3, ncharATPase61, ncharATPase62, ncharATPase63, ncharATPase81, ncharATPase82, ncharATPase83], [COIseqs1, COIseqs2, COIseqs3, COIIseqs1, COIIseqs2, COIIseqs3, ATPase6seqs1, ATPase6seqs2, ATPase6seqs3, ATPase8seqs1, ATPase8seqs2, ATPase8seqs3])
    COI1_boundaries     = writeNexusFile(os.path.join(byboth_data_dir, 'COI-1st.nex'),     taxa, None, [ncharCOI1],     [COIseqs1])
    COI2_boundaries     = writeNexusFile(os.path.join(byboth_data_dir, 'COI-2nd.nex'),     taxa, None, [ncharCOI2],     [COIseqs2])
    COI3_boundaries     = writeNexusFile(os.path.join(byboth_data_dir, 'COI-3rd.nex'),     taxa, None, [ncharCOI3],     [COIseqs3])
    COII1_boundaries    = writeNexusFile(os.path.join(byboth_data_dir, 'COII-1st.nex'),    taxa, None, [ncharCOII1],    [COIIseqs1])
    COII2_boundaries    = writeNexusFile(os.path.join(byboth_data_dir, 'COII-2nd.nex'),    taxa, None, [ncharCOII2],    [COIIseqs2])
    COII3_boundaries    = writeNexusFile(os.path.join(byboth_data_dir, 'COII-3rd.nex'),    taxa, None, [ncharCOII3],    [COIIseqs3])
    ATPase61_boundaries = writeNexusFile(os.path.join(byboth_data_dir, 'ATPase6-1st.nex'), taxa, None, [ncharATPase61], [ATPase6seqs1])
    ATPase62_boundaries = writeNexusFile(os.path.join(byboth_data_dir, 'ATPase6-2nd.nex'), taxa, None, [ncharATPase62], [ATPase6seqs2])
    ATPase63_boundaries = writeNexusFile(os.path.join(byboth_data_dir, 'ATPase6-3rd.nex'), taxa, None, [ncharATPase63], [ATPase6seqs3])
    ATPase81_boundaries = writeNexusFile(os.path.join(byboth_data_dir, 'ATPase8-1st.nex'), taxa, None, [ncharATPase81], [ATPase8seqs1])
    ATPase82_boundaries = writeNexusFile(os.path.join(byboth_data_dir, 'ATPase8-2nd.nex'), taxa, None, [ncharATPase82], [ATPase8seqs2])
    ATPase83_boundaries = writeNexusFile(os.path.join(byboth_data_dir, 'ATPase8-3rd.nex'), taxa, None, [ncharATPase83], [ATPase8seqs3])
    print('byboth:')
    print('  ncharCOI1     = %d' % ncharCOI1)
    print('  ncharCOI2     = %d' % ncharCOI2)
    print('  ncharCOI3     = %d' % ncharCOI3)
    print('  ncharCOII1    = %d' % ncharCOII1)
    print('  ncharCOII2    = %d' % ncharCOII2)
    print('  ncharCOII3    = %d' % ncharCOII3)
    print('  ncharATPase61 = %d' % ncharATPase61)
    print('  ncharATPase62 = %d' % ncharATPase62)
    print('  ncharATPase63 = %d' % ncharATPase63)
    print('  ncharATPase81 = %d' % ncharATPase81)
    print('  ncharATPase82 = %d' % ncharATPase82)
    print('  ncharATPase83 = %d' % ncharATPase83)

########################
# Create slurm scripts #
########################

if include_lorad or include_ghm:
    submit_lorad = '#!/bin/bash\n\n'
    slurm_lorad_script_template = open('slurm-lorad-template.txt','r').read()

if include_gss:
    submit_gss = '#!/bin/bash\n\n'
    slurm_gss_script_template = open('slurm-gss-template.txt','r').read()

if include_rev:
    submit_rev = '#!/bin/bash\n\n'
    slurm_rev_script_template = open('slurm-rev-template.txt','r').read()

#if include_lorad or include_ghm or include_gss or include_rev:
#    submit_all = '#!/bin/bash\n\n'
    
# These two variables determine whether the user is emailed when each run ends
# Set notify = True to do this or notify = False to avoid lots of email notifications
sbatch_mail_type = ''
sbatch_mail_user = ''
if notify:
    sbatch_mail_type = '#SBATCH --mail-type=END'
    sbatch_mail_user = '#SBATCH --mail-user=%s' % email[userid]

#########################################
# Create slurm file for unpart analyses #
#########################################

if include_unpart:
    if not fan_etal_2011 and (include_lorad or include_ghm):
        unpart_lorad_slurm_filename = os.path.join(unpart_lorad_dir,'s.sh')
        unpart_lorad_slurm_contents = re.sub('__JOBNAME__',             'lrd0%d' % dest_dir_index, slurm_lorad_script_template, re.M | re.S)
        unpart_lorad_slurm_contents = re.sub('__SBATCH_MAIL_TYPE__',    sbatch_mail_type,          unpart_lorad_slurm_contents, re.M | re.S)
        unpart_lorad_slurm_contents = re.sub('__SBATCH_MAIL_USER__',    sbatch_mail_user,          unpart_lorad_slurm_contents, re.M | re.S)
        unpart_lorad_slurm_contents = re.sub('__USERID__',              userid,                    unpart_lorad_slurm_contents, re.M | re.S)
        unpart_lorad_slurm_contents = re.sub('__FNPREFIX__',            "unpart-lorad-",           unpart_lorad_slurm_contents, re.M | re.S)
        submit_lorad += 'cd %s; sbatch s.sh; cd ../..\n' % unpart_lorad_dir_short
        #submit_all   += 'cd %s; sbatch s.sh; cd ../..\n' % unpart_lorad_dir_short
        f = open(unpart_lorad_slurm_filename,'w')
        f.write(unpart_lorad_slurm_contents)
        f.close()

    if not fan_etal_2011 and include_gss:
        unpart_gss_slurm_filename = os.path.join(unpart_gss_dir,'s.sh')
        unpart_gss_slurm_contents = re.sub('__JOBNAME__',            'gss0%d' % dest_dir_index, slurm_gss_script_template, re.M | re.S)
        unpart_gss_slurm_contents = re.sub('__SBATCH_MAIL_TYPE__',   sbatch_mail_type,          unpart_gss_slurm_contents, re.M | re.S)
        unpart_gss_slurm_contents = re.sub('__SBATCH_MAIL_USER__',   sbatch_mail_user,          unpart_gss_slurm_contents, re.M | re.S)
        unpart_gss_slurm_contents = re.sub('__USERID__',             userid,                    unpart_gss_slurm_contents, re.M | re.S)
        unpart_gss_slurm_contents = re.sub('__FNPREFIX__',           "unpart-gss-",             unpart_gss_slurm_contents, re.M | re.S)
        submit_gss += 'cd %s; sbatch s.sh; cd ../..\n' % unpart_gss_dir_short
        #submit_all += 'cd %s; sbatch s.sh; cd ../..\n' % unpart_gss_dir_short
        f = open(unpart_gss_slurm_filename,'w')
        f.write(unpart_gss_slurm_contents)
        f.close()

    if include_rev:
        unpart_rev_slurm_filename = os.path.join(unpart_rev_dir,'s.sh')
        unpart_rev_slurm_contents = re.sub('__JOBNAME__',             'rev0%d' % dest_dir_index, slurm_rev_script_template, re.M | re.S)
        unpart_rev_slurm_contents = re.sub('__SBATCH_MAIL_TYPE__',    sbatch_mail_type,          unpart_rev_slurm_contents, re.M | re.S)
        unpart_rev_slurm_contents = re.sub('__SBATCH_MAIL_USER__',    sbatch_mail_user,          unpart_rev_slurm_contents, re.M | re.S)
        unpart_rev_slurm_contents = re.sub('__USERID__',              userid,                    unpart_rev_slurm_contents, re.M | re.S)
        unpart_rev_slurm_contents = re.sub('__FNPREFIX__',            "unpart-rev-",             unpart_rev_slurm_contents, re.M | re.S)
        submit_rev += 'cd %s; sbatch s.sh; cd ../..\n' % unpart_rev_dir_short
        #submit_all += 'cd %s; sbatch s.sh; cd ../..\n' % unpart_rev_dir_short
        f = open(unpart_rev_slurm_filename,'w')
        f.write(unpart_rev_slurm_contents)
        f.close()

##########################################
# Create slurm file for bycodon analyses #
##########################################

if include_bycodon:
    if not fan_etal_2011 and (include_lorad or include_ghm):
        bycodon_lorad_slurm_filename = os.path.join(bycodon_lorad_dir,'s.sh')
        bycodon_lorad_slurm_contents = re.sub('__JOBNAME__',             'lrd1%d' % dest_dir_index, slurm_lorad_script_template, re.M | re.S)
        bycodon_lorad_slurm_contents = re.sub('__SBATCH_MAIL_TYPE__',    sbatch_mail_type,          bycodon_lorad_slurm_contents, re.M | re.S)
        bycodon_lorad_slurm_contents = re.sub('__SBATCH_MAIL_USER__',    sbatch_mail_user,          bycodon_lorad_slurm_contents, re.M | re.S)
        bycodon_lorad_slurm_contents = re.sub('__USERID__',              userid,                    bycodon_lorad_slurm_contents, re.M | re.S)
        bycodon_lorad_slurm_contents = re.sub('__FNPREFIX__',            "bycodon-lorad-",          bycodon_lorad_slurm_contents, re.M | re.S)
        submit_lorad += 'cd %s; sbatch s.sh; cd ../..\n' % bycodon_lorad_dir_short
        #submit_all   += 'cd %s; sbatch s.sh; cd ../..\n' % bycodon_lorad_dir_short
        f = open(bycodon_lorad_slurm_filename,'w')
        f.write(bycodon_lorad_slurm_contents)
        f.close()

    if not fan_etal_2011 and include_gss:
        bycodon_gss_slurm_filename = os.path.join(bycodon_gss_dir,'s.sh')
        bycodon_gss_slurm_contents = re.sub('__JOBNAME__',             'gss1%d' % dest_dir_index, slurm_gss_script_template, re.M | re.S)
        bycodon_gss_slurm_contents = re.sub('__SBATCH_MAIL_TYPE__',    sbatch_mail_type,          bycodon_gss_slurm_contents, re.M | re.S)
        bycodon_gss_slurm_contents = re.sub('__SBATCH_MAIL_USER__',    sbatch_mail_user,          bycodon_gss_slurm_contents, re.M | re.S)
        bycodon_gss_slurm_contents = re.sub('__USERID__',              userid,                    bycodon_gss_slurm_contents, re.M | re.S)
        bycodon_gss_slurm_contents = re.sub('__FNPREFIX__',            "bycodon-gss-",            bycodon_gss_slurm_contents, re.M | re.S)
        submit_gss += 'cd %s; sbatch s.sh; cd ../..\n' % bycodon_gss_dir_short
        #submit_all += 'cd %s; sbatch s.sh; cd ../..\n' % bycodon_gss_dir_short
        f = open(bycodon_gss_slurm_filename,'w')
        f.write(bycodon_gss_slurm_contents)
        f.close()

    if include_rev:
        bycodon_rev_slurm_filename = os.path.join(bycodon_rev_dir,'s.sh')
        bycodon_rev_slurm_contents = re.sub('__JOBNAME__',             'rev1%d' % dest_dir_index, slurm_rev_script_template, re.M | re.S)
        bycodon_rev_slurm_contents = re.sub('__SBATCH_MAIL_TYPE__',    sbatch_mail_type,          bycodon_rev_slurm_contents, re.M | re.S)
        bycodon_rev_slurm_contents = re.sub('__SBATCH_MAIL_USER__',    sbatch_mail_user,          bycodon_rev_slurm_contents, re.M | re.S)
        bycodon_rev_slurm_contents = re.sub('__USERID__',              userid,                    bycodon_rev_slurm_contents, re.M | re.S)
        bycodon_rev_slurm_contents = re.sub('__FNPREFIX__',            "bycodon-rev-",            bycodon_rev_slurm_contents, re.M | re.S)
        submit_rev += 'cd %s; sbatch s.sh; cd ../..\n' % bycodon_rev_dir_short
        #submit_all += 'cd %s; sbatch s.sh; cd ../..\n' % bycodon_rev_dir_short
        f = open(bycodon_rev_slurm_filename,'w')
        f.write(bycodon_rev_slurm_contents)
        f.close()
    
#########################################
# Create slurm file for bygene analyses #
#########################################

if include_bygene:
    if not fan_etal_2011 and (include_lorad or include_ghm):
        bygene_lorad_slurm_filename = os.path.join(bygene_lorad_dir,'s.sh')
        bygene_lorad_slurm_contents = re.sub('__JOBNAME__',             'lrd2%d' % dest_dir_index, slurm_lorad_script_template, re.M | re.S)
        bygene_lorad_slurm_contents = re.sub('__SBATCH_MAIL_TYPE__',    sbatch_mail_type,          bygene_lorad_slurm_contents, re.M | re.S)
        bygene_lorad_slurm_contents = re.sub('__SBATCH_MAIL_USER__',    sbatch_mail_user,          bygene_lorad_slurm_contents, re.M | re.S)
        bygene_lorad_slurm_contents = re.sub('__USERID__',              userid,                    bygene_lorad_slurm_contents, re.M | re.S)
        bygene_lorad_slurm_contents = re.sub('__FNPREFIX__',            "bygene-lorad-",           bygene_lorad_slurm_contents, re.M | re.S)
        submit_lorad += 'cd %s; sbatch s.sh; cd ../..\n' % bygene_lorad_dir_short
        #submit_all   += 'cd %s; sbatch s.sh; cd ../..\n' % bygene_lorad_dir_short
        f = open(bygene_lorad_slurm_filename,'w')
        f.write(bygene_lorad_slurm_contents)
        f.close()

    if not fan_etal_2011 and include_gss:
        bygene_gss_slurm_filename = os.path.join(bygene_gss_dir,'s.sh')
        bygene_gss_slurm_contents = re.sub('__JOBNAME__',             'gss2%d' % dest_dir_index, slurm_gss_script_template, re.M | re.S)
        bygene_gss_slurm_contents = re.sub('__SBATCH_MAIL_TYPE__',    sbatch_mail_type,          bygene_gss_slurm_contents, re.M | re.S)
        bygene_gss_slurm_contents = re.sub('__SBATCH_MAIL_USER__',    sbatch_mail_user,          bygene_gss_slurm_contents, re.M | re.S)
        bygene_gss_slurm_contents = re.sub('__USERID__',              userid,                    bygene_gss_slurm_contents, re.M | re.S)
        bygene_gss_slurm_contents = re.sub('__FNPREFIX__',            "bygene-gss-",             bygene_gss_slurm_contents, re.M | re.S)
        submit_gss += 'cd %s; sbatch s.sh; cd ../..\n' % bygene_gss_dir_short
        #submit_all += 'cd %s; sbatch s.sh; cd ../..\n' % bygene_gss_dir_short
        f = open(bygene_gss_slurm_filename,'w')
        f.write(bygene_gss_slurm_contents)
        f.close()

    if include_rev:
        bygene_rev_slurm_filename = os.path.join(bygene_rev_dir,'s.sh')
        bygene_rev_slurm_contents = re.sub('__JOBNAME__',            'rev2%d' % dest_dir_index, slurm_rev_script_template, re.M | re.S)
        bygene_rev_slurm_contents = re.sub('__SBATCH_MAIL_TYPE__',    sbatch_mail_type,         bygene_rev_slurm_contents, re.M | re.S)
        bygene_rev_slurm_contents = re.sub('__SBATCH_MAIL_USER__',    sbatch_mail_user,         bygene_rev_slurm_contents, re.M | re.S)
        bygene_rev_slurm_contents = re.sub('__USERID__',              userid,                   bygene_rev_slurm_contents, re.M | re.S)
        bygene_rev_slurm_contents = re.sub('__FNPREFIX__',            "bygene-rev-",            bygene_rev_slurm_contents, re.M | re.S)
        submit_rev += 'cd %s; sbatch s.sh; cd ../..\n' % bygene_rev_dir_short
        #submit_all += 'cd %s; sbatch s.sh; cd ../..\n' % bygene_rev_dir_short
        f = open(bygene_rev_slurm_filename,'w')
        f.write(bygene_rev_slurm_contents)
        f.close()

#########################################
# Create slurm file for byboth analyses #
#########################################

if include_byboth:
    if not fan_etal_2011 and (include_lorad or include_ghm):
        byboth_lorad_slurm_filename = os.path.join(byboth_lorad_dir,'s.sh')
        byboth_lorad_slurm_contents = re.sub('__JOBNAME__',             'lrd3%d' % dest_dir_index, slurm_lorad_script_template, re.M | re.S)
        byboth_lorad_slurm_contents = re.sub('__SBATCH_MAIL_TYPE__',    sbatch_mail_type,          byboth_lorad_slurm_contents, re.M | re.S)
        byboth_lorad_slurm_contents = re.sub('__SBATCH_MAIL_USER__',    sbatch_mail_user,          byboth_lorad_slurm_contents, re.M | re.S)
        byboth_lorad_slurm_contents = re.sub('__USERID__',              userid,                    byboth_lorad_slurm_contents, re.M | re.S)
        byboth_lorad_slurm_contents = re.sub('__FNPREFIX__',            "byboth-lorad-",           byboth_lorad_slurm_contents, re.M | re.S)
        submit_lorad += 'cd %s; sbatch s.sh; cd ../..\n' % byboth_lorad_dir_short
        #submit_all   += 'cd %s; sbatch s.sh; cd ../..\n' % byboth_lorad_dir_short
        f = open(byboth_lorad_slurm_filename,'w')
        f.write(byboth_lorad_slurm_contents)
        f.close()

    if not fan_etal_2011 and include_gss:
        byboth_gss_slurm_filename = os.path.join(byboth_gss_dir,'s.sh')
        byboth_gss_slurm_contents = re.sub('__JOBNAME__',             'gss3%d' % dest_dir_index, slurm_gss_script_template, re.M | re.S)
        byboth_gss_slurm_contents = re.sub('__SBATCH_MAIL_TYPE__',    sbatch_mail_type,          byboth_gss_slurm_contents, re.M | re.S)
        byboth_gss_slurm_contents = re.sub('__SBATCH_MAIL_USER__',    sbatch_mail_user,          byboth_gss_slurm_contents, re.M | re.S)
        byboth_gss_slurm_contents = re.sub('__USERID__',              userid,                    byboth_gss_slurm_contents, re.M | re.S)
        byboth_gss_slurm_contents = re.sub('__FNPREFIX__',            "byboth-gss-",             byboth_gss_slurm_contents, re.M | re.S)
        submit_gss += 'cd %s; sbatch s.sh; cd ../..\n' % byboth_gss_dir_short
        #submit_all += 'cd %s; sbatch s.sh; cd ../..\n' % byboth_gss_dir_short
        f = open(byboth_gss_slurm_filename,'w')
        f.write(byboth_gss_slurm_contents)
        f.close()

    if include_rev:
        byboth_rev_slurm_filename = os.path.join(byboth_rev_dir,'s.sh')
        byboth_rev_slurm_contents = re.sub('__JOBNAME__',             'rev3%d' % dest_dir_index, slurm_rev_script_template, re.M | re.S)
        byboth_rev_slurm_contents = re.sub('__SBATCH_MAIL_TYPE__',    sbatch_mail_type,          byboth_rev_slurm_contents, re.M | re.S)
        byboth_rev_slurm_contents = re.sub('__SBATCH_MAIL_USER__',    sbatch_mail_user,          byboth_rev_slurm_contents, re.M | re.S)
        byboth_rev_slurm_contents = re.sub('__USERID__',              userid,                    byboth_rev_slurm_contents, re.M | re.S)
        byboth_rev_slurm_contents = re.sub('__FNPREFIX__',            "byboth-rev-",             byboth_rev_slurm_contents, re.M | re.S)
        submit_rev += 'cd %s; sbatch s.sh; cd ../..\n' % byboth_rev_dir_short
        #submit_all += 'cd %s; sbatch s.sh; cd ../..\n' % byboth_rev_dir_short
        f = open(byboth_rev_slurm_filename,'w')
        f.write(byboth_rev_slurm_contents)
        f.close()
    
#########################
# Create submit scripts #
#########################

if include_lorad or include_ghm:
    submit_lorad_filename = os.path.join(dest_dir,'submit-lorad.sh')
    f = open(submit_lorad_filename,'w')
    f.write(submit_lorad)
    f.close()

if include_gss:
    submit_gss_filename = os.path.join(dest_dir,'submit-gss.sh')
    f = open(submit_gss_filename,'w')
    f.write(submit_gss)
    f.close()

if include_rev:
    submit_rev_filename = os.path.join(dest_dir,'submit-rev.sh')
    f = open(submit_rev_filename,'w')
    f.write(submit_rev)
    f.close()

# Commented out because it makes no sense to start gss before lorad is finished    
#if include_lorad or include_ghm or include_gss or include_rev:
#    submit_all_filename = os.path.join(dest_dir,'submit-all.sh')
#    f = open(submit_all_filename,'w')
#    f.write(submit_all)
#    f.close()

#############################
# Create lorad.conf scripts #
#############################

#################################################
# Create lorad.conf scripts for unpart analyses #
#################################################

if include_unpart:
    if not fan_etal_2011 and (include_lorad or include_ghm):
        unpart_lorad_conf_template = open('conf-unpart-lorad-template.txt','r').read()
        unpart_lorad_conf_filename = os.path.join(unpart_lorad_dir,'lorad.conf')
        unpart_lorad_conf_contents = re.sub('__LAST_SITE__',         str(unpart_boundaries[0]),  unpart_lorad_conf_template, re.M | re.S)
        unpart_lorad_conf_contents = re.sub('__BURNIN__',            lorad_burnin,               unpart_lorad_conf_contents, re.M | re.S)
        unpart_lorad_conf_contents = re.sub('__NITER__',             lorad_niter,                unpart_lorad_conf_contents, re.M | re.S)
        unpart_lorad_conf_contents = re.sub('__SAMPLEFREQ__',        lorad_samplefreq,           unpart_lorad_conf_contents, re.M | re.S)
        unpart_lorad_conf_contents = re.sub('__PRINTFREQ__',         lorad_printfreq,            unpart_lorad_conf_contents, re.M | re.S)
        unpart_lorad_conf_contents = re.sub('__COVERAGE10__',        lorad_coverage10,           unpart_lorad_conf_contents, re.M | re.S)
        unpart_lorad_conf_contents = re.sub('__COVERAGE20__',        lorad_coverage20,           unpart_lorad_conf_contents, re.M | re.S)
        unpart_lorad_conf_contents = re.sub('__COVERAGE30__',        lorad_coverage30,           unpart_lorad_conf_contents, re.M | re.S)
        unpart_lorad_conf_contents = re.sub('__COVERAGE40__',        lorad_coverage40,           unpart_lorad_conf_contents, re.M | re.S)
        unpart_lorad_conf_contents = re.sub('__COVERAGE50__',        lorad_coverage50,           unpart_lorad_conf_contents, re.M | re.S)
        unpart_lorad_conf_contents = re.sub('__COVERAGE60__',        lorad_coverage60,           unpart_lorad_conf_contents, re.M | re.S)
        unpart_lorad_conf_contents = re.sub('__COVERAGE70__',        lorad_coverage70,           unpart_lorad_conf_contents, re.M | re.S)
        unpart_lorad_conf_contents = re.sub('__COVERAGE80__',        lorad_coverage80,           unpart_lorad_conf_contents, re.M | re.S)
        unpart_lorad_conf_contents = re.sub('__COVERAGE90__',        lorad_coverage90,           unpart_lorad_conf_contents, re.M | re.S)
        unpart_lorad_conf_contents = re.sub('__COVERAGE95__',        lorad_coverage95,           unpart_lorad_conf_contents, re.M | re.S)
        unpart_lorad_conf_contents = re.sub('__COVERAGE99__',        lorad_coverage99,           unpart_lorad_conf_contents, re.M | re.S)
        unpart_lorad_conf_contents = re.sub('__RNSEED__',            rnseed,                     unpart_lorad_conf_contents, re.M | re.S)
        unpart_lorad_conf_contents = re.sub('__TREEFILE__',          tree_file_name,             unpart_lorad_conf_contents, re.M | re.S)
        f = open(unpart_lorad_conf_filename,'w')
        f.write(unpart_lorad_conf_contents)
        f.close()
    
    if not fan_etal_2011 and include_gss:
        # This part commented out because, in the paper, GSS used reference distributions computed from the posterior sample created for LoRaD.
        ## GSS requires a reference distribution estimated from an initial MCMC run, so this part is identical to the include_lorad
        ## part above except that gss_burnin, gss_niter, gss_samplefreq, and gss_printfreq are used instead of lorad_burnin, lorad_niter,
        ## lorad_samplefreq, and lorad_printfreq (also, the file generated is <unpart_gss_dir>/lorad-mcmc.conf, not <unpart_lorad_dir>/lorad.conf).
        ##unpart_lorad_conf_template = open('conf-unpart-lorad-template.txt','r').read()
        #unpart_lorad_conf_filename = os.path.join(unpart_gss_dir,'lorad-mcmc.conf')
        #unpart_lorad_conf_contents = re.sub('__LAST_SITE__',         str(unpart_boundaries[0]),  unpart_lorad_conf_template, re.M | re.S)
        #unpart_lorad_conf_contents = re.sub('__BURNIN__',            gss_burnin,                 unpart_lorad_conf_contents, re.M | re.S)
        #unpart_lorad_conf_contents = re.sub('__NITER__',             gss_niter,                  unpart_lorad_conf_contents, re.M | re.S)
        #unpart_lorad_conf_contents = re.sub('__SAMPLEFREQ__',        gss_samplefreq,             unpart_lorad_conf_contents, re.M | re.S)
        #unpart_lorad_conf_contents = re.sub('__PRINTFREQ__',         gss_printfreq,              unpart_lorad_conf_contents, re.M | re.S)
        #unpart_lorad_conf_contents = re.sub('__COVERAGE10__',        lorad_coverage10,           unpart_lorad_conf_contents, re.M | re.S)
        #unpart_lorad_conf_contents = re.sub('__COVERAGE20__',        lorad_coverage20,           unpart_lorad_conf_contents, re.M | re.S)
        #unpart_lorad_conf_contents = re.sub('__COVERAGE30__',        lorad_coverage30,           unpart_lorad_conf_contents, re.M | re.S)
        #unpart_lorad_conf_contents = re.sub('__COVERAGE40__',        lorad_coverage40,           unpart_lorad_conf_contents, re.M | re.S)
        #unpart_lorad_conf_contents = re.sub('__COVERAGE50__',        lorad_coverage50,           unpart_lorad_conf_contents, re.M | re.S)
        #unpart_lorad_conf_contents = re.sub('__COVERAGE60__',        lorad_coverage60,           unpart_lorad_conf_contents, re.M | re.S)
        #unpart_lorad_conf_contents = re.sub('__COVERAGE70__',        lorad_coverage70,           unpart_lorad_conf_contents, re.M | re.S)
        #unpart_lorad_conf_contents = re.sub('__COVERAGE80__',        lorad_coverage80,           unpart_lorad_conf_contents, re.M | re.S)
        #unpart_lorad_conf_contents = re.sub('__COVERAGE90__',        lorad_coverage90,           unpart_lorad_conf_contents, re.M | re.S)
        #unpart_lorad_conf_contents = re.sub('__COVERAGE95__',        lorad_coverage95,           unpart_lorad_conf_contents, re.M | re.S)
        #unpart_lorad_conf_contents = re.sub('__COVERAGE99__',        lorad_coverage99,           unpart_lorad_conf_contents, re.M | re.S)
        #unpart_lorad_conf_contents = re.sub('__USE_REGRESSION__',    lorad_regression,           unpart_lorad_conf_contents, re.M | re.S)
        #unpart_lorad_conf_contents = re.sub('__LINEAR_REGRESSION__', lorad_linearregression,     unpart_lorad_conf_contents, re.M | re.S)
        #unpart_lorad_conf_contents = re.sub('__RNSEED__',            rnseed,                     unpart_lorad_conf_contents, re.M | re.S)
        #unpart_lorad_conf_contents = re.sub('__TREEFILE__',          tree_file_name,             unpart_lorad_conf_contents, re.M | re.S)
        #f = open(unpart_lorad_conf_filename,'w')
        #f.write(unpart_lorad_conf_contents)
        #f.close()

        unpart_gss_conf_template = open('conf-unpart-gss-template.txt','r').read()
        unpart_gss_conf_filename = os.path.join(unpart_gss_dir,'lorad-gss.conf')
        unpart_gss_conf_contents = re.sub('__LAST_SITE__',  str(unpart_boundaries[0]), unpart_gss_conf_template, re.M | re.S)
        unpart_gss_conf_contents = re.sub('__BURNIN__',     gss_burnin,                unpart_gss_conf_contents, re.M | re.S)
        unpart_gss_conf_contents = re.sub('__NITER__',      gss_niter,                 unpart_gss_conf_contents, re.M | re.S)
        unpart_gss_conf_contents = re.sub('__SAMPLEFREQ__', gss_samplefreq,            unpart_gss_conf_contents, re.M | re.S)
        unpart_gss_conf_contents = re.sub('__PRINTFREQ__',  gss_printfreq,             unpart_gss_conf_contents, re.M | re.S)
        unpart_gss_conf_contents = re.sub('__NSTONES__',    gss_nstones,               unpart_gss_conf_contents, re.M | re.S)
        unpart_gss_conf_contents = re.sub('__ALPHA__',      gss_alpha,                 unpart_gss_conf_contents, re.M | re.S)
        unpart_gss_conf_contents = re.sub('__RNSEED__',     rnseed,                    unpart_gss_conf_contents, re.M | re.S)
        unpart_gss_conf_contents = re.sub('__TREEFILE__',   tree_file_name,            unpart_gss_conf_contents, re.M | re.S)
        f = open(unpart_gss_conf_filename,'w')
        f.write(unpart_gss_conf_contents)
        f.close()

    if include_rev:
        unpart_rev_template = open('rev-unpart-template.txt','r').read()
        unpart_rev_filename = os.path.join(unpart_rev_dir,'ss.Rev')
        unpart_rev_contents = re.sub('__RNSEED__',         rnseed,              unpart_rev_template, re.M | re.S)
        unpart_rev_contents = re.sub('__BURNIN__',         rev_burnin,          unpart_rev_contents, re.M | re.S)
        unpart_rev_contents = re.sub('__NITER__',          rev_niter,           unpart_rev_contents, re.M | re.S)
        unpart_rev_contents = re.sub('__SAMPLEFREQ__',     rev_samplefreq,      unpart_rev_contents, re.M | re.S)
        unpart_rev_contents = re.sub('__PRINTFREQ__',      rev_printfreq,       unpart_rev_contents, re.M | re.S)
        unpart_rev_contents = re.sub('__TUNINGINTERVAL__', rev_tuning_interval, unpart_rev_contents, re.M | re.S)
        unpart_rev_contents = re.sub('__NSTONES__',        rev_nstones,         unpart_rev_contents, re.M | re.S)
        unpart_rev_contents = re.sub('__ALPHA__',          rev_alpha,           unpart_rev_contents, re.M | re.S)
        unpart_rev_contents = re.sub('__TREEFILE__',       tree_file_name,      unpart_rev_contents, re.M | re.S)
        if fan_etal_2011:
            unpart_rev_contents = re.sub('__FAN_ETAL_2011__',  'TRUE',          unpart_rev_contents, re.M | re.S)
        else:
            unpart_rev_contents = re.sub('__FAN_ETAL_2011__',  'FALSE',         unpart_rev_contents, re.M | re.S)
        f = open(unpart_rev_filename,'w')
        f.write(unpart_rev_contents)
        f.close()

##################################################
# Create lorad.conf scripts for bycodon analyses #
##################################################

if include_bycodon:
    if not fan_etal_2011 and (include_lorad or include_ghm):
        bycodon_lorad_conf_template = open('conf-bycodon-lorad-template.txt','r').read()
        bycodon_lorad_conf_filename = os.path.join(bycodon_lorad_dir,'lorad.conf')
        bycodon_lorad_conf_contents = re.sub('__FIRST_SITE_1ST_CODON__', '1',                            bycodon_lorad_conf_template, re.M | re.S)
        bycodon_lorad_conf_contents = re.sub('__LAST_SITE_1ST_CODON__',  str(bycodon_boundaries[0]),     bycodon_lorad_conf_contents, re.M | re.S)
        bycodon_lorad_conf_contents = re.sub('__FIRST_SITE_2ND_CODON__', str(bycodon_boundaries[0] + 1), bycodon_lorad_conf_contents, re.M | re.S)
        bycodon_lorad_conf_contents = re.sub('__LAST_SITE_2ND_CODON__',  str(bycodon_boundaries[1]),     bycodon_lorad_conf_contents, re.M | re.S)
        bycodon_lorad_conf_contents = re.sub('__FIRST_SITE_3RD_CODON__', str(bycodon_boundaries[1] + 1), bycodon_lorad_conf_contents, re.M | re.S)
        bycodon_lorad_conf_contents = re.sub('__LAST_SITE_3RD_CODON__',  str(bycodon_boundaries[2]),     bycodon_lorad_conf_contents, re.M | re.S)
        bycodon_lorad_conf_contents = re.sub('__BURNIN__',               lorad_burnin,                   bycodon_lorad_conf_contents, re.M | re.S)
        bycodon_lorad_conf_contents = re.sub('__NITER__',                lorad_niter,                    bycodon_lorad_conf_contents, re.M | re.S)
        bycodon_lorad_conf_contents = re.sub('__SAMPLEFREQ__',           lorad_samplefreq,               bycodon_lorad_conf_contents, re.M | re.S)
        bycodon_lorad_conf_contents = re.sub('__PRINTFREQ__',            lorad_printfreq,                bycodon_lorad_conf_contents, re.M | re.S)
        bycodon_lorad_conf_contents = re.sub('__COVERAGE10__',           lorad_coverage10,               bycodon_lorad_conf_contents, re.M | re.S)
        bycodon_lorad_conf_contents = re.sub('__COVERAGE20__',           lorad_coverage20,               bycodon_lorad_conf_contents, re.M | re.S)
        bycodon_lorad_conf_contents = re.sub('__COVERAGE30__',           lorad_coverage30,               bycodon_lorad_conf_contents, re.M | re.S)
        bycodon_lorad_conf_contents = re.sub('__COVERAGE40__',           lorad_coverage40,               bycodon_lorad_conf_contents, re.M | re.S)
        bycodon_lorad_conf_contents = re.sub('__COVERAGE50__',           lorad_coverage50,               bycodon_lorad_conf_contents, re.M | re.S)
        bycodon_lorad_conf_contents = re.sub('__COVERAGE60__',           lorad_coverage60,               bycodon_lorad_conf_contents, re.M | re.S)
        bycodon_lorad_conf_contents = re.sub('__COVERAGE70__',           lorad_coverage70,               bycodon_lorad_conf_contents, re.M | re.S)
        bycodon_lorad_conf_contents = re.sub('__COVERAGE80__',           lorad_coverage80,               bycodon_lorad_conf_contents, re.M | re.S)
        bycodon_lorad_conf_contents = re.sub('__COVERAGE90__',           lorad_coverage90,               bycodon_lorad_conf_contents, re.M | re.S)
        bycodon_lorad_conf_contents = re.sub('__COVERAGE95__',           lorad_coverage95,               bycodon_lorad_conf_contents, re.M | re.S)
        bycodon_lorad_conf_contents = re.sub('__COVERAGE99__',           lorad_coverage99,               bycodon_lorad_conf_contents, re.M | re.S)
        bycodon_lorad_conf_contents = re.sub('__RNSEED__',               rnseed,                         bycodon_lorad_conf_contents, re.M | re.S)
        bycodon_lorad_conf_contents = re.sub('__TREEFILE__',             tree_file_name,                 bycodon_lorad_conf_contents, re.M | re.S)
        f = open(bycodon_lorad_conf_filename,'w')
        f.write(bycodon_lorad_conf_contents)
        f.close()

    if not fan_etal_2011 and include_gss:
        # This part commented out because, in the paper, GSS used reference distributions computed from the posterior sample created for LoRaD.
        ## GSS requires a reference distribution estimated from an initial MCMC run, so this part is identical to the include_lorad
        ## part above except that gss_burnin, gss_niter, gss_samplefreq, and gss_printfreq are used instead of lorad_burnin, lorad_niter,
        ## lorad_samplefreq, and lorad_printfreq (also, the file generated is <bycodon_gss_dir>/lorad-mcmc.conf, not <bycodon_lorad_dir>/lorad.conf).
        #bycodon_lorad_conf_template = open('conf-bycodon-lorad-template.txt','r').read()
        #bycodon_lorad_conf_filename = os.path.join(bycodon_gss_dir,'lorad-mcmc.conf')
        #bycodon_lorad_conf_contents = re.sub('__FIRST_SITE_1ST_CODON__', '1',                            bycodon_lorad_conf_template, re.M | re.S)
        #bycodon_lorad_conf_contents = re.sub('__LAST_SITE_1ST_CODON__',  str(bycodon_boundaries[0]),     bycodon_lorad_conf_contents, re.M | re.S)
        #bycodon_lorad_conf_contents = re.sub('__FIRST_SITE_2ND_CODON__', str(bycodon_boundaries[0] + 1), bycodon_lorad_conf_contents, re.M | re.S)
        #bycodon_lorad_conf_contents = re.sub('__LAST_SITE_2ND_CODON__',  str(bycodon_boundaries[1]),     bycodon_lorad_conf_contents, re.M | re.S)
        #bycodon_lorad_conf_contents = re.sub('__FIRST_SITE_3RD_CODON__', str(bycodon_boundaries[1] + 1), bycodon_lorad_conf_contents, re.M | re.S)
        #bycodon_lorad_conf_contents = re.sub('__LAST_SITE_3RD_CODON__',  str(bycodon_boundaries[2]),     bycodon_lorad_conf_contents, re.M | re.S)
        #bycodon_lorad_conf_contents = re.sub('__BURNIN__',               gss_burnin,                     bycodon_lorad_conf_contents, re.M | re.S)
        #bycodon_lorad_conf_contents = re.sub('__NITER__',                gss_niter,                      bycodon_lorad_conf_contents, re.M | re.S)
        #bycodon_lorad_conf_contents = re.sub('__SAMPLEFREQ__',           gss_samplefreq,                 bycodon_lorad_conf_contents, re.M | re.S)
        #bycodon_lorad_conf_contents = re.sub('__PRINTFREQ__',            gss_printfreq,                  bycodon_lorad_conf_contents, re.M | re.S)
        #bycodon_lorad_conf_contents = re.sub('__COVERAGE10__',           lorad_coverage10,               bycodon_lorad_conf_contents, re.M | re.S)
        #bycodon_lorad_conf_contents = re.sub('__COVERAGE20__',           lorad_coverage20,               bycodon_lorad_conf_contents, re.M | re.S)
        #bycodon_lorad_conf_contents = re.sub('__COVERAGE30__',           lorad_coverage30,               bycodon_lorad_conf_contents, re.M | re.S)
        #bycodon_lorad_conf_contents = re.sub('__COVERAGE40__',           lorad_coverage40,               bycodon_lorad_conf_contents, re.M | re.S)
        #bycodon_lorad_conf_contents = re.sub('__COVERAGE50__',           lorad_coverage50,               bycodon_lorad_conf_contents, re.M | re.S)
        #bycodon_lorad_conf_contents = re.sub('__COVERAGE60__',           lorad_coverage60,               bycodon_lorad_conf_contents, re.M | re.S)
        #bycodon_lorad_conf_contents = re.sub('__COVERAGE70__',           lorad_coverage70,               bycodon_lorad_conf_contents, re.M | re.S)
        #bycodon_lorad_conf_contents = re.sub('__COVERAGE80__',           lorad_coverage80,               bycodon_lorad_conf_contents, re.M | re.S)
        #bycodon_lorad_conf_contents = re.sub('__COVERAGE90__',           lorad_coverage90,               bycodon_lorad_conf_contents, re.M | re.S)
        #bycodon_lorad_conf_contents = re.sub('__COVERAGE95__',           lorad_coverage95,               bycodon_lorad_conf_contents, re.M | re.S)
        #bycodon_lorad_conf_contents = re.sub('__COVERAGE99__',           lorad_coverage99,               bycodon_lorad_conf_contents, re.M | re.S)
        #bycodon_lorad_conf_contents = re.sub('__USE_REGRESSION__',       lorad_regression,               bycodon_lorad_conf_contents, re.M | re.S)
        #bycodon_lorad_conf_contents = re.sub('__LINEAR_REGRESSION__',    lorad_linearregression,         bycodon_lorad_conf_contents, re.M | re.S)
        #bycodon_lorad_conf_contents = re.sub('__RNSEED__',               rnseed,                         bycodon_lorad_conf_contents, re.M | re.S)
        #bycodon_lorad_conf_contents = re.sub('__TREEFILE__',             tree_file_name,                 bycodon_lorad_conf_contents, re.M | re.S)
        #f = open(bycodon_lorad_conf_filename,'w')
        #f.write(bycodon_lorad_conf_contents)
        #f.close()

        bycodon_gss_conf_template = open('conf-bycodon-gss-template.txt','r').read()
        bycodon_gss_conf_filename = os.path.join(bycodon_gss_dir,'lorad-gss.conf')
        bycodon_gss_conf_contents = re.sub('__FIRST_SITE_1ST_CODON__', '1',                            bycodon_gss_conf_template, re.M | re.S)
        bycodon_gss_conf_contents = re.sub('__LAST_SITE_1ST_CODON__',  str(bycodon_boundaries[0]),     bycodon_gss_conf_contents, re.M | re.S)
        bycodon_gss_conf_contents = re.sub('__FIRST_SITE_2ND_CODON__', str(bycodon_boundaries[0] + 1), bycodon_gss_conf_contents, re.M | re.S)
        bycodon_gss_conf_contents = re.sub('__LAST_SITE_2ND_CODON__',  str(bycodon_boundaries[1]),     bycodon_gss_conf_contents, re.M | re.S)
        bycodon_gss_conf_contents = re.sub('__FIRST_SITE_3RD_CODON__', str(bycodon_boundaries[1] + 1), bycodon_gss_conf_contents, re.M | re.S)
        bycodon_gss_conf_contents = re.sub('__LAST_SITE_3RD_CODON__',  str(bycodon_boundaries[2]),     bycodon_gss_conf_contents, re.M | re.S)
        bycodon_gss_conf_contents = re.sub('__BURNIN__',               gss_burnin,                     bycodon_gss_conf_contents, re.M | re.S)
        bycodon_gss_conf_contents = re.sub('__NITER__',                gss_niter,                      bycodon_gss_conf_contents, re.M | re.S)
        bycodon_gss_conf_contents = re.sub('__SAMPLEFREQ__',           gss_samplefreq,                 bycodon_gss_conf_contents, re.M | re.S)
        bycodon_gss_conf_contents = re.sub('__PRINTFREQ__',            gss_printfreq,                  bycodon_gss_conf_contents, re.M | re.S)
        bycodon_gss_conf_contents = re.sub('__NSTONES__',              gss_nstones,                    bycodon_gss_conf_contents, re.M | re.S)
        bycodon_gss_conf_contents = re.sub('__ALPHA__',                gss_alpha,                      bycodon_gss_conf_contents, re.M | re.S)
        bycodon_gss_conf_contents = re.sub('__RNSEED__',               rnseed,                         bycodon_gss_conf_contents, re.M | re.S)
        bycodon_gss_conf_contents = re.sub('__TREEFILE__',             tree_file_name,                 bycodon_gss_conf_contents, re.M | re.S)
        f = open(bycodon_gss_conf_filename,'w')
        f.write(bycodon_gss_conf_contents)
        f.close()

    if include_rev:
        bycodon_rev_template     = open('rev-bycodon-template.txt','r').read()
        bycodon_rev_filename     = os.path.join(bycodon_rev_dir,'ss.Rev')
        bycodon_rev_contents     = re.sub('__RNSEED__',         rnseed,              bycodon_rev_template, re.M | re.S)
        bycodon_rev_contents     = re.sub('__BURNIN__',         rev_burnin,          bycodon_rev_contents, re.M | re.S)
        bycodon_rev_contents     = re.sub('__NITER__',          rev_niter,           bycodon_rev_contents, re.M | re.S)
        bycodon_rev_contents     = re.sub('__SAMPLEFREQ__',     rev_samplefreq,      bycodon_rev_contents, re.M | re.S)
        bycodon_rev_contents     = re.sub('__PRINTFREQ__',      rev_printfreq,       bycodon_rev_contents, re.M | re.S)
        bycodon_rev_contents     = re.sub('__TUNINGINTERVAL__', rev_tuning_interval, bycodon_rev_contents, re.M | re.S)
        bycodon_rev_contents     = re.sub('__NSTONES__',        rev_nstones,         bycodon_rev_contents, re.M | re.S)
        bycodon_rev_contents     = re.sub('__ALPHA__',          rev_alpha,           bycodon_rev_contents, re.M | re.S)
        bycodon_rev_contents     = re.sub('__TREEFILE__',       tree_file_name,      bycodon_rev_contents, re.M | re.S)
        if fan_etal_2011:
            bycodon_rev_contents = re.sub('__FAN_ETAL_2011__',  'TRUE',              bycodon_rev_contents, re.M | re.S)
        else:
            bycodon_rev_contents = re.sub('__FAN_ETAL_2011__',  'FALSE',             bycodon_rev_contents, re.M | re.S)
        f = open(bycodon_rev_filename,'w')
        f.write(bycodon_rev_contents)
        f.close()

#################################################
# Create lorad.conf scripts for bygene analyses #
#################################################

if include_bygene:
    if not fan_etal_2011 and (include_lorad or include_ghm):
        bygene_lorad_conf_template = open('conf-bygene-lorad-template.txt','r').read()
        bygene_lorad_conf_filename = os.path.join(bygene_lorad_dir,'lorad.conf')
        bygene_lorad_conf_contents = re.sub('__FIRST_SITE_COI__',     '1',                           bygene_lorad_conf_template, re.M | re.S)
        bygene_lorad_conf_contents = re.sub('__LAST_SITE_COI__',      str(bygene_boundaries[0]),     bygene_lorad_conf_contents, re.M | re.S)
        bygene_lorad_conf_contents = re.sub('__FIRST_SITE_COII__',    str(bygene_boundaries[0] + 1), bygene_lorad_conf_contents, re.M | re.S)
        bygene_lorad_conf_contents = re.sub('__LAST_SITE_COII__',     str(bygene_boundaries[1]),     bygene_lorad_conf_contents, re.M | re.S)
        bygene_lorad_conf_contents = re.sub('__FIRST_SITE_ATPASE6__', str(bygene_boundaries[1] + 1), bygene_lorad_conf_contents, re.M | re.S)
        bygene_lorad_conf_contents = re.sub('__LAST_SITE_ATPASE6__',  str(bygene_boundaries[2]),     bygene_lorad_conf_contents, re.M | re.S)
        bygene_lorad_conf_contents = re.sub('__FIRST_SITE_ATPASE8__', str(bygene_boundaries[2] + 1), bygene_lorad_conf_contents, re.M | re.S)
        bygene_lorad_conf_contents = re.sub('__LAST_SITE_ATPASE8__',  str(bygene_boundaries[3]),     bygene_lorad_conf_contents, re.M | re.S)
        bygene_lorad_conf_contents = re.sub('__BURNIN__',             lorad_burnin,                  bygene_lorad_conf_contents, re.M | re.S)
        bygene_lorad_conf_contents = re.sub('__NITER__',              lorad_niter,                   bygene_lorad_conf_contents, re.M | re.S)
        bygene_lorad_conf_contents = re.sub('__SAMPLEFREQ__',         lorad_samplefreq,              bygene_lorad_conf_contents, re.M | re.S)
        bygene_lorad_conf_contents = re.sub('__PRINTFREQ__',          lorad_printfreq,               bygene_lorad_conf_contents, re.M | re.S)
        bygene_lorad_conf_contents = re.sub('__COVERAGE10__',         lorad_coverage10,              bygene_lorad_conf_contents, re.M | re.S)
        bygene_lorad_conf_contents = re.sub('__COVERAGE20__',         lorad_coverage20,              bygene_lorad_conf_contents, re.M | re.S)
        bygene_lorad_conf_contents = re.sub('__COVERAGE30__',         lorad_coverage30,              bygene_lorad_conf_contents, re.M | re.S)
        bygene_lorad_conf_contents = re.sub('__COVERAGE40__',         lorad_coverage40,              bygene_lorad_conf_contents, re.M | re.S)
        bygene_lorad_conf_contents = re.sub('__COVERAGE50__',         lorad_coverage50,              bygene_lorad_conf_contents, re.M | re.S)
        bygene_lorad_conf_contents = re.sub('__COVERAGE60__',         lorad_coverage60,              bygene_lorad_conf_contents, re.M | re.S)
        bygene_lorad_conf_contents = re.sub('__COVERAGE70__',         lorad_coverage70,              bygene_lorad_conf_contents, re.M | re.S)
        bygene_lorad_conf_contents = re.sub('__COVERAGE80__',         lorad_coverage80,              bygene_lorad_conf_contents, re.M | re.S)
        bygene_lorad_conf_contents = re.sub('__COVERAGE90__',         lorad_coverage90,              bygene_lorad_conf_contents, re.M | re.S)
        bygene_lorad_conf_contents = re.sub('__COVERAGE95__',         lorad_coverage95,              bygene_lorad_conf_contents, re.M | re.S)
        bygene_lorad_conf_contents = re.sub('__COVERAGE99__',         lorad_coverage99,              bygene_lorad_conf_contents, re.M | re.S)
        bygene_lorad_conf_contents = re.sub('__RNSEED__',             rnseed,                        bygene_lorad_conf_contents, re.M | re.S)
        bygene_lorad_conf_contents = re.sub('__TREEFILE__',           tree_file_name,                bygene_lorad_conf_contents, re.M | re.S)
        f = open(bygene_lorad_conf_filename,'w')
        f.write(bygene_lorad_conf_contents)
        f.close()

    if not fan_etal_2011 and include_gss:
        # This part commented out because, in the paper, GSS used reference distributions computed from the posterior sample created for LoRaD.
        ## GSS requires a reference distribution estimated from an initial MCMC run, so this part is identical to the include_lorad
        ## part above except that gss_burnin, gss_niter, gss_samplefreq, and gss_printfreq are used instead of lorad_burnin, lorad_niter,
        ## lorad_samplefreq, and lorad_printfreq (also, the file generated is <bygene_gss_dir>/lorad-mcmc.conf, not <bygene_lorad_dir>/lorad.conf).
        #bygene_lorad_conf_template = open('conf-bygene-lorad-template.txt','r').read()
        #bygene_lorad_conf_filename = os.path.join(bygene_gss_dir,'lorad-mcmc.conf')
        #bygene_lorad_conf_contents = re.sub('__FIRST_SITE_COI__',     '1',                           bygene_lorad_conf_template, re.M | re.S)
        #bygene_lorad_conf_contents = re.sub('__LAST_SITE_COI__',      str(bygene_boundaries[0]),     bygene_lorad_conf_contents, re.M | re.S)
        #bygene_lorad_conf_contents = re.sub('__FIRST_SITE_COII__',    str(bygene_boundaries[0] + 1), bygene_lorad_conf_contents, re.M | re.S)
        #bygene_lorad_conf_contents = re.sub('__LAST_SITE_COII__',     str(bygene_boundaries[1]),     bygene_lorad_conf_contents, re.M | re.S)
        #bygene_lorad_conf_contents = re.sub('__FIRST_SITE_ATPASE6__', str(bygene_boundaries[1] + 1), bygene_lorad_conf_contents, re.M | re.S)
        #bygene_lorad_conf_contents = re.sub('__LAST_SITE_ATPASE6__',  str(bygene_boundaries[2]),     bygene_lorad_conf_contents, re.M | re.S)
        #bygene_lorad_conf_contents = re.sub('__FIRST_SITE_ATPASE8__', str(bygene_boundaries[2] + 1), bygene_lorad_conf_contents, re.M | re.S)
        #bygene_lorad_conf_contents = re.sub('__LAST_SITE_ATPASE8__',  str(bygene_boundaries[3]),     bygene_lorad_conf_contents, re.M | re.S)
        #bygene_lorad_conf_contents = re.sub('__BURNIN__',             gss_burnin,                    bygene_lorad_conf_contents, re.M | re.S)
        #bygene_lorad_conf_contents = re.sub('__NITER__',              gss_niter,                     bygene_lorad_conf_contents, re.M | re.S)
        #bygene_lorad_conf_contents = re.sub('__SAMPLEFREQ__',         gss_samplefreq,                bygene_lorad_conf_contents, re.M | re.S)
        #bygene_lorad_conf_contents = re.sub('__PRINTFREQ__',          gss_printfreq,                 bygene_lorad_conf_contents, re.M | re.S)
        #bygene_lorad_conf_contents = re.sub('__COVERAGE10__',         lorad_coverage10,              bygene_lorad_conf_contents, re.M | re.S)
        #bygene_lorad_conf_contents = re.sub('__COVERAGE20__',         lorad_coverage20,              bygene_lorad_conf_contents, re.M | re.S)
        #bygene_lorad_conf_contents = re.sub('__COVERAGE30__',         lorad_coverage30,              bygene_lorad_conf_contents, re.M | re.S)
        #bygene_lorad_conf_contents = re.sub('__COVERAGE40__',         lorad_coverage40,              bygene_lorad_conf_contents, re.M | re.S)
        #bygene_lorad_conf_contents = re.sub('__COVERAGE50__',         lorad_coverage50,              bygene_lorad_conf_contents, re.M | re.S)
        #bygene_lorad_conf_contents = re.sub('__COVERAGE60__',         lorad_coverage60,              bygene_lorad_conf_contents, re.M | re.S)
        #bygene_lorad_conf_contents = re.sub('__COVERAGE70__',         lorad_coverage70,              bygene_lorad_conf_contents, re.M | re.S)
        #bygene_lorad_conf_contents = re.sub('__COVERAGE80__',         lorad_coverage80,              bygene_lorad_conf_contents, re.M | re.S)
        #bygene_lorad_conf_contents = re.sub('__COVERAGE90__',         lorad_coverage90,              bygene_lorad_conf_contents, re.M | re.S)
        #bygene_lorad_conf_contents = re.sub('__COVERAGE95__',         lorad_coverage95,              bygene_lorad_conf_contents, re.M | re.S)
        #bygene_lorad_conf_contents = re.sub('__COVERAGE99__',         lorad_coverage99,              bygene_lorad_conf_contents, re.M | re.S)
        #bygene_lorad_conf_contents = re.sub('__USE_REGRESSION__',     lorad_regression,              bygene_lorad_conf_contents, re.M | re.S)
        #bygene_lorad_conf_contents = re.sub('__LINEAR_REGRESSION__',  lorad_linearregression,        bygene_lorad_conf_contents, re.M | re.S)
        #bygene_lorad_conf_contents = re.sub('__RNSEED__',             rnseed,                        bygene_lorad_conf_contents, re.M | re.S)
        #bygene_lorad_conf_contents = re.sub('__TREEFILE__',           tree_file_name,                bygene_lorad_conf_contents, re.M | re.S)
        #f = open(bygene_lorad_conf_filename,'w')
        #f.write(bygene_lorad_conf_contents)
        #f.close()

        bygene_gss_conf_template = open('conf-bygene-gss-template.txt','r').read()
        bygene_gss_conf_filename = os.path.join(bygene_gss_dir,'lorad-gss.conf')
        bygene_gss_conf_contents = re.sub('__FIRST_SITE_COI__',     '1',                           bygene_gss_conf_template, re.M | re.S)
        bygene_gss_conf_contents = re.sub('__LAST_SITE_COI__',      str(bygene_boundaries[0]),     bygene_gss_conf_contents, re.M | re.S)
        bygene_gss_conf_contents = re.sub('__FIRST_SITE_COII__',    str(bygene_boundaries[0] + 1), bygene_gss_conf_contents, re.M | re.S)
        bygene_gss_conf_contents = re.sub('__LAST_SITE_COII__',     str(bygene_boundaries[1]),     bygene_gss_conf_contents, re.M | re.S)
        bygene_gss_conf_contents = re.sub('__FIRST_SITE_ATPASE6__', str(bygene_boundaries[1] + 1), bygene_gss_conf_contents, re.M | re.S)
        bygene_gss_conf_contents = re.sub('__LAST_SITE_ATPASE6__',  str(bygene_boundaries[2]),     bygene_gss_conf_contents, re.M | re.S)
        bygene_gss_conf_contents = re.sub('__FIRST_SITE_ATPASE8__', str(bygene_boundaries[2] + 1), bygene_gss_conf_contents, re.M | re.S)
        bygene_gss_conf_contents = re.sub('__LAST_SITE_ATPASE8__',  str(bygene_boundaries[3]),     bygene_gss_conf_contents, re.M | re.S)
        bygene_gss_conf_contents = re.sub('__BURNIN__',             gss_burnin,                    bygene_gss_conf_contents, re.M | re.S)
        bygene_gss_conf_contents = re.sub('__NITER__',              gss_niter,                     bygene_gss_conf_contents, re.M | re.S)
        bygene_gss_conf_contents = re.sub('__SAMPLEFREQ__',         gss_samplefreq,                bygene_gss_conf_contents, re.M | re.S)
        bygene_gss_conf_contents = re.sub('__PRINTFREQ__',          gss_printfreq,                 bygene_gss_conf_contents, re.M | re.S)
        bygene_gss_conf_contents = re.sub('__NSTONES__',            gss_nstones,                   bygene_gss_conf_contents, re.M | re.S)
        bygene_gss_conf_contents = re.sub('__ALPHA__',              gss_alpha,                     bygene_gss_conf_contents, re.M | re.S)
        bygene_gss_conf_contents = re.sub('__RNSEED__',             rnseed,                        bygene_gss_conf_contents, re.M | re.S)
        bygene_gss_conf_contents = re.sub('__TREEFILE__',           tree_file_name,                bygene_gss_conf_contents, re.M | re.S)
        f = open(bygene_gss_conf_filename,'w')
        f.write(bygene_gss_conf_contents)
        f.close()

    if include_rev:
        bygene_rev_template     = open('rev-bygene-template.txt','r').read()
        bygene_rev_filename     = os.path.join(bygene_rev_dir,'ss.Rev')
        bygene_rev_contents     = re.sub('__RNSEED__',         rnseed,              bygene_rev_template, re.M | re.S)
        bygene_rev_contents     = re.sub('__BURNIN__',         rev_burnin,          bygene_rev_contents, re.M | re.S)
        bygene_rev_contents     = re.sub('__NITER__',          rev_niter,           bygene_rev_contents, re.M | re.S)
        bygene_rev_contents     = re.sub('__SAMPLEFREQ__',     rev_samplefreq,      bygene_rev_contents, re.M | re.S)
        bygene_rev_contents     = re.sub('__PRINTFREQ__',      rev_printfreq,       bygene_rev_contents, re.M | re.S)
        bygene_rev_contents     = re.sub('__TUNINGINTERVAL__', rev_tuning_interval, bygene_rev_contents, re.M | re.S)
        bygene_rev_contents     = re.sub('__NSTONES__',        rev_nstones,         bygene_rev_contents, re.M | re.S)
        bygene_rev_contents     = re.sub('__ALPHA__',          rev_alpha,           bygene_rev_contents, re.M | re.S)
        bygene_rev_contents     = re.sub('__TREEFILE__',       tree_file_name,      bygene_rev_contents, re.M | re.S)
        if fan_etal_2011:
            bygene_rev_contents = re.sub('__FAN_ETAL_2011__',  'TRUE',              bygene_rev_contents, re.M | re.S)
        else:
            bygene_rev_contents = re.sub('__FAN_ETAL_2011__',  'FALSE',             bygene_rev_contents, re.M | re.S)
        f = open(bygene_rev_filename,'w')
        f.write(bygene_rev_contents)
        f.close()

#################################################
# Create lorad.conf scripts for byboth analyses #
#################################################

if include_byboth:
    if not fan_etal_2011 and (include_lorad or include_ghm):
        byboth_lorad_conf_template = open('conf-byboth-lorad-template.txt','r').read()
        byboth_lorad_conf_filename = os.path.join(byboth_lorad_dir,'lorad.conf')
        byboth_lorad_conf_contents = re.sub('__FIRST_SITE_COI1__',     '1',                            byboth_lorad_conf_template, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__LAST_SITE_COI1__',      str(byboth_boundaries[0]),      byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__FIRST_SITE_COI2__',     str(byboth_boundaries[0] + 1),  byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__LAST_SITE_COI2__',      str(byboth_boundaries[1]),      byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__FIRST_SITE_COI3__',     str(byboth_boundaries[1] + 1),  byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__LAST_SITE_COI3__',      str(byboth_boundaries[2]),      byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__FIRST_SITE_COII1__',    str(byboth_boundaries[2] + 1),  byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__LAST_SITE_COII1__',     str(byboth_boundaries[3]),      byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__FIRST_SITE_COII2__',    str(byboth_boundaries[3] + 1),  byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__LAST_SITE_COII2__',     str(byboth_boundaries[4]),      byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__FIRST_SITE_COII3__',    str(byboth_boundaries[4] + 1),  byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__LAST_SITE_COII3__',     str(byboth_boundaries[5]),      byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__FIRST_SITE_ATPASE61__', str(byboth_boundaries[5] + 1),  byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__LAST_SITE_ATPASE61__',  str(byboth_boundaries[6]),      byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__FIRST_SITE_ATPASE62__', str(byboth_boundaries[6] + 1),  byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__LAST_SITE_ATPASE62__',  str(byboth_boundaries[7]),      byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__FIRST_SITE_ATPASE63__', str(byboth_boundaries[7] + 1),  byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__LAST_SITE_ATPASE63__',  str(byboth_boundaries[8]),      byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__FIRST_SITE_ATPASE81__', str(byboth_boundaries[8] + 1),  byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__LAST_SITE_ATPASE81__',  str(byboth_boundaries[9]),      byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__FIRST_SITE_ATPASE82__', str(byboth_boundaries[9] + 1),  byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__LAST_SITE_ATPASE82__',  str(byboth_boundaries[10]),     byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__FIRST_SITE_ATPASE83__', str(byboth_boundaries[10] + 1), byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__LAST_SITE_ATPASE83__',  str(byboth_boundaries[11]),     byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__BURNIN__',              lorad_burnin,                   byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__NITER__',               lorad_niter,                    byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__SAMPLEFREQ__',          lorad_samplefreq,               byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__PRINTFREQ__',           lorad_printfreq,                byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__COVERAGE10__',          lorad_coverage10,               byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__COVERAGE20__',          lorad_coverage20,               byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__COVERAGE30__',          lorad_coverage30,               byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__COVERAGE40__',          lorad_coverage40,               byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__COVERAGE50__',          lorad_coverage50,               byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__COVERAGE60__',          lorad_coverage60,               byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__COVERAGE70__',          lorad_coverage70,               byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__COVERAGE80__',          lorad_coverage80,               byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__COVERAGE90__',          lorad_coverage90,               byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__COVERAGE95__',          lorad_coverage95,               byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__COVERAGE99__',          lorad_coverage99,               byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__RNSEED__',              rnseed,                         byboth_lorad_conf_contents, re.M | re.S)
        byboth_lorad_conf_contents = re.sub('__TREEFILE__',            tree_file_name,                 byboth_lorad_conf_contents, re.M | re.S)
        f = open(byboth_lorad_conf_filename,'w')
        f.write(byboth_lorad_conf_contents)
        f.close()

    if not fan_etal_2011 and include_gss:
        # This part commented out because, in the paper, GSS used reference distributions computed from the posterior sample created for LoRaD.
        ## GSS requires a reference distribution estimated from an initial MCMC run, so this part is identical to the include_lorad
        ## part above except that gss_burnin, gss_niter, gss_samplefreq, and gss_printfreq are used instead of lorad_burnin, lorad_niter,
        ## lorad_samplefreq, and lorad_printfreq (also, the file generated is <byboth_gss_dir>/lorad-mcmc.conf, not <byboth_lorad_dir>/lorad.conf).
        #byboth_lorad_conf_template = open('conf-byboth-lorad-template.txt','r').read()
        #byboth_lorad_conf_filename = os.path.join(byboth_gss_dir,'lorad-mcmc.conf')
        #byboth_lorad_conf_contents = re.sub('__FIRST_SITE_COI1__',     '1',                            byboth_lorad_conf_template, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__LAST_SITE_COI1__',      str(byboth_boundaries[0]),      byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__FIRST_SITE_COI2__',     str(byboth_boundaries[0] + 1),  byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__LAST_SITE_COI2__',      str(byboth_boundaries[1]),      byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__FIRST_SITE_COI3__',     str(byboth_boundaries[1] + 1),  byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__LAST_SITE_COI3__',      str(byboth_boundaries[2]),      byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__FIRST_SITE_COII1__',    str(byboth_boundaries[2] + 1),  byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__LAST_SITE_COII1__',     str(byboth_boundaries[3]),      byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__FIRST_SITE_COII2__',    str(byboth_boundaries[3] + 1),  byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__LAST_SITE_COII2__',     str(byboth_boundaries[4]),      byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__FIRST_SITE_COII3__',    str(byboth_boundaries[4] + 1),  byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__LAST_SITE_COII3__',     str(byboth_boundaries[5]),      byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__FIRST_SITE_ATPASE61__', str(byboth_boundaries[5] + 1),  byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__LAST_SITE_ATPASE61__',  str(byboth_boundaries[6]),      byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__FIRST_SITE_ATPASE62__', str(byboth_boundaries[6] + 1),  byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__LAST_SITE_ATPASE62__',  str(byboth_boundaries[7]),      byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__FIRST_SITE_ATPASE63__', str(byboth_boundaries[7] + 1),  byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__LAST_SITE_ATPASE63__',  str(byboth_boundaries[8]),      byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__FIRST_SITE_ATPASE81__', str(byboth_boundaries[8] + 1),  byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__LAST_SITE_ATPASE81__',  str(byboth_boundaries[9]),      byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__FIRST_SITE_ATPASE82__', str(byboth_boundaries[9] + 1),  byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__LAST_SITE_ATPASE82__',  str(byboth_boundaries[10]),     byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__FIRST_SITE_ATPASE83__', str(byboth_boundaries[10] + 1), byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__LAST_SITE_ATPASE83__',  str(byboth_boundaries[11]),     byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__BURNIN__',              gss_burnin,                     byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__NITER__',               gss_niter,                      byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__SAMPLEFREQ__',          gss_samplefreq,                 byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__PRINTFREQ__',           gss_printfreq,                  byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__COVERAGE10__',          lorad_coverage10,               byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__COVERAGE20__',          lorad_coverage20,               byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__COVERAGE30__',          lorad_coverage30,               byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__COVERAGE40__',          lorad_coverage40,               byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__COVERAGE50__',          lorad_coverage50,               byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__COVERAGE60__',          lorad_coverage60,               byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__COVERAGE70__',          lorad_coverage70,               byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__COVERAGE80__',          lorad_coverage80,               byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__COVERAGE90__',          lorad_coverage90,               byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__COVERAGE95__',          lorad_coverage95,               byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__COVERAGE99__',          lorad_coverage99,               byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__RNSEED__',              rnseed,                         byboth_lorad_conf_contents, re.M | re.S)
        #byboth_lorad_conf_contents = re.sub('__TREEFILE__',            tree_file_name,                 byboth_lorad_conf_contents, re.M | re.S)
        #f = open(byboth_lorad_conf_filename,'w')
        #f.write(byboth_lorad_conf_contents)
        #f.close()

        byboth_gss_conf_template = open('conf-byboth-gss-template.txt','r').read()
        byboth_gss_conf_filename = os.path.join(byboth_gss_dir,'lorad-gss.conf')
        byboth_gss_conf_contents = re.sub('__FIRST_SITE_COI1__',     '1',                            byboth_gss_conf_template, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__LAST_SITE_COI1__',      str(byboth_boundaries[0]),      byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__FIRST_SITE_COI2__',     str(byboth_boundaries[0] + 1),  byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__LAST_SITE_COI2__',      str(byboth_boundaries[1]),      byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__FIRST_SITE_COI3__',     str(byboth_boundaries[1] + 1),  byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__LAST_SITE_COI3__',      str(byboth_boundaries[2]),      byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__FIRST_SITE_COII1__',    str(byboth_boundaries[2] + 1),  byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__LAST_SITE_COII1__',     str(byboth_boundaries[3]),      byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__FIRST_SITE_COII2__',    str(byboth_boundaries[3] + 1),  byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__LAST_SITE_COII2__',     str(byboth_boundaries[4]),      byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__FIRST_SITE_COII3__',    str(byboth_boundaries[4] + 1),  byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__LAST_SITE_COII3__',     str(byboth_boundaries[5]),      byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__FIRST_SITE_ATPASE61__', str(byboth_boundaries[5] + 1),  byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__LAST_SITE_ATPASE61__',  str(byboth_boundaries[6]),      byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__FIRST_SITE_ATPASE62__', str(byboth_boundaries[6] + 1),  byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__LAST_SITE_ATPASE62__',  str(byboth_boundaries[7]),      byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__FIRST_SITE_ATPASE63__', str(byboth_boundaries[7] + 1),  byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__LAST_SITE_ATPASE63__',  str(byboth_boundaries[8]),      byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__FIRST_SITE_ATPASE81__', str(byboth_boundaries[8] + 1),  byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__LAST_SITE_ATPASE81__',  str(byboth_boundaries[9]),      byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__FIRST_SITE_ATPASE82__', str(byboth_boundaries[9] + 1),  byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__LAST_SITE_ATPASE82__',  str(byboth_boundaries[10]),     byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__FIRST_SITE_ATPASE83__', str(byboth_boundaries[10] + 1), byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__LAST_SITE_ATPASE83__',  str(byboth_boundaries[11]),     byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__BURNIN__',              gss_burnin,                     byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__NITER__',               gss_niter,                      byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__SAMPLEFREQ__',          gss_samplefreq,                 byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__PRINTFREQ__',           gss_printfreq,                  byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__NSTONES__',             gss_nstones,                    byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__ALPHA__',               gss_alpha,                      byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__RNSEED__',              rnseed,                         byboth_gss_conf_contents, re.M | re.S)
        byboth_gss_conf_contents = re.sub('__TREEFILE__',            tree_file_name,                 byboth_gss_conf_contents, re.M | re.S)
        f = open(byboth_gss_conf_filename,'w')
        f.write(byboth_gss_conf_contents)
        f.close()

    if include_rev:
        byboth_rev_template     = open('rev-byboth-template.txt','r').read()
        byboth_rev_filename     = os.path.join(byboth_rev_dir,'ss.Rev')
        byboth_rev_contents     = re.sub('__RNSEED__',         rnseed,              byboth_rev_template, re.M | re.S)
        byboth_rev_contents     = re.sub('__BURNIN__',         rev_burnin,          byboth_rev_contents, re.M | re.S)
        byboth_rev_contents     = re.sub('__NITER__',          rev_niter,           byboth_rev_contents, re.M | re.S)
        byboth_rev_contents     = re.sub('__SAMPLEFREQ__',     rev_samplefreq,      byboth_rev_contents, re.M | re.S)
        byboth_rev_contents     = re.sub('__PRINTFREQ__',      rev_printfreq,       byboth_rev_contents, re.M | re.S)
        byboth_rev_contents     = re.sub('__TUNINGINTERVAL__', rev_tuning_interval, byboth_rev_contents, re.M | re.S)
        byboth_rev_contents     = re.sub('__NSTONES__',        rev_nstones,         byboth_rev_contents, re.M | re.S)
        byboth_rev_contents     = re.sub('__ALPHA__',          rev_alpha,           byboth_rev_contents, re.M | re.S)
        byboth_rev_contents     = re.sub('__TREEFILE__',       tree_file_name,      byboth_rev_contents, re.M | re.S)
        if fan_etal_2011:
            byboth_rev_contents = re.sub('__FAN_ETAL_2011__',  'TRUE',              byboth_rev_contents, re.M | re.S)
        else:
            byboth_rev_contents = re.sub('__FAN_ETAL_2011__',  'FALSE',             byboth_rev_contents, re.M | re.S)
        f = open(byboth_rev_filename,'w')
        f.write(byboth_rev_contents)
        f.close()

#################################################################################
#################################################################################
### Create loradml-untransformed.conf and loradml-logtransformed.conf scripts ###
#################################################################################
#################################################################################

#####################################################
# Create loradml-*.conf scripts for unpart analyses #
#####################################################

if include_unpart and (include_lorad or include_ghm):
    ####################################################
    # Copy the loradml-unpart-logtransformed.conf file #
    ####################################################
    file_contents = open('loradml-unpart-logtransformed.conf', 'r').read()
    file_name = os.path.join(unpart_lorad_dir,'loradml-logtransformed.conf')
    f = open(file_name,'w')
    f.write(file_contents)
    f.close()
    
    ###################################################
    # Copy the loradml-unpart-untransformed.conf file #
    ###################################################
    file_contents = open('loradml-unpart-untransformed.conf', 'r').read()
    file_name = os.path.join(unpart_lorad_dir,'loradml-untransformed.conf')
    f = open(file_name,'w')
    f.write(file_contents)
    f.close()

#####################################################
# Create loradml-*.conf scripts for bycodon analyses #
#####################################################

if include_bycodon and (include_lorad or include_ghm):
    ####################################################
    # Copy the loradml-bycodon-logtransformed.conf file #
    ####################################################
    file_contents = open('loradml-bycodon-logtransformed.conf', 'r').read()
    file_name = os.path.join(bycodon_lorad_dir,'loradml-logtransformed.conf')
    f = open(file_name,'w')
    f.write(file_contents)
    f.close()
    
    ###################################################
    # Copy the loradml-bycodon-untransformed.conf file #
    ###################################################
    file_contents = open('loradml-bycodon-untransformed.conf', 'r').read()
    file_name = os.path.join(bycodon_lorad_dir,'loradml-untransformed.conf')
    f = open(file_name,'w')
    f.write(file_contents)
    f.close()

#####################################################
# Create loradml-*.conf scripts for bygene analyses #
#####################################################

if include_bygene and (include_lorad or include_ghm):
    ####################################################
    # Copy the loradml-bygene-logtransformed.conf file #
    ####################################################
    file_contents = open('loradml-bygene-logtransformed.conf', 'r').read()
    file_name = os.path.join(bygene_lorad_dir,'loradml-logtransformed.conf')
    f = open(file_name,'w')
    f.write(file_contents)
    f.close()
    
    ###################################################
    # Copy the loradml-bygene-untransformed.conf file #
    ###################################################
    file_contents = open('loradml-bygene-untransformed.conf', 'r').read()
    file_name = os.path.join(bygene_lorad_dir,'loradml-untransformed.conf')
    f = open(file_name,'w')
    f.write(file_contents)
    f.close()

#####################################################
# Create loradml-*.conf scripts for byboth analyses #
#####################################################

if include_byboth and (include_lorad or include_ghm):
    ####################################################
    # Copy the loradml-byboth-logtransformed.conf file #
    ####################################################
    file_contents = open('loradml-byboth-logtransformed.conf', 'r').read()
    file_name = os.path.join(byboth_lorad_dir,'loradml-logtransformed.conf')
    f = open(file_name,'w')
    f.write(file_contents)
    f.close()
    
    ###################################################
    # Copy the loradml-byboth-untransformed.conf file #
    ###################################################
    file_contents = open('loradml-byboth-untransformed.conf', 'r').read()
    file_name = os.path.join(byboth_lorad_dir,'loradml-untransformed.conf')
    f = open(file_name,'w')
    f.write(file_contents)
    f.close()

##########################
##########################
### Copy the tree file ###
##########################
##########################

if thirtytwo:
    #################################
    # Copy the gtrg-32taxa.tre file #
    #################################
    gtrg_32taxa_tre_contents = open('gtrg-32taxa.tre', 'r').read()
    gtrg_32taxa_tre_filename = os.path.join(dest_dir,'gtrg-32taxa.tre')
    f = open(gtrg_32taxa_tre_filename,'w')
    f.write(gtrg_32taxa_tre_contents)
    f.close()
else:
    #################################
    # Copy the gtrg-31taxa.tre file #
    #################################
    gtrg_31taxa_tre_contents = open('gtrg-31taxa.tre', 'r').read()
    gtrg_31taxa_tre_filename = os.path.join(dest_dir,'gtrg-31taxa.tre')
    f = open(gtrg_31taxa_tre_filename,'w')
    f.write(gtrg_31taxa_tre_contents)
    f.close()
