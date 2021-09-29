# Reads S1679.nex and creates from it these 23 files:
# Unartitioned:
#    1. unpart.nex
# Partitioned by codon:
#    2. bycodon.nex
#    3. codon1st.nex
#    4. codon2nd.nex
#    5. codon3rd.nex
# Partitioned by gene:
#    6. bygene.nex
#    7. COI.nex
#    8. COII.nex
#    9. ATPase6.nex
#   10. ATPase8.nex
# Partitioned by both gene and codon position:
#   11. byboth.nex
#   12. COI-1st.nex
#   13. COI-2nd.nex
#   14. COI-3rd.nex, 
#   15. COII-1st.nex
#   16. COII-2nd.nex
#   17. COII-3rd.nex, 
#   18. ATPase6-1st.nex
#   19. ATPase6-2nd.nex
#   20. ATPase6-3rd.nex
#   21. ATPase8-1st.nex
#   22. ATPase8-2nd.nex
#   23. ATPase8-3rd.nex

# True produces data used in Fan et al. 2011 paper
# False produces data used in Wang et al. 2021 paper
fan_etal_2011 = True

# True uses all 32 taxa
# False excludes Kikihia muta east, whick lacks data for ATPase6 and ATPase8
thirtytwo = True

import sys,os,re

rnseed     = '13579'                   # the pseudorandom number seed to use for all analyses
userid     = 'pol02003'                # the home directory will be assumed to be /home/<userid>
email      = 'paul.o.lewis@gmail.com'  # the email address for notifications
executable = 'yubo/hpdml'              # the relative path to the file to execute with respect to /home/<userid>
revbayes   = 'bin/rb ss.Rev'           # the relative path to the revbayes executable with respect to /home/<userid> plus script to run
dest_dir   = 'g32'                     # directory under which entire directory structure below will be created

excluded_taxa = ['Kikihia muta east']
if thirtytwo:
    excluded_taxa = []
    
tree_file_name = 'gtrg-31taxa.tre'
if thirtytwo:
    tree_file_name = 'gtrg-32taxa.tre'

# Directory structure:
# <dest_dir>
#    unpart 
#       data
#           unpart.nex
#       hpd 
#           hpdml.conf
#           s.sh
#       ss  
#           hpdml.conf
#           s.sh
#       rbss
#           rb.Rev
#           s.sh
#    bycodon
#       data
#           bycodon.nex
#           codon1st.nex
#           codon2nd.nex
#           codon3rd.nex
#       hpd 
#           hpdml.conf
#           s.sh
#       ss  
#           hpdml.conf
#           s.sh
#       rbss
#           rb.Rev
#           s.sh
#    bygene 
#       data
#           bygene.nex
#           COI.nex
#           COII.nex
#           ATPase6.nex
#           ATPase8.nex
#       hpd 
#           hpdml.conf
#           s.sh
#       ss  
#           hpdml.conf
#           s.sh
#       rbss
#           rb.Rev
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
#       hpd 
#           hpdml.conf
#           s.sh
#       ss
#           hpdml.conf
#           s.sh
#       rbss
#           rb.Rev
#           s.sh
#    gtrg-31taxa.tre
#    submit-all.sh  

# This taxon ordering is from the translate statement in gtrg_ml.tre
taxa_in_order = [
  'Kikihia acoustica',
  'Kikihia angusta',
  'Kikihia aotea east',
  'Kikihia aotea west',
  'Kikihia astragali',
  'Kikihia balaena',
  'Kikihia cauta',
  'Kikihia convicta',
  'Kikihia cutora cumberi',
  'Kikihia cutora cutora',
  'Kikihia dugdalei',
  'Kikihia cutora exulis',
  'Kikihia horologium',
  'Kikihia laneorum',
  'Kikihia longula',
  'Kikihia murihikua',
  'Kikihia muta east',
  'Kikihia muta',
  'Kikihia nelsonensis',
  'Kikihia westlandica north',
  'Kikihia ochrina',
  'Kikihia paxillulae',
  'Kikihia peninsularis',
  'Kikihia rosea',
  'Kikihia scutellaris',
  'Kikihia subalpina',
  'Kikihia flemingi',
  'Kikihia westlandica south',
  'Kikihia tasmani',
  'Kikihia tuta',
  'Rhodopsalta microdora',
  'Maoricicada cassiope'
]

genes_used = ['COI','COII','ATPase8','ATPase6']
gene_boundaries = {}
gene_boundaries['COI']     = (1,774)
gene_boundaries['COII']    = (775,1476)
gene_boundaries['tRNA']    = (1477,1538) # 1538 - 1476 = 62
if fan_etal_2011:
    gene_boundaries['ATPase8'] = (1539,1689) # (1477,1627) <-- includes AT from beginning of ATPase6 at end of ATPase8
    gene_boundaries['ATPase6'] = (1690,2152) # (1628,2090) <-- ATPase6 starts with 3rd position G
else:
    gene_boundaries['ATPase8'] = (1539,1687) # (1477,1625)
    gene_boundaries['ATPase6'] = (1691,2152) # (1629,2090) <-- note AT from ATPase6 and G from ATPase6 deleted

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
    # fn is filename to use for the nexus file saved
    # taxa is a list of taxon names serving as indices into the maps stored in the sequences vector
    # masks_vect is a vector holding masks for each element of sequences (a mask is a string containing
    #   '-' for every site that should be included and '*' for every site that should be excluded.
    # nchar_vect is a vector holding the number of sites for each element of sequences
    # sequences_vect is a vector of maps, each with taxa as keys and sequences as values
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
        boundaries.append(ntotal)
    
    ntax = len(taxa) - len(excluded_taxa)
    nchar = sum(nchar_vect)
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
            f.write('\n')
    f.write('  ;\n')
    f.write('end;\n')
    f.close()
    return boundaries

ntax, nchar, mask, taxa, sequences = readNexusFile('S1679.nex')

#for t in taxa:
#    print(t)
#sys.exit('debug stop')

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

unpart_dir       = os.path.join(dest_dir, 'unpart')
unpart_data_dir  = os.path.join(unpart_dir, 'data')
unpart_hpd_dir   = os.path.join(unpart_dir, 'hpd')
unpart_ss_dir    = os.path.join(unpart_dir, 'ss')
unpart_rb_dir    = os.path.join(unpart_dir, 'rb')

bycodon_dir       = os.path.join(dest_dir, 'bycodon')
bycodon_data_dir  = os.path.join(bycodon_dir, 'data')
bycodon_hpd_dir   = os.path.join(bycodon_dir, 'hpd')
bycodon_ss_dir    = os.path.join(bycodon_dir, 'ss')
bycodon_rb_dir    = os.path.join(bycodon_dir, 'rb')
       
bygene_dir       = os.path.join(dest_dir, 'bygene')
bygene_data_dir  = os.path.join(bygene_dir, 'data')
bygene_hpd_dir   = os.path.join(bygene_dir, 'hpd')
bygene_ss_dir    = os.path.join(bygene_dir, 'ss')
bygene_rb_dir    = os.path.join(bygene_dir, 'rb')
       
byboth_dir       = os.path.join(dest_dir, 'byboth')
byboth_data_dir  = os.path.join(byboth_dir, 'data')
byboth_hpd_dir   = os.path.join(byboth_dir, 'hpd')
byboth_ss_dir    = os.path.join(byboth_dir, 'ss')
byboth_rb_dir    = os.path.join(byboth_dir, 'rb')

######################
# Create directories #
######################

if os.path.exists(dest_dir):
    sys.exit('destination directory (%s) exists; please rename, delete, or move it and try again' % dest_dir)
os.mkdir(dest_dir)

os.mkdir(unpart_dir     )
os.mkdir(unpart_data_dir)
os.mkdir(unpart_hpd_dir )
os.mkdir(unpart_ss_dir  )
os.mkdir(unpart_rb_dir  )

os.mkdir(bycodon_dir     )
os.mkdir(bycodon_data_dir)
os.mkdir(bycodon_hpd_dir )
os.mkdir(bycodon_ss_dir  )
os.mkdir(bycodon_rb_dir  )

os.mkdir(bygene_dir     )
os.mkdir(bygene_data_dir)
os.mkdir(bygene_hpd_dir )
os.mkdir(bygene_ss_dir  )
os.mkdir(bygene_rb_dir  )

os.mkdir(byboth_dir     )
os.mkdir(byboth_data_dir)
os.mkdir(byboth_hpd_dir )
os.mkdir(byboth_ss_dir  )
os.mkdir(byboth_rb_dir  )
       
####################
# Write data files #
####################

unpart_boundaries  = writeNexusFile(os.path.join(unpart_data_dir, 'unpart.nex'),  taxa, None, [nchar0], [unpartseqs])

bycodon_boundaries  = writeNexusFile(os.path.join(bycodon_data_dir,'bycodon.nex'), taxa, None, [nchar1,nchar2,nchar3], [codon1seqs,codon2seqs,codon3seqs])
codon1st_boundaries = writeNexusFile(os.path.join(bycodon_data_dir,'codon1st.nex'), taxa, None, [nchar1], [codon1seqs])
codon2nd_boundaries = writeNexusFile(os.path.join(bycodon_data_dir,'codon2nd.nex'), taxa, None, [nchar2], [codon2seqs])
codon3rd_boundaries = writeNexusFile(os.path.join(bycodon_data_dir,'codon3rd.nex'), taxa, None, [nchar3], [codon3seqs])

bygene_boundaries  = writeNexusFile(os.path.join(bygene_data_dir, 'bygene.nex'),  taxa, None, [ncharCOI,ncharCOII,ncharATPase6,ncharATPase8], [COIseqs,COIIseqs,ATPase6seqs,ATPase8seqs])
COI_boundaries     = writeNexusFile(os.path.join(bygene_data_dir, 'COI.nex'),     taxa, None, [ncharCOI],     [COIseqs])
COII_boundaries    = writeNexusFile(os.path.join(bygene_data_dir, 'COII.nex'),    taxa, None, [ncharCOII],    [COIIseqs])
ATPase6_boundaries = writeNexusFile(os.path.join(bygene_data_dir, 'ATPase6.nex'), taxa, None, [ncharATPase6], [ATPase6seqs])
ATPase8_boundaries = writeNexusFile(os.path.join(bygene_data_dir, 'ATPase8.nex'), taxa, None, [ncharATPase8], [ATPase8seqs])

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

########################
# Create slurm scripts #
########################

submitall = '#!/bin/bash\n\n'
submitall += 'cd ..\n\n'
slurm_script_template = open('slurm-template.txt','r').read()

unpart_hpd_slurm_filename = os.path.join(unpart_hpd_dir,'s.sh')
unpart_hpd_slurm_contents = re.sub('__JOBNAME__',    'none-hpd', slurm_script_template, re.M | re.S)
unpart_hpd_slurm_contents = re.sub('__EMAIL__',      email,      unpart_hpd_slurm_contents, re.M | re.S)
unpart_hpd_slurm_contents = re.sub('__USERID__',     userid,     unpart_hpd_slurm_contents, re.M | re.S)
unpart_hpd_slurm_contents = re.sub('__EXECUTABLE__', executable, unpart_hpd_slurm_contents, re.M | re.S)
submitall += 'cd %s; sbatch s.sh; cd ../../..\n' % unpart_hpd_dir
f = open(unpart_hpd_slurm_filename,'w')
f.write(unpart_hpd_slurm_contents)
f.close()

unpart_ss_slurm_filename = os.path.join(unpart_ss_dir,'s.sh')
unpart_ss_slurm_contents = re.sub('__JOBNAME__',     'none-ss',  slurm_script_template, re.M | re.S)
unpart_ss_slurm_contents = re.sub('__EMAIL__',       email,      unpart_ss_slurm_contents, re.M | re.S)
unpart_ss_slurm_contents = re.sub('__USERID__',      userid,     unpart_ss_slurm_contents, re.M | re.S)
unpart_ss_slurm_contents = re.sub('__EXECUTABLE__',  executable, unpart_ss_slurm_contents, re.M | re.S)
submitall += 'cd %s; sbatch s.sh; cd ../../..\n' % unpart_ss_dir
f = open(unpart_ss_slurm_filename,'w')
f.write(unpart_ss_slurm_contents)
f.close()

unpart_rb_slurm_filename = os.path.join(unpart_rb_dir,'s.sh')
unpart_rb_slurm_contents = re.sub('__JOBNAME__',     'none-rb',  slurm_script_template, re.M | re.S)
unpart_rb_slurm_contents = re.sub('__EMAIL__',       email,      unpart_rb_slurm_contents, re.M | re.S)
unpart_rb_slurm_contents = re.sub('__USERID__',      userid,     unpart_rb_slurm_contents, re.M | re.S)
unpart_rb_slurm_contents = re.sub('__EXECUTABLE__',  revbayes,   unpart_rb_slurm_contents, re.M | re.S)
submitall += 'cd %s; sbatch s.sh; cd ../../..\n' % unpart_rb_dir
f = open(unpart_rb_slurm_filename,'w')
f.write(unpart_rb_slurm_contents)
f.close()

bycodon_hpd_slurm_filename = os.path.join(bycodon_hpd_dir,'s.sh')
bycodon_hpd_slurm_contents = re.sub('__JOBNAME__',   'codon-hpd', slurm_script_template, re.M | re.S)
bycodon_hpd_slurm_contents = re.sub('__EMAIL__',     email,       bycodon_hpd_slurm_contents, re.M | re.S)
bycodon_hpd_slurm_contents = re.sub('__USERID__',    userid,      bycodon_hpd_slurm_contents, re.M | re.S)
bycodon_hpd_slurm_contents = re.sub('__EXECUTABLE__', executable, bycodon_hpd_slurm_contents, re.M | re.S)
submitall += 'cd %s; sbatch s.sh; cd ../../..\n' % bycodon_hpd_dir
f = open(bycodon_hpd_slurm_filename,'w')
f.write(bycodon_hpd_slurm_contents)
f.close()

bycodon_ss_slurm_filename = os.path.join(bycodon_ss_dir,'s.sh')
bycodon_ss_slurm_contents = re.sub('__JOBNAME__',    'codon-ss', slurm_script_template, re.M | re.S)
bycodon_ss_slurm_contents = re.sub('__EMAIL__',      email,      bycodon_ss_slurm_contents, re.M | re.S)
bycodon_ss_slurm_contents = re.sub('__USERID__',     userid,     bycodon_ss_slurm_contents, re.M | re.S)
bycodon_ss_slurm_contents = re.sub('__EXECUTABLE__', executable, bycodon_ss_slurm_contents, re.M | re.S)
submitall += 'cd %s; sbatch s.sh; cd ../../..\n' % bycodon_ss_dir
f = open(bycodon_ss_slurm_filename,'w')
f.write(bycodon_ss_slurm_contents)
f.close()

bycodon_rb_slurm_filename = os.path.join(bycodon_rb_dir,'s.sh')
bycodon_rb_slurm_contents = re.sub('__JOBNAME__',    'codon-rb', slurm_script_template, re.M | re.S)
bycodon_rb_slurm_contents = re.sub('__EMAIL__',      email,      bycodon_rb_slurm_contents, re.M | re.S)
bycodon_rb_slurm_contents = re.sub('__USERID__',     userid,     bycodon_rb_slurm_contents, re.M | re.S)
bycodon_rb_slurm_contents = re.sub('__EXECUTABLE__', revbayes, bycodon_rb_slurm_contents, re.M | re.S)
submitall += 'cd %s; sbatch s.sh; cd ../../..\n' % bycodon_rb_dir
f = open(bycodon_rb_slurm_filename,'w')
f.write(bycodon_rb_slurm_contents)
f.close()

bygene_hpd_slurm_filename = os.path.join(bygene_hpd_dir,'s.sh')
bygene_hpd_slurm_contents = re.sub('__JOBNAME__',    'gene-hpd', slurm_script_template, re.M | re.S)
bygene_hpd_slurm_contents = re.sub('__EMAIL__',      email,      bygene_hpd_slurm_contents, re.M | re.S)
bygene_hpd_slurm_contents = re.sub('__USERID__',     userid,     bygene_hpd_slurm_contents, re.M | re.S)
bygene_hpd_slurm_contents = re.sub('__EXECUTABLE__', executable, bygene_hpd_slurm_contents, re.M | re.S)
submitall += 'cd %s; sbatch s.sh; cd ../../..\n' % bygene_hpd_dir
f = open(bygene_hpd_slurm_filename,'w')
f.write(bygene_hpd_slurm_contents)
f.close()

bygene_ss_slurm_filename = os.path.join(bygene_ss_dir,'s.sh')
bygene_ss_slurm_contents = re.sub('__JOBNAME__',    'gene-ss',  slurm_script_template, re.M | re.S)
bygene_ss_slurm_contents = re.sub('__EMAIL__',      email,      bygene_ss_slurm_contents, re.M | re.S)
bygene_ss_slurm_contents = re.sub('__USERID__',     userid,     bygene_ss_slurm_contents, re.M | re.S)
bygene_ss_slurm_contents = re.sub('__EXECUTABLE__', executable, bygene_ss_slurm_contents, re.M | re.S)
submitall += 'cd %s; sbatch s.sh; cd ../../..\n' % bygene_ss_dir
f = open(bygene_ss_slurm_filename,'w')
f.write(bygene_ss_slurm_contents)
f.close()

bygene_rb_slurm_filename = os.path.join(bygene_rb_dir,'s.sh')
bygene_rb_slurm_contents = re.sub('__JOBNAME__',    'gene-rb',  slurm_script_template, re.M | re.S)
bygene_rb_slurm_contents = re.sub('__EMAIL__',      email,      bygene_rb_slurm_contents, re.M | re.S)
bygene_rb_slurm_contents = re.sub('__USERID__',     userid,     bygene_rb_slurm_contents, re.M | re.S)
bygene_rb_slurm_contents = re.sub('__EXECUTABLE__', revbayes, bygene_rb_slurm_contents, re.M | re.S)
submitall += 'cd %s; sbatch s.sh; cd ../../..\n' % bygene_rb_dir
f = open(bygene_rb_slurm_filename,'w')
f.write(bygene_rb_slurm_contents)
f.close()

byboth_hpd_slurm_filename = os.path.join(byboth_hpd_dir,'s.sh')
byboth_hpd_slurm_contents = re.sub('__JOBNAME__',    'both-hpd', slurm_script_template, re.M | re.S)
byboth_hpd_slurm_contents = re.sub('__EMAIL__',      email,     byboth_hpd_slurm_contents, re.M | re.S)
byboth_hpd_slurm_contents = re.sub('__USERID__',     userid,     byboth_hpd_slurm_contents, re.M | re.S)
byboth_hpd_slurm_contents = re.sub('__EXECUTABLE__', executable, byboth_hpd_slurm_contents, re.M | re.S)
submitall += 'cd %s; sbatch s.sh; cd ../../..\n' % byboth_hpd_dir
f = open(byboth_hpd_slurm_filename,'w')
f.write(byboth_hpd_slurm_contents)
f.close()

byboth_ss_slurm_filename = os.path.join(byboth_ss_dir,'s.sh')
byboth_ss_slurm_contents = re.sub('__JOBNAME__',    'both-ss',  slurm_script_template, re.M | re.S)
byboth_ss_slurm_contents = re.sub('__EMAIL__',      email,     byboth_ss_slurm_contents, re.M | re.S)
byboth_ss_slurm_contents = re.sub('__USERID__',     userid,     byboth_ss_slurm_contents, re.M | re.S)
byboth_ss_slurm_contents = re.sub('__EXECUTABLE__', executable, byboth_ss_slurm_contents, re.M | re.S)
submitall += 'cd %s; sbatch s.sh; cd ../../..\n' % byboth_ss_dir
f = open(byboth_ss_slurm_filename,'w')
f.write(byboth_ss_slurm_contents)
f.close()

byboth_rb_slurm_filename = os.path.join(byboth_rb_dir,'s.sh')
byboth_rb_slurm_contents = re.sub('__JOBNAME__',    'both-rb',  slurm_script_template, re.M | re.S)
byboth_rb_slurm_contents = re.sub('__EMAIL__',      email,     byboth_rb_slurm_contents, re.M | re.S)
byboth_rb_slurm_contents = re.sub('__USERID__',     userid,     byboth_rb_slurm_contents, re.M | re.S)
byboth_rb_slurm_contents = re.sub('__EXECUTABLE__', revbayes, byboth_rb_slurm_contents, re.M | re.S)
submitall += 'cd %s; sbatch s.sh; cd ../../..\n' % byboth_rb_dir
f = open(byboth_rb_slurm_filename,'w')
f.write(byboth_rb_slurm_contents)
f.close()

submit_all_filename = os.path.join(dest_dir,'submit-all.sh')
f = open(submit_all_filename,'w')
f.write(submitall)
f.close()

#############################
# Create hpdml.conf scripts #
#############################

unpart_hpd_conf_template = open('conf-unpart-hpd-template.txt','r').read()
unpart_hpd_conf_filename = os.path.join(unpart_hpd_dir,'hpdml.conf')
unpart_hpd_conf_contents = re.sub('__LAST_SITE__', str(unpart_boundaries[0]),  unpart_hpd_conf_template, re.M | re.S)
unpart_hpd_conf_contents = re.sub('__RNSEED__',    rnseed,                     unpart_hpd_conf_contents, re.M | re.S)
unpart_hpd_conf_contents = re.sub('__TREEFILE__',  tree_file_name,             unpart_hpd_conf_contents, re.M | re.S)
f = open(unpart_hpd_conf_filename,'w')
f.write(unpart_hpd_conf_contents)
f.close()

unpart_ss_conf_template = open('conf-unpart-ss-template.txt','r').read()
unpart_ss_conf_filename = os.path.join(unpart_ss_dir,'hpdml.conf')
unpart_ss_conf_contents = re.sub('__LAST_SITE__', str(unpart_boundaries[0]), unpart_ss_conf_template, re.M | re.S)
unpart_ss_conf_contents = re.sub('__RNSEED__',    rnseed,                    unpart_ss_conf_contents, re.M | re.S)
unpart_ss_conf_contents = re.sub('__TREEFILE__',  tree_file_name,            unpart_ss_conf_contents, re.M | re.S)
f = open(unpart_ss_conf_filename,'w')
f.write(unpart_ss_conf_contents)
f.close()

unpart_rev_template = open('rev-unpart-template.txt','r').read()
unpart_rev_filename = os.path.join(unpart_rb_dir,'ss.Rev')
unpart_rev_contents = re.sub('__RNSEED__',    rnseed,                    unpart_rev_template, re.M | re.S)
unpart_rev_contents = re.sub('__TREEFILE__',  tree_file_name,            unpart_rev_contents, re.M | re.S)
f = open(unpart_rev_filename,'w')
f.write(unpart_rev_contents)
f.close()

bycodon_hpd_conf_template = open('conf-bycodon-hpd-template.txt','r').read()
bycodon_hpd_conf_filename = os.path.join(bycodon_hpd_dir,'hpdml.conf')
bycodon_hpd_conf_contents = re.sub('__FIRST_SITE_1ST_CODON__', '1',                            bycodon_hpd_conf_template, re.M | re.S)
bycodon_hpd_conf_contents = re.sub('__LAST_SITE_1ST_CODON__',  str(bycodon_boundaries[0]),     bycodon_hpd_conf_contents, re.M | re.S)
bycodon_hpd_conf_contents = re.sub('__FIRST_SITE_2ND_CODON__', str(bycodon_boundaries[0] + 1), bycodon_hpd_conf_contents, re.M | re.S)
bycodon_hpd_conf_contents = re.sub('__LAST_SITE_2ND_CODON__',  str(bycodon_boundaries[1]),     bycodon_hpd_conf_contents, re.M | re.S)
bycodon_hpd_conf_contents = re.sub('__FIRST_SITE_3RD_CODON__', str(bycodon_boundaries[1] + 1), bycodon_hpd_conf_contents, re.M | re.S)
bycodon_hpd_conf_contents = re.sub('__LAST_SITE_3RD_CODON__',  str(bycodon_boundaries[2]),     bycodon_hpd_conf_contents, re.M | re.S)
bycodon_hpd_conf_contents = re.sub('__RNSEED__',               rnseed,                         bycodon_hpd_conf_contents, re.M | re.S)
bycodon_hpd_conf_contents = re.sub('__TREEFILE__',             tree_file_name,                 bycodon_hpd_conf_contents, re.M | re.S)
f = open(bycodon_hpd_conf_filename,'w')
f.write(bycodon_hpd_conf_contents)
f.close()

bycodon_ss_conf_template = open('conf-bycodon-ss-template.txt','r').read()
bycodon_ss_conf_filename = os.path.join(bycodon_ss_dir,'hpdml.conf')
bycodon_ss_conf_contents = re.sub('__FIRST_SITE_1ST_CODON__', '1',                            bycodon_ss_conf_template, re.M | re.S)
bycodon_ss_conf_contents = re.sub('__LAST_SITE_1ST_CODON__',  str(bycodon_boundaries[0]),     bycodon_ss_conf_contents, re.M | re.S)
bycodon_ss_conf_contents = re.sub('__FIRST_SITE_2ND_CODON__', str(bycodon_boundaries[0] + 1), bycodon_ss_conf_contents, re.M | re.S)
bycodon_ss_conf_contents = re.sub('__LAST_SITE_2ND_CODON__',  str(bycodon_boundaries[1]),     bycodon_ss_conf_contents, re.M | re.S)
bycodon_ss_conf_contents = re.sub('__FIRST_SITE_3RD_CODON__', str(bycodon_boundaries[1] + 1), bycodon_ss_conf_contents, re.M | re.S)
bycodon_ss_conf_contents = re.sub('__LAST_SITE_3RD_CODON__',  str(bycodon_boundaries[2]),     bycodon_ss_conf_contents, re.M | re.S)
bycodon_ss_conf_contents = re.sub('__RNSEED__',               rnseed,                         bycodon_ss_conf_contents, re.M | re.S)
bycodon_ss_conf_contents = re.sub('__TREEFILE__',             tree_file_name,                 bycodon_ss_conf_contents, re.M | re.S)
f = open(bycodon_ss_conf_filename,'w')
f.write(bycodon_ss_conf_contents)
f.close()

bycodon_rev_template = open('rev-bycodon-template.txt','r').read()
bycodon_rev_filename = os.path.join(bycodon_rb_dir,'ss.Rev')
bycodon_rev_contents = re.sub('__RNSEED__',               rnseed,                         bycodon_rev_template, re.M | re.S)
bycodon_rev_contents = re.sub('__TREEFILE__',             tree_file_name,                 bycodon_rev_contents, re.M | re.S)
f = open(bycodon_rev_filename,'w')
f.write(bycodon_rev_contents)
f.close()

bygene_hpd_conf_template = open('conf-bygene-hpd-template.txt','r').read()
bygene_hpd_conf_filename = os.path.join(bygene_hpd_dir,'hpdml.conf')
bygene_hpd_conf_contents = re.sub('__FIRST_SITE_COI__',     '1',                           bygene_hpd_conf_template, re.M | re.S)
bygene_hpd_conf_contents = re.sub('__LAST_SITE_COI__',      str(bygene_boundaries[0]),     bygene_hpd_conf_contents, re.M | re.S)
bygene_hpd_conf_contents = re.sub('__FIRST_SITE_COII__',    str(bygene_boundaries[0] + 1), bygene_hpd_conf_contents, re.M | re.S)
bygene_hpd_conf_contents = re.sub('__LAST_SITE_COII__',     str(bygene_boundaries[1]),     bygene_hpd_conf_contents, re.M | re.S)
bygene_hpd_conf_contents = re.sub('__FIRST_SITE_ATPASE6__', str(bygene_boundaries[1] + 1), bygene_hpd_conf_contents, re.M | re.S)
bygene_hpd_conf_contents = re.sub('__LAST_SITE_ATPASE6__',  str(bygene_boundaries[2]),     bygene_hpd_conf_contents, re.M | re.S)
bygene_hpd_conf_contents = re.sub('__FIRST_SITE_ATPASE8__', str(bygene_boundaries[2] + 1), bygene_hpd_conf_contents, re.M | re.S)
bygene_hpd_conf_contents = re.sub('__LAST_SITE_ATPASE8__',  str(bygene_boundaries[3]),     bygene_hpd_conf_contents, re.M | re.S)
bygene_hpd_conf_contents = re.sub('__RNSEED__',             rnseed,                        bygene_hpd_conf_contents, re.M | re.S)
bygene_hpd_conf_contents = re.sub('__TREEFILE__',           tree_file_name,                 bygene_hpd_conf_contents, re.M | re.S)
f = open(bygene_hpd_conf_filename,'w')
f.write(bygene_hpd_conf_contents)
f.close()

bygene_ss_conf_template = open('conf-bygene-ss-template.txt','r').read()
bygene_ss_conf_filename = os.path.join(bygene_ss_dir,'hpdml.conf')
bygene_ss_conf_contents = re.sub('__FIRST_SITE_COI__',     '1',                           bygene_ss_conf_template, re.M | re.S)
bygene_ss_conf_contents = re.sub('__LAST_SITE_COI__',      str(bygene_boundaries[0]),     bygene_ss_conf_contents, re.M | re.S)
bygene_ss_conf_contents = re.sub('__FIRST_SITE_COII__',    str(bygene_boundaries[0] + 1), bygene_ss_conf_contents, re.M | re.S)
bygene_ss_conf_contents = re.sub('__LAST_SITE_COII__',     str(bygene_boundaries[1]),     bygene_ss_conf_contents, re.M | re.S)
bygene_ss_conf_contents = re.sub('__FIRST_SITE_ATPASE6__', str(bygene_boundaries[1] + 1), bygene_ss_conf_contents, re.M | re.S)
bygene_ss_conf_contents = re.sub('__LAST_SITE_ATPASE6__',  str(bygene_boundaries[2]),     bygene_ss_conf_contents, re.M | re.S)
bygene_ss_conf_contents = re.sub('__FIRST_SITE_ATPASE8__', str(bygene_boundaries[2] + 1), bygene_ss_conf_contents, re.M | re.S)
bygene_ss_conf_contents = re.sub('__LAST_SITE_ATPASE8__',  str(bygene_boundaries[3]),     bygene_ss_conf_contents, re.M | re.S)
bygene_ss_conf_contents = re.sub('__RNSEED__',             rnseed,                        bygene_ss_conf_contents, re.M | re.S)
bygene_ss_conf_contents = re.sub('__TREEFILE__',           tree_file_name,                bygene_ss_conf_contents, re.M | re.S)
f = open(bygene_ss_conf_filename,'w')
f.write(bygene_ss_conf_contents)
f.close()

bygene_rev_template = open('rev-bygene-template.txt','r').read()
bygene_rev_filename = os.path.join(bygene_rb_dir,'ss.Rev')
bygene_rev_contents = re.sub('__RNSEED__',             rnseed,                        bygene_rev_template, re.M | re.S)
bygene_rev_contents = re.sub('__TREEFILE__',           tree_file_name,                bygene_rev_contents, re.M | re.S)
f = open(bygene_rev_filename,'w')
f.write(bygene_rev_contents)
f.close()

byboth_hpd_conf_template = open('conf-byboth-hpd-template.txt','r').read()
byboth_hpd_conf_filename = os.path.join(byboth_hpd_dir,'hpdml.conf')
byboth_hpd_conf_contents = re.sub('__FIRST_SITE_COI1__',     '1',                            byboth_hpd_conf_template, re.M | re.S)
byboth_hpd_conf_contents = re.sub('__LAST_SITE_COI`__',      str(byboth_boundaries[0]),      byboth_hpd_conf_contents, re.M | re.S)
byboth_hpd_conf_contents = re.sub('__FIRST_SITE_COI2__',     str(byboth_boundaries[0] + 1),  byboth_hpd_conf_template, re.M | re.S)
byboth_hpd_conf_contents = re.sub('__LAST_SITE_COI2__',      str(byboth_boundaries[1]),      byboth_hpd_conf_contents, re.M | re.S)
byboth_hpd_conf_contents = re.sub('__FIRST_SITE_COI3__',     str(byboth_boundaries[1] + 1),  byboth_hpd_conf_template, re.M | re.S)
byboth_hpd_conf_contents = re.sub('__LAST_SITE_COI3__',      str(byboth_boundaries[2]),      byboth_hpd_conf_contents, re.M | re.S)
byboth_hpd_conf_contents = re.sub('__FIRST_SITE_COII1__',    str(byboth_boundaries[2] + 1),  byboth_hpd_conf_contents, re.M | re.S)
byboth_hpd_conf_contents = re.sub('__LAST_SITE_COII1__',     str(byboth_boundaries[3]),      byboth_hpd_conf_contents, re.M | re.S)
byboth_hpd_conf_contents = re.sub('__FIRST_SITE_COII2__',    str(byboth_boundaries[3] + 1),  byboth_hpd_conf_contents, re.M | re.S)
byboth_hpd_conf_contents = re.sub('__LAST_SITE_COII2__',     str(byboth_boundaries[4]),      byboth_hpd_conf_contents, re.M | re.S)
byboth_hpd_conf_contents = re.sub('__FIRST_SITE_COII3__',    str(byboth_boundaries[4] + 1),  byboth_hpd_conf_contents, re.M | re.S)
byboth_hpd_conf_contents = re.sub('__LAST_SITE_COII3__',     str(byboth_boundaries[5]),      byboth_hpd_conf_contents, re.M | re.S)
byboth_hpd_conf_contents = re.sub('__FIRST_SITE_ATPASE61__', str(byboth_boundaries[5] + 1),  byboth_hpd_conf_contents, re.M | re.S)
byboth_hpd_conf_contents = re.sub('__LAST_SITE_ATPASE61__',  str(byboth_boundaries[6]),      byboth_hpd_conf_contents, re.M | re.S)
byboth_hpd_conf_contents = re.sub('__FIRST_SITE_ATPASE62__', str(byboth_boundaries[6] + 1),  byboth_hpd_conf_contents, re.M | re.S)
byboth_hpd_conf_contents = re.sub('__LAST_SITE_ATPASE62__',  str(byboth_boundaries[7]),      byboth_hpd_conf_contents, re.M | re.S)
byboth_hpd_conf_contents = re.sub('__FIRST_SITE_ATPASE63__', str(byboth_boundaries[7] + 1),  byboth_hpd_conf_contents, re.M | re.S)
byboth_hpd_conf_contents = re.sub('__LAST_SITE_ATPASE63__',  str(byboth_boundaries[8]),      byboth_hpd_conf_contents, re.M | re.S)
byboth_hpd_conf_contents = re.sub('__FIRST_SITE_ATPASE81__', str(byboth_boundaries[8] + 1),  byboth_hpd_conf_contents, re.M | re.S)
byboth_hpd_conf_contents = re.sub('__LAST_SITE_ATPASE81__',  str(byboth_boundaries[9]),      byboth_hpd_conf_contents, re.M | re.S)
byboth_hpd_conf_contents = re.sub('__FIRST_SITE_ATPASE82__', str(byboth_boundaries[9] + 1),  byboth_hpd_conf_contents, re.M | re.S)
byboth_hpd_conf_contents = re.sub('__LAST_SITE_ATPASE82__',  str(byboth_boundaries[10]),     byboth_hpd_conf_contents, re.M | re.S)
byboth_hpd_conf_contents = re.sub('__FIRST_SITE_ATPASE83__', str(byboth_boundaries[10] + 1), byboth_hpd_conf_contents, re.M | re.S)
byboth_hpd_conf_contents = re.sub('__LAST_SITE_ATPASE83__',  str(byboth_boundaries[11]),     byboth_hpd_conf_contents, re.M | re.S)
byboth_hpd_conf_contents = re.sub('__RNSEED__',              rnseed,                         byboth_hpd_conf_contents, re.M | re.S)
byboth_hpd_conf_contents = re.sub('__TREEFILE__',            tree_file_name,                 byboth_hpd_conf_contents, re.M | re.S)
f = open(byboth_hpd_conf_filename,'w')
f.write(byboth_hpd_conf_contents)
f.close()

byboth_ss_conf_template = open('conf-byboth-ss-template.txt','r').read()
byboth_ss_conf_filename = os.path.join(byboth_ss_dir,'hpdml.conf')
byboth_ss_conf_contents = re.sub('__FIRST_SITE_COI1__',     '1',                            byboth_ss_conf_template, re.M | re.S)
byboth_ss_conf_contents = re.sub('__LAST_SITE_COI`__',      str(byboth_boundaries[0]),      byboth_ss_conf_contents, re.M | re.S)
byboth_ss_conf_contents = re.sub('__FIRST_SITE_COI2__',     str(byboth_boundaries[0] + 1),  byboth_ss_conf_template, re.M | re.S)
byboth_ss_conf_contents = re.sub('__LAST_SITE_COI2__',      str(byboth_boundaries[1]),      byboth_ss_conf_contents, re.M | re.S)
byboth_ss_conf_contents = re.sub('__FIRST_SITE_COI3__',     str(byboth_boundaries[1] + 1),  byboth_ss_conf_template, re.M | re.S)
byboth_ss_conf_contents = re.sub('__LAST_SITE_COI3__',      str(byboth_boundaries[2]),      byboth_ss_conf_contents, re.M | re.S)
byboth_ss_conf_contents = re.sub('__FIRST_SITE_COII1__',    str(byboth_boundaries[2] + 1),  byboth_ss_conf_contents, re.M | re.S)
byboth_ss_conf_contents = re.sub('__LAST_SITE_COII1__',     str(byboth_boundaries[3]),      byboth_ss_conf_contents, re.M | re.S)
byboth_ss_conf_contents = re.sub('__FIRST_SITE_COII2__',    str(byboth_boundaries[3] + 1),  byboth_ss_conf_contents, re.M | re.S)
byboth_ss_conf_contents = re.sub('__LAST_SITE_COII2__',     str(byboth_boundaries[4]),      byboth_ss_conf_contents, re.M | re.S)
byboth_ss_conf_contents = re.sub('__FIRST_SITE_COII3__',    str(byboth_boundaries[4] + 1),  byboth_ss_conf_contents, re.M | re.S)
byboth_ss_conf_contents = re.sub('__LAST_SITE_COII3__',     str(byboth_boundaries[5]),      byboth_ss_conf_contents, re.M | re.S)
byboth_ss_conf_contents = re.sub('__FIRST_SITE_ATPASE61__', str(byboth_boundaries[5] + 1),  byboth_ss_conf_contents, re.M | re.S)
byboth_ss_conf_contents = re.sub('__LAST_SITE_ATPASE61__',  str(byboth_boundaries[6]),      byboth_ss_conf_contents, re.M | re.S)
byboth_ss_conf_contents = re.sub('__FIRST_SITE_ATPASE62__', str(byboth_boundaries[6] + 1),  byboth_ss_conf_contents, re.M | re.S)
byboth_ss_conf_contents = re.sub('__LAST_SITE_ATPASE62__',  str(byboth_boundaries[7]),      byboth_ss_conf_contents, re.M | re.S)
byboth_ss_conf_contents = re.sub('__FIRST_SITE_ATPASE63__', str(byboth_boundaries[7] + 1),  byboth_ss_conf_contents, re.M | re.S)
byboth_ss_conf_contents = re.sub('__LAST_SITE_ATPASE63__',  str(byboth_boundaries[8]),      byboth_ss_conf_contents, re.M | re.S)
byboth_ss_conf_contents = re.sub('__FIRST_SITE_ATPASE81__', str(byboth_boundaries[8] + 1),  byboth_ss_conf_contents, re.M | re.S)
byboth_ss_conf_contents = re.sub('__LAST_SITE_ATPASE81__',  str(byboth_boundaries[9]),      byboth_ss_conf_contents, re.M | re.S)
byboth_ss_conf_contents = re.sub('__FIRST_SITE_ATPASE82__', str(byboth_boundaries[9] + 1),  byboth_ss_conf_contents, re.M | re.S)
byboth_ss_conf_contents = re.sub('__LAST_SITE_ATPASE82__',  str(byboth_boundaries[10]),     byboth_ss_conf_contents, re.M | re.S)
byboth_ss_conf_contents = re.sub('__FIRST_SITE_ATPASE83__', str(byboth_boundaries[10] + 1), byboth_ss_conf_contents, re.M | re.S)
byboth_ss_conf_contents = re.sub('__LAST_SITE_ATPASE83__',  str(byboth_boundaries[11]),     byboth_ss_conf_contents, re.M | re.S)
byboth_ss_conf_contents = re.sub('__RNSEED__',              rnseed,                         byboth_ss_conf_contents, re.M | re.S)
byboth_ss_conf_contents = re.sub('__TREEFILE__',            tree_file_name,                 byboth_ss_conf_contents, re.M | re.S)
f = open(byboth_ss_conf_filename,'w')
f.write(byboth_ss_conf_contents)
f.close()

byboth_rev_template = open('rev-byboth-template.txt','r').read()
byboth_rev_filename = os.path.join(byboth_rb_dir,'ss.Rev')
byboth_rev_contents = re.sub('__RNSEED__',              rnseed,                         byboth_rev_template, re.M | re.S)
byboth_rev_contents = re.sub('__TREEFILE__',            tree_file_name,                 byboth_rev_contents, re.M | re.S)
f = open(byboth_rev_filename,'w')
f.write(byboth_rev_contents)
f.close()

#################################
# Copy the gtrg-31taxa.tre file #
#################################
gtrg_31taxa_tre_contents = open('gtrg-31taxa.tre', 'r').read()
gtrg_31taxa_tre_filename = os.path.join(dest_dir,'gtrg-31taxa.tre')
f = open(gtrg_31taxa_tre_filename,'w')
f.write(gtrg_31taxa_tre_contents)
f.close()

#################################
# Copy the gtrg-32taxa.tre file #
#################################
gtrg_32taxa_tre_contents = open('gtrg-32taxa.tre', 'r').read()
gtrg_32taxa_tre_filename = os.path.join(dest_dir,'gtrg-32taxa.tre')
f = open(gtrg_32taxa_tre_filename,'w')
f.write(gtrg_32taxa_tre_contents)
f.close()


