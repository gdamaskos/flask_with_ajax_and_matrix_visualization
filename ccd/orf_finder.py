import os
import tempfile
import re
import io
import warnings

from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio import BiopythonParserWarning
from Bio.SeqUtils.ProtParam import ProteinAnalysis as PA

from ccd import vector_sequence
from ccd import app

#silencing biopython's warning about malformed locus
#TODO: fix malformed locus
warnings.simplefilter('ignore', BiopythonParserWarning)

def trim_frame(seq, mode='rev'):
    '''
    Makes sure that a seqobject's length is multiple of 3 to preserve frame
    seq: a biopython Seq or Seqrecord object or string
    mode: forward ("fwd") or reverse ("rev", default).
    If mode == fwd, the string will be trimmed from the end. If mode == rev, 
    from the beginning.
    e.g. ATGcccGATc -> ATGcccGAT if mode == fwd; 
         cATGcccGAT -> ATGcccGAT if mode == rev
    '''
    extra = len(seq) % 3
    if extra:
        if mode == 'rev':
            seq = seq[extra:]
        if mode == 'fwd':
            seq = seq[:-(extra)]
    return seq

def extend_to_Ct(seq, insert_end):
    '''
    Finds and returns the part of the ORF that extends C-terminal to the insert
    Input:
    seq = Bio.Seq.Seq object, typically a plasmid map
    insert_end: where does the cloned insert end
    Return:
    Bio.Seq.Seq object corresponding to the Ct extension of the ORF (translated)
    '''
    #todo: check for insert_end sanity
    #todo: circularize the plasmid in case we are at the end of the sequence
    Ct_ext = trim_frame(seq[insert_end:], mode='fwd')
    assert len(Ct_ext) % 3 == 0
    Ct_tag = Ct_ext.seq.translate().split('*')[0]
    return Ct_tag

def extend_to_Nt(seq, insert_start):
    Nt_ext = trim_frame(seq[:insert_start], mode='rev')
    Nt_ext = Nt_ext.seq.translate().split('*')[-1]
    Nt_tag = Nt_ext[Nt_ext.find('M'):]
    return Nt_tag[:insert_start]

def find_cleavage_sites(tagged_orf):
    '''
    Finds cleavage sites in a protein sequence
    Input: a Seq object having alphabet IUPAC.IUPACProtein 
    (e.g. from a Seqrecord translate() method, but not a SeqRecord object!, 
    because of incompatible APIs (thanks Biopython) ))
    Output: a dictionary in the form 
    {protease: [pos1, pos2,pos3]}
    where site_type is the kind of protease to be used for the cleavage
    and posx is the position in the ORF where the cleavage would occur
    multiple cleavage sites are possible, and they are enumerated Nt to Ct
    '''
    #TODO: put sites in a more sensible position for production code
    '''sites format: {protease_name: (recognition sequence, 
                                      cleavage_position_relative_to_rec._sequence)}
        and protease_name is each of ['3C', 'TEV', 'Thrombin', 'Senp2']
    '''
    sites = {'TEV':('ENLYFQ[GS]', 6),
             '3C':('VLFQGP', 4),
             'Thrombin':('LVPRGS', 4),
             'Senp2': ('[QE]QTGG', 4)}
    matches = {}
    for key, value in list(sites.items()):
        offset = sites[key][1]
        a = str(tagged_orf)
        x = re.findall(value[0], str(tagged_orf))
        matches[key] = [s.start() + offset  
                        for s in re.finditer(value[0], str(tagged_orf))]
           
    return matches

def perform_cleavage(tag, cleavage_sites, side='Nt'):
    '''
    Returns the remainder of the tag after digestion by all proteases in ccd.
    tags are returned as Seq objects
    
    Input: 
    tagged_orf = a Seq object having alphabet IUPAC.IUPACProtein
    cleavage_sites = the position of all cleavage sites (returned by find_cleavage_sites)
    
    Output:
    {protease_name: [orf0, orf1, ... orfx]}
    where orf1 is the digestion of tagged_orf at cleavage_sites[protease_name][0],
    etc.
    and protease_name is each of ['3C', 'TEV', 'Thrombin', 'Senp2'] 
    if no cleavage sites are present for a protease type, the uncleaved tag is returned
    within a list
    '''
    cleaved_orfs = {}
    for protease, match_positions in list(cleavage_sites.items()):
        cleaved_orfs[protease] = []
        if match_positions:
            if side == 'Nt':
                tags = [tag[m:] for m in match_positions]
            if side == 'Ct':
                tags = [tag[:m] for m in match_positions]
            #in case we have multiple tags, the shortest is the one resulting
            #from complete digestion
            cleaved_orfs[protease] = sorted(tags, key=lambda x: len(x))[0]
        else:
            cleaved_orfs[protease] = ''
    return cleaved_orfs

def build_cleaved_orfs(insert_translation, Nt_tag, Ct_tag):
    orfs = {}
    for protease in list(Nt_tag.keys()):
        if Nt_tag[protease] or Ct_tag[protease]: #we only return a cleaved orf if there is cleaveage site
            orfs[protease] = Nt_tag[protease] + insert_translation + Ct_tag[protease]
        else:
            orfs[protease] = []
    return orfs
    
def digest(insert_object, plasmid_map):
    '''
    Given a plasmid map and a insert (Bio.SeqFeature object), returns the tagged
    ORF and the possible ORFs resulting from protease digestion
    
    Input: 
    map = the plasmid map that has been generated and contains that insert
    insert_object = the insert itself, as a SeqFeature Biopython object 
                    (i.e. the same object that you sent to ligate_LIC to obtain 
                    the plasmid map)
    Output:
    tagged_orf = the insert with all tags (N- and C-terminal)
    cleaved_orfs = dict with format {protease_name: digested_orf}
        where protease_name is one of '3C', 'TEV', 'Thrombin' or 'Senp2'
        and digested_orf is a Seq object (IUPACProtein alphabet) containing the
        sequence resulting from digestion with that protease.
        If the protease does not cleave the ORF, digested_orf is an empty list
    '''
#     if run_tests:
#         test()
    #get some attributes from the insert (SeqRecord object)
    insert_dna_sequence = insert_object.extract(plasmid_map)
    insert_start =  insert_object.location.start
    insert_end = insert_object.location.end
    insert_translation = insert_dna_sequence.seq.translate()
    # Determining the N- and C-terminal tag if present
    Nt_tag = extend_to_Nt(plasmid_map, insert_start)
    Ct_tag = extend_to_Ct(plasmid_map, insert_end)
    # Determining what would happen to the tags if digested with protease
    Nt_tag_cleaved = perform_cleavage(Nt_tag, find_cleavage_sites(Nt_tag), side='Nt')
    Ct_tag_cleaved = perform_cleavage(Ct_tag, find_cleavage_sites(Ct_tag), side='Ct')

    tagged_orf = Nt_tag + insert_translation + Ct_tag
    cleaved_orfs = build_cleaved_orfs(insert_translation, Nt_tag_cleaved, 
                                      Ct_tag_cleaved)
    return tagged_orf, cleaved_orfs
    
def tag_cleave_pps(inserts_dna_seq, vname):
    vectors = vector_sequence.VectorLoader().make_vectors()
    zbuff = io.StringIO()
    tagged_cleaved = []
    
    
    for insert in inserts_dna_seq:
        cleaved_seqs = []
        proteases = []
        cleaved_MWs = []
        cleaved_pIs = []
        cleaved_epsilons = []
        
        vector_name = vectors[vname].name
        vector_short_name = vectors[vname].short_name
        #turn the insert into a SeqRecord & ligate
        annotated_insert = vector_sequence.seqrecord_from_string(insert)
        final_construct = vectors[vname].ligate_LIC(annotated_insert)
        #some annotation
        insert_name = insert[0]
        final_construct.description = '{}_in_{}'.format(insert_name, vector_name)
        short_construct_name = '{}_in_{}'.format(insert_name, vector_short_name)
        filename = '{}.gb'.format(short_construct_name)
        #preparing the map in .gb format
        with tempfile.NamedTemporaryFile(mode='r+') as temp1:
            with tempfile.TemporaryFile(mode='r+') as temp2:
                #output the map as genbank file via Biopython
                SeqIO.write(final_construct, temp2, 'gb')
                temp2.seek(0)
                old_header = temp2.readline()
                body = temp2.readlines()
            #extract and annotate the header
            new_header = vector_sequence.annotate_header(old_header, short_construct_name)
            #combine new header and body
            temp1.write(new_header)
            temp1.writelines(body)
            temp1.seek(0)
            #opening plasmid map
            plasmid_file = temp1.name
            with open(plasmid_file) as f:
                plasmid_map = SeqIO.read(f, 'gb', IUPAC.IUPACAmbiguousDNA())
#             #to George -> this is where the bug is... you need to add a check to
#             #see which one has the sequence of the insert (annotated_insert above)
#             insert_objects = [f for f in plasmid_map.features if f.type == 'CDS']
#             insert_object = insert_objects[1]
#             tagged_orf, cleaved_orfs = digest(insert_object, plasmid_map)
            # The following code will make is more clear when integrated! THe current solution is suboptimal.
            orfs_on_map = [f for f in plasmid_map.features if f.type == 'CDS']
            for orf in orfs_on_map:
                if orf.extract(plasmid_map).seq == annotated_insert.seq:
                    good_orf = orf
                    break
            tagged_orf, cleaved_orfs = digest(good_orf, plasmid_map)
            for p, s in list(cleaved_orfs.items()):
                if s:
                    cleaved_seqs.append(str(s))
                    proteases.append(p)
#             if cleaved_seq == '': # this is the case that we already have the antibiotic gene, so we pick the other element of the list
#                 insert_object = insert_objects[0]
#                 tagged_orf, cleaved_orfs = digest(insert_object, plasmid_map)
#                 for p,s in cleaved_orfs.items():
#                     if s:
#                         cleaved_seq = str(s)
#                         protease = p
                        
        tagged_seq = str(tagged_orf)
        #getting protparam values
        tagged_MW, tagged_pI, tagged_epsilon = get_protparam(tagged_seq)
        
        for cleaved_seq in cleaved_seqs:
            cleaved_MW, cleaved_pI, cleaved_epsilon = get_protparam(cleaved_seq)
            cleaved_MWs.append(cleaved_MW)
            cleaved_pIs.append(cleaved_pI)
            cleaved_epsilons.append(cleaved_epsilon)
            
        tagged_cleaved.append([insert_name, tagged_MW, tagged_pI, tagged_epsilon,
                               tagged_seq, cleaved_MWs, cleaved_pIs, 
                               cleaved_epsilons, cleaved_seqs, proteases])
    return tagged_cleaved


def get_protparam(sequence):
        
        s = PA(sequence)
        MW = round(s.molecular_weight()/1000, 1) #in kDa
        pI = round(s.isoelectric_point(), 1)
        aa_dict = s.count_amino_acids()
        # To calculate the epsilon, we use this formula from protparam (web.expasy.org/protparam)
        # Epsilon (Prot) = N(Tyr)*Ext(Tyr) + N(Trp)*Ext(Trp) + N(Cystine)*Ext(Cystine) / MW in Dalton
        epsilon = round((aa_dict['Y']*1490 + aa_dict['W']*5500 + aa_dict['C']*125)/(MW*1000), 2)          
        return MW, pI, epsilon
    
if __name__ == '__main__':
    #example usage of tag cleavage
    map_folder = os.path.abspath('../tests/unit/static')
    my_map = 'ex_1_20_in_NKI_1_1.gb' #3C cleavable, N-terminal tag
    my_map2 = 'ex_1_20_in_NKI_2_10.gb' #3C cleavable, N-terminal  and C-terminal tags
    map = os.path.join(map_folder, my_map)
    #opening map for a mock run
    with open(map) as f:
        plasmid_map = SeqIO.read(f, 'gb', IUPAC.IUPACAmbiguousDNA())
    '''
    We grab the insert from the map & various attributes
    Note: there are two orfs in the plasmid. One is the antibiotic resistance gene,
    the other is the ORF we want to cleave. Supposing the ORF of interest is 
    the second one:
    '''
    insert_object = [f for f in plasmid_map.features if f.type == 'CDS'][1]
    tagged_orf, cleaved_orfs = digest(insert_object, plasmid_map)
    
    print('The peptide with tags has the sequence: {}'.format(str(tagged_orf)))
    for p,s in cleaved_orfs.items():
        if s:
            print('After digestion with {protease}, '
                  'the peptide will have sequence: {sequence}'.format(protease = p,
                                                              sequence = s))
