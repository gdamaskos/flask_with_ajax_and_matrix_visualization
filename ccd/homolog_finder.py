import re
import requests
import urllib.request, urllib.error, urllib.parse
import sys
import logging
from io import BytesIO

from bs4 import BeautifulSoup

from ccd import app
from ccd import dna_finder_async as dba
from ccd.sequence_classes import MuscleAligner, HomologSequence
from ccd.db_utilities_async import UniprotMapper, UniprotFetcher
from ccd.custom_exceptions import NoAlignmentCached, IdNotFoundError, IGiveUpError

model_orgs = ["saccharomyces_cerevisiae",       #Baker's yeast
"caenorhabditis_elegans",                       #Roundworm  
"drosophila_melanogaster",                      #Fruit fly
"danio_rerio",                                  #Zebrafish
"oryzias_latipes",                              #Medaka fish
"xenopus_laevis",                               #African clawed frog
"mus_musculus",                                 #Mouse
"homo_sapiens",                                 #You
"gallus_gallus",                                #Chicken
"anolis_carolinensis",                          #Carolina anole lizard
"meleagris_gallopavo",                          #Turkey
"sarcophilus_harrisii"]                         #Tasmanian devil

default_flags = {'ortholog_one2one': True,
                 'ortholog_many2many': False,
                 'ortholog_one2many': True,
                 'within_species_paralog' : False,
                 'model' : True,
                 'other_paralog': False} #we ignore other paralogs anyway

class Evolutionator(object):
    
    def __init__(self):
        super(Evolutionator, self).__init__()
        
    def fill_data(self, id_, seq, species, flags=default_flags):
        correct_init = seq or id_
        if not correct_init:
            msg = 'Please provide either a protein sequence or a valid uniprot id'
            raise ValueError(msg)
        self.flags = flags
        self.seq = seq
        self.id_ = id_
        self.species = species
        self.configure_servers()
        self.headers = {}
        self.hits_by_flag = {}
        logging.basicConfig(stream=sys.stdout, level=logging.INFO)
        self.logger = logging.getLogger(__name__)
        
    def search(self):
        #case 1: user did not provide an id: we map using local blast to the best match
        if not self.id_: #we try to map by sequence because the user gave no id
            blast_results = self.run_local_blast(self.seq, 'swissprot',
                                                 'csv', fast=True)
            parsed = self.parse_blast_results(blast_results, 'csv')
            high_identity = self.filter_blast_results(parsed, 'evalue', 0, 1e-10)
            try:
                best_matches = high_identity[:5] #take the best 5
                to_match = [match.uniprot_id for match in best_matches]
            except IndexError: #all matches are lousy
                raise IdNotFoundError 
        else: # case 2: user provided an id
            if not isinstance(self.id_, list): #always the case for CCD, might not be the case if called from outside
                self.id_ = [self.id_] 
            to_match =  self.id_ #a list of one in this case 
            self.uniprot_id = self.id_[0]
        # let's map the uniprot id  or the 5 best matches to ENSEMBL. 
        # Not all ids map out of the box         # will raise IdNotFoundError
        mapper = dba.UniprotMapper()  
        mapping_to_ensembl = mapper.map_(to_match, 'ACC+ID', 'ENSEMBL_ID')
        # mapping to ensembl is a dictionary of uniprot_id : ensembl_id key : value pairs
        # in python 3.6 the dict is ordered by insertion
        # We choose here not to trust that and do some more jogging
        # to make sure that we always choose the highest scoring match
        for uniprot_id in to_match:
            try:
                self.ensembl_id = mapping_to_ensembl[uniprot_id]
                break #we found the ensembl_id of the best identit match
            except KeyError:
                continue
        else: #none of the best matches has an ensembl ID. we give up
            raise IdNotFoundError
        ali = self.get_ensembl_alignment(self.ensembl_id)
        if ali:
            self.logger.info("Found an alignment for {}".format(self.ensembl_id))
        sequences = self.parse_ensembl_alignment(ali)
        self.cached_sequences = self.rename_sequences(sequences)
        #counting how many hits we have for each flag and storing it in the main object
        counts = {}
        for flag in default_flags:
            counts[flag] = 0
        for s in sequences: 
            counts[s.homology_type] += 1
            if s.species in model_orgs:
                counts['model'] += 1
        self.hits_by_flag = counts
        #counting how many sequences are from model organisms
        return self.cached_sequences
    
    def get_organism(self, id_):
        a = dba.DatabaseCrawler(id_)
        try:
            a.search_uniprot()
        except IGiveUpError:
            msg = 'Cannot find source species for {}'.format(self.id_)
            raise IdNotFoundError(msg)
        return a.protein.source_organism.decode('utf-8').replace(' ', '_').lower()
    
    def configure_servers(self):
        self.ensembl_address = 'http://rest.ensembl.org'
        self.uniprot_id_map_address = 'https://www.uniprot.org/uploadlists/'
    
    def get_ensembl_alignment(self, id_):
        '''given a ENSEMBL ID, retrieves a precalculated alignment of orthologs''' 
        server = self.ensembl_address
        ext = "/homology/id/{}".format(id_)
        ali = requests.get(server+ext, headers={ "Content-Type" : "application/json"})#application/json"})
        if not ali.ok:
            ali.raise_for_status()
            sys.exit()
        self.logger.info('Retrieved alignment for {}'.format(id_))
        return ali.json()
    
    def run_local_blast(self, query, database, return_format, fast=True):
        blaster = dba.BlastRunner()
        return blaster.run_blast(query, database, return_format,
                                 fast=fast)
    
    def parse_blast_results(self, blast_result, result_format):
        if not result_format in ['xml', 'csv']:
            raise ValueError
        parser = dba.BlastParser()
        return parser.parse_csv(blast_result)
    
    def filter_blast_results(self, parsed_results, filter_by, 
                             min_value, max_value):
        parser = dba.BlastParser()
        filtered_results = parser.filter_results(parsed_results, 
                                                 filter_by, 
                                                 min_value, 
                                                 max_value)
        sorted_results = sorted(filtered_results, 
                                key= lambda x: getattr(x, filter_by))
        return sorted_results
#         id_ = ''
#         if id_:
#             return [i.accession.text for i in id_]
#         else:
#             msg = 'Cannot map the protein sequence to a valid Uniprot id'
#             raise IdNotFoundError(msg)
    
    def parse_ensembl_alignment(self, alignment):
        sequences = []
#         with open('/home/andrea/Desktop/abraxas.xml','w') as f:
#             json.dump(alignment, f)
        for hit in alignment['data'][0]['homologies']:
            species = hit['target']['species'] 
            protein_seq = hit['target']['align_seq'].replace('-', '')
            ensembl_id_ = hit['target']['id']
            homology_type = hit['type']
            sequences.append(HomologSequence(ensembl_id_, protein_seq, species,
                                             homology_type=homology_type))
        return sequences
    
    def align(self, to_align):
        aligner = MuscleAligner(to_align)
        muscle_output = aligner.align()
        return muscle_output
    
    def select_sequences_from_ali(self, flags, sequences):
        #select by flag
        to_align = []
        for flag in list(flags.keys()):
            if self.flags[flag]: 
                for seq in sequences:
                    if seq.homology_type == flag:
                        to_align.append(seq)
        to_align = list(set(to_align)) #There should be no duplicates, but let's be sure
        #select by mode (model organism only vs full alignment)
        if flags['model']:
            for seq in list(to_align):
                if seq.species not in model_orgs:
                    to_align.remove(seq)
        return flags, to_align   
    
    def rename_sequences(self, sequences):
        #mapping ensembl id to iuniprot Id_, and then uniprot_id_ to gene name
        uniprot_ids = self.uniprot_id_map('ENSEMBL_ID', 'ID', 
                                          ' '.join([i.ensembl_id_ 
                                                    for i in sequences]))
        gene_names = self.uniprot_id_map('ID', 'GENENAME', 
                                          ' '.join([i for i in list(uniprot_ids.values())]))
        #updating our sequence objects
        for seq in sequences:
            try:
                seq.uniprot_id_ = uniprot_ids[seq.ensembl_id_]
            except KeyError:
                seq.uniprot_id_ = ''
            try:
                seq.gene_name = gene_names[seq.uniprot_id_]
            except KeyError:
                seq.gene_name = ''
        return sequences
        
    def get_default_alignment(self, flags, sequences):
        #selecting which sequences to align
        _, to_align = self.select_sequences_from_ali(flags, sequences)
        #checking if the query is present in the alignment retrieved from the database 
        self.query = HomologSequence(self.id_, self.seq, self.species, gene_name='Query')
        #if an identical sequence is present, we rename it to query
        for s in list(to_align):
            if s.protein_seq == self.query.protein_seq:
                s.gene_name = 'Query'
                break
        #otherwise we add our sequence to the alignment queue
        else:
            to_align.append(self.query)
        if len(to_align) > 1:
            self.logger.info('Starting alignment')
            alignment = self.align(to_align)
            self.logger.info('Alignment completed')
            #hit_by_type -> how many hits are present under each flag. type: dict {flag : #hits}
        return (flags, self.serialize_alignment(alignment), self.hits_by_flag)
    
    def serialize_alignment(self, alignment):
        return [s.serialize() for s in alignment]
    
    def fasta_to_list(self, aligned):
        # algorithm to parse fasta multisequence
        lines = aligned.splitlines() 
    
        aligned_seqs = []
        name = ''
        seq = ''
        for line in lines:
            seq = seq.replace('-', ' ')
            if line.startswith('>'):
                if seq != '':
                    aligned_seqs.append([name[1:], seq])
                    seq = ''
                    name = line
                else:
                    name = line
            else:
                seq = seq + line
        seq = seq.replace('-', ' ')
        aligned_seqs.append([name[1:], seq])
        return aligned_seqs
    
    def uniprot_id_map(self, from_, to, ids):
        params = {
        'from':from_,
        'to':to,
        'format':'tab',
        'query': ids
        }
        data = urllib.parse.urlencode(params).encode('utf-8')
        request = urllib.request.Request(self.uniprot_id_map_address, data)
        #TODO add ccd email / url
        contact = app.config['EMAIL'] # Please set your email address here to help us debug in case of problems.
        request.add_header('User-Agent', 'Python %s' % contact)
        response = urllib.request.urlopen(request)
        page = response.read(200000).decode()
        accs = dict([tuple((i.split('\t'))) 
                      for i in page.split('\n')[1:]
                      if i])
        return accs
    
    def realign(self, flags, sequences):
        _, to_align = self.select_sequences_from_ali(flags, sequences)
        try:
            to_align.append(self.query)
        except AttributeError:
            msg = 'Realign can only be used after align()'
            raise NoAlignmentCached(msg)
        alignment = self.align(to_align)
        self.logger.info('Aligning with flags {}'.format(flags))
        return (flags, self.serialize_alignment(alignment))
        
    def get_fasta(self, alignment):
        # (alignment[0], alignment[1]) == (flags, [AlignedSequence1...AlignedSequencen]
        bytes_IO = BytesIO()
        #AlignedSequence objects __repr__ to '>id_\n\protein_seq\n (aka fasta format)'
        fasta = ''.join([str(a) for a in alignment[1]])
        bytes_IO.write(fasta.encode('utf-8')) 
        bytes_IO.seek(0)
        return bytes_IO 
        
if __name__ == '__main__':
    usage = '''
    Usage: 
        - can search ENSEMBL for a list of homolog sequences to a given protein
        - can align such proteins using muscle
    Initialize: 
        Evolutionator(uniprot_id, 
                      protein_seq, 
                      source_organism)
    Required:
        - a valid protein sequence OR
        - a valid uniprot id_ AND a protein sequence (to account for different isoforms)
        
    initialization via id only would be possible (and trivial), but would require making 
    implicit assumptions about what the user wants - we prefer to let the user choose
    Error codes:
        - raises IdNotFoundError if no mapping to ENSEMBL is possible
    '''
    print(usage)

