import os
import platform
import getpass
import string
import random
import logging

from subprocess import Popen, PIPE

from ccd import app
from ccd.custom_exceptions import (NoIsoformsPresentError, NoNeedToAlignError)

class PdbCrossref(object):
    
    def __init__(self, parent_entry_id: str, pdb_code: str, start: int, 
                 stop: int, chains: list):
        self.parent_entry_id = parent_entry_id
        self.pdb_code = pdb_code
        self.start = start
        self.stop = stop
        if not isinstance(chains, list):
            raise ValueError('Please pass a list of chains, even for a single chain')
        self.chains = chains
        
    def __repr__(self):
        return f'{self.pdb_code}_{"".join(self.chains)}: {self.start}-{self.stop}'
    
    def __lt__(self, other):
        #compare in order: pdb_code, start, stop (i.e. same start, shorter first), chains (should never change)) 
        return self.__key() < other.__key()
    
    def __gt__(self, other):
        return self.__key() > other.__key()
    
    def __eq__(self, other):
        #equal if everything is equal except parent entry. 
        #The same pdb can be referenced by different parent entries (different uniprot versions)
        return isinstance(self, type(other)) and self.__key() == other.__key()
    
    def __key(self):
        return self.pdb_code, self.start, self.stop, self.chains
    
    def __hash__(self):
        return hash(self.__key())
                                   
class AlignedProtein(object):
    
    def __init__(self, id_, protein_seq):
        super(AlignedProtein, self).__init__()
        self.id_ = id_
        self.protein_seq = protein_seq
        
    def __repr__(self):
#         return ':'.join((self.id_, self.protein_seq[:10] + ' ... ' + self.protein_seq[-10:]))
        return '>{}\n{}\n'.format(self.id_, self.protein_seq)
    
    def serialize(self):
        return[self.id_, self.protein_seq]
        
class HomologSequence(object):
    
    def __init__(self, ensembl_id_, protein_seq, species, gene_name='', 
                 homology_type=None):
        self.ensembl_id_ = ensembl_id_
        self.protein_seq = protein_seq
        self.species = species
        self.homology_type = homology_type 
        if gene_name:
            self.gene_name = gene_name
        else:
            self.gene_name = ''
        
    def __repr__(self):
        try:
            if self.uniprot_id_ and self.gene_name:
                return '{}_{}_({})'.format(self.get_short_species(),
                                            self.gene_name, 
                                            self.uniprot_id_)
            elif self.gene_name:
                return '{}_{}'.format(self.get_short_species(),
                                      self.gene_name)
            elif self.uniprot_id_:
                return '{}_({})'.format(self.get_short_species(),
                                        self.uniprot_id_)
            else:
                return '{}_{}'.format(self.get_short_species(),
                                      self.ensembl_id_)
            
        except AttributeError:
            return 'Protein {} from {}'.format(self.ensembl_id_,
                                               self.species)
    
    def __str__(self):
        if self.gene_name:
            return '{}_{}'.format(self.get_short_species(), self.gene_name)
        else:
            return self.__repr__()
    
    def get_short_species(self):
        try:
            genus = self.species.split('_')[0][0].upper()
            species = self.species.split('_')[1]
            return '{}.{}'.format(genus, species)
        except IndexError:
            return ('Unknown')

class UniprotProtein(object):
    
    def __init__(self, id_):
        super(UniprotProtein, self).__init__()
        self.id_ = id_
        self.uniprot_html_soup = None
        self.uniprot_xml_soup = None
        self.gene_xml_soup = None
        self.isoforms = []
        self.db_crossrefs = []
        self.cached_cds = []
        self.source_organism = ''
        self.source_organism_common = ''
        self.short_name = ''
        self.full_name = ''
        
    def __repr__(self):
        matched = [i for i in self.isoforms if i.matched_with]
        r = '''
        Protein: {}
        Isoforms: {}
        Matched_isoforms: {}
        '''.format(self.id_, len(self.isoforms), len(matched))
        return r
    
    def serialize_protein(self):
        try:
            id_ = self.uniprot_id
        except AttributeError:
            id_ = self.id_
        return {
                'protein_id': id_,
                'protein_isoforms': str(len(self.isoforms)),
                'protein_isoforms_list': [isoform.serialize_isoform() for isoform in self.isoforms],
                'protein_source_organism': self.source_organism,
                'protein_source_organism_common': self.source_organism_common,
            }
    
    def isoforms_are_present(self):
        if self.isoforms:
            return 1
        else:
            return 0
    
    def align_isoforms(self):
        if not self.isoforms_are_present():
            msg = 'No protein sequences to align. Please use the run() method before calling align_isoforms()'
            raise NoIsoformsPresentError
        elif len(self.isoforms) == 1:
            msg = 'Only one isoform is present, no need to align'
            raise NoNeedToAlignError(msg)
        
        aligner = MuscleAligner(self.isoforms)
        muscle_output = aligner.align()
        return muscle_output
    
    def fasta_to_list(self, fasta):
        lines = fasta.splitlines() 
        aligned_seqs = []
        name = ''
        seq = ''
        for line in lines:
            if line.startswith('>'):
                if seq != '':
                    aligned_seqs.append([name[1:], seq])
                    seq = ''
                    name = line
                else:
                    name = line
            else:
                seq = seq + line
        aligned_seqs.append([name[1:], seq])
        return aligned_seqs
    
    def get_fasta(self):
        if not self.isoforms_are_present():
            msg = 'No protein sequences to align. Please use the run() method before calling align_isoforms()'
            raise NoIsoformsPresentError(msg)
        seq = []
        for i in self.isoforms:
            seq.append('>{}\n{}\n'.format(i.id_, i.protein_seq))
        return ''.join(seq)
    
    def get_url(self):
        return 'https://www.uniprot.org/uniprot/{}'.format(self.id_)
    
class UniprotIsoform(object):
    
    def __init__(self, id_, protein_seq):
        super(UniprotIsoform, self).__init__()
        self.id_ = id_
        self.protein_seq = protein_seq
        self.dna_seq = None
        self.codingsequences = None
        self.matched_with = []
        self.mismatches = []
        self.description = ''
        
    def serialize_isoform(self):
        return {
            'isoform_id': self.id_,
            'isoform_matches': str(len(self.matched_with)),
            'isoform_matched_with': [cds.serialize_cds() for cds in self.matched_with],
            'isoform_mismatches': str(len(self.mismatches)),
            'isoform_mismatched_with': self.mismatches,
            'isoform_short_dna': str(self.get_short_dna_sequence()),
            'isoform_short_protein': str(self.get_short_protein_sequence()),
            'isoform_protein_length': str(len(self.protein_seq)),
            'isoform_description': self.description, #already unicode
            'isoform_dna': str(self.dna_seq)
        }
    
    def __repr__(self):
        r = '''
        Protein {id_}:
           DNA sequence = {dna_seq};
           Protein_sequence = {protein_seq}
           Description = {description}'''.format(id_ = self.id_,
                               dna_seq = self.get_short_dna_sequence(),
                               protein_seq = self.get_short_protein_sequence(),
                               description = self.description)
        return r
    
    def get_short_protein_sequence(self):
#         l = ' Length: {}'.format(len(self.protein_seq))
        return '{}...{}'.format(self.protein_seq[:10],
                                  self.protein_seq[-10:])
    
    def get_short_dna_sequence(self):
        if self.dna_seq:
            return '{}...{}'.format(self.dna_seq[:10],
                                      self.dna_seq[-10:])
        else: 
            return 'Undefined'
        
    def get_length(self):
        return len(self.protein_seq)
    
    def get_url(self):
        return 'https://www.uniprot.org/uniprot/{}#{}'.format(self.id_.split('-')[0],
                                                             self.id_)
        
class CodingSequence(object):
    
    def __init__(self, id_, database, dna_seq, protein_seq):
        self.id_ = id_
        self.database = database
        self.dna_seq = dna_seq
        self.protein_seq = protein_seq
        self.url = self.get_url()
    
    def __repr__(self):
        r = '''
        Coding Sequence {id_}:
           DNA sequence = {dna_seq};
           DNA source: {source}
           Protein_sequence = {protein_seq}'''.format(id_ = self.id_,
                               source = self.get_source(),
                               dna_seq = self.get_short_dna_sequence(),
                               protein_seq = self.get_short_protein_sequence())
                               
        return r
    
    def __hash__(self):
        try:
            return(hash(self.id_) ^ hash(str(self.dna_seq)))
        except AttributeError: #no dna_seq yet
            return(hash(self.id_))
    
    def __eq__(self, other):
        try:
            return self.__hash__() == other.__hash__()
        except TypeError:
            print(f'self is {self}')
            print(f'other is {other}')
        
    def serialize_cds(self):
        return {
            'cds_id': self.id_,
            'cds_database': self.database,
            'cds_dna': str(self.dna_seq),
            'cds_protein': str(self.protein_seq),
            'cds_url': self.get_url()
        }
    
    def get_short_protein_sequence(self):
#         l = ' Length: {}'.format(len(self.protein_seq))
        return (self.protein_seq[:10] + '...' + self.protein_seq[-10:]).upper() #+ l
    
    def get_short_dna_sequence(self):
        if self.dna_seq:
            return self.dna_seq[:10] + '...' + self.dna_seq[-10:]
        else: 
            return 'Undefined'
    
    def get_source(self):
        return '{id_} ({database})'.format(id_=self.id_, 
                                           database=self.database)
        
    def get_url(self):
        if self.database in ['EMBL', 'RefSeq']:
            return 'https://www.ncbi.nlm.nih.gov/nuccore/{}'.format(self.id_)
        elif self.database == 'CCDS':
            base = 'https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&DATA={}'
            return base.format(self.id_)
        
class CandidateORF(CodingSequence):
    
    def __init__(self, id_, database, parent_dna_seq, start, stop):
        super(CandidateORF, self).__init__(id_, database, parent_dna_seq[start:stop+1], '')
        self.parent_dna_seq = parent_dna_seq
        self.start = start
        self.stop = stop
    
    def serialize_cds(self):
        return {
            'cds_id': self.id_.decode('utf-8'),
            'cds_database': self.database.decode('utf-8'),
            'cds_dna': str(self.dna_seq).decode('utf-8'),
            'cds_protein': str(self.protein_seq).decode('utf-8'),
            'cds_url': self.get_url().decode('utf-8'),
            'cds_start': str(self.start).decode('utf-8'),
            'cds_stop' : str(self.stop).decode('utf-8')
            }
        
class MuscleAligner(object):
    
    '''
    Initializes with a list of sequence-like objects with the following attributes:
             protein_seq -> a valid protein sequence
             either:
                 id_ -> name of the sequence
             or __str__ if id_not present
             
    Will return a list of aligned sequence objects upon running the align()
    method. 
    '''
    
    def __init__(self, sequences):
        super(MuscleAligner, self).__init__()
        self.sequences = sequences
        self.muscle_exe = app.config['MUSCLE'] #TODO make this configurable on instantiating the class
        self.logger = logging.getLogger(__name__)
    
    def align(self):
        try:
            if not len(self.sequences) > 1:
                raise NoNeedToAlignError
        except TypeError: #somehow someway only a single homologsequence is passed
            raise NoNeedToAlignError
        fasta = self.sequences_to_fasta(self.sequences)
        muscle_out = self.run_muscle(fasta) #returns fasta format
        return self.fasta_to_sequences(muscle_out)
    
    def id_generator(self, size=6, chars=string.ascii_uppercase + string.digits):
        return ''.join(random.choice(chars) for _ in list(range(size)))
    
    def run_muscle(self, fasta):
        #NamedTemporaryFiles don't work on windows. Therefore this.
        if platform.system() == 'Windows':
            temp_path = 'C:\\Users\\{}\\AppData\\Local\\Temp\\'.format(getpass.getuser())
        elif platform.system() == 'Linux':
            temp_path = '/tmp'
        #to avoid race conditions when parallelizing we need random file names
        #the NamedTempFile module would not work under Windows, and I am using Windows
        #to code this.
        muscle_in = os.path.join(temp_path, self.id_generator())
        muscle_out = os.path.join(temp_path, self.id_generator())
        with open(muscle_in, 'w') as f:
            f.write(''.join(fasta))
            f.close()
        muscle_command = '{muscle} -in {input} -out {output}'.format(
                                                        muscle = self.muscle_exe,
                                                        input=muscle_in,
                                                        output=muscle_out)
        muscle_command = muscle_command.split()
        #python 2.7 has no context manager with Popen, unfortunately...
        p =  Popen(muscle_command, stdout=PIPE, stderr = PIPE)
        ali_length = len([line for line in fasta if line.startswith('>')])
        self.logger.info('Performing alignment of {} sequences'.format(ali_length))
        _, err = p.communicate() #muscle writes to stderr for whatever reason
        err = err.decode('utf-8')
        if platform.system() == 'Windows':
            p.kill()
        if 'ERROR' in err:
            #print('Errors: {}'.format(err))
            msg = 'Error in aligning'
            raise IOError(msg)
        if err:
            with open(muscle_out, 'r') as f:
                muscle_output = f.read()
        try:
            os.unlink(muscle_in)
            os.unlink(muscle_out)
        except OSError:
            pass
        return muscle_output
           
    def sequences_to_fasta(self, sequences):
        fasta = []
        for s in sequences:
            try:
                name = s.id_
            except AttributeError:
                name = str(s)
            fasta.append('>{}\n'.format(name))
            fasta.append('{}\n'.format(s.protein_seq))
        return fasta
        
    def fasta_to_sequences(self, fasta):
        proteins_list = []
        fasta_proteins = ''.join(fasta).split('>')
        #make list of (seq_name, sequence) tuples
        for fasta_protein in fasta_proteins:
            if fasta_protein.strip():
                name = fasta_protein.split('\n')[0]
                sequence = ''.join(fasta_protein.split('\n')[1:]).replace(' ', '')  # don't use first line which is the nam
                proteins_list.append(AlignedProtein(name, sequence))
        return proteins_list

if __name__ == '__main__':
    raise NotImplementedError('Nothing to run, move along')
    
# class AlignmentHelper(object):
#         '''
#         Helper class to launche multiple muscle alignments at once
#         Input: 
#         -> to_align: [ [ali1], [ali2], ... [aliN]] 
#         where aliN is a list of (name, sequence) tuples to be aligned with muscle
#         -> processes: how many threads should be spawned (def. 4)
#         -> muscle_exe: location of the muscle executable (default '/usr/bin/muscle')
#         Output:
#         a list of alignments [[ali1], [ali2], ... [aliN]]
#         where aliN is a list of (name, sequences) 
#         
#         '''
#         
#         def __init__(self, to_align, processes = 4, muscle_exe=None):
#             super(AlignmentHelper, self).__init__()
#             self.to_align = to_align
#             self.processes = processes
#             self.muscle_exe = muscle_exe #Windows debug monkeypatch
#         
#         def __call__(self, sequences):
#             return self.ali(sequences)
#         
#         def ali(self, sequences):
#             
#             aligner = MuscleAligner(sequences)
#             if self.muscle_exe: #Windows debug monkeypatch
#                 aligner.muscle_exe = self.muscle_exe 
#             return aligner.align()
#         
#         def parallel(self):
#             from multiprocessing import Pool
#             pool = Pool(processes = self.processes)
#             #Pool.map wants a callable that can be serialized.
#             #Class methods are not serializable, and thus don't fit the bill 
#             #defining the __call__method turns the class into a callable. 
#             #if __call__ returns the function we want to map to, we can pass the class
#             #to pool.map, and obtain parallelization
#             result = pool.map(self, self.to_align) # self == self.ali()
#             return [r for r in result if r]
    