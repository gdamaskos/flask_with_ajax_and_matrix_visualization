
import sys
import logging
import os
import urllib.request
import argparse

from subprocess import Popen, PIPE
from requests.exceptions import HTTPError

try:
    import ccd
except ImportError:
    sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from ccd import app

pdb_seqres_url = 'ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt'



class Entry(tuple):
    '''
    convenient way of ordering pdb_entries.
    t = tuple
    t[0] = entry id in format XXXX_A, where A indicates chain
    t[1] = fasta sequence from seqres
    pdb entries are considered identical if they have the same pdb_id and same
    sequence, regardless of chain, i.e. 
    Entry(1a2x_A, seq) == Entry(1a2x_B, seq) 
    
    >>> a = Entry(('1a2x_1','AAA'))
    >>> b = Entry(('1a2x_2','AAA'))
    >>> c = Entry(('1a2x_3','AAB'))
    >>> l = [a,b,c]
    >>> u = list(set(l))
    >>> l
    [('1a2x_1', 'AAA'), ('1a2x_2', 'AAA'), ('1a2x_3', 'AAB')]
    >>> u
    [('1a2x_3', 'AAB'), ('1a2x_1', 'AAA')]
    >>> (a.__hash__() == b.__hash__())
    True
    >>> (a.__hash__() == c.__hash__())
    False
    >>> a is b
    False
    >>> a == b
    True
    >>> a == c
    False
    Note that this "breaks" the in operator:
    >>>x = set([a,b,c])
    {a,c}
    >>>a in x
    True
    >>> b in x
    True # Set works on hash, and a and b hash equal, even though a is not b
    I don't care because as far as I am concerned, I need only one of either a or b,
    and I don't care which one, 
    '''
    
    def __init__(self, t):
        self.t = t
        
    def __eq__(self, other):
        return (self.t[0][:4], self.t[1]) == (other.t[0][:4], other.t[1])
     
    def __ne__(self, other):
        return not self.__eq__(self, other) 
        
    def __hash__(self):
        return hash(
            (self.t[0][:4], self.t[1])
            )
    
    def __len__(self):
        return len(self.t[1])

def read_fasta(fasta_file):
    seqs = []
    count = 1
    with open (fasta_file, 'r') as f:
        while True:
            print('Processing entry {}'.format(count), end='\r')
            fasta_header = f.readline()[1:] # >id_chain_number - remove '>'
            if not fasta_header:
                break
            if 'mol:protein' in fasta_header:
                seq = f.readline()
                '''some entries have XXXX_A and XXXX_a as separate chains. 
                makeblastdb with parse seq id would complain of duplicate entries
                because it uppercases everything. 
                we add the counter as a unique identifier''' 
                pdb_id = f'{fasta_header.split("mol:protein")[0].strip()}_{count}'
                seqs.append(Entry(
                    (pdb_id, seq.upper().rstrip())
                    ))
            count += 1
    print('\n')
    return seqs

def to_fasta(entries):
    out = []
    for entry in entries:
        out.append ('>{}\n{}\n'.format(entry[0], entry[1]))
    return ''.join(out)
                
def remove_redundancy(seqs):
    return list(set(seqs))

def _store(processed, destination):
    import pickle
    with open(destination, 'w+') as f:
        pickle.dump(processed, f)
        
def write_fasta(fasta, destination):
    with open(destination, 'w+') as f:
        f.write(fasta)
        
def retrieve(destination):
    import pickle
    with open(destination, 'r+') as f:
        return pickle.load(f)

def fetch_sequences_from_pdb(destination_file):
    
    req = urllib.request.Request(pdb_seqres_url)
    try:
        response = urllib.request.urlopen(req)
    except HTTPError as e:
        logging.warn( 'The server couldn\'t fulfill the request.')
        logging.warn( 'Error code: ', e.code)
    with open(destination_file, 'w') as f:
        f.write(response.read().decode('utf-8'))
        
def convert_database(fasta_in, db_out):
    cmd = f'makeblastdb -in {fasta_in} -out {db_out} -dbtype prot -parse_seqids'
    p =  Popen(cmd.split(), stderr=PIPE, stdout=PIPE)
    out, err = p.communicate()
    logging.info('Makeblastdb output:\n')
    for msg in out.decode('utf-8').split('\n'):
        if msg:
            logging.info(f'{msg}')
    if err:
        logging.warn(('Makeblastdb reported errors\n {}'.format(err)))

def process(destination_folder, update_database, keep_intermediate_files):
    
    PDB_FASTA = os.path.join(destination_folder, 'pdb_seqres.txt')
    PDB_TEMP_FASTA_CLEANED = os.path.join(destination_folder, 'pdb_processed.fasta')
    PDB_BLAST_DB = os.path.join(destination_folder, 'pdb_blastdb')
    
    logging.basicConfig(level=logging.INFO)
    if not os.path.isdir(destination_folder):
        os.makedirs(destination_folder, exist_ok=True)
    if update_database:
        logging.info(f'Fetching sequences from wwwPDB to {PDB_FASTA}')
        logging.info('This might take a while')
        fetch_sequences_from_pdb(PDB_FASTA)
        logging.info('Download completed')
    seqs = read_fasta(PDB_FASTA)
    logging.info('Cleaning non unique entries')
    unique = remove_redundancy(seqs)
    logging.info('Converting to Fasta')
    fasta_string =  to_fasta(unique)
    logging.info(f'Writing fasta file in {PDB_TEMP_FASTA_CLEANED}')
    write_fasta(fasta_string, PDB_TEMP_FASTA_CLEANED)
    logging.info('Converting database for usage with blast')
    convert_database(PDB_TEMP_FASTA_CLEANED, PDB_BLAST_DB)
    if not keep_intermediate_files:
        logging.info('Removing temporary files')
        try:
            os.remove(PDB_TEMP_FASTA_CLEANED)
            os.remove(PDB_FASTA)
        except OSError:
            pass
    logging.info('Done!')
    
def get_cmd_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--destination_folder',
                        help='The folder the database files will be stored in. It will be created if necessary. Default: ./ccd/static/blastdb',
                        )
    parser.add_argument('--do_not_update_pdb', 
                        help='Will not fetch the latest version of the PDB database from the web.',
                        action='store_false')
    parser.add_argument('--keep_intermediate_files', 
                        help='The script will clean intermediate files (e.g. the raw pdb database). [True/False] default = True. ',
                        action='store_true')
    parser = parser.parse_args() 
    return parser

def check_arguments(parser):
    update_database = True 
    keep_intermediate_files = False
    if parser.do_not_update_pdb is not None:
        update_database = parser.do_not_update_pdb
    if parser.keep_intermediate_files is not None:
        keep_intermediate_files = parser.keep_intermediate_files
    if parser.destination_folder is not None:
        destination_folder = os.path.abspath(parser.destination_folder)
    else:
        destination_folder = os.path.join(app.root_path, 'static/blast_db')
    return destination_folder, update_database, keep_intermediate_files 

if __name__ == '__main__':
    args = get_cmd_arguments()
    destination_folder, update_database, keep_intermediate_files  = check_arguments(args)
    seqres_path = os.path.join(destination_folder, 'pdb_seqres.txt')
    if not update_database and not os.path.isfile(seqres_path):
        import sys
        msg = f'I cannot find a raw fasta file at {seqres_path}. ' +\
              'Please allow update from the web, or point the script to an appropriate fasta file'
        sys.exit(msg)
    process(destination_folder, update_database, keep_intermediate_files)
    
    