import json
import os
import types
import tempfile
import zipfile
import datetime
import uuid
import io

from Bio.Seq import Seq
from io import BytesIO, StringIO

from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ccd import app

maps_folder =  os.path.join(app.static_folder, 'plasmid_maps')

class BaseVector(object):

    def __init__(self):
        super(BaseVector, self).__init__()
        with open(os.path.join(maps_folder, self.filename), 'r',
                  encoding='utf-8', errors='ignore') as f:
            #Bio.SeqIO stumbles a bit on utf-8 conversion
            #let's have python do it properly beforehand
            data = StringIO(''.join(f.readlines()))  
        self.seq = SeqIO.read(data, 'gb')

    def __repr__(self):
        if self.Ct_tag:
            t = 'Ct tag vector'
        else:
            t = 'Nt tag vector'
        if self.is_LIC:
            t += ', LIC'
        return '{}; {}'.format(self.short_name, t)
        
    def ligate_LIC(self, insert):
        '''
        Accepts an insert, i.e. your ORF (Seq or SeqRecord object); 
        returns the insert correctly ligated. 
        Note: do not add primer overhangs.
        Note: we don't check the insert at all (i.e. Starts, stops, etc.)
        We assume that the insert is correct.
        '''
        lig = self.seq[:self.cut_left] + self.add_left + insert + \
              self.add_right + self.seq[self.cut_right:]
        return lig
    
class Ct_tag_vector(BaseVector):

    def __init__(self):
        super(Ct_tag_vector, self).__init__()

    def ligate_LIC(self, insert):
        '''
        We need to override this in case we have a stop codon in the ORF
        '''
        insert_nostop = insert
        for stop in ['TAA', 'TGA', 'TAG']:
            if insert.seq.upper().endswith(stop):
                insert_nostop = insert[:-3]
                qualifiers = insert.features[0].qualifiers
                feature = SeqFeature(FeatureLocation(0, len(insert)), 
                         type="CDS", strand=1, qualifiers=qualifiers)
                insert_nostop.features.append(feature)
                break
        lig = self.seq[:self.cut_left] + self.add_left + insert_nostop + \
              self.add_right + self.seq[self.cut_right:]
        return lig


class VectorLoader(object):
    '''
    Helper class to instantiate a dictionary of {vector.short_name : vector_object}
    From a correctly formatted json file
    vectors are read from two configuration files (json)
    1 - base vectors from static/base_vectors.json
    2 - custom vectors from static/custom_vectors.json
    See vectors_help.txt for documentation on how to add custom vectors to ccd
    '''
    
    basepath = app.static_folder
    BASE_VECTOR_FILE = os.path.join(basepath, 'base_vectors.json')
    CUSTOM_VECTOR_FILE = os.path.join(basepath, 'custom_vectors.json')
    MAP_PATH = os.path.join(basepath, 'plasmid_maps')

    def __init__(self):
        super(VectorLoader, self).__init__()
        
    def load_vectors_from_file(self, json_file):
        try:
            with open(json_file, 'r') as f:
                return json.load(f)
        except json.decoder.JSONDecodeError as e: #Let's help the user out...
            msg = f'Error in decoding custom vector file {json_file}.\nCheck your syntax'
            raise type(e)(f'{str(e)}\n{msg}', e.doc, e.pos)
            
    def instantiate_vectors(self, vectors):
        c = {}
        for i in vectors:
            if not i['Ct_tag']:
            # this creates the vector object from the json dictionary element
                vector_object = type(i['class_name'], (BaseVector,), i)() # note the () - we directly instantiate each vector object
                    #this simply creates a dictionary entry for accessing the objects later
                c[vector_object.short_name] = vector_object
            else:
                vector_object = type(i['class_name'], (Ct_tag_vector,), i)()
                c[vector_object.short_name] = vector_object
        return c
        
    def make_vectors(self):
        base_vectors = self.load_vectors_from_file(self.BASE_VECTOR_FILE)
        container = {}
        container.update(self.instantiate_vectors(base_vectors))
        try: #load custom vectors if present
            custom_vectors = self.load_vectors_from_file(self.CUSTOM_VECTOR_FILE)
            if custom_vectors[0]: #empty file = [{}]
                container.update(self.instantiate_vectors(custom_vectors))
        except OSError:
            pass 
        return container
    
def vector_names():
    # make sure maps and vector are in this folder
    vectors = os.path.join(app.static_folder, 'base_vectors.json')
    if not os.path.isfile(vectors):
        raise OSError(f'The file containing plasmid data is not found at {vectors}')
    vectors = VectorLoader().make_vectors()
    vnames = [i.short_name for _, i in list(vectors.items())]
    vnames = sorted(vnames, #sort by vector type and by number, i.e. NKI 1_14 goes before NKI_2_1
                    key = lambda x: tuple(int(num) for num in x.split('_')[1:3]))
    return vnames


def seqrecord_from_string(insert):
    # insert[0]: name, insert[1]: sequence
    # Making a suitable SeqRecord from the annotated_sequence
    # 1 create a Seq object
    ins = Seq(insert[1], IUPAC.unambiguous_dna)
    # 2 turn the Seq into a SeqRecord
    annotated_sequence = SeqRecord(ins)
    # 3 create a feature for the annotated_sequence
    qualifiers = {'label' : insert[0]}
    feature = SeqFeature(FeatureLocation(0, len(annotated_sequence)), 
                         type="CDS", strand=1, qualifiers=qualifiers)
    # 4 append the feature to the SeqRecord object
    annotated_sequence.features.append(feature)
    return annotated_sequence

def annotate_header(old_header, construct_name):
    '''
    Annotates the Biopython-created header to:
    1 - indicate that the plasmid is circular and synthetic DNA
    2 - adds locus information
    3 - adds creation date
    
    Input: 
    old_header : the old header (whole line)
    insert_name: a short name for the final construct, e.g. "insert_in_vectorname"
    Output: properly annotated header
    
    Genbank format specification: 
    https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html
    '''
    header = old_header.split()
    #exactly 7 spaces are needed after 'LOCUS' for the header to 
    #be read correctly. Three now, four in the join below
    header[0] = header[0].replace ('LOCUS', 'LOCUS   ')
    #add construct description
    header[1] = construct_name
    #indicate circular DNA
    header[4] = header[4].replace('DNA', 'DNA\tcircular') 
    #indicate this is a synthetic construct
    header[5] = 'SYN'
    #Add correct creation date 
    date_now = datetime.datetime.now()
    date_formatted = date_now.strftime('%d-%b-%Y')
    header[6] = '{}\n'.format(date_formatted)
    return '    '.join(header)

def make_plasmid_maps(inserts_dna_seq, vname):
    vectors = VectorLoader().make_vectors()
    zbuff = io.BytesIO()
    zip_file = zipfile.ZipFile(zbuff, 'w')
    for insert in inserts_dna_seq:
        vector_name = vectors[vname].name
        vector_short_name = vectors[vname].short_name
        #turn the insert into a SeqRecord & ligate
        annotated_insert = seqrecord_from_string(insert)
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
                new_header = annotate_header(old_header, short_construct_name)
                #combine new header and body
                temp1.write(new_header)
                temp1.writelines(body)
                temp1.seek(0)
            zip_file.write(temp1.name, arcname = filename, 
                          compress_type = zipfile.ZIP_DEFLATED)
    zip_file.close()
    zbuff.seek(0)
    return zbuff

def poly2dna_seq(dna_orf, pname, polypeptides_starts, polypeptides_stops):
    inserts_dna_seq = []
    for index, _ in enumerate(polypeptides_starts):
        mapped_start = (int(polypeptides_starts[index]) - 1) * 3
        mapped_stop = int(polypeptides_stops[index]) * 3
        ppname = '{}_{}_{}'.format(pname, polypeptides_starts[index], polypeptides_stops[index])
        inserts_dna_seq.append([ppname, dna_orf[mapped_start: mapped_stop]])
    return inserts_dna_seq # format name, sequence
