# coding: utf8
import re 
import requests
import urllib.request, urllib.error, urllib.parse

from Bio import Entrez
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA
from Bio.Data.CodonTable import TranslationError
from bs4 import BeautifulSoup
from requests.exceptions import HTTPError

from ccd import app
from ccd.custom_exceptions import (IdNotFoundError, SequenceNotFoundError, 
                                   NotAnORF, IGiveUpError)

HEADERS = {'User agent': 'CrystallizationConstructDesigner_{}'.format(app.config['EMAIL'])}

class PicrMapper(object):
    def __init__(self):       
        raise NotImplementedError('PICR has been obsoleted by EBI. Use UniprotMapper instead')

class UniprotFetcher(object):
    '''
    Object to retrieve entries from Uniprot
    '''
    def __init__(self):
        super(UniprotFetcher, self).__init__()
        self.uniprot_address = 'https://www.uniprot.org/uniprot/'

    def fetch(self, protein_id, format_):
        '''
        Retrieves the uniprot entry corresponding to uniprot_id.
        Can fetch the entry in 'xml' or 'html' or 'sequence' format.
        Returns a Beautifulsoup of the webpage ('xml', 'html') or the protein
        sequence ('sequence')
        '''
        if format_ not in ['xml', 'html', 'sequence']:
            msg = '{} is not a supported format. Choose either "xml" or "html" or "sequence"'
            raise ValueError(msg)
        if format_ == 'sequence':   #"sequence" is easier to remember,
            format_ = 'fasta'       #but Uniprot uses the word "fasta"
        url = '{}{}.{}'.format(self.uniprot_address, protein_id, format_)
        try:
            response = RestLauncher().run_query(url, headers={})
        except HTTPError as e:
            raise HTTPError (f'Error in fetching the Uniprot entry {protein_id} from {url}') from e
        if format_ == 'xml': 
            return BeautifulSoup(response.text, 'xml')
        elif format_ == 'html':
            return BeautifulSoup(response.text, 'html.parser')
        elif format_ == 'fasta':
            # first line is sequence name
            return ''.join(response.text.split('\n')[1:])  

class UniprotParser(object):
    
    def __init__(self, xml_soup, html_soup):
        super(UniprotParser, self).__init__()
        self.xml_soup = xml_soup
        self.html_soup = html_soup    
    
    def get_scientific_name(self):
        return str(self.xml_soup.find_all('name', {'type':'scientific'})[0].string)
    
    def get_common_name(self): 
        try:
            return str(self.xml_soup.find_all('name', {'type':'common'})[0].string)
        except IndexError:  # no common name in entry
            return ''
         
    def get_protein_names(self):
        names = self.xml_soup.find_all('recommendedName')
        if names:  # well annotated entry
            for n in names:
                try:
                    short_name = str(n.find('shortName').string)  # the ways of BeautifulSoup are inscrutable
                except AttributeError:  # no short name present in entry
                    short_name = None 
                full_name = str(n.find('fullName').string)
                return full_name, short_name
        else:  # not very well annotated entry
            names = self.xml_soup.find_all('fullName')
            for n in names:
                full_name = str(n.next)  # BeautifulSoup is a bit weird sometimes...
            return full_name, None
    
    def get_kingdom(self):
        taxonomy = self.xml_soup.find_all('taxon')
        if taxonomy:
            for branch in taxonomy:
                if str(branch.next) in ['Archaea', 'Eukaryota', 'Bacteria']:
                    return str(branch.next)
        else:
            raise IGiveUpError #Kingdom is necessary for dna_finder to work.
    
    def get_isoform_descriptions(self):
        '''
        Input: a uniprot entry, html format, parsed by beautifulsoup
        Output: dict in format {isoform:description}
        where isoform is the uniprot isoform identifier
        and description is the description of the splicing variants the isoform has with
        respect to the canonical entry.
        Returns None if no isoform identifiers/descriptions are found
        '''
        entries = self.html_soup.find_all('div', {'class':'sequence-isoform'})
        parsed = {}
        for entry in entries:
            isoform_descr = []
            ''' the isoform identifier is present in the html as (identifier: XXXXXX-X)
            However, if only one isoform is present, no identifier is present,
            and re.findall is an empty list. In that case, we return None and let
            the caller figure it out'''
            try:
                isoform_name = re.findall(r'(identifier: .{6}-.)', entry.text)[0].split(': ')[1]
            except IndexError:
                return None
            '''the description of each isoform is as follows:
             first line is 'The sequence of this isoform differs from the canonical sequence as follows:'
            afterwards, missing exons are reported as e.g. u'1-138: Missing.'
            alternative exons are reported as e.g. u'1-8: MVVVTGRE u'\\u2192 u'MMASSYHE'
            the end of description contains the word ".« Hide".
            '''
            start = entry.text.find('The sequence of this isoform')
            stop = entry.text.find('.« Hide')
            description = entry.text[start:stop]
            # \xab and \xa0 are used instead of newlines - we replace and split
            cleaned = description.replace('\xab', '\n').replace('\xa0','\n').split('\n')
            # the canonical isoform is described as in standard text, and will result in cleaned == [''] above
            standard_txt = 'This isoform has been chosen by Uniprot as the "canonical" sequence.\n'
            if not cleaned[0]: #this is isoform 1, canonical
                parsed[isoform_name] = [standard_txt]
                continue
            #splitting the description for parsing
            txt = iter([t for t in cleaned if t])
            #after splitting, alternative exons will produce three separate elements in cleaned
            #e.g. ['1-8: MVVVTGRE, u'\\u2192, 'MMASSYHE']
            #here we rejoin these three elements in one single line
            while True:
                try:
                    t = txt.__next__()
                    #The sequence of this isoform differs from the canonical sequence as follows:
                    # or #1-#2: Missing
                    if t.startswith('The sequence') or 'Missing' in t:
                        isoform_descr.append(f'{t}\n')
                    elif t[0].isdigit(): #alternative exon
                        isoform_descr.append(f'{t} {txt.__next__()} {txt.__next__()}\n')
                    elif 'Hide' in t: #end of description
                        raise StopIteration
                except StopIteration:
                    parsed[isoform_name]=isoform_descr
                    break
        return parsed
    
    def get_crossreferences(self):
        '''
        takes a uniprot entry and returns all the crossrefrences to useful DNA
        databases in the form:
        {<database>:[ref1,ref2,...refn]
        for database in ['CCDS', 'EMBL', 'RefSeq'] 
        Reference XML schema:
        <...>
        <dbReference type="EMBL" id="BC015330">
            <property type="protein sequence ID" value="AAH15330.1"/>
            <property type="molecule type" value="mRNA"/>
        </dbReference>
        </...>
        '''
        dbreferences = {}
        # grabbing references from uniprot
        #crossreferences to CCDS and EMBL
        for database in ['CCDS', 'EMBL']:
            refs = self.xml_soup.find_all('dbReference', {'type': database})
            #remove references to genomic DNA - we don't need them, and they take ages to download
            if database == 'EMBL':
                for ref in list(refs):
                    is_eukaryote = self.get_kingdom() == 'Eukaryota'
                    for content in ref.contents:
                        try:
                            is_genomic_dna = 'Genomic_DNA' in content['value'] 
                            if is_genomic_dna and is_eukaryote:
                                refs.remove(ref)
                                break
                        except TypeError:
                            pass
            refs = [i['id'] for i in refs]
            dbreferences[database] = refs
        #crossreferences to RefSeq
        refs = []
        for ref in self.xml_soup.find_all('dbReference', {'type': 'RefSeq'}):
            a = ref.property.attrs['value']
            refs.append(a) 
        dbreferences['RefSeq'] = refs
        return dbreferences
    
class GeneFetcher(object):
    
    def __init__(self):
        super(GeneFetcher, self).__init__()
        self.ncbi_address = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
        
    def fetch(self, entrez_id, format_='xml'):
        if format_ not in ['xml', 'gb']:
            msg = 'Please specify a valid return format. Accepted: "xml" (default) or "gb"'
            raise ValueError(msg)
        req = 'efetch.fcgi?db=gene&id={}&retmode={}'.format(entrez_id, format_)
        launcher = UrllibLauncher()
        entry = launcher.run_query(self.ncbi_address + req, params={},
                                         headers={'tool':'CrystallizationConstructDesigner',
                                                  'email':app.config['EMAIL']})
        if format_ == 'xml':
            return BeautifulSoup(entry, 'xml')
        elif format_ == 'gb':
            raise NotImplementedError #we'll need this in the future

class GeneParser(object):
    
    def __init__(self):
        super(GeneParser, self).__init__()
    
    def get_crossreferences(self, gene_xml_soup):
        ccds_entries = list(set(
            [c for c in gene_xml_soup.find_all(text=re.compile(
                                                      r'CCDS[0-9]+[.]?[0-9]*'))]))
        return ccds_entries

class GenBankFetcher(object):
    '''
    Given a list of Genbank identifier, returns their xml entries as an iterable
    '''
    def __init__(self):
        super(GenBankFetcher, self).__init__()
    
    def check_size(self, id_):
        Entrez.email = app.config['EMAIL']
        handle = Entrez.esummary(db='nuccore', id=id_)
        try:
            summary = Entrez.read(handle)
        except RuntimeError: #some error in the esummary utility
            return -1
        filesize = summary[0]['Length'] #Biopython, a mistery wrapped in an enigma...
        return int(filesize)
        
    def fetch(self, id_list):
        Entrez.email = app.config['EMAIL']
        id_list = [id_ for id_ in id_list] #otherwise Efetch will panic because these are bytes in python 3
        try:
            handle = Entrez.efetch(db='nuccore', id=id_list, rettype="gb", retmode="xml")
        #defensive - if all the elements of the list of ids are invalid, HTTPError 404
        #is raised. If at least one Id in the list is valid, the invalid entries are ignored 
        except urllib.error.HTTPError: 
            return [] 
        soup = BeautifulSoup(handle, 'xml')
        entries = soup.find_all('GBSeq')
        return entries

class GenBankParser(object):
    
    def __init__(self):
        super(GenBankParser, self).__init__()
    
    def parse(self, xml_soup):
        '''
        Input: Beautifulsoup(xml) of Genbank entry
        Output:
        entry_id, DNA sequence, Protein sequence
        Raises NotAnORF if DNA does not contain a coding sequence (CDS) 
        '''
        is_mrna = bool(xml_soup.find_all('GBSeq_moltype')[0].text.strip() == 'mRNA')
        is_dna = bool(xml_soup.find_all('GBSeq_moltype')[0].text.strip() == 'cDNA')
        if not (is_dna or is_mrna):
            raise NotAnORF
        id_ = xml_soup.find_all('GBSeq_locus')[0].text.strip()
        features = xml_soup.find_all('GBFeature')
        if not features:  # entirely unannotated entry - ultra rare
            raise SequenceNotFoundError
        found = False
        for f in features: 
            if f.GBFeature_key.text.strip() == 'CDS':
                loc = f.find_all('GBFeature_location')[0].text
                loc = loc
                #sometimes format is start...pos2,pos3..stop);
                #if pos2 != pos3, or more intervals are indicated, we ignore this entry   
                if loc.startswith('join('):
                    loc = loc[5:-1].split('..')
                    positions = []
                    for pos in loc:
                        positions = positions + pos.split(',') #pos2,pos2
                    positions = list(set(positions))
                    try:
                        assert len(positions) == 3
                        start, stop = positions[0], positions[2]
                    except (AssertionError, ValueError):
                        raise NotAnORF 
                #mostly feature location is simply encoded as start..stop;
                else:
                    start, stop = loc.split('..')
                if '<' in start or '>' in stop:  # start or stop codon not known
                    raise NotAnORF
                start, stop = int(start), int(stop)
                found = True
                break
        if not found:  # not sure this ever happens
            raise SequenceNotFoundError
        dna_seq = xml_soup.GBSeq_sequence.text.strip().upper()
        orf = dna_seq[start-1:stop] 
        try: 
            assert orf.startswith('ATG') 
            assert orf.endswith(('TAA', 'TGA', 'TAG')) 
            assert len(orf)%3 == 0 
        except AssertionError: 
            raise NotAnORF 
        cds = Seq(orf, IUPACUnambiguousDNA()) 
        return id_, cds, cds.translate(cds=True)
    
    def parse_non_eukaryotes(self, xml_soup):
        id_ = str(xml_soup.find_all('GBSeq_primary-accession')[0].text).strip()
        features = xml_soup.find_all('GBFeature')
        if not features:  # entirely unannotated entry - ultra rare
            raise SequenceNotFoundError
        found = False
        #unlike eukaryotes, there are usually multiple CDS per entry, and they
        #might be on complementary strands...
        coding_sequences = []
        for f in features: 
            if f.GBFeature_key.text.strip() == 'CDS':
                try:
                    dna_seq = Seq(xml_soup.GBSeq_sequence.text.strip().upper(),
                                  IUPACUnambiguousDNA()) #don't move out of loop...
                except AttributeError: #entry does not actually have a normal sequence (e.g. HOPD_ECOLX)
                    raise SequenceNotFoundError 
                loc = f.find_all('GBFeature_location')[0].text
                start, stop = loc.split('..')
                if '<' in start or '>' in stop:  # start or stop codon not known
                    continue
                try:
                    if not 'complement('.upper() in start.upper(): #cds on sense strand
                        start, stop = int(start), int(stop)
                    elif 'complement('.upper() in start.upper(): #cds on other strand
                        #complement([start]..[stop])
                        start = int(start.split('(')[-1])
                        stop = int(stop.replace(')', '')) 
                        #reverse complement dna and remap
                        dna_seq = dna_seq.reverse_complement()
                        temp = start
                        start = len(dna_seq) - stop +1
                        stop = len(dna_seq) - temp +1
                except ValueError: #some other abstruse way of indicating starts and stops
                    continue 
                orf = dna_seq[start-1:stop]
                try: 
                    protein_seq = orf.translate(table=11, cds=True) #note that we use bacterial codon table
                    coding_sequences.append((id_, orf, protein_seq))
                    found = True
                except TranslationError: 
                    continue #not a good CDS
        if not found:  # not sure this ever happens
            raise NotAnORF
        return coding_sequences

class CCDSFetcher(object):
    
    def __init__(self):
        super(CCDSFetcher, self).__init__()
        
    def fetch(self, id_):
        '''
        Fetches the page corresponding to a ccds entry. Only HTML, because 
        there is no obvious API
        '''
        server = 'https://www.ncbi.nlm.nih.gov/projects/CCDS/'
        req = 'CcdsBrowse.cgi?REQUEST=CCDS&DATA={id_}&ORGANISM=0&BUILDS=CURRENTBUILDS'.format(id_=id_)
        launcher = UrllibLauncher()
        page = launcher.run_query(server + req, params={}, headers={})
        return BeautifulSoup(str(page), 'html.parser')     

class CCDSParser(object):
    
    def __init__(self):
        super(CCDSParser, self).__init__()
        
    def parse(self, html_soup): 
        # title: Report for CCDS[id].[version] (current version)
        #.[version] is optional.
        #" (current version)" might not be present
        titlematcher = re.compile(r'Report for CCDS[0-9]*(?:\.[0-9]*)(?:\ \(current version\))?')
        id_ = html_soup.find_all(string=titlematcher)[0] #find() does not take kwargs
        idmatcher = r'CCDS[0-9]*(?:.[0-9]*)?'
        id_ = re.search(idmatcher, id_).group(0)
        nucleotides = html_soup.find_all('span', {'id':re.compile('n[0-9]+')})
        aminoacids = html_soup.find_all('span', {'id':re.compile('p[0-9]+')})
        dna_seq = ''.join([nt.text for nt in nucleotides])
        aa_seq = ''.join([aa.text for aa in aminoacids])
        return id_, dna_seq, aa_seq