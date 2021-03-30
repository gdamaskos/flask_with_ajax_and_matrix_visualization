# coding: utf8
import os
import re 
import requests
import urllib.request, urllib.parse
import asyncio
import aiohttp
import glob

from subprocess import PIPE, Popen
from tempfile import NamedTemporaryFile

from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA, ExtendedIUPACProtein
from Bio.Data.CodonTable import TranslationError
from bs4 import BeautifulSoup
from requests.exceptions import HTTPError

from ccd import app
from ccd.sequence_classes import CodingSequence, PdbCrossref
from ccd.custom_exceptions import (IdNotFoundError, SequenceNotFoundError, 
                                   NotAnORF, IGiveUpError)
from collections import namedtuple

HEADERS = {'User agent': 'CrystallizationConstructDesigner_{}'.format(app.config['EMAIL'])}

class PicrMapper(object):
    def __init__(self):       
        raise NotImplementedError('PICR has been obsoleted by EBI. Use UniprotMapper instead')

class UniprotMapper(object):
    '''
    USAGE:
    mapper = UniprotMapper()
    mapped_id = mapper.map_(old_id, old_database, new_database)
    raises IdNotFounderror if not found
    '''
    
    def __init__(self):
        super(UniprotMapper, self).__init__()
        self.uniprot_address = 'https://www.uniprot.org/uploadlists/'
        #see https://www.uniprot.org/help/api_idmapping for more databases
        #here we use only protein- and dna-centric databases
        from_dbs = ['ACC+ID', 'ACC', 'ID', 'UPARC', 'NF50', 'NF90', 'NF100', 'GENENAME',
               'CRC64', 'EMBL_ID', 'EMBL', 'P_ENTREZGENEID', 'P_GI', 'PIR',
               'REFSEQ_NT_ID', 'P_REFSEQ_AC', 'UNIGENE_ID', 'PDB_ID', 'ENSEMBL_ID']
        to_dbs = list(from_dbs)
        to_dbs.remove('ACC+ID')
        self.allowed = {'from': from_dbs, 'to': to_dbs}
        
    def map_(self, query, from_, to):
        self.validate_input(from_, to) #will raise ValueError if invalid input
        if isinstance(query, list):
            q = ' '.join(query) #Uniprot wants id=id1 id2 id3 ... idn format
        else:
            q= query
        params = {'from': from_,
                  'to': to,
                  'format':'tab',
                  'query':q}
        launcher = UrllibLauncher()
        #get page
        page = launcher.run_query(self.uniprot_address, params)  #uniprot wants no headers
        #parse page
        if page == 'From\tTo\n':
            msg = '{} not found in database {}'.format(query, to)
            raise IdNotFoundError(msg)
        elif isinstance(query, str):
            gene_id = page.split('\t')[-1].strip()
        elif isinstance(query, list):
            gene_id = [(i.split()[0], i.split()[1]) #(uniprot_id, ensembl_id)
                       for i in page.splitlines()[1:] #discard first line
                       ]
            gene_id = dict(gene_id)
        else:
            raise ValueError
        return gene_id
        
    def validate_input(self, from_, to):
        if not from_ in self.allowed['from']:
            msg = f'Unknown value for from_ parameter: {from_}'
            msg += f'Allowed values: {self.allowed["from"]}'
            raise ValueError(msg)
        if not to in self.allowed['to']:
            msg = f'Unknown value for from_ parameter: {to}'
            msg += f'Allowed values: {self.allowed["to"]}'
            raise ValueError(msg)

class UrllibLauncher(object):
    
    def __init__(self):
        super(UrllibLauncher, self).__init__()
        
    def run_query(self, url, params=None, headers=None):
        if params is not None:
            data = urllib.parse.urlencode(params).encode('utf-8')
        request = urllib.request.Request(url, data)
        if headers is not None:
            for key, value in list(headers.items()):
                request.add_header(key, value)
        response = urllib.request.urlopen(request)
        page = response.read(200000).decode('utf-8')
        return page

class RestLauncher(object):
    '''
    Runs GET requests and returns the response
    Raises error on 4xx or 5xx status 
    '''
    def __init__(self):
        super(RestLauncher, self).__init__()
    
    def run_query(self, req, headers=None):
        if headers is None:
            r = requests.get(req, timeout=20)
        else:
            r = requests.get(req, headers=headers, timeout=15)
        if not r.ok:
            r.raise_for_status()  # raises the appropriate errors if 4xx or 5xx 
        return r   

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
        if format_ not in ['xml', 'html', 'sequence', 'fasta']:
            msg = '{} is not a supported format. Choose either "xml" or "html" or "sequence"/"fasta"'
            raise ValueError(msg)
        if format_ == 'sequence':   #"sequence" is easier to remember,
            format_ = 'fasta'       #but Uniprot uses the word "fasta"
        url = '{}{}.{}'.format(self.uniprot_address, protein_id, format_)
        try:
            response = RestLauncher().run_query(url)
        except HTTPError as e:
            raise IdNotFoundError(f'Error in fetching the Uniprot entry {protein_id} from {url}') from e
        if format_ == 'xml': 
            return BeautifulSoup(response.text, 'xml')
        elif format_ == 'html':
            return BeautifulSoup(response.text, 'html.parser')
        elif format_ == 'fasta':
            # first line is sequence name
            return ''.join(response.text.split('\n')[1:])  

class UrlFormatter(object):
    '''Create url(s) to fetch nucleic acid sequences from databases'''
    
    def __init__(self):
        super().__init__()
        self.ncbi_address = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
        self.CCDS_address = 'https://www.ncbi.nlm.nih.gov/projects/CCDS/'
        self.ncbi_esummary_address = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
        self.uniprot_address = 'https://www.uniprot.org/uniprot/'
        self.uniprot_batch_address = 'https://www.uniprot.org/uploadlists'
        self.email = app.config['EMAIL']
        
    def format(self, database: str, id_list, format_=None) -> list:
        '''Produce url(s) to retrieve data from databases
        Input: 
        id_list: one or more valid gene identifiers 
        database: one of ['EMBL', 'RefSeq', 'CCDS', 'Gene', 'Uniprot']
        format_: required id database == 'Uniprot', ignored otherwise.
                 allowed_formats = ['fasta','xml','html', 'batch_xml', 'batch_fasta']
                 batch retrieval uses uploadlists
        Output:
        a list of tuples (database, url).
        If database == CCDS: len(output) == len(id_list)
        Else: len(output) == 1 because EMBL, RefSeq and Gene allow batch retrieval
        '''
        ## some input checking
        if not isinstance(id_list, list):
            id_list = [id_list]
        allowed_databases = ['EMBL', 'RefSeq', 'CCDS', 'Gene', 'Uniprot']
        allowed_formats = ['fasta','xml','html', 'batch_xml', 'batch_fasta', 'batch_html', None] #uniprot only
        if database not in allowed_databases:
            raise ValueError(f'Supported database types: {allowed_databases}')
        elif format_ not in allowed_formats:
            raise ValueError(f'Supported formats: {allowed_formats}')
        elif database == 'Uniprot' and format_ is None:
            raise ValueError(f'Uniprot requires to specify a format: {allowed_formats}')
        urls = []
        #make the actual URLs
        if database in ['EMBL', 'RefSeq']:
            url = (f'{self.ncbi_address}efetch.fcgi?db=nuccore&id={",".join(id_list)}&retmode=xml&rettype=gb&email={self.email}&tool=ccd.rhpc.nki.nl')
            urls.append((database, url))
        elif database == 'CCDS':
            for id_ in id_list:
                url = f'{self.CCDS_address}CcdsBrowse.cgi?REQUEST=CCDS&DATA={id_}&ORGANISM=0&BUILDS=CURRENTBUILDS'
                urls.append((database, url))
        elif database == 'Gene':
            url = f'{self.ncbi_address}efetch.fcgi?db=Gene&id={",".join(id_list)}&retmode=xml&email={self.email}&tool=ccd.rhpc.nki.nl'
            urls.append((database, url))
        elif database == 'Uniprot':
            if 'batch' not in format_:
                for id_ in id_list:
                    url = f'{self.uniprot_address}{id_}.{format_}'
                    urls.append((database, url))
            else:
                f = format_.split('_')[1]
                url = f'{self.uniprot_batch_address}?query={" ".join(id_list)}&from=ACC+ID&to=ACC&format={f}'
                urls.append((database, url))
        return urls
    
    def format_for_summary(self, database: str, id_list) -> list:
        if not isinstance(id_list, list):
            id_list = [id_list]
        allowed = ['EMBL', 'RefSeq']
        if database not in allowed:
            raise ValueError(f'Supported database types: {allowed}')
        url = f'{self.ncbi_esummary_address}?db=nuccore&id={",".join(id_list)}'
        return [('Summary', url)]
        
class Entry_fetcher(object):
    '''Given database and url, fetch entries asynchronously'''
    
    
    #solution for retry on HTTP errors from 
    #http://allyouneedisbackend.com/blog/2017/09/15/how-backend-software-should-retry-on-failures/    
    def retry(self, *exceptions, retries=3, cooldown=1, verbose=True):
        import logging
        from functools import wraps
        log = logging.getLogger(__name__)
        """Decorate an async function to execute it a few times before giving up.
        Hopes that problem is resolved by another side shortly.
    
        Args:
            exceptions (Tuple[Exception]) : The exceptions expected during function execution
            retries (int): Number of retries of function execution.
            cooldown (int): Seconds to wait before retry.
            verbose (bool): Specifies if we should log about not successful attempts.
        """
    
        def wrap(func):
            @wraps(func)
            async def inner(*args, **kwargs):
                retries_count = 0
    
                while True:
                    try:
                        result = await func(*args, **kwargs)
                    except exceptions as err:
                        retries_count += 1
                        message = f"Exception during {func} execution. {retries_count} of {retries} retries attempted"
                        if retries_count > retries:
                            verbose and log.exception(message)
                            raise (
                                func.__qualname__, args, kwargs) from err
                        else:
                            verbose and log.warning(message)
    
                        if cooldown:
                            await asyncio.sleep(cooldown)
                    else:
                        return result
            return inner
        return wrap

    @retry(aiohttp.web.HTTPRequestTimeout, aiohttp.web.HTTPTooManyRequests, aiohttp.client_exceptions.ClientResponseError)
    async def _fetch(self, session, database : str, url : str) -> tuple:
        '''Fetch a webpage (coroutine)
        Inout: 
        session = aiohttp session
        database =  a string, that is simply returned
        url = url to fetch
        
        Output:
        (database, response text)
        '''
        async with session.get(url) as response:
            if response.status != 200:
                response.raise_for_status()
            return (database, await response.text())
    
    async def fetch_all(self, queries):
        '''Fetch DNA from NCBI or CCDS asynchronously
        Input:
        list of (database, url) tuples 
        database in ['EMBL', 'Refseq', 'CCDS', 'Gene']
        Return:
        list of (database, response.text) tuples. 
        If response.text of EMBL or Gene or RefSeq contains multiple entries, these are
        returned as separate (database, entry.text) tuples. 
        '''
        #fetching all results
        timeout = aiohttp.ClientTimeout(total=60)
        async with aiohttp.ClientSession(timeout=timeout) as session:
#             for query in queries:
#                 tasks.append(asyncio.ensure_future(self._fetch(session, query[0], query[1])))
            tasks = [asyncio.ensure_future(self._fetch(session, query[0], query[1]))
                     for query in queries]
            return await asyncio.gather(*tasks)
            
class EntrySplitter(object):
    '''
    Takes a list of (database, xml_blob) and returns a list of 
    (database, single_xml_entry) tuples
    '''
    
    def split(self, results):
        #refseq and EMBL entries are concatenated. Let's break them up
        entries = []
        for r in results:
            database = r[0]
            entry_blob = r[1]
            entries += (self._split_entries(database, entry_blob))
        return entries
    
    def _split_entries(self, database, xml_response):
        soup = BeautifulSoup(xml_response, 'xml')
        if database == 'CCDS':
            return [(database, xml_response)]
        elif database == 'Gene':
            return [(database, str(entry)) for entry in soup.find_all('Entrezgene')]
        elif database in ['EMBL', 'RefSeq']:
            return [(database, str(entry)) for entry in soup.find_all('GBSeq')]
        elif database == 'Summary':
            return [(database, str(entry)) for entry in soup.find_all('DocSum')]
        elif database == 'Uniprot':
            return [(database, str(entry)) for entry in soup.find_all('entry')]
        else:
            raise ValueError(f'Database {database} is not recognized')

class DnaParser(object):
    
    def parse(self, entries):
        parsed = []
        CCDS_Parser = CCDSParser()
        GenBank_Parser = GenBankParser()
        Gene_Parser = GeneParser()
        for entry in entries:
            database = entry[0]
            text = entry[1]
            if database in ['EMBL', 'RefSeq']:
                soup = BeautifulSoup(text, 'xml')
                try:
                    parsed.append(GenBank_Parser.parse(database, soup))
                except (NotAnORF, SequenceNotFoundError): #some entries have no usable DNA sequence
                    pass
            elif database == 'CCDS':
                soup = BeautifulSoup(text, 'html.parser')
                parsed.append(CCDS_Parser.parse(soup))
            elif database == 'Gene':
                soup = BeautifulSoup(text, 'xml')
                parsed.append(Gene_Parser.get_crossreferences(soup))
            else:
                raise ValueError(f'Unknown database: {database}')
        return parsed

class SummaryParser(object):
    
    def get_id_and_size(self, summary_page):
        soup = BeautifulSoup(summary_page, 'xml')
        id_ = soup.find_all('Item', {'Name': 'AccessionVersion'})[0].text
        id_ = id_.split('.')[0] #sometime it's ABCD1234.1, we want only ABCD1234
        size = soup.find_all('Item', {'Name': 'Length'})[0].text
        return id_, int(size)

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
    
    def match_entry_to_id(self, id_list):
        '''finds out which id_ of id_list matches with the current uniprot entry'''
        protein_ids = self.xml_soup.find_all('accession')
        protein_ids = set([x.text for x in protein_ids])
        id_list = [i.upper() for i in id_list] #just in case
        return list(set(id_list).intersection(protein_ids))
    
    def get_protein_sequence(self):
        sequence = self.xml_soup.find_all('sequence', 
                                          {'checksum': re.compile('^[A-Z0-9]+$')}) #checksum is alphanumeric uppercase
        return sequence[0].text.replace('\n','')
    
    def get_pdb_constructs(self, parent_id=None):
        if parent_id is None:
            parent_id = self.get_uniprot_id()
        pdbs = self.xml_soup.find_all('dbReference', {'type': 'PDB'})
        pdb_crossrefs = []
        for crossref in pdbs:
            '''crossref format example:
            <dbReference type="PDB" id="2HS7">
               <property type="method" value="X-ray" />
                <property type="chains" value="A=19-147" />
            </dbReference>
            also possible:
            ...
            <property type="chains" value="A=19-88, B=90-147" />
            ...'''
            try:
                chains_token = crossref.find_all('property', {'type': 'chains'})[0].attrs['value']
            except IndexError: 
                '''PDB codes are also present in the evidence section
                but lack start, stop, chains definitions. These are redundant, we ignore.
                Format is:
                <evidence ...>
                  <source><dbReference id="1IEE" type="PDB"/></source>
                </evidence>'''
                if 'evidence' in str(crossref.parent.parent):
                    continue
                else:
                    raise
            pdb_code = crossref.attrs['id']
            try:
                if ',' not in chains_token: #i.e. "A=19-147"
                    pdb_chains = chains_token.split('=')[0].split('/')
                    start, stop = [int(i) for i in chains_token.split('=')[1].split('-')]
                    pdb_crossrefs.append(PdbCrossref(parent_id, pdb_code, 
                                                     start, stop, pdb_chains))
                elif ',' in chains_token: #i.e. "A=19-88, B=90-147"
                    for chain in chains_token.replace(' ', '').split(','):
                        pdb_chains = chain.split('=')[0].split('/')
                        start, stop =  [int(i) for i in chain.split('=')[1].split('-')]
                        pdb_crossrefs.append(PdbCrossref(parent_id, pdb_code, 
                                                     start, stop, pdb_chains))
            except ValueError:
                print(f'Offending token is {chains_token} for pdb_code '
                      f'{pdb_code} in entry {parent_id}')
                continue
        return pdb_crossrefs
    
    def get_uniprot_id(self):
        #uniprot id is at beginning of page, within <accession> tags
        #top level tag indicates the latest (non deprecated) id
        ids = self.xml_soup.find_all('accession')
        return str(ids[0].next_element)
    
    def get_ensembl_references(self, isoform_number=1):
        '''
        If multiple isoforms:
        <dbReference type="Ensembl" id="ENST00000359435">
            <molecule id="Q9NWV8-1"/>
            <property type="protein sequence ID" value="ENSP00000352408"/>
            <property type="gene ID" value="ENSG00000105393"/>
        </dbReference>
        If single isoform:
        <dbReference type="Ensembl" id="ENSMFAT00000005779">
            <property type="protein sequence ID" value="ENSMFAP00000031561"/>
            <property type="gene ID" value="ENSMFAG00000036845"/>
        </dbReference>
        '''
        regexp = re.compile(r'Ensembl+')
        ensembl_references = self.xml_soup.find_all('dbReference', {'type': regexp})
        if not len(ensembl_references):
            return None #no references found!
        if '<isoform>' in str(self.xml_soup): #multiple isoforms present:
            single_isoform = False 
            isoform = f'{self.get_uniprot_id()}-{isoform_number}'
        else:
            single_isoform = True
        #this for loop will succeed if there is only one isoform, or if there is 
        #an ensembl reference for our specific isoform
        for r in ensembl_references:
            soup2 = BeautifulSoup(str(r), 'xml')
            if single_isoform or isoform in str(r): #if multiple isoforms, we only want ours
                gene_id = soup2.find_all('property', {'type': "gene ID"})
                return gene_id[0].attrs['value']
        #there might be ensembl references, but not for this isoform; we try again
        if not single_isoform: 
            for r in ensembl_references:
                soup2 = BeautifulSoup(str(r), 'xml')
                gene_id = soup2.find_all('property', {'type': "gene ID"})
                return f'*{gene_id[0].attrs["value"]}' #we mark with asterisk because match not exact
            else:
                raise IGiveUpError #something wrong with the entry?
    
class GeneParser(object):
     
    def __init__(self):
        super(GeneParser, self).__init__()
     
    def get_crossreferences(self, gene_xml_soup):
        ccds_entries = list(set(
            [c for c in gene_xml_soup.find_all(text=re.compile(
                                                      r'CCDS[0-9]+[.]?[0-9]*'))]))
        return ccds_entries
 
class GenBankParser(object):
     
    def __init__(self):
        super(GenBankParser, self).__init__()
     
    def parse(self, database, xml_soup):
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
        return CodingSequence(id_, database, cds, cds.translate(cds=True))
     
    def parse_non_eukaryotes(self, database, xml):
        xml_soup = BeautifulSoup(xml, 'xml')
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
                    coding_sequences.append(CodingSequence(id_, database, orf, protein_seq))
                    found = True
                except TranslationError: 
                    continue #not a good CDS
        if not found:  # not sure this ever happens
            raise NotAnORF
        return coding_sequences

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
        dna_seq = Seq(''.join([nt.text for nt in nucleotides]), 
                      IUPACUnambiguousDNA())
        aa_seq = Seq(''.join([aa.text for aa in aminoacids]), ExtendedIUPACProtein())
        assert aa_seq == dna_seq.translate(cds=True) 
        return CodingSequence(id_, 'CCDS', dna_seq, aa_seq)

class BlastRunner(object):
    
    def __init__(self):
        super().__init__()
        self.dbs = {'swissprot' : app.config['LOCAL_SWISSPROT_DATABASE'],
                          'trembl' : app.config['LOCAL_TREMBL_DATABASE'],
                          'pdb': app.config['LOCAL_PDB_DATABASE']}
        self.return_formats = {'xml' : '5',
                              'csv': '10'}
        self.blastp_exe = app.config['BLAST']
        self.check_blast_executable_exists()
        
    def check_blast_executable_exists(self):
        if not os.path.isfile(self.blastp_exe):
            raise OSError('Cannot find the blast executable (blastp). Please install blast')
    
    def check_search_parameters(self, query, dbname, return_format):
        if not query:
            raise ValueError('Please provide a query')
        try:
            format_code = self.return_formats[return_format]
        except KeyError:
            raise ValueError(f'Unknown format: {return_format}') 
        try:
            dbfile = self.dbs[dbname]
        except KeyError:
            raise ValueError(f'Unknown database {dbname}') 
        #check database exists at path
        if not glob.glob(f'{dbfile}.*'):
            raise OSError(f'The specified database: {dbfile} does not exist')
        #check that database is compiled with makeblastdb
        phr = glob.glob(f'{dbfile}*.phr') #files have format database_name.01.phr or database_name.phr
        pin = glob.glob(f'{dbfile}*.pin')
        psq = glob.glob(f'{dbfile}*.psq')
        if not (phr and pin and psq):
            raise OSError(f'Compiled database files {dbfile}.* not found. Did you compile the database?')
        return format_code, dbfile 
    
    def run_blast(self, query, dbname, return_format, fast=True):
        format_code, dbfile = self.check_search_parameters(query, dbname, return_format)
        if not query.startswith('>'):
            query = '>Unknown_query\n' + query
        #blastp works with file inputs. delete=False is for windows to work properly
        with NamedTemporaryFile(delete=False) as temp_in, NamedTemporaryFile(delete=False) as temp_out:
            #write query
            temp_in.write(query.encode('utf-8'))
            temp_in.seek(0)
            #run blast
            cmd = (f'{self.blastp_exe} -query {temp_in.name}'
                   f' -db {dbfile} -out {temp_out.name}'
                   f' -outfmt {format_code}')
            if fast:
                cmd += ' -task blastp-fast'
            p = Popen(cmd.split(), stdout=PIPE, stderr=PIPE)
            o, e = p.communicate() #blastp uses stdout to output errors
            if e or o:
                raise OSError(e)
            return temp_out.read().decode('utf-8')

class BlastParser(object):
    
    def __init__(self):
        super().__init__()
        self.filter_criteria = ['query_start','query_stop',
                                'identity', 'evalue']
        
    def parse_xml(self, blast_result):
        BlastHit = namedtuple('BlastHit', ['uniprot_id','query_start','query_stop',
                                           'identity', 'evalue'])
        soup = BeautifulSoup(blast_result, 'xml')
        hits = soup.find_all('Hit')
        results = []
        for hit in hits:
            entry = BeautifulSoup(str(hit), 'xml')
            #uniprot id is formatted as sp|Q9UJ41|RABX5_HUMAN bla bla bla
            uniprot_id = entry.find('Hit_accession').text 
            query_start = int(entry.find('Hsp_query-from').text)
            query_stop = int(entry.find('Hsp_query-to').text)
            identical = int(entry.find('Hsp_identity').text)
            aligned_length = int(entry.find('Hsp_align-len').text)
            identity = float(identical/aligned_length)*100
            evalue = float(entry.find('Hsp_evalue').text)
            results.append(BlastHit(uniprot_id, query_start, query_stop,
                                    identity, evalue ))
        return results
    
    def parse_csv(self, blast_result):
        BlastHit = namedtuple('BlastHit', ['uniprot_id','query_start','query_stop',
                                           'identity', 'evalue'])
        hits = [line for line in blast_result.split('\n') if line]
        results = []
        for line in hits:
            fields = line.split(',')
            uniprot_id = fields[1].split('|')[1]
            query_start = int(fields[6])
            query_stop = int(fields[7])
            identity = float(fields[2])
            evalue = float(fields[10])
            results.append(BlastHit(uniprot_id, query_start, query_stop,
                                    identity, evalue ))
        return results
    
    def filter_results(self, parsed_results, filter_by, min_val, max_val):
        if not filter_by in self.filter_criteria:
            msg = (f'Unknown filtering value {filter_by}. Valid choices are: '
                   f'{" ".join(self.filter_criteria)}')
            raise ValueError(msg)
        filtered = [res for res in parsed_results 
                    if getattr(res, filter_by) >= min_val
                    if getattr(res, filter_by) <= max_val]
        return filtered
    
    
            
        
                
        
             
        
    