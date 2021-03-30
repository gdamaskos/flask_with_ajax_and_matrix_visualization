# coding: utf8
import re 
import urllib.request, urllib.error, urllib.parse
import logging
import sys

from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA

from ccd import app
from ccd.sequence_classes import (UniprotProtein, UniprotIsoform, 
                                  CodingSequence, CandidateORF)
from ccd.db_utilities_async import (UniprotMapper, UniprotFetcher, UniprotParser,
                              GeneFetcher, GeneParser,
                              GenBankFetcher, GenBankParser,
                              CCDSFetcher, CCDSParser)
#all exceptions are simply pass, nothing special
from ccd.custom_exceptions import (IdNotFoundError, SequenceNotFoundError, 
                                   NotAnORF, IGiveUpError, NoSequenceFound)  

#TODO: check that an email is properly set in ccd_settings. We want this compulsory.
HEADERS = {'User agent': 'CrystallizationConstructDesigner_{}'.format(app.config['EMAIL'])}

class ORFMatcher(object):
    
    def __init__(self):
        super(ORFMatcher, self).__init__()
    
    def match_exact(self, uniprot_isoform, cached_cds):
        for c in cached_cds:
            if uniprot_isoform.protein_seq == c.protein_seq:
                uniprot_isoform.dna_seq = c.dna_seq  # TODO : account for possible different cds with same protein_seq
                uniprot_isoform.matched_with.append(c)
        return uniprot_isoform #annotated with all matched cds

    def match_fuzzy(self, uniprot_isoform, cached_cds, max_mismatches=3):
        '''
        Compares the protein sequences of uniprot_isoform and candidate_cds
        position by position. 
        Returns None if more than max_err differences are found.
        If less than max_error differences are found, the isoforms are considered
        matched, and uniprot_isoform is updated to reflect the match 
        '''
        n_errors = 0
        i = uniprot_isoform #I'm a lazy typist
        for c in cached_cds:
            mismatches = []
            # different length (plus/minus tolerance) implies mismatch and is faster to check
            length_mismatch = len(i.protein_seq) - len(c.protein_seq)
            if abs(length_mismatch) > max_mismatches:
                continue 
            # otherwise we just match position by position until we find more 
            #than mismatches than max_mismatches or we run out of positions => match
            for pos, aa in enumerate(i.protein_seq):
                if n_errors > max_mismatches:
                    break
                if aa != c.protein_seq[pos]:
                    mismatch = ['{aa}{pos}'.format(aa=aa, pos=pos),
                                '{aa2}'.format(aa2=c.protein_seq[pos])]
                    mismatches.append(mismatch) 
                    n_errors += 1
            i.dna_seq = c.dna_seq.upper()
            i.matched_with.append(c)
            i.mismatches = mismatches
        return i
    
class ORFFinder(object):
    
    def __init__(self):
        super(ORFFinder, self).__init__()
    
    def _find_boundaries(self, genomic_dna_seq):
        # find all possible starts, stops
        try:
            _ = genomic_dna_seq.alphabet #we want a Biopython Seq object or similar
        except AttributeError:
            genomic_dna_seq = Seq(genomic_dna_seq, IUPACUnambiguousDNA()) 
        start = re.compile('ATG')
        stop = re.compile('TAA|TGA|TAG')
        stops = re.finditer(stop, str(genomic_dna_seq).upper())
        stops = [s.start() + 2 for s in stops]  # stop codon ends at +3 from regexp match
        starts = re.finditer(start, str(genomic_dna_seq).upper())
        starts = [s.start() for s in starts]
        candidate_boundaries = []
        for start in starts:
            for stop in stops:
                #ORFs shorter than 20 aminoacids == 60 base pairs don't make sense
                sensible = stop - start >= 60   
                not_in_frame = (stop - start) % 3 != 2  # expect 2 if start & stop are in frame
                if not_in_frame or not sensible:
                    continue
                else:
                    candidate_boundaries.append((start, stop))
        return candidate_boundaries
        
    def find_orfs(self, dna_seq):
        possible_orfs = []
        candidate_boundaries = self._find_boundaries(dna_seq) #list of tuples
        for start, stop in candidate_boundaries:
            orf_sequence = Seq(dna_seq[start-1, stop+1], IUPACUnambiguousDNA())
            possible_orfs.append(CandidateORF(None, None, orf_sequence, start, stop))        
        return possible_orfs 
    
class DatabaseCrawler(object):
    # TODO: write proper docstring
    '''
    given a database id or a dna sequence, fetches all info necessary to start CCD 
    '''

    def __init__(self, search_id, max_mismatches=3):
        super(DatabaseCrawler, self).__init__()
        self.search_id = search_id
        self.max_mismatches = max_mismatches
        self.protein = None
        logging.basicConfig(stream=sys.stdout, level=logging.INFO)
        self.logger = logging.getLogger(__name__)
        
    def search(self):
        #mapping the given id to uniprot primary accession
        #we could use the given id directly, but secondary accession are not good for
        #downstream applications
        #raises IfNotFoundError if the id does not exist in Uniprot
        try:
            self.uniprot_id = self.map_id_to_uniprot(self.search_id)
        except IdNotFoundError:
            msg = 'Failed in retrieving a valid uniprot id using {}'.format(
                                                                    self.search_id)
            self.logger.debug(msg)
            raise
        msg = 'Mapped user-provided id {} to {} primary accession number'.format(
                                                                self.search_id,
                                                                self.uniprot_id) 
        self.logger.info(msg)
        # fetching uniprot entry and identifying protein isoforms
        self.protein = self.fetch_and_parse_uniprot_data(self.uniprot_id)
        # fetching the corresponding entry from the NCBI Gene database
        # returns None if none exists
        self.protein.gene_xml_soup = self.search_gene_database()
        #fetch all crossreferenced database entries from uniprot and gene entries
        self.protein.db_crossrefs = self.scrape_crossreferences()
        #fetch and cache all entries corresponding to cDNAs
        if self.protein.kingdom == 'Eukaryota':
            self.protein.cached_cds = self.fetch_crossreferences(
                                                        self.protein.db_crossrefs)
            self.protein.isoforms = self.match_eukaryotes()
        elif self.protein.kingdom in ['Archaea', 'Bacteria']:
            self.protein.cached_cds = self.fetch_crossreferences_non_eukaryota(
                                                        self.protein.db_crossrefs)
            self.protein.isoforms = self.match_non_eukaryotes()
        # checking the results of the matching operations
        if not self.protein.isoforms:
            msg = 'No crossreferences to DNA sequences found in Uniprot entry {}'.format(self.protein.id_)
            self.logger.info(msg)
            raise IGiveUpError
        unmatched_isoforms = [i for i in self.protein.isoforms if not i.matched_with]
        if unmatched_isoforms:
            msg = 'Warning: {} isoforms could not be matched to any DNA within {} substitutions'.format(
                                                                    len(unmatched_isoforms),
                                                                    self.max_mismatches)
            self.logger.info(msg)
        if len(unmatched_isoforms) == len(self.protein.isoforms):
            msg = 'None of the {} known isoforms could be matched'.format(
                                                        len(self.protein.isoforms))
        return self.protein
    
    def fetch_and_parse_uniprot_data(self, uniprot_id):
        p = UniprotProtein(uniprot_id)
        # 1 Let's retrieve the uniprot entry in html and xml form
        self.logger.info('Fetching Uniprot data for {}'.format(uniprot_id))
        fetcher = UniprotFetcher()
        p.uniprot_xml_soup = fetcher.fetch(uniprot_id, 'xml')
        p.uniprot_html_soup = fetcher.fetch(uniprot_id, 'html')
        # 2 Grab some basic protein data from uniprot
        self.logger.info('Parsing Uniprot data for {}'.format(uniprot_id))
        parser = UniprotParser(p.uniprot_xml_soup, p.uniprot_html_soup)
        p.source_organism = parser.get_scientific_name()
        p.source_organism_common = parser.get_common_name()
        p.full_name, p.short_name = parser.get_protein_names()
        try:
            p.kingdom = parser.get_kingdom()
        except IGiveUpError: 
            msg = 'Cannot identify kingdom for {}'.format(self.uniprot_id)
            raise IGiveUpError(msg) #We need kingdom to proceed with cds_matching
        # 3 let's parse the entries for isoform ids, descriptions and sequences 
        isoform_descriptions = parser.get_isoform_descriptions()
        if not isoform_descriptions: #none found in the page
            isoform_descriptions = {uniprot_id: ['This is the only known isoform']}
        # 4let's get the sequences of the isoforms from uniprot 
        #(more robust than parsing the entry ourselves, i think)
        isoforms = []
        for id_ in isoform_descriptions:
            protein_seq = fetcher.fetch(id_, 'sequence')
            isoforms.append(UniprotIsoform(id_=id_,
                                           protein_seq=protein_seq))
        # 5 let's match descriptions and sequences
        for i in isoforms:
            i.description = isoform_descriptions[i.id_]
        p.isoforms = isoforms
        return p
    
    def search_gene_database(self):
        # fetching gene entry and parsing data
        mapper = UniprotMapper()
        try:
            self.logger.info('Searching for entries in Gene')
            gene_id = mapper.map_(self.search_id, from_='ACC+ID', to='P_ENTREZGENEID')
            gene_fetcher = GeneFetcher()
            return  gene_fetcher.fetch(gene_id)
        except IdNotFoundError:
            msg = 'Id {} not found in Entrez Gene; Skipping.'.format(self.protein.id_)
            self.logger.info(msg)
            return None
    
    def scrape_crossreferences(self):
        #start from the uniprot page
        uniprot_parser = UniprotParser(self.protein.uniprot_xml_soup, 
                               self.protein.uniprot_html_soup)
        dbreferences = uniprot_parser.get_crossreferences()
        #get more crossreferences to CCDS from the NCBI gene page
        if self.protein.gene_xml_soup:
            try:
                gene_parser = GeneParser()
                ccds_refs = gene_parser.get_crossreferences(self.protein.gene_xml_soup)
                #some duplicates might be present 
                dbreferences['CCDS'] = list(set(dbreferences['CCDS'] + ccds_refs))
            except AttributeError:  # no gene entry present
                pass
        # If any crossreference has been found, we are good. Otherwise we give up
        has_entries = [True for i in list(dbreferences.values()) if i] #any will do
        if has_entries:
            return dbreferences
        else:
            msg = 'No crossreferences to DNA sequences found in Uniprot entry {}'.format(self.protein.id_)
            self.logger.info(msg)
            raise IGiveUpError
    
    def fetch_crossreferences(self, crossrefs):
        cached_entries = {}
        #Refseq and EMBL support batch
        genbank_fetcher = GenBankFetcher()
        for database in ['EMBL', 'RefSeq']: #we can batch fetch these
            self.logger.info('Fetching {} entries from {}'.format(
                                                    len(crossrefs[database]),
                                                    database))
            cached_entries[database]  = genbank_fetcher.fetch(crossrefs[database]) 
        #CCDS cross references are fetched one by one
        CCDS_fetcher = CCDSFetcher()
        cached_entries['CCDS'] = []
        self.logger.info('Fetching {} entries from CCDS'.format(
                                                    len(crossrefs['CCDS'])))
        for id_ in crossrefs['CCDS']: 
            cached_entries['CCDS'].append(CCDS_fetcher.fetch(id_))
        return cached_entries
    
    def fetch_crossreferences_non_eukaryota(self, crossrefs):
        genbank_fetcher = GenBankFetcher()
        #first we check the size of these things - we don't want to fetch anything too big
        #we cannot really parse entire genomes anyway...
        self.logger.info('Checking that the DNA references are not too big to fetch')
        reasonable = {'EMBL':[],'RefSeq':[]}
        for database in ['RefSeq', 'EMBL']:
            for ref in crossrefs[database]:
                if genbank_fetcher.check_size(ref) < 150000: #arbitrary 150 kBases
                    reasonable[database].append(ref)
        reasonable['CCDS'] = crossrefs['CCDS'] #always empty unless mouse or human
        self.logger.info('Fetching DNA references')
        return self.fetch_crossreferences(reasonable)
        
    def match_eukaryotes(self):
        matcher = ORFMatcher()
        #parse the database entries to a uniform format
        self.protein.cached_cds = self.parse_crossrefs_eukaryota(self.protein.cached_cds)
        if not self.protein.cached_cds:
            return []
        #Attempt to match protein sequences exactly
        for isoform in self.protein.isoforms:
            isoform = matcher.match_exact(isoform, self.protein.cached_cds)
        exactly_matched = sum([1 for i in self.protein.isoforms if i.matched_with])
        self.logger.info('{} out of {} isoforms match a DNA sequence exactly'.format(
                    exactly_matched, len(self.protein.isoforms)))
        #Match any unmatched isoform within <max_mismatches> substitutions
        unmatched_isoforms = [i for i in self.protein.isoforms if not i.matched_with]
       
        if unmatched_isoforms: 
            for isoform in unmatched_isoforms:
                isoform = matcher.match_fuzzy(isoform, self.protein.cached_cds,
                                              max_mismatches=self.max_mismatches)
            fuzzy_matched = sum([1 for i in self.protein.isoforms if i.mismatches])
            self.logger.info('{} out of {} isoforms match a DNA sequence within {} substitutions'.format(
                fuzzy_matched, len(self.protein.isoforms), self.max_mismatches))
        return self.protein.isoforms
    
    def match_non_eukaryotes(self):
        matcher = ORFMatcher()
        annotated_orfs, non_annotated_orfs = self.parse_crossrefs_non_eukaryota(
                                                        self.protein.cached_cds)
        #annotated isoforms will 99% of times be a good match
        if annotated_orfs:
            for isoform in self.protein.isoforms:
                isoform = matcher.match_exact(isoform, annotated_orfs)
        exactly_matched = sum([1 for i in self.protein.isoforms if i.matched_with])
        unmatched_isoforms = [i for i in self.protein.isoforms if not i.matched_with]
        if not unmatched_isoforms:
            self.logger.info('{} out of {} isoforms match a DNA sequence exactly'.format(
                    exactly_matched, len(self.protein.isoforms)))
            return self.protein.isoforms
        else:
            raise IGiveUpError 

    def parse_crossrefs_eukaryota(self, cached_crossrefs):
        parsed_entries = []
        for database, entry_list in list(cached_crossrefs.items()):
            if not entry_list:
                self.logger.info('No entries found in {}'.format(database))
            if database in ['EMBL', 'RefSeq']:
                parser = GenBankParser()
            elif database == 'CCDS':
                parser = CCDSParser()
            for e in entry_list:
                try:
                    id_, dna_seq, protein_seq = parser.parse(e)
                    parsed_entries.append(CodingSequence(id_, database, dna_seq, protein_seq))
                except NotAnORF:
                    try:
                        self.logger.info('{} is not a complete ORF'.format(id_))
                    except UnboundLocalError:
                        self.logger.info('Cannot parse id of this {} entry'.format(database))
                    pass
                
        '''
        Sometime identical entries are present in both RefSeq and EMBL databases.
        We only need one copy. I don't care which one, they are identical anyway.
        I have overridden __hash__ in CodingSequence so that set can operate
        based on CodingSequence.id_ XOR str(CodingSequence.dna_seq)
        '''
        parsed_entries = list(set(parsed_entries))
        return parsed_entries

    def parse_crossrefs_non_eukaryota(self, cached_crossrefs):
        annotated_orfs = []
        non_annotated_orfs = []
        parser = GenBankParser()
        for database in ['RefSeq', 'EMBL']:
            for ref in cached_crossrefs[database]:
                try:
                    cds_coordinates = parser.parse_non_eukaryotes(ref)
                    for cds in cds_coordinates:
                        annotated_orfs.append(CodingSequence(cds[0],
                                                             database,
                                                             cds[1], cds[2]))
                except NotAnORF:
                    dna_seq = str(ref.find_all('GBSeq_sequence')[0].string).upper()
                    id_ = str(ref.find_all('GBSeq_primary-accession')[0].text).strip()
                    non_annotated_orfs.append((id_, dna_seq))
                except SequenceNotFoundError:
                    continue
        return annotated_orfs, non_annotated_orfs
        
    def cache_cds_prokaryotes(self, crossref):
        id_, database = crossref
        fetcher, parser, get_sequence_method = self.disambiguate(database)
#         fetch_method, _, get_sequence_method = self.disambiguate(database)        
        self.logger.info('{} has no annotated cds; looking for open reading frames'.format(
                        id_))
        possible_orfs = []
        try:
            entry = fetcher.fetch(id_)
        except urllib.error.HTTPError as err:
            if err.code == 400 or err.code == 404:
                self.logger.info('Cannot find entry {}'.format(id_))
                return None
        try:
            dna_seq = get_sequence_method(entry)
            possible_boundaries = self.find_boundaries(dna_seq)
            possible_orfs = []
            for i in possible_boundaries:
                start, stop = i[0], i[1]
                possible_orfs.append(CandidateORF(id_,
                                                  database,
                                                  dna_seq,
                                                  start,
                                                  stop))
            return possible_orfs 
        except NoSequenceFound:
            msg = 'No usable sequence found for entry {}'.format(id_)
            self.logger.info(msg)
            return None
        except UnboundLocalError as e:
            msg = f'Error encountered with id = {id_}'
            raise UnboundLocalError(msg) from e

    def match_orfs(self, isoforms, possible_orfs):
        for isoform in isoforms:
            # we check that the orf has the correct length to code for protein_seq
            pruned_orfs = [orf for orf in possible_orfs 
                           if len(orf.dna_seq) == len(isoform.protein_seq) * 3 + 3]
            for orf in pruned_orfs:
                t = str(orf.dna_seq.translate()).replace('*', '')  # stop codon = '*'
                if t == isoform.protein_seq:
                    orf.protein_seq = t
                    #we expect only one, but we return a list for coherence with 
                    #multiple matches that can happen on Codingsequences  
                    isoform.matched_with = [orf]
                    isoform.dna_seq = orf.dna_seq 
                    break
        return
    
    def map_id_to_uniprot(self, id_):
        '''
        Maps any protein id to a UNIPROT or TREMBL protein id, and returns it.
        Raises IdNotFoundError if no id is found
        '''
        mapper = UniprotMapper()
        try:
            new_id = mapper.map_(id_, 'ACC+ID', 'ACC') 
            return new_id
        except IdNotFoundError:
            raise
#             try:
#                 new_id = mapper.map_(id_, 'TREMBL')
#                 return new_id
#             except IdNotFoundError as e:
#                 msg = ('Id {} not found'.format(id_))
#                 raise IdNotFoundError(msg) from e #OK we give up and raise
    
if __name__ == '__main__':
#    the following code generates json files for black box testing in test_dna_finder
    import getpass, platform, os, json
    generate_test_jsons = True # creates references for blackbox testing of this unit
    if platform.system() == 'Windows':
        app.config['MUSCLE'] = r'D:\\Users\\{}\\Downloads\\muscle\\muscle3.8.31_i86win32.exe'.format(getpass.getuser())#
    #the following code will generate json files for testing
    max_errors = 3
    id_ = 'UIMC1_HUMAN'
    static_path = os.path.normpath(os.getcwd() + '/../tests/unit/static')
    n = DatabaseCrawler(id_)
    f = os.path.join(static_path, id_.upper() + '.json')
    protein = n.search()    
    
    if generate_test_jsons:
        with open(f, 'w') as f:
            json.dump(protein.serialize_protein(), f)
     
