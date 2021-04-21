# coding: utf8
# import urllib.request, urllib.error, urllib.parse
import logging
import sys
import time

from ccd import app
from ccd.sequence_classes import (UniprotProtein, UniprotIsoform, CandidateORF)
from ccd.db_utilities_async import *
from ccd.custom_exceptions import (IdNotFoundError, SequenceNotFoundError, 
                                   NotAnORF, IGiveUpError, NoSequenceFound)  

HEADERS = {'User agent': f'CrystallizationConstructDesigner_{app.config["EMAIL"]}'}
MAX_SIZE = 10_000 # 10 kBases of DNA should be sufficient 
class ORFMatcher(object):
    
    def __init__(self):
        super(ORFMatcher, self).__init__()
    
    def match_exact(self, uniprot_isoform, cached_cds):
        for c in cached_cds:
            if uniprot_isoform.protein_seq == c.protein_seq:
                uniprot_isoform.dna_seq = c.dna_seq
                uniprot_isoform.matched_with.append(c)
        return uniprot_isoform # annotated with all matched cds

    def match_fuzzy(self, uniprot_isoform, cached_cds, max_mismatches=3):
        '''
        Compares the protein sequences of uniprot_isoform and candidate_cds
        position by position. 
        Returns None if more than max_err differences are found.
        If less than max_error differences are found, the isoforms are considered
        matched, and uniprot_isoform is updated to reflect the match 
        '''
        n_errors = 0
        i = uniprot_isoform
        for c in cached_cds:
            mismatches = []
            # different length (plus/minus tolerance) implies mismatch and is faster to check
            length_mismatch = len(i.protein_seq) - len(c.protein_seq)
            if abs(length_mismatch) > max_mismatches:
                continue 
            # otherwise we just match position by position until we find more 
            # than mismatches than max_mismatches or we run out of positions => match
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
    
class DatabaseCrawler(object):
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
        self.loop = asyncio.get_event_loop()
        self.fetcher = Entry_fetcher()
        
    def search(self):
        '''
        first we assume the user has provided a valid uniprot entry, and we try 
        to fetch the data. If 404 (=> user supplied an invalid id), we try
        to map using the Uniprot uploadlists service. If that also fails, we raise
        IdNotFoudError
        All parsed data is stored in a UniprotProtein object 
        '''
        try:
            self.protein = self.fetch_uniprot_entry(self.search_id)
            self.uniprot_id = self.search_id
        except IdNotFoundError: # maybe obsoleted entry? let's try mapping it to a new version
            try:
                self.uniprot_id = self.map_id_to_uniprot(self.search_id)
                self.protein = self.fetch_uniprot_entry(self.uniprot_id)
                msg = 'Mapped user-provided id {} to {} primary accession number'.format(
                                                                self.search_id,
                                                                self.uniprot_id) 
                
                self.logger.info(msg)
            except IdNotFoundError:
                msg = 'Failed in retrieving a valid Uniprot id using {}'.format(
                                                                    self.search_id)
                self.logger.info(msg)
                raise
        self.protein = self.parse_uniprot_data(self.protein)
        # fetching the corresponding entry from the NCBI Gene database
        # returns None if none exists
        self.protein.gene_xml_soup = self.search_gene_database()
        # fetch all crossreferenced database entries from uniprot and gene entries
        self.protein.db_crossrefs = self.scrape_crossreferences()
        # in some rare instances, the gene database gets queried too fast and we get a 429
        # let's wait one second to be sure

        # fetch and cache all entries corresponding to cDNAs
        if self.protein.kingdom == 'Eukaryota':
            self.protein.cached_cds = self.fetch_crossreferences(self.protein.db_crossrefs)
            self.logger.info('Done fetching crossreferences')
            parser = DnaParser()
            parsed_entries = parser.parse(self.protein.cached_cds)
            self.protein.isoforms = self.match_eukaryotes(parsed_entries)
        elif self.protein.kingdom in ['Archaea', 'Bacteria']:
            self.protein.db_crossrefs = self.exclude_big_entries(self.protein.db_crossrefs)
            self.protein.cached_cds = self.fetch_crossreferences(self.protein.db_crossrefs)
            self.logger.info('Done fetching crossreferences')
            parsed_entries = self.parse_crossrefs_non_eukaryota(self.protein.cached_cds)
            self.protein.isoforms = self.match_non_eukaryotes(parsed_entries)
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
            self.logger.info(msg)
            raise IGiveUpError(msg)
        return self.protein
    
    def fetch_uniprot_entry(self, uniprot_id):
        '''
        Given a valid uniprot_id, fetches and returns the xml and html pages for the entry,
        as a beautifulsoup.
        Returns a UniprotPtorein object
        '''
        p = UniprotProtein(uniprot_id)
        # 1 let's retrieve the uniprot entry in html and xml form
        self.logger.info('Fetching Uniprot data for {}'.format(uniprot_id))
        fetcher = UniprotFetcher()
        p.uniprot_xml_soup = fetcher.fetch(uniprot_id, 'xml')
        p.uniprot_html_soup = fetcher.fetch(uniprot_id, 'html')
        p.uniprot_id = uniprot_id
        return p
    
    def parse_uniprot_data(self, protein_object):
        # 2 grab some basic protein data from uniprot
        p = protein_object
        self.logger.info('Parsing Uniprot data for {}'.format(p.uniprot_id))
        parser = UniprotParser(p.uniprot_xml_soup, p.uniprot_html_soup)
        p.uniprot_id = parser.get_uniprot_id()
        p.source_organism = parser.get_scientific_name()
        p.source_organism_common = parser.get_common_name()
        p.full_name, p.short_name = parser.get_protein_names()
        try:
            p.kingdom = parser.get_kingdom()
        except IGiveUpError: 
            msg = 'Cannot identify kingdom for {}'.format(self.uniprot_id)
            raise IGiveUpError(msg) # We need kingdom to proceed with cds_matching
        # 3 let's parse the entries for isoform ids, descriptions and sequences 
        isoform_descriptions = parser.get_isoform_descriptions()
        if not isoform_descriptions: #none found in the page
            isoform_descriptions = {p.uniprot_id: ['This is the only known isoform']}
        # 4 let's get the sequences of the isoforms from uniprot 
        # (more robust than parsing the entry ourselves, i think)
        isoforms = []
        fetcher = UniprotFetcher()
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
            try:
                gene_id = mapper.map_(self.search_id, from_='ACC+ID', to='P_ENTREZGENEID')
            except urllib.error.HTTPError:
                return None
            query = UrlFormatter().format('Gene', gene_id)
            res = self.loop.run_until_complete(self.fetcher.fetch_all(query))
            return BeautifulSoup(res[0][1], 'xml')
        except IdNotFoundError:
            msg = 'Id {} not found in Entrez Gene; Skipping.'.format(self.protein.id_)
            self.logger.info(msg)
            return None
    
    def scrape_crossreferences(self):
        # start from the uniprot page
        uniprot_parser = UniprotParser(self.protein.uniprot_xml_soup, 
                               self.protein.uniprot_html_soup)
        dbreferences = uniprot_parser.get_crossreferences()
        # get more crossreferences to CCDS from the NCBI gene page
        if self.protein.gene_xml_soup:
            try:
                gene_parser = GeneParser()
                ccds_refs = gene_parser.get_crossreferences(self.protein.gene_xml_soup)
                # some duplicates might be present 
                dbreferences['CCDS'] = list(set(dbreferences['CCDS'] + ccds_refs))
            except AttributeError:  # no gene entry present
                pass
        # if any crossreference has been found we are good, therwise we give up
        has_entries = [True for i in list(dbreferences.values()) if i] #any will do
        if has_entries:
            return dbreferences
        else:
            msg = 'No crossreferences to DNA sequences found in Uniprot entry {}'.format(self.protein.id_)
            self.logger.info(msg)
            raise IGiveUpError
    
    def exclude_big_entries(self, db_crossrefs, max_size=MAX_SIZE):
        formatter = UrlFormatter()
        splitter = EntrySplitter()
        summary_parser = SummaryParser()
        to_fetch = []
        query = []
        # generate urls for checking summary of all entries 
        for database, id_list in db_crossrefs.items():
            if database == 'CCDS':
                to_fetch += id_list
            elif database in ['EMBL','RefSeq']:
                query += (formatter.format_for_summary(database, id_list))
        # fetching and parsing summaries
        entry_summaries = self.loop.run_until_complete(self.fetcher.fetch_all(query))
        entry_summaries = splitter.split(entry_summaries) # list of (entry_id, entry_summary) tuples
        # removing all entries with size > max_size
        for s in entry_summaries:
            id_, size = summary_parser.get_id_and_size(s[1])
            if size <= max_size:
                to_fetch.append(id_)
        for database, entries in db_crossrefs.items(): 
            for id_ in list(entries):
                if id_ not in to_fetch:
                    entries.remove(id_)
        return db_crossrefs
        
    def prepare_queries(self, crossrefs):
        formatter = UrlFormatter()
        queries = []
        for database, id_list in crossrefs.items():
            queries += formatter.format(database, id_list)
            self.logger.info(f'Fetching {len(id_list)} entries from {database}')
        return queries
        
        
    def fetch_crossreferences(self, crossrefs):
        # preparing URLs to fetch asynchronously
        formatter = UrlFormatter()
        queries = []
        for database, id_list in crossrefs.items():
            queries += formatter.format(database, id_list)
            self.logger.info(f'Fetching {len(id_list)} entries from {database}')
        # fetching asynchronously
        entries = self.loop.run_until_complete(self.fetcher.fetch_all(queries))
        # EMBL, RefSeq are returned as a single page and need to be splitted.
        splitter = EntrySplitter()
        splitted = splitter.split(entries)
        return splitted
    
    def match_eukaryotes(self, parsed_entries):
        matcher = ORFMatcher()
        # parse the database entries to a uniform format
        if not self.protein.cached_cds:
            return []
        # attempt to match protein sequences exactly
        for isoform in self.protein.isoforms:
            isoform = matcher.match_exact(isoform, parsed_entries)
        exactly_matched = sum([1 for i in self.protein.isoforms if i.matched_with])
        self.logger.info('{} out of {} isoforms match a DNA sequence exactly'.format(
                    exactly_matched, len(self.protein.isoforms)))
        # match any unmatched isoform within <max_mismatches> substitutions
        unmatched_isoforms = [i for i in self.protein.isoforms if not i.matched_with]
       
        if unmatched_isoforms: 
            for isoform in unmatched_isoforms:
                isoform = matcher.match_fuzzy(isoform, parsed_entries,
                                              max_mismatches=self.max_mismatches)
            fuzzy_matched = sum([1 for i in self.protein.isoforms if i.mismatches])
            self.logger.info('{} out of {} isoforms match a DNA sequence within {} substitutions'.format(
                fuzzy_matched, len(self.protein.isoforms), self.max_mismatches))
        return self.protein.isoforms
    
    def match_non_eukaryotes(self, parsed_entries):
        matcher = ORFMatcher()
        for isoform in self.protein.isoforms:
            isoform = matcher.match_exact(isoform, parsed_entries)
        exactly_matched = sum([1 for i in self.protein.isoforms if i.matched_with])
        unmatched_isoforms = [i for i in self.protein.isoforms if not i.matched_with]
        if not unmatched_isoforms:
            self.logger.info('{} out of {} isoforms match a DNA sequence exactly'.format(
                    exactly_matched, len(self.protein.isoforms)))
            return self.protein.isoforms
        else:
            raise IGiveUpError

    def parse_crossrefs_non_eukaryota(self, cached_crossrefs):
        annotated_orfs = []
        parser = GenBankParser()
        for ref in cached_crossrefs:
            try:
                annotated_orfs += parser.parse_non_eukaryotes(*ref)
            except NotAnORF:
                continue
            except SequenceNotFoundError:
                continue
        return annotated_orfs

    def cache_cds_prokaryotes(self, crossref):
        id_, database = crossref
        fetcher, parser, get_sequence_method = self.disambiguate(database)
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
                    # we expect only one, but we return a list for coherence with 
                    # multiple matches that can happen on Codingsequences  
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
            return mapper.map_(id_, 'ACC+ID', 'ID')
        except IdNotFoundError:
            raise
    
if __name__ == '__main__':
    # the following code generates json files for black box testing in test_dna_finder
    import getpass, platform, os, json
    generate_test_jsons = True # creates references for blackbox testing of this unit
    if platform.system() == 'Windows':
        app.config['MUSCLE'] = r'D:\\Users\\{}\\Downloads\\muscle\\muscle3.8.31_i86win32.exe'.format(getpass.getuser())
    # the following code will generate json files for testing
    max_errors = 3
    id_ = 'EFTU2_ECOLI'
    static_path = os.path.normpath(os.getcwd() + '/../tests/unit/static')
    n = DatabaseCrawler(id_)
    f = os.path.join(static_path, id_.upper() + '.json')
    protein = n.search()    
    
    if generate_test_jsons:
        with open(f, 'w') as f:
            json.dump(protein.serialize_protein(), f)
     
