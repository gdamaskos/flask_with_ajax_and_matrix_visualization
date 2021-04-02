import asyncio
import unittest

from pprint import pprint
from requests.exceptions import ConnectionError, HTTPError

from bs4 import BeautifulSoup

from ccd import app
import ccd.db_utilities_async as dba
from ccd.custom_exceptions import (IdNotFoundError, SequenceNotFoundError,
                                   NotAnORF) 
from ccd.sequence_classes import CodingSequence, PdbCrossref
from _ast import Name

class TestBlastParser(unittest.TestCase):
    
    def setUp(self):
        self.parser = dba.BlastParser()
    
    def test_parse_csv(self):
        query = f'>\n{q9uj41_1_protein_sequence}'
        blast_results = dba.BlastRunner().run_blast(query, 'swissprot', 'csv')
        res = self.parser.parse_csv(blast_results)
        self.assertEqual(res[0].uniprot_id, 'Q9UJ41')
        self.assertEqual(res[0].query_start, 139)
        self.assertEqual(res[0].query_stop, 708)
        self.assertAlmostEqual(res[0].identity, 86.14)
        self.assertAlmostEqual(res[0].evalue, 0.0)
        
    def test_parse_xml(self):
        query = f'>\n{q9uj41_1_protein_sequence}'
        blast_results = dba.BlastRunner().run_blast(query, 'swissprot', 'xml')
        res = self.parser.parse_xml(blast_results)
        self.assertEqual(res[0].uniprot_id, 'Q9UJ41')
        self.assertEqual(res[0].query_start, 139)
        self.assertEqual(res[0].query_stop, 708)
        self.assertAlmostEqual(res[0].identity, 86.14, 2)
        self.assertAlmostEqual(res[0].evalue, 0, 2)
    
    def test_filter_results(self):
        query = f'>\n{q9uj41_1_protein_sequence}'
        blast_results = dba.BlastRunner().run_blast(query, 'swissprot', 'csv')
        parsed = self.parser.parse_csv(blast_results)
        for filter_by in [('evalue', 0, 0.001, 48), 
                          ('query_start', 1, 139, 8),
                          ('query_stop', 350, 491, 25),
                          ('identity', 70, 100, 16)]:
            filtered = self.parser.filter_results(parsed, *filter_by[:3])
            self.assertEqual(filter_by[3], len(filtered), 
                             f'Failed to filter on {filter_by[0]}')
        with self.assertRaises(ValueError):
            filtered = self.parser.filter_results(parsed, 'invalid', 1, 2)
            
        
    
class toDO(unittest.TestCase):
    
    def test(self):
        msg = 'fix the GenBank Parser tests, you couch potato!'
        raise NotImplementedError(msg)
    
class Test_entry_splitter(unittest.TestCase):
    
    def test_split_entries_Refseq(self):
        xml_response = '<GBSet><GBSeq>a</GBseq><GBSeq>b</GBseq></GBSet>'
        database = 'RefSeq'
        splitter = dba.EntrySplitter()
        exp = [(database, '<GBSeq>a</GBSeq>'),
               (database, '<GBSeq>b</GBSeq>')]
        res = splitter.split(
            [(database, xml_response)])
        self.assertEqual(res, exp)
        
    def test_split_entries_EMBL(self):
        xml_response = '<GBSet><GBSeq>a</GBseq><GBSeq>b</GBseq></GBSet>'
        database = 'EMBL'
        splitter = dba.EntrySplitter()
        exp = [(database, '<GBSeq>a</GBSeq>'),
               (database, '<GBSeq>b</GBSeq>')]
        res = splitter.split(
            [(database, xml_response)])
        self.assertEqual(res, exp)
        
    def test_split_entries_Gene(self):
        xml_response = '<Entrezgene-Set><Entrezgene>a</Entrezgene><Entrezgene>b</Entrezgene></Entrezgene-Set>'
        database = 'Gene'
        splitter = dba.EntrySplitter()
        exp = [(database, '<Entrezgene>a</Entrezgene>'),
               (database, '<Entrezgene>b</Entrezgene>')]
        res = splitter.split(
            [(database, xml_response)])
        self.assertEqual(res, exp)
    
    def test_split_entries_summary(self):
        xml_response = '<eSummaryResult><DocSum>a</DocSum><DocSum>b</DocSum></eSummaryResult>'
        database = 'Summary'
        splitter = dba.EntrySplitter()
        exp = [(database, '<DocSum>a</DocSum>'),
               (database, '<DocSum>b</DocSum>')]
        res = splitter.split([(database, xml_response)])
        self.assertEqual(res, exp)
    
    def test_nodatabase(self):
        with self.assertRaises(ValueError):
            splitter = dba.EntrySplitter()
            splitter.split(
                [('NotADatabase', '')])
    
    def test_split_entries_Uniprot(self):
        xml_response = '''<uniprot>
        <entry><accession>P13368</accession></entry>
        <entry><accession>P13456</accession></entry>
        </uniprot>'''
        database = 'Uniprot'
        splitter = dba.EntrySplitter()
        exp = [(database, '<entry><accession>P13368</accession></entry>'),
               (database, '<entry><accession>P13456</accession></entry>')]
        res = splitter.split([(database, xml_response)])
        self.assertEqual(res, exp)
        
class Test_url_formatter(unittest.TestCase):
    
    def setUp(self):
        app.config['EMAIL'] = 'proteinccd@gmail.com'
        self.email = app.config['EMAIL']
    
    def test_EMBL_format_single_id(self):
        id_ = 'CR456855.1'
        exp = [('EMBL',
               f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={id_}&retmode=xml&rettype=gb&email={self.email}&tool=ccd.rhpc.nki.nl')]
        formatter = dba.UrlFormatter()
        res = formatter.format('EMBL', id_)
        self.assertEqual(res, exp)
    
    def test_EMBL_format_list_of_ids(self):
        id_ = ['CR456855.1', 'GBYX01232236.1']
        url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={",".join(id_)}&retmode=xml&rettype=gb&email={self.email}&tool=ccd.rhpc.nki.nl'
        exp = [('EMBL', url)]
        formatter = dba.UrlFormatter()
        res = formatter.format('EMBL', id_)
        self.assertEqual(res, exp)
        
    def test_RefSeq_single_id(self):
        id_ = 'NM_001270952.1'
        exp = [('RefSeq',
               f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={id_}&retmode=xml&rettype=gb&email={self.email}&tool=ccd.rhpc.nki.nl')]
        formatter = dba.UrlFormatter()
        res = formatter.format('RefSeq', id_)
        self.assertEqual(res, exp)
        
    def test_RefSeq_format_list_of_ids(self):
        id_ = ['NM_001270952.1', 'NM_006002.4']
        url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={",".join(id_)}&retmode=xml&rettype=gb&email={self.email}&tool=ccd.rhpc.nki.nl'
        exp = [('RefSeq', url)]
        formatter = dba.UrlFormatter()
        res = formatter.format('RefSeq', id_)
        self.assertEqual(res, exp)
        
    def test_CCDS_single_id(self):
        id_ = 'CCDS73586.1'
        exp = [('CCDS',
               f'https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&DATA={id_}&ORGANISM=0&BUILDS=CURRENTBUILDS')]
        formatter = dba.UrlFormatter()
        res = formatter.format('CCDS', id_)
        self.assertEqual(res, exp)
        
    def test_CCDS_format_list_of_ids(self):
        id_list = ['CCDS73586.1', 'CCDS86041.1']
        urls = [f'https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&DATA={id_}&ORGANISM=0&BUILDS=CURRENTBUILDS'
                for id_ in id_list]
        exp = [('CCDS', url) for url in urls]
        formatter = dba.UrlFormatter()
        res = formatter.format('CCDS', id_list)
        self.assertEqual(res, exp)
        
    def test_Gene_format_single_id(self):
        id_ = '7347'
        exp = [('Gene',
               f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=Gene&id={id_}&retmode=xml&email={self.email}&tool=ccd.rhpc.nki.nl')]
        formatter = dba.UrlFormatter()
        res = formatter.format('Gene', id_)
        self.assertEqual(res, exp)
    
    def test_Gene_format_list_of_ids(self):
        id_ = ['7347','50933']
        url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=Gene&id={",".join(id_)}&retmode=xml&email={self.email}&tool=ccd.rhpc.nki.nl'
        exp = [('Gene', url)]
        formatter = dba.UrlFormatter()
        res = formatter.format('Gene', id_)
        self.assertEqual(res, exp)
        
    def test_invalid_database(self):
        formatter = dba.UrlFormatter()
        with self.assertRaises(ValueError):
            formatter.format('NotADatabase', ['id1','id2'])
        with self.assertRaises(ValueError):
            formatter.format_for_summary('NotADatabase', ['id1','id2'])
            
    def test_format_for_summary(self):
        query1 = ('EMBL', ['CR456855.1', 'DQ917642.1']) #list of ids
        query2 = ('RefSeq', 'NM_001270952.1') #single id
        exp = [('Summary', 
                'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id=CR456855.1,DQ917642.1'), 
               ('Summary', 
                'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id=NM_001270952.1')]
        formatter = dba.UrlFormatter()
        res = formatter.format_for_summary(*query1)
        res += formatter.format_for_summary(*query2)
        self.assertEqual(exp, res)
    
    def test_uniprot_all_formats(self):
        id_ = 'q9uj41'
        formatter = dba.UrlFormatter()
        for format_ in ['xml','fasta','html']:
            url = f'https://www.uniprot.org/uniprot/{id_}.{format_}'
            exp = [('Uniprot', url)]
            res = formatter.format('Uniprot', id_, format_=format_)
            self.assertEqual(res, exp)
            
    def test_uniprot_no_format(self):
        formatter = dba.UrlFormatter()
        with self.assertRaises(ValueError):
            formatter.format('Uniprot', 'q9uj41') #no format
        with self.assertRaises(ValueError):
            formatter.format('Uniprot', 'q9uj41', format_='notvalid') #unsupported format
            
    def test_uniprot_batch(self):
        formatter = dba.UrlFormatter()
        id_list = ['q9uj41','Q15287','Q96RL1']
        formats = ['xml', 'html', 'fasta'] #batch added below
        for f in formats: 
            exp = [('Uniprot',
                    f'https://www.uniprot.org/uploadlists?query=q9uj41 Q15287 Q96RL1&from=ACC+ID&to=ACC&format={f}')]
            res = formatter.format('Uniprot', id_list, format_=f'batch_{f}')
            self.assertEqual(res, exp)
    
class Test_summary_parser(unittest.TestCase):
    
    def test_it(self):
        summary = '''<DocSum>
                     <Item Name="AccessionVersion" Type="String">CR456855.1</Item>
                     <Item Name="Length" Type="Integer">693</Item>
                     </DocSum>'''
        exp = ('CR456855', 693)
        parser = dba.SummaryParser()
        res = parser.get_id_and_size(summary)
        self.assertEqual(res, exp)

class TestRestLauncher(unittest.TestCase):
    
    def __init__(self, *args, **kwargs):
        super(TestRestLauncher, self).__init__(*args, **kwargs)
    
    def test_HTTPError(self):
        launcher = dba.RestLauncher()
        fake_url = 'https://ThisURLDoesNot.exist'
        with self.assertRaises(ConnectionError):
            launcher.run_query(fake_url, headers={}) 
    
class TestUniprotFetcher(unittest.TestCase):
    
    def __init__(self, *args, **kwargs):
        super(TestUniprotFetcher, self).__init__(*args, **kwargs)
        self.formats = ['xml', 'html', 'sequence']
    
    def test_wrong_format(self):
        fetcher = dba.UniprotFetcher()
        with self.assertRaises(ValueError):
            fetcher.fetch('Q9UJ41', 'NonSupportedFormat')
    
    def test_nonexistent_id(self):
        fetcher = dba.UniprotFetcher()
        nonexistent_id = 'ThisDoesNotExist'
        for format_ in self.formats:
            with self.assertRaises(IdNotFoundError):
                fetcher.fetch(nonexistent_id, format_)
            
    def test_returns_beautifulsoup(self):
        fetcher = dba.UniprotFetcher()
        result = fetcher.fetch('q9uj41', 'xml')
        self.assertTrue(type(BeautifulSoup('<a></a>', 'xml') == type(result)))
        result = fetcher.fetch('q9uj41', 'html')
        self.assertTrue(type(BeautifulSoup('<a></a>', 'html.parser') == type(result)))
    
    def test_returns_fasta(self):
        fetcher = dba.UniprotFetcher()
        result = fetcher.fetch('q9uj41', 'sequence')
        expected = q9uj41_2_protein_seq #sequence at end of file
        self.assertEqual(result, expected)

class TestUniprotParser(unittest.TestCase):
    
    def __init__(self, *args, **kwargs):
        super(TestUniprotParser, self).__init__(*args, **kwargs)
        fetcher = dba.UniprotFetcher()
        self.xml_soup = fetcher.fetch('Q9NWV8', 'xml')
        self.html_soup = fetcher.fetch('Q9NWV8', 'html')
        self.protein_seq = fetcher.fetch('Q9NWV8', 'sequence')
        self.maxDiff = None
        
    def setUp(self):
        self.parser = dba.UniprotParser(self.xml_soup, self.html_soup)
        
    def test_get_scientific_name(self):
        self.assertEqual(self.parser.get_scientific_name(),
                         'Homo sapiens' )
        
    def test_get_common_name(self):
        self.assertEqual(self.parser.get_common_name(), 
                        'Human')
        
    def test_common_name_not_in_entry(self):
        fetcher = dba.UniprotFetcher()
        self.parser.xml_soup = fetcher.fetch('Q7JKC3', 'xml')
        self.assertEqual(self.parser.get_common_name(), '')
    
    def test_get_protein_names(self):
        expected = ('BRISC and BRCA1-A complex member 1', None) 
        self.assertEqual(self.parser.get_protein_names(), 
                        expected)
    
    def test_get_kingdom(self):
        self.assertEqual(self.parser.get_kingdom(), 
                        'Eukaryota')
        
    def test_get_isoform_descriptions_malformed_uniprot_in_entry(self):
        #sometimes the uniprot entry has a difference that causes the parsing 
        #of the descriptions to glitch. Q9NWV8 is one such offenders
        exp = {u'Q9NWV8-2': [u'The sequence of this isoform differs from the canonical sequence as follows:\n', 
                             u'1-115: MEVAEPSSPT...SLPKLESFNG \u2192 MMGASTLQEP...HHGSTVQRKC\n', 
                             u'263-329: Missing\n'], 
               u'Q9NWV8-3': [u'The sequence of this isoform differs from the canonical sequence as follows:\n', 
                             u'116-190: Missing.\n', 
                             u'234-329: KMFQCPYFFF...DEAIEVEATV \u2192 VGEAWEREGCLQP\n'], 
               u'Q9NWV8-1': [u'This isoform has been chosen by Uniprot as the "canonical" sequence.\n']}
        res = self.parser.get_isoform_descriptions()
        self.assertEqual(res, exp)
    
    def test_get_isoform_descriptions_normal(self):
        fetcher = dba.UniprotFetcher()
        self.parser.html_soup = fetcher.fetch('Q9uj41', 'html')
        exp = {u'Q9UJ41-4': [u'The sequence of this isoform differs from the canonical sequence as follows:\n',
                             u'1-8: MVVVTGRE \u2192 MMASSYHE\n',
                             u'9-132: Missing.\n',
                             u'199-238: Missing.\n',
                             u'377-415: Missing\n'],
               u'Q9UJ41-1': [u'This isoform has been chosen by Uniprot as the "canonical" sequence.\n'],
               u'Q9UJ41-3': [u'The sequence of this isoform differs from the canonical sequence as follows:\n',
                             u'1-138: Missing.\n',
                             u'377-415: Missing\n'],
               u'Q9UJ41-2': [u'The sequence of this isoform differs from the canonical sequence as follows:\n',
                             u'1-138: Missing.\n',
                             u'199-238: Missing.\n',
                             u'377-415: Missing\n']}
        self.assertEqual(self.parser.get_isoform_descriptions(), exp)
        
    def test_get_crossreferences(self):
        parser = dba.UniprotParser(self.xml_soup, self.html_soup)
        exp = {'CCDS': ['CCDS46012.1'],
               'RefSeq': ['NM_001033549.2', 'NM_001288756.1',
                          'NM_001288757.1', 'NM_014173.3'],
               'EMBL': ['AF161491', 'AL136692', 'AK000578', 'AK299493', 'AK301193',
                        'CR533526', 'BC000788', 'BC006244', 'BC091491']}
#                'Translated':['AAH06244.1', 'CAB66627.1', 'BAG61451.1', 
#                              'CAG38557.1', 'AAH91491.1', 'BAG62773.1', 
#                              'BAA91268.1', 'AAH00788.1', 'AAF29106.1']}
        self.assertEqual(exp, parser.get_crossreferences())
    
    def test_match_entry_to_id(self):
        id_list = ['Q9NWV8', 'notreally', 'noteither']
        parser = dba.UniprotParser(self.xml_soup, self.html_soup)
        exp = ['Q9NWV8']
        res = parser.match_entry_to_id(id_list)
        self.assertEqual(res, exp)
        id_list = ['notreally', 'noteither']
        exp = []
        res = parser.match_entry_to_id(id_list)
        self.assertEqual(res, exp)
    
    def test_get_protein_sequence(self):
        parser = dba.UniprotParser(self.xml_soup, self.html_soup)
        exp = self.protein_seq
        res = parser.get_protein_sequence()
        self.assertEqual(res, exp)
    
    def test_get_uniprot_id(self):
        entry = '''
        <entry dataset="x" created="x" modified="x" version="1">
        <accession>P00698</accession>
        <accession>Q90884</accession>
        </entry>
        '''    
        parser = dba.UniprotParser(BeautifulSoup(entry, 'xml'), None)
        exp = 'P00698'
        self.assertEqual(exp, parser.get_uniprot_id())

    def test_get_pdb_constructs(self):
        entry = '''
        <entry>
            <accession>P00698</accession>
            <dbReference type="PDB" id="132L">
                <property type="method" value="X-ray" />
                <property type="resolution" value="1.80" />
                <property type="chains" value="A=19-147" />
            </dbReference>
            <dbReference type="PDB" id="4WL6">
                <property type="method" value="X-ray" />
                <property type="resolution" value="1.85" />
                <property type="chains" value="A=1-147" />
            </dbReference>
            <dbReference type="PDB" id="2X89">
                <property type="method" value="X-ray" />
                <property type="resolution" value="2.16" />
                <property type="chains" value="D/E/F/G=27-119" />
            </dbReference>
            <evidence ...>
                <source><dbReference id="1IEE" type="PDB"/></source>
            </evidence>   
        </entry>
        '''
        #the pdb code in evidence must be ignored
        exp = [PdbCrossref('P00698', '132L', 19, 147, ['A']),
               PdbCrossref('P00698', '4WL6', 1, 147, ['A']),
               PdbCrossref('P00698', '2X89', 27, 119, ['D','E','F','G'])]
        parser = dba.UniprotParser(BeautifulSoup(entry, 'xml'), None)
        res = parser.get_pdb_constructs()
        self.assertEqual(sorted(res), sorted(exp))
        
    def test_get_pdb_construct_none_present(self):
        exp = [] #should return empty list
        entry = '''<entry><accession>P18181</accession></entry>'''
        parser = dba.UniprotParser(BeautifulSoup(entry, 'xml'), None)
        res = parser.get_pdb_constructs()
        self.assertEqual(res, exp)
        
    def test_get_pdb_construct_uniprot_id_given(self):
        entry = '''
        <entry>
            <accession>P00698</accession>
            <dbReference type="PDB" id="132L">
                <property type="method" value="X-ray" />
                <property type="resolution" value="1.80" />
                <property type="chains" value="A=19-147" />
            </dbReference>
        </entry>
        '''
        exp = [PdbCrossref('override', '132L', 19, 147, ['A'])]
        parser = dba.UniprotParser(BeautifulSoup(entry, 'xml'), None)
        res = parser.get_pdb_constructs('override')
        self.assertEqual(sorted(res), sorted(exp))
        self.assertEqual(res[0].parent_entry_id, 'override')
    
    def test_get_pdb_construct_double_chain_notation(self):
        entry = '''
        <entry>
            <accession>P00698</accession>
            <dbReference type="PDB" id="132L">
                <property type="method" value="X-ray" />
                <property type="resolution" value="1.80" />
                <property type="chains" value="A=19-88, B=90-147" />
            </dbReference>
        </entry>
        '''    
        exp = [PdbCrossref('P00698', '132L', 19, 88, ['A']),
               PdbCrossref('P00698', '132L', 90, 147, ['B'])]
        parser = dba.UniprotParser(BeautifulSoup(entry, 'xml'), None)
        res = parser.get_pdb_constructs()
        self.assertEqual(sorted(res), sorted(exp))
        
    def test_get_ensembl_id_isoform1(self):
        fetcher = dba.UniprotFetcher()
        xml_soup = fetcher.fetch('Q9NWV8', 'xml')
        parser = dba.UniprotParser(xml_soup, None)
        exp = 'ENSG00000105393'
        res = parser.get_ensembl_references()
        self.assertEqual(res, exp)
        
    def test_get_ensembl_id_isoform_not1(self):
        fetcher = dba.UniprotFetcher()
        xml_soup = fetcher.fetch('Q9UJ41', 'xml')
        parser = dba.UniprotParser(xml_soup, None)
        exp = 'ENSG00000154710'
        res = parser.get_ensembl_references(isoform_number=2)
        self.assertEqual(res, exp)
        res = parser.get_ensembl_references(isoform_number=4)
        self.assertEqual(res, exp)
    
    def test_get_ensembl_id_single_isoform(self):
        fetcher = dba.UniprotFetcher()
        xml_soup = fetcher.fetch('P61258', 'xml')
        parser = dba.UniprotParser(xml_soup, None)
        exp = 'ENSMFAG00000036845'
        res = parser.get_ensembl_references()
        self.assertEqual(res, exp)
        
    def test_get_ensembl_id_non_matching_isoform(self):
        fetcher = dba.UniprotFetcher()
        xml_soup = fetcher.fetch('Q9UJ41', 'xml')
        parser = dba.UniprotParser(xml_soup, None)
        exp = '*ENSG00000154710'
        res = parser.get_ensembl_references(isoform_number=3)
        self.assertEqual(res, exp)       

class TestUniprotMapper(unittest.TestCase):
    
    def __init__(self, *args, **kwargs):
        super(TestUniprotMapper, self).__init__(*args, **kwargs)
        
    def test_nonexisting_id(self):
        mapper = dba.UniprotMapper()
        with self.assertRaises(IdNotFoundError):
            mapper.map_('NonexistingId', 'ID', 'P_ENTREZGENEID')
    
    def test_nonexisting_from_database(self):
        mapper = dba.UniprotMapper()
        with self.assertRaises(ValueError):
            mapper.map_('38144', 'ID', 'NotADatabase')
    
    def test_nonexisting_to_database(self):
        mapper = dba.UniprotMapper()
        with self.assertRaises(ValueError):
            mapper.map_('38144', 'NotADatabase', 'P_ENTREZGENEID')
            
    def test_uniprot_to_gene(self):
        mapper = dba.UniprotMapper()
        result = mapper.map_('Q9UJ41', 'ID', 'P_ENTREZGENEID')
        exp = '27342'
        self.assertEqual(result, exp)
        
    def test_list_of_ids(self):
        mapper = dba.UniprotMapper()
        result = mapper.map_(['P00698','P00918'], 'ID', 'ENSEMBL_ID')
        exp = {'P00698': 'ENSGALG00000009963', 
               'P00918': 'ENSG00000104267'}
        self.assertEqual(result, exp)

class TestGeneParser(unittest.TestCase):
    
    def __init__(self, *args, **kwargs):
        super(TestGeneParser, self).__init__(*args, **kwargs)
    
    def test_correct_functioning(self):
        loop = asyncio.get_event_loop()
        parser = dba.GeneParser()
        uf = dba.UrlFormatter()
        fetcher = dba.Entry_fetcher()
        query = uf.format('Gene','27342')
        res = loop.run_until_complete(fetcher.fetch_all(query))
        xml_soup = BeautifulSoup(res[0][1], 'xml')
        exp = ['CCDS69308.1', 'CCDS5535.1', 'CCDS75610.1']
        self.assertEqual(exp.sort(), 
                         parser.get_crossreferences(xml_soup).sort())

class TestGenBankParser(unittest.TestCase):
    '''
    test_euk* tests check parser.parse()
    test_bact* tests check parser.parse_non_eukaryotes()
    they are very different in expected behaviour...
    '''
    def __init__(self, *args, **kwargs):
        super(TestGenBankParser, self).__init__(*args, **kwargs)
        self.parser = dba.GenBankParser()
        
    def _fetch_data(self, database, id_list):
        loop = asyncio.get_event_loop()
        uf = dba.UrlFormatter()
        fetcher = dba.Entry_fetcher()
        query = uf.format(database,id_list)
        return loop.run_until_complete(fetcher.fetch_all(query))
        
    def test_euk_correct_functioning(self):
        id_list = ['AJ250042']
        database = 'EMBL'
        exp = CodingSequence(id_list[0], database,
                             q9uj41_2_dna_seq, q9uj41_2_protein_seq)
        entry = self._fetch_data(database, id_list)[0][1]
        res = self.parser.parse(database , BeautifulSoup(entry, 'xml')) 
        self.assertEqual(res, exp)
    
    def test_bact_normal_function(self):
        id_list = ['AH002539']
        database = 'EMBL'
        exp = CodingSequence(id_list[0], database,
                             rpsL_ecoli_dna_seq, rpsL_ecoli_protein_seq)
        entry = self._fetch_data(database, id_list)[0][1]
        res = self.parser.parse_non_eukaryotes(database, entry) 
        self.assertIn(exp, res) #others will be there - multiple cds
        
    def test_bact_non_standard_start_codon(self):
        id_list = ['AH002539']
        database = 'EMBL'
        exp = CodingSequence(id_list[0], database,
                             eftu1_ecoli_dna_seq, eftu1_ecoli_protein_seq)
        entry = self._fetch_data(database, id_list)[0][1]
        res = self.parser.parse_non_eukaryotes(database, entry) 
        self.assertIn(exp, res) #others will be there - multiple cds
        
    def test_bact_orf_in_complement(self):
        id_list = ['AF058450']
        database = 'EMBL'
        exp = CodingSequence(id_list[0], database, 
                             hopd_ecolx_dna_seq, hopd_ecolx_protein_seq)
        entry = self._fetch_data(database, id_list)[0][1]
        res = self.parser.parse_non_eukaryotes(database, entry)
        self.assertIn(exp, res) #others will be there - multiple cds
    
    def test_euk_not_mrna_or_dna_entry(self):
        #eukaryote
        entry = self.fetcher.fetch(['AC027644'])
        with self.assertRaises(NotAnORF):
            self.parser.parse(entry[0])
        #prokaryote - fails if run with eukaryotic parser
        entry = self.fetcher.fetch(['Z37500']) 
        with self.assertRaises(NotAnORF):
            self.parser.parse(entry[0])

    def test_euk_incomplete_orf(self):
        entry = self.fetcher.fetch(['Z14117'])
        with self.assertRaises(NotAnORF):
            self.parser.parse(entry[0])
            
    def test_bact_incomplete_orf(self):
        entry = self.fetcher.fetch(['M10459'])
        with self.assertRaises(NotAnORF):
            self.parser.parse(entry[0])
    
    def test_euk_badly_defined_cds(self): #missing atg despite being declared a complete cds
        entry = self.fetcher.fetch(['M11848'])
        with self.assertRaises(NotAnORF):
            self.parser.parse(entry[0])
    
    def test_bact_entry_with_no_sequence(self): 
        #sequence reconstructed programmatically from other sources but not in entry
        entry = self.fetcher.fetch(['NZ_NMWW01000089.1'])
        with self.assertRaises(SequenceNotFoundError):
            self.parser.parse_non_eukaryotes(entry[0])
                        
class TestBlastRunner(unittest.TestCase):
    
    def __init__(self, *args, **kwargs):
        super(TestBlastRunner, self).__init__(*args, **kwargs)
        self.runner = dba.BlastRunner()
    
    def test_wrong_arguments(self):
        query = ''
        with self.assertRaises(ValueError): 
            self.runner.run_blast(query, 'swissprot', 'xml')
        query = f'>\n{q9uj41_1_protein_sequence}'
        with self.assertRaises(ValueError): 
            self.runner.run_blast(query, 'wrongdb', 'xml')
        with self.assertRaises(ValueError): 
            self.runner.run_blast(query, 'swissprot', 'wrongformat')
    
    def test_no_database(self):
        query = f'>\n{q9uj41_1_protein_sequence}'
        with self.assertRaises(OSError):
            self.runner.dbs['swissprot'] = 'notafile'
            self.runner.run_blast(query, 'swissprot', 'xml')
        from tempfile import NamedTemporaryFile
        db= NamedTemporaryFile()
        with self.assertRaises(OSError):
            self.runner.dbs['swissprot'] = db.name
            self.runner.run_blast(query, 'swissprot', 'xml')
            
    def test_csv_format(self):
        query = f'>\n{q9uj41_1_protein_sequence}'
        res = self.runner.run_blast(query, 'swissprot', 'csv')
        exp = 'unnamed,sp|Q9UJ41|RABX5_HUMAN,86.14,570,0,2,139,708,1,491,0.0,984'
        self.assertEqual(exp, res.split('\n')[0])
        
    def test_xml_format(self):
        query = f'>\n{q9uj41_1_protein_sequence}'
        res = self.runner.run_blast(query, 'swissprot', 'xml')
        exp = 'sp|Q9UJ41|RABX5_HUMAN Rab5 GDP/GTP exchange factor OS=Homo sapiens OX=9606 GN=RABGEF1 PE=1 SV=3'
        soup = BeautifulSoup(res, 'xml')
        hits = soup.find_all('Hit_def')
        self.assertEqual(exp, hits[0].text)
        
        


q9uj41_1_protein_sequence = '''MVVVTGREPDSRRQDGAMSSSDAEDDFLEPATPTATQAGHALPLLPQERCAEFPALRGPP
TQGACSSCVQRGPVLCHRAPPGAAGEHAATEGREGAPSVSGTHALLQRPLGADCGDRPAA
CGPAEGPLCQAQVVSRKKMSLKSERRGIHVDQSDLLCKKGCGYYGNPAWQGFCSKCWREE
YHKARQKQIQEDWELAERVLLCCPGWSAMVQFQLTATSASWAQVILLLQPPKWLGLQKLQ
REEEEAFASSQSSQGAQSLTFSKFEEKKTNEKTRKVTTVKKFFSASSRVGSKKEIQEAKA
PSPSINRQTSIETDRVSKEFIEFLKTFHKTGQEIYKQTKLFLEGMHYKRDLSIEEQSECA
QDFYHNVAERMQTRGKERRFHHVGQAGLELLTSGDPPASASQSAGNTGVEPPHPAVPPER
VEKIMDQIEKYIMTRLYKYVFCPETTDDEKKDLAIQKRIRALRWVTPQMLCVPVNEDIPE
VSDMVVKAITDIIEMDSKRVPRDKLACITKCSKHIFNAIKITKNEPASADDFLPTLIYIV
LKGNPPRLQSNIQYITRFCNPSRLMTGEDGYYFTNLCCAVAFIEKLDAQSLNLSQEDFDR
YMSGQTSPRKQEAESWSPDACLGVKQMYKNLDLLSQLNERQERIMNEAKKLEKDLIDWTD
GIAREVQDIVEKYPLEIKPPNQPLAAIDSENVENDKLPPPLQPQVYAG'''.replace('\n','')

q9uj41_2_dna_seq = 'ATGAGCCTTAAGTCTGAACGCCGAGGAATTCATGTGGATCAATCGGATCTCCTGTGCAAGAAAGGATGTGGTTACTACGGCAACCCTGCCTGGCAGGGTTTCTGCTCCAAGTGCTGGAGGGAAGAGTACCACAAAGCCAGGCAGAAGCAGATTCAGGAGGACTGGGAGCTGGCGGAGCGACTCCAGCGGGAGGAAGAAGAGGCCTTTGCCAGCAGTCAGAGCAGCCAAGGGGCCCAATCCCTCACATTCTCCAAGTTTGAAGAAAAGAAAACCAACGAGAAGACCCGCAAGGTTACCACAGTGAAGAAATTCTTCAGTGCATCTTCCAGGGTCGGATCAAAGAAGGAAATTCAGGAAGCAAAAGCTCCCAGTCCTTCCATAAACCGGCAAACCAGCATTGAAACGGATAGAGTGTCTAAGGAGTTCATAGAATTTCTCAAGACCTTCCACAAGACAGGCCAAGAAATCTATAAACAGACCAAGCTGTTTTTGGAAGGAATGCATTACAAAAGGGATCTAAGCATTGAAGAACAGTCAGAGTGTGCTCAGGATTTCTACCACAATGTGGCCGAAAGGATGCAAACTCGTGGGAAAGTGCCTCCAGAAAGAGTCGAGAAGATAATGGATCAGATTGAAAAGTACATCATGACTCGTCTCTATAAATATGTATTCTGTCCAGAAACTACTGATGATGAGAAGAAAGATCTTGCCATTCAAAAGAGAATCAGAGCCCTGCGCTGGGTTACGCCTCAGATGCTGTGTGTCCCTGTTAATGAAGACATCCCAGAAGTGTCTGATATGGTGGTGAAGGCGATCACAGATATCATTGAAATGGATTCCAAGCGTGTGCCTCGAGACAAGCTGGCCTGCATCACCAAGTGCAGCAAGCACATCTTCAATGCCATCAAGATCACCAAGAATGAGCCGGCGTCAGCGGATGACTTCCTCCCCACCCTCATCTACATTGTTTTGAAGGGCAACCCCCCACGCCTTCAGTCTAATATCCAGTATATCACGCGCTTCTGCAATCCAAGCCGACTGATGACTGGAGAGGATGGCTACTATTTCACCAATCTGTGCTGTGCTGTGGCTTTCATTGAGAAGCTAGACGCCCAGTCTTTGAATCTAAGTCAGGAGGATTTTGATCGCTACATGTCTGGCCAGACCTCTCCCAGGAAGCAAGAAGCTGAGAGTTGGTCTCCTGATGCTTGCTTAGGCGTCAAGCAAATGTATAAGAACTTGGATCTCTTGTCTCAGTTGAATGAACGACAAGAAAGGATCATGAATGAAGCCAAGAAACTGGAAAAAGACCTCATAGATTGGACAGATGGAATTGCAAGAGAAGTTCAAGACATCGTTGAGAAATACCCACTGGAAATTAAGCCTCCGAATCAACCGTTAGCAGCTATTGACTCTGAAAACGTTGAAAATGATAAACTTCCTCCACCACTGCAACCTCAAGTTTATGCAGGATGA'
q9uj41_2_protein_seq = '''MSLKSERRGIHVDQSDLLCKKGCGYYGNPAWQGFCSKCWREEYHKARQKQIQEDWELAER
LQREEEEAFASSQSSQGAQSLTFSKFEEKKTNEKTRKVTTVKKFFSASSRVGSKKEIQEA
KAPSPSINRQTSIETDRVSKEFIEFLKTFHKTGQEIYKQTKLFLEGMHYKRDLSIEEQSE
CAQDFYHNVAERMQTRGKVPPERVEKIMDQIEKYIMTRLYKYVFCPETTDDEKKDLAIQK
RIRALRWVTPQMLCVPVNEDIPEVSDMVVKAITDIIEMDSKRVPRDKLACITKCSKHIFN
AIKITKNEPASADDFLPTLIYIVLKGNPPRLQSNIQYITRFCNPSRLMTGEDGYYFTNLC
CAVAFIEKLDAQSLNLSQEDFDRYMSGQTSPRKQEAESWSPDACLGVKQMYKNLDLLSQL
NERQERIMNEAKKLEKDLIDWTDGIAREVQDIVEKYPLEIKPPNQPLAAIDSENVENDKL
PPPLQPQVYAG'''.replace('\n','')

hopd_ecolx_dna_seq='''ttaacg gcctgccagc
aaactatccc aacctacaac gaaacccgca gccgccagaa atggcccaaa cggcagcggg
ttttttaatg actcttttcc tctcatgagc agaccaacaa ctacggctcc gcaagcgaat
gaggcagcaa ggaaaaccag tcgtggcaag aaggtccagg tatgccaggc accaagggcg
gcgagaaact tcacatcacc atagcctaat ccttctttat gacgcagtat gcgatacccc
cagtaaatga cagcgaatgt gccgtaacca ataatcgccc cccataacgc atcggctaaa
cagtcaggat tacaaacctg agagaacaaa agacctgacc agagtaacgg gcaggtaaat
ctgtcaggta acaagccatg ctttgcatcc cataaaaaaa gtaggacaga taaacaagcg
tacaaaatca gaaagggaag ggtagcgtcc at'''.replace('\n','').replace(' ','').upper()
hopd_ecolx_protein_seq='''MDATLPFLILYACLSVLLFLWDAKHGLLPDRFTCPLLWSGLLFSQVCNPDCLADALWGAI
IGYGTFAVIYWGYRILRHKEGLGYGDVKFLAALGAWHTWTFLPRLVFLAASFACGAVVVG
LLMRGKESLKNPLPFGPFLAAAGFVVGWDSLLAGR'''.replace('\n','')

eftu1_ecoli_protein_seq='''MSKEKFERTKPHVNVGTIGHVDHGKTTLTAAITTVLAKTYGGAA
RAFDQIDNAPEEKARGITINTSHVEYDTPTRHYAHVDCPGHADYVKNMITGAAQMDGA
ILVVAATDGPMPQTREHILLGRQVGVPYIIVFLNKCDMVDDEELLELVEMEVRELLSQ
YDFPGDDTPIVRGSALKALEGDAEWEAKILELAGFLDSYIPEPERAIDKPFLLPIEDV
FSISGRGTVVTGRVERGIIKVGEEVEIVGIKETQKSTCTGVEMFRKLLDEGRAGENVG
VLLRGIKREEIERGQVLAKPGTIKPHTKFESEVYILSKDEGGRHTPFFKGYRPQFYFR
TTDVTGTIELPEGVEMVMPGDNIKMVVTLIHPIAMDDGLRFAIREGGRTVGAGVVAKV
LG'''.replace(' ','').replace('\n','')
eftu1_ecoli_dna_seq='''gtgtctaaag aaaaatttga acgtacaaaa ccgcacgtta
acgttggtac tatcggccac gttgaccacg gtaaaactac tctgaccgct gcaatcacca
ccgtactggc taaaacctac ggcggtgctg ctcgtgcatt cgaccagatc gataacgcgc
cggaagaaaa agctcgtggt atcaccatca acacttctca cgttgaatac gacaccccga
cccgtcacta cgcacacgta gactgcccgg ggcacgccga ctatgttaaa aacatgatca
ccggtgctgc tcagatggac ggcgcgatcc tggtagttgc tgcgactgac ggcccgatgc
cgcagactcg tgagcacatc ctgctgggtc gtcaggtagg cgttccgtac atcatcgtgt
tcctgaacaa atgcgacatg gttgatgacg aagagctgct ggaactggtt gaaatggaag
ttcgtgaact tctgtctcag tacgacttcc cgggcgacga cactccgatc gttcgtggtt
ctgctctgaa agcgctggaa ggcgacgcag agtgggaagc gaaaatcctg gaactggctg
gcttcctgga ttcttatatt ccggaaccag agcgtgcgat tgacaagccg ttcctgctgc
cgatcgaaga cgtattctcc atctccggtc gtggtaccgt tgttaccggt cgtgtagaac
gcggtatcat caaagttggt gaagaagttg aaatcgttgg tatcaaagag actcagaagt
ctacctgtac tggcgttgaa atgttccgca aactgctgga cgaaggccgt gctggtgaga
acgtaggtgt tctgctgcgt ggtatcaaac gtgaagaaat cgaacgtggt caggtactgg
ctaagccggg caccatcaag ccgcacacca agttcgaatc tgaagtgtac attctgtcca
aagatgaagg cggccgtcat actccgttct tcaaaggcta ccgtccgcag ttctacttcc
gtactactga cgtgactggt accatcgaac tgccggaagg cgtagagatg gtaatgccgg
gcgacaacat caaaatggtt gttaccctga tccacccgat cgcgatggac gacggtctgc
gtttcgcaat ccgtgaaggc ggccgtaccg ttggcgcggg cgttgttgct aaagttctgg
gctaa'''.replace('\n','').replace(' ','').upper()

rpsL_ecoli_protein_seq='''MATVNQLVRKPRARKVAKSNVPALEACPQKRGVCTRVYTTTPKK
PNSALRKVCRVRLTNGFEVTSYIGGEGHNLQEHSVILIRGGRVKDLPGVRYHTVRGAL
DCSGVKDRKQARSKYGVKRPKA'''.replace(' ','').replace('\n','')
rpsL_ecoli_dna_seq='''atggcaac agttaaccag ctggtacgca aaccacgtgc tcgcaaagtt gcgaaaagca
acgtgcctgc gctggaagca tgcccgcaaa aacgtggcgt atgtactcgt gtatatacta
ccactcctaa aaaaccgaac tccgcgctgc gtaaagtatg ccgtgttcgt ctgactaacg
gtttcgaagt gacttcctac atcggtggtg aaggtcacaa cctgcaggag cactccgtga
tcctgatccg tggcggtcgt gttaaagacc tcccgggtgt tcgttaccac accgtacgtg
gtgcgcttga ctgctccggc gttaaagacc gtaagcaggc tcgttccaag tatggcgtga
agcgtcctaa ggcttaa'''.replace('\n','').replace(' ','').upper()

 
 
 