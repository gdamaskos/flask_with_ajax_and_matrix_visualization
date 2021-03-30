
import unittest
import json
import os
import logging
import requests

from ccd.dna_finder_async import DatabaseCrawler 
from ccd.custom_exceptions import IGiveUpError, IdNotFoundError

           
class TestBlackbox(unittest.TestCase):
    '''
    Tests dna_finder performance on actual protein entries.
    Results are manually compiled and are valid against Uniprot release 2017_07,
    Genbank release 220.0 (15/07) and CCDS release 21 (Oct 2016)
    The results of tests that use Uniprot ids are not expected to change over time
    The tests that use long uniprot accession (e..g baba2_human) might change. 
    Should a test fail, use your judgement!
    '''
    
    def __init__(self, *args, **kwargs):
        super(TestBlackbox, self).__init__(*args, **kwargs)
        self.static_path = os.path.abspath('tests/unit/static')
        if not os.path.isdir(self.static_path): #testing via eclipse vs. tesing command line
            self.static_path = os.path.abspath('static')
        self.debug = False
        if not self.debug:
            logging.basicConfig(level=logging.CRITICAL) #let the thing shut up
        self.test_logger = logging.getLogger('test_logger')
        self.test_logger.setLevel(level=logging.INFO)
        self.maxDiff = 100000
        
            
    def open_results_and_search(self, id_):
        f = os.path.join(self.static_path, id_ + '.json')
        with open(f, 'r') as f:
            expected = json.load(f)
        a = DatabaseCrawler(id_)
        return a.search().serialize_protein(), expected
    
    def cleanup_results(self, res):
        #results are returned as big dict, a.k.a. hard to parse any difference
        #also, depending on retrieval order, some crossreferences might be marked
        #as EMbL or Refseq, with no functional differences, but would fail comparison 
        #Here we extract and sort the meaningful parts
        isoforms = {}
        isoforms['isoform_ids'] = sorted([i['isoform_id'] 
                                   for i in res['protein_isoforms_list']])
        isoforms['isoforms_dna'] = [i['isoform_dna']
                                    for i in res['protein_isoforms_list']]
        #depending on random retrieval order, the dna might have a different stop 
        #codon to the one specified in the static reference
        for i in isoforms['isoforms_dna']:
            i = i[:-3] + 'TAA'
        isoforms['isoforms_dna'] = sorted(isoforms['isoforms_dna'])
        isoforms['isoforms_dna'] = isoforms['isoforms_dna'].sort() 
        isoforms['isoforms_short_protein'] = sorted([i['isoform_short_protein'] 
                                         for i in res['protein_isoforms_list']])
        matches = {}
        for i in res['protein_isoforms_list']:
            matches[i['isoform_id']] = {}
            matches[i['isoform_id']]['id_'] = sorted([m['cds_id'] 
                                               for m in i['isoform_matched_with']])
            matches[i['isoform_id']]['dna_seq'] = sorted([m['cds_dna'] 
                                               for m in i['isoform_matched_with']])
            matches[i['isoform_id']]['url'] = sorted([m['cds_url'] 
                                               for m in i['isoform_matched_with']])
        return isoforms, matches
            
    def test_archeal_sequence(self):
        self.test_logger.info('Testing an archaeal sequence')
        id_ = 'P61875' #Pyrococcus furiosus DNA pol
        res, exp = self.open_results_and_search(id_)
        self.assertEqual(self.cleanup_results(res), self.cleanup_results(exp))
     
    def test_multiple_isoforms_fuzzy_and_exact_matches(self):
        self.test_logger.info('Testing entry with multiple isoforms and exact + fuzzy matches')
        id_ = 'Q9UJ41' #2 exact, 2 fuzzy matches
        res, exp = self.open_results_and_search(id_)
        self.assertEqual(self.cleanup_results(res), self.cleanup_results(exp))
      
    def test_no_crossreferences_found(self):
        self.test_logger.info('Testing entry with no crossreferences')
        id_ = 'G1MX40'
        a = DatabaseCrawler(id_)
        with self.assertRaises(IGiveUpError):
            a.search()
    
    def test_no_crossreferences_match_protein(self):
        self.test_logger.info('Testing entry with no crossreferences matching the protein sequence')
        id_ = 'F4I443'
        a = DatabaseCrawler(id_)
        with self.assertRaises(IGiveUpError):
            a.search()
      
    def test_bad_uniprot_unicode_in_description(self):
        #this entry would fail at correctly parsing the description of isoform 3
        #because the html of uniprot would have the << character encoded as \xab
        #instead of the usual \xa0
        self.test_logger.info('Testing entry with inconsistent unicode in uniprot entry')
        id_ = 'Q9NWV8' #3 exact matches; entry point uniprot id
        res, exp = self.open_results_and_search(id_)
        self.assertEqual(self.cleanup_results(res), self.cleanup_results(exp))
 
    def test_entry_via_accession_instead_of_id(self):
        self.test_logger.info('Testing entry with accession UIMC1_HUMAN instead of ID Q96RL1')
        id_ = 'UIMC1_HUMAN' #3 exact matches; entry point uniprot id
        res, exp = self.open_results_and_search(id_)
        self.assertEqual(self.cleanup_results(res), self.cleanup_results(exp))
        
    def test_trembl_match_via_accession(self):
        self.test_logger.info('Testing entry with TREMBL accession H3AMY9_LATCH')
        id_ = 'H3AMY9_LATCH' 
        res, exp = self.open_results_and_search(id_)
        self.assertEqual(self.cleanup_results(res), self.cleanup_results(exp))
    
    def test_trembl_match_via_id(self):
        self.test_logger.info('Testing entry with TREMBL id A8DPD6')
        id_ = 'A8DPD6' 
        res, exp = self.open_results_and_search(id_)
        self.assertEqual(self.cleanup_results(res), self.cleanup_results(exp))
        
    def test_match_Refseq_only(self):
        self.test_logger.info('Testing entry that has no CCDS or Genbank matches')
        id_ = 'A7SD85'   
        res, exp = self.open_results_and_search(id_)
        self.assertEqual(self.cleanup_results(res), self.cleanup_results(exp))
    
    def test_match_Refseq_only_and_hallucinates_beautifulsoup(self):
        #For reasons I cannot fathom, this entry hallucinates Uniprot parser (Beautifulsoup)
        #the dbReference tag gets split in two and a crossreference to a protein 
        #slips through. This should cause a HTTP 404 in Genbank fetcher, that should
        #be silenced by the fetcher as a last resort. 
        self.test_logger.info('Testing entry that has no CCDS or Genbank matches, and is problematic')
        id_ = 'A7S327'
        res, exp = self.open_results_and_search(id_)
        self.assertEqual(self.cleanup_results(res), self.cleanup_results(exp))
        
#     def test_no_gene_entry(self): #uniprot changed; now I cannot find another entry that would satisfy test criteria
#         self.test_logger.info('Testing entry that has no Gene entry')
#         id_ = 'F4I443' 
#         res, exp = self.open_results_and_search(id_)
#         self.assertEqual(self.cleanup_results(res), self.cleanup_results(exp))
    
    def test_nonexisting_id(self):
        self.test_logger.info('Testing entry with no crossreferences')
        id_ = 'thismostdefinitelydoesnotexist'
        a = DatabaseCrawler(id_)
        with self.assertRaises(IdNotFoundError):
            a.search()    
    
    def test_bacterial_matching_polycystronic_mRNA(self):
        self.test_logger.info('Testing E. coli entry that matches one ORF from a multicistronic gene')
        id_ = 'Q47223' 
        res, exp = self.open_results_and_search(id_)
        self.assertEqual(self.cleanup_results(res), self.cleanup_results(exp))
        id_ = 'P0ACF8' 
        self.assertEqual(self.cleanup_results(res), self.cleanup_results(exp))
        self.assertEqual(res, exp)
    
    def test_HTTPError500(self): 
        id_ = '' 
        with self.assertRaises(requests.HTTPError):
            a = DatabaseCrawler(id_)
            a.search()
            
    def test_join_in_CCDS_start_stop(self):
        #some entries in CCDS can have the start and stop of the CDS
        #encoded as join(start.pos2,pos2..stop) e.g. uniprot id: Q68EN5
        self.test_logger.info('Testing Uniprot entry with non-joined CDS boundaries in CCDS')
        id_ = 'Q68EN5' 
        res, exp = self.open_results_and_search(id_)
        self.assertEqual(self.cleanup_results(res), self.cleanup_results(exp))

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
GIAREVQDIVEKYPLEIKPPNQPLAAIDSENVENDKLPPPLQPQVYAG'''.replace('\n', '')

q9uj41_2_dna_seq = 'ATGAGCCTTAAGTCTGAACGCCGAGGAATTCATGTGGATCAATCGGATCTCCTGTGCAAGAAAGGATGTGGTTACTACGGCAACCCTGCCTGGCAGGGTTTCTGCTCCAAGTGCTGGAGGGAAGAGTACCACAAAGCCAGGCAGAAGCAGATTCAGGAGGACTGGGAGCTGGCGGAGCGACTCCAGCGGGAGGAAGAAGAGGCCTTTGCCAGCAGTCAGAGCAGCCAAGGGGCCCAATCCCTCACATTCTCCAAGTTTGAAGAAAAGAAAACCAACGAGAAGACCCGCAAGGTTACCACAGTGAAGAAATTCTTCAGTGCATCTTCCAGGGTCGGATCAAAGAAGGAAATTCAGGAAGCAAAAGCTCCCAGTCCTTCCATAAACCGGCAAACCAGCATTGAAACGGATAGAGTGTCTAAGGAGTTCATAGAATTTCTCAAGACCTTCCACAAGACAGGCCAAGAAATCTATAAACAGACCAAGCTGTTTTTGGAAGGAATGCATTACAAAAGGGATCTAAGCATTGAAGAACAGTCAGAGTGTGCTCAGGATTTCTACCACAATGTGGCCGAAAGGATGCAAACTCGTGGGAAAGTGCCTCCAGAAAGAGTCGAGAAGATAATGGATCAGATTGAAAAGTACATCATGACTCGTCTCTATAAATATGTATTCTGTCCAGAAACTACTGATGATGAGAAGAAAGATCTTGCCATTCAAAAGAGAATCAGAGCCCTGCGCTGGGTTACGCCTCAGATGCTGTGTGTCCCTGTTAATGAAGACATCCCAGAAGTGTCTGATATGGTGGTGAAGGCGATCACAGATATCATTGAAATGGATTCCAAGCGTGTGCCTCGAGACAAGCTGGCCTGCATCACCAAGTGCAGCAAGCACATCTTCAATGCCATCAAGATCACCAAGAATGAGCCGGCGTCAGCGGATGACTTCCTCCCCACCCTCATCTACATTGTTTTGAAGGGCAACCCCCCACGCCTTCAGTCTAATATCCAGTATATCACGCGCTTCTGCAATCCAAGCCGACTGATGACTGGAGAGGATGGCTACTATTTCACCAATCTGTGCTGTGCTGTGGCTTTCATTGAGAAGCTAGACGCCCAGTCTTTGAATCTAAGTCAGGAGGATTTTGATCGCTACATGTCTGGCCAGACCTCTCCCAGGAAGCAAGAAGCTGAGAGTTGGTCTCCTGATGCTTGCTTAGGCGTCAAGCAAATGTATAAGAACTTGGATCTCTTGTCTCAGTTGAATGAACGACAAGAAAGGATCATGAATGAAGCCAAGAAACTGGAAAAAGACCTCATAGATTGGACAGATGGAATTGCAAGAGAAGTTCAAGACATCGTTGAGAAATACCCACTGGAAATTAAGCCTCCGAATCAACCGTTAGCAGCTATTGACTCTGAAAACGTTGAAAATGATAAACTTCCTCCACCACTGCAACCTCAAGTTTATGCAGGATGA'
q9uj41_2_protein_seq = '''MSLKSERRGIHVDQSDLLCKKGCGYYGNPAWQGFCSKCWREEYHKARQKQIQEDWELAER
LQREEEEAFASSQSSQGAQSLTFSKFEEKKTNEKTRKVTTVKKFFSASSRVGSKKEIQEA
KAPSPSINRQTSIETDRVSKEFIEFLKTFHKTGQEIYKQTKLFLEGMHYKRDLSIEEQSE
CAQDFYHNVAERMQTRGKVPPERVEKIMDQIEKYIMTRLYKYVFCPETTDDEKKDLAIQK
RIRALRWVTPQMLCVPVNEDIPEVSDMVVKAITDIIEMDSKRVPRDKLACITKCSKHIFN
AIKITKNEPASADDFLPTLIYIVLKGNPPRLQSNIQYITRFCNPSRLMTGEDGYYFTNLC
CAVAFIEKLDAQSLNLSQEDFDRYMSGQTSPRKQEAESWSPDACLGVKQMYKNLDLLSQL
NERQERIMNEAKKLEKDLIDWTDGIAREVQDIVEKYPLEIKPPNQPLAAIDSENVENDKL
PPPLQPQVYAG'''.replace('\n', '')

hopd_ecolx_dna_seq='''ttaacg gcctgccagc
aaactatccc aacctacaac gaaacccgca gccgccagaa atggcccaaa cggcagcggg
ttttttaatg actcttttcc tctcatgagc agaccaacaa ctacggctcc gcaagcgaat
gaggcagcaa ggaaaaccag tcgtggcaag aaggtccagg tatgccaggc accaagggcg
gcgagaaact tcacatcacc atagcctaat ccttctttat gacgcagtat gcgatacccc
cagtaaatga cagcgaatgt gccgtaacca ataatcgccc cccataacgc atcggctaaa
cagtcaggat tacaaacctg agagaacaaa agacctgacc agagtaacgg gcaggtaaat
ctgtcaggta acaagccatg ctttgcatcc cataaaaaaa gtaggacaga taaacaagcg
tacaaaatca gaaagggaag ggtagcgtcc at'''.replace('\n', '').replace(' ', '').upper()
hopd_ecolx_protein_seq='''MDATLPFLILYACLSVLLFLWDAKHGLLPDRFTCPLLWSGLLFSQVCNPDCLADALWGAI
IGYGTFAVIYWGYRILRHKEGLGYGDVKFLAALGAWHTWTFLPRLVFLAASFACGAVVVG
LLMRGKESLKNPLPFGPFLAAAGFVVGWDSLLAGR'''.replace('\n', '')

eftu1_ecoli_protein_seq='''MSKEKFERTKPHVNVGTIGHVDHGKTTLTAAITTVLAKTYGGAA
RAFDQIDNAPEEKARGITINTSHVEYDTPTRHYAHVDCPGHADYVKNMITGAAQMDGA
ILVVAATDGPMPQTREHILLGRQVGVPYIIVFLNKCDMVDDEELLELVEMEVRELLSQ
YDFPGDDTPIVRGSALKALEGDAEWEAKILELAGFLDSYIPEPERAIDKPFLLPIEDV
FSISGRGTVVTGRVERGIIKVGEEVEIVGIKETQKSTCTGVEMFRKLLDEGRAGENVG
VLLRGIKREEIERGQVLAKPGTIKPHTKFESEVYILSKDEGGRHTPFFKGYRPQFYFR
TTDVTGTIELPEGVEMVMPGDNIKMVVTLIHPIAMDDGLRFAIREGGRTVGAGVVAKV
LG'''.replace(' ', '').replace('\n', '')
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
gctaa'''.replace('\n', '').replace(' ', '').upper()

rpsL_ecoli_protein_seq='''MATVNQLVRKPRARKVAKSNVPALEACPQKRGVCTRVYTTTPKK
PNSALRKVCRVRLTNGFEVTSYIGGEGHNLQEHSVILIRGGRVKDLPGVRYHTVRGAL
DCSGVKDRKQARSKYGVKRPKA'''.replace(' ', '').replace('\n', '')
rpsL_ecoli_dna_seq='''atggcaac agttaaccag ctggtacgca aaccacgtgc tcgcaaagtt gcgaaaagca
acgtgcctgc gctggaagca tgcccgcaaa aacgtggcgt atgtactcgt gtatatacta
ccactcctaa aaaaccgaac tccgcgctgc gtaaagtatg ccgtgttcgt ctgactaacg
gtttcgaagt gacttcctac atcggtggtg aaggtcacaa cctgcaggag cactccgtga
tcctgatccg tggcggtcgt gttaaagacc tcccgggtgt tcgttaccac accgtacgtg
gtgcgcttga ctgctccggc gttaaagacc gtaagcaggc tcgttccaag tatggcgtga
agcgtcctaa ggcttaa'''.replace('\n', '').replace(' ', '').upper()


if __name__ == '__main__':
    unittest.main()