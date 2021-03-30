import unittest

from mock import MagicMock
from pprint import pprint

from ccd.homolog_finder import Evolutionator
from ccd.dna_finder_async import DatabaseCrawler

class test_homolog_finder_black_box(unittest.TestCase):
    
    def test_start_from_id(self):
        e = Evolutionator()
        e.fill_data('q9uj41', q9uj41_seq, '')
        e.query_picr_by_sequence = MagicMock()
        e.rename_sequences = MagicMock(return_value = []) #saves test time
        e.search()
        e.query_picr_by_sequence.assert_not_called()
        self.assertEqual('ENSG00000154710', e.ensembl_id)
         
    def test_start_from_sequence(self):
        e = Evolutionator()
        e.fill_data('', q9uj41_seq, '')
        e.rename_sequences = MagicMock(return_value = []) #saves test time
        e.search()
        self.assertEqual('ENSG00000154710', e.ensembl_id)
    
    def test_ambiguous_mapping_from_sequence(self):
        e = Evolutionator()
        e.fill_data('', pcna_seq, '')
        e.rename_sequences = MagicMock(return_value = []) #saves test time
        e.search()
        self.assertTrue(e.ensembl_id in ['ENSMFAG00000036845',
                                         'ENSPTRG00000013222']) #they are actually identical and returned at random by ENSEMBL
        
    
    def test_ambiguous_mapping_from_id(self):
        e = Evolutionator()
        e.fill_data('P61258', pcna_seq, '')
#         e.extensive_map = MagicMock(return_value='ENSG00000132646')
        e.rename_sequences = MagicMock(return_value = []) #saves test time
        e.search()
        self.assertEqual('ENSMFAG00000036845', e.ensembl_id)
        
    def test_more_sequences(self):
        e = Evolutionator()
        e.rename_sequences = MagicMock(return_value = [])
        uniprot_ids = ['G3VK03', 'SAP18_HUMAN', 'RNPS1_HUMAN', 'G3WX24_SARHA']
        ensembl_ids = ['ENSSHAG00000003081', 'ENSG00000150459',
                       'ENSG00000205937', 'ENSSHAG00000016951']
         
        for index, uid in enumerate(uniprot_ids):
            finder = DatabaseCrawler(uid)
            protein = finder.fetch_uniprot_entry(uid)
            protein = finder.parse_uniprot_data(protein)
            seq = protein.isoforms[0].protein_seq
            e.fill_data(uid, seq, '')
            e.search()
            self.assertEqual(e.ensembl_id, ensembl_ids[index])
            
    def test_get_ensembl_alignment(self):
        e = Evolutionator()
        query = 'ENSG00000105393'
        e.fill_data('P61258', pcna_seq, '')
        res = e.get_ensembl_alignment('ENSG00000105393')
        res = res['data'][0]['homologies'] #ENSEMBL has some interesting JSON schemes....
        ids = [i['target']['id'] for i in res]
        self.assertIn('ENSPTRG00000010664', ids)
        
q9uj41_seq = """
MVVVTGREPDSRRQDGAMSSSDAEDDFLEPATPTATQAGHALPLLPQERCAEFPALRGPP
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
GIAREVQDIVEKYPLEIKPPNQPLAAIDSENVENDKLPPPLQPQVYAG
""".replace('\n', '')
pcna_seq = """
MFEARLVQGSILKKVLEALKDLINEACWDISSSGVNLQSMDSSHVSLVQLTLRSEGFDTY
RCDRNLAMGVNLTSMSKILKCAGNEDIITLRAEDNADTLALVFEAPNQEKVSDYEMKLMD
LDVEQLGIPEQEYSCVVKMPSGEFARICRDLSHIGDAVVISCAKDGVKFSASGELGNGNI
KLSQTSNVDKEEEAVTIEMNEPVQLTFALRYLNFFTKATPLSSTVTLSMSADVPLVVEYK
IADMGHLKYYLAPKIEDEEGS
""".replace('\n', '')

if __name__ == '__main__':
    unittest.main()