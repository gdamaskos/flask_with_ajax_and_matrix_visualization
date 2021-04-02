import unittest

from subprocess import Popen
from mock import Mock
from pprint import pprint

from ccd.pdb_crawler import (PdbCrawler, BlastFailError, BlastNoResult)

class testPdbCrawlerOffline(unittest.TestCase):
    
    def __init__(self, *args, **kwargs):
        super(testPdbCrawlerOffline, self).__init__(*args, **kwargs)
        
    def test_checks_for_protein_sequence(self):
        with self.assertRaises(ValueError):
            p = PdbCrawler(id_='Q9uj41')
            p.run(online=False)
    
    def test_raises_on_generic_error(self):
        try:
            x = Popen.communicate
            p = PdbCrawler(protein_seq='NTASEQ')
            Popen.communicate = Mock(return_value=('Error message', ''))
            with self.assertRaises(BlastFailError):
                p.run(online=False)
        except:
            raise
        finally: #otherwise the following tests will fail horribly
            Popen.communicate = x
            
    def test_raises_on_blank_result(self):
        seq = 'NTASEQ'
        p = PdbCrawler(protein_seq=seq)
        with self.assertRaises(BlastNoResult):
            p.run(online=False)
        
    def test_correct_run(self):
        self.maxDiff = None
        p = PdbCrawler(protein_seq=rnps1_seq)
        res = p.run(online=False)
        exp = ('------------------------------------------------------------------------------------------------------------------------------------------------------------>--------------------------------------------------------------------------------------<-------------------------------------------------------------',
               '------------------------------------------------------------------------------------------------------------------------------------------------------------1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111-------------------------------------------------------------')
        assert res[0][0][156] == '>'
        assert res[0][0][243] == '<'
        assert len(res[0][0]) == len(rnps1_seq)
        assert len(res[0][0]) == len(res[0][1])
        assert set(res[0][1][156:244]) == {'1'} #unless someone solves another structure ;-)
        self.assertEqual(exp, res[0])
        
#test stuff#
rnps1_seq = '''
    MDLSGVKKKSLLGVKENNKKSSTRAPSPTKRKDRSDEKSKDRSKDKGATKESSEKDRGRD
    KTRKRRSASSGSSSTRSRSSSTSSSGSSTSTGSSSGSSSSSASSRSGSSSTSRSSSSSSS
    SGSPSPSRRRHDNRRRSRSKSKPPKRDEKERKRRSPSPKPTKVHIGRLTRNVTKDHIMEI
    FSTYGKIKMIDMPVERMHPHLSKGYAYVEFENPDEAEKALKHMDGGQIDGQEITATAVLA
    PWPRPPPRRFSPPRRMLPPPPMWRRSPPRMRRRSRSPRRRSPVRRRSRSPGRRRHRSRSS
    SNSSR
    '''.replace('\n', '').replace('\t','').replace(' ','')
rnps1_rrm = '''
    SPKPTKVHIGRLTRNVTKDHIMEIFSTYGKIKMIDMPVERMHPHLSKGYAYVEFENPDEAE
    KALKHMDGGQIDGQEITATAVLAPWPR
    '''.replace('\n', '').replace('\t','').replace(' ','')
        
