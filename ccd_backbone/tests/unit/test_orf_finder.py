import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

import ccd.orf_finder as of

class TestTrimFrame(unittest.TestCase):
    
    def test_trim_frame_fwd_str_object(self):
        self.assertEqual(of.trim_frame('ATGcccG', mode='fwd'), 'ATGccc')
        self.assertEqual(of.trim_frame('ATGcccGC', mode='fwd'), 'ATGccc')
        self.assertEqual(of.trim_frame('ATGcccGCT', mode='fwd'), 'ATGcccGCT')
    
    def test_trim_frame_rev_str_object(self):
        self.assertEqual(of.trim_frame('aATGcccGCT', mode='rev'), 'ATGcccGCT') 
        self.assertEqual(of.trim_frame('aaATGcccGCT', mode='rev'), 'ATGcccGCT')
        self.assertEqual(of.trim_frame('aaaATGcccGCT', mode='rev'), 'aaaATGcccGCT')
        
    def test_trim_frame_fwd_Seq_object(self):
        a = SeqRecord(Seq('ATGcccG', IUPAC.IUPACUnambiguousDNA))
        self.assertEqual('ATGccc', 
                         of.trim_frame(a, mode='fwd').seq)
        b = SeqRecord(Seq('ATGcccGC', IUPAC.IUPACUnambiguousDNA))
        self.assertEqual('ATGccc', 
                         of.trim_frame(b, mode='fwd').seq)
        c = SeqRecord(Seq('ATGcccGCT', IUPAC.IUPACUnambiguousDNA))
        self.assertEqual('ATGcccGCT', 
                         of.trim_frame(c, mode='fwd').seq)
    def test_trim_frame_rev_Seq_object(self):
        a = SeqRecord(Seq('aATGcccGCT', IUPAC.IUPACUnambiguousDNA))
        self.assertEqual('ATGcccGCT', 
                         of.trim_frame(a, mode='rev').seq)
        b = SeqRecord(Seq('aaATGcccGCT', IUPAC.IUPACUnambiguousDNA))
        self.assertEqual('ATGcccGCT', 
                         of.trim_frame(b, mode='rev').seq)
        c = SeqRecord(Seq('aaaATGcccGCT', IUPAC.IUPACUnambiguousDNA))
        self.assertEqual('aaaATGcccGCT', 
                         of.trim_frame(c, mode='rev').seq)
        
class TestExtendFunctions(unittest.TestCase):
     
    def extend_to_Ct(self):
        a = SeqRecord(Seq('aaacccATGATGTAA', IUPAC.IUPACAmbiguousDNA()))
        self.assertEqual(of.extend_to_Ct(a, 6), 'MM')
    
    def extend_to_Nt(self):
        a = SeqRecord(Seq('TAAATGCCCCAT', IUPAC.IUPACUnambiguousDNA()))
        self.assertEqual(of.extend_to_Nt(a, 9), 'MP')
        
class TestCleavage(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        super(TestCleavage, cls).setUpClass()
        cls.matches = {'TEV':[6, 13],
                        '3C':[6],
                        'Thrombin':[4],
                        'Senp2': [4, 9]}
        cls.ORFs = {}
    
    def test_TEV_cleavage(self):    
        TEV_seq = Seq('ENLYFQGENLYFQST', IUPAC.IUPACProtein())
        self.ORFs['TEV'] = TEV_seq[-2:]
        sites = of.find_cleavage_sites(TEV_seq)
        cleaved = of.perform_cleavage(TEV_seq, sites) 
        self.assertEqual(self.matches['TEV'], sites['TEV'])
        self.assertEqual([str(i) for i in self.ORFs['TEV']], 
                         [str(i) for i in cleaved['TEV']])
        
    def test_3C_cleavage(self):
        Psc_seq = Seq('LEVLFQGPT', IUPAC.IUPACProtein())
        self.ORFs['3C'] = Psc_seq[6:]
        sites = of.find_cleavage_sites(Psc_seq)
        cleaved = of.perform_cleavage(Psc_seq, sites)
        self.assertEqual(self.matches['3C'], sites['3C'])
        self.assertEqual([str(i) for i in self.ORFs['3C']], 
                         [str(i) for i in cleaved['3C']])   

    def test_thrombin_cleavage(self):
        thrombin_seq = Seq('LVPRGST', IUPAC.IUPACProtein())
        self.ORFs['Thrombin'] = thrombin_seq[4:]
        sites = of.find_cleavage_sites(thrombin_seq)
        cleaved = of.perform_cleavage(thrombin_seq, sites)
        self.assertEqual(self.matches['Thrombin'], sites['Thrombin'])
        self.assertEqual([str(i) for i in self.ORFs['Thrombin']], 
                         [str(i) for i in cleaved['Thrombin']])
        
    def test_SUMO_cleavage(self):
        SUMO_seq = Seq('QQTGGEQTGGT', IUPAC.IUPACProtein())
        self.ORFs['Senp2'] = SUMO_seq[9:]
        sites = of.find_cleavage_sites(SUMO_seq)
        cleaved = of.perform_cleavage(SUMO_seq, sites)
        self.assertEqual(self.matches['Senp2'], sites['Senp2'])
        self.assertEqual([str(i) for i in self.ORFs['Senp2']], 
                         [str(i) for i in cleaved['Senp2']])
        
class Test_build_cleaved_orfs(unittest.TestCase):
    
    def test_a_bunch(self):
        Nt_tags = {'TEV' : 'GG', '3C' : 'PP', 'Thrombin' : 'NN', 'SUMO' : 'SS'}
        Ct_tags =  {'TEV' : 'CGG', '3C' : 'CPP', 'Thrombin' : 'CNN', 'SUMO' : 'CSS'}
        orf = 'AAAA'
        expected =  {'TEV' : 'GGAAAACGG', '3C' : 'PPAAAACPP',
                   'Thrombin' : 'NNAAAACNN', 'SUMO' : 'SSAAAACSS'}
        result = of.build_cleaved_orfs(orf, Nt_tags, Ct_tags)
        self.assertEqual(result, expected)
    
if __name__ == '__main__':
    unittest.main()