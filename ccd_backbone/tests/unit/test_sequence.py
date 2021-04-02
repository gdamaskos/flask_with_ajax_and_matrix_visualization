import unittest

from Bio.Seq import Seq
from Bio.Alphabet import HasStopCodon, IUPAC


from ccd.sequence import translate, build_primers_tm, build_primers_bs, dna_validator, over_validator

class test_stuff(unittest.TestCase):
    
    def test_dna_validator(self):
        """Test that the dna validator works."""
        
        dna          = 'TGGAGACGGAAACATCCGAGGACATCCGGAGGAACCCGGGGAGTTCTGAGTGGTAAT'
        # test invalid characters
        invalid_dna1 = 'EETGGAGACGGAAACASTCCGAGGACATCCGGAGGAACCCGGGGAGTZVTHHCTGAGTGGTAAT'
        # test invalid length
        invalid_dna2 = 'GGAGACGGAAACATCCGAGGACATCCGGAGGAACCCGGGGAGTTCTGAGTGGTAAT'
        # test for invalid internal stop
        invalid_dna3 = 'TGGAGACGGAAACATCCGAGGACATCCGGAGGAACCCGGGGAGTTTGAGTGGTAATC'
        expected_validationT = True
        expected_validationF = False
        result_validation1 = dna_validator(dna)
        self.assertEqual(result_validation1, expected_validationT)
        result_validation2 = dna_validator(invalid_dna1)
        self.assertEqual(result_validation2, expected_validationF)
        result_validation3 = dna_validator(invalid_dna2)
        self.assertEqual(result_validation3, expected_validationF)
        result_validation4 = dna_validator(invalid_dna3)
        self.assertEqual(result_validation4, expected_validationF)
        
        
    def test_over_validator(self):
        """Test that the overhang validator works. It should accept any letter."""
        overhang = "awscxzfrtyuiopklkjhgfdswqcvbnmklGHYTTFDESVNJIJ"
        invalid_overhang = "123gggcccccaaaatttttTTTAAGGGCCCC"
        expected_validationT = True
        expected_validationF = False
        result_validation1 = over_validator(overhang)
        self.assertEqual(result_validation1, expected_validationT)
        result_validation2 = over_validator(invalid_overhang)
        self.assertEqual(result_validation2, expected_validationF)
    
    def test_translate(self):
        """Test that a dna open reading frame is translated correctly."""
        
        orf = 'ATGTGGAGACGGAAACATCCGAGGACATCCGGAGGAACCCGGGGAGTTCTGAGTGGTAATTAG'
        expected_primers = Seq('MWRRKHPRTSGGTRGVLSGN*', HasStopCodon(IUPAC.ExtendedIUPACProtein(), '*'))
        result_primers = translate(orf)
        self.assertEqual(result_primers, expected_primers)
        self.assertEqual(len(result_primers), 21)
        self.assertEqual(isinstance(result_primers, Seq), True)
    
    
    def test_build_primers_tm(self):
        """Test that primers are built correctly with melting temperature method.
        Also test that the index errors are correctly raised and the valid starts and stops are returned."""
    
        starts = [5, 9, 62]
        stops = [2, 13, 17]
        forward_overhang = 'cagggacccggt'
        reverse_overhang = 'cgaggagaagcccggtta'
        dna_orf = 'ATGTGGAGACGGAAACATCCGAGGACATCCGGAGGAACCCGGGGAGTTCTGAGTGGTAATTAG'
        expected_primers = [['Fw_5', 'cagggacccggtAAACATCCGAGGACATCCGGAGGAACCCG'],
                            ['Fw_9', 'cagggacccggtACATCCGGAGGAACCCGGGGAGTTCTG'],
                            ['Rv_13', 'cgaggagaagcccggttaGGTTCCTCCGGATGTCCTCGGATGTTTCC'],
                            ['Rv_17', 'cgaggagaagcccggttaCAGAACTCCCCGGGTTCCTCCGGATG']]
        expected_errors = [['Fw_62', ' Not enough bases in that direction'], ['Rv_2', ' Not enough bases in that direction']]
        expected_vstarts = [5, 9]
        expected_vstops = [13, 17]
        result_primers, result_errors, result_vstarts, result_vstops = build_primers_tm(dna_orf, starts, stops,
                                       forward_overhang, reverse_overhang, 65)
        self.assertEqual(result_primers, expected_primers)
        self.assertEqual(result_errors, expected_errors)
        self.assertEqual(expected_vstarts, result_vstarts)
        self.assertEqual(expected_vstops, result_vstops)
        
        
    def test_build_primers_bs(self):
        """Test that primers are built correctly with base length method.
        Also test that the index errors are correctly raised and the valid starts and stops are returned."""
    
        starts = [5, 9, 62]
        stops = [2, 13, 17]
        forward_overhang = 'cagggacccggt'
        reverse_overhang = 'cgaggagaagcccggtta'
        dna_orf = 'ATGTGGAGACGGAAACATCCGAGGACATCCGGAGGAACCCGGGGAGTTCTGAGTGGTAATTAG'
        expected_primers = [['Fw_5', 'cagggacccggtAAACATCCGA'],
                            ['Fw_9', 'cagggacccggtACATCCGGAG'],
                            ['Rv_13', 'cgaggagaagcccggttaGGTTCCTCCG'],
                            ['Rv_17', 'cgaggagaagcccggttaCAGAACTCCC']]
        expected_errors = [['Fw_62', ' Not enough bases in that direction'], ['Rv_2', ' Not enough bases in that direction']]
        expected_vstarts = [5, 9]
        expected_vstops = [13, 17]
        result_primers, result_errors, result_vstarts, result_vstops = build_primers_bs(dna_orf, starts, stops,
                                       forward_overhang, reverse_overhang, 10)
        self.assertEqual(result_primers, expected_primers)
        self.assertEqual(result_errors, expected_errors)
        self.assertEqual(expected_vstarts, result_vstarts)
        self.assertEqual(expected_vstops, result_vstops)
        
if __name__ == '__main__':
    unittest.main()