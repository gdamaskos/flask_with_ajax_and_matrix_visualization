import json
import os
import unittest

import ccd.vector_sequence as vs

class testVectorContainer(unittest.TestCase):
    
    #testing load_vectors_from_file 
    def test_load_vectors_from_file(self):
        f = os.path.abspath('static/test_vectors.json')
        c = vs.VectorLoader()
        res = c.load_vectors_from_file(f)
        res_tags = sorted([r['short_name'] for r in res])
        exp_tags = sorted(['NKI_1_1', 'NKI_1_2', 'NKI_1_13', 'NKI_3_11'])
        self.assertEqual(res_tags, exp_tags, 'Vectors were not properly loaded')
        
    def test_load_empty_vector_file(self):
        c = vs.VectorLoader()
        #NamedTempFile does not work on windows...
        res = c.load_vectors_from_file(
            os.path.abspath('static/empty_vector_file.json')) 
        self.assertEqual(res, [{}]) # list with empty dictionary
        
    def test_custom_vector_file_absent(self):
        #should proceed as usual
        c = vs.VectorLoader()
        c.CUSTOM_VECTOR_FILE = 'notafile'
        a = c.make_vectors()
        self.assertTrue('NKI_1_1' in a.keys())
        
    def test_custom_vector_syntax_error(self):
        c = vs.VectorLoader()
        #NamedTempFile does not work on windows...
        with self.assertRaises(json.decoder.JSONDecodeError):
            c.load_vectors_from_file( 
                os.path.abspath('static/test_vectors_invalid.json')) 
            
    def test_make_vectors(self):
        c = vs.VectorLoader()
        a = c.make_vectors()
        #test that two random vectors have been loaded correctly
        self.assertTrue('NKI_1_1' in a.keys())
        self.assertTrue(hasattr(a['NKI_1_1'], 'seq'))
        self.assertTrue(hasattr(a['NKI_1_1'], 'Nt_tag'))
        self.assertTrue(hasattr(a['NKI_1_1'], 'Ct_tag'))
        self.assertEqual(a['NKI_1_1'].Nt_tag, True)
        self.assertEqual(a['NKI_1_1'].Ct_tag, False)
        self.assertTrue('NKI_1_13' in a.keys())
        self.assertTrue(hasattr(a['NKI_1_13'], 'Nt_tag'))
        self.assertTrue(hasattr(a['NKI_1_13'], 'Ct_tag'))
        self.assertEqual(a['NKI_1_13'].Nt_tag, True)
        self.assertEqual(a['NKI_1_13'].Ct_tag, True)
    
class testBaseVector(unittest.TestCase):
    
    def test_ligate_LIC(self):
        pass

class testCtTagVector(unittest.TestCase):
    
    def test_ligate_LIC(self):
        pass