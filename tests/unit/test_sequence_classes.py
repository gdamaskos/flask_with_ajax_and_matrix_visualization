import unittest

import ccd.sequence_classes as seq

class TestPdbCrossref(unittest.TestCase):
    
    def test_lt_gt(self):
        a = seq.PdbCrossref('a','b',1,5,['A','B'])
        b = seq.PdbCrossref('a','b',1,5,['A'])
        self.assertTrue(b<a)
        self.assertTrue(a>b)
        b = seq.PdbCrossref('a','b',1,4,['A','B'])
        self.assertTrue(b<a)
        self.assertTrue(a>b)
        b = seq.PdbCrossref('a','b',0,5,['A','B'])
        self.assertTrue(b<a)
        self.assertTrue(a>b)
        b = seq.PdbCrossref('a','a',1,5,['A','B'])
        self.assertTrue(b<a)
        self.assertTrue(a>b)
        
    def test_eq(self):
        #parent_entry_id does not count in equality comparisons 
        a = seq.PdbCrossref('a','b',1,5,['A','B'])
        b = seq.PdbCrossref('b','b',1,5,['A','B'])
        self.assertTrue(a==b)