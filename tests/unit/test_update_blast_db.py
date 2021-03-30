import unittest

from subprocess import Popen
from mock import Mock
from pprint import pprint

from ccd.update_blast_db import Entry

class test_Entry(unittest.TestCase):
    
    def test_equality(self):
        a = Entry(('2C7N_A', 'ABC')) 
        b = Entry(('2C7N_B', 'ABC')) #different chain id
        c = Entry(('2C7N_A', 'XYZ')) # different sequence
        e = Entry(('2C7M_A', 'XYZ')) # different PDB_ID
        
        self.assertTrue(a == b)
        self.assertFalse(a == c)
        self.assertFalse(a == e)
        
    def test_set_operation(self):
        from operator import xor
        a = Entry(('2C7N_A', 'ABC'))
        b = Entry(('2C7N_B', 'ABC')) #different chain id
        c = Entry(('2C7N_A', 'XYZ')) # different sequence
        d = Entry(('2C7M_A', 'XYZ')) # different PDB_ID
        res = set([a,b,c,d])
        #note that since in compares on hash, a in res and b in res are both true
        # even if only either a or b are actually present, and a is not b.
        # I don't think this will cause problems, therefore I leave the loophole
        #open and move on 
        self.assertTrue(a in res or b in res) 
        self.assertTrue(c in res)
        self.assertTrue(d in res)
        
    def test_hash(self):
        a = Entry(('2C7N_A', 'ABC'))
        b = Entry(('2C7N_B', 'ABC')) #different chain id
        c = Entry(('2C7N_A', 'XYZ')) # different sequence
        d = Entry(('2C7M_A', 'XYZ')) # different PDB_ID
        self.assertEqual(a.__hash__(), b.__hash__())
        self.assertNotEqual(a.__hash__(), c.__hash__())
        self.assertNotEqual(a.__hash__(), d.__hash__())

    def test_real_entries(self):
        a = Entry(('2c7n_K_70174',
                   'MSLKSERRGIHVDQSDLLCKKGCGYYGNPAWQGFCSKCWREEYHKARQKQIQEDWELAERLQREEEEAFASSQS'))
        a2 = Entry(('2c7n_L_70174',
                   'MSLKSERRGIHVDQSDLLCKKGCGYYGNPAWQGFCSKCWREEYHKARQKQIQEDWELAERLQREEEEAFASSQS'))
        b = Entry(('2c7n_B_70175',
                   'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'))
        c = Entry(('2c7m_A_70167',
                   'MSLKSERRGIHVDQSDLLCKKGCGYYGNPAWQGFCSKCWREEYHKARQKQIQEDWELAERLQREEEEAFASSQS'))
        d = Entry(('2c7m_B_70168',
                   'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'))
        self.assertFalse(a is b)
        self.assertFalse(a is a2)
        self.assertFalse(c is d)
        self.assertFalse(a == b) #differ in sequence
        self.assertTrue(a == a2) #same PDB ID, same sequence, only different chain
        self.assertFalse(a == c) #same sequence, different PDB id
        
        