# flask_testing/test_base.py
import json
import sys
import io
import logging
import unittest

from mock import MagicMock
from flask import url_for
from flask_testing import TestCase as ft_TestCase

from ccd.dna_finder import DatabaseCrawler

class BaseTestCase(ft_TestCase):
    
    #Standard flask-testing procedure
    def create_app(self):
        from ccd import app
        self.app = app
        return app

class TestExceptionsAreReportedToFrontend(BaseTestCase):
    
    def setUp(self):
        self.silence_sublogger = True
        if self.silence_sublogger:
            self.old_out = sys.stdout
            sys.stdout = io.StringIO()
        
    def tearDown(self):
        if self.silence_sublogger:
            sys.stdout = self.old_out
      
    def test_IdNotFoundError(self):
        uniprot_id = 'notAnEntry'
        url = url_for('search_uniprot', 
                      uniprot_accession_or_entry_name=uniprot_id)
        exp_response = json.loads(json.dumps({"reason": "ID not found.",
                                              "status": "ERROR"}))
        response = json.loads(self.client.get(url).data)
        self.assertEqual(response, exp_response)
        
    def test_IGiveUpError(self):
        uniprot_id = 'H2SBB3_TAKRU'
        url = url_for('search_uniprot', 
                      uniprot_accession_or_entry_name=uniprot_id)
        exp_response = json.loads(json.dumps({'status': 'ERROR',
                                              'reason': 'Unable to map entry to DNA.'}))
        response = json.loads(self.client.get(url).data)
        self.assertEqual(response, exp_response)
        
    def test_503(self):
        import requests
        uniprot_id = 'H2SBB3_TAKRU'
        url = url_for('search_uniprot', 
                      uniprot_accession_or_entry_name=uniprot_id)
        DatabaseCrawler.search = MagicMock()
        DatabaseCrawler.search.side_effect = requests.get('http://httpbin.org/status/503')
        response = json.loads(self.client.get(url).data)
        exp = json.loads(json.dumps({'status': 'ERROR', 
               'reason': 'Could not complete the multiple alignment. Remote server error.'}))
        self.assertEqual(response, exp)
        
        
    uid = 'H2SBB3'

#following code adapted from https://eli.thegreenplace.net/2011/08/02/python-unit-testing-parametrized-test-cases
#it's a quick factory class for repeated tests with parameter sets
#e.g. to quickly test multiple entry - response pairs  
#see main() for usage
class ParametrizedTestCase(unittest.TestCase):
    """ TestCase classes that want to be parametrized should
        inherit from this class.
    """
    def __init__(self, methodName='runTest', param=None):
        super(ParametrizedTestCase, self).__init__(methodName)
        self.param = param

    @staticmethod
    def parametrize(testcase_klass, param=None):
        """ Create a suite containing all tests taken from the given
            subclass, passing them the parameter 'param'.
        """
        testloader = unittest.TestLoader()
        testnames = testloader.getTestCaseNames(testcase_klass)
        suite = unittest.TestSuite()
        for name in testnames:
            suite.addTest(testcase_klass(name, param=param))
        return suite
    
class TestNormalOperation(ParametrizedTestCase, BaseTestCase):
    
    def setUp(self):
        self.logger = logging.getLogger('self')
        submodule_logger = logging.getLogger('ccd.dna_finder')
        submodule_logger.level = logging.ERROR
            
    def test_protein(self):
        uniprot_id, organism, exp_isoforms = self.param
#         organism = self.param[1]
        msg = 'Testing {}'.format(uniprot_id)
        self.logger.info(msg)
        url = url_for('search_uniprot', 
                      uniprot_accession_or_entry_name=uniprot_id)
        response = json.loads(self.client.get(url).data)
        self.assertEqual(len(response['protein_isoforms']),
                         exp_isoforms) 
        self.assertEqual(response['protein_source_organism'], organism)

def main(test_set):
            
    suite = unittest.TestSuite()
    for t in test_set:
        suite.addTest(ParametrizedTestCase.parametrize(TestNormalOperation, param=t))
    unittest.TextTestRunner(verbosity=1).run(suite)


if __name__ == '__main__':
    test_set = [('q9uj41', 'Homo sapiens', 4),
               ('PCNA_HUMAN', 'Homo sapiens', 1),
               ('P02990', 'Escherichia coli', 1)]
    main()


