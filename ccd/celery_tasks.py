from ccd import celery

import logging
from urllib.parse import urlparse
from http.client import HTTPConnection 
from bs4 import BeautifulSoup
from Bio.SeqUtils.ProtParam import ProteinAnalysis as PA

_log = logging.getLogger(__name__)


@celery.task
def crysol_service(service_name, pp_seqs):
    """Takes as arguments the service name and the amino acid sequence for predictions input.
    Uses the pythonic way to create a programmatic switch between the services.
    A dictionary is created with the function names and then the corresponding function is called,
    with the amino acids as input string for predictions.
    """
    options = {
               'SECRET': run_secret,
               'PROSO2': run_proso2,
               'RPSP': run_rpsp,
               }
    return options[service_name](pp_seqs)


def run_rpsp(pp_seqs):
    """Makes several web requests to the remote predictions server.
    Processes the reply (HTML) with beautiful soup and extracts the predictions list. 
    Returns the predictions in the end.
    """

    try:
        rpsp_predictions = []
        for pp_seq in pp_seqs:
            conn = HTTPConnection('www.biotech.ou.edu', timeout = 180)
            params = PA(pp_seq) # works with string or Seq objects
            MW = str(round(params.molecular_weight(), 1)) # in Dalton, rounded to 1 decimal
            pI = str(round(params.isoelectric_point(), 1))
            conn.request("GET", '/cgi-bin/solcalc4.cgi?AvgPI=' + pI + '&MW=' + MW + '&amino=' + pp_seq + '&Submit=Submit')
            body = conn.getresponse().read()
            soup = BeautifulSoup(body, 'html.parser')
            prediction = soup.find_all('font')[1].getText()[:-9]
            rpsp_predictions.append(prediction)
        return rpsp_predictions
    except Exception as e:
        error_string = '{} {}'.format('RPSP', str(e))
        _log.error(error_string)
        return 'ERROR'


def run_proso2(pp_seqs):
    """Makes several web requests to the remote predictions server.
    Processes the reply (HTML) with beautiful soup and extracts the predictions list. 
    Returns the predictions in the end.
    """
 
    try:   
        proso2_predictions = []
        url = "http://mbiljj45.bio.med.uni-muenchen.de:8888/prosoII/prosoII.seam"
        urlparts = urlparse(url)
        conn = HTTPConnection(urlparts.netloc)
        conn.request("GET", urlparts.path)
        body = conn.getresponse().read()
        
        soup = BeautifulSoup(body, 'html.parser')
        action = soup.find('form', id='form1').get('action')
        view_state = soup.find('input', id='javax.faces.ViewState').get('value')
        
        for pp_seq in pp_seqs:
            conn = HTTPConnection('mbiljj45.bio.med.uni-muenchen.de:8888', timeout = 180)
            conn.request("GET", action + '?javax.portlet.faces.DirectLink=true&AJAXREQUEST=_viewRoot&form1=form1&form1%3AinText=' + pp_seq + '&form1_link_hidden_=form1%3Aj_id28&javax.faces.ViewState=' + view_state + '&form1%3Aj_id27=form1%3Aj_id27')
            body = conn.getresponse().read()
            soup = BeautifulSoup(body, 'html.parser')
            prediction = soup.find_all('td')[2].getText()[-5:]
            proso2_predictions.append(prediction)
        return proso2_predictions
    except Exception as e:
        error_string = '{} {}'.format('PROSOII', str(e))
        _log.error(error_string)
        return 'ERROR'


def run_secret(pp_seqs):
    """Makes several web requests to the remote predictions server.
    Processes the reply (HTML) with beautiful soup and extracts the predictions list. 
    Returns the predictions in the end.
    """

    try:    
        secret_predictions = []
        url = "http://mbiljj45.bio.med.uni-muenchen.de:8888/secret/secret.seam"
        urlparts = urlparse(url)
        conn = HTTPConnection(urlparts.netloc)
        conn.request("GET", urlparts.path)
        body = conn.getresponse().read()
        
        soup = BeautifulSoup(body, 'html.parser')
        action = soup.find('form', id='form1').get('action')
        view_state = soup.find('input', id='javax.faces.ViewState').get('value')
        # parallelization can take place here
        for pp_seq in pp_seqs:
            conn = HTTPConnection('mbiljj45.bio.med.uni-muenchen.de:8888', timeout = 180)
            conn.request("GET", action + '?javax.portlet.faces.DirectLink=true&AJAXREQUEST=_viewRoot&form1=form1&form1%3AinText=' + pp_seq + '&form1_link_hidden_=form1%3Aj_id17&javax.faces.ViewState=' + view_state + '&form1%3Aj_id16=form1%3Aj_id16')
            body = conn.getresponse().read()
            soup = BeautifulSoup(body, 'html.parser')
            prediction = soup.find_all('td')[2].getText()[4:]
            secret_predictions.append(prediction)
        return secret_predictions
    except Exception as e:
        error_string = '{} {}'.format('SECRET', str(e))
        _log.error(error_string)
        return 'ERROR'
