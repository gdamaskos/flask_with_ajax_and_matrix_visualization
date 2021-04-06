import logging
import os
import time
import json
import re
import math
import requests
import tempfile

from requests.compat import urljoin
from subprocess import Popen, PIPE
from urllib.parse import urlparse
from http.client import HTTPSConnection, HTTPConnection 
from bs4 import BeautifulSoup


# for the log variables to be visible
from ccd import app
from ccd.pdb_crawler import PdbCrawler
from ccd.utils import env_var
from ccd.pdb_crawler import BlastNoResult, BlastFailError, BlastTimeoutError
# from Bio.Blast.Record import Blast
# from errno import ECHRNG

_log = logging.getLogger(__name__)

crawler = PdbCrawler(protein_seq='TEMP') #global object to deal with idiosncracy of frontend - see comment in service()

class CoilsError(RuntimeError):
    pass

class IupredError(RuntimeError):
    pass

class PredatorError(RuntimeError):
    pass


def run_coils(aa_seq):
    """Creates a virtual temp file, adds a fasta header and the amino acids.
    Then feeds this file to the executable and gets the output.
    It checks for errors from the executable's run and raises an error if there is any.
    If not the case, it processes and returns the predictions.
    """
    aa_fasta = '>CCD_protein\n{}'.format(aa_seq)
    with tempfile.NamedTemporaryFile() as th:
        th.write(aa_fasta.encode('utf-8'))
        th.seek(0)
        _log.info('Running ncoils..')
        coils_pipe = Popen(app.config['NCOILS'], stdout=PIPE, 
                           stderr=PIPE, stdin=th)
        coils_output, coils_error = coils_pipe.communicate()
        if coils_error:
            raise CoilsError(coils_error)
        coils_predictions = ''
        for line in coils_output.splitlines():
            current_value = float(line.split()[4])
            norm_value = int(round(current_value*10))
            if norm_value == 10:
                out_value = 0
            elif norm_value <= 4:
                out_value = '-'
            else:
                out_value = norm_value
            coils_predictions = '{}{}'.format(coils_predictions, 
                                              str(out_value))
    return ['COILS', coils_predictions]

def run_iupred(aa_seq):
    """Creates a virtual temp file, adds a fasta header and the amino acids.
    Then passes this file name as argument to the executable and gets the output.
    It checks for errors from the executable's run and raises an error if there is any.
    If not the case, it processes and returns the predictions.
    """
    aa_fasta = f'>CCD_protein\n{aa_seq}'.encode('utf-8')
    iupred_exe = app.config['IUPRED']
    iupred_dir = os.path.dirname(iupred_exe) 
    with tempfile.NamedTemporaryFile() as th, \
         env_var('IUPred_PATH', iupred_dir):
        _log.debug(f'Set IUPRed_PATH environment variable to {iupred_dir}')
        #create fasta file for input
        th.write(aa_fasta)
        th.seek(0)
        _log.info('Running iupred..')
        with Popen([iupred_exe, th.name, 'short'], 
                   stdout=PIPE, stderr=PIPE) as p:
            out, err = p.communicate()
        iupred_prediction = ''
        if err:
            raise IupredError(err)
        for line in out.decode('utf-8').splitlines()[9:]:
            if float(line.split()[2]) > 0.5: #cutoff suggested by iupred
                iupred_prediction += 'd'
            else:
                iupred_prediction += '-'
        assert len(aa_seq) == len(iupred_prediction) 
    return ['IUPRED', iupred_prediction]

def run_predator(aa_seq):
    """Creates a virtual temp file, adds a fasta header and the amino acids.
    Then passes this file name as argument to the executable and gets the output.
    It checks for errors from the executable's run and raises an error if there is any.
    If not the case, it processes and returns the predictions.
    """
    # First FASTA line is a header
    fasta = f'>CCD_protein\n{aa_seq}'.encode('utf-8')
    predator_exe = app.config['PREDATOR']
    datafile = os.path.join(os.path.dirname(predator_exe), 'stride.dat')
    
    with tempfile.NamedTemporaryFile() as th:
        th.write(fasta)
        th.seek(0)
        _log.info('Running predator..')
        cmd = f'{predator_exe} -l -b{datafile} {th.name}' #no space after -b !
        with Popen(cmd.split(), stdout=PIPE, stderr=PIPE) as p:
            out, err = p.communicate()
        if 'ERROR' in err.decode('utf-8'):
            raise PredatorError(err.splitlines()[1].decode('utf-8'))
        out = out.decode('utf-8').splitlines()[3:]
        pred = []
        for line in out:
            pred.append(line.split()[4].replace('c', '-'))
        pred = ''.join(pred)
    return ['PREDATOR', pred]


def run_globplot(aa_seq):
    """Makes a web request to the remote predictions server.
    Processes the reply (HTML) with beautiful soup and extracts the predictions list. 
    In the current precessing two predictions strings are merged into one, and priority is given in the predictions in the second string.
    Returns the predictions in the end.
    """
    url = 'http://globplot.embl.de/cgiDict.py'
    data = 'key=process&sequence_string={}'.format(aa_seq)
    _log.info('Web request to GLOBPLOT..')
    reply = query_server(url, data).decode('utf-8')
    soup1 = BeautifulSoup(reply, 'html.parser')
    res = soup1.find_all('tt')
    annotated = '-'*len(aa_seq)
    #find regions classified as disordered and mark them with 'd'
    #disordered intervals are marked as "start1-stop1, start2-stop2, ..." etc
    disordered = res[1].text.splitlines()[2].split(',')  
    if disordered != ['none']: #comes from globplot when no disorder is found
        #cast intervals to tuples
        disordered = [(int(i.split('-')[0]),
                       int(i.split('-')[1])) for i in disordered]
        annotated = replace_range(annotated, disordered, 'd')
    #find regions classified as globular and mark them with 'D'
    globular = res[2].text.splitlines()[2].split(',')
    if globular != ['none']:
        globular = [(int(i.split('-')[0]),
                     int(i.split('-')[1])) for i in globular]
        annotated = replace_range(annotated, globular, 'G')
    return ['GLOBPLOT', annotated]

#helper function for run_globplot

def replace_range(seq, ranges, filler_symbol): #ranges are counted as 1-based
    annotated = seq
    for rng in ranges:
            start, stop = rng
            head = annotated[:start-1] 
            filler = filler_symbol*(stop-start+1)
            tail = annotated[stop:]
            annotated = f'{head}{filler}{tail}'
    return annotated

def run_smart(aa_seq):
    """Makes a web request to the remote predictions server.
    Processes the reply (HTML) with beautiful soup and extracts the predictions list. 
    Returns the predictions in the end.
    """
    
    base_url = "http://smart.embl-heidelberg.de/smart/"
    smart_script_url = urljoin(base_url, 'show_motifs.pl')
    data = 'SEQUENCE={}'.format(aa_seq)
    _log.info('Web request to SMART..')
    reply = query_server(smart_script_url, data)
    soup = BeautifulSoup(reply, 'html.parser')
    #if the sequence is not precalculated, it might be put in queue
    if 'Sequence submitted into the processing queue' in soup.title.string:
        tags = soup.find_all('a')
        for t in tags:
            if 'job_status.pl?' in t['href']:
                job_url = t['href']
        newreq = urljoin(base_url, job_url)
        #we keep refreshing until done
        while True:
            conn = HTTPConnection("smart.embl-heidelberg.de")
            conn.request("GET", newreq)
            res = conn.getresponse()
            nbody = res.read()
            conn.close()
            soup = BeautifulSoup(nbody, 'html.parser')
            if soup.title.string == 'SMART: Sequence analysis results':
                break
            time.sleep(10)
            soup = BeautifulSoup(reply, 'html.parser')
    # If there are multiple identical sequences in the database (e.g. PCNA_HUMAN
    # SMART will ask the user to choose which one to display 
    if "Multiple matching sequences" in soup.title.string:
        protein_url = soup.find_all('td', {'class': 'prot'})[0].next['href']
        reply = requests.get(urljoin(base_url, protein_url))
        soup = BeautifulSoup(reply.content, 'html.parser') 
    scripts = soup('script', {'type' : 'text/javascript'})
    for s in scripts:
        if 'var domNfo' in s.text:
            script_text = s.text
            break
    else:
        return ['SMART','-'*len(aa_seq),[]] #return dashes only
    domain_info = re.compile(r'var domNfo={((.|\n)*)};')
    match = domain_info.search(script_text)
    if not match: #no domains found in protein
        return ['SMART','-'*len(aa_seq),[]] #return dashes only
    #remove var declaration and final ; to get clean json
    match = match.group().replace('var domNfo=','').replace('};','}') 
    domains = json.loads(match)
    annotated = ['-']*len(aa_seq)
    pdomains = []
    for domain in domains.values():
        start = int(domain['st'])
        stop = int(domain['en'])
        dname = domain['n']
        if dname == 'Low complexity region':
            dname = 'low complexity'
        start = int(domain['st'])
        stop = int(domain['en'])
        dlen = stop - start
        displayed = f'{dname} * * '
        dnlen = len(dname)
        if dlen <= dnlen:
            mult = 1
        else:
            mult = (dlen // dnlen) + 1
        dspan = displayed * mult
        annotated[start-1: stop] = dspan[:dlen+1]
        pdomains.append([start, stop, dname])
    annotated = ''.join(annotated)
    return ['SMART', annotated, pdomains]


def secure_query_server(url, data):
    """Creates and issues a secure POST web request to the specified url.
    A timeout of 3 minutes is set for server reply. 
    Returns the reply in the end.
    """
    urlparts = urlparse(url)
    conn = HTTPSConnection(urlparts.netloc, timeout = 180)
    conn.request("POST", urlparts.path, data)
    resp = conn.getresponse()
    return resp.read()

def query_server(url, data):
    """Creates and issues a POST web request to the specified url.
    A timeout of 3 minutes is set for server reply. 
    Returns the reply in the end.
    """
    urlparts = urlparse(url)
    conn = HTTPConnection(urlparts.netloc, timeout = 180)
    conn.request("POST", urlparts.path, data)
    resp = conn.getresponse()
    return resp.read()

def service(service_name, aa_seq):
    """Takes as arguments the service name and the amino acid sequence for predictions input.
    Unfortunately (#Georgewhatwereyouthinking) creating a display line in the frontend 
    is coupled to launching one prediction. So for the PDB, either I run the same prediction 3 times 
    or I have some hidden things that fake returns stuff from the same persistent object...
    and since I don't really want to refactor the frontend...
    """
    options = {
        'PREDATOR': run_predator,
        'IUPRED': run_iupred,
        'GLOBPLOT': run_globplot,
        'COILS': run_coils,
        'SMART': run_smart,
        }
    return options[service_name](aa_seq)

if __name__ == '__main__':
#     raise NotImplementedError('This module should not be run as main')
    crawler.protein_seq = 'RKRRSPSPKPTKVHIGRLTRNVTKDHIMEIFSTYGKIKMIDMPVERMHPH'
    try:
        res = crawler.run(online=False)
    except BlastNoResult:
        print('no result')
    print(res)
    
