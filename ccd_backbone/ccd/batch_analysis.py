import asyncio
import logging
from subprocess import Popen, PIPE
import tempfile

from bs4 import BeautifulSoup


from ccd import app
from ccd.services.predict import CoilsError
from ccd.homolog_finder import Evolutionator
from ccd.sequence_classes import UniprotProtein
import ccd.db_utilities_async as dba
import ccd.services.predict as pred 
from builtins import str

class Result (object):
    
    def __init__(self, service: str, id_: str, result: str, ):
        self.service = service
        self.id_= id_
        self.result = result
        
    def __repr__(self):
        return f'Prediction result of {self.service} for protein {self.id_}'
    

def get_uniprot_data(id_list):
    '''Gather data necessary for starting Evolutionator in batch, rather than from 
    the GUI'''
    formatter = dba.UrlFormatter()
    queries = formatter.format('Uniprot', id_list, format_='batch_xml')
    fetcher = dba.Entry_fetcher()
    loop = asyncio.get_event_loop()
    uniprot_xml = loop.run_until_complete(fetcher.fetch_all(queries))
    splitter = dba.EntrySplitter()
    entries = splitter.split(uniprot_xml)
    proteins = []
    for entry in list(entries):
        soup = BeautifulSoup(entry[1], 'xml')
        parser = dba.UniprotParser(soup, None)
        protein_id = parser.match_entry_to_id(id_list)
        p = UniprotProtein(protein_id)
        p.uniprot_xml_soup = soup
        p.source_organism = parser.get_scientific_name()
        p.source_organism_common = parser.get_common_name()
        p.full_name, p.short_name = parser.get_protein_names()
        p.kingdom = parser.get_kingdom()
        p.sequence = parser.get_protein_sequence()
        proteins.append(p)
    return proteins

async def get_homologs(protein):
    data = (protein.id_, protein.sequence, protein.source_organism)
    e = Evolutionator()
    e.configure_servers()
    e.fill_data(*data)
    return (protein.id_, await e.search())

@asyncio.coroutine
def run_coils_coro(id_, aa_seq):
    return Result('COILS', id_[0], pred.run_coils(aa_seq)[1])

@asyncio.coroutine
def run_iupred_coro(id_, aa_seq):
    return Result('IUPRED', id_[0], pred.run_iupred(aa_seq)[1])

@asyncio.coroutine
def run_predator_coro(id_, aa_seq):
    return Result('PREDATOR', id_[0], pred.run_predator(aa_seq)[1])

def main():
    from pprint import pprint
    id_list = ['Q9UJ41','Q15287','Q96RL1']
    url = 'https://www.uniprot.org/uniprot/?query=method:X-ray&format=xml&force=true&sort=score&compress=yes'
    data = get_uniprot_data(id_list)
    loop = asyncio.get_event_loop()
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    coros = []
#     coros = [asyncio.ensure_future(get_homologs(protein)) for protein in data]
    coros += [asyncio.ensure_future(run_coils_coro(protein.id_, protein.sequence))
                        for protein in data]
    coros += [asyncio.ensure_future(run_iupred_coro(protein.id_, protein.sequence))
                        for protein in data]
    coros += [asyncio.ensure_future(run_predator_coro(protein.id_, protein.sequence))
                        for protein in data]
    tasks = asyncio.gather(*coros)
    res = loop.run_until_complete(tasks)
    pprint(res)
    
    '''
    To limit aiohttp connection pool:
    connector = aiohttp.TCPConnector(limit=50) #global
    connector = aiohttp.TCPConnector(limit_per_host=50) #per host
    client = aiohttp.ClientSession(connector=connector)
    '''
    
        
if __name__ == '__main__':
    main()
