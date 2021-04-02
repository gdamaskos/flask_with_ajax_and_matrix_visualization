# coding: utf8
import asyncio
import aiohttp

from bs4 import BeautifulSoup
from ccd import app, db_utilities

class UrlFormatter(object):
    '''Create url(s) to fetch nucleic acid sequences from databases'''
    
    def __init__(self):
        super().__init__()
        self.ncbi_address = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
        self.CCDS_address = 'https://www.ncbi.nlm.nih.gov/projects/CCDS/'
        self.ncbi_esummary_address = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
        self.email = app.config['EMAIL']
        
    def format(self, database: str, id_list) -> list:
        '''Produce url(s) to retrieve list of identifiers from database
        Input: 
        id_list: one or more valid gene identifiers 
        database: one of ['EMBL', 'RefSeq', 'CCDS', 'Gene']
        Output:
        a list of tuples (database, url).
        If database == CCDS: len(output) == len(id_list)
        Else: len(output) == 1 because EMBL, RefSeq and Gene allow batch retrieval
        '''
        if not isinstance(id_list, list):
            id_list = [id_list]
        allowed = ['EMBL', 'RefSeq', 'CCDS', 'Gene']
        if database not in allowed:
            raise ValueError(f'Supported database types: {allowed}')
        if database in ['EMBL', 'RefSeq']:
            url = f'{self.ncbi_address}efetch.fcgi?db=nuccore&id={",".join(id_list)}&retmode=xml&rettype=gb&email={self.email}'
            return [(database, url)]
        elif database == 'CCDS':
            urls = []
            for id_ in id_list:
                url = f'{self.CCDS_address}CcdsBrowse.cgi?REQUEST=CCDS&DATA={id_}&ORGANISM=0&BUILDS=CURRENTBUILDS'
                urls.append((database, url))
            return urls
        elif database == 'Gene':
            url = f'{self.ncbi_address}efetch.fcgi?db=Gene&id={",".join(id_list)}&retmode=xml&email={self.email}'
            return [(database, url)]
    
    def format_for_summary(self, database: str, id_list) -> list:
        if not isinstance(id_list, list):
            id_list = [id_list]
        allowed = ['EMBL', 'RefSeq']
        if database not in allowed:
            raise ValueError(f'Supported database types: {allowed}')
        url = f'{self.ncbi_esummary_address}?db=nuccore&id={",".join(id_list)}'
        return [('Summary', url)]
        
class Entry_fetcher(object):
    '''Given database and url, fetch entries asynchronously'''
    
    async def _fetch(self, session, database : str, url : str) -> tuple:
        '''Fetch a webpage (coroutine)
        Inout: 
        session = aiohttp session
        database =  a string, that is simply returned
        url = url to fetch
        
        Output:
        (database, response text)
        '''
        async with session.get(url) as response:
            if response.status != 200:
                response.raise_for_status()
            print(f'getting {database}')
            return (database, await response.text())
    
    async def fetch_all(self, queries):
        '''Fetch DNA from NCBI or CCDS asynchronously
        Input:
        list of (database, url) tuples 
        database in ['EMBL', 'Refseq', 'CCDS', 'Gene']
        Return:
        list of (database, response.text) tuples. 
        If response.text of EMBL or Gene or RefSeq contains multiple entries, these are
        returned as separate (database, entry.text) tuples. 
        '''
        #fetching all results
        async with aiohttp.ClientSession() as session:
            tasks = [asyncio.ensure_future(self._fetch(session, query[0], query[1]))
                     for query in queries]
            return await asyncio.gather(*tasks)
            
class EntrySplitter(object):
    
    def split(self, results):
        #refseq and EMBL entries are concatenated. Let's break them up
        entries = []
        for r in results:
            entries += (self._split_entries(r[0],r[1]))
        return entries
    
    def _split_entries(self, database, xml_response):
        soup = BeautifulSoup(xml_response, 'xml')
        if database == 'CCDS':
            return [(database, xml_response)]
        elif database == 'Gene':
            return [(database, str(entry)) for entry in soup.find_all('Entrezgene')]
        elif database in ['EMBL', 'RefSeq']:
            return [(database, str(entry)) for entry in soup.find_all('GBSeq')]
        elif database == 'Summary':
            return [(database, str(entry)) for entry in soup.find_all('DocSum')]
        else:
            raise ValueError

class DnaParser(object):
    
    def parse(self, entries):
        parsed = []
        CCDSParser = db_utilities.CCDSParser()
        GenBankParser = db_utilities.GenBankParser()
        GeneParser = db_utilities.GeneParser()
        for entry in entries:
            database = entry[0]
            text = entry[1]
            if database in ['EMBL', 'RefSeq']:
                parsed.append((database, #tuple
                              GenBankParser.parse(BeautifulSoup(text, 'xml'))))
            elif database == 'CCDS':
                parsed.append((database, 
                              CCDSParser.parse(BeautifulSoup(text, 'html.parser'))))
            elif database == 'Gene':
                parsed.append((database,
                              GeneParser.get_crossreferences(BeautifulSoup(text, 'xml'))))
            else:
                raise ValueError(f'Unknown database: {database}')
        return parsed

class SummaryParser(object):
    
    def get_id_and_size(self, summary_page):
        soup = BeautifulSoup(summary_page, 'xml')
        id_ = soup.find_all('Item', {'Name': 'AccessionVersion'})[0].text
        size = soup.find_all('Item', {'Name': 'Length'})[0].text
        return id_, int(size)
    
if __name__ == '__main__':
    MAX_SIZE = 1000
    loop = asyncio.get_event_loop()
    urlform = UrlFormatter()
    fetcher = Entry_fetcher()
    summary_parser = SummaryParser()
    ids = [('EMBL', ['CR456855.1', 'DQ917642.1']),
           ('RefSeq', ['NM_001270952.1']),
           ('CCDS', ['CCDS73586.1', 'CCDS86041.1']),
           ('Gene', ['7347','50933'])]
    formatter = UrlFormatter()
    splitter = EntrySplitter()
    summaries = []
    #check for size
    for id_list in ids:
        if id_list[0] in ['EMBL', 'RefSeq']:
            summaries += formatter.format_for_summary(id_list[0], id_list[1])
#     for id_list in ids:
#         queries += formatter.format_for_summary(id_list[0], id_list[1])
    entries = loop.run_until_complete(fetcher.fetch_all(summaries))
    res = splitter.split(entries)
    retrieve = []
    for page in res:
        id_, size = summary_parser.get_id_and_size(page[1])
        print(id_, size)
        if size < MAX_SIZE:
            retrieve.append(id_)
    for x in ids:
        for id_ in list(x[1]):
            if id_ not in retrieve:
                x[1].remove(id_)
    print(ids) 
#     parser = DnaParser()
#     print(parser.parse(entries))

    loop.close()
    