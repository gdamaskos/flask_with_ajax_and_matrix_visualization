# -*- coding: utf-8 -*-

import time
import os
import sys
import requests
import json

from tempfile import NamedTemporaryFile
from subprocess import Popen, PIPE

from bs4 import BeautifulSoup

try:
    import ccd
except ImportError:
    sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from ccd import app

BASE_URL = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'
PDB_URL = 'https://www.ebi.ac.uk/pdbe/entry/pdb'
CONTACT = '{}'.format(app.config['EMAIL'])
EVALUE_THRESHOLD = 1e-4
HIGH_IDENTITY_THRESHOLD = 95
MEDIUM_IDENTITY_THRESHOLD = 50
LOW_IDENTITY_THRESHOLD = 30
#inspect.getfile(ccd) gets the correct app path no matter where the script is run from
BLAST_DB = os.path.join(app.root_path, 'static/blast_db/pdb_blastdb') 
blast_exe = app.config['BLAST']

class BlastFailError(Exception):
    pass

class BlastTimeoutError(Exception):
    pass

class BlastNoResult(Exception):
    pass

class Hsp(object):
    '''
    Stores the data of a high scoring pair of a blast result. 
    I wish I had data classes available in 2.7 ;-)
    '''
    def __init__(self, attrs):
        super(Hsp, self).__init__()
        for key in list(attrs.keys()):
            setattr(self, key, attrs[key])  
            try:
                self.percent_identity = float(self.identity) / float(self.length) *100
            except AttributeError: #makes testing easier I am lazy
                pass
        #remove the unique numeric identifier e.g. >2C7M1_A_12345 -> 2C7M_A
        self.id_ = '_'.join(self.id_.split('_')[:-1])
    
    def __repr__(self):
        return '<Hsp: {}>'.format(self.__dict__)
    
    def __len__(self):
        return len(self.span())
    
    def span(self):
        return set(range(self.start-1, self.stop)) #start and stop count from 1
    
    def spanned_sequence(self, seq):
        return seq[self.start-1:self.stop]
        
class Hsp_strainer(object):
    '''
    Given a bunch of hsps from blast, it finds the smallest possible subset
    that covers the biggest sequence area.
    It also marks them on the sequence.
    '''
    def __init__(self):
        super(Hsp_strainer, self).__init__()
        
    def count_structures_per_residue(self, sequence, hsps):
        #dictionaries cannot have numbers as keys, so initialize a pseudo counting
        #list with structure to count how many hits per residue
        #shape: [ [count1, [structure list], [count2, [structure_list], ... ] 
        #of length == len(sequence);  
        count_list = []
        for _ in range(len(sequence)):
            count_list.append([0,[]])
        assert len(count_list) == len(sequence)
        for hsp in hsps: #add the span of each hsp to the relative portion of count_list
            for index in range(hsp.start-1, hsp.stop): #hsps count starting from 1, internally we use 0
                try:
                    count_list[index][0] += 1
                    count_list[index][1].append(hsp.id_)
                except IndexError:
                    print(index, len(sequence), len(count_list))
        return count_list 
    
    def decorate_by_counts(self, sequence, counts): 
        decorated = [] 
        for pos in counts:
            count = pos[0] 
            if count and count <=9: 
                decorated.append(str(pos[0])) 
            elif count>9: 
                decorated.append('0') 
            else: 
                decorated.append('-') 
        assert len(decorated) == len(sequence)
        return ''.join(decorated)
    
    def count_boundaries(self, sequence, hsps):
        counts = {}
        for pos in range(len(sequence)):
            counts[str(pos)] = [[], []] #hsps starting, # hsps stopping
        for h in hsps:
            counts[str(h.start-1)][0].append(h.id_) #BLAST counts from 1, lists from 0
            counts[str(h.stop-1)][1].append(h.id_)
        return counts
     
    def decorate_by_boundaries(self, sequence, boundaries):
        seq = ''
        for pos in range(len(sequence)):
            has_start = bool(boundaries[str(pos)][0])
            has_stop = bool(boundaries[str(pos)][1])
            if has_start and has_stop: 
                seq += '\u25CA'
            elif has_start: 
                seq += '>'
            elif has_stop:
                seq += '<'
            else: 
                seq += '-'
        return seq
    
    def get_boundaries_by_id(self, hsps):
        boundaries_by_id = {}
        for h in hsps:
            try:
                boundaries_by_id[str(h.start-1)].append(h.id_.upper())#-1 because frontend counts from 0
            except KeyError:
                boundaries_by_id[str(h.start-1)] = [h.id_.upper()]
            try:
                boundaries_by_id[str(h.stop-1)].append(h.id_.upper())
            except KeyError:
                boundaries_by_id[str(h.stop-1)] = [h.id_.upper()]
        return boundaries_by_id
    
class PdbCrawler(object):
    
    def __init__(self, protein_seq='', id_=''):
        super(PdbCrawler, self).__init__()
        self.base_url = BASE_URL
        if not (protein_seq or id_):
            msg = 'One of either protein_seq or id_ must not be empty'
            raise ValueError(msg)
        if '_' in id_:
            msg = 'Only primary uniprot accession number are supported (e.g. Q15287, NOT RNPS1_HUMAN)'
            raise ValueError(msg)
        self.protein_seq = protein_seq
        self.id_ = id_
        self.pdb90, self.pdb50, self.pdb30 = None, None, None
    
    def run_online_search(self):
        params = {'QUERY': '',
                  'PROGRAM': 'blastp',
                  'DATABASE': 'pdbaa',
                  'CMD': 'Put'}
        if self.id_:
            params['QUERY'] = self.id_.upper()
        else:
            params['QUERY'] = self.protein_seq.upper()
        #funny: requests already puts User Agent; sending it twice causes a 400
        response = requests.get(self.base_url, params=params)
        search_id, exp_time_to_completion = self.get_search_id(response)
        status = 'WAITING'
        iterations = 1
        while status == 'WAITING':
            print('Iteration {}: status is {}. ETA is {}'.format(iterations, 
                                                                 status,
                                                                 exp_time_to_completion))
            time.sleep(int(exp_time_to_completion)/2)
            status, response = self.check_status(search_id)
            if status == 'READY':
                break
            if iterations > 10:
                raise BlastTimeoutError
            iterations += 1
        if not r'ThereAreHits=yes' in response:
            return None
        else:
            blast_result = self.get_results(search_id)
            hsps = self.parse_result(blast_result)
        return hsps
    
    def run_offline_search(self):
        #check if blast installed locally
        try:
            p=Popen(blast_exe, stderr=PIPE, stdout=PIPE)
            p.communicate()
        except OSError:
            raise OSError('Cannot find the blast executable (blastp). Please install blast')
        #Offline cannot search by id
        if not self.protein_seq:
            raise ValueError('Offline searches by sequence only')
        #blastp works with file inputs. delete=False is for windows to work properly
        with NamedTemporaryFile(delete=False) as query, NamedTemporaryFile(delete=False) as out:
            #write query
            query.write(self.protein_seq.encode('utf-8'))
            query.seek(0)
            #run blast
            cmd = f'{blast_exe} -query {query.name} -db {BLAST_DB} -out {out.name} -outfmt 5'
            p = Popen(cmd.split(), stdout=PIPE, stderr=PIPE)
            o, e = p.communicate() #blastp uses stdout to output errors
            if e or o:
                raise BlastFailError(e)
            result = out.read().decode('utf-8')
        if r'No hits found' in result:
            return None
        else:
            hsps = self.parse_result(result)
        return hsps
    
    def run(self, online=False):
        if online:
            hsps = self.run_online_search() #query NCBI (note: max 5 queries hour)
        else:
            hsps = self.run_offline_search() #query local database
        if not hsps:
            raise BlastNoResult
        pdb_90, pdb_50, pdb_30 = self.bin_hsps_by_identity(hsps)
        strainer = Hsp_strainer()
        to_frontend = []
        for s in (pdb_90, pdb_50, pdb_30):
            counts = strainer.count_structures_per_residue(self.protein_seq, s) 
            decorated = strainer.decorate_by_counts(self.protein_seq, counts) 
            boundaries = strainer.count_boundaries(self.protein_seq, s)
            bounded = strainer.decorate_by_boundaries(self.protein_seq, boundaries)
            to_frontend.append([bounded, decorated])
        #each result is: (string_to_display, string_with_color_info)
        res_pdb_90, res_pdb_50, res_pdb_30 = to_frontend 
        #for pdb_90, we also add extra data detailing which structures start/stop 
        #at which boundary. for the other predictions. we leave this as emtpy dictionary
        boundaries_by_id = strainer.get_boundaries_by_id(pdb_90)
        res_pdb_90.append(json.dumps(boundaries_by_id))
        res_pdb_50.append(json.dumps({"0":['empty']}))
        res_pdb_30.append(json.dumps({"0":['empty']}))
        return res_pdb_90, res_pdb_50, res_pdb_30
     
    def __run_debug(self): 
        import os, pickle  
        with open(os.path.abspath('../tests/unit/static/BLASTRes.pickle')) as f:
            return pickle.load(f) 
     
    def get_search_id(self, response):
        r = response.text
        start = r.find('<!--QBlastInfoBegin')
        stop = r.find('QBlastInfoEnd\n-->') + 17
        # RID block = <!--QBlastInfoBegin\n     RID = [request_id]\n     RTOE = 18\n QBlastInfoEnd\n-->
        RID_block = r[start:stop]
        x = [x.strip(' \t\n\r') for x in RID_block.split('\n')]
        request_id = [i for i in x if 'RID' in i][0].split(' =')[1].strip()
        exp_time_to_completion =  [i for i in x if 'RTOE' in i][0].split(' =')[1].strip() # in seconds
        return request_id, exp_time_to_completion
    
    def parse_result(self, blast_result):
        soup = BeautifulSoup(blast_result, 'lxml')
#         dest = '/tmp/tmpresult.html'
#         with open(dest, 'w+') as f:
#             f.write(blast_result)
        #warning, lxml turns the original "Hit" tag to "hit" (lowercase),
        #because reasons... it took me forever to find this out...
        hits = soup.findAll('hit')
        hsps = [] 
        for hit in hits:
            attrs = {'id_': hit.hit_accession.text,
                     'length': hit.hit_len.text}
            for h in hit.hit_hsps:
                h = str(h).replace('-','_') #'-' in tags gives invalid attribute names in beautifulsoup
                
                try:
                    soup = BeautifulSoup(str(h), 'lxml')
                    s = soup.findAll('hsp')[0]
                except IndexError: #h == '\n'
                    continue
                try:
                    data = {'score': float(s.hsp_score.text),
                            'bit_score': float(s.hsp_bit_score.text),
                           'evalue': float(s.hsp_evalue.text.replace('_','-')),
                           'start': int(s.hsp_query_from.text),
                           'stop': int(s.hsp_query_to.text),
                           'id_': attrs['id_'],
                           'length': int(s.hsp_align_len.text),
                           'identity': int(s.hsp_identity.text),
                           'similarity': int(s.hsp_positive.text)}
                except:
                    pass
                hsps.append(Hsp(data))
        return hsps 

    def check_status(self, request_id):
        params = {'RID': request_id,
                  'FORMAT_OBJECT': 'SearchInfo',
                  'CMD': 'Get'}
        result = requests.get(self.base_url, params=params)
        r = result.text
        start = r.find('QBlastInfoBegin')
        stop = r.find('QBlastInfoEnd') + 13
        # Status block = QBlastInfoBegin\n     Status=[READY|FAILED|UNKNOWN]\n QBlastInfoEnd\n
        status = [x.strip(' \t\n\r') for x in r[start:stop].split('\n') 
                  if 'Status' in x][0].split('=')[1]
        if status in ['UNNOWN', 'FAILED']:
            msg = 'Blast failed with error status {}'.format(status)
            raise BlastFailError(msg)
        return status, r
    
    def get_results(self, request_id):
        params = {'RID': request_id,
                  'CMD': 'Get',
                  'FORMAT_TYPE': 'xml'}
        r = requests.get(self.base_url, params=params)
        return r.text
    
    def bin_hsps_by_identity(self, hsps):
        identical = [hsp for hsp in hsps if hsp.percent_identity >= HIGH_IDENTITY_THRESHOLD]
        similar_50 = [hsp for hsp in hsps 
                   if hsp.percent_identity >= MEDIUM_IDENTITY_THRESHOLD
                   if hsp.percent_identity < HIGH_IDENTITY_THRESHOLD
                   if hsp.evalue < 1e-4]
        similar_30 = [hsp for hsp in hsps 
                   if hsp.percent_identity < MEDIUM_IDENTITY_THRESHOLD
                   if hsp.percent_identity >= LOW_IDENTITY_THRESHOLD
                   if hsp.evalue < EVALUE_THRESHOLD]
        return identical, similar_50, similar_30
        
    def requests_query(self, url, params):
        r = requests.get(url, params=params)
        r.raise_for_status()
        return r.text

if __name__ == '__main__':
    from pprint import pprint
    rnps1_id ='Q15287' 
    rnps1_seq = '''
    MDLSGVKKKSLLGVKENNKKSSTRAPSPTKRKDRSDEKSKDRSKDKGATKESSEKDRGRD
    KTRKRRSASSGSSSTRSRSSSTSSSGSSTSTGSSSGSSSSSASSRSGSSSTSRSSSSSSS
    SGSPSPSRRRHDNRRRSRSKSKPPKRDEKERKRRSPSPKPTKVHIGRLTRNVTKDHIMEI
    FSTYGKIKMIDMPVERMHPHLSKGYAYVEFENPDEAEKALKHMDGGQIDGQEITATAVLA
    PWPRPPPRRFSPPRRMLPPPPMWRRSPPRMRRRSRSPRRRSPVRRRSRSPGRRRHRSRSS
    SNSSR
    '''.replace('\n', '').replace('\t','').replace(' ','')
    usp7_id = 'Q93009'
    usp7_seq = '''
    MNHQQQQQQQKAGEQQLSEPEDMEMEAGDTDDPPRITQNPVINGNVALSDGHNTAEEDME
    DDTSWRSEATFQFTVERFSRLSESVLSPPCFVRNLPWKIMVMPRFYPDRPHQKSVGFFLQ
    CNAESDSTSWSCHAQAVLKIINYRDDEKSFSRRISHLFFHKENDWGFSNFMAWSEVTDPE
    KGFIDDDKVTFEVFVQADAPHGVAWDSKKHTGYVGLKNQGATCYMNSLLQTLFFTNQLRK
    AVYMMPTEGDDSSKSVPLALQRVFYELQHSDKPVGTKKLTKSFGWETLDSFMQHDVQELC
    RVLLDNVENKMKGTCVEGTIPKLFRGKMVSYIQCKEVDYRSDRREDYYDIQLSIKGKKNI
    FESFVDYVAVEQLDGDNKYDAGEHGLQEAEKGVKFLTLPPVLHLQLMRFMYDPQTDQNIK
    INDRFEFPEQLPLDEFLQKTDPKDPANYILHAVLVHSGDNHGGHYVVYLNPKGDGKWCKF
    DDDVVSRCTKEEAIEHNYGGHDDDLSVRHCTNAYMLVYIRESKLSEVLQAVTDHDIPQQL
    VERLQEEKRIEAQKRKERQEAHLYMQVQIVAEDQFCGHQGNDMYDEEKVKYTVFKVLKNS
    SLAEFVQSLSQTMGFPQDQIRLWPMQARSNGTKRPAMLDNEADGNKTMIELSDNENPWTI
    FLETVDPELAASGATLPKFDKDHDVMLFLKMYDPKTRSLNYCGHIYTPISCKIRDLLPVM
    CDRAGFIQDTSLILYEEVKPNLTERIQDYDVSLDKALDELMDGDIIVFQKDDPENDNSEL
    PTAKEYFRDLYHRVDVIFCDKTIPNDPGFVVTLSNRMNYFQVAKTVAQRLNTDPMLLQFF
    KSQGYRDGPGNPLRHNYEGTLRDLLQFFKPRQPKKLYYQQLKMKITDFENRRSFKCIWLN
    SQFREEEITLYPDKHGCVRDLLEECKKAVELGEKASGKLRLLEIVSYKIIGVHQEDELLE
    CLSPATSRTFRIEEIPLDQVDIDKENEMLVTVAHFHKEVFGTFGIPFLLRIHQGEHFREV
    MKRIQSLLDIQEKEFEKFKFAIVMMGRHQYINEDEYEVNLKDFEPQPGNMSHPRPWLGLD
    HFNKAPKRSRYTYLEKAIKIHN'''.replace('\n', '').replace('\t','').replace(' ','')
    rbx5_id = 'Q9UJ41'
    rbx5_seq = '''
    MSLKSERRGIHVDQSDLLCKKGCGYYGNPAWQGFCSKCWREEYHKARQKQIQEDWELAER
    LQREEEEAFASSQSSQGAQSLTFSKFEEKKTNEKTRKVTTVKKFFSASSRVGSKKEIQEA
    KAPSPSINRQTSIETDRVSKEFIEFLKTFHKTGQEIYKQTKLFLEGMHYKRDLSIEEQSE
    CAQDFYHNVAERMQTRGKVPPERVEKIMDQIEKYIMTRLYKYVFCPETTDDEKKDLAIQK
    RIRALRWVTPQMLCVPVNEDIPEVSDMVVKAITDIIEMDSKRVPRDKLACITKCSKHIFN
    AIKITKNEPASADDFLPTLIYIVLKGNPPRLQSNIQYITRFCNPSRLMTGEDGYYFTNLC
    CAVAFIEKLDAQSLNLSQEDFDRYMSGQTSPRKQEAESWSPDACLGVKQMYKNLDLLSQL
    NERQERIMNEAKKLEKDLIDWTDGIAREVQDIVEKYPLEIKPPNQPLAAIDSENVENDKL
    PPPLQPQVYAG'''.replace('\n', '').replace('\t','').replace(' ','')
    acinus_id = 'Q9UKV3'
    acinus_sequence = '''
    MWRRKHPRTSGGTRGVLSGNRGVEYGSGRGHLGTFEGRWRKLPKMPEAVGTDPSTSRKMA
    ELEEVTLDGKPLQALRVTDLKAALEQRGLAKSGQKSALVKRLKGALMLENLQKHSTPHAA
    FQPNSQIGEEMSQNSFIKQYLEKQQELLRQRLEREAREAAELEEASAESEDEMIHPEGVA
    SLLPPDFQSSLERPELELSRHSPRKSSSISEEKGDSDDEKPRKGERRSSRVRQARAAKLS
    EGSQPAEEEEDQETPSRNLRVRADRNLKTEEEEEEEEEEEEDDEEEEGDDEGQKSREAPI
    LKEFKEEGEEIPRVKPEEMMDERPKTRSQEQEVLERGGRFTRSQEEARKSHLARQQQEKE
    MKTTSPLEEEEREIKSSQGLKEKSKSPSPPRLTEDRKKASLVALPEQTASEEETPPPLLT
    KEASSPPPHPQLHSEEEIEPMEGPAPAVLIQLSPPNTDADTRELLVSQHTVQLVGGLSPL
    SSPSDTKAESPAEKVPEESVLPLVQKSTLADYSAQKDLEPESDRSAQPLPLKIEELALAK
    GITEECLKQPSLEQKEGRRASHTLLPSHRLKQSADSSSSRSSSSSSSSSRSRSRSPDSSG
    SRSHSPLRSKQRDVAQARTHANPRGRPKMGSRSTSESRSRSRSRSRSASSNSRKSLSPGV
    SRDSSTSYTETKDPSSGQEVATPPVPQLQVCEPKERTSTSSSSVQARRLSQPESAEKHVT
    QRLQPERGSPKKCEAEEAEPPAATQPQTSETQTSHLPESERIHHTVEEKEEVTMDTSENR
    PENDVPEPPMPIADQVSNDDRPEGSVEDEEKKESSLPKSFKRKISVVSATKGVPAGNSDT
    EGGQPGRKRRWGASTATTQKKPSISITTESLKSLIPDIKPLAGQEAVVDLHADDSRISED
    ETERNGDDGTHDKGLKICRTVTQVVPAEGQENGQREEEEEEKEPEAEPPVPPQVSVEVAL
    PPPAEHEVKKVTLGDTLTRRSISQQKSGVSITIDDPVRTAQVPSPPRGKISNIVHISNLV
    RPFTLGQLKELLGRTGTLVEEAFWIDKIKSHCFVTYSTVEEAVATRTALHGVKWPQSNPK
    FLCADYAEQDELDYHRGLLVDRPSETKTEEQGIPRPLHPPPPPPVQPPQHPRAEQREQER
    AVREQWAEREREMERRERTRSEREWDRDKVREGPRSRSRSRDRRRKERAKSKEKKSEKKE
    KAQEEPPAKLLDDLFRKTKAAPCIYWLPLTDSQIVQKEAERAERAKEREKRRKEQEEEEQ
    KEREKEAERERNRQLEREKRREHSRERDRERERERERDRGDRDRDRERDRERGRERDRRD
    TKRHSRSRSRSTPVRDRGGRR
    '''.replace('\n', '').replace('\t','').replace(' ','')
    rnps1_short_seq = 'PSRRRHDNRRRSRSKSKPPKRDEKERKRRSPSPKPTKVHIGRLTRNVTKDHIMEIFSTYGKIKMIDMPVERMHPHLSKGYAYVEFENPDEAEKALKHMDGGQIDGQEITATAVLAPWPRPPPRRFSPP'
    crawler = PdbCrawler(protein_seq=acinus_sequence)
    pdb_90,pdb_50,pdb_30 = crawler.run(online=False)
    pprint(str(acinus_sequence))
    pprint(pdb_90)
#     for line in (pdb_90,pdb_50,pdb_30):
#         print((str(line)))
