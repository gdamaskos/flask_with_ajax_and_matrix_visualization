import os

import unittest

from pprint import pprint

from ccd.services import predict
from ccd import app

class Test_smart(unittest.TestCase):
    
    def test_multiple_proteins_matched(self):
        res = predict.run_smart(sap18_sequence)
        exp = ['SMART',
               '---------low complex-------------------------------------------------------------------------------------------------------------------------------------', 
               [[10, 20, 'low complexity']]]
        self.assertEqual(res[:2], exp[:2])
        self.assertEqual(sorted(res[2]), sorted(exp[2]))
        print(res[2])
        
    def test_no_domain_returned_and_multiple_matches(self):
        res = predict.run_smart(PCNA_sequence)
        exp = ['SMART', '-'*len(PCNA_sequence), []]
        self.assertEqual(res[:2], exp[:2])
        self.assertEqual(sorted(res[2]), sorted(exp[2]))
        
    def test_search_queried(self):
        res = predict.run_smart(obscure_sequence)
        exp = ['SMART',
               '-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------low complexity * *-----------------',
               [[236, 253, 'low complexity']]]
        self.assertEqual(res[:2], exp[:2])
        self.assertEqual(sorted(res[2]), sorted(exp[2]))
        
    def test_normal_operation(self):
        res = predict.run_smart(rnps1_sequence)
        exp = ['SMART',
               '-----------------------------low complexity * --------low complexity * * low complexity * * low complexity * * low complexity * * low complexity * * low comp----RRM * * RRM * * RRM * * RRM * * RRM * * RRM * * RRM * * RRM * * RRM * * RRM------low complexity * * low complexity * * low complexity * * low co',
               [[30, 46, 'low complexity'], [55, 157, 'low complexity'], [162, 236, 'RRM'], [243, 305, 'low complexity']]]
        self.assertEqual(res[:2], exp[:2])
        self.assertEqual(sorted(res[2]), sorted(exp[2]))
class Test_predator(unittest.TestCase):
    
    def test_executable_is_present(self):
        exe = app.config['PREDATOR']
        path = os.path.dirname(exe)
        self.assertTrue(os.path.isdir(path), 'Predator folder not found')
        self.assertTrue(os.path.isfile(exe), 'Predator executable not found')
        
    def test_normal_run(self):
        res = predict.run_predator(rnps1_sequence)
        exp = ['PREDATOR', predator_exp]
        self.assertEqual(len(res[1]), len(exp[1]), 'Wrong length of prediction')
        self.assertEqual(res, exp, 'Predator prediction not correct')
        
class Test_coils(unittest.TestCase):
    
    def test_executable_is_present(self):
        exe = app.config['NCOILS']
        path = os.path.dirname(exe)
        self.assertTrue(os.path.isdir(path), 'Ncoils folder not found')
        self.assertTrue(os.path.isfile(exe), 'Ncoils executable not found')
    
    def test_normal_run_no_coils(self):
        res = predict.run_coils(rnps1_sequence)
        exp = ['COILS', '-'*len(rnps1_sequence)]
        self.assertEqual(len(res[1]), len(exp[1]), 'Wrong length of prediction')
        self.assertEqual(res, exp)
        
    def test_normal_run_coils_predicted(self):
        res = predict.run_coils(acinus_sequence)
        exp = ['COILS', coils_exp]
        self.assertEqual(len(res[1]), len(exp[1]), 'Wrong length of prediction')
        self.assertEqual(res, exp)
        
class Test_iupred(unittest.TestCase):
    
    def test_executable_is_present(self):
        exe = app.config['IUPRED']
        path = os.path.dirname(exe)
        self.assertTrue(os.path.isdir(path), 'Iupred folder not found')
        self.assertTrue(os.path.isfile(exe), 'Iupred executable not found')
    
    def test_normal_run(self):
        res = predict.run_iupred(sap18_sequence)
        exp = ['IUPRED', iupred_exp]
        self.assertEqual(len(res[1]), len(exp[1]), 'Wrong length of prediction')
        self.assertEqual(res, exp)
        
class Test_globplot(unittest.TestCase):
    
    def test_normal_run(self):
        res = predict.run_globplot(rnps1_sequence)
        exp = ['GLOBPLOT', globplot_exp]
        self.assertEqual(len(res[1]), len(exp[1]), 'Wrong length of prediction')
        self.assertEqual(res, exp)
        
    def test_no_disorder(self):
        res = predict.run_globplot(ubiquitin_sequence)
        exp = ['GLOBPLOT', 'G'*len(ubiquitin_sequence)]
        self.assertEqual(len(res[1]), len(exp[1]), 'Wrong length of prediction')
        self.assertEqual(res, exp)
        
    def test_no_globular(self):
        seq = 'NKKSSTRAPSPTKRKDRSDEKSKDRSKDKGATKESSEKDRGRD'
        res = predict.run_globplot(seq)
        exp = ['GLOBPLOT', '--ddddddddddddddddddddddddddddddddddddddddd']
        self.assertEqual(len(res[1]), len(exp[1]), 'Wrong length of prediction')
        self.assertEqual(res, exp)
        
class Test_replace_range(unittest.TestCase):
    
    def test_no_range(self):
        res = predict.replace_range('aaaa', [], 'x')
        self.assertEqual(res, 'aaaa')
    
    def test_one_range(self):
        res = predict.replace_range('aaaa', [(2,3)], 'x')
        self.assertEqual(res, 'axxa')
        
    def test_multiple_ranges(self):
        res = predict.replace_range('aaaabbbb', [(2,3), (5,7)], 'x')
        exp = 'axxaxxxb'
        self.assertEqual(res, exp)
    
    def test_overlapping_ranges(self):
        res = predict.replace_range('aaaabbbb', [(2,3), (2,7)], 'x')
        exp = 'axxxxxxb'
        self.assertEqual(res, exp)
        
PCNA_sequence ='''
MFEARLVQGSILKKVLEALKDLINEACWDISSSGVNLQSMDSSHVSLVQLTLRSEGFDTY
RCDRNLAMGVNLTSMSKILKCAGNEDIITLRAEDNADTLALVFEAPNQEKVSDYEMKLMD
LDVEQLGIPEQEYSCVVKMPSGEFARICRDLSHIGDAVVISCAKDGVKFSASGELGNGNI
KLSQTSNVDKEEEAVTIEMNEPVQLTFALRYLNFFTKATPLSSTVTLSMSADVPLVVEYK
IADMGHLKYYLAPKIEDEEGS
'''.replace('\n', '') 
sap18_sequence='''
MAVESRVTQEEIKKEPEKPIDREKTCPLLLRVFTTNNGRHHRMDEFSRGNVPSSELQIYT
WMDATLKELTSLVKEVYPEARKKGTHFNFAIVFTDVKRPGYRVKEIGSTMSGRKGTDDSM
TLQSQKFQIGDYLDIAITPPNRAPPPSGRMRPY
'''.replace('\n', '')
#tr|A0A0U8TRG9|A0A0U8TRG9_TETTH - unlikely it will ever be mapped by id
obscure_sequence ='''
MSNRVQGGFDNNSGNNQSAQKQQAEKIPQITVPLNCFMINQIVKAAKENPQAHSGNHYEW
YGAFENAIITAKFEFLQSINDSPKIMGKLSDSTGCIEVVIQKSKMSDELPEFVQAYEIEL
QNNGNRHKYVRAMLKMRKNAQIQLLYFSIVNDANEISRHGLDLCLRYLQRKHGIEDFMHM
TNDKAHNNHNASAQKVHYQIDRNQQPKEQVLELMRQILKHNPNDQIPKSKIIEFFQSQLN
QVQINQILQQLVSANEIFSVGSDNYLLNV 
'''.replace('\n', '')
rnps1_sequence = '''
MDLSGVKKKSLLGVKENNKKSSTRAPSPTKRKDRSDEKSKDRSKDKGATKESSEKDRGRD
KTRKRRSASSGSSSTRSRSSSTSSSGSSTSTGSSSGSSSSSASSRSGSSSTSRSSSSSSS
SGSPSPSRRRHDNRRRSRSKSKPPKRDEKERKRRSPSPKPTKVHIGRLTRNVTKDHIMEI
FSTYGKIKMIDMPVERMHPHLSKGYAYVEFENPDEAEKALKHMDGGQIDGQEITATAVLA
PWPRPPPRRFSPPRRMLPPPPMWRRSPPRMRRRSRSPRRRSPVRRRSRSPGRRRHRSRSS
SNSSR
'''.replace('\n', '')
acinus_sequence = '''ITIDDPVRTA 
QVPSPPRGKI SNIVHISNLV RPFTLGQLKE LLGRTGTLVE EAFWIDKIKS 
HCFVTYSTVE EAVATRTALH GVKWPQSNPK FLCADYAEQD ELDYHRGLLV 
DRPSETKTEE QGIPRPLHPP PPPPVQPPQH PRAEQREQER AVREQWAERE 
REMERRERTR SEREWDRDKV REGPRSRSRS RDRRRKERAK SKEKKSEKKE 
KAQEEPPAKL LDDLFRKTKA APCIYWLPLT DSQIVQKEAE RAERAKEREK 
RRKEQEEEEQ KEREKEAERE RNRQLEREKR REHSRERDRE
'''.replace(' ', '').replace('\n', '')
ubiquitin_sequence = '''MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYN
IQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLI
FAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKA
KIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKT
ITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLR
LRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTL
SDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQ
QRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTLTGKTITLEVEPSDTIE
NVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGMQIFVKTL
TGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLH
LVLRLRGGMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLED
GRTLSDYNIQKESTLHLVLRLRGGV'''.replace('\n','')
acinus_full_sequence = '''MWRRKHPRTSGGTRGVLSGNRGVEYGSGRGHLGTFEGRWRKLPKMPEAVGTDPS
TSRKMAELEEVTLDGKPLQALRVTDLKAALEQRGLAKSGQKSALVKRLKGALMLENLQKHSTPHAAFQPNSQIGEEMSQN
SFIKQYLEKQQELLRQRLEREAREAAELEEASAESEDEMIHPEGVASLLPPDFQSSLERPELELSRHSPRKSSSISEEKG
DSDDEKPRKGERRSSRVRQARAAKLSEGSQPAEEEEDQETPSRNLRVRADRNLKTEEEEEEEEEEEEDDEEEEGDDEGQK
SREAPILKEFKEEGEEIPRVKPEEMMDERPKTRSQEQEVLERGGRFTRSQEEARKSHLARQQQEKEMKTTSPLEEEEREI
KSSQGLKEKSKSPSPPRLTEDRKKASLVALPEQTASEEETPPPLLTKEASSPPPHPQLHSEEEIEPMEGPAPPVLIQLSP
PNTDADTRELLVSQHTVQLVGGLSPLSSPSDTKAESPAEKVPEESVLPLVQKSTLADYSAQKDLEPESDRSAQPLPLKIE
ELALAKGITEECLKQPSLEQKEGRRASHTLLPSHRLKQSADSSSSRSSSSSSSSSRSRSRSPDSSGSRSHSPLRSKQRDV
AQARTHANPRGRPKMGSRSTSESRSRSRSRSRSASSNSRKSLSPGVSRDSSTSYTETKDPSSGQEVATPPVPQLQVCEPK
ERTSTSSSSVQARRLSQPESAEKHVTQRLQPERGSPKKCEAEEAEPPAATQPQTSETQTSHLPESERIHHTVEEKEEVTM
DTSENRPENDVPEPPMPIADQVSNDDRPEGSVEDEEKKESSLPKSFKRKISVVSATKGVPAGNSDTEGGQPGRKRRWGAS
TATTQKKPSISITTESLKSLIPDIKPLAGQEAVVDLHADDSRISEDETERNGDDGTHDKGLKICRTVTQVVPAEGQENGQ
REEEEEEKEPEAEPPVPPQVSVEVALPPPAEHEVKKVTLGDTLTRRSISQQKSGVSITIDDPVRTAQVPSPPRGKISNIV
HISNLVRPFTLGQLKELLGRTGTLVEEAFWIDKIKSHCFVTYSTVEEAVATRTALHGVKWPQSNPKFLCADYAEQDELDY
HRGLLVDRPSETKTEEQGIPRPLHPPPPPPVQPPQHPRAEQREQERAVREQWAEREREMERRERTRSEREWDRDKVREGP
RSRSRSRDRRRKERAKSKEKKSEKKEKAQEEPPAKLLDDLFRKTKAAPCIYWLPLTDSQIVQKEAERAERAKEREKRRKE
QEEEEQKEREKEAERERNRQLEREKRREHSRERDRERERERERDRGDRDRDRERDRERGRERDRRDTKRHSRSRSRSTPV
RDRGGRR'''.replace('\n', '')
predator_exp = '''
------------------------------------hhhhh--------------------hhhhh--------------
------------------------------------------------------------------hhhhhh--------
-eeeeee-------hhhhhhhh----eeee--------------eeeeee----hhhhhhhhh--------hhhhhhh--
-----------------------------------------------------------------
'''.replace('\n', '')
mlr_exp = '''-------eeeeeeee----------------------------------------------------
--------------------------------------------------------------------------------
---hh----------eeeee------hhhhhhhhhhh--eeeee-------------eeeeee---hhhhhhhhhh----
----heeehhhh------------------------------------------------------------------
'''.replace('\n', '')
hnn_exp='''-----h---eeeeee------------------------------------------------------
-------------------------------------------------------------------------------h
hh-----------eeeee------hhhhhhhhhh---eeee---hh---------eeeeee----hhhhhhhhh------
--eeeeeeeh------------------------------------------------------------------
'''.replace('\n', '')
dpm_exp='''--e--e----ee-e-tttttt-----t---t-ttt-t-t-t-t-t---t-ht-ttt-t-t-hhhhhttt
ttttt-tttttt-ttttttttttttttttt-tttttttttt-tttttttttt-t-ttt-hhhttthht-t---tt---hh
hhhh-t-t-t--eeeeeeeeeeee--hehheeeee--eheehh-ehhhh--tt--eheehhht--hhhhhhhhht-----
-hheehehehh----------t---hh-------h-t---hhhht-ttt--t-ehhht-t-t-hhhhttttt----
'''.replace('\n', '')
coils_exp = '''-----------------------------------------------------------------
--------------------------------------------------------------------------------
-------------------------------------------------555555555555555555555----------
----------------000000000000000000000000000000000000000000000009999999-----
'''.replace('\n', '')
iupred_exp = '''
dddddddddddddddddddd--d--------------------ddd----------------------------------
--------------------------ddd--d-dddddddddddd------dddddddddddddddddddddd
'''.replace('\n', '')
globplot_exp = '''-----------------dddddddddddddddddddddddddddddddddddddd-------
-ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
dddddddddddddddddddddGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
GGGGGGGGGGGGGGGGGdddddddddddddddddddddddddddddddddddd-----dddddddddddddddddddddd
d--'''.replace('\n', '')
nls_exp = '''-------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
---------------------------------------------------------777777777777777777777--
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--------------------------------------------------------555555555555555555555555
555555555-----------------------------------------------------------------
'''.replace('\n', '')

if __name__ == '__main__':
    unittest.main()