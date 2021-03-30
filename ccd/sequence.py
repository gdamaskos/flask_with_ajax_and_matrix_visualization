import logging, re
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, generic_dna
from Bio.SeqUtils.ProtParam import ProteinAnalysis as PA
from numpy import append


_log = logging.getLogger(__name__)


def dna_validator(dna_seq):
    """Take a DNA sequence and check its validity as dna.
    Return False if invalid and True if valid.
    This validation is mainly for security reasons since another validation takes place in client's side. 
    """
    # Validate against an alphabet
    alphabet = 'ACTG'
    for base in dna_seq:
        if base not in alphabet:
            return False
    # Validate length as a multiple of 3
    if len(dna_seq) % 3 != 0:
        return False
    # Check for internal stop
    for base1 in range(0, len(dna_seq) - 3, 3):
        if dna_seq[base1: base1 + 3] == 'TAG' or dna_seq[base1: base1 + 3] == 'TGA' or dna_seq[base1: base1 + 3] == 'TAA':
            return False
    return True


def over_validator(overhang):
    """Take a overhang sequence and check its validity.
    Return False if invalid and True if valid.
    This validation is mainly for security reasons since another validation takes place in client's side. 
    """
    # Validate as character
    p = re.compile(r'^[a-zA-Z]+$')
    searchObj = p.search(overhang)
    if searchObj:
        return True
    return False

def translate(dna_orf):
    """Translate a DNA sequence to an amino acid sequence.
    Take an input string dna_orf and use the IUPAC alphabet for translation.
    Return a biopython Seq object.
    """
    _log.debug('Translating the input DNA open reading frame..')
    dna_orf = Seq(dna_orf, IUPAC.unambiguous_dna)
    aa_seq = dna_orf.translate()
    return aa_seq


def build_primers_tm(dna_orf, starts, stops, forward_overhang, reverse_overhang, user_tm):
    """Build primers using melting temperature method and validate starts/stops.
    Return forward and reverse primers, possible errors (like index error), valid starts and stops.
    """
    # starts part
    starts_primers = []
    valid_starts = []
    errors = []
    for start in starts:
        # keep in mind that start number/index begin from 1
        primer_name = 'Fw_{}'.format(int(start))
        try:
            dna_start_mapped = (int(start)-1)*3    
            curr_base = dna_start_mapped
            length = 0
            curr_tm = 0
            gc_no = 0
            while curr_tm < user_tm:
                if dna_orf[curr_base] == 'G' or dna_orf[curr_base] == 'C':
                    gc_no += 1
                length += 1
                curr_base += 1
                curr_tm = 64.9 + 41 * (gc_no - 16.4) / length  
            primer_sequence = '{}{}'.format(forward_overhang, dna_orf[dna_start_mapped : dna_start_mapped + length])
            starts_primers.append([primer_name, primer_sequence])
            valid_starts.append(start)
        except IndexError:
            errors.append([primer_name, ' Not enough bases in that direction'])
    # stops part
    stops_primers = []
    valid_stops = []
    reverse_stops = []
    reverse_complementary_dna = Seq(dna_orf, generic_dna)
    reverse_complementary_dna = str(reverse_complementary_dna.reverse_complement())
    reverse_stops[:] = [len(reverse_complementary_dna) / 3 - int(x) for x in stops]
    for stop in reverse_stops:
        primer_name = 'Rv_{}'.format(len(reverse_complementary_dna) / 3 - int(stop)) #Stop Number/Index begin from 1
        try:
            dna_stop_mapped = int(stop)*3    
            curr_base = dna_stop_mapped
            length = 0
            curr_tm = 0
            gc_no = 0
            while curr_tm < user_tm:
                if reverse_complementary_dna[curr_base] == 'G' or reverse_complementary_dna[curr_base] == 'C':
                    gc_no += 1
                length += 1
                curr_base += 1
                curr_tm = 64.9 + 41 * (gc_no - 16.4) / length  
            primer_sequence = '{}{}'.format(reverse_overhang, reverse_complementary_dna[dna_stop_mapped : dna_stop_mapped + length])
            stops_primers.append([primer_name, primer_sequence])
            valid_stops.append(len(reverse_complementary_dna) / 3 - int(stop))
        except IndexError:
            errors.append([primer_name, ' Not enough bases in that direction'])
    return starts_primers + stops_primers, errors, valid_starts, valid_stops


def build_primers_bs(dna_orf, starts, stops, forward_overhang, reverse_overhang, bs_no):
    """Build primers using given bases number method and validate starts/stops.
    Return forward and reverse primers, possible errors (like index error), valid starts and stops.
    """
    # starts part
    starts_primers = []
    valid_starts = []
    errors = []
    for start in starts:
        primer_name = 'Fw_{}'.format(int(start)) #Start Number/Index begin from 1
        try:
            dna_start_mapped = (int(start)-1)*3
            dna_orf[dna_start_mapped + bs_no] # here it throws index error
            primer_sequence = '{}{}'.format(forward_overhang, 
                                            dna_orf[dna_start_mapped : dna_start_mapped + bs_no])
            starts_primers.append([primer_name, primer_sequence])
            valid_starts.append(start)
        except IndexError:
            errors.append([primer_name, ' Not enough bases in that direction'])
    stops_primers = []
    reverse_stops = []
    valid_stops = []
    reverse_complementary_dna = Seq(dna_orf, generic_dna)
    reverse_complementary_dna = str(reverse_complementary_dna.reverse_complement())
    reverse_stops[:] = [len(reverse_complementary_dna) / 3 - int(x) for x in stops]
    for stop in reverse_stops:
        primer_name = 'Rv_{}'.format(len(reverse_complementary_dna) / 3 - int(stop)) #Stop Number/Index begin from 1
        try:
            dna_stop_mapped = int(stop)*3
            reverse_complementary_dna[dna_stop_mapped + bs_no]
            primer_sequence = '{}{}'.format(reverse_overhang, reverse_complementary_dna[dna_stop_mapped : dna_stop_mapped + bs_no])
            stops_primers.append([primer_name, primer_sequence])
            valid_stops.append(len(reverse_complementary_dna) / 3 - int(stop))
        except IndexError:
            errors.append([primer_name, ' Not enough bases in that direction'])
    return starts_primers + stops_primers, errors, valid_starts, valid_stops

def protparams(aa_seq, vstarts, vstops):
    """Compute a set of parameters for a polypepeptide,
    which would helps assess the potenial of this peptide as a crystalization candidate.
    """
    MWs = []
    pIs = []
    epsilons = []
    for start in vstarts:
        for stop in vstops:
            if int(start) < int(stop):
                params = PA(aa_seq[int(start) : int(stop)]) # works with string or Seq objects
                MW = params.molecular_weight()
                MW = round(MW/1000, 1) # in kiloDalton, rounded to 1 decimal
                pI = round(params.isoelectric_point(), 1)
                # To calculate the epsilon, we use this formula from protparam (web.expasy.org/protparam)
                # Epsilon (Prot) = N(Tyr)*Ext(Tyr) + N(Trp)*Ext(Trp) + N(Cystine)*Ext(Cystine) / MW in Dalton
                aa_dict = params.count_amino_acids() # returns a dict {'aa' : count } where aa is one letter code for the aminoacid
                epsilon = round((aa_dict['Y']*1490 + aa_dict['W']*5500 + aa_dict['C']*125)/(MW*1000), 2)
                MWs.append(MW)
                pIs.append(pI)
                epsilons.append(epsilon)
    return MWs, pIs, epsilons
