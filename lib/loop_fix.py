# standard library
import io
import sys
import re
from typing import Dict, Literal
import os
import time
from math import isnan
import gzip

# third-party libraries
from liftover import ChainFile as get_lifter_from_ChainFile # type: ignore # pylance mistakenly doesn't recognize ChainFile

# local
from math_utils import normal_p_area_two_tailed, normal_z_score_two_tailed
from standard_column_order import STANDARD_COLUMN_ORDER
from validate_utils import read_report_from_dir
from env import GWASSS_BUILD_NUMBER_ENV, get_build, set_build


# # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                 #
#                      INPUT                      #
#                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # #

if len(sys.argv) < 8:  # the very first 0th arg is the name of this script
    print("ERROR: you should specify args:")
    print("  #1 GWAS summary statistics file in the internal \"standard\" tsv format")
    print("  #2 directory with the report about the GWAS summary statistics file")
    print("  #3 OUTPUT: filename for GWAS summary statistics with fixes")
    print("  #4 preprocessed dbSNP1 file, or \"None\"")
    print("  #5 preprocessed dbSNP2 file, or \"None\"")
    print("  #6 chain file for liftover from build 36 or 37 to build 38, or \"None\"")
    print("  #7 frequency database slug (e.g.: \"GnomAD\", \"dbGaP_PopFreq\", \"TOMMO\"), or \"None\"")
    print("  #8 (optional) Either \"rsID\" or \"ChrBP\". Denotes the sorting of the input GWAS SS file")
    exit(1)

# GWAS_FILE has to be in the internal "standard" tsv format
GWAS_FILE = sys.argv[1]
REPORT_DIR = sys.argv[2]
OUTPUT_GWAS_FILE = sys.argv[3]
SNPs_FILE = sys.argv[4]
SNPs_rsID_FILE = sys.argv[5]
CHAIN_FILE = sys.argv[6]
FREQ_DATABASE_SLUG = sys.argv[7].lower() if sys.argv[7] != 'None' else None

GWAS_SORTING: Literal[None, 'rsID', 'ChrBP'] = None
if len(sys.argv) > 8:
    GWAS_SORTING = sys.argv[8] if sys.argv[8] in ('rsID', 'ChrBP') else None # type: ignore # pylance doesn't collapse types properly atm



def file_exists(path: str):
    return os.path.isfile(path)

if not file_exists(GWAS_FILE):
    print(f"ERROR: provided gwas file doesn't exist: {GWAS_FILE}")
    exit(1)


GWAS_FILE_o = open(GWAS_FILE, 'r')
OUTPUT_GWAS_FILE_o = open(OUTPUT_GWAS_FILE, 'w')
line_i=0



cols_i: Dict[str, int] = {STANDARD_COLUMN_ORDER[i]:i for i in range(len(STANDARD_COLUMN_ORDER))}


# # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                 #
#                    CONSTANTS                    #
#                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # #

NUCLEOTIDES = ['a', 't', 'c', 'g']
NO_NUCLEOTIDE = '.'

ALLOW_MULTI_NUCLEOTIDE_POLYMORPHISMS = True

CATEGORY_CHR = [
'1', '01', '2', '02', '3', '03', '4', '04', '5', '05', '6', '06', '7', '07', '8', '08', '9', '09',
'10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
'21', '22', '23', 'X', 'x', 'Y', 'y', 'M', 'm']

class kukdefaultdict(dict):
    def __missing__(self, key):
        return key
CHR_ORDER = kukdefaultdict() # if unknown key was passed, returns the key itself
CHR_ORDER['1']  = 1
CHR_ORDER['01'] = 1 #
CHR_ORDER['2']  = 2
CHR_ORDER['02'] = 2 #
CHR_ORDER['3']  = 3
CHR_ORDER['03'] = 3 #
CHR_ORDER['4']  = 4
CHR_ORDER['04'] = 4 #
CHR_ORDER['5']  = 5
CHR_ORDER['05'] = 5 #
CHR_ORDER['6']  = 6
CHR_ORDER['06'] = 6 #
CHR_ORDER['7']  = 7
CHR_ORDER['07'] = 7 #
CHR_ORDER['8']  = 8
CHR_ORDER['08'] = 8 #
CHR_ORDER['9']  = 9
CHR_ORDER['09'] = 9 #
CHR_ORDER['10'] = 10
CHR_ORDER['11'] = 11
CHR_ORDER['12'] = 12
CHR_ORDER['13'] = 13
CHR_ORDER['14'] = 14
CHR_ORDER['15'] = 15
CHR_ORDER['16'] = 16
CHR_ORDER['17'] = 17
CHR_ORDER['18'] = 18
CHR_ORDER['19'] = 19
CHR_ORDER['20'] = 20
CHR_ORDER['21'] = 21
CHR_ORDER['22'] = 22
CHR_ORDER['23'] = 23

CHR_ORDER['X'] =  25
CHR_ORDER['x'] =  25
CHR_ORDER['Y'] =  26
CHR_ORDER['y'] =  26
CHR_ORDER['M'] =  27
CHR_ORDER['m'] =  27

CHR_LIFTOVER = kukdefaultdict() 
CHR_LIFTOVER['1']  = '1'
CHR_LIFTOVER['01'] = '1' #
CHR_LIFTOVER['2']  = '2'
CHR_LIFTOVER['02'] = '2' #
CHR_LIFTOVER['3']  = '3'
CHR_LIFTOVER['03'] = '3' #
CHR_LIFTOVER['4']  = '4'
CHR_LIFTOVER['04'] = '4' #
CHR_LIFTOVER['5']  = '5'
CHR_LIFTOVER['05'] = '5' #
CHR_LIFTOVER['6']  = '6'
CHR_LIFTOVER['06'] = '6' #
CHR_LIFTOVER['7']  = '7'
CHR_LIFTOVER['07'] = '7' #
CHR_LIFTOVER['8']  = '8'
CHR_LIFTOVER['08'] = '8' #
CHR_LIFTOVER['9']  = '9'
CHR_LIFTOVER['09'] = '9' #
CHR_LIFTOVER['10'] = '10'
CHR_LIFTOVER['11'] = '11'
CHR_LIFTOVER['12'] = '12'
CHR_LIFTOVER['13'] = '13'
CHR_LIFTOVER['14'] = '14'
CHR_LIFTOVER['15'] = '15'
CHR_LIFTOVER['16'] = '16'
CHR_LIFTOVER['17'] = '17'
CHR_LIFTOVER['18'] = '18'
CHR_LIFTOVER['19'] = '19'
CHR_LIFTOVER['20'] = '20'
CHR_LIFTOVER['21'] = '21'
CHR_LIFTOVER['22'] = '22'
CHR_LIFTOVER['23'] = '23'

CHR_LIFTOVER['X'] =  'X'
CHR_LIFTOVER['x'] =  'X'
CHR_LIFTOVER['Y'] =  'Y'
CHR_LIFTOVER['y'] =  'Y'
CHR_LIFTOVER['M'] =  'M'
CHR_LIFTOVER['m'] =  'M'


# # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                 #
#                    FUNCTIONS                    #
#                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # #

def copy_line(line_i):
    OUTPUT_GWAS_FILE_o.write(GWAS_FILE_o.readline())
    return line_i + 1


# # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                 #
#         FUNCTIONS THAT RUN MANY TIMES           #
#                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # #

##### read/write #####

def get_next_line_in_GWASSS():
    line = GWAS_FILE_o.readline()
    if line == '': raise EOFError('attempt to read beyond the end of GWAS SS file')
    return line.split("\t")

def write_line_to_GWASSS(fields):
    OUTPUT_GWAS_FILE_o.write("\t".join(fields))


def read_dbSNP1_data_row(FILE_o: io.TextIOWrapper):
    """Reads a row from the preprocessed dbSNP file 1"""
    line = FILE_o.readline()
    words = line.split()
    return (
        words[0], # Chr
        int(words[1]), # BP
        words[2], # rsID
        words[3], # REF
        words[4], # ALT
        words[5], # freq
    )

def read_dbSNP2_data_row(FILE_o: io.TextIOWrapper):
    """Reads a row from the preprocessed dbSNP file 2 (which is sorted by rsID)"""
    line = FILE_o.readline()
    words = line.split()
    return (
        words[1], # Chr
        words[2], # BP
        words[0], # rsID
        words[3], # REF
        words[4], # ALT
        words[5], # freq
    )

def gt(val1, val2):
    """
    A safe "greater than" operator. Accepts int and str for both args.

    It has the following 3 important features:
      - if both values are numbers (like chromosome order number), it compares numerically
      - if both values are strings (like unrecognized chromosome numbers), compares alphabetically
      - if the second value is string (like unrecognized chromosome number), and
        the first value is number (like chromosome order number),
        then it assumes the number is always no bigger than the string.
        This is because (here) "unrecognized" chromosome numbers in dbSNP always go after
        the common chromosome numbers (known in the CHR_ORDER dictionary).

    With the assumption that the case where the first value is string
    and the second value is number doesn't occur,
    this order relation defines a *totally ordered set* which is mathematically defined as:
    1. First go chromosome numbers defined as keys in CHR_ORDER dictionary,
    starting with '1' = '01', ordered by the corresponding values in the dictionary.
    2. After that goes the set of all other strings
       (the complement set of strings to the set of known chromosome numbers)
       in the ascending alphabetical order

    This matches the way GWAS SS is sorted by Chr and BP here.
    """
    try:
        return val1 > val2
    except:
        return False


##### Math & Stats functions #####

"""
These functions resolve each of the three values (p-value, beta, standard error) from the other two

                    s = b/z,
where:
    s - standard error,
    b - beta-value,
    z - the normal z-score that corresponds to the p-value

e.g. for the two-tailed test (which is used in GWAS):
                z = qnorm(1 - p/2),
where:
    qnorm - inverse cumulative function for normal distribution
"""
def get_StdErr_from_beta_pval(beta, p):
    z = normal_z_score_two_tailed(p)
    return abs(beta)/z if z != 0 else 'nan'

def get_beta_from_StdErr_pval(se, p):
    z = normal_z_score_two_tailed(p)
    return se*z

def get_pval_from_beta_StdErr(beta, se):
    z = abs(beta)/se if se != 0 else 'nan'
    return normal_p_area_two_tailed(z)


##### Field validators #####
"""
These boolean functions:
    accept the list of fields read from a line;
    return a boolean value answering whether or not a particular field in a row is valid or not
"""

def is_valid_rsID(fields):
    try:
        rsid = fields[cols_i["rsID"]]
        if not re.match("^rs\d+$", rsid):
            return False
    except:
        return False
    return True

def is_valid_Chr(fields):
    try:
        chr = fields[cols_i["Chr"]]
        if chr not in CATEGORY_CHR:
            return False
    except:
        return False
    return True

def is_valid_BP(fields):
    try:
        bp = int(float(fields[cols_i["BP"]])) # using float allows sci notation string
        if bp < 0:
            return False
    except:
        return False
    return True

def is_valid_EA_allowMNP(fields):
    try:
        ea = fields[cols_i["EA"]]
        if ea == '':
            return False
        if ea == NO_NUCLEOTIDE:
            return True
        for char in ea.lower():
            if char not in NUCLEOTIDES:
                return False
    except:
        return False
    return True

def is_valid_OA_allowMNP(fields):
    try:
        oa = fields[cols_i["OA"]]
        if oa == '':
            return False
        if oa == NO_NUCLEOTIDE:
            return True
        for char in oa.lower():
            if char not in NUCLEOTIDES:
                return False
    except:
        return False
    return True

def is_valid_EA_dontallowMNP(fields):
    try:
        ea = fields[cols_i["EA"]]
        if ea == '':
            return False
        if ea == NO_NUCLEOTIDE:
            return True
        if ea.lower() not in NUCLEOTIDES:
            return False
    except:
        return False
    return True

def is_valid_OA_dontallowMNP(fields):
    try:
        oa = fields[cols_i["OA"]]
        if oa == '':
            return False
        if oa == NO_NUCLEOTIDE:
            return True
        if oa.lower() not in NUCLEOTIDES:
            return False
    except:
        return False
    return True

is_valid_EA = None
is_valid_OA = None
if ALLOW_MULTI_NUCLEOTIDE_POLYMORPHISMS:
    is_valid_EA = is_valid_EA_allowMNP
    is_valid_OA = is_valid_OA_allowMNP
else:
    is_valid_EA = is_valid_EA_dontallowMNP
    is_valid_OA = is_valid_OA_dontallowMNP


def is_valid_EAF(fields):
    try:
        af = float(fields[cols_i["EAF"]])
        if not (0 <= af <= 1):
            return False
    except:
        return False
    return True

def is_valid_SE(fields):
    try:
        se = float(fields[cols_i["SE"]])
        if isnan(se):
            return False
    except:
        return False
    return True

def is_valid_beta(fields):
    try:
        beta = float(fields[cols_i["beta"]])
        if isnan(beta):
            return False
    except:
        return False
    return True

def is_valid_pval(fields):
    try:
        pval = float(fields[cols_i["pval"]])
        if isnan(pval) or not (0 <= pval <= 1):
            return False
    except:
        return False
    return True



##### RESOLVERS helpers #####
"""
These void functions accept the list of fields read from a line and may mutate it

`fields`: list[str]
    list of entries in a row
"""

def resolve_allele(fields, REF, ALT):
    """
    Runs iff exactly one allele entry is missing.
    Depending on the REF and ALT (from the dbSNP), decides which one is the allele that's present,
    and restores the other one in accord.
    """
    if   is_valid_EA(fields) and not is_valid_OA(fields):
        MA = 'OA'  # missing allele
        PA = 'EA'  # present allele
    elif is_valid_OA(fields) and not is_valid_EA(fields):
        MA = 'EA'  # missing allele
        PA = 'OA'  # present allele
    else:
        return

    PA_val = fields[cols_i[PA]]

    if PA_val == REF:
        fields[cols_i[MA]] = ALT.split(',')[0]
        return
    for a in ALT.split(','):
        if PA_val == a:
            fields[cols_i[MA]] = REF
            return


def resolve_EAF(fields, REF, ALT, SNP_freq_field):
    """
    1. Tries to find the allele by exact match in REF or ALT;
    2. Checks if frequency entry from the given database slug is present;
    3. Takes the allele frequency for the corresponding allele
    """

    """
    example value for `SNP_freq_field`:
        "freq=1000Genomes:0.9988,.,0.001198|GnomAD:0.9943,0.005747,."

    This describes allele frequencies for 3 alleles,
    from 2 databases: 1000Genomes, GnomAD

    Having 3 alleles, means an allele in REF field and two alleles in ALT separated with comma.
    For this example, REF and ALT values are the following:
        "TA"
        "T,TAA"
    and `FREQ_DATABASE_SLUG`:
        "GnomAD"
    """

    if not is_valid_EAF(fields) and is_valid_EA(fields):
        try:
            freqs = SNP_freq_field.replace('freq=','').replace('|',':').split(':') # ["1000Genomes", "0.9988,.,0.001198", "GnomAD", "0.9943,0.005747,."]
            freqs = [f.lower() for f in freqs]  # ["1000genomes", "0.9988,.,0.001198", "gnomad", "0.9943,0.005747,."]
            alleles = (REF+','+ALT).split(',') # ["TA", "T", "TAA"]
            EA = fields[cols_i['EA']] # "T"

            the_freq_db_i = freqs.index(FREQ_DATABASE_SLUG) # "2"
            SNP_freqs = freqs[the_freq_db_i+1].split(',') # ["0.9943", "0.005747", "."]
            allele_i = alleles.index(EA) # 1
            fields[cols_i['EAF']] = SNP_freqs[allele_i]  # "0.005747"
        except:
            fields[cols_i['EAF']] = '.'


##### RESOLVERS #####
"""
These void functions accept the list of fields read from a line and may mutate it

`fields`: list[str]
    list of entries in a row
"""

def resolve_build38(fields, converter):
    """
    Will use the input converter dictionary to liftover
    from the build specified by the user to build38 (with 'chr' prefix)
    """
    if is_valid_Chr(fields) and is_valid_BP(fields):
        chr_gwas = CHR_LIFTOVER[fields[cols_i['Chr']]]
        bp_gwas  = int(float(fields[cols_i['BP']])) # using float allows sci notation string
        try:
            new_chr, new_bp, _ = converter[chr_gwas][bp_gwas][0]
            fields[cols_i["Chr"]] = new_chr.replace('chr', '')
            fields[cols_i["BP"]] = str(new_bp)
        # if it can't liftover
        except:
            fields[cols_i["Chr"]] = '.'
            fields[cols_i["BP"]] = '.'


def resolve_rsID(fields, SNPs_FILE_o):
    """
    Loops through the SNPs file entries until it finds the current SNP in GWAS SS file
    Current SNP is defined by Chr and BP from the passed `fields`, which is one row of GWAS SS
    Therefore, rsID is restored by Chr and BP, which won't always work for biallelic sites

    Assumes given GWAS SS file is sorted by Chr and BP, in the same way SNPs file is.

    `SNPs_FILE_o`
        opened SNPs file object
    """
    if is_valid_Chr(fields) and is_valid_BP(fields) and not all([
        is_valid_rsID(fields),
        is_valid_OA(fields), is_valid_EA(fields),
        is_valid_EAF(fields),
    ]):
        try:
            chr_gwas = fields[cols_i['Chr']]
            bp_gwas  = int(float(fields[cols_i['BP']]))

            while True:
                chr_snps, bp_snps, rsid, ref, alt, freq = read_dbSNP1_data_row(SNPs_FILE_o)
                # SNPs_FILE_line_i += 1

                if CHR_ORDER[chr_gwas] == CHR_ORDER[chr_snps]:
                    if bp_snps < bp_gwas:
                        continue
                    elif bp_gwas == bp_snps:
                        fields[cols_i['rsID']] = rsid
                        resolve_allele(fields, ref, alt)
                        resolve_EAF(fields, ref, alt, freq)
                        break # after this a new line of GWAS SS should be read and index incremented
                    else: #bp_snps > bp_gwas:
                        fields[cols_i['rsID']] = '.'
                        break # after this a new line of GWAS SS should be read and index incremented
                elif gt(CHR_ORDER[chr_snps], CHR_ORDER[chr_gwas]):
                    fields[cols_i['rsID']] = '.'
                    break # after this a new line of GWAS SS should be read and index incremented

        except Exception as e:
            if isinstance(e, IndexError) or isinstance(e, EOFError):
                # it reached the end of an either file
                pass
            else:
                # print(f'An error occured on line {SNPs_FILE_line_i} of the SNPs file (see below)')
                print(f'An error occured while looping through the SNPs file (see below)')
                raise e


def resolve_ChrBP(fields, SNPs_rsID_FILE_o):
    """
    Loops through the SNPs file entries until it finds the current locus
    Current locus is defined by rsID from the passed `fields`, which is one row of GWAS SS

    Assumes given GWAS SS file is sorted by rsID, in the same way this processed SNPs file is.

    `SNPs_rsID_FILE_o`
        opened SNPs file object
    """
    if is_valid_rsID(fields) and not all([
        is_valid_Chr(fields), is_valid_BP(fields),
        is_valid_OA(fields), is_valid_EA(fields),
        is_valid_EAF(fields),
    ]):
        try:
            rsID_gwas = fields[cols_i['rsID']]

            while True:
                chr_snps, bp_snps, rsid, ref, alt, freq = read_dbSNP2_data_row(SNPs_rsID_FILE_o)
                # SNPs_FILE_line_i += 1

                if rsid < rsID_gwas:
                    continue
                elif rsID_gwas == rsid:
                    fields[cols_i['Chr']] = chr_snps
                    fields[cols_i['BP']] = bp_snps
                    resolve_allele(fields, ref, alt)
                    resolve_EAF(fields, ref, alt, freq)
                    break # after this a new line of GWAS SS should be read and index incremented
                else: #rsid > bp_gwas:
                    fields[cols_i['Chr']] = '.'
                    fields[cols_i['BP']] = '.'
                    break # after this a new line of GWAS SS should be read and index incremented

        except Exception as e:
            if isinstance(e, IndexError) or isinstance(e, EOFError):
                # it reached the end of an either file
                pass
            else:
                # print(f'An error occured on line {SNPs_FILE_line_i} of the SNPs file (see below)')
                print(f'An error occured while looping through the SNPs file (see below)')
                raise e


def resolve_SE(fields):
    if not is_valid_SE(fields) and is_valid_beta(fields) and is_valid_pval(fields):
        fields[cols_i["SE"]] = str(get_StdErr_from_beta_pval(
            float(fields[cols_i["beta"]]), float(fields[cols_i["pval"]])
        ))
def resolve_beta(fields):
    if not is_valid_beta(fields) and is_valid_SE(fields) and is_valid_pval(fields):
        fields[cols_i["beta"]] = str(get_beta_from_StdErr_pval(
            float(fields[cols_i["SE"]]), float(fields[cols_i["pval"]])
        ))
def resolve_pval(fields):
    if not is_valid_pval(fields) and is_valid_beta(fields) and is_valid_SE(fields):
        fields[cols_i["pval"]] = str(get_pval_from_beta_StdErr(
            float(fields[cols_i["beta"]]), float(fields[cols_i["SE"]])
        ))


##### FULL RESOLVER #####
"""
This function will be called for each row of the input GWAS SS file
It calls resolvers one by one from the list of resolvers.
Each resolver attempts to resolve one or many values for the given row.

Each resolver has `fields: list[str]` as the first argument,
and may have other args defined as a list under the corresponding index in `resolvers_args` list
"""
resolvers = [] # list of functions
resolvers_args = [] # list of lists of arguments for these functions

def run_all(resolvers, fields, args):
    for res_i in range(len(resolvers)):
        resolvers[res_i](fields, *args[res_i])



# # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                 #
#                      MAIN                       #
#                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # #

MAIN_start_time = STEP1_start_time = time.time()

#
# STEP #1
#     Assemble the full resolver function in accord to present issues.
#
#     E.g.:
#       - if all values from a pacticular column are valid,
#         there's no need to add a resolver function for it.
#       - if all values from a particular column are invalid,
#         there's no need to add a resolver for another column that depends on the first,
#         as it would always run uselessly
#
#

issues, total_entries = read_report_from_dir(REPORT_DIR)

DOING_LIFTOVER: bool = False

current_build = get_build()
converter = None
if current_build != 'hg38' and file_exists(CHAIN_FILE):
    DOING_LIFTOVER = True

if DOING_LIFTOVER:
    converter = get_lifter_from_ChainFile(CHAIN_FILE, current_build, 'hg38')
    set_build('hg38')
    resolvers.append(resolve_build38)
    resolvers_args.append([converter])


if not DOING_LIFTOVER and GWAS_SORTING == 'rsID' and (issues['Chr'] or issues['BP'] or issues['OA'] or issues['EA'] or issues['EAF']) and file_exists(SNPs_rsID_FILE):
    """
    This ChrBP resolver assumes GWAS SS file is sorted by rsID
    """
    # open files here
    SNPs_rsID_FILE_o_gz: io.RawIOBase = gzip.open(SNPs_rsID_FILE, 'r')  # type: ignore # GzipFile and RawIOBase _are_ in fact compatible
    SNPs_rsID_FILE_o = io.TextIOWrapper(io.BufferedReader(SNPs_rsID_FILE_o_gz))

    resolvers.append(resolve_ChrBP)
    resolvers_args.append([SNPs_rsID_FILE_o])


if not DOING_LIFTOVER and GWAS_SORTING == 'ChrBP' and (issues['rsID'] or issues['OA'] or issues['EA'] or issues['EAF']) and file_exists(SNPs_FILE):
    """
    These resolvers assumes GWAS SS file is sorted by Chr and BP in accord to the SNPs file
    """
    SNPs_FILE_o_gz: io.RawIOBase = gzip.open(SNPs_FILE, 'r')  # type: ignore # GzipFile and RawIOBase _are_ in fact compatible
    SNPs_FILE_o = io.TextIOWrapper(io.BufferedReader(SNPs_FILE_o_gz))

    resolvers.append(resolve_rsID)
    resolvers_args.append([SNPs_FILE_o])


if issues['SE'] and issues['beta']<total_entries and issues['pval']<total_entries:
    resolvers.append(resolve_SE)
    resolvers_args.append([])

if issues['beta'] and issues['SE']<total_entries and issues['pval']<total_entries:
    resolvers.append(resolve_beta)
    resolvers_args.append([])

if issues['pval'] and issues['beta']<total_entries and issues['SE']<total_entries:
    resolvers.append(resolve_pval)
    resolvers_args.append([])


#
# STEP #2
#     Loop through the GWAS SS file,
#     For each row, run all assembled resolvers, and FINALLY save the row to the output file
#

# copy the first line that is the header
line_i = copy_line(line_i)

try:
    while True:
        fields = get_next_line_in_GWASSS()
        run_all(resolvers, fields, resolvers_args)
        write_line_to_GWASSS(fields)

except Exception as e:
    if isinstance(e, IndexError) or isinstance(e, EOFError):
        # it reached the end of the file
        pass
    else:
        print(f'An error occured on line {line_i} of the GWAS SS file (see below)')
        raise e
 

GWAS_FILE_o.close()
OUTPUT_GWAS_FILE_o.close()

