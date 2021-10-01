# standard library
import io
import sys
import re
from typing import Dict
import os
import time
from math import isnan
import gzip

# local
from math_utils import normal_p_area_two_tailed, normal_z_score, normal_z_score_two_tailed
from standard_column_order import STANDARD_COLUMN_ORDER
from validate_utils import read_report_from_dir


# # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                 #
#                      INPUT                      #
#                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # #

if len(sys.argv) < 5:
    print("ERROR: you should specify args:")
    print("  #1 GWAS summary statistics file in the internal \"standard\" tsv format")
    print("  #2 directory with the report about the GWAS summary statistics file")
    print("  #3 output: filename for GWAS summary statistics with fixes")
    print("  #4 dbSNP file path")
    exit(1)

# GWAS_FILE has to be in the internal "standard" tsv format
GWAS_FILE = sys.argv[1]
REPORT_DIR = sys.argv[2]
OUTPUT_GWAS_FILE = sys.argv[3]
SNPs_FILE = sys.argv[4]

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

ALLOW_MULTI_NUCLEOTIDE_POLYMORPHISMS = True

CATEGORY_CHR = [
'1', '01', '2', '02', '3', '03', '4', '04', '5', '05', '6', '06', '7', '07', '8', '08', '9', '09',
'10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
'21', '22', '23', 'X', 'x', 'Y', 'y', 'M', 'm']
CHR_ORDER = {
    '1':  1,
    '01': 1,
    '2':  2,
    '02': 2,
    '3':  3,
    '03': 4,
    '4':  4,
    '04': 4,
    '5':  5,
    '05': 5,
    '6':  6,
    '06': 6,
    '7':  7,
    '07': 7,
    '8':  8,
    '08': 8,
    '9':  9,
    '09': 9,
    '10': 10,
    '11': 11,
    '12': 12,
    '13': 13,
    '14': 14,
    '15': 15,
    '16': 16,
    '17': 17,
    '18': 18,
    '19': 19,
    '20': 20,
    '21': 21,
    '22': 22,
    '23': 23,

    'X':  25,
    'x':  25,
    'Y':  26,
    'y':  26,
    'M':  27,
    'm':  27,
}



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


def read_dbSNPs_data_row(FILE_o: io.TextIOWrapper):
    line = FILE_o.readline()
    words = line.split()
    return (
        words[0][3:], # bc in SNPs_file it chromosome number is prepended with "chr"
        int(words[1]), # BP
        words[2] # rsID
    )


##### Math & Stats functions #####

def get_StdErr_from_beta_pval(beta, p):
    z = normal_z_score_two_tailed(p)
    return abs(beta)/z

def get_beta_from_StdErr_pval(se, p):
    z = normal_z_score_two_tailed(p)
    return se*z

def get_pval_from_beta_StdErr(beta, se):
    z = abs(beta)/se
    return normal_p_area_two_tailed(z)


##### Field validators #####
"""
These boolean functions accept the list of fields read from a line
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
        bp = int(fields[cols_i["BP"]])
        if bp < 0:
            return False
    except:
        return False
    return True

def is_valid_EA_allowMNP(fields):
    try:
        ea = fields[cols_i["EA"]]
        for char in ea.lower():
            if char not in NUCLEOTIDES:
                return False
    except:
        return False
    return True

def is_valid_OA_allowMNP(fields):
    try:
        oa = fields[cols_i["OA"]]
        for char in oa.lower():
            if char not in NUCLEOTIDES:
                return False
    except:
        return False
    return True

def is_valid_EA_dontallowMNP(fields):
    try:
        ea = fields[cols_i["EA"]]
        if ea.lower() not in NUCLEOTIDES:
            return False
    except:
        return False
    return True

def is_valid_OA_dontallowMNP(fields):
    try:
        oa = fields[cols_i["OA"]]
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


##### RESOLVERS #####
"""
These void functions accept the list of fields read from a line and may mutate it
"""

def resolve_rsID(fields, SNPs_FILE_o):
    """
    Loops through the SNPs file entries until it finds the current locus
    Current locus is defined by Chr and BP from the passed `fields`, which is one row of GWAS SS

    Assumes given GWAS SS file is sorted by Chr and BP, in the same way SNPs file is.

    `SNPs_FILE_o`
        opened SNPs file object
    """
    if not is_valid_rsID(fields) and is_valid_Chr(fields) and is_valid_BP(fields):
        try:
            chr_gwas = fields[cols_i['Chr']]
            bp_gwas  = int(fields[cols_i['BP']])

            while True:
                chr_snps, bp_snps, rsid = read_dbSNPs_data_row(SNPs_FILE_o)
                # SNPs_FILE_line_i += 1

                if CHR_ORDER[chr_gwas] == CHR_ORDER[chr_snps]:
                    if bp_snps < bp_gwas:
                        continue
                    elif bp_gwas == bp_snps:
                        fields[cols_i['rsID']] = rsid
                        break # after this a new line of GWAS SS should be read and index incremented
                    else: #bp_snps > bp_gwas:
                        fields[cols_i['rsID']] = '.'
                        break # after this a new line of GWAS SS should be read and index incremented
                elif CHR_ORDER[chr_snps] > CHR_ORDER[chr_gwas]:
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


def resolve_SE(fields):
    if not is_valid_SE(fields):
        fields[cols_i["SE"]] = str(get_StdErr_from_beta_pval(
            float(fields[cols_i["beta"]]), float(fields[cols_i["pval"]])
        ))
def resolve_beta(fields):
    if not is_valid_beta(fields):
        fields[cols_i["beta"]] = str(get_beta_from_StdErr_pval(
            float(fields[cols_i["SE"]]), float(fields[cols_i["pval"]])
        ))
def resolve_pval(fields):
    if not is_valid_pval(fields):
        fields[cols_i["pval"]] = str(get_pval_from_beta_StdErr(
            float(fields[cols_i["beta"]]), float(fields[cols_i["SE"]])
        ))


##### FULL RESOLVER #####
"""
This function will be called for each row.
It calls resolvers one by one from the list of resolvers (hope this helps)
"""
resolve_fields = None
resolvers = []
resolvers_args = []

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
#     Assemble the full resolver function in accord to present issues
#

issues = read_report_from_dir(REPORT_DIR)

if issues['rsID']:
    """
    This rsID resolver assumes GWAS SS file is sorted by Chr and BP in accord to the SNPs file
    """
    SNPs_FILE_o_gz: io.RawIOBase = gzip.open(SNPs_FILE, 'r')  # type: ignore # GzipFile and RawIOBase _are_ in fact compatible
    SNPs_FILE_o = io.TextIOWrapper(io.BufferedReader(SNPs_FILE_o_gz))
    # SNPs_FILE_line_i = 0

    """
    SNPs_FILE, being a VCF file, starts with comment lines at the beginning. Comments start with ##
    Then comes a header, starts with #
    after the header comes the dataframe table

    THIS codeblock skips all the comment lines,
    then reads the header line before going to the next loop.
    this prevents putting more conditions into the loop
    """
    SNPs_line = SNPs_FILE_o.readline()
    # SNPs_FILE_line_i += 1
    while SNPs_line.startswith('##'):
        SNPs_line = SNPs_FILE_o.readline()
        # SNPs_FILE_line_i += 1

    resolvers.append(resolve_rsID)
    resolvers_args.append([SNPs_FILE_o])


if issues['SE']:
    resolvers.append(resolve_SE)
    resolvers_args.append([])

if issues['beta']:
    resolvers.append(resolve_beta)
    resolvers_args.append([])

if issues['pval']:
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

