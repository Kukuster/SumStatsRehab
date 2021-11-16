# standard library
import sys
import os
from subprocess import call
import inspect
import time
from typing import Dict, Literal, Union
import json
import argparse
import pathlib

# local
from lib.validate_utils import read_report_from_dir
from lib.standard_column_order import STANDARD_COLUMN_ORDER
from lib.env import get_build, set_build




# # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                 #
#                    CONSTANTS                    #
#                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # #

# define paths to libs
the_dir = os.path.dirname(inspect.getfile(inspect.currentframe())) or os.getcwd()  # type: ignore
lib_dir = the_dir + "/lib"
prepare_two_dbSNPs = lib_dir+"/prepare_two_dbSNPs.py"
prepare_GWASSS_columns = lib_dir+"/prepare_GWASSS_columns.py"
validate_GWASSS_entries = lib_dir+"/validate_GWASSS_entries.py"
sort_GWASSS_by_ChrBP = lib_dir+"/sort_GWASSS_by_ChrBP.py"
sort_GWASSS_by_rsID = lib_dir+"/sort_GWASSS_by_rsID.py"
loop_fix = lib_dir+"/loop_fix.py"



# # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                 #
#                    FUNCTIONS                    #
#                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # #


def remove_last_ext(filename: str):
    return filename.rsplit(".", 1)[0] # passed 1 means do max 1 split; _rightmost_ splits first


def perc(x, total):
    if x/total < 0.0001:
        return "<0.01%"
    else:
        return str(round((x/total)*100, 2)) + "%"



# # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                 #
#                    COMMANDS                     #
#                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # #



def fix(INPUT_GWAS_FILE: str, OUTPUT_FILE: str, dbSNP_FILE: str, dbSNP2_FILE: str, CHAIN_FILE: str, FREQ_DATABASE_SLUG: str):

    ### PROCESS INPUT ###
    INPUT_GWAS_FILE = str(INPUT_GWAS_FILE)
    JSON_CONFIG = INPUT_GWAS_FILE + '.json'
    OUTPUT_FILE = str(OUTPUT_FILE)
    dbSNP_FILE = str(dbSNP_FILE)
    dbSNP2_FILE = str(dbSNP2_FILE)
    CHAIN_FILE = str(CHAIN_FILE)
    FREQ_DATABASE_SLUG = str(FREQ_DATABASE_SLUG)


    ### Set the environment varaible for the build ###
    config: Dict[str, Union[int, str]] = json.load(open(JSON_CONFIG,))
    try:
        input_build = config['build']
    except:
        raise ValueError('config file that corresponds to the GWAS SS file has to have "build" key set')
    set_build(input_build)
    input_build = get_build()
    print(f'input build: {input_build}')


    i_step = 0

    ##### 1 #####
    i_step += 1
    print(f'=== Step {i_step}: Format the GWAS SS file ===')
    start_time = time.time()

    INPUT_GWAS_FILE_standard = remove_last_ext(INPUT_GWAS_FILE) + "_standard.tsv"
    ec = call(["python3",
            prepare_GWASSS_columns,
            INPUT_GWAS_FILE,
            INPUT_GWAS_FILE_standard,
            ])
    print(f"  Step {i_step} finished in {(time.time() - start_time)} seconds\n")

    if ec != 0:
        print(f"ERROR: prepare_GWASSS_columns script finished with exit code: {ec}")
        exit(11)



    ##### 2 #####
    i_step += 1
    print(f'=== Step {i_step}: Validate entries in the formatted GWAS SS file and save the report ===')
    start_time = time.time()

    input_validation_report_dir = INPUT_GWAS_FILE + "_input-report"
    ec = call(["python3",
            validate_GWASSS_entries,
            INPUT_GWAS_FILE_standard,
            "standard",
            input_validation_report_dir,
            ])
    print(f"  Step {i_step} finished in {(time.time() - start_time)} seconds\n")

    if ec != 0:
        print(f"ERROR: validate_GWASSS_entries script finished with exit code: {ec}")
        exit(12)



    ##### 3 #####
    i_step += 1
    print(f'=== Step {i_step}: Analyze the report and prepare for REHAB ===')
    start_time = time.time()

    issues, total_entries = read_report_from_dir(input_validation_report_dir)

    required_sorting: bool = False
    required_liftover: bool = False
    sorted_by: Literal[None, 'rsID', 'ChrBP'] = None
    INPUT_GWAS_FILE_standard_sorted = remove_last_ext(INPUT_GWAS_FILE) + "_standard_sorted.tsv"
    INPUT_GWAS_FILE_standard_lifted = remove_last_ext(INPUT_GWAS_FILE) + "_standard_lifted.tsv"


    if get_build() != 'hg38' and CHAIN_FILE and CHAIN_FILE != "None" and issues['BP']<total_entries and issues['Chr']<total_entries:
        # if either BP or Chr is fully missing,
        # there's no need for liftover since those will be restored with dbSNPs in the target build
        required_liftover = True
        ec = call(["python3",
                loop_fix,
                INPUT_GWAS_FILE_standard,
                input_validation_report_dir,
                INPUT_GWAS_FILE_standard_lifted,
                dbSNP_FILE,
                dbSNP2_FILE,
                CHAIN_FILE,
                FREQ_DATABASE_SLUG if FREQ_DATABASE_SLUG else 'None',
                sorted_by if sorted_by else ''
                ])

        if ec != 0:
            print(f"ERROR: loop_fix (liftover) script finished with exit code: {ec}")
            exit(13)

        print("finished liftover to hg38 (saved report)")
        set_build('hg38')

        input_lifted_validation_report_dir = INPUT_GWAS_FILE + "_input-lifted-report"
        ec = call(["python3",
            validate_GWASSS_entries,
            INPUT_GWAS_FILE_standard_lifted,
            "standard",
            input_lifted_validation_report_dir,
            ])
        if ec != 0:
            print(f"ERROR: validate_GWASSS_entries script finished with exit code: {ec}")
            exit(13)
    
    if (issues['BP'] or issues['Chr']) and issues['rsID']<total_entries:
        # here if some alleles are invalid, they will be attempted to be restored by rsID, which is better then by ChrBP
        required_sorting = True
        sorted_by = 'rsID'
        for col in ('Chr', 'BP'):
            if issues[col]:
                print(f"{issues[col]}/{total_entries} entries are missing {col}")
        print(f"Going to sort the GWAS SS file by rsID")
        ec = call(["python3",
            sort_GWASSS_by_rsID,
            INPUT_GWAS_FILE_standard_lifted if required_liftover else INPUT_GWAS_FILE_standard,
            INPUT_GWAS_FILE_standard_sorted,
            ])
        if ec == 0:
            print(f"Sorted by rsID")
        else:
            print(f"ERROR: sort_GWASSS_by_rsID script finished with exit code: {ec}")

    elif (issues['rsID'] or issues['OA'] or issues['EA']) and issues['Chr']<total_entries and issues['BP']<total_entries:
        required_sorting = True
        sorted_by = 'ChrBP'
        for col in ('rsID', 'OA', 'EA'):
            if issues[col]:
                print(f"{issues[col]}/{total_entries} entries are missing {col}")
        print(f"Going to sort the GWAS SS file by Chr and BP")
        ec = call(["python3",
            sort_GWASSS_by_ChrBP,
            INPUT_GWAS_FILE_standard_lifted if required_liftover else INPUT_GWAS_FILE_standard,
            INPUT_GWAS_FILE_standard_sorted,
            ])
        if ec == 0:
            print(f"Sorted by Chr and BP")
        else:
            print(f"ERROR: sort_GWASSS_by_ChrBP script finished with exit code: {ec}")

    print(f"  Step {i_step} finished in {(time.time() - start_time)} seconds\n")

    if not any(issues.values()):
        print(f"The input summary statistics file has not been identified to have any issues!")
        print(f"all {total_entries} SNPs are good")
        exit(0)

    if required_sorting and ec != 0:
        exit(13)



    ##### 4 #####
    i_step += 1
    print(f'=== Step {i_step}: REHAB: loopping through the GWAS SS file and fixing entries ===')
    start_time = time.time()

    FILE_FOR_FIXING = INPUT_GWAS_FILE_standard_sorted if required_sorting else INPUT_GWAS_FILE_standard
    REHAB_OUTPUT_FILE = OUTPUT_FILE + '.SSrehabed.tsv'
    ec = call(["python3",
            loop_fix,
            FILE_FOR_FIXING,
            input_validation_report_dir,
            REHAB_OUTPUT_FILE,
            dbSNP_FILE,
            dbSNP2_FILE,
            'None',
            FREQ_DATABASE_SLUG if FREQ_DATABASE_SLUG else 'None',
            sorted_by if sorted_by else ''
            ])
    print(f"  Step {i_step} finished in {(time.time() - start_time)} seconds\n")

    if ec != 0:
        print(f"ERROR: loop_fix script finished with exit code: {ec}")
        exit(14)
    else:
        print(f"see fixed file at: \"{REHAB_OUTPUT_FILE}\"")



    ##### 5 #####
    i_step += 1
    print(f'=== Step {i_step}: Validate entries in the fixed GWAS SS file and save the report ===')
    start_time = time.time()

    REHABed_validation_report_dir = INPUT_GWAS_FILE + "_REHABed-report"
    ec = call(["python3",
            validate_GWASSS_entries,
            REHAB_OUTPUT_FILE,
            "standard",
            REHABed_validation_report_dir,
            ])
    print(f"  Step {i_step} finished in {(time.time() - start_time)} seconds\n")

    if ec != 0:
        print(f"ERROR: validate_GWASSS_entries script finished with exit code: {ec}")
        exit(15)



    ##### 6 #####
    i_step += 1
    print(f'=== Step {i_step}: Analyze the report after REHAB ===')
    start_time = time.time()

    issues_REHABed, total_entries = read_report_from_dir(REHABed_validation_report_dir)
    issues_solved = {c: issues[c]-issues_REHABed[c] for c in issues}
    for col in STANDARD_COLUMN_ORDER:
        if col not in ('N', 'INFO') and issues_solved[col]:
            if issues_solved[col] < 0:
                if col in ('Chr', 'BP') and get_build() != input_build:
                    print(f"lost {-issues_solved[col]} ({perc(-issues_solved[col], total_entries)}) \"{col}\" fields after liftover")
                else:
                    print(f"lost {-issues_solved[col]} ({perc(-issues_solved[col], total_entries)}) \"{col}\" fields")
            else:
                print(f"restored {issues_solved[col]} ({perc(issues_solved[col], total_entries)}) \"{col}\" fields")

    INPUT_GWAS_FILE_standard_sorted2 = remove_last_ext(INPUT_GWAS_FILE) + "_standard_sorted2.tsv"

    required_sorting2: bool = False
    if sorted_by != 'ChrBP' and (issues_REHABed['rsID'] or issues_REHABed['OA'] or issues_REHABed['EA']):
        required_sorting2 = True
        sorted_by = 'ChrBP'
        for col in ('rsID', 'OA', 'EA'):
            if issues_REHABed[col]:
                print(f"{issues_REHABed[col]}/{total_entries} entries are missing {col}")
        print(f"Going to sort the GWAS SS file by Chr and BP")
        ec = call(["python3",
            sort_GWASSS_by_ChrBP,
            REHAB_OUTPUT_FILE,
            INPUT_GWAS_FILE_standard_sorted2,
            ])
        if ec == 0:
            print(f"Sorted by Chr and BP")
        else:
            print(f"ERROR: sort_GWASSS_by_ChrBP script finished with exit code: {ec}")

    print(f"  Step {i_step} finished in {(time.time() - start_time)} seconds\n")

    if not any(issues_REHABed.values()):
        print(f"The REHABed summary statistics file has not been identified to have any issues!")
        print(f"all {total_entries} SNPs are good")
        exit(0)

    if required_sorting2 and ec != 0:
        exit(16)
    if not required_sorting2:
        # if the file doesn't require any other sorting,
        # then with the current sorting everything that could be restored is already restored
        exit(0)



    ##### 7 #####
    i_step += 1
    print(f'=== Step {i_step}: REHAB: loopping through the GWAS SS file again and fixing entries ===')
    start_time = time.time()

    FILE_FOR_FIXING = INPUT_GWAS_FILE_standard_sorted2 if required_sorting2 else REHAB_OUTPUT_FILE
    REHAB2_OUTPUT_FILE = OUTPUT_FILE + '.SSrehabed-twice.tsv'
    ec = call(["python3",
            loop_fix,
            FILE_FOR_FIXING,
            REHABed_validation_report_dir,
            REHAB2_OUTPUT_FILE,
            dbSNP_FILE,
            dbSNP2_FILE,
            'None', # setting to None suppresses liftover second time
            FREQ_DATABASE_SLUG if FREQ_DATABASE_SLUG else 'None',
            sorted_by if sorted_by else ''
            ])
    print(f"  Step {i_step} finished in {(time.time() - start_time)} seconds\n")

    if ec != 0:
        print(f"ERROR: loop_fix script finished with exit code: {ec}")
        exit(17)
    else:
        print(f"see fixed file at: \"{REHAB2_OUTPUT_FILE}\"")



    ##### 8 #####
    i_step += 1
    print(f'=== Step {i_step}: Validate entries in the twice REHABed GWAS SS file and save the report ===')
    start_time = time.time()

    REHABed_twice_validation_report_dir = INPUT_GWAS_FILE + "_REHABed-twice-report"
    ec = call(["python3",
            validate_GWASSS_entries,
            REHAB2_OUTPUT_FILE,
            "standard",
            REHABed_twice_validation_report_dir,
            ])
    print(f"  Step {i_step} finished in {(time.time() - start_time)} seconds\n")

    if ec != 0:
        print(f"ERROR: validate_GWASSS_entries script finished with exit code: {ec}")
        exit(18)



    ##### 9 #####
    i_step += 1
    print(f'=== Step {i_step}: Analyze the report after the second REHAB ===')
    start_time = time.time()

    issues_REHABed_twice, total_entries = read_report_from_dir(REHABed_twice_validation_report_dir)
    issues_solved = {c: issues[c]-issues_REHABed_twice[c] for c in issues}
    for col in STANDARD_COLUMN_ORDER:
        if col not in ('N', 'INFO') and issues_solved[col]:
            if issues_solved[col] < 0:
                print(f"lost {-issues_solved[col]} ({perc(-issues_solved[col], total_entries)}) \"{col}\" fields")
            else:
                print(f"restored {issues_solved[col]} ({perc(issues_solved[col], total_entries)}) \"{col}\" fields")


    print(f"  Step {i_step} finished in {(time.time() - start_time)} seconds\n")

    if not any(issues.values()):
        print(f"The twice REHABed summary statistics file has not been identified to have any issues!")
        print(f"all {total_entries} SNPs are good")
        exit(0)

    if required_sorting and ec != 0:
        exit(19)




def diagnose(INPUT_GWAS_FILE: str, REPORT_DIR: str):

    ### PROCESS INPUT ###
    INPUT_GWAS_FILE = str(INPUT_GWAS_FILE)
    REPORT_DIR = str(REPORT_DIR)

    JSON_CONFIG = INPUT_GWAS_FILE + '.json'
    if not os.path.isfile(JSON_CONFIG):
        JSON_CONFIG = 'standard'


    ### RUN ###
    print(f'=== Diagnosis ===')
    start_time = time.time()

    ec = call(["python3",
            validate_GWASSS_entries,
            INPUT_GWAS_FILE,
            JSON_CONFIG,
            REPORT_DIR if REPORT_DIR and REPORT_DIR != 'None' else '',
            ])
    print(f"  Diagnosis finished in {(time.time() - start_time)} seconds\n")

    if ec != 0:
        print(f"ERROR: validate_GWASSS_entries script finished with exit code: {ec}")
        exit(11)

    return None


def sort(INPUT_GWAS_FILE: str, OUTPUT_FILE: str, SORT_BY: str):

    ### PROCESS INPUT ###
    INPUT_GWAS_FILE = str(INPUT_GWAS_FILE)
    OUTPUT_FILE = str(OUTPUT_FILE)
    SORT_BY = str(SORT_BY)

    FILE_TO_SORT = INPUT_GWAS_FILE

    if os.path.isfile(INPUT_GWAS_FILE+".json"):
        ##### 1 #####
        print(f'=== Format the GWAS SS file ===')
        start_time = time.time()

        FILE_TO_SORT = OUTPUT_FILE + ".unsorted.tsv"
        ec = call(["python3",
                prepare_GWASSS_columns,
                INPUT_GWAS_FILE,
                FILE_TO_SORT,
                ])
        print(f"  Formatting finished in {(time.time() - start_time)} seconds\n")

        if ec != 0:
            print(f"ERROR: prepare_GWASSS_columns script finished with exit code: {ec}")
            exit(11)
    else:
        print("there's no corresponding .json file, so STANDARD_COLUMN_ORDER is assumed")


    ##### 2 #####
    print(f'=== Sorting the GWAS SS file by {SORT_BY} ===')
    start_time = time.time()

    if SORT_BY == 'rsID':
        ec = call(["python3",
            sort_GWASSS_by_rsID,
            FILE_TO_SORT,
            OUTPUT_FILE,
            ])
        if ec != 0:
            print(f"ERROR: sort_GWASSS_by_rsID script finished with exit code: {ec}")
    elif SORT_BY == 'ChrBP':
        ec = call(["python3",
            sort_GWASSS_by_ChrBP,
            FILE_TO_SORT,
            OUTPUT_FILE,
            ])
        if ec != 0:
            print(f"ERROR: sort_GWASSS_by_ChrBP script finished with exit code: {ec}")

    print(f"  Sorting finished in {(time.time() - start_time)} seconds\n")

    return None


def prepare_dbSNPs(SNPs_FILE: str, OUTPUT_FILE: str, gzsort: str, bcftools: str, buffer_size: str):

    ### PROCESS INPUT ###
    SNPs_FILE = str(SNPs_FILE)
    OUTPUT_FILE = str(OUTPUT_FILE)
    gzsort = str(gzsort)
    bcftools = str(bcftools)
    buffer_size = str(buffer_size)

    ### RUN ###
    ec = call(["python3",
            prepare_two_dbSNPs,
            SNPs_FILE,
            gzsort,
            bcftools,
            buffer_size,
            OUTPUT_FILE,
            ])

    if ec != 0:
        print(f"ERROR: validate_GWASSS_entries script finished with exit code: {ec}")
        exit(11)

    return None





# # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                 #
#                      MAIN                       #
#                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # #

##### Custom arguments types #####

def GWASSS_path_type(string: str):
    """
    A GWAS SS file must have a corresponding .json file:
    a JSON object where column indices are specified and build
    """
    if not os.path.isfile(string):
        raise ValueError("No such file:", string)
    if not os.path.isfile(string + '.json'):
        raise ValueError("No such file:", string + '.json')
    return string

def maybe_dir_type(string: str):
    if not os.path.isdir(string) and os.path.exists(string):
        raise ValueError("Unable to create directory at: ", string)
    return string

def file_path_type(string: str):
    if not os.path.isfile(string):
        raise ValueError("No such file:", string)
    return string

def maybe_file_path_type(string: str):
    if not string or string.lower() in ('none', 'na', 'null'):
        return None
    if not os.path.isfile(string):
        raise ValueError("No such file:", string)
    return string


def main():
    # version = "0.1.0"

    p = argparse.ArgumentParser(description='GWAS summary statistics QC tool')
    subparser = p.add_subparsers(dest='command')
    FIX_PARSER = subparser.add_parser('fix', help="diagnoses and tries to fix the file")
    PREPARE_DBSNPS_PARSER = subparser.add_parser('prepare_dbSNPs', help="prepares two DBs from the given dbSNP database. These two DBs are required for restoring rsID, chr, BP, alleles, and allele frequencies")
    DIAGNOSE_PARSER = subparser.add_parser('diagnose', help="only diagnosis. Produce report to a directory or just pop up plots")
    SORT_PARSER = subparser.add_parser('sort', help="sort GWAS SS file either by Chr:BP or rsID")


    # fix.add_argument('-v', '--version', action='version', version='%(prog)s {}'.format(version))
    FIX_PARSER.add_argument('--INPUT', dest='INPUT_GWAS_FILE', type=GWASSS_path_type, required=True,
        help='Path to GWAS summary stats in tab-separated format (.tsv, .tsv.gz, .tsv.zip), with a config file at the same path with .json suffix')
    FIX_PARSER.add_argument('--OUTPUT', dest='OUTPUT_FILE', type=pathlib.Path, required=True,
        help='Output path. This name will be used as a base for output file(s)')
    FIX_PARSER.add_argument('--dbsnp-1', dest='dbSNP1_FILE', type=maybe_file_path_type, required=False, default='None',
        help='Path to prepared dbSNP file #1 for the target build')
    FIX_PARSER.add_argument('--dbsnp-2', dest='dbSNP2_FILE', type=maybe_file_path_type, required=False, default='None',
        help='Path to prepared dbSNP file #2 for the target build')
    FIX_PARSER.add_argument('--chain-file', dest='CHAIN_FILE', type=maybe_file_path_type, required=False, default='None',
        help='Path to the chain file for liftover from given build to GrCh38')
    FIX_PARSER.add_argument('--freq-db', dest='FREQ_DATABASE_SLUG', type=str, required=False, default='dbGaP_PopFreq',
        help='Population slug from frequency database in dbSNP (e.g.: "GnomAD", "dbGaP_PopFreq", "TOMMO", "1000Genomes", etc.). Default: "dbGaP_PopFreq"')


    DIAGNOSE_PARSER.add_argument('--INPUT', dest='INPUT_GWAS_FILE', type=file_path_type, required=True,
        help='Path to GWAS summary stats in tab-separated format (.tsv, .tsv.gz), with a config file at the same path with .json suffix. If the config file is absent, internal "STANDARD_COLUMN_ORDER" is assumed.')
    DIAGNOSE_PARSER.add_argument('--REPORT-DIR', dest='REPORT_DIR', type=maybe_dir_type, required=False, default='None',
        help='A directory where the report should be saved. If not specified, doesn\'t save the full report and pops up matplotlib plots instead')


    SORT_PARSER.add_argument('--INPUT', dest='INPUT_GWAS_FILE', type=file_path_type, required=True,
        help='Path to GWAS summary stats in tab-separated format (.tsv, .tsv.gz, .tsv.zip), with a config file at the same path with .json suffix. If the config file is absent, internal "STANDARD_COLUMN_ORDER" is assumed.')
    SORT_PARSER.add_argument('--OUTPUT', dest='OUTPUT_FILE', type=pathlib.Path, required=True,
        help='Output path. This name will be used as a base for output file(s)')
    SORT_PARSER.add_argument('--by', dest='SORT_BY', choices=['rsID', 'ChrBP'], required=False, default='ChrBP',
        help='How to sort. Default: by Chr and BP')


    PREPARE_DBSNPS_PARSER.add_argument('--dbsnp', dest='DBSNP', type=file_path_type, required=True,
        help='Path to dbSNP file (.vcf, .vcf.gz, .bcf, bcf.gz)')
    PREPARE_DBSNPS_PARSER.add_argument('--OUTPUT', dest='OUTPUT', type=pathlib.Path, required=True,
        help='Base name for the two prepared DBs')
    PREPARE_DBSNPS_PARSER.add_argument('--gz-sort', dest='GZ_SORT', type=file_path_type, required=True,
        help='Path to gz-sort executable. Get from: https://github.com/keenerd/gz-sort/')
    PREPARE_DBSNPS_PARSER.add_argument('--bcftools', dest='BCFTOOLS', type=file_path_type, required=True,
        help='Path to bcftools executable. Get from: http://samtools.github.io/bcftools/ (recommended version - 1.11)')
    PREPARE_DBSNPS_PARSER.add_argument('--buffer', dest='BUFFER', type=str, required=False, default="1G",
        help='Buffer size for sorting (size of presort), supports k/M/G suffix. Default: 1G. Recommended: at least 200M, ideally 4G or more')


    args = p.parse_args()

    if args.command == 'fix':
        fix(args.INPUT_GWAS_FILE, args.OUTPUT_FILE,
            args.dbSNP1_FILE, args.dbSNP2_FILE, args.CHAIN_FILE, args.FREQ_DATABASE_SLUG)

    elif args.command == 'diagnose':
        diagnose(args.INPUT_GWAS_FILE, args.REPORT_DIR)

    elif args.command == 'sort':
        sort(args.INPUT_GWAS_FILE, args.OUTPUT_FILE,
             args.SORT_BY)

    elif args.command == 'prepare_dbSNPs':
        prepare_dbSNPs(args.DBSNP, args.OUTPUT,
                       args.GZ_SORT, args.BCFTOOLS, args.BUFFER)

    else:
        p.print_help(sys.stderr)
        exit(1)




if __name__ == "__main__":
    main()

