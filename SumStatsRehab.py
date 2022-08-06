# standard library
import sys
import os
from subprocess import call
import inspect
import time
from typing import Dict, List, Literal, Union
import json
import argparse
import pathlib

# local
from lib.prepare_two_dbSNPs import prepare_two_dbSNPs
from lib.prepare_GWASSS_columns import prepare_GWASSS_columns
from lib.validate_GWASSS_entries import validate_GWASSS_entries
from lib.sort_GWASSS_by_ChrBP import sort_GWASSS_by_ChrBP
from lib.sort_GWASSS_by_rsID import sort_GWASSS_by_rsID
from lib.loop_fix import ResolverName, resolvers_names, loop_fix, ActivatedResolvers
from lib.report_utils import read_report_from_dir
from lib.standard_column_order import STANDARD_COLUMN_ORDER
from lib.env import get_build, set_build
from lib.utils import mv, rm, rm_r, rm_rf




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


class brag:
    def see_formatted_file(self, OUTPUT_FILE: str):
        print(f" see formatted file at: \"{OUTPUT_FILE}\"")

    def see_fixed_file(self, OUTPUT_FILE: str):
        print(f" see fixed file at: \"{OUTPUT_FILE}\"")


    def all_issues_resolved(self, total_entries: int):
        print(f"All issues with all data points have been resolved!")
        print(f"all {total_entries} SNP entries are good")

    def all_which_can_be_solved_is_resolved(self):
        print(f"Those issues which were possible to resolve have been resolved")

    def the_input_file_has_no_issues(self, total_entries: int):
        print(f"The input summary statistics file doesn't seem to have any issues!")
        print(f"all {total_entries} SNPs are good")

    def the_input_file_has_nothing_to_resolve(self):
        print(f"The input file has nothing to resolve")

BRAG = brag()

class inform:
    def this_number_of_ChrBP_lost_after_liftover(self, Chr_lost: int, BP_lost: int, total_entries: int):
        if Chr_lost > 0:
            print(f"lost {Chr_lost} ({perc(Chr_lost,total_entries)}) \"Chr\" fields")
        if BP_lost > 0:
            print(f"lost {BP_lost}  ({perc(BP_lost, total_entries)}) \"BP\" fields")

    def no_need_for_liftover_since_ChrBP_will_be_restored(self):
        print(f"There's no need for liftover since all Chr and BP will be attempted to be restored in the target build")

INFORM = inform()
        
class warn:
    def impossible_to_liftover_since_all_Chr_or_BP_are_invalid(self):
        print(f"WARNING: Impossible to perform liftover. Liftover requires at least 1 entry with both valid Chr and BP. dbSNP2 wasn't passed, so Chr and BP will not be restored.")

    def beta_resolver_was_enabled(self):
        print(f"WARNING: if standard error field is provided without a sign, then restored beta will be unsigned.")

WARN = warn()


# # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                 #
#                    COMMANDS                     #
#                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # #



def fix(INPUT_GWAS_FILE: str,
        OUTPUT_FILE: str,
        dbSNP_FILE: str,
        dbSNP2_FILE: str,
        CHAIN_FILE: str,
        FREQ_DATABASE_SLUG: str,
        ACTIVATED_RESOLVERS: Dict[ResolverName, bool] = {
            "ChrBP": True,
            "rsID":  True,
            "OA":    True,
            "EA":    True,
            "EAF":   True,
            "beta":  False,
            "SE":    True,
            "pval":  True,
        },
        VERBOSE: bool = False,
    ):

    ### PROCESS INPUT ###
    INPUT_GWAS_FILE = str(INPUT_GWAS_FILE)
    JSON_CONFIG = INPUT_GWAS_FILE + '.json'
    OUTPUT_FILE = str(OUTPUT_FILE)
    dbSNP_FILE = str(dbSNP_FILE)
    dbSNP2_FILE = str(dbSNP2_FILE)
    CHAIN_FILE = str(CHAIN_FILE)
    FREQ_DATABASE_SLUG = str(FREQ_DATABASE_SLUG)
    ACTIVATED_RESOLVERS = ActivatedResolvers(ACTIVATED_RESOLVERS)


    ### declare shortcut functions ###
    def gonna_resolve(field: str, issues: Dict[str, int]) -> bool:
        if field in ('Chr', 'BP'):
            return bool(issues[field] and dbSNP2_FILE != 'None' and ACTIVATED_RESOLVERS['ChrBP']) # Chr and BP are always resolved together
        elif field == 'rsID':
            return bool(issues[field] and dbSNP_FILE  != 'None' and ACTIVATED_RESOLVERS['rsID'])
        elif field in ('OA','EA','EAF'):
            return bool(issues[field] and (dbSNP_FILE != 'None' or dbSNP2_FILE != 'None') and ACTIVATED_RESOLVERS[field])
        else:
            return bool(issues[field] and ACTIVATED_RESOLVERS[field])
    
    def any_issues_to_resolve(issues: Dict[str, int]):
        issues_to_resolve: Dict[str, int] = {}
        for field, issue_count in issues.items():
            if issue_count > 0 and gonna_resolve(field, issues):
                issues_to_resolve[field] = True
        return any(issues_to_resolve.values())

    intermediate_files: List[str] = []
    def present_output(result_file: str):
        if VERBOSE:
            return result_file
        else:
            for junk in intermediate_files:
                rm_rf(junk)
            mv(result_file, OUTPUT_FILE)
            return OUTPUT_FILE



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
    prepare_GWASSS_columns(
        INPUT_GWAS_FILE,
        INPUT_GWAS_FILE_standard,
    )
    print(f"  Step {i_step} finished in {(time.time() - start_time)} seconds\n")


    ##### 2 #####
    i_step += 1
    print(f'=== Step {i_step}: Validate entries in the formatted GWAS SS file and save the report ===')
    start_time = time.time()

    input_validation_report_dir = INPUT_GWAS_FILE + "_input-report"
    validate_GWASSS_entries(
        INPUT_GWAS_FILE_standard,
        "standard",
        input_validation_report_dir,   
    )
    intermediate_files.append(input_validation_report_dir)
    print(f"  Step {i_step} finished in {(time.time() - start_time)} seconds\n")


    ##### 3 #####
    i_step += 1
    print(f'=== Step {i_step}: Analyze the report and prepare for REHAB ===')
    start_time = time.time()

    issues, total_entries = read_report_from_dir(input_validation_report_dir)

    required_sorting: bool = False
    required_liftover: bool = False
    ChrBP_lost_because_of_liftover: int = 0
    sorted_by: Literal[None, 'rsID', 'ChrBP'] = None
    INPUT_GWAS_FILE_standard_sorted = remove_last_ext(INPUT_GWAS_FILE) + "_standard_sorted.tsv"
    INPUT_GWAS_FILE_standard_lifted = remove_last_ext(INPUT_GWAS_FILE) + "_standard_lifted.tsv"
    INPUT_GWAS_FILE_prepared: str = INPUT_GWAS_FILE_standard


    if get_build() != 'hg38' and CHAIN_FILE and CHAIN_FILE != "None":
        if issues['BP']<total_entries and issues['Chr']<total_entries:
            required_liftover = True
            ChrBP_lost_because_of_liftover = loop_fix(
                INPUT_GWAS_FILE_standard,
                input_validation_report_dir,
                INPUT_GWAS_FILE_standard_lifted,
                dbSNP_FILE,
                dbSNP2_FILE,
                CHAIN_FILE,
                FREQ_DATABASE_SLUG if FREQ_DATABASE_SLUG else 'None',
                sorted_by if sorted_by else None,
            )
            INPUT_GWAS_FILE_prepared = INPUT_GWAS_FILE_standard_lifted
            intermediate_files.append(INPUT_GWAS_FILE_standard)

            print("finished liftover to hg38 (saved report)")
            set_build('hg38')

            input_lifted_validation_report_dir = INPUT_GWAS_FILE + "_input-lifted-report"
            validate_GWASSS_entries(
                INPUT_GWAS_FILE_standard_lifted,
                "standard",
                input_lifted_validation_report_dir,
            )
            intermediate_files.append(input_lifted_validation_report_dir)
        elif dbSNP2_FILE != 'None':
            # if either Chr or BP is fully missing, there's no need for liftover.
            # Because Chr and BP will be attempted to be restored with dbSNPs in the target build
            set_build('hg38')
            INFORM.no_need_for_liftover_since_ChrBP_will_be_restored()
        else:
            WARN.impossible_to_liftover_since_all_Chr_or_BP_are_invalid()
    
    if (
            gonna_resolve('BP',issues) or 
            gonna_resolve('Chr',issues) or 
            gonna_resolve('EAF',issues)
        ) and issues['rsID']<total_entries:
        # here if some alleles are invalid, they will be attempted to be restored by rsID,
        # which is better then by ChrBP
        required_sorting = True
        sorted_by = 'rsID'
        for col in ('Chr', 'BP'):
            if issues[col]:
                print(f"{issues[col]}/{total_entries} entries are missing {col}")
        print(f"Going to sort the GWAS SS file by rsID")
        sort_GWASSS_by_rsID(
            INPUT_GWAS_FILE_standard_lifted if required_liftover else INPUT_GWAS_FILE_standard,
            INPUT_GWAS_FILE_standard_sorted,
        )
        intermediate_files.append(INPUT_GWAS_FILE_standard_lifted)
        intermediate_files.append(INPUT_GWAS_FILE_standard)
        
        INPUT_GWAS_FILE_prepared = INPUT_GWAS_FILE_standard_sorted
        print(f"Sorted by rsID")

    elif (
            gonna_resolve('rsID',issues) or 
            gonna_resolve('OA',issues) or 
            gonna_resolve('EA',issues) or 
            gonna_resolve('EAF',issues)
        ) and issues['Chr']<total_entries and issues['BP']<total_entries:
        # however if anyway going to restore any alleles and either:
        #  - all rsIDs are missing, or 
        #  - if going to restore rsID too
        # then sort by ChrBP and those will be restored from Chr and BP
        required_sorting = True
        sorted_by = 'ChrBP'
        for col in ('rsID', 'OA', 'EA'):
            if issues[col]:
                print(f"{issues[col]}/{total_entries} entries are missing {col}")
        print(f"Going to sort the GWAS SS file by Chr and BP")
        sort_GWASSS_by_ChrBP(
            INPUT_GWAS_FILE_standard_lifted if required_liftover else INPUT_GWAS_FILE_standard,
            INPUT_GWAS_FILE_standard_sorted,
        )
        intermediate_files.append(INPUT_GWAS_FILE_standard_lifted)
        intermediate_files.append(INPUT_GWAS_FILE_standard)
        INPUT_GWAS_FILE_prepared = INPUT_GWAS_FILE_standard_sorted
        print(f"Sorted by Chr and BP")

    print(f"  Step {i_step} finished in {(time.time() - start_time)} seconds\n")


    if not any(issues.values()) and not ChrBP_lost_because_of_liftover:
        BRAG.the_input_file_has_no_issues(total_entries)
        BRAG.see_formatted_file(
            present_output(INPUT_GWAS_FILE_prepared)
        )
        return
    elif not any(issues.values()) and ChrBP_lost_because_of_liftover:
        INFORM.this_number_of_ChrBP_lost_after_liftover(
            Chr_lost=ChrBP_lost_because_of_liftover,
             BP_lost=ChrBP_lost_because_of_liftover,
            total_entries=total_entries,
        )
        BRAG.the_input_file_has_nothing_to_resolve()
        BRAG.see_formatted_file(
            present_output(INPUT_GWAS_FILE_prepared)
        )
        return
    elif not any_issues_to_resolve(issues):
        BRAG.the_input_file_has_nothing_to_resolve()
        BRAG.see_formatted_file(
            present_output(INPUT_GWAS_FILE_prepared)
        )
        return



    ##### 4 #####
    i_step += 1
    print(f'=== Step {i_step}: REHAB: loopping through the GWAS SS file and fixing entries ===')
    start_time = time.time()

    FILE_FOR_FIXING = INPUT_GWAS_FILE_prepared
    REHAB_OUTPUT_FILE = OUTPUT_FILE + '.rehabed.tsv'
    loop_fix(
        FILE_FOR_FIXING,
        input_validation_report_dir,
        REHAB_OUTPUT_FILE,
        dbSNP_FILE,
        dbSNP2_FILE,
        'None',
        FREQ_DATABASE_SLUG if FREQ_DATABASE_SLUG else 'None',
        sorted_by if sorted_by else None,
        ACTIVATED_RESOLVERS,
    )
    intermediate_files.append(FILE_FOR_FIXING)
    print(f"  Step {i_step} finished in {(time.time() - start_time)} seconds\n")



    ##### 5 #####
    i_step += 1
    print(f'=== Step {i_step}: Validate entries in the fixed GWAS SS file and save the report ===')
    start_time = time.time()

    REHABed_validation_report_dir = INPUT_GWAS_FILE + "_REHABed-report"
    validate_GWASSS_entries(
        REHAB_OUTPUT_FILE,
        "standard",
        REHABed_validation_report_dir,   
    )
    intermediate_files.append(REHABed_validation_report_dir)
    print(f"  Step {i_step} finished in {(time.time() - start_time)} seconds\n")



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
            else: # issues_solved[col] > 0:
                print(f"restored {issues_solved[col]} ({perc(issues_solved[col], total_entries)}) \"{col}\" fields")

    INPUT_GWAS_FILE_standard_sorted2 = remove_last_ext(INPUT_GWAS_FILE) + "_standard_sorted2.tsv"

    required_sorting2: bool = False
    if sorted_by != 'ChrBP' and (
            gonna_resolve('rsID',issues_REHABed) or
            gonna_resolve('OA',issues_REHABed) or
            gonna_resolve('EA',issues_REHABed) or
            gonna_resolve('EAF',issues_REHABed)
        ) and \
        (not issues['Chr']==total_entries and not issues['BP']==total_entries):
        # if either rsID, OA, or EA are invalid, we can try sorting by ChrBP to restore them.
        # But if Chr or BP was totally missing at first, it won't help to try this

        required_sorting2 = True
        sorted_by = 'ChrBP'
        for col in ('rsID', 'OA', 'EA'):
            if issues_REHABed[col]:
                print(f"{issues_REHABed[col]}/{total_entries} entries are missing {col}")
        print(f"Going to sort the GWAS SS file by Chr and BP")
        sort_GWASSS_by_ChrBP(
            REHAB_OUTPUT_FILE,
            INPUT_GWAS_FILE_standard_sorted2,
        )
        intermediate_files.append(REHAB_OUTPUT_FILE)

        print(f"Sorted by Chr and BP")

    print(f"  Step {i_step} finished in {(time.time() - start_time)} seconds\n")

    if not any(issues_REHABed.values()):
        BRAG.all_issues_resolved(total_entries)
        BRAG.see_fixed_file(
            present_output(REHAB_OUTPUT_FILE)
        )
        return

    if not required_sorting2:
        # if the file doesn't require any other sorting,
        # then with the current sorting everything that could be restored is already restored
        BRAG.all_which_can_be_solved_is_resolved()
        BRAG.see_fixed_file(
            present_output(REHAB_OUTPUT_FILE)
        )
        return



    ##### 7 #####
    i_step += 1
    print(f'=== Step {i_step}: REHAB: loopping through the GWAS SS file again and fixing entries ===')
    start_time = time.time()

    FILE_FOR_FIXING = INPUT_GWAS_FILE_standard_sorted2 if required_sorting2 else REHAB_OUTPUT_FILE
    REHAB2_OUTPUT_FILE = OUTPUT_FILE + '.rehabed-twice.tsv'
    loop_fix(
        FILE_FOR_FIXING,
        REHABed_validation_report_dir,
        REHAB2_OUTPUT_FILE,
        dbSNP_FILE,
        dbSNP2_FILE,
        'None',  # setting to None suppresses liftover second time
        FREQ_DATABASE_SLUG if FREQ_DATABASE_SLUG else 'None',
        sorted_by if sorted_by else None,
        ACTIVATED_RESOLVERS,
    )
    intermediate_files.append(FILE_FOR_FIXING)
    print(f"  Step {i_step} finished in {(time.time() - start_time)} seconds\n")



    ##### 8 #####
    i_step += 1
    print(f'=== Step {i_step}: Validate entries in the twice REHABed GWAS SS file and save the report ===')
    start_time = time.time()

    REHABed_twice_validation_report_dir = INPUT_GWAS_FILE + "_REHABed-twice-report"
    validate_GWASSS_entries(
        REHAB2_OUTPUT_FILE,
        "standard",
        REHABed_twice_validation_report_dir,
    )
    intermediate_files.append(REHABed_twice_validation_report_dir)
    print(f"  Step {i_step} finished in {(time.time() - start_time)} seconds\n")



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

    if not any(issues_REHABed_twice.values()):
        BRAG.all_issues_resolved(total_entries)
        BRAG.see_fixed_file(
            present_output(REHAB2_OUTPUT_FILE)
        )
        return
    else:
        BRAG.all_which_can_be_solved_is_resolved()
        BRAG.see_fixed_file(
            present_output(REHAB2_OUTPUT_FILE)
        )
        return




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

    validate_GWASSS_entries(
        INPUT_GWAS_FILE,
        JSON_CONFIG,
        REPORT_DIR if REPORT_DIR and REPORT_DIR != 'None' else None,
    )
    print(f"  Diagnosis finished in {(time.time() - start_time)} seconds\n")

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
        prepare_GWASSS_columns(
            INPUT_GWAS_FILE,
            FILE_TO_SORT,
        )
        print(f"  Formatting finished in {(time.time() - start_time)} seconds\n")

    else:
        print("there's no corresponding .json file, so STANDARD_COLUMN_ORDER is assumed")


    ##### 2 #####
    print(f'=== Sorting the GWAS SS file by {SORT_BY} ===')
    start_time = time.time()

    if SORT_BY == 'rsID':
        sort_GWASSS_by_rsID(
            FILE_TO_SORT,
            OUTPUT_FILE,
        )
    elif SORT_BY == 'ChrBP':
        sort_GWASSS_by_ChrBP(
            FILE_TO_SORT,
            OUTPUT_FILE,
        )

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
    prepare_two_dbSNPs(
        SNPs_FILE,
        gzsort,
        bcftools,
        buffer_size,
        OUTPUT_FILE,
    )

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
    version = "1.2.0"

    p = argparse.ArgumentParser(description='GWAS summary statistics QC tool')
    p.prog = 'SumStatsRehab'
    subparser = p.add_subparsers(dest='command')
    FIX_PARSER = subparser.add_parser('fix', help="diagnoses and tries to fix the file")
    PREPARE_DBSNPS_PARSER = subparser.add_parser('prepare_dbSNPs', help="prepares two DBs from the given dbSNP database. These two DBs are required for restoring rsID, chr, BP, alleles, and allele frequencies")
    DIAGNOSE_PARSER = subparser.add_parser('diagnose', help="only diagnosis. Produce report to a directory or just pop up plots")
    SORT_PARSER = subparser.add_parser('sort', help="sort GWAS SS file either by Chr:BP or rsID")


    # fix.add_argument('-v', '--version', action='version', version='%(prog)s {}'.format(version))
    FIX_PARSER.add_argument('--INPUT', dest='INPUT_GWAS_FILE', type=GWASSS_path_type, required=True,
        help='Path to GWAS summary stats in tab-separated format (.tsv, .tsv.gz, .tsv.zip), with a config file at the same path with .json suffix')
    FIX_PARSER.add_argument('--OUTPUT', dest='OUTPUT_FILE', type=pathlib.Path, required=True,
        help='Output path for the final fixed file.\nIf --verbose key is set, then this name will be used as a base (prefix) for output file(s)')
    FIX_PARSER.add_argument('--dbsnp-1', dest='dbSNP1_FILE', type=maybe_file_path_type, required=False, default='None',
        help='Path to prepared dbSNP file #1 for the target build')
    FIX_PARSER.add_argument('--dbsnp-2', dest='dbSNP2_FILE', type=maybe_file_path_type, required=False, default='None',
        help='Path to prepared dbSNP file #2 for the target build')
    FIX_PARSER.add_argument('--chain-file', dest='CHAIN_FILE', type=maybe_file_path_type, required=False, default='None',
        help='Path to the chain file for liftover from the given build to GrCh38')
    FIX_PARSER.add_argument('--freq-db', dest='FREQ_DATABASE_SLUG', type=str, required=False, default='dbGaP_PopFreq',
        help='Population slug from frequency database in dbSNP (e.g.: "GnomAD", "dbGaP_PopFreq", "TOMMO", "1000Genomes", etc.). Default: "dbGaP_PopFreq"')
    FIX_PARSER.add_argument('--restore', choices=resolvers_names, nargs='+',
        help=f'Enable resotration of particular fields', required=False)
    FIX_PARSER.add_argument('--do-not-restore', choices=resolvers_names, nargs='+',
        help=f'Disable resotration of particular fields. By default, everything is enabled except "beta". This key takes priority over --restore', required=False)
    FIX_PARSER.add_argument('--verbose', dest='VERBOSE', action='store_true',
        help=f"If set, preserves all intermediate files and their diagnoses on all stages of fixing.\n" +
        "The --OUTPUT key path will be used as a base name for the resulting files, and will not be the exact path to a file\n" + 
        "This key doesn't affect logging."
        , required=False)


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

    if args.command:
        print(f"{p.prog} v{version} - {args.command} command")
    else:
        print(f"{p.prog} v{version}")

    if args.command == 'fix':
        chosen_resolvers: Dict[ResolverName, bool] = {
            "ChrBP": True,
            "rsID":  True,
            "OA":    True,
            "EA":    True,
            "EAF":   True,
            "beta":  False,
            "SE":    True,
            "pval":  True,
        }
        if args.restore:
            for field in args.restore:
                chosen_resolvers[field] = True
        if args.do_not_restore:
            for field in args.do_not_restore:
                chosen_resolvers[field] = False
        if chosen_resolvers['beta']:
            WARN.beta_resolver_was_enabled()

        fix(args.INPUT_GWAS_FILE, args.OUTPUT_FILE,
            args.dbSNP1_FILE, args.dbSNP2_FILE, args.CHAIN_FILE, args.FREQ_DATABASE_SLUG,
            chosen_resolvers, args.VERBOSE)

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

