# standard library
import sys
import os
from subprocess import call
import inspect
import time
from typing import Dict, Literal, Union
import json

# local
from lib.validate_utils import read_report_from_dir
from lib.standard_column_order import STANDARD_COLUMN_ORDER
from lib.env import get_build, set_build




# # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                 #
#                      INPUT                      #
#                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # #

if len(sys.argv) < 4:  # the very first 0th arg is the name of this script
    print("ERROR: you should specify args:")
    print("  #1 GWAS summary statistics file in tsv format (bare, zipped, or gzipped), that has a corresponding config file (suffixed \".json\") with column indices")
    print("  #2 output file name for a fixed file")
    print("  #3 dbSNP file")
    print("  #4 dbSNP file sorted by rsID")
    print("  #5 (optional) chain file")
    exit(1)

# INPUT_GWAS_FILE has to be in a tabular tab-sep format with a header on the first line
INPUT_GWAS_FILE = sys.argv[1]
JSON_CONFIG = sys.argv[1] + '.json'
OUTPUT_FILE = sys.argv[2]
dbSNP_FILE = sys.argv[3]
dbSNP2_FILE = sys.argv[4]

CHAIN_FILE: Union[str, None] = None
if len(sys.argv) > 5:
    CHAIN_FILE = sys.argv[5]


if not os.path.isfile(INPUT_GWAS_FILE):
    print(f"ERROR: provided GWAS SS file doesn't exist: {INPUT_GWAS_FILE}")
    exit(2)

if not os.path.isfile(JSON_CONFIG):
    print(f"ERROR: there's no corresponding json config file: {JSON_CONFIG}. Please create one based on config.example.json file")
    exit(2)

if not os.path.isfile(dbSNP_FILE):
    print(f"ERROR: there's no dbSNP file at the path: {dbSNP_FILE}")
    exit(2)

if not os.path.isfile(dbSNP2_FILE):
    print(f"ERROR: there's no preprocessed dbSNP file at the path: {dbSNP2_FILE}")
    exit(2)

if CHAIN_FILE and not os.path.isfile(CHAIN_FILE):
    print(f"ERROR: there's no chain file at the path: {CHAIN_FILE}")


# define paths to libs
the_dir = os.path.dirname(inspect.getfile(inspect.currentframe())) or os.getcwd()  # type: ignore
lib_dir = the_dir + "/lib"
prepare_GWASSS_columns = lib_dir+"/prepare_GWASSS_columns.py"
validate_GWASSS_entries = lib_dir+"/validate_GWASSS_entries.py"
sort_GWASSS_by_ChrBP = lib_dir+"/sort_GWASSS_by_ChrBP.py"
sort_GWASSS_by_rsID = lib_dir+"/sort_GWASSS_by_rsID.py"
loop_fix = lib_dir+"/loop_fix.py"


def remove_last_ext(filename: str):
    return filename.rsplit(".", 1)[0] # passed 1 means do max 1 split; _rightmost_ splits first


# # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                 #
#                      MAIN                       #
#                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # #


### Set the environment varaible for the build ###
config: Dict[str, Union[int, str]] = json.load(open(JSON_CONFIG,))
try:
    input_build = config['build']
except:
    raise ValueError('config file that corresponds to the GWAS SS file has to have "build" key set')
set_build(input_build)
print(f'build of the GWAS SS file: {get_build()}')


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

input_validation_report_dir = INPUT_GWAS_FILE.split('/')[-1] + "_input-report"
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
sorted_by: Literal[None, 'rsID', 'ChrBP'] = None
INPUT_GWAS_FILE_standard_sorted = remove_last_ext(INPUT_GWAS_FILE) + "_standard_sorted.tsv"

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
           INPUT_GWAS_FILE_standard,
           INPUT_GWAS_FILE_standard_sorted,
           ])
    if ec == 0:
        print(f"Sorted by rsID")
    else:
        print(f"ERROR: sort_GWASSS_by_rsID script finished with exit code: {ec}")

elif issues['rsID'] or issues['OA'] or issues['EA']:
    required_sorting = True
    sorted_by = 'ChrBP'
    for col in ('rsID', 'OA', 'EA'):
        if issues[col]:
            print(f"{issues[col]}/{total_entries} entries are missing {col}")
    print(f"Going to sort the GWAS SS file by Chr and BP")
    ec = call(["python3",
           sort_GWASSS_by_ChrBP,
           INPUT_GWAS_FILE_standard,
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
REHAB_OUTPUT_FILE = remove_last_ext(OUTPUT_FILE) + '_SSrehabed.tsv'
ec = call(["python3",
           loop_fix,
           FILE_FOR_FIXING,
           input_validation_report_dir,
           REHAB_OUTPUT_FILE,
           dbSNP_FILE,
           dbSNP2_FILE,
           CHAIN_FILE if CHAIN_FILE else 'None',
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

REHABed_validation_report_dir = INPUT_GWAS_FILE.split('/')[-1] + "_REHABed-report"
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
        if issues_solved[col] < 0 and col in ('Chr', 'BP') and get_build() != input_build:
            print(f"lost {issues_solved[col]} entries for \"{col}\" column after liftover")
        else:
            print(f"restored {issues_solved[col]} entries for \"{col}\" column")

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
REHAB2_OUTPUT_FILE = remove_last_ext(OUTPUT_FILE) + '_SSrehabed-twice.tsv'
ec = call(["python3",
           loop_fix,
           FILE_FOR_FIXING,
           REHABed_validation_report_dir,
           REHAB2_OUTPUT_FILE,
           dbSNP_FILE,
           dbSNP2_FILE,
           CHAIN_FILE if CHAIN_FILE else 'None',
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

REHABed_twice_validation_report_dir = INPUT_GWAS_FILE.split('/')[-1] + "_REHABed-twice-report"
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
        if issues_solved[col] < 0 and col in ('Chr', 'BP') and get_build() != input_build:
            print(f"lost {issues_solved[col]} entries for \"{col}\" column after liftover")
        else:
            print(f"restored {issues_solved[col]} entries for \"{col}\" column")


print(f"  Step {i_step} finished in {(time.time() - start_time)} seconds\n")

if not any(issues.values()):
    print(f"The twice REHABed summary statistics file has not been identified to have any issues!")
    print(f"all {total_entries} SNPs are good")
    exit(0)

if required_sorting and ec != 0:
    exit(19)
