# standard library
import sys
import os
from subprocess import call
import inspect
import time

# local
from lib.validate_utils import read_report_from_dir




# # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                 #
#                      INPUT                      #
#                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # #

if len(sys.argv) < 3:  # the very first 0th arg is the name of this script
    print("ERROR: you should specify args:")
    print("  #1 GWAS summary statistics file in tsv format (bare, zipped, or gzipped), that has a corresponding config file (suffixed \".json\") with column indices")
    print("  #2 output file name for a fixed file")
    print("  #3 dbSNP file")
    exit(1)

# INPUT_GWAS_FILE has to be in a tabular tab-sep format with a header on the first line
INPUT_GWAS_FILE = sys.argv[1]
JSON_CONFIG = sys.argv[1] + '.json'
OUTPUT_FILE = sys.argv[2]
dbSNP_FILE = sys.argv[3]

if not os.path.isfile(INPUT_GWAS_FILE):
    print(f"ERROR: provided GWAS SS file doesn't exist: {INPUT_GWAS_FILE}")
    exit(2)

if not os.path.isfile(JSON_CONFIG):
    print(f"ERROR: there's no corresponding json config file: {JSON_CONFIG}. Please create one based on config.example.json file")
    exit(2)

if not os.path.isfile(dbSNP_FILE):
    print(f"ERROR: there's no dbSNP file at the path: {dbSNP_FILE}")
    exit(2)


# define paths to libs
the_dir = os.path.dirname(inspect.getfile(inspect.currentframe())) or os.getcwd()  # type: ignore
lib_dir = the_dir + "/lib"
prepare_GWASSS_columns = lib_dir+"/prepare_GWASSS_columns.py"
validate_GWASSS_entries = lib_dir+"/validate_GWASSS_entries.py"
sort_GWASSS_by_ChrBP = lib_dir+"/sort_GWASSS_by_ChrBP.py"
loop_fix = lib_dir+"/loop_fix.py"




# # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                 #
#                      MAIN                       #
#                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # #

numsteps: int = 5

##########
print(f'=== Step 1/{numsteps}: Format the GWAS SS file ===')
start_time = time.time()

INPUT_GWAS_FILE_standard = INPUT_GWAS_FILE + "_standard.tsv"
ec = call(["python3",
           prepare_GWASSS_columns,
           INPUT_GWAS_FILE,
           INPUT_GWAS_FILE_standard,
           ])
print(f"  Step 1/{numsteps} finished in {(time.time() - start_time)} seconds\n")

if ec != 0:
    print(f"ERROR: prepare_GWASSS_columns script finished with exit code: {ec}")
    exit(11)



##########
print(f'=== Step 2/{numsteps}: Validate entries in the formatted GWAS SS file and save the report ===')
start_time = time.time()

input_validation_report_dir = INPUT_GWAS_FILE.split('/')[-1] + "_input-report"
ec = call(["python3",
           validate_GWASSS_entries,
           INPUT_GWAS_FILE_standard,
           "standard",
           input_validation_report_dir,
           ])
print(f"  Step 2/{numsteps} finished in {(time.time() - start_time)} seconds\n")

if ec != 0:
    print(f"ERROR: validate_GWASSS_entries script finished with exit code: {ec}")
    exit(12)



##########
print(f'=== Step 3/{numsteps}: Analyze the report and prepare for REHAB ===')
start_time = time.time()

issues = read_report_from_dir(input_validation_report_dir)

required_sorting: bool = False
INPUT_GWAS_FILE_standard_sorted = INPUT_GWAS_FILE + "_standard_sorted.tsv"

if issues['rsID']:
    required_sorting = True
    print(f"{issues['rsID']}/{issues['total_entries']} entries are missing rsID")
    print(f"Going to sort the GWAS SS file by Chr and BP in order to restore rsIDs")
    ec = call(["python3",
           sort_GWASSS_by_ChrBP,
           INPUT_GWAS_FILE_standard,
           INPUT_GWAS_FILE_standard_sorted,
           ])
    if ec != 0:
        print(f"ERROR: sort_GWASSS_by_ChrBP script finished with exit code: {ec}")
    else:
        print(f"sorting finished")

elif issues['BP']:
    # need to sort by rsID
    pass

print(f"  Step 3/{numsteps} finished in {(time.time() - start_time)} seconds\n")

if not any(issues.values()):
    print(f"The input summary statistics file has not been identified to have any issues!")
    print(f"all {issues['total_entries']} SNPs are good")
    exit(0)

if required_sorting and ec != 0:
    exit(13)



##########
print(f'=== Step 4/{numsteps}: REHAB: loopping through the GWAS SS file and fixing entries ===')
start_time = time.time()

FILE_FOR_FIXING = INPUT_GWAS_FILE_standard_sorted if required_sorting else INPUT_GWAS_FILE_standard
ec = call(["python3",
           loop_fix,
           FILE_FOR_FIXING,
           input_validation_report_dir,
           OUTPUT_FILE,
           dbSNP_FILE,
           ])
print(f"  Step 4/{numsteps} finished in {(time.time() - start_time)} seconds\n")

if ec != 0:
    print(f"ERROR: loop_fix script finished with exit code: {ec}")
    exit(14)
else:
    print(f"see fixed file at: \"{OUTPUT_FILE}\"")


##########
print(f'=== Step 5/{numsteps}: Validate entries in the fixed GWAS SS file and save the report ===')
start_time = time.time()

fixed_validation_report_dir = INPUT_GWAS_FILE.split('/')[-1] + "_FIXED-report"
ec = call(["python3",
           validate_GWASSS_entries,
           OUTPUT_FILE,
           "standard",
           fixed_validation_report_dir,
           ])
print(f"  Step 5/{numsteps} finished in {(time.time() - start_time)} seconds\n")

if ec != 0:
    print(f"ERROR: validate_GWASSS_entries script finished with exit code: {ec}")
    exit(15)


##########
