# standard library
import sys
import os

# local
from utils import run_bash
from standard_column_order import STANDARD_COLUMN_ORDER


# # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                 #
#                      INPUT                      #
#                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # #

if len(sys.argv) < 3:  # the very first 0th arg is the name of this script
    print("ERROR: you should specify args:")
    print("  #1 GWAS summary statistics file in the internal \"standard\" tsv format")
    print("  #2 output file name, prepared GWAS summary statistics file with all the columns in the standard order")
    exit(1)

# GWAS_FILE has to be in the internal "standard" tsv format
GWAS_FILE = sys.argv[1]
OUTPUT_FILE = sys.argv[2]


if not os.path.isfile(GWAS_FILE):
    print(f"ERROR: provided GWAS SS file doesn't exist: {GWAS_FILE}")
    exit(1)



# # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                 #
#                      MAIN                       #
#                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # #

rsID_col_i = STANDARD_COLUMN_ORDER.index('rsID') + 1
run_bash(f'(head -n 1 "{GWAS_FILE}" && tail -n +2 "{GWAS_FILE}" | LC_ALL=C sort -t $\'\\t\' -k{rsID_col_i},{rsID_col_i}) > \"{OUTPUT_FILE}\"')
