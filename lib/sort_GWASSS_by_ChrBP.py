# standard library
import sys
from typing import Dict
import json
import os

# local
from utils import run_bash


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

#
# STEP #1
#    Prepend a column that represents the alphabetical order of the chromosomes using awk
#

GWAS_FILE_w_col = GWAS_FILE + "_chr-order-key.tsv"

run_bash(f"""paste -d$'\t' <( echo "chr_order_key" && tail -n +2 "{GWAS_FILE}" | awk 'BEGIN {{ 
    chrs["1"]  = "Aa";
    chrs["01"] = "Aa";
    chrs["2"]  = "Ab";
    chrs["02"] = "Ab";
    chrs["3"]  = "Ac";
    chrs["03"] = "Ac";
    chrs["4"]  = "Ad";
    chrs["04"] = "Ad";
    chrs["5"]  = "Ae";
    chrs["05"] = "Ae";
    chrs["6"]  = "Af";
    chrs["06"] = "Af";
    chrs["7"]  = "Ag";
    chrs["07"] = "Ag";
    chrs["8"]  = "Ah";
    chrs["08"] = "Ah";
    chrs["9"]  = "Ai";
    chrs["09"] = "Ai";
    chrs["10"] = "Aj";
    chrs["11"] = "Ak";
    chrs["12"] = "Al";
    chrs["13"] = "Am";
    chrs["14"] = "An";
    chrs["15"] = "Ao";
    chrs["16"] = "Ap";
    chrs["17"] = "Aq";
    chrs["18"] = "Ar";
    chrs["19"] = "As";
    chrs["20"] = "At";
    chrs["21"] = "Au";
    chrs["22"] = "Av";
    chrs["23"] = "Aw";
    chrs["24"] = "Ax";
    chrs["25"] = "Ay";
    chrs["26"] = "Az";
    chrs["X"]  = "Ba";
    chrs["x"]  = "Ba";
    chrs["Y"]  = "Bb";
    chrs["y"]  = "Bb";
    chrs["M"]  = "Bc";
    chrs["m"]  = "Bc";
 }}
 # awk has associative arrays
 # *in* operator checks if the string is in the array keys, which is probably O(1)
 # if the given chr entry is not recognized, the chr entry is used with "ZZ" prefix
 {{ print !($1 in chrs) ? "ZZ" $1 : chrs[$1] ; }}
' ) "{GWAS_FILE}" > "{GWAS_FILE_w_col}"
""")


#
# STEP #2
#    Sort the file by Chr and BP.
#
#    Specifically, sorts:
#        1. alphabetically by the first column (chromosome order key)
#        2. and numerically by the third column (BP)
#    and preserves the header
#

sorted_GWAS_FILE_w_col = GWAS_FILE + "_chr-order-key_sorted.tsv"
run_bash(f'(head -n 1 "{GWAS_FILE_w_col}" && tail -n +2 "{GWAS_FILE_w_col}" | sort -k1,1 -k3,3n) > "{sorted_GWAS_FILE_w_col}"')


#
# STEP #3
#    FINALLY, cuts the temporary column after the sorting, leaving properly sorted file
#

run_bash(f"cut -d$'\t' -f2-20 \"{sorted_GWAS_FILE_w_col}\" > \"{OUTPUT_FILE}\"")


exit(0)
