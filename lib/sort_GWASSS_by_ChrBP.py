# standard library
import sys
import os

# local
from lib.utils import run_bash
from lib.standard_column_order import STANDARD_COLUMN_ORDER



def sort_GWASSS_by_ChrBP(GWAS_FILE: str, OUTPUT_FILE: str):
    """
    Sorts formatted GWAS summary stats file by two columns: chromosome and base pair position.

    Parameters
    ----------
    GWAS_FILE : str
        GWAS summary statistics file in the internal \"standard\" tsv format

    OUTPUT_FILE : str
        output file name, GWAS summary statistics file in the internal \"standard\" tsv format sorted by Chr and BP
    """

    if not os.path.isfile(GWAS_FILE):
        raise ValueError(f"passed GWAS SS file doesn't exist at path {GWAS_FILE}")


    #
    # STEP #1
    #    Prepend a column that represents the alphabetical order of the chromosomes using awk
    #

    GWAS_FILE_w_col = GWAS_FILE + "_chr-order-key.tsv"

    run_bash(f"""paste -d$'\t' <( echo "chr_order_key" && tail -n +2 "{GWAS_FILE}" | awk 'BEGIN {{ 
        chrs["1"]  = "Aa";
        chrs["01"] = "Aa";#
        chrs["2"]  = "Ab";
        chrs["02"] = "Ab";#
        chrs["3"]  = "Ac";
        chrs["03"] = "Ac";#
        chrs["4"]  = "Ad";
        chrs["04"] = "Ad";#
        chrs["5"]  = "Ae";
        chrs["05"] = "Ae";#
        chrs["6"]  = "Af";
        chrs["06"] = "Af";#
        chrs["7"]  = "Ag";
        chrs["07"] = "Ag";#
        chrs["8"]  = "Ah";
        chrs["08"] = "Ah";#
        chrs["9"]  = "Ai";
        chrs["09"] = "Ai";#
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
    # each chr number is mapped into the corresponding key, representing the relative order
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

    # once to shift to 1-indexed, once to account for prepended chr_order column
    BP_col_i = STANDARD_COLUMN_ORDER.index('BP') + 1 + 1
    run_bash(f'(head -n 1 "{GWAS_FILE_w_col}" && tail -n +2 "{GWAS_FILE_w_col}" | LC_ALL=C sort -t $\'\\t\' -k1,1 -k{BP_col_i},{BP_col_i}n) > "{sorted_GWAS_FILE_w_col}"')


    #
    # STEP #3
    #    FINALLY, cuts the temporary column after the sorting, leaving properly sorted file.
    #     and removes the intermediate files
    #
    run_bash(f"cut -d$'\t' -f2- \"{sorted_GWAS_FILE_w_col}\" > \"{OUTPUT_FILE}\"")
    run_bash(f'rm "{GWAS_FILE_w_col}" "{sorted_GWAS_FILE_w_col}"')




if __name__ == "__main__":
    GWAS_FILE = sys.argv[1]
    OUTPUT_FILE = sys.argv[2]

    sort_GWASSS_by_ChrBP(GWAS_FILE, OUTPUT_FILE)


