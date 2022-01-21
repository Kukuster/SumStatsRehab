# standard library
import sys
import os

# local
from lib.utils import run_bash
from lib.standard_column_order import STANDARD_COLUMN_ORDER



def sort_GWASSS_by_rsID(GWAS_FILE: str, OUTPUT_FILE: str):
    """
    Sorts formatted GWAS summary stats file by two columns: chromosome and base pair position.

    Parameters
    ----------
    GWAS_FILE : str
        GWAS summary statistics file in the internal \"standard\" tsv format

    OUTPUT_FILE : str
        output file name, GWAS summary statistics file in the internal \"standard\" tsv format sorted by rsID
    """

    if not os.path.isfile(GWAS_FILE):
        raise ValueError(f"passed GWAS SS file doesn't exist at path {GWAS_FILE}")

    rsID_col_i = STANDARD_COLUMN_ORDER.index('rsID') + 1
    run_bash(f'(head -n 1 "{GWAS_FILE}" && tail -n +2 "{GWAS_FILE}" | LC_ALL=C sort -t $\'\\t\' -k{rsID_col_i},{rsID_col_i}) > \"{OUTPUT_FILE}\"')



if __name__ == "__main__":
    GWAS_FILE = sys.argv[1]
    OUTPUT_FILE = sys.argv[2]

    sort_GWASSS_by_rsID(GWAS_FILE, OUTPUT_FILE)
