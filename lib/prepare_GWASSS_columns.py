# standard library
import sys
from typing import Dict, Union
import json
import os

# local
from lib.utils import run_bash
from lib.file import resolve_bare_text_file
from lib.standard_column_order import STANDARD_COLUMN_ORDER



def prepare_GWASSS_columns(INPUT_GWAS_FILE: str, OUTPUT_FILE: str):
    """
    Preprocesses the input GWAS summary statistics file with the .json config file with into the internal standardized format.

    Parameters
    ----------
    INPUT_GWAS_FILE : str
        GWAS summary statistics file in tsv format (bare, zipped, or gzipped), that has a corresponding config file (suffixed ".json") with column indices and build
    
    OUTPUT_FILE : str
        output file name, prepared GWAS summary statistics file with all the columns in the standard order

    """

    JSON_CONFIG = INPUT_GWAS_FILE + '.json'
    config: Dict[str, Union[int,str]] = json.load(open(JSON_CONFIG,))

    if not os.path.isfile(INPUT_GWAS_FILE):
        raise ValueError(f"passed GWAS SS file doesn't exist at path {INPUT_GWAS_FILE}")

    if not os.path.isfile(JSON_CONFIG):
        raise ValueError(f"passed GWAS SS file doesn't have the corresponding json config file at path: {JSON_CONFIG}. Please create one based on the config.example.json file")




    #
    # STEP #1
    #    Parse config file
    #

    cols_i: Dict[str, int] = {}
    for key, value in config.items():
        if isinstance(key, str) and isinstance(value, int):
            cols_i[key] = value


    #
    # STEP #2
    #    Unpack if the input file is an archive
    #
    BARE_GWAS_FILE = resolve_bare_text_file(INPUT_GWAS_FILE, f"{INPUT_GWAS_FILE}.tsv")


    #
    # STEP #3
    #    Reorder the columns using paste(1),
    #    while cutting with cut(1) on the fly using bash process substitution,
    #    and FINALLY save to the output filename specified by user
    #

    BASH_CMD = ["paste", "-d$'\\t'"]

    for i in range(len(STANDARD_COLUMN_ORDER)):
        col_name = STANDARD_COLUMN_ORDER[i]

        if col_name in cols_i.keys():
            # user specifies column indices as starting with 0,
            # whereas Unix cut(1) and awk count columns starting with 1
            c_i = current_col_index = cols_i[col_name] + 1

            # for any relevant column that's present, cut it.
            # if this is a chromosome column, make sure there's no "chr" prefix
            if col_name == 'Chr':
                chrom_col_fd = f"<( "
                chrom_col_fd += f"head -n 1 \"{BARE_GWAS_FILE}\" | cut -d$'\\t' -f{c_i} ; "
                chrom_col_fd += f"tail -n +2 \"{BARE_GWAS_FILE}\" | awk -F $'\\t' '{{if (tolower(${c_i}) ~ /^chr/) {{print substr(${c_i},4)}} else {{print ${c_i}}} }}' ; "
                chrom_col_fd += ")"
                BASH_CMD.append(chrom_col_fd)
            else:
                BASH_CMD.append(f"<(cut -d$'\\t' -f{c_i} \"{BARE_GWAS_FILE}\")")
        else:
            # if user didn't specify index for the column, a template column is added (header only)
            # in this case, paste(1) will leave such columns empty, i.e. values are the empty string
            BASH_CMD.append(f"<(echo {col_name}_rehab)")

    BASH_CMD.append(f">\"{OUTPUT_FILE}\"")


    bash_code = ' '.join(BASH_CMD)
    run_bash(bash_code)




if __name__ == "__main__":
    INPUT_GWAS_FILE = sys.argv[1]
    OUTPUT_FILE = sys.argv[2]

    prepare_GWASSS_columns(INPUT_GWAS_FILE, OUTPUT_FILE)
