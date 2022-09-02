# standard library
import sys
from typing import Dict, List, Union
from collections import OrderedDict
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

    # usually user specifies an integer for each column: a 0-indexed column index from the table
    cols_i: Dict[str, int] = {}

    # for EAF, a user may specify not a column index (int), but an object,
    #       where keys are column indices, and values are weights for those 
    input_cols_i_with_avg: Dict[str, Dict[Union[int, str], Union[int, float]]] = {}
    parsed_cols_i_with_avg: Dict[str, Dict[int, float]] = {}

    # another parameter that user may have specified is an array of additional columns to pick
    #   these columns will not participate in processing, but remain in the output as they are.
    # This may especially be useful if the SS file undegoes sorting at some stage.
    input_readonly_cols_i: Dict[str, List[Union[int, str]]] = {}
    parsed_readonly_cols_i: Dict[str, List[int]] = {}

    for key, value in config.items():
        if isinstance(key, str):
            if isinstance(value, int):
                cols_i[key] = value
            elif key == 'EAF' and isinstance(value, dict):
                input_cols_i_with_avg[key] = value
            elif key == 'other' and isinstance(value, list):
                input_readonly_cols_i[key] = value

    # validate json object format and parse numbers at the same time
    for key, value in input_cols_i_with_avg.items():
        try:
            for col, weight in value.items():
                if key not in parsed_cols_i_with_avg.keys():
                    parsed_cols_i_with_avg[key] = {}
                parsed_cols_i_with_avg[key][int(col)] = float(weight)
        except ValueError as e:
            print(e)
            raise ValueError(f"Invalid json configuration for column: {key}."+
            "Weighted average configuration has to be an object mapping column indices (int) to weights (numbers)")

    # validate json array format and parse at the same time
    for key, value in input_readonly_cols_i.items():
        try:
            for col in value:
                if key not in parsed_readonly_cols_i.keys():
                    parsed_readonly_cols_i[key] = []
                parsed_readonly_cols_i[key].append(int(float(col)))
        except ValueError as e:
            print(e)
            raise ValueError(f"Invalid json configuration for column: {key}."+
            "\"other\" configuration has to be an array of column indices (int)")


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

    awk_function_o_class = f'''
        # see answer https://unix.stackexchange.com/a/363471/387925
        # on question https://unix.stackexchange.com/questions/281271/can-i-determine-type-of-an-awk-variable
        function o_class(obj,   q, x, z){{
            q = CONVFMT
            CONVFMT = "% g"
                split(" " obj "\1" obj, x, "\1")
                x[1] = obj == x[1]
                x[2] = obj == x[2]
                x[3] = obj == 0
                x[4] = obj "" == +obj
            CONVFMT = q
            z["0001"] = z["1101"] = z["1111"] = "number"
            z["0100"] = z["0101"] = z["0111"] = "string"
            z["1100"] = z["1110"] = "strnum"
            z["0110"] = "undefined"
            return z[x[1] x[2] x[3] x[4]]
        }}
    '''
    # NOTE: '\1' in python will be interpreted into the same character as '\1' in awk, so no need to additioanlly escape it in python

    for i in range(len(STANDARD_COLUMN_ORDER)):
        col_name = STANDARD_COLUMN_ORDER[i]

        if col_name in cols_i.keys():
            # user specifies column indices as starting with 0,
            # whereas Unix cut(1) and awk count columns starting with 1
            c_i = cols_i[col_name] + 1

            # for any relevant column that's present, cut it.
            # if this is a chromosome column, make sure there's no "chr" prefix
            # if this is a BP column, make sure all numerical values are in integer format (not in float, or sci notation, etc.)
            # if this is an allele column, make all letters uppercase
            if col_name == 'Chr':
                chrom_col_fd = f"<( "
                chrom_col_fd += f"head -n 1 \"{BARE_GWAS_FILE}\" | cut -d$'\\t' -f{c_i} | sed -e 's/\r$//' ; "
                chrom_col_fd += f"tail -n +2 \"{BARE_GWAS_FILE}\" | awk -F $'\\t' '{{if (tolower(${c_i}) ~ /^chr/) {{print substr(${c_i},4)}} else {{print ${c_i}}} }}' | sed -e 's/\r$//' ; "
                chrom_col_fd += f" )"
                BASH_CMD.append(chrom_col_fd)
            elif col_name == 'BP':
                bp_col_fd = f"<( "
                bp_col_fd += f"head -n 1 \"{BARE_GWAS_FILE}\" | cut -d$'\\t' -f{c_i} ; "
                bp_col_fd += f"tail -n +2 \"{BARE_GWAS_FILE}\" | "
                bp_col_fd += f"""awk -F $'\\t' '
                {awk_function_o_class}
                {{
                    datum_type = o_class(${c_i})
                    if (datum_type == "number" || datum_type == "strnum"){{
                        printf("%.0f\\n", ${c_i})
                    }} else {{
                        print ${c_i}
                    }}
                }}' | """
                bp_col_fd += f" sed -e 's/\r$//' "
                bp_col_fd += f" )"
                BASH_CMD.append(bp_col_fd)
            elif col_name in ['OA', 'EA']:
                chrom_col_fd = f"<( "
                chrom_col_fd += f"head -n 1 \"{BARE_GWAS_FILE}\" | cut -d$'\\t' -f{c_i} | sed -e 's/\r$//' ; "
                chrom_col_fd += f"tail -n +2 \"{BARE_GWAS_FILE}\" | awk -F $'\\t' '{{print toupper(${c_i})}}' | sed -e 's/\r$//' ; "
                chrom_col_fd += f" )"
                BASH_CMD.append(chrom_col_fd)
            else:
                BASH_CMD.append(f"<(cut -d$'\\t' -f{c_i} \"{BARE_GWAS_FILE}\" | sed -e 's/\r$//')")

        elif col_name in parsed_cols_i_with_avg.keys():
            cols_obj = parsed_cols_i_with_avg[col_name]
            # user may have specified multiple columns with corresponding weights using syntax:
            # "EAF": {
            #    4: 35653,
            #    5: 25624
            # }
            # which denotes that column EAF should be a weighted average of columns 4 and 5,
            # with corresponding the weights

            shifted_cols_obj: OrderedDict[int, float] = OrderedDict({})
            for col_i, weight in cols_obj.items():
                shifted_cols_obj[col_i+1] = weight

            col_indices_for_cut = ",".join(map(str, shifted_cols_obj.keys()))
            col_weights_for_awk = " ".join(map(str, shifted_cols_obj.values()))

            avg_col_df = f"<( "
            avg_col_df += f"""echo {col_name}_rehab ;
            tail -n +2 \"{BARE_GWAS_FILE}\" | cut -d$'\t' -f{col_indices_for_cut} | sed -e 's/\r//' | \\
                awk -F$'\t' '
                    {awk_function_o_class}
                    BEGIN{{
                        split("{col_weights_for_awk}", weights, " ")

                        divisor = 0
                        for(i=1; i<=length(weights); i++)
                            divisor += weights[i];
                    }}
                    {{
                        rowsum=0
                        for(i=1;i<=NF;i++){{
                            datum_type = o_class($i)
                            if (datum_type == "number" || datum_type == "strnum"){{
                                rowsum += $i * weights[i]
                            }} else {{
                                print "."
                                next
                            }}
                        }}
                        print rowsum/divisor
                    }}
                '
            """
            avg_col_df += f" )"
            BASH_CMD.append(avg_col_df)

        else:
            # if user didn't specify index for the column, a template column is added (header only)
            # in this case, paste(1) will leave such columns empty, i.e. values are the empty string
            BASH_CMD.append(f"<(echo {col_name}_rehab)")


    if "other" in parsed_readonly_cols_i.keys():
        # user specifies column indices as starting with 0,
        # whereas Unix cut(1) and awk count columns starting with 1
        for col_i in parsed_readonly_cols_i["other"]:
            c_i = col_i + 1
            BASH_CMD.append(f"<(cut -d$'\\t' -f{c_i} \"{BARE_GWAS_FILE}\")")


    BASH_CMD.append(f">\"{OUTPUT_FILE}\"")


    bash_code = ' '.join(BASH_CMD)
    run_bash(bash_code)




if __name__ == "__main__":
    INPUT_GWAS_FILE = sys.argv[1]
    OUTPUT_FILE = sys.argv[2]

    prepare_GWASSS_columns(INPUT_GWAS_FILE, OUTPUT_FILE)
