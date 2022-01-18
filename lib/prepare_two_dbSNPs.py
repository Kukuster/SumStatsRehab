# standard library
import sys
import os
import time

# local
from utils import run_bash



def remove_last_ext(filename: str):
    return filename.rsplit(".", 1)[0] # passed 1 means do max 1 split; _rightmost_ splits first


def prepare_two_dbSNPs(SNPs_FILE: str, gzsort: str, bcftools: str, buffer_size: str, OUTPUT_FILE: str):
    """
    Takes a dbSNP dataset and turns it into two datasets for later use in FIX command.

    Parameters
    ----------
    SNPs_FILE : str
        path to a SNPs file in vcf, vcf.gz, bcf, bcf.gz format

    gzsort : str
        path to the gz-sort executable

    bcftools : str
        path to bcftools executable

    buffer_size : str
        buffer size for sorting (size of presort), supports k/M/G suffix

    OUTPUT_FILE : str
        base path for output file names of two prepared DBs
    """

    if not os.path.isfile(SNPs_FILE):
        print(f"ERROR: provided SNPs file doesn't exist: {SNPs_FILE}")
        exit(1)


    buffer_size = buffer_size.strip().replace(' ', '')
    SNPs_FILE_DATA = OUTPUT_FILE + ".1.tsv.gz"
    SNPs_FILE_DATA_RSID_SORTED_TMP = OUTPUT_FILE + ".2.unsorted.tsv.gz"
    SNPs_FILE_DATA_RSID_SORTED = OUTPUT_FILE + ".2.tsv.gz"


    #
    # STEP #1
    #    Get formatted table data from dbSNP
    #
    print("=== Preparing DB1 ===")
    start_time = time.time()

    get_header = f"echo 'CHROM\tPOS\tID\tREF\tALT\tFREQ'"
    query_fields = f"\"{bcftools}\" query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\tfreq=%FREQ\n' \"{SNPs_FILE}\""
    format_fields = f"""awk -F $'\t' 'BEGIN {{
        # GRCh38.p13, GenBank sequence
        chrs["CM000663.2"] = "1"
        chrs["CM000664.2"] = "2"
        chrs["CM000665.2"] = "3"
        chrs["CM000666.2"] = "4"
        chrs["CM000667.2"] = "5"
        chrs["CM000668.2"] = "6"
        chrs["CM000669.2"] = "7"
        chrs["CM000670.2"] = "8"
        chrs["CM000671.2"] = "9"
        chrs["CM000672.2"] = "10"
        chrs["CM000673.2"] = "11"
        chrs["CM000674.2"] = "12"
        chrs["CM000675.2"] = "13"
        chrs["CM000676.2"] = "14"
        chrs["CM000677.2"] = "15"
        chrs["CM000678.2"] = "16"
        chrs["CM000679.2"] = "17"
        chrs["CM000680.2"] = "18"
        chrs["CM000681.2"] = "19"
        chrs["CM000682.2"] = "20"
        chrs["CM000683.2"] = "21"
        chrs["CM000684.2"] = "22"
        chrs["CM000685.2"] = "X"
        chrs["CM000686.2"] = "Y"

        # GRCh38.p13, RefSeq sequence
        chrs["NC_000001.11"] = "1"
        chrs["NC_000002.12"] = "2"
        chrs["NC_000003.12"] = "3"
        chrs["NC_000004.12"] = "4"
        chrs["NC_000005.10"] = "5"
        chrs["NC_000006.12"] = "6"
        chrs["NC_000007.14"] = "7"
        chrs["NC_000008.11"] = "8"
        chrs["NC_000009.12"] = "9"
        chrs["NC_000010.11"] = "10"
        chrs["NC_000011.10"] = "11"
        chrs["NC_000012.12"] = "12"
        chrs["NC_000013.11"] = "13"
        chrs["NC_000014.9"]  = "14"
        chrs["NC_000015.10"] = "15"
        chrs["NC_000016.10"] = "16"
        chrs["NC_000017.11"] = "17"
        chrs["NC_000018.10"] = "18"
        chrs["NC_000019.10"] = "19"
        chrs["NC_000020.11"] = "20"
        chrs["NC_000021.9"]  = "21"
        chrs["NC_000022.11"] = "22"
        chrs["NC_000023.11"] = "X"
        chrs["NC_000024.10"] = "Y"

        # GRCh37, GenBank sequence
        chrs["CM000663.1"] = "1"
        chrs["CM000664.1"] = "2"
        chrs["CM000665.1"] = "3"
        chrs["CM000666.1"] = "4"
        chrs["CM000667.1"] = "5"
        chrs["CM000668.1"] = "6"
        chrs["CM000669.1"] = "7"
        chrs["CM000670.1"] = "8"
        chrs["CM000671.1"] = "9"
        chrs["CM000672.1"] = "10"
        chrs["CM000673.1"] = "11"
        chrs["CM000674.1"] = "12"
        chrs["CM000675.1"] = "13"
        chrs["CM000676.1"] = "14"
        chrs["CM000677.1"] = "15"
        chrs["CM000678.1"] = "16"
        chrs["CM000679.1"] = "17"
        chrs["CM000680.1"] = "18"
        chrs["CM000681.1"] = "19"
        chrs["CM000682.1"] = "20"
        chrs["CM000683.1"] = "21"
        chrs["CM000684.1"] = "22"
        chrs["CM000685.1"] = "X"
        chrs["CM000686.1"] = "Y"

        # GRCh37, RefSeq sequence
        chrs["NC_000001.10"] = "1"
        chrs["NC_000002.11"] = "2"
        chrs["NC_000003.11"] = "3"
        chrs["NC_000004.11"] = "4"
        chrs["NC_000005.9"]  = "5"
        chrs["NC_000006.11"] = "6"
        chrs["NC_000007.13"] = "7"
        chrs["NC_000008.10"] = "8"
        chrs["NC_000009.11"] = "9"
        chrs["NC_000010.10"] = "10"
        chrs["NC_000011.9"]  = "11"
        chrs["NC_000012.11"] = "12"
        chrs["NC_000013.10"] = "13"
        chrs["NC_000014.8"]  = "14"
        chrs["NC_000015.9"]  = "15"
        chrs["NC_000016.9"]  = "16"
        chrs["NC_000017.10"] = "17"
        chrs["NC_000018.9"]  = "18"
        chrs["NC_000019.9"]  = "19"
        chrs["NC_000020.10"] = "20"
        chrs["NC_000021.8"]  = "21"
        chrs["NC_000022.10"] = "22"
        chrs["NC_000023.10"] = "X"
        chrs["NC_000024.9"]  = "Y"

        # NCBI36, RefSeq sequence
        chrs["NC_000001.9"]  = "1"
        chrs["NC_000002.10"] = "2"
        chrs["NC_000003.10"] = "3"
        chrs["NC_000004.10"] = "4"
        chrs["NC_000005.8"]  = "5"
        chrs["NC_000006.10"] = "6"
        chrs["NC_000007.12"] = "7"
        chrs["NC_000008.9"]  = "8"
        chrs["NC_000009.10"] = "9"
        chrs["NC_000010.9"]  = "10"
        chrs["NC_000011.8"]  = "11"
        chrs["NC_000012.10"] = "12"
        chrs["NC_000013.9"]  = "13"
        chrs["NC_000014.7"]  = "14"
        chrs["NC_000015.8"]  = "15"
        chrs["NC_000016.8"]  = "16"
        chrs["NC_000017.9"]  = "17"
        chrs["NC_000018.8"]  = "18"
        chrs["NC_000019.8"]  = "19"
        chrs["NC_000020.9"]  = "20"
        chrs["NC_000021.7"]  = "21"
        chrs["NC_000022.9"]  = "22"
        chrs["NC_000023.9"]  = "X"
        chrs["NC_000024.8"]  = "Y"

        chrs["01"] = "1";
        chrs["02"] = "2";
        chrs["03"] = "3";
        chrs["04"] = "4";
        chrs["05"] = "5";
        chrs["06"] = "6";
        chrs["07"] = "7";
        chrs["08"] = "8";
        chrs["09"] = "9";
        chrs["x"]  = "X";
        chrs["y"]  = "Y";
        chrs["m"]  = "M";

    }}
    # note: *in* operator checks if the string is in the array keys, which is probably O(1)
    # the next block is executed on each row and ensures that:
    #   - the chr prefix is removed wherever it is present
    #   - if a field is one of the keys in the *chrs* array, it is mapped to the corresponding value
    {{
        if ($0 ~ /^chr/) {{
            print (substr($1,4) in chrs) ? chrs[substr($1,4)]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6 : substr($1,4)"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6
        }} else {{
            print ($1 in chrs) ? chrs[$1]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6 : $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6
        }}
    }}'"""
    gzip = f"gzip"
    run_bash(f"{query_fields} | {format_fields} | {gzip} > \"{SNPs_FILE_DATA}\"")

    print(f"  Preparing DB1 finished in {(time.time() - start_time)} seconds\n")


    #
    # STEP #2
    #    FINALLY, sort the formatted table data by rsID, and remove the intermediate file
    #
    print("=== Preparing DB2 ===")
    start_time = time.time()

    get_rsID_col = f"gunzip -c \"{SNPs_FILE_DATA}\" | cut -d$'\t' -f3"
    get_other_cols = f"gunzip -c \"{SNPs_FILE_DATA}\" | cut -d$'\t' -f1-2,4-6"
    run_bash(f"paste -d$'\t' <({get_rsID_col}) <({get_other_cols}) | gzip > \"{SNPs_FILE_DATA_RSID_SORTED_TMP}\"")
    run_bash(f"\"{gzsort}\" -S {buffer_size} \"{SNPs_FILE_DATA_RSID_SORTED_TMP}\" \"{SNPs_FILE_DATA_RSID_SORTED}\"")
    run_bash(f"rm \"{SNPs_FILE_DATA_RSID_SORTED_TMP}\"")

    print(f"  Preparing DB2 finished in {(time.time() - start_time)} seconds\n")



if __name__ == "__main__":
    SNPs_FILE = sys.argv[1]
    gzsort = sys.argv[2]
    bcftools = sys.argv[3]
    buffer_size = sys.argv[4]
    OUTPUT_FILE = sys.argv[5]

    prepare_two_dbSNPs(SNPs_FILE, gzsort, bcftools, buffer_size, OUTPUT_FILE)
