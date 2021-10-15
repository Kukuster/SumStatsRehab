# standard library
import sys
import os
import time

# local
from utils import run_bash


# # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                 #
#                      INPUT                      #
#                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # #

if len(sys.argv) < 4:  # the very first 0th arg is the name of this script
    print("ERROR: you should specify args:")
    print("  #1 SNPs file in vcf, vcf.gz, bcf, bcf.gz format")
    print("  #2 path to gz-sort executable")
    print("  #3 path to bcftools executable")
    print("  #4 buffer size for sorting (size of presort), supports k/M/G suffix")
    print("  #5 base for output file names of two prepared DBs")
    exit(1)

# GWAS_FILE has to be in the internal "standard" tsv format
SNPs_FILE = sys.argv[1]
gzsort = sys.argv[2]
bcftools = sys.argv[3]
buffer_size = sys.argv[4].strip().replace(' ', '')
OUTPUT_FILE = sys.argv[5]


if not os.path.isfile(SNPs_FILE):
    print(f"ERROR: provided SNPs file doesn't exist: {SNPs_FILE}")
    exit(1)


def remove_last_ext(filename: str):
    return filename.rsplit(".", 1)[0] # passed 1 means do max 1 split; _rightmost_ splits first


SNPs_FILE_DATA = OUTPUT_FILE + ".1.tsv.gz"
SNPs_FILE_DATA_RSID_SORTED_TMP = OUTPUT_FILE + ".2.unsorted.tsv.gz"
SNPs_FILE_DATA_RSID_SORTED = OUTPUT_FILE + ".2.tsv.gz"


# # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                 #
#                      MAIN                       #
#                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # #


#
# STEP #1
#    Get formatted table data from dbSNP
#
print("=== Preparing DB1 ===")
start_time = time.time()

get_header = f"echo 'CHROM\tPOS\tID\tREF\tALT\tFREQ'"
query_fields = f"\"{bcftools}\" query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\tfreq=%FREQ\n' \"{SNPs_FILE}\""
format_fields = f"awk -F $'\t' '{{if ($0 ~ /^chr/) {{print substr($0,4)}} else {{print $0}} }}'"
gzip = f"gzip"
run_bash(f"{query_fields} | {format_fields} | {gzip} > \"{SNPs_FILE_DATA}\"")

print(f"  Preparing DB1 finished in {(time.time() - start_time)} seconds\n")


#
# STEP #2
#    Sort the formatted table data by rsID
#
print("=== Preparing DB2 ===")
start_time = time.time()

get_rsID_col = f"gunzip -c \"{SNPs_FILE_DATA}\" | cut -d$'\t' -f3"
get_other_cols = f"gunzip -c \"{SNPs_FILE_DATA}\" | cut -d$'\t' -f1-2,4-6"
run_bash(f"paste -d$'\t' <({get_rsID_col}) <({get_other_cols}) | gzip > \"{SNPs_FILE_DATA_RSID_SORTED_TMP}\"")
run_bash(f"\"{gzsort}\" -S {buffer_size} \"{SNPs_FILE_DATA_RSID_SORTED_TMP}\" \"{SNPs_FILE_DATA_RSID_SORTED}\"")
run_bash(f"rm \"{SNPs_FILE_DATA_RSID_SORTED_TMP}\"")

print(f"  Preparing DB2 finished in {(time.time() - start_time)} seconds\n")
