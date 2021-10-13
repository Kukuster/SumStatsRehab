# standard library
import sys
import os

# local
from utils import run_bash


# # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                 #
#                      INPUT                      #
#                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # #

if len(sys.argv) < 4:  # the very first 0th arg is the name of this script
    print("ERROR: you should specify args:")
    print("  #1 SNPs file in vcf format")
    print("  #2 buffer size: size of presort, supports k/M/G suffix")
    print("  #3 output file name, SNPs file in vcf format, sorted by rsID")
    exit(1)

# GWAS_FILE has to be in the internal "standard" tsv format
SNPs_FILE = sys.argv[1]
buffer_size = sys.argv[2].strip().replace(' ', '')
OUTPUT_FILE = sys.argv[3]


if not os.path.isfile(SNPs_FILE):
    print(f"ERROR: provided SNPs file doesn't exist: {SNPs_FILE}")
    exit(1)


def remove_last_ext(filename: str):
    return filename.rsplit(".", 1)[0] # passed 1 means do max 1 split; _rightmost_ splits first


SNPs_FILE_nohead = remove_last_ext(OUTPUT_FILE) + "_original_nohead.vcf.gz"


# # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                 #
#                      MAIN                       #
#                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # #



#
# STEP #1
#    Get the copy of the gz file without the header
#
first_data_line_i = run_bash(f"gunzip -c \"{SNPs_FILE}\" | awk '!/^#/ {{print FNR; exit;}}'")

print("copying original SNPs file without the header, with columns reordered so rsID goes first")

get_raw_data=f"gunzip -c \"{SNPs_FILE}\" | tail -n +{first_data_line_i}"

run_bash(f"paste -d$'\\t' <({get_raw_data} | cut -d$'\\t' -f3) <({get_raw_data} | cut -d$'\\t' -f1,2,4,5) | gzip -c > \"{SNPs_FILE_nohead}\"")

#
# STEP #2
#    Sort the big gzipped file using gz-sort
#
print("sorting with gz-sort")
run_bash(f"\"$HOME/gz-sort/gz-sort\" -S {buffer_size} \"{SNPs_FILE_nohead}\" \"{OUTPUT_FILE}\"")


# # This way, sort(1) uses an unpredictably big amount of cache and fails silently with exit code 0 when there's not enough memory
# first_data_line_i = run_bash(f"gunzip -c \"{SNPs_FILE}\" | awk '!/^#/ {{print FNR; exit;}}'")
# run_bash(f'gunzip -c \"{SNPs_FILE}\" | tail -n +{first_data_line_i} | cut -d$\'\\t\' -f1-5 | LC_ALL=C sort -k3 | gzip -c > \"{OUTPUT_FILE}\"')

