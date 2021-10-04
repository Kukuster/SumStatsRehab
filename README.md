# SSrehab

## how to run
run `SSrehab.py` with `python3`. It accepts 3 arguments:
 1. Path to GWAS summary statistics file in tsv format, that has a corresponding config file (suffixed \".json\") with column indices. Use `sample/config.example.json` as a template. No columns are required to specify. If you don't specify a column, `SSrehab` will attempt to fully restore it, when possible.
 2. Path to the output fixed GWAS SS file
 3. Path to a dbSNP file that corresponds to the input GWAS SS file

e.g.:
```python
python3 SSrehab.py "sample/29559693_randlines50000.tsv" "sample/29559693_randlines50000_SSREHAB-FIXED.tsv" "/media/$USER/exFAT_share/SelfDecode/dbSNP151_GRCh37.vcf.gz"
# this implies that config file exists at: "sample/29559693_randlines50000.tsv.json"
```

## main dependencies:
 - python3
 - Linux system with bash v4 or 5


## NOTES

### "standard" format
 - file is in tsv format, i.e. tabular tab-separated format (bare, zipped, or gzipped)
 - there's one-line header in the file on the first line. All other lines are the data entries
 - the file has precisely columns defined as `STANDARD_COLUMN_ORDER` in `lib/standard_column_order.py`.
    - file has exactly these columns, exactly this number of columns, and no other columns
    - columns are in this exact order
    - if the original file was missing a column, an empty column should be taking its place (entries are *the empty string*)




### TODO:
 - add requirements.txt and things to easily recreate the environment
 - improve interface. Use argparse or something like that for convenient CLI arguments
 - add conditional statements that will check cases when all entries in particular columns are missing. If vital columns are missing for restoring something, don't attempt to restore it and provide a warning message
 - add resolvers for Chr&BP, alleles, and MAF
 - improve comments in `lib/loop_fix.py`

