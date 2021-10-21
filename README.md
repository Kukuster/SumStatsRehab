# SSrehab

## how to run
...

## main dependencies:
 - python 3.8+
 - Linux system with bash v4 or 5
 - python packages in `requirements.txt`


## NOTES

### "standard" format
 - file is in tsv format, i.e. tabular tab-separated format (bare, zipped, or gzipped)
 - there's one-line header in the file on the first line. All other lines are the data entries
 - the file has precisely columns defined as `STANDARD_COLUMN_ORDER` in `lib/standard_column_order.py`.
    - file has exactly these columns, exactly this number of columns, and no other columns
    - columns are in this exact order
    - if the original file was missing a column, an empty column should be taking its place (entries are *the empty string*)




### BACKLOG
 - a config file has to be generated with all the names of the intermediary files (or does it). This will improve refactoring into the actual pipeline.
 - (maybe) improve restoring alleles by adding checks for exact match of flipped alleles, if other checks didn't help. This requires having all SNPs for a particular ChrBP in the memory, and is relevant only for restoring alleles by looping through file sorted by Chr and BP.
 - add ability to specify additional columns from the GWAS SS file that user wants to include in the end file. This would be an array in the the json config file for the input GWAS SS file.
 - improve code in the main file: `SSrehab.py`
 - improve resolver architecture in `loop_fix.py`: make a separate function loopDB1 and loopDB2 that will loop through enough entries in a DB before every resolver and rewrite a "global" object with properties to be fields from the DB: rsID, Chr, BP, alleles, EAF. So resolvers for rsID and ChrBP will be similar to ones for alleles and EAF. Resolvers for these fields then should operate on `fields` and that object with fields from a DB. This way a really strong optimization, flexibility, and modularity of resolvers will be achieved. `run_all` doesn't have to have resolvers and resolvers_args object to be passed, it can just use the global ones.

