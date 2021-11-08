# SSrehab

## dependencies:
 - python 3.8+
 - Linux system with bash v4 or 5
 - python packages in `requirements.txt`
 - [bcftools](https://github.com/samtools/bcftools) (only for `prepare_dbSNPs`)
 - [gz-sort](http://kmkeen.com/gz-sort/) (only for `prepare_dbSNPs`)

## Installation and basics
1. clone this repo
```bash
git clone https://github.com/Kukuster/SSrehab.git
```
2. install requirements
```bash
pip install -r requirements.txt
```

3. run the eponymous script in the cloned directory using the following syntax:
```bash
python3 SSrehab.py <command> [keys]
```

Use `diagnose` to check the validity of entries in the GWAS SS file.

Use `fix` to restore missing/invalid data in the GWAS SS file.

Use `prepare_dbSNPs` to preprocess a given dbSNP dataset into 2 datasets, which are used in the `fix` command.

Use `sort` to format the input GWAS SS file and sort either by Chr and BP or by rsID.

To use the `fix` command to the its fullest, a user needs: 
 - SNPs datasets in the target build, preprocessed with the `prepare_dbSNPs` command.
 - chain file, if the GWAS SS file is provided in build different from the target build 


## Tutorial
### 1. Download dbSNP dataset
Download dbSNP datasets from ncbi, in the target build, in vcf, vcf.gz, bcf, or bcf.gz format. Latest versions are recommended.
dbSNP datasets are used to restore the following data: Chr, BP, rsID, OA, EA, EAF. Although only builds 37 and 38 are explicitly supported, build 36 may work as well.

For example, curently latest datasets for [build 38](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39/) and [build 37](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.25/) can be downloaded here:

https://ftp.ncbi.nih.gov/snp/latest_release/VCF/


### 2. Download the chain file
Chain file is necessary to perform liftover. If GWAS SS file is provided in the target build, then a chain file is not used.

### 3. Preprocess dbSNPs datasets

#### 3.1 Download and install [bcftools](https://github.com/samtools/bcftools) and [gz-sort](http://kmkeen.com/gz-sort/)
see instructions on their websites and/or githubs

recommended version for bcftools - 1.11

NOTE: after preprocessing of the necessary dbSNPs is finished, these tools are no longer needed

#### 3.2 Run preprocessing
Run `prepare_dbSNPs` using the following syntax:
```bash
python3 SSrehab.py prepare_dbSNPs --dbsnp DBSNP --OUTPUT OUTPUT --gz-sort GZ_SORT --bcftools BCFTOOLS
                                  [--buffer BUFFER]
```
where:
 - `DBSNP` is the dbSNP dataset of intented 
 - `OUTPUT` is the base name for the two output dbSNPs datasets
 - `GZ_SORT` is a path to the gz-sort executable
 - `BCFTOOLS` is a path to the bcftools executable
 - `BUFFER` is buffer size for sorting (size of presort), supports k/M/G suffix. Defaults to 1G. Recommended: at least 200M, ideally: 4G or more

Depending on the size of the dataset, specified buffer size, and specs of the machine, preprocessing may take somewhere from 30 minutes to 6 hours.

After preprocessing is finished, steps 4 and 5 may be repeated ad lib.

### 4. Create a config file for your GWAS SS file
Config file is used as meta data for GWAS SS file, and contains:
 1) columns indices starting with 0
 2) input build slug.

This config file has to have the same file name as the GWAS SS file but with additional `.json` extension.

For example, if your GWAS SS file is named `WojcikG_PMID_htn.gz`, and the first 5 lines in the unpacked file are:
```
Chr     Position_hg19   SNP     Other-allele    Effect-allele   Effect-allele-frequency Sample-size     Effect-allele-frequency-cases   Sample-size-cases       Beta    SE      P-val    INFO-score      rsid
1       10539   1:10539:C:A     C       A       0.004378065     49141   0.003603676     27123   -0.1041663      0.1686092       0.5367087       0.46    rs537182016
1       10616   rs376342519:10616:CCGCCGTTGCAAAGGCGCGCCG:C      CCGCCGTTGCAAAGGCGCGCCG  C       0.9916342       49141   0.9901789       27123   -0.1738814      0.109543        0.1124369        0.604   rs376342519
1       10642   1:10642:G:A     G       A       0.006042409     49141   0.007277901     27123   0.1794246       0.1482529       0.226179        0.441   rs558604819
1       11008   1:11008:C:G     C       G       0.1054568       49141   0.1042446       27123   -0.007140072    0.03613677      0.84337 0.5     rs575272151
```

your config file should have the name `WojcikG_PMID_htn.gz.json` and the following contents:
```json
{
    "Chr": 0,
    "BP": 1,
    "rsID": 13,
    "OA": 3,
    "EA": 4,
    "EAF": 5,
    "beta": 9,
    "SE": 10,
    "pval": 11,
    "INFO": 12,

    "build": "grch37"
}
```

Notes:
 - SSrehab will only consider data in the columns which indices are specified in the config file.
 - In this example, all the 10 columns from the list of supported columns are present. But none of the columns above are mandatory. If certain columns are missing, `fix` command will attempt to restore them if possible.


### 5. Run the `fix` command
When the config file is created, and dbSNP datasets are preprocessed, chain file is downloaded if necessary, then `fix` command may be used utilizing all its features.

Although it is normally a part of execution of `fix` command, user may choose to manually run `diagnose` beforehand.

If `diagnose` is ran without additional arguments, it is "readonly", i.e. doesn't write into the file system.

Run `diagnose` as follows:

```bash
python3 SSrehab.py diagnose --INPUT INPUT_GWAS_FILE
```
where `INPUT_GWAS_FILE` is the path to the GWAS SS file with the corresponding config file at `*.json`

as a result it will generate the main plot: stacked histogram plot, and an addtional bar chart plot for each of the bins in the stacked histogram plot.

This plots will popup in a new matplotlib window.

The stacked histogram maps the number of invalid SNPs against p-value, allowing assessment of the distribution of invalid SNPs by significance. On the histogram, valid SNPs are shown as blue, and SNPs that have issues are shown as red. Height of the red plot over each bin with the red caption represents the proportion of invalid SNPs in the corresponding bin.

![WojcikG_PMID_htn gz](https://user-images.githubusercontent.com/12045236/140248528-8ed3bc3f-b53e-4cef-af46-6d60d9899bcc.png)

A bar chart is generated for each bin of the stacked histogram plot and reports the number of issues that invalidate the SNP entries in a particular bin.

![bin_3__1e-5—1e-3](https://user-images.githubusercontent.com/12045236/140248537-28c4c287-ef84-4bb6-a9ee-30c4c879b5cd.png)

If the linux system runs without GUI, report should be saved on the file system. For this, run the command as follows:

```bash
python3 SSrehab.py diagnose --INPUT INPUT_GWAS_FILE --REPORT-DIR REPORT_DIR
```
where `REPORT_DIR` is existing or not existing directory under which the generated report will be contained. When saved onto a disk, report also includes a small table with exact numbers of invalid fields and other issues in the GWAS SS file.


Finally, user may want to decide to run `fix` command.

User should run the `fix` command as follows:
```bash
python3 SSrehab.py fix --INPUT INPUT_GWAS_FILE --OUTPUT OUTPUT_FILE
                       [--dbsnp-1 DBSNP1_FILE] [--dbsnp-2 DBSNP2_FILE]
                       [--chain-file CHAIN_FILE]
                       [--freq-db FREQ_DATABASE_SLUG]
```
where:
 - `INPUT_GWAS_FILE` is the input GWAS SS file with the corresponding `.json` config file create at step 4
 - `OUTPUT_FILE` is the base name for the fixed file(s) 
 - `DBSNP1_FILE` is a path to the preprocessed dbSNP #1
 - `DBSNP2_FILE` is a path to the preprocessed dbSNP #2
 - `CHAIN_FILE` is a path to the chain file
 - `FREQ_DATABASE_SLUG` is a population slug from a frequency database in dbSNP

example:

```bash
python3 SSrehab.py fix --INPUT "29559693.tsv" --OUTPUT "SSrehab_fixed/29559693" --dbsnp-1 "dbSNP_155_b38.1.tsv.gz" --dbsnp-2 "dbSNP_155_b38.2.tsv.gz" --chain-file "hg19_to_hg38.chain" --freq-db TOPMED
```

As the normal process of `fix`, a report will be generated for the input file, as well as for the file after each step of processing. Depending on the availability of invalid/missing data in the GWAS SS file and the input arguments, a different number of steps may be required for a complete run of `fix` command, with 1 or 2 _loops_ ran on the GWAS SS file. All steps are performed automatically without prompt. The process of `fix`ing is represented in logging to the stanard output, and may take anywhere from 5 minutes to 1.5 hours, depending on the size of the file and the number of steps.

As a result, if 1 loop was required to fix the file, then the resulting file will be available with the suffix `.SSrehabed.tsv`. If 2 loops were required, then the resulting file is available with the suffix `.SSrehabed-twice.tsv`.

The report made with a `diagnose` command will be available in a separate directory for:
 - the input file
 - for the file after 1 loop of fixing
 - for the file after 2 loops of fixing (applicable only if 2 loops were required)

<hr>

## Manual

Please refer to the instructions by running
```bash
python3 SSrehab.py -h
```
or
```bash
python3 SSrehab.py <command> -h
```


## NOTES

### "standard" format
 - file is in tsv format, i.e. tabular tab-separated format (bare, zipped, or gzipped)
 - there's one-line header in the file on the first line. All other lines are the data entries
 - the file has precisely columns defined as `STANDARD_COLUMN_ORDER` in `lib/standard_column_order.py`.
    - file has exactly these columns, exactly this number of columns, and no other columns
    - columns are in this exact order
    - if the original file was missing a column, an empty column should be taking its place (entries are *the empty string*)


### BACKLOG
 - upon execution of `fix` command, a config file has to be generated with all the names of the intermediary files. This will improve refactoring into the actual pipeline.
 - (maybe) improve restoring alleles by adding checks for exact match of flipped alleles, if other checks didn't help. This requires having all SNPs for a particular ChrBP in the memory, and is relevant only for restoring alleles by looping through file sorted by Chr and BP.
 - add ability to specify additional columns from the GWAS SS file that user wants to include in the end file. This would be an array of integers in the the json config file for the input GWAS SS file.
 - improve code in the main file: `SSrehab.py`
 - improve resolver architecture in `loop_fix.py`: make a separate function loopDB1 and loopDB2 that will loop through enough entries in a DB before every resolver and rewrite a "global" object with properties to be fields from the DB: rsID, Chr, BP, alleles, EAF. So resolvers for rsID and ChrBP will be similar to ones for alleles and EAF. Resolvers for these fields then should operate on `fields` and that object with fields from a DB. This way a really strong optimization, flexibility, and modularity of resolvers will be achieved. `run_all` doesn't have to have resolvers and resolvers_args object to be passed, it can just use the global ones.
 - improve interface for liftover. SSrehab fix should work for all sorts of liftover between builds 36, 37, and 38, including back liftover. If user omits the preprocessed dbSNP databases as input, but specifies the chain file, it can perform liftover only.
 - introduce dependency on the STANDARD_COLUMN_ORDER in the `validate_GWASSS_entries.py` script
 - add support for OR, and, maybe, restoring OR from beta or vice versa.
 - add a keyword argument that will cause SSrehab fix to clean up all intermediate files and leave only the last resulting file after the processing.
 - add a keyword argument that specifies tmp directory for intermediate files. GWAS SS files are usually 1-4 Gigs unpacked.
 - set alleles column to uppercase during preparation (in `prepare_GWASSS_columns.py` script).
 - add a keyword argument that will cause SSrehab `fix` to prompt for fixing after the diagnosis step
