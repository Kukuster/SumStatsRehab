# SumStatsRehab

SumStatsRehab is a universal GWAS SumStats [pre-processing](https://en.wikipedia.org/wiki/Data_preparation) tool. SumStatsRehab takes care of each of the original data points to maximize statistical power of downstream calculations. Currently, the only supported processing which may result in a loss of the original data is liftover, which is a common task, and is optional to the user.

Examples of what the tool does:
 - data validation (e.g. `diagnose` command),
 - data enrichment (e.g. restoration of a missing Chr and BP fields from rsID in the input GWAS SumStats file by rsID lookup in the input dbSNP dataset),
 - data correction (restoration of invariably erroneous values),
 - data formatting (e.g. sorting),
 - data restoration (e.g. calculation of a StdErr field from present p-value and beta fields).

Example of what the tool does not:
 - data cleaning/reduction

SumStatsRehab aims to be a production-grade software, but you may find it to be not a complete data preparation solution for sumstats-ingesting pipelines. Yet, the tool streamlines the development of the data preparation part. This comes out from focusing on solving more complex problems that span many different use-cases, instead of covering only specific use-cases.

## dependencies:
 - python 3.8+
 - a GNU/Linux with bash v4 or 5.
 - several python packages in `requirements.txt`
 - [bcftools](https://github.com/samtools/bcftools) (only for the `prepare_dbSNPs` command)
 - [gz-sort](http://kmkeen.com/gz-sort/) (only for the `prepare_dbSNPs` command)

## Installation and basics
1. clone this repo
```bash
git clone https://github.com/Kukuster/SumStatsRehab.git && cd SumStatsRehab
```

2. install requirements
```bash
pip install -r requirements.txt
```

3. install SumStatsRehab as a package
```bash
python3 setup.py build
python3 setup.py install
```

4. run the command using the following syntax:
```bash
SumStatsRehab <command> [keys]
```

Use `diagnose` to check the validity of entries in the GWAS SS file.

Use `fix` to restore missing/invalid data in the GWAS SS file.

Use `prepare_dbSNPs` to preprocess a given dbSNP dataset into 2 datasets, which are used in the `fix` command.

Use `sort` to format the input GWAS SS file and sort either by Chr and BP or by rsID.

To use the `fix` command to its fullest, a user needs: 
 - SNPs datasets in the target build, preprocessed with the `prepare_dbSNPs` command.
 - chain file, if the GWAS SS file is provided in build different from the target build 


## Deprecated installation method (since March 15, 2022):
```
pip install git+https://github.com/Kukuster/SumStatsRehab.git
```
or, for specific version:
```
pip install git+https://github.com/Kukuster/SumStatsRehab.git@v1.2.0 --upgrade
```

This installation method doesn't work with the currently upcoming git protocol security update on github:
 - https://github.blog/2021-09-01-improving-git-protocol-security-github/




## Tutorial
### 1. Download dbSNP dataset
Download dbSNP datasets from NCBI, in the target build, in vcf, vcf.gz, bcf, or bcf.gz format. The latest versions are recommended.
dbSNP datasets are used to restore the following data: Chr, BP, rsID, OA, EA, EAF. Although only builds 37 and 38 are explicitly supported, build 36 may work as well.

For example, curently latest datasets for [build 38](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39/) and [build 37](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.25/) can be downloaded here:

https://ftp.ncbi.nih.gov/snp/latest_release/VCF/


### 2. Download the chain file
A chain file is necessary to perform liftover. If a GWAS SS file is provided in the target build, then a chain file is not used.

### 3. Preprocess dbSNPs datasets

#### 3.1 Download and install [bcftools](https://github.com/samtools/bcftools) and [gz-sort](http://kmkeen.com/gz-sort/)
see instructions on their websites and/or githubs

recommended bcftools version: 1.11

NOTE: after preprocessing of the necessary dbSNPs is finished, these tools are no longer needed

#### 3.2 Run preprocessing
Run `prepare_dbSNPs` using the following syntax:
```bash
SumStatsRehab prepare_dbSNPs --dbsnp DBSNP --OUTPUT OUTPUT --gz-sort GZ_SORT --bcftools BCFTOOLS
                                  [--buffer BUFFER]
```
where:
 - `DBSNP` is the dbSNP dataset in vcf, vcf.gz, bcf, or bcf.gz format referencing build 38 or 37
 - `OUTPUT` is the base name for the two output dbSNPs datasets
 - `GZ_SORT` is a path to the gz-sort executable
 - `BCFTOOLS` is a path to the bcftools executable
 - `BUFFER` is buffer size for sorting (size of presort), supports k/M/G suffix. Defaults to 1G. Recommended: at least 200M, ideally: 4G or more

Depending on the size of the dataset, specified buffer size, and specs of the machine, preprocessing may take somewhere from 30 minutes to 6 hours.

After preprocessing, steps 4 and 5 may be repeated ad-lib.

### 4. Create a config file for your GWAS SS file
Config file is used as meta data for GWAS SS file, and contains:
 1) columns' indices (indices start from 0)
 2) input build slug (such as "GRCh38", "GRCh37", "hg18", "hg19")

This config file has to have the same file name as the GWAS SS file but with an additional `.json` extension.

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
    "N": 6,

    "build": "grch37"
}
```

Notes:
 - SumStatsRehab will only consider data from the columns which indices are specified in the config file. If one of the above columns is present in the SS file but wasn't specified in the config file, then SumStatsRehab treats the column as missing.
 - In this example, all the 11 columns from the list of supported columns are present. Yet, none of the columns above are mandatory. If certain columns are missing, the `fix` command will attempt to restore them.


### 5. Run the `fix` command
When the config file is created, and dbSNP datasets are preprocessed, the chain file is downloaded if necessary, then the `fix` command can use all its features.

Although it is normally a part of the execution of the `fix` command, a user may choose to manually run `diagnose` beforehand.

If `diagnose` is ran without additional arguments, it is "read-only", i.e. doesn't write into the file system.

Run `diagnose` as follows:

```bash
SumStatsRehab diagnose --INPUT INPUT_GWAS_FILE
```
where `INPUT_GWAS_FILE` is the path to the GWAS SS file with the corresponding config file at `*.json`

as a result, it will generate the main plot: stacked histogram plot, and an additional bar chart plot for each of the bins in the stacked histogram plot.

These plots will pop up in a new matplotlib window.

The stacked histogram maps the number of invalid SNPs against p-value, allowing assessment of the distribution of invalid SNPs by significance. On the histogram, valid SNPs are shown as blue, and SNPs that have issues are shown as red. The height of the red plot over each bin with the red caption represents the proportion of invalid SNPs in the corresponding bin.

![WojcikG_PMID_htn gz](https://user-images.githubusercontent.com/12045236/140248528-8ed3bc3f-b53e-4cef-af46-6d60d9899bcc.png)

A bar chart is generated for each bin of the stacked histogram plot and reports the number of issues that invalidate the SNP entries in a particular bin.

![bin_3__1e-5—1e-3](https://user-images.githubusercontent.com/12045236/140248537-28c4c287-ef84-4bb6-a9ee-30c4c879b5cd.png)

If a Linux system runs without GUI, the report should be saved on the file system. For this, run the command as follows:

```bash
SumStatsRehab diagnose --INPUT INPUT_GWAS_FILE --REPORT-DIR REPORT_DIR
```
where `REPORT_DIR` is an existing or not existing directory under which the generated report will be contained. When saved onto a disk, the report also includes a small table with exact numbers of invalid fields and other issues in the GWAS SS file.


Finally, a user may want to decide to run the `fix` command.

A user should run the `fix` command as follows:
```bash
SumStatsRehab fix --INPUT INPUT_GWAS_FILE --OUTPUT OUTPUT_FILE
                       [--dbsnp-1 DBSNP1_FILE] [--dbsnp-2 DBSNP2_FILE]
                       [--chain-file CHAIN_FILE]
                       [--freq-db FREQ_DATABASE_SLUG]
                       [{--restore,--do-not-restore} {ChrBP,rsID,OA,EA,EAF,beta,SE,pval}+]
```
where:
 - `INPUT_GWAS_FILE` is the input GWAS SS file with the corresponding `.json` config file create at step 4
 - `OUTPUT_FILE` is the base name for the fixed file(s) 
 - `DBSNP1_FILE` is a path to the preprocessed dbSNP #1
 - `DBSNP2_FILE` is a path to the preprocessed dbSNP #2
 - `CHAIN_FILE` is a path to the chain file
 - `FREQ_DATABASE_SLUG` is a slug of a frequency database contained in the dbSNP

example:

```bash
SumStatsRehab fix --INPUT "29559693.tsv" --OUTPUT "SumStatsRehab_fixed/29559693" --dbsnp-1 "dbSNP_155_b38.1.tsv.gz" --dbsnp-2 "dbSNP_155_b38.2.tsv.gz" --chain-file "hg19_to_hg38.chain" --freq-db TOPMED --do-not-restore OA EA
```

As the normal process of `fix`, a report will be generated for the input file, as well as for the file after each step of processing. Depending on the availability of invalid/missing data in the GWAS SS file and the input arguments, a different number of steps may be required for a complete run of the `fix` command, with 1 or 2 _loops_ performed on the GWAS SS file. All steps are performed automatically without prompt. The process of `fix`ing is represented in logging to the standard output and may take anywhere from 5 minutes to 1.5 hours, depending on the size of the file and the number of steps.

As a result, if 1 loop was required to fix the file, then the resulting file will be available with the suffix `.rehabed.tsv`. If 2 loops were required, then the resulting file is available with the suffix `.rehabed-twice.tsv`.

The report made with a `diagnose` command will be available in a separate directory for:
 - the input file
 - for the file after 1 loop of fixing
 - for the file after 2 loops of fixing (applicable only if 2 loops were required)

<hr>

## Manual

Please refer to the instructions by running
```bash
SumStatsRehab -h
```
or
```bash
SumStatsRehab <command> -h
```


## NOTES

### config file
Config file is a json object which supports the following properties:
 - `"Chr"`
 - `"BP"`
 - `"rsID"`
 - `"EA"`
 - `"OA"`
 - `"EAF"`
 - `"beta"`
 - `"SE"`
 - `"pval"`
 - `"N"`
 - `"INFO"`
should be evaluated to an integer: the corresponding column index starting from 0.

And 
 - `"build"`
should be evaluated to either one of the following (case insensitive): `'hg38'`, `'GRCh38'`, `'hg19'`, `'GRCh37'`, `'hg18'`, `'ncbi36'`.
 - `"other"`
should be evaluated to an array of integers: indicies of other columns to include


It is also possible to set `EAF` to a weighted average of multiple colums. E.g. if there are separate freq. columns for case and control groups, and average freq. is needed, number of participants in each group will serve as weights for the two columns:
```json
{
    ...
    "EAF": {
        "4": 1001,
        "5": 2500
    },
    ...
}
```


During the `fix` command, the input sumstats file may undergo sorting. If you want any other columns to be included in the resulting fixed file, add 0-indexed column indices in an array as the `"other"` parameter in the config.

E.g.:
```json
{
    ...
    "other": [7, 8, 2],
    ...
}
```

With this entry, the resulting file will have 3 additional columns at the end in the order as their indices appear in this config entry.




### "standard" format
Commands `fix` and `sort` always output files in this format. Internally, input files are always converted before undergoing any processing. `lib/prepare_GWASSS_columns.py` is responsible for formatting the input.

 - file is in the tsv format, i.e. tabular tab-separated format (bare, zipped, or gzipped)
 - there's a one-line header in the file on the first line. All other lines are the data entries
 - the file has columns defined as `STANDARD_COLUMN_ORDER` in `lib/standard_column_order.py`.
    - the file always has these columns
    - these columns necessarily go first and appear in this order, therefore have exactly these indices
    - if the original file was missing a column, an empty column should be taking its place (entries are *the empty string*)
 - the file may have any number of any other columns going after the columns defined in the `STANDARD_COLUMN_ORDER`


## Supplementary Information
The two key functions of SumStatsRehab are validation and restoration, implemented for 9 data categories: chromosome, base pair position, rsID, effect allele, other allele, allele frequency, standard error, beta, p-value. Each column is validated independently of the others, and is regarded to have two possible states: valid and invalid. With the exception of the case where only one of either the chromosome or base pair position is a valid entry, valid entries are always kept and invalid entries are subject to restoration. 
1. Value in the chromosome field is considered to be valid if and only if it includes one of the following values: '1', '01', '2', '02', '3', '03', '4', '04', '5', '05', '6', '06', '7', '07', '8', '08', '9', '09','10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', 'X', 'x', 'Y', 'y', 'M', 'm'. In practice, ‘chr1’, ‘chr01’, ‘Chr1’,‘Chr01’, ‘1’, and ‘01’ are all similarly recognized as referring to chromosome 1. Any entries which do not include a specific chromosome number reference are subject to restoration.
2. Base pair position entry is considered to be valid if and only if it is a non-negative integer.
3. rsID entry is considered valid if and only if it is an arbitrary integer value with ‘rs’ prefix.
4. Effect allele and other allele entries are considered to be valid if and only if the value is either a dash (representing a deleted allele), or any non-empty combination of letters A, T, C, and G.
5. Any floating-point value between 0 and 1 inclusively is considered to be a valid allele frequency entry and p-value entry.
6. Any floating-point value is considered to be valid standard error entry and beta entry.

While most of the columns are assessed independently, the restoration process often involves multiple columns.
1. Chromosome and base pair position entries are only restored together for any particular row. If any of the two values are invalid or missing, and rsID entry in the row is present and valid, then the two columns are subject to restoration. This is the only case where a single potentially valid entry can be rewritten or removed as a result of restoration. Chromosome and base pair position entries are restored by a lookup in the preprocessed DB2 for a matching rsID entry. If a match wasn’t found, both chromosome and base pair position entries are replaced with a dot ‘.’, representing a missing value.
2. Similar to the current implementation of MungeSumstats, SumStatsRehab restores rsID entry by chromosome and base pair position. If rsID entry is missing or invalid, but both chromosome and base pair position are valid, then SumStatsRehab attempts to restore rsID by a lookup in the preprocessed DB1 for matching chromosome and base pair position entries. If a match wasn’t found, rsID is set to a dot ‘.’ This requires chromosome and base pair position entries in the row to be present and valid.
3. EA and OA entries are restored in the same way. If one of the two alleles is invalid or missing and the other one is present, and either rsID or both chromosome and base pair position are valid, then the missing allele can be restored. First, if rsID is present, then the missing allele is restored by a lookup in the preprocessed DB2 for matching rsID entry and the present allele. When rsID is not present, but Chr and BP are present, then the missing allele is restored by a lookup in the preprocessed DB1 for matching Chr and BP entries and the present allele. Otherwise, or if a lookup fails, then no entries are rewritten.
4. Allele frequency is restored as the effect allele frequency, and it’s done in a similar way: first by a lookup in DB1 for matching rsID and EA, if rsID is present; second by a lookup into DB2 for matching Chr and BP and EA. Restoring allele frequency requires frequency entries in the INFO field in the SNPs dataset that was preprocessed. By default, frequencies are taken from the database of Genotypes and Phenotypes (dbGaP), but an alternative database contained in the dbSNP dataset may be used as an argument in the FIX command.
5. Each of the columns, standard error, beta, and p-value, are restored from the other two, using the following relation: *__s = β/z__*, where *__s__* is the standard error, *__β__* is beta, and *__z__* is z-score that corresponds to the p-value in the two-tailed test [(Zhu et al. 2016)](https://www.nature.com/articles/ng.3538). If standard error is signed, then beta can be restored accurately, but if standard error is provided as an absolute value, the restored beta will be correct only as an absolute value, with an unknown sign. Files restored in this manner will have little utility in downstream applications, but may be useful in comparing relative effect sizes. 

Warning: Restored betas with an unknown sign should not be utilized in any downstream application. 


## BACKLOG
 - **add switches to the `fix` command to turn on and off the restoration of the ambiguous effect alleles (in cases when the dbSNPs present multiple options for ALT)**
 - (maybe) improve restoring alleles by adding checks for an exact match of flipped alleles if other checks didn't help. This requires having all SNPs for a particular ChrBP in the memory and is relevant only for restoring alleles by looping through the file sorted by Chr and BP.
 - **improve code in the main file: `SumStatsRehab.py`**
 - improve resolver architecture in `loop_fix.py`: make a separate function loopDB1 and loopDB2 that will loop through enough entries in a DB before every resolver and rewrite a "global" object with properties to be fields from the DB: rsID, Chr, BP, alleles, EAF. So resolvers for rsID and ChrBP will be similar to ones for alleles and EAF. Resolvers for these fields then should operate on `fields` and that object with fields from a DB. This way a really strong optimization, flexibility, and modularity of resolvers will be achieved. `run_all` doesn't have to have resolvers and resolvers_args object to be passed, it can just use the global ones.
 - improve the interface for liftover. SumStatsRehab fix should work for all sorts of liftovers between builds 36, 37, and 38, including back liftover. If the user omits the preprocessed dbSNP databases as input but specifies the chain file, it can perform liftover only.
 - add support for:
   - OR, and maybe restoration of OR from beta and vice versa
   - Z-score, and maybe restoration of z-score from p-value and vice versa
   - standard deviation, and maybe restoration of std.err., std.dev., and N from each other in accord to the relation
 - add a keyword argument that specifies a temp directory for intermediate files. GWAS SS files are usually 1-4 Gigs unpacked.
 - feature: save a human-readable textual report about the overall results of restoration (e.g. "performed a liftover, n rsIDs restored, n Chrs lost, ...")
 - at the moment of 2021.11.14, the following executables are assumed to be available in PATH: `bash`, `cut`, `paste`, `sort`, `awk`, `gzip`, `gunzip`, `head`, `tail`, `rm`, `wc`. Need to do more tests of SumStatsRehab with a different versions of `bash`, `awk` (including `gawk`, `nawk`, `mawk`. E.g. even though `gawk` is default for GNU/Linux, Ubuntu has `mawk` by default).
 - **make SumStatsRehab installable via `pip`**
 - Study what is a better approach to restoring EAF from other dbs. Bc for other populations there are not a lot of snps having frequency data. When you try to restore eaf for a more specific populations, it will miss a lot of snps in the ss files, therefore reducing overall accuracy. Idea for workaround: ability to specify multiple dbs, so the each next one in a list will be a lower priority.
 - `diagnose` command script should also generate a bar chart for the bin with missing p-value
 - **add data to the csv report of issues: rsID, rsID_restorable, Chr, Chr_restorable, ... etc. Use this report in the FIX command when deciding on the workflow.**
 - add feature: `amputate` command. Runs `diagnosis` and removes all the invalid rows.
 - improve try & catch clauses in `lib/loop_fix.py`: it has to catch speicifc Exceptions everywhere
 - catch KeyboardInterrupt exception in main, and, first, don't show the python traceback, second, call some kind of "destruct method", removing temp files. E.g. if for the `fix` command `--verbose` key was not set, then remove the intermediate files.
 - catch JSONDecodeError exception and provide a message with something like "sorry your config file is not a valid json", and provide a link to the github or the JSON docs
 - catch exceptions raise because of other wrong user input, e.g. wrong `gz-sort` or `bcftools` executable.
 - change the `fix` command interface: if `--verbose` is set, then forbid the `--OUTPUT` argument and allow `--OUTPUT-PREFIX`.
 - add feature: if the `fix` command was not set to `--verbose`, then if the `--OUTPUT` path ends with `.gz` or `.zip`, then compress the resulting file on the output.



## I accept donations!

### Paypal

<p>
<!--   <a href="https://www.paypal.com/donate/?hosted_button_id=485PXFAM75G4E">
      <img src="https://www.paypalobjects.com/en_US/i/btn/btn_donateCC_LG.gif" alt="paypal">
  </a> -->
  <a href="https://www.paypal.com/donate/?hosted_button_id=485PXFAM75G4E">
      <img src="https://www.paypalobjects.com/en_US/i/btn/btn_donate_SM.gif" alt="paypal">
  </a>
</p>

### Cryptocurrency

You can add a transaction message with the name of a project or a custom message if your wallet and the blockchain support this

Preferred blockchains:

blockchain | address |  
--- | --- | ---
<a href="javascript:void(0)" style="cursor: default;" alt="Donate via Bitcoin"><img src="https://img.shields.io/badge/-Bitcoin-402607?logo=data:image/svg%2bxml;base64,PHN2ZyBmaWxsPSIjRjc5MzFBIiByb2xlPSJpbWciIHZpZXdCb3g9IjAgMCAyNCAyNCIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj48dGl0bGU+Qml0Y29pbjwvdGl0bGU+PHBhdGggZD0iTTIzLjYzOCAxNC45MDRjLTEuNjAyIDYuNDMtOC4xMTMgMTAuMzQtMTQuNTQyIDguNzM2QzIuNjcgMjIuMDUtMS4yNDQgMTUuNTI1LjM2MiA5LjEwNSAxLjk2MiAyLjY3IDguNDc1LTEuMjQzIDE0LjkuMzU4YzYuNDMgMS42MDUgMTAuMzQyIDguMTE1IDguNzM4IDE0LjU0OHYtLjAwMnptLTYuMzUtNC42MTNjLjI0LTEuNTktLjk3NC0yLjQ1LTIuNjQtMy4wM2wuNTQtMi4xNTMtMS4zMTUtLjMzLS41MjUgMi4xMDdjLS4zNDUtLjA4Ny0uNzA1LS4xNjctMS4wNjQtLjI1bC41MjYtMi4xMjctMS4zMi0uMzMtLjU0IDIuMTY1Yy0uMjg1LS4wNjctLjU2NS0uMTMyLS44NC0uMmwtMS44MTUtLjQ1LS4zNSAxLjQwN3MuOTc1LjIyNS45NTUuMjM2Yy41MzUuMTM2LjYzLjQ4Ni42MTUuNzY2bC0xLjQ3NyA1LjkyYy0uMDc1LjE2Ni0uMjQuNDA2LS42MTQuMzE0LjAxNS4wMi0uOTYtLjI0LS45Ni0uMjRsLS42NiAxLjUxIDEuNzEuNDI2LjkzLjI0Mi0uNTQgMi4xOSAxLjMyLjMyNy41NC0yLjE3Yy4zNi4xLjcwNS4xOSAxLjA1LjI3M2wtLjUxIDIuMTU0IDEuMzIuMzMuNTQ1LTIuMTljMi4yNC40MjcgMy45My4yNTcgNC42NC0xLjc3NC41Ny0xLjYzNy0uMDMtMi41OC0xLjIxNy0zLjE5Ni44NTQtLjE5MyAxLjUtLjc2IDEuNjgtMS45M2guMDF6bS0zLjAxIDQuMjJjLS40MDQgMS42NC0zLjE1Ny43NS00LjA1LjUzbC43Mi0yLjljLjg5Ni4yMyAzLjc1Ny42NyAzLjMzIDIuMzd6bS40MS00LjI0Yy0uMzcgMS40OS0yLjY2Mi43MzUtMy40MDUuNTVsLjY1NC0yLjY0Yy43NDQuMTggMy4xMzcuNTI0IDIuNzUgMi4wODR2LjAwNnoiLz48L3N2Zz4=" /></a> |  `bc1pjd2c4xcgq978979htc9admycue4nqqhda3vwsc38agked8yya50qz454xc` | 
<a href="javascript:void(0)" style="cursor: default;" alt="Donate via Ethereum"><img src="https://img.shields.io/badge/-Ethereum-6784c7?logo=data:image/svg%2bxml;base64,PHN2ZyBmaWxsPSIjM0MzQzNEIiByb2xlPSJpbWciIHZpZXdCb3g9IjAgMCAyNCAyNCIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj48dGl0bGU+RXRoZXJldW08L3RpdGxlPjxwYXRoIGQ9Ik0xMS45NDQgMTcuOTdMNC41OCAxMy42MiAxMS45NDMgMjRsNy4zNy0xMC4zOC03LjM3MiA0LjM1aC4wMDN6TTEyLjA1NiAwTDQuNjkgMTIuMjIzbDcuMzY1IDQuMzU0IDcuMzY1LTQuMzVMMTIuMDU2IDB6Ii8+PC9zdmc+" /></a> |  `0x176D1b6c3Fc1db5f7f967Fdc735f8267cCe741F3` | <span>![Tether](https://raw.githubusercontent.com/Kukuster/Kukuster/refs/heads/master/tether_20x20.svg)</span> supports USDT ERC-20
<a href="javascript:void(0)" style="cursor: default;" alt="Donate via TRON"><img src="https://img.shields.io/badge/-TRON-5C0E0E?logo=data:image/svg%2bxml;base64,PHN2ZyBmaWxsPSIjRkYwNjBBIiBpZD0iQ2FscXVlXzEiIGRhdGEtbmFtZT0iQ2FscXVlIDEiIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyIgdmlld0JveD0iMCAwIDY0IDY0Ij48ZGVmcz48c3R5bGU+LmNscy0xe2ZpbGw6I2ZmMDYwYTt9PC9zdHlsZT48L2RlZnM+PHRpdGxlPnRyb248L3RpdGxlPjxnIGlkPSJ0cm9uIj48cGF0aCBjbGFzcz0iY2xzLTEiIGQ9Ik02MS41NSwxOS4yOGMtMy0yLjc3LTcuMTUtNy0xMC41My0xMGwtLjItLjE0YTMuODIsMy44MiwwLDAsMC0xLjExLS42MmwwLDBDNDEuNTYsNywzLjYzLS4wOSwyLjg5LDBhMS40LDEuNCwwLDAsMC0uNTguMjJMMi4xMi4zN2EyLjIzLDIuMjMsMCwwLDAtLjUyLjg0bC0uMDUuMTN2LjcxbDAsLjExQzUuODIsMTQuMDUsMjIuNjgsNTMsMjYsNjIuMTRjLjIuNjIuNTgsMS44LDEuMjksMS44NmguMTZjLjM4LDAsMi0yLjE0LDItMi4xNFM1OC40MSwyNi43NCw2MS4zNCwyM2E5LjQ2LDkuNDYsMCwwLDAsMS0xLjQ4QTIuNDEsMi40MSwwLDAsMCw2MS41NSwxOS4yOFpNMzYuODgsMjMuMzcsNDkuMjQsMTMuMTJsNy4yNSw2LjY4Wm0tNC44LS42N0wxMC44LDUuMjZsMzQuNDMsNi4zNVpNMzQsMjcuMjdsMjEuNzgtMy41MS0yNC45LDMwWk03LjkxLDcsMzAuMywyNiwyNy4wNiw1My43OFoiLz48L2c+PC9zdmc+" /></a> | `TMuNqEgEeBQ2GseWsqgaSdbtqasnJi8ePw` | <span>![Tether](https://raw.githubusercontent.com/Kukuster/Kukuster/refs/heads/master/tether_20x20.svg)</span> supports USDT TRC-20



<details>
  <summary>Alternative options (Ethereum L2, LN, EVM)</summary>
  
  blockchain | address
  --- | ---
  <a href="javascript:void(0)" style="cursor: default;" alt="Donate via Polygon"><img src="https://img.shields.io/badge/-Polygon-2a0c60?logo=data:image/svg%2bxml;base64,PHN2ZyBmaWxsPSIjN0IzRkU0IiByb2xlPSJpbWciIHZpZXdCb3g9IjAgMCAyNCAyNCIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj48dGl0bGU+UG9seWdvbjwvdGl0bGU+PHBhdGggZD0ibTE3LjgyIDE2LjM0MiA1LjY5Mi0zLjI4N0EuOTguOTggMCAwIDAgMjQgMTIuMjFWNS42MzVhLjk4Ljk4IDAgMCAwLS40ODgtLjg0NmwtNS42OTMtMy4yODZhLjk4Ljk4IDAgMCAwLS45NzcgMEwxMS4xNSA0Ljc4OWEuOTguOTggMCAwIDAtLjQ4OS44NDZ2MTEuNzQ3TDYuNjcgMTkuNjg2bC0zLjk5Mi0yLjMwNHYtNC42MWwzLjk5Mi0yLjMwNCAyLjYzMyAxLjUyVjguODk2TDcuMTU4IDcuNjU4YS45OC45OCAwIDAgMC0uOTc3IDBMLjQ4OCAxMC45NDVhLjk4Ljk4IDAgMCAwLS40ODguODQ2djYuNTczYS45OC45OCAwIDAgMCAuNDg4Ljg0N2w1LjY5MyAzLjI4NmEuOTgxLjk4MSAwIDAgMCAuOTc3IDBsNS42OTItMy4yODZhLjk4Ljk4IDAgMCAwIC40ODktLjg0NlY2LjYxOGwuMDcyLS4wNDEgMy45Mi0yLjI2MyAzLjk5IDIuMzA1djQuNjA5bC0zLjk5IDIuMzA0LTIuNjMtMS41MTd2My4wOTJsMi4xNCAxLjIzNmEuOTgxLjk4MSAwIDAgMCAuOTc4IDB2LS4wMDFaIi8+PC9zdmc+" /></a> |  `0x176D1b6c3Fc1db5f7f967Fdc735f8267cCe741F3`
  <a href="javascript:void(0)" style="cursor: default;" alt="Donate via Base"><img src="https://img.shields.io/badge/-Base-152846?logo=data:image/svg%2bxml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCA2NCA2NCIgeG1sOnNwYWNlPSJwcmVzZXJ2ZSIgd2lkdGg9IjMwcHgiIGhlaWdodD0iMzBweCI+PHBhdGggZmlsbD0iI0ZGRkZGRiIgZD0iTTYzLjYgMzJjMCAxNy40LTE0LjIgMzEuNi0zMS42IDMxLjZDMTUuNSA2My42IDEuOSA1MC45LjUgMzQuN2g0MS43di01LjNILjVDMS45IDEzLjEgMTUuNS40IDMyIC40IDQ5LjUuNCA2My42IDE0LjYgNjMuNiAzMnoiPjwvcGF0aD48L3N2Zz4=" /></a> |  `0x176D1b6c3Fc1db5f7f967Fdc735f8267cCe741F3`
  <a href="javascript:void(0)" style="cursor: default;" alt="Donate via Arbitrum"><img src="https://img.shields.io/badge/-Arbitrum-3F3F3F?logo=data:image/svg%2bxml;base64,PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0iVVRGLTgiPz4KPHN2ZyB4bWxuczp4b2RtPSJodHRwOi8vd3d3LmNvcmVsLmNvbS9jb3JlbGRyYXcvb2RtLzIwMDMiIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyIgeG1sbnM6eGxpbms9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkveGxpbmsiIHZlcnNpb249IjEuMSIgaWQ9IkxheWVyXzEiIHg9IjBweCIgeT0iMHB4IiB2aWV3Qm94PSIwIDAgMjUwMCAyNTAwIiBzdHlsZT0iZW5hYmxlLWJhY2tncm91bmQ6bmV3IDAgMCAyNTAwIDI1MDA7IiB4bWw6c3BhY2U9InByZXNlcnZlIj4KPHN0eWxlIHR5cGU9InRleHQvY3NzIj4KCS5zdDB7ZmlsbDpub25lO30KCS5zdDF7ZmlsbDojMjEzMTQ3O30KCS5zdDJ7ZmlsbDojMTJBQUZGO30KCS5zdDN7ZmlsbDojOURDQ0VEO30KCS5zdDR7ZmlsbDojRkZGRkZGO30KPC9zdHlsZT4KPGcgaWQ9IkxheWVyX3gwMDIwXzEiPgoJPGcgaWQ9Il8yNDA1NTg4NDc3MjMyIj4KCQk8cmVjdCBjbGFzcz0ic3QwIiB3aWR0aD0iMjUwMCIgaGVpZ2h0PSIyNTAwIj48L3JlY3Q+CgkJPGc+CgkJCTxnPgoJCQkJPHBhdGggY2xhc3M9InN0MSIgZD0iTTIyNiw3NjB2OTgwYzAsNjMsMzMsMTIwLDg4LDE1Mmw4NDksNDkwYzU0LDMxLDEyMSwzMSwxNzUsMGw4NDktNDkwYzU0LTMxLDg4LTg5LDg4LTE1MlY3NjAgICAgICBjMC02My0zMy0xMjAtODgtMTUybC04NDktNDkwYy01NC0zMS0xMjEtMzEtMTc1LDBMMzE0LDYwOGMtNTQsMzEtODcsODktODcsMTUySDIyNnoiPjwvcGF0aD4KCQkJCTxnPgoJCQkJCTxnPgoJCQkJCQk8Zz4KCQkJCQkJCTxwYXRoIGNsYXNzPSJzdDIiIGQ9Ik0xNDM1LDE0NDBsLTEyMSwzMzJjLTMsOS0zLDE5LDAsMjlsMjA4LDU3MWwyNDEtMTM5bC0yODktNzkzQzE0NjcsMTQyMiwxNDQyLDE0MjIsMTQzNSwxNDQweiI+PC9wYXRoPgoJCQkJCQk8L2c+CgkJCQkJCTxnPgoJCQkJCQkJPHBhdGggY2xhc3M9InN0MiIgZD0iTTE2NzgsODgyYy03LTE4LTMyLTE4LTM5LDBsLTEyMSwzMzJjLTMsOS0zLDE5LDAsMjlsMzQxLDkzNWwyNDEtMTM5TDE2NzgsODgzVjg4MnoiPjwvcGF0aD4KCQkJCQkJPC9nPgoJCQkJCTwvZz4KCQkJCTwvZz4KCQkJCTxnPgoJCQkJCTxwYXRoIGNsYXNzPSJzdDMiIGQ9Ik0xMjUwLDE1NWM2LDAsMTIsMiwxNyw1bDkxOCw1MzBjMTEsNiwxNywxOCwxNywzMHYxMDYwYzAsMTItNywyNC0xNywzMGwtOTE4LDUzMGMtNSwzLTExLDUtMTcsNSAgICAgICBzLTEyLTItMTctNWwtOTE4LTUzMGMtMTEtNi0xNy0xOC0xNy0zMFY3MTljMC0xMiw3LTI0LDE3LTMwbDkxOC01MzBjNS0zLDExLTUsMTctNWwwLDBWMTU1eiBNMTI1MCwwYy0zMywwLTY1LDgtOTUsMjVMMjM3LDU1NSAgICAgICBjLTU5LDM0LTk1LDk2LTk1LDE2NHYxMDYwYzAsNjgsMzYsMTMwLDk1LDE2NGw5MTgsNTMwYzI5LDE3LDYyLDI1LDk1LDI1czY1LTgsOTUtMjVsOTE4LTUzMGM1OS0zNCw5NS05Niw5NS0xNjRWNzE5ICAgICAgIGMwLTY4LTM2LTEzMC05NS0xNjRMMTM0NCwyNWMtMjktMTctNjItMjUtOTUtMjVsMCwwSDEyNTB6Ij48L3BhdGg+CgkJCQk8L2c+CgkJCQk8cG9seWdvbiBjbGFzcz0ic3QxIiBwb2ludHM9IjY0MiwyMTc5IDcyNywxOTQ3IDg5NywyMDg4IDczOCwyMjM0ICAgICAiPjwvcG9seWdvbj4KCQkJCTxnPgoJCQkJCTxwYXRoIGNsYXNzPSJzdDQiIGQ9Ik0xMTcyLDY0NEg5MzljLTE3LDAtMzMsMTEtMzksMjdMNDAxLDIwMzlsMjQxLDEzOWw1NTAtMTUwN2M1LTE0LTUtMjgtMTktMjhMMTE3Miw2NDR6Ij48L3BhdGg+CgkJCQkJPHBhdGggY2xhc3M9InN0NCIgZD0iTTE1ODAsNjQ0aC0yMzNjLTE3LDAtMzMsMTEtMzksMjdMNzM4LDIyMzNsMjQxLDEzOWw2MjAtMTcwMWM1LTE0LTUtMjgtMTktMjhWNjQ0eiI+PC9wYXRoPgoJCQkJPC9nPgoJCQk8L2c+CgkJPC9nPgoJPC9nPgo8L2c+Cjwvc3ZnPgo=" /></a> |  `0x176D1b6c3Fc1db5f7f967Fdc735f8267cCe741F3`
  <a href="javascript:void(0)" style="cursor: default;" alt="Donate via Avalanche"><img src="https://img.shields.io/badge/-Avalanche-4B2224?logo=data:image/svg%2bxml;base64,PHN2ZyB3aWR0aD0iMTUwMyIgaGVpZ2h0PSIxNTA0IiB2aWV3Qm94PSIwIDAgMTUwMyAxNTA0IiBmaWxsPSJub25lIiB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciPgo8cmVjdCB4PSIyODciIHk9IjI1OCIgd2lkdGg9IjkyOCIgaGVpZ2h0PSI4NDQiIGZpbGw9IndoaXRlIi8+CjxwYXRoIGZpbGwtcnVsZT0iZXZlbm9kZCIgY2xpcC1ydWxlPSJldmVub2RkIiBkPSJNMTUwMi41IDc1MkMxNTAyLjUgMTE2Ni43NyAxMTY2LjI3IDE1MDMgNzUxLjUgMTUwM0MzMzYuNzM0IDE1MDMgMC41IDExNjYuNzcgMC41IDc1MkMwLjUgMzM3LjIzNCAzMzYuNzM0IDEgNzUxLjUgMUMxMTY2LjI3IDEgMTUwMi41IDMzNy4yMzQgMTUwMi41IDc1MlpNNTM4LjY4OCAxMDUwLjg2SDM5Mi45NEMzNjIuMzE0IDEwNTAuODYgMzQ3LjE4NiAxMDUwLjg2IDMzNy45NjIgMTA0NC45NkMzMjcuOTk5IDEwMzguNSAzMjEuOTExIDEwMjcuOCAzMjEuMTczIDEwMTUuOTlDMzIwLjYxOSAxMDA1LjExIDMyOC4xODQgOTkxLjgyMiAzNDMuMzEyIDk2NS4yNTVMNzAzLjE4MiAzMzAuOTM1QzcxOC40OTUgMzAzLjk5OSA3MjYuMjQzIDI5MC41MzEgNzM2LjAyMSAyODUuNTVDNzQ2LjUzNyAyODAuMiA3NTkuMDgzIDI4MC4yIDc2OS41OTkgMjg1LjU1Qzc3OS4zNzcgMjkwLjUzMSA3ODcuMTI2IDMwMy45OTkgODAyLjQzOCAzMzAuOTM1TDg3Ni40MiA0NjAuMDc5TDg3Ni43OTcgNDYwLjczOEM4OTMuMzM2IDQ4OS42MzUgOTAxLjcyMyA1MDQuMjg5IDkwNS4zODUgNTE5LjY2OUM5MDkuNDQzIDUzNi40NTggOTA5LjQ0MyA1NTQuMTY5IDkwNS4zODUgNTcwLjk1OEM5MDEuNjk1IDU4Ni40NTUgODkzLjM5MyA2MDEuMjE1IDg3Ni42MDQgNjMwLjU0OUw2ODcuNTczIDk2NC43MDJMNjg3LjA4NCA5NjUuNTU4QzY3MC40MzYgOTk0LjY5MyA2NjEuOTk5IDEwMDkuNDYgNjUwLjMwNiAxMDIwLjZDNjM3LjU3NiAxMDMyLjc4IDYyMi4yNjMgMTA0MS42MyA2MDUuNDc0IDEwNDYuNjJDNTkwLjE2MSAxMDUwLjg2IDU3My4wMDQgMTA1MC44NiA1MzguNjg4IDEwNTAuODZaTTkwNi43NSAxMDUwLjg2SDExMTUuNTlDMTE0Ni40IDEwNTAuODYgMTE2MS45IDEwNTAuODYgMTE3MS4xMyAxMDQ0Ljc4QzExODEuMDkgMTAzOC4zMiAxMTg3LjM2IDEwMjcuNDMgMTE4Ny45MiAxMDE1LjYzQzExODguNDUgMTAwNS4xIDExODEuMDUgOTkyLjMzIDExNjYuNTUgOTY3LjMwN0MxMTY2LjA1IDk2Ni40NTUgMTE2NS41NSA5NjUuNTg4IDExNjUuMDQgOTY0LjcwNkwxMDYwLjQzIDc4NS43NUwxMDU5LjI0IDc4My43MzVDMTA0NC41NCA3NTguODc3IDEwMzcuMTIgNzQ2LjMyNCAxMDI3LjU5IDc0MS40NzJDMTAxNy4wOCA3MzYuMTIxIDEwMDQuNzEgNzM2LjEyMSA5OTQuMTk5IDc0MS40NzJDOTg0LjYwNSA3NDYuNDUzIDk3Ni44NTcgNzU5LjU1MiA5NjEuNTQ0IDc4NS45MzRMODU3LjMwNiA5NjQuODkxTDg1Ni45NDkgOTY1LjUwN0M4NDEuNjkgOTkxLjg0NyA4MzQuMDY0IDEwMDUuMDEgODM0LjYxNCAxMDE1LjgxQzgzNS4zNTIgMTAyNy42MiA4NDEuNDQgMTAzOC41IDg1MS40MDIgMTA0NC45NkM4NjAuNDQzIDEwNTAuODYgODc1Ljk0IDEwNTAuODYgOTA2Ljc1IDEwNTAuODZaIiBmaWxsPSIjRTg0MTQyIi8+Cjwvc3ZnPgo=" /></a> |  `0x176D1b6c3Fc1db5f7f967Fdc735f8267cCe741F3`
</details>

