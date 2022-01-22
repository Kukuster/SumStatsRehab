# SumStatsRehab

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
pip install git+https://github.com/Kukuster/SumStatsRehab.git@v1.1.1 --upgrade
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

    "build": "grch37"
}
```

Notes:
 - SumStatsRehab will only consider data from the columns which indices are specified in the config file. If one of the above columns is present in the SS file but wasn't specified in the config file, then SumStatsRehab treats the column as missing.
 - In this example, all the 10 columns from the list of supported columns are present. But none of the columns above are mandatory. If certain columns are missing, the `fix` command will attempt to restore them if possible.


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
 - `FREQ_DATABASE_SLUG` is a population slug from a frequency database in dbSNP

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

It is also possible to set `EAF` to a weighted average of multiple colums. E.g. if there are separate freq. columns for case and control groups, and average freq. is needed, number of participants in each group will serve as weights for the two columns:
```
{
    "EAF": {
        "4": 1001,
        "5": 2500
    },
}
```



### "standard" format
When `fix`ing, file is first formatted into this internal format. Output file is also in this format.

 - file is in the tsv format, i.e. tabular tab-separated format (bare, zipped, or gzipped)
 - there's a one-line header in the file on the first line. All other lines are the data entries
 - the file has precisely columns defined as `STANDARD_COLUMN_ORDER` in `lib/standard_column_order.py`.
    - file has exactly these columns, exactly this number of columns, and no other columns
    - columns are in this exact order
    - if the original file was missing a column, an empty column should be taking its place (entries are *the empty string*)


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
4. Allele frequency is restored as the effect allele frequency, and it’s done in a similar way: first by a lookup in DB1 for matching rsID and EA, if rsID is present; second by a lookup into DB2 for matching Chr and BP and EA. Restoring allele frequency requires frequency entries in the INFO field in the SNPs dataset that was preprocessed. By default, frequencies are taken from the database of Genotypes and Phenotypes (dbGaP), but a different population can be specified as an argument to the FIX command.
5. Each of the columns, standard error, beta, and p-value, are restored from the other two, using the following relation: *__s = β/z__*, where *__s__* is the standard error, *__β__* is beta, and *__z__* is z-score that corresponds to the p-value in the two-tailed test [(Zhu et al. 2016)](https://www.nature.com/articles/ng.3538). If standard error is signed, then beta can be restored accurately, but if standard error is provided as an absolute value, the restored beta will be correct only as an absolute value, with an unknown sign. Files restored in this manner will have little utility in downstream applications, but may be useful in comparing relative effect sizes. 

Warning: Restored betas with an unknown sign should not be utilized in any downstream application. 


## BACKLOG
 - (maybe) improve restoring alleles by adding checks for an exact match of flipped alleles if other checks didn't help. This requires having all SNPs for a particular ChrBP in the memory and is relevant only for restoring alleles by looping through the file sorted by Chr and BP.
 - add the ability to specify additional columns from the GWAS SS file that the user wants to include in the end file. This would be an array of integers in the json config file for the input GWAS SS file.
 - **improve code in the main file: `SumStatsRehab.py`**
 - improve resolver architecture in `loop_fix.py`: make a separate function loopDB1 and loopDB2 that will loop through enough entries in a DB before every resolver and rewrite a "global" object with properties to be fields from the DB: rsID, Chr, BP, alleles, EAF. So resolvers for rsID and ChrBP will be similar to ones for alleles and EAF. Resolvers for these fields then should operate on `fields` and that object with fields from a DB. This way a really strong optimization, flexibility, and modularity of resolvers will be achieved. `run_all` doesn't have to have resolvers and resolvers_args object to be passed, it can just use the global ones.
 - improve the interface for liftover. SumStatsRehab fix should work for all sorts of liftovers between builds 36, 37, and 38, including back liftover. If the user omits the preprocessed dbSNP databases as input but specifies the chain file, it can perform liftover only.
 - add support for OR, and, maybe, restoration of OR from beta or vice versa.
 - add a keyword argument that specifies a temp directory for intermediate files. GWAS SS files are usually 1-4 Gigs unpacked.
 - set alleles column to uppercase during preparation (in `prepare_GWASSS_columns.py` script).
 - feature: save a human-readable textual report about the overall results of restoration (e.g. "performed a liftover, n rsIDs restored, n Chrs lost, ...")
 - at the moment of 2021.11.14, the following executables are assumed to be available in PATH: `bash`, `cut`, `paste`, `sort`, `awk`, `gzip`, `gunzip`, `head`, `tail`, `rm`, `wc`. Need to test SumStatsRehab with a different versions of `bash`, `awk` (including `gawk`, `nawk`, `mawk`. E.g. even though `gawk` is default for GNU/Linux, Ubuntu has `mawk` by default).
 - **make SumStatsRehab installable via `pip`**
 - Study what is a better approach to restoring EAF from other dbs. Bc for other populations there are not a lot of snps having frequency data. When you try to restore eaf for a more specific populations, it will miss a lot of snps in the ss files, therefore reducing overall accuracy. Idea for workaround: ability to specify multiple dbs, so the each next one in a list will be a lower priority.
 - `diagnose` command script should also generate a bar chart for the bin with missing p-value
 - **add data to the csv report of issues: rsID, rsID_restorable, Chr, Chr_restorable, ... etc. Use this report in the FIX command when deciding on the workflow.**
 - **add the following logic for restoring MAF: if the DB specifies allele frequency for all alleles except one (dot '.'), then it should calculated as 1 minus frequency of other alleles instead of leaving a dot**
