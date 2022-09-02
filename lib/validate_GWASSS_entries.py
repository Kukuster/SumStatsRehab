# standard library
import sys
import re
from typing import Any, Dict, List, Literal, Tuple, Union
import os
import json
import subprocess
import time
import gzip
import io

# third-party libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import magic
from tqdm import tqdm

# local
from lib.utils import run_bash
from lib.report_utils import write_report_to_dir
from lib.standard_column_order import STANDARD_COLUMN_ORDER



def file_exists(path: str):
    return os.path.isfile(path)
def dir_exists(path: str):
    return os.path.isdir(path)

def perc(x, total):
    if x/total < 0.0001:
        return "<0.01%"
    else:
        return str(round((x/total)*100, 2)) + "%"



def validate_GWASSS_entries(
    GWAS_FILE: str,
    FORMAT_OR_CONFIG_FILE: str = "standard",
    REPORT_DIR: Union[str, None] = None,
    TICK_LABELS: List[str] = ["0", "1e-8", "1e-5", "1e-3", ".03", ".3", "1"],
    TICKS_WIDTH_RULE: Literal['even', 'log10'] = 'log10',
):
    """
    Loops through the GWAS summary stats file and analyses which data points are missing or invalid.
    Then, creates a report, showing the proportions of missing or invalid data points with respect to the statistical significance (p-value)

    Parameters
    ----------
    GWAS_FILE : str
        GWAS summary statistics file in tsv or tsv.gz format

    JSON_CONFIG : str | "standard"
        Either:
         - path to json config file specifying which columns are which
         - string "standard", which means the input file is in the internal "standard" format (e.g. formatted with fix or sort commands)

    REPORT_DIR : str | None
        If set, is a dir name where textual and graphical reports will be saved.
        Otherwise, doesn't save any reports into the fs but only pops up a graphical report in matplotlib plots

    TICK_LABELS : List[str]
        p-value interval tick points for the stacked histogram plot.

        List of strings with numerical values between 0 and 1 in ascending order. Sets the tick labels for the stacked histogram plot.
        Actual ticks' numerical values will be inferred from the strings.

    TICKS_WIDTH_RULE : 'even' | 'log10'
        Sets the rule for calculation of width of bins in the stacked histogram.

        With 'even', all bars will have unit width, regardless of the numerical difference between ticks.

        With 'log10', all bars will have widths adjusted in accord to log10 scale.
        If the very first tick is zero, then the bin size from zero to the next tick equals to distance between 1e-4 and 1.
    """

    if REPORT_DIR is None:
        REPORT_ABS_DIR = None
    else:
        REPORT_ABS_DIR = os.path.abspath(REPORT_DIR)

    if not file_exists(GWAS_FILE):
        raise ValueError(f"passed GWAS SS file doesn't exist at path {GWAS_FILE}")

    mime: str = magic.from_file(GWAS_FILE, mime=True)
    if mime == 'application/gzip' or mime == 'application/x-gzip':
        GWAS_FILE_o_gz: io.RawIOBase = gzip.open(GWAS_FILE, 'r')  # type: ignore # GzipFile and RawIOBase _are_ in fact compatible
        GWAS_FILE_o = io.TextIOWrapper(io.BufferedReader(GWAS_FILE_o_gz))
    elif mime == 'text/plain':
        GWAS_FILE_o = open(GWAS_FILE, 'r')
    elif mime == 'inode/x-empty':
        raise FileNotFoundError('The provided file is empty!')
    else:
        raise ValueError(f"Got unexpected type of file: {mime}")


    if FORMAT_OR_CONFIG_FILE == "standard":
        cols_i: Dict[str, int] = {STANDARD_COLUMN_ORDER[i]: i for i in range(len(STANDARD_COLUMN_ORDER))}
    else:
        if not file_exists(FORMAT_OR_CONFIG_FILE):
            raise ValueError(f"passed GWAS SS file doesn't have the corresponding json config file at path: {FORMAT_OR_CONFIG_FILE}. Please create one based on the config.example.json file")
        cols_i: Dict[str, int] = json.load(open(FORMAT_OR_CONFIG_FILE,))


    separator = '\t'

    x = TICK_LABELS

    # p-value interval points:
    ticks = [float(x_label) for x_label in x]

    if not (  ticks == sorted(ticks) and len(ticks) == len(set(ticks))   ):
        raise ValueError("Ticks have to be in strictly ascending order")
    
    if not (  ticks[0] >= 0 and ticks[-1] <= 1   ):
        raise ValueError("Ticks have to be in range from 0 to 1")



    # # # # # # # # # # # # # # # # # # # # # # # # # #
    #                                                 #
    #                    CONSTANTS                    #
    #                                                 #
    # # # # # # # # # # # # # # # # # # # # # # # # # #

    GOOD_ENTRY = 0
    MISSING_P_VALUE = 1
    INVALID_ENTRY = 2


    # indices for boolean values in a list of issues for each SNP
    INVALID_ROW = 0
    INVALID_RSID = 1
    INVALID_CHR = 2
    INVALID_BP = 3
    INVALID_EA = 4
    INVALID_OA = 5
    INVALID_EAF = 6
    INVALID_SE = 7
    INVALID_ES = 8
    ISSUES=[
        INVALID_ROW,
        INVALID_RSID,
        INVALID_CHR,
        INVALID_BP,
        INVALID_EA,
        INVALID_OA,
        INVALID_EAF,
        INVALID_SE,
        INVALID_ES,
    ]
    ISSUES_LABELS = [
        "format",
        "rsID",
        "Chr",
        "BP",
        "EA",
        "OA",
        "EAF",
        "SE",
        "beta",
    ]
    ISSUES_COLORS=[
        "#ff0000", # format
        "#777ae5", # rsID
        "#cf44a1", # Chr
        "#ff4481", # BP
        "#ffa121", # EA
        "#ff9191", # OA
        "#fdbc64", # EAF
        "#563E3E", # std. err.
        "#175a63", # beta

    ]

    NUCLEOTIDES = ['a', 't', 'c', 'g']
    NO_NUCLEOTIDE = '-'

    ALLOW_MULTI_NUCLEOTIDE_POLYMORPHISMS = True

    CATEGORY_CHR = [
    '1', '01', '2', '02', '3', '03', '4', '04', '5', '05', '6', '06', '7', '07', '8', '08', '9', '09',
    '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
    '21', '22', '23', 'X', 'x', 'Y', 'y', 'M', 'm']
    # CATEGORY_CHR = [
    # '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
    # '21', '22', '23', 'X', 'x', 'Y', 'y', 'M', 'm']



    # # # # # # # # # # # # # # # # # # # # # # # # # #
    #                                                 #
    #                    FUNCTIONS                    #
    #                                                 #
    # # # # # # # # # # # # # # # # # # # # # # # # # #

    # https://stackoverflow.com/a/850962/6041933  # a comment to the answer
    # https://gist.github.com/zed/0ac760859e614cd03652
    def wccount(filename: str):
        """counts the number of lines in the file"""
        if mime == 'application/gzip' or mime == 'application/x-gzip':
            out = run_bash(f'gunzip -c "{filename}" | wc -l')
        else:
            out = run_bash(f'wc -l "{filename}"')
        return int(out.partition(' ')[0])

    def is_null(val: str) -> bool:
        return val.lower() in ["", " ", ".", "-", "na", "nan"]




    # # # # # # # # # # # # # # # # # # # # # # # # # #
    #                                                 #
    #       FUNCTION THAT CHECKS EACH SNP ENTRY       #
    #                                                 #
    # # # # # # # # # # # # # # # # # # # # # # # # # #

    def check_row(line_cols: List[str]) -> Union[
            # "good" entry
            Tuple[float, Literal[0], List[bool]],
            # "missing/invalid p-value" entry
            Tuple[None,  Literal[1], List[bool]],
            # "invalid" entry, having some issues (listed in the list)
            Tuple[float, Literal[2], List[bool]],
        ]:
        """
        This function runs for EVERY LINE of the input file,
        which may be MILLIONS OF TIMES per script execution

        Returns:
            - a p-value of a SNP,
            - a report of whether a SNP entry is valid,
            - and a list of issues with it
        """

        issues = [False] * len(ISSUES)
        missing_pvalue = False

        ### First check if p-value itself is present ###
        pval = None
        try:
            pval = line_cols[cols_i["pval"]]
            if is_null(pval) or not (0 <= float(pval) <= 1):
                missing_pvalue = True
            pval = float(pval)
        except:
            pval = None
            missing_pvalue = True


        ### Try getting all columns. If some not present, will throw ###
        try:
            rsid  = line_cols[cols_i["rsID"]]
            chrom = line_cols[cols_i["Chr"]]
            bp    = line_cols[cols_i["BP"]]
            ea    = line_cols[cols_i["EA"]]
            oa    = line_cols[cols_i["OA"]]
            af    = line_cols[cols_i["EAF"]]
            se    = line_cols[cols_i["SE"]]
            es    = line_cols[cols_i["beta"]]
            # n     = line_cols[cols_i["N"]]

        except:
            issues[INVALID_ROW] = True
            if missing_pvalue:
                return None, MISSING_P_VALUE, issues
            else:
                assert pval is not None
                return pval, INVALID_ENTRY, issues # pval is None

        ### Check any reasons this SNP will be discarded later ###

        # 1. rsID
        try:
            if not re.match("^rs\d+$", rsid):
                issues[INVALID_RSID] = True
        except:
            issues[INVALID_RSID] = True

        # 2. chromosome
        try:
            if chrom not in CATEGORY_CHR and chrom[3:] not in CATEGORY_CHR:
                issues[INVALID_CHR] = True
        except:
            issues[INVALID_CHR] = True

        # 3. base pair position
        try:
            bp = int(float(bp)) # using float allows sci notation string
            if bp < 0:
                issues[INVALID_BP] = True
        except:
            issues[INVALID_BP] = True

        # 4. effect allele
        try:
            if ea == '':
                issues[INVALID_EA] = True
            elif ea == NO_NUCLEOTIDE:
                issues[INVALID_EA] = False
            elif ALLOW_MULTI_NUCLEOTIDE_POLYMORPHISMS:
                for char in ea.lower():
                    if char not in NUCLEOTIDES:
                        issues[INVALID_EA] = True
            else:
                if ea.lower() not in NUCLEOTIDES:
                    issues[INVALID_EA] = True
        except:
            issues[INVALID_EA] = True

        # 5. other allele
        try:
            if oa == '':
                issues[INVALID_OA] = True
            elif oa == NO_NUCLEOTIDE:
                issues[INVALID_OA] = False
            elif ALLOW_MULTI_NUCLEOTIDE_POLYMORPHISMS:
                for char in oa.lower():
                    if char not in NUCLEOTIDES:
                        issues[INVALID_OA] = True
            else:
                if oa.lower() not in NUCLEOTIDES:
                    issues[INVALID_OA] = True
        except:
            issues[INVALID_OA] = True

        # 6. effect allele frequency or minor allele frequency
        try:
            if not (0 <= float(af) <= 1):
                issues[INVALID_EAF] = True
        except:
            issues[INVALID_EAF] = True

        # 7. standard error
        try:
            float(se) # will throw if not float
            if is_null(se):
                issues[INVALID_SE] = True
        except:
            issues[INVALID_SE] = True

        # 8. effect size (odds ratio or beta-value)
        try:
            float(es) # will throw if not float
            if is_null(es):
                issues[INVALID_ES] = True
        except:
            issues[INVALID_ES] = True

        # # 9. n - sample size
        # #sometimes sample size is fractional
        # if null_entry(n) or not (0 < float(n)):
        #     return INVALID_ENTRY, pval

        if missing_pvalue:
            return None, MISSING_P_VALUE, issues
        else:
            assert pval is not None
            if any(issues):
                return pval, INVALID_ENTRY, issues
            else:
                # all good?
                return pval, GOOD_ENTRY, issues




    # # # # # # # # # # # # # # # # # # # # # # # # # #
    #                                                 #
    #                      MAIN                       #
    #                                                 #
    # # # # # # # # # # # # # # # # # # # # # # # # # #

    MAIN_start_time = STEP1_start_time = time.time()


    #
    # STEP #1
    #    read the file line by line and check the validity of each SNP entry,
    #    and save the report and the p-value if present
    #

    num_of_lines = wccount(GWAS_FILE)
    num_of_snps = num_of_lines - 1
    print(f"number of lines in the file: {num_of_lines}")

    line_i=0
    # skip the first line that is the header
    GWAS_FILE_o.readline()
    line_i+=1

    SNPs_pval = np.zeros(num_of_snps, dtype=np.float64)
    SNPs_report = np.zeros(num_of_snps, dtype=np.int8)
    SNPs_issues = np.zeros((num_of_snps, len(ISSUES)), dtype=np.bool_)

    ### populate the allocated array with report for each SNP as well as its p-value ###
    pbar = tqdm(total=num_of_snps, desc="validating entries ")
    try:
        snp_i = 0
        while True:
            SNPs_pval[snp_i], SNPs_report[snp_i], SNPs_issues[snp_i] = check_row(GWAS_FILE_o.readline().replace('\n','').split(separator))
            snp_i += 1
            pbar.update(1)

    except Exception as e:
        if isinstance(e, IndexError) or isinstance(e, EOFError):
            # it reached the end of the file
            pass
        else:
            print(f'An error occured on line {line_i} of the GWAS SS file (see below)')
            raise e
    pbar.close()
    ### ###


    GWAS_FILE_o.close()
    # print("--- STEP1: %s seconds ---" % (time.time() - STEP1_start_time))

    # result: SNPs_report, SNPs_pval


    #
    # STEP #2
    #    sum up the reports and calculate the parameters before plotting
    #
    STEP2_start_time = time.time()


    # 2.1
    """
    Calculate bars widths and coordinates.
    
    Bar widths start from the right, i.e. the last bar is the range between the last two ticks.
    The very first bar is "no p-value" bar, it has unit width and goes to the left of the first tick.

    If the very first tick is zero, then the bin size from zero to the next tick is 2 units.

    User may have choosen the rule of ticks width to either 'even' or 'log10':
    • With 'even', all bars will have unit width, regardless of the numerical difference between ticks
    • With 'log10', all bars will have widths adjusted in accord to log10 scale,
    with a constant unit_width defined here
    """
    bars_widths = [1.]

    if TICKS_WIDTH_RULE == 'even':
        for i in range(1, len(ticks)):
            bars_widths.append(1.)

    elif TICKS_WIDTH_RULE == 'log10':
        unit_width = np.log10(1) - np.log10(1e-2)
        for i in range(1, len(ticks)):
            if ticks[i-1] == 0:
                bars_widths.append(2.)
            else:
                bars_widths.append(
                    (np.log10(ticks[i]) - np.log10(ticks[i-1]))
                                    / unit_width
                )

    else:
        raise ValueError(f'unknown ticks_width_rule: {TICKS_WIDTH_RULE}')

    negative_bars_widths = list(-np.array(bars_widths))


    # 2.2
    """
    Ticks location (ticks_loc) equals to cumulative of bars_widths,
    shifted by the very first bar width to the left, so the tick #1 equals 0
    """
    ticks_loc: List[float] = []
    cumulative = -bars_widths[0] # starting such that the first loc is 0
    for width in bars_widths:
        cumulative += width
        ticks_loc.append(cumulative)
    del cumulative


    assert len(ticks) == len(bars_widths) == len(ticks_loc), "lists: `ticks`, `ticks_loc`, and `bars_widths` should have the same lengths"


    # 2.3
    ### Counting how many entries are valid, don't have p-value, or invalid for other reasons ###
    missing_pval_bins = [0]*len(ticks)
    good_entry_bins = [0]*len(ticks)
    invalid_entry_bins = [0]*len(ticks)
    # besides total, for each of the bin we'll store the number of invalid entries for each type
    invalid_entry_bins_reason_bins = np.zeros((len(ticks), max(ISSUES)+1)).astype(int)

    pbar = tqdm(total=len(SNPs_report), desc="calculating reports")
    for line_i in range(len(SNPs_report)):

        if SNPs_report[line_i] == MISSING_P_VALUE:
            missing_pval_bins[0] += 1

        elif SNPs_report[line_i] == GOOD_ENTRY:
            for j in range(1,len(ticks)):
                if SNPs_pval[line_i] <= ticks[j]:
                    good_entry_bins[j] += 1
                    break

        elif SNPs_report[line_i] == INVALID_ENTRY:
            for j in range(1, len(ticks)):
                if SNPs_pval[line_i] <= ticks[j]:
                    invalid_entry_bins[j] += 1
                    invalid_entry_bins_reason_bins[j] += SNPs_issues[line_i]
                    break

        pbar.update(1)
    pbar.close()

    ### ###

    # print("--- STEP2: %s seconds ---" % (time.time() - STEP2_start_time))
    # print("=== MAIN: %s seconds ===" % (time.time() - MAIN_start_time)) # plotting doesn't count
    print("generating reports")



    #
    # STEP #3
    #     save csv file with report for each issue
    #

    if REPORT_ABS_DIR:
        if not dir_exists(REPORT_ABS_DIR):
            os.makedirs(REPORT_ABS_DIR)


    issues_count_arr = np.sum(SNPs_issues, axis=0)
    issues_count: Dict[str, int] = {}

    for issue_i in range(0, len(ISSUES)):
        issues_count[ISSUES_LABELS[issue_i]] = issues_count_arr[issue_i]

    issues_count["pval"] = sum(missing_pval_bins)

    if any(issues_count.values()):
        print("found issues:")
        for issue, count in issues_count.items():
            if count:
                print(f"    {issue}: {count}/{num_of_snps} ({perc(count,num_of_snps)})")

    issues_count["total_entries"] = num_of_snps


    if REPORT_ABS_DIR: 
        write_report_to_dir(issues_count, REPORT_ABS_DIR)



    #
    # STEP #4
    #    plot
    #

    ### CALC: proportion of invalid entries in total ###

    invalid_entries_totally = sum(invalid_entry_bins) + sum(missing_pval_bins)
    proportion_of_invalid_entries_totally = invalid_entries_totally / num_of_snps
    percentage_of_invalid_entries_totally = proportion_of_invalid_entries_totally * 100
    percentage_of_invalid_entries_totally_str = str(np.round(percentage_of_invalid_entries_totally, 1)) + "%"


    ### PLOT: the figure, labels, ticks, bars ###

    fig, ax = plt.subplots(num="valid/invalid SNPs")

    image_name = GWAS_FILE.split('/')[-1]
    fig.canvas.set_window_title(image_name) # sets the window title to the filename

    ax.set_title(
        f"invalid SNPs: {invalid_entries_totally}/{num_of_snps} ({percentage_of_invalid_entries_totally_str})")

    # # Hide the right and top spines
    # ax.spines['right'].set_visible(False)
    # ax.spines['top'].set_visible(False)

    ax.tick_params(axis='x', labelsize=9)

    ax.set_xticks(ticks_loc)
    ax.set_xticklabels(x)
    ax.bar(ticks_loc, missing_pval_bins, negative_bars_widths, color='#7f7f7f', align='edge')
    ax.bar(ticks_loc, good_entry_bins, negative_bars_widths, color='#0000ff', align='edge', label="valid SNPs")
    ax.bar(ticks_loc, invalid_entry_bins, negative_bars_widths, color='#ff0000', align='edge', bottom=good_entry_bins, label="invalid SNPs")
    ax.set_xlabel("p-value", fontweight="heavy", fontsize=14)
    ax.set_ylabel("N of SNPs", fontweight="heavy", fontsize=14)
    ax.set_xlim([ticks_loc[0]+negative_bars_widths[0], ticks_loc[-1]])

    max_bar_height = max(
        np.array(missing_pval_bins) + np.array(good_entry_bins) + np.array(invalid_entry_bins) # np arrays add element-wise
    )

    plt_bottom, plt_top = ax.set_ylim(0, max_bar_height*1.15 if max_bar_height else 1)
    plt_height = plt_top - plt_bottom


    ### CALC: points right at the middle of each bin ###
    bins_mid_points = list(np.array(ticks_loc) - np.array(bars_widths)/2)


    ### PLOT: caption for "no p-value" bin ###
    ax.text(x=bins_mid_points[0], y=plt_top*-0.08, s="no\np-value",
        horizontalalignment='center',
    )


    ### CALC: total entries, proportion and percentage of invalid entries  ###

    total_p_value_entries_bins = [good_entry_bins[i]+invalid_entry_bins[i] for i in range(len(good_entry_bins))]

    proportion_of_invalid_entries_bins = [0.] + [
        invalid_entry_bins[i]/total_p_value_entries_bins[i] if total_p_value_entries_bins[i] else 0.
                    for i in range(1, len(good_entry_bins))]

    percentage_of_invalid_entries_bins = np.round(np.array(proportion_of_invalid_entries_bins)*100).astype(int)
    percentage_of_invalid_entries_bins_str = np.char.array(percentage_of_invalid_entries_bins) + "%" # type: ignore # pylance mistakenly doesn't recognize np.char



    ### PLOT: representation of the percentage of invalid entries for each bin ###
    # the bottom and the top spines of the plot represent 0% and 100%
    # skipping the first, "no p-value" bin

    ax.plot(
        # X: points right at the mid of bins (except the no p-value bin)
        bins_mid_points[1:],

        # Y: how much proportion -> that much height within the plot
        np.array(proportion_of_invalid_entries_bins[1:]) * plt_height + plt_bottom,

        linestyle='-',
        color="red",
        alpha=0.5,
        linewidth=2,
    )


    ### PLOT: Captions above the points of the plot ###
    # shows percentage of invalid entries for each bin in text
    # skipping the first, "no p-value" bin

    # points right at the mid of bins
    X = bins_mid_points
    # how much proportion -> that much height within the plot, also lifted 5%
    Y = (np.array(proportion_of_invalid_entries_bins) * plt_height) + plt_height*0.05 + plt_bottom

    for i in range(1, len(percentage_of_invalid_entries_bins_str)):
        if proportion_of_invalid_entries_bins[i] > 0.15:
            text = ax.text(X[i], Y[i], s=percentage_of_invalid_entries_bins_str[i],
                horizontalalignment='center',
                color="#bf3f3f", # caption may overlap with the red stuff from stacked bar
                fontsize=10,
                fontweight="demibold",
            )
        else:
            text = ax.text(X[i], Y[i], s=percentage_of_invalid_entries_bins_str[i],
                horizontalalignment='center',
                color="#ff0000",
                fontsize=10,
            )
        # text.set_path_effects([path_effects.Stroke(linewidth=1, foreground='#000000'),
        #                        path_effects.Normal()])


    if REPORT_ABS_DIR: fig.savefig(os.path.join(REPORT_ABS_DIR, image_name+'.png'))


    ### PLOT: bar chart for issues in each of the p-value bins ###
    for i in range(1, len(invalid_entry_bins)):
        image_name = f'bin_{i}__{x[i-1]}—{x[i]}'
        plot_i, ax = plt.subplots(num=image_name)

        issues_proportions = [0] * len(ISSUES)
        total_invalids = invalid_entry_bins[i]
        total_snps = invalid_entry_bins[i] + good_entry_bins[i]

        proportion_of_invalids = total_invalids / total_snps if total_snps else 0
        percentage_of_invalids = proportion_of_invalids * 100
        percentage_of_invalids_str = str(np.round(percentage_of_invalids, 1)) + "%"

        ax.set_title(f'issues in p-value: {x[i-1]} — {x[i]}\ninvalid SNPs: {total_invalids}/{total_snps} ({percentage_of_invalids_str})')
        ax.set_ylim(0, total_invalids if total_invalids > 0 else 1)
        ax.set_ylabel("N of invalid SNPs", fontweight="demibold", fontsize=14)

        if total_invalids == 0:
            ax.bar(ISSUES_LABELS, height=issues_proportions, width=1, color=ISSUES_COLORS)
        else:
            for issue in range(len(ISSUES)): # for each ISSUE
                issues_proportions[issue] = invalid_entry_bins_reason_bins[i][issue] / total_invalids
            ax.bar(ISSUES_LABELS, height=invalid_entry_bins_reason_bins[i], width=1, color=ISSUES_COLORS)

        if REPORT_ABS_DIR: plot_i.savefig(os.path.join(REPORT_ABS_DIR, image_name+'.png'))



    ### FINALLY: display results and the figure (if report directory was not set) ###

    # print(f'missing_pval_bins = {missing_pval_bins}')
    # print(f'good_entry_bins = {good_entry_bins}')
    # print(f'invalid_entry_bins = {invalid_entry_bins}')

    # print(f'proportion_of_invalid_entries_bins = {proportion_of_invalid_entries_bins}')

    # for i in range(1, len(invalid_entry_bins)):
    #     print(f"{x[i-1]} — {x[i]}: {[invalid_entry_bins_reason_bins[i][issue] for issue in ISSUES]}")



    if not REPORT_DIR:
        plt.show()
        input("")
        input("")



if __name__ == "__main__":
    GWAS_FILE = sys.argv[1]
    FORMAT_OR_CONFIG_FILE = sys.argv[2]
    REPORT_DIR = sys.argv[3]

    validate_GWASSS_entries(GWAS_FILE, FORMAT_OR_CONFIG_FILE, REPORT_DIR)

