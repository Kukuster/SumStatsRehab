import sys

import pandas as pd

ss_1 = sys.argv[1]
ss_2 = sys.argv[2]

ss_1_df = pd.read_csv(ss_1, sep='\t')
ss_2_df = pd.read_csv(ss_2, sep='\t')

print("Processing ...")

def compare_equality(col_1, col_2):
    """
    This function compare the equality of two columns
    and output an accuracy percentage
    """

    return float(col_1.eq(col_2.values).mean())

def compare_numeric(col_1, col_2):
    """
    This function compare two numeric values by a distance measure
    d(a, b) = |a - b|
    """

    df['distance'] = col_1 - col_2
    df['distance'].mean()
    return 1 - df['distance'].mean()


cols = {
    # "rsID":compare_equality,
    "CHR":compare_equality,
    "POS":compare_equality,
    "EA":compare_equality,
    "NEA":compare_equality,
    "EAF":compare_numeric,
    "OR":compare_numeric,
    "beta":compare_numeric,
    "SE":compare_numeric,
    "pval":compare_numeric,
    "N":compare_numeric,
    "INFO":compare_numeric,
}


df = pd.merge(ss_1_df, ss_2_df, on="rsID")

df = df[df.rsID != "."]
columns = list(ss_1_df.columns)

columns.remove("rsID")

accuracy = {}
for col in columns:
    accuracy[col] = cols[col](df[f"{col}_x"], df[f"{col}_y"])

# These accuracy metrics are relative ones, you should compare two file
# accuracy tables. And the closer to one in all fields, the more accurate. 
print("accuracy_table")
print(accuracy)
