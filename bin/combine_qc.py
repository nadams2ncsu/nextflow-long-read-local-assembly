#!/usr/bin/env python3

import argparse
import pandas as pd

#########################################################
#  Combine coverage qc results into 1 summary text file
#  
# - input: sample coverage qc files
# - output: summary of all sample coverage qc files


#######################
# combine files 


## required arguments
parser = argparse.ArgumentParser()
parser.add_argument("--input", nargs="+", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

## store all values in a dataframe
dfs = [
    pd.read_csv(f, sep="\t", keep_default_na=False)
    for f in args.input
]

## save summary text file
pd.concat(dfs, ignore_index=True).to_csv(
    args.output,
    sep="\t",
    index=False
)
