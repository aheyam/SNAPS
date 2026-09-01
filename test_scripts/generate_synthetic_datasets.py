"""
Script to generate synthetic observed and predicted shifts for testing SNAPS
"""

import numpy as np
import pandas as pd
from scipy.stats import norm, multivariate_normal
from random import choices
# from Bio import SeqIO, Align
import logging
import yaml
import sys
from pathlib import Path

import argparse

#### Import arguments
parser = argparse.ArgumentParser(
        description="Generate synthetic ovserved and predicted shifts")

# Mandatory arguments
parser.add_argument("output_dir",
                    help="The directory the output of the script will be stored in")
parser.add_argument("config_file", 
                    help="Where to read information such as atom mean shifts, expected error etc.")

# Optional arguments
parser.add_argument("-N", default=1, type=int, help="The number of synthetic datasets to generate")
parser.add_argument("--seq_length", default=200, type=int, help="The number of amino acids in the sequence")
parser.add_argument("--seed", default=0, type=int, help="Random seed")

# Parameters controlling missing data

# Parameters controlling wrong predictions

args = parser.parse_args(sys.argv[1:])

#### Import config file
f = open(args.config_file, 'r')

pars = yaml.safe_load(f)

# Get the full paths of extra parameter files
par_dir = Path(args.config_file).parent
for par_name, par in pars.items():
    if par_name.endswith('_file'):
        pars[par_name] = str(par_dir / par)

# Load the mean and standard deviations of the observed shifts
obs_mean = pd.read_csv(pars["generic_shift_file"], index_col=0)
obs_stdev = pd.read_csv(pars["generic_shift_stdev_file"], index_col=0)

# Calculate stdev of predicted shifts
pred_stdev = obs_stdev.copy()
for atom in pars["atom_list"]:
    pred_stdev.loc[:,atom] = np.sqrt(obs_stdev.loc[:,atom]**2 - pars["atom_prediction_error"][atom]**2)
    pred_stdev[pred_stdev<0] = 0    # Set stdev to 0 if the calculated stdev is -ve

#### Generate the synthetic data
rng = np.random.default_rng(seed=args.seed)   # Set up numpy random number generator

synthetic_datasets = []
for x in range(args.N):
    # Generate the sequence
    # Temporary bad sequence generator
    sequence_list = [rng.choice(list("ACDEFGHIKLMNPQRSTVWY")) for i in range(args.seq_length)]
    # TODO: use [rng.choice(aa_list, p=list_of_aa_probabilities) for i in range(args.seq_length)]
    # this will allow you to have probabilities for each amino acid type.
    # Need to calculate the probabilities from the testset proteins
    # Possibly also extendable to probabilities conditional on previous amino acid,
    # but may require a loop and selecting different columns of the transition matrix

    df = pd.DataFrame({"Res_N":range(1,args.seq_length+1), "Res_type":sequence_list})
    df["Res_name"] = df.Res_N.astype(str) + df.Res_type
    df.Res_name = df.Res_name.str.rjust(5)
    df.index = df.Res_name
    df.index.name = None

    # Generate the predicted shifts
    # (this has to be done first to get a distribution of errors to match ShiftX2)
    for atom in pars["atom_list"]:
        df[atom+"_pred"] = np.nan
        for res in df.Res_type.unique():
            mask = df.Res_type == res
            df.loc[mask, atom+"_pred"] = rng.normal(loc=obs_mean.loc[res, atom], scale=pred_stdev.loc[res, atom], size=mask.sum())
        df.loc[:, atom+"_pred"] = df.loc[:, atom+"_pred"].round(decimals=2)

    # Generate the observed shifts
    for atom in pars["atom_list"]:
        df[atom] = df[atom+"_pred"] + rng.normal(loc=0, scale=pars["atom_prediction_error"][atom], size=len(df.index))
        
        df.loc[:, atom] = df.loc[:, atom].round(decimals=2)
    
    # Get rid of proline H and glycine CB
    if "H" in pars["atom_list"]:
        df.loc[df.Res_type=="P", "H"] = np.nan
        df.loc[df.Res_type=="P", "H_pred"] = np.nan
    if "CB" in pars["atom_list"]:
        df.loc[df.Res_type=="G", "CB"] = np.nan
        df.loc[df.Res_type=="G", "CB_pred"] = np.nan

    synthetic_datasets += [df.copy()]

#### Output the results
output_dir = Path(args.output_dir)
output_dir.mkdir(parents=True, exist_ok=True)

for i in range(args.N):
    target_dir = output_dir/str(i+1)
    target_dir.mkdir(parents=True, exist_ok=True)

    # Print sequence
    sequence = "".join(list(synthetic_datasets[i].Res_type))
    with open(target_dir/"sequence.txt", 'w') as fp:
        print(sequence, file=fp)
    
    # Print observed data
    columns = ["Res_name"]+pars["atom_list"]
    obs_df = synthetic_datasets[i].loc[:, columns]
    obs_df.to_csv(target_dir/"observed_shifts.csv", float_format="%.2f", index=False)

    # Print predictions
    columns = ["Res_N", "Res_type"] + [atom+"_pred" for atom in pars["atom_list"]]
    preds_df = synthetic_datasets[i].loc[:, columns]
    preds_df.to_csv(target_dir/"predicted_shifts.csv", float_format="%.2f", index=False)

    



