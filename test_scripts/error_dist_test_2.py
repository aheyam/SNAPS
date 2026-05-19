#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to analyse the error distribution of predicted shifts

Created 12/5/2026

@author: Alex Heyam
"""

import pandas as pd
import numpy as np
from scipy.stats import norm, multivariate_normal, linregress
from Bio.SeqUtils import seq1
import plotnine
from plotnine import *
from plotnine.ggplot import save_as_pdf_pages
from pathlib import Path
import argparse

def import_testset_shifts(filename, remove_Pro=True, 
                          short_aa_names=True, SS_class=None, SS_class_m1=None):
    """ Import observed chemical shifts from testset data
    
    This function is intended for use with test data only, and is unlikely 
    to work well on 'real' data.
    
    filename: The simplified BMRB file containing observed shift info.
    remove_Pro: If True, remove proline residues from output
    short_aa_names: If True, single letter aa codes are used, otherwise 3 
        letter codes are used
    SS_class: Either None or a list of strings, each of which is a list of 
        amino acids (eg. ["VIA","G","S","T","DN","FHYWC","REKPQML"] would 
        give the HADAMAC classes). If not None, a column SS_class will be 
        created which gives the class containing the residue type.
    SS_class_m1: as above, but for the i-1 residue.
    
    """
    #### Import the observed chemical shifts
    obs_long = pd.read_table(filename)
    obs_long = obs_long[["Residue_PDB_seq_code","Residue_label",
                            "Atom_name","Chem_shift_value"]]
    obs_long.columns = ["Res_N","Res_type","Atom_type","Shift"]
    # Convert residue type to single-letter code
    if short_aa_names: 
        obs_long["Res_type"] = obs_long["Res_type"].apply(seq1)
        obs_long["SS_name"] = (obs_long["Res_N"].astype(str) + 
                obs_long["Res_type"])
        obs_long["SS_name"] = obs_long["SS_name"].str.rjust(5)
    else:
        obs_long["SS_name"] = (obs_long["Res_N"].astype(str) + 
                                obs_long["Res_type"])
        obs_long["SS_name"] = obs_long["SS_name"].str.rjust(7)
        obs_long["Res_type"] = obs_long["Res_type"].apply(seq1)
    obs_long = obs_long.reindex(columns=["Res_N","Res_type","SS_name",
                                            "Atom_type","Shift"])
    
    # Convert from long to wide
    obs = obs_long.pivot(index="Res_N", columns="Atom_type", 
                            values="Shift")
    
    # Add the other columns back in
    tmp = obs_long[["Res_N","Res_type","SS_name"]]
    tmp = tmp.drop_duplicates(subset="SS_name")
    tmp.index = tmp["Res_N"]
    obs = pd.concat([tmp, obs], axis=1)
    
    # Make columns for the i-1 observed shifts of C, CA and CB
    obs_m1 = obs[list({"C","CA","CB","Res_type"}.intersection(obs.columns))]
    obs_m1.index = obs_m1.index+1
    obs_m1.columns = obs_m1.columns + "_m1"
    obs = pd.merge(obs, obs_m1, how="left", left_index=True, 
                    right_index=True)
    
    # Restrict to specific atom types
    atom_set = {"H","N","C","CA","CB","C_m1","CA_m1","CB_m1","HA"}
    obs = obs[["Res_N","Res_type","Res_type_m1","SS_name"]+
                list(atom_set.intersection(obs.columns))]
    
    # Add SS_class information
    if SS_class is not None:
        obs["SS_class"]=obs["Res_type"]
        for g in SS_class:
            obs["SS_class"] = obs["SS_class"].str.replace("["+g+"]", g)
    if SS_class_m1 is not None:
        obs["SS_class_m1"]=obs["Res_type_m1"]
        for g in SS_class_m1:
            obs["SS_class_m1"] = obs["SS_class_m1"].str.replace("["+g+"]", g)
    
    obs.index = obs["SS_name"]
    obs.index.name = None
    
    if remove_Pro:
        # Remove prolines, as they wouldn't be observed in a real spectrum
        obs = obs.drop(obs.index[obs["Res_type"].isin(["PRO","P"])]) 
    
    return(obs)

def import_pred_shifts(filename, filetype, offset=0):
    """ Import predicted chemical shifts from a ShiftX2 results file.

    Returns
    A DataFrame containing the predicted shifts, or None if the import failed.

    Parameters
    filename: path to file containing predicted shifts
    filetype: either "shiftx2" or "sparta+"
    offset: an optional integer to add to the residue number.
    """

    #### Import the raw data
    if filetype == "shiftx2":
        preds_long = pd.read_csv(filename)
        if any(preds_long.columns == "CHAIN"):
            if len(preds_long["CHAIN"].unique())>1:
                print(
                        """Chain identifier dropped - if multiple chains are
                        present in the predictions, they will be merged.""")
            preds_long = preds_long.drop("CHAIN", axis=1)
        preds_long = preds_long.reindex(columns=["NUM","RES","ATOMNAME",
                                                    "SHIFT"])
        preds_long.columns = ["Res_N","Res_type","Atom_type","Shift"]
    elif filetype == "sparta+":
        # Work out where the column names and data are
        with open(filename, 'r') as f:
            for num, line in enumerate(f, 1):
                if line.find("VARS")>-1:
                    colnames_line = num
                    colnames = line.split()[1:]
                    break

        preds_long = pd.read_table(filename, sep="\s+", names=colnames,
                                    skiprows=colnames_line+1)
        preds_long = preds_long.reindex(columns=["RESID","RESNAME",
                                                    "ATOMNAME","SHIFT"])
        preds_long.columns = ["Res_N","Res_type","Atom_type","Shift"]

        # Sparta+ uses HN for backbone amide proton - convert to H
        preds_long.loc[preds_long["Atom_type"]=="HN", "Atom_type"] = "H"
    else:
        print("""Invalid predicted shift type: '%s'. Allowed
                            options are 'shiftx2' or 'sparta+'""" % (filetype))
        return(None)

    print("Imported %d predicted chemical shifts from %s"
                        % (len(preds_long.index), filename))

    #### Initial processing and conversion from long to wide
    # Add sequence number offset and create residue names
    preds_long["Res_N"] = preds_long["Res_N"] + offset
    preds_long.insert(1, "Res_name", (preds_long["Res_N"].astype(str) +
                preds_long["Res_type"]))
    # Left pad with spaces to a constant length (helps with sorting)
    preds_long["Res_name"] = preds_long["Res_name"].str.rjust(5)

    # Convert from long to wide format
    preds = preds_long.pivot(index="Res_N", columns="Atom_type",
                                values="Shift")
    preds.index.name = None

    # Add the other residue name and type back in
    tmp = preds_long[["Res_N","Res_type","Res_name"]]
    tmp = tmp.drop_duplicates(subset="Res_name")
    tmp.index = tmp["Res_N"]
    tmp.index.name = None
    preds = pd.concat([tmp, preds], axis=1)

    #### Make consistent with seq_df (and create if it doesn't already exist)
    # TODO: Maybe this should be split off into a separate function?
    # If seq_df is missing, create it based on preds
    seq_df = None

    if seq_df is None:
        seq_df = preds.copy()[["Res_N","Res_type","Res_name"]]

        # If there are any missing residue numbers, create them
        min_N = seq_df["Res_N"].min()
        max_N = seq_df["Res_N"].max()
        missing_residue_numbers = (set(range(min_N, max_N+1)).
                                    difference(seq_df["Res_N"]))
        if len(missing_residue_numbers)>0:
            tmp = pd.DataFrame({"Res_N":list(missing_residue_numbers),
                                "Res_type":"X","Res_name":np.nan})
            tmp["Res_name"] = tmp["Res_N"].astype(str) + tmp["Res_type"]
            tmp["Res_name"] = tmp["Res_name"].str.rjust(5)
            tmp.index = tmp["Res_N"]
            tmp.index.name = None
            #seq_df = seq_df.append(tmp).sort_index()
            seq_df = pd.concat([seq_df, tmp]).sort_index()

    # Add/delete residues from preds so it matches seq_df
    # Only keep predictions that are in seq_df
    tmp = len(preds.index)
    preds = preds[preds["Res_N"].isin(seq_df["Res_N"])]
    tmp2 = tmp - len(preds.index)
    if tmp2>0:
        print(("Predictions for %d residues were discarded "
                            "because they were not present in the imported "
                            "sequence file") % tmp2)

    # Add predictions for any residue number that is in seq_df but not preds
    tmp = pd.DataFrame(data=seq_df[~seq_df["Res_N"].isin(preds["Res_N"])],columns=preds.columns)
    #preds = preds.append(tmp)
    preds = pd.concat([preds, tmp])
    if len(tmp.index)>0:
        print(("%d residues from the sequence were missing "
                            "from the predictions") % len(tmp.index))

    # If Res_name is inconsistent between seq_df and preds, make a compromise
    preds.index = preds["Res_N"]
    preds = preds.sort_index()
    seq_df.index = seq_df["Res_N"]
    seq_df = seq_df.sort_index()
    mask = preds["Res_name"] != seq_df["Res_name"]

    if any(mask):
        ambiguous_res_names = (seq_df.loc[mask, "Res_name"] + "(" +
                                preds.loc[mask, "Res_type"] + "?)")
        preds.loc[mask, "Res_name"] = ambiguous_res_names
        seq_df.loc[mask, "Res_name"] = ambiguous_res_names
        # Note we don't update the Res_type in preds, because this is potentially
        # used for correcting the chemical shift
        print(
                ("There were inconsistencies between the provided sequence "
                "file and the predicted shifts. The following residues had "
                "inconsistent amino acid types: %s") %
                ", ".join(ambiguous_res_names))

#        preds.index = preds["Res_name"]
#        preds.index.name = None
    seq_df.index = seq_df["Res_name"]
    seq_df.index.name = None

    #### Add the chemical shift info back in

    # Make columns for the i-1 predicted shifts of C, CA and CB
    preds_m1 = preds[list({"C","CA","CB","Res_type","Res_name"}.
                            intersection(preds.columns))].copy()
    preds_m1.index = preds_m1.index+1
    preds_m1.columns = preds_m1.columns + "_m1"
    preds = pd.merge(preds, preds_m1, how="left",
                        left_index=True, right_index=True)

    # Make column for the i+1 Res_name
    preds_p1 = preds[["Res_name"]].copy()
    preds_p1.index = preds_p1.index-1
    preds_p1.columns = ["Res_name_p1"]
    preds = pd.merge(preds, preds_p1, how="left",
                        left_index=True, right_index=True)

    # Set index to Res_name
    preds.index = preds["Res_name"]
    preds.index.name = None

    # Restrict to only certain atom types
    atom_set = {"H","N","C","CA","CB","C_m1","CA_m1","CB_m1","HA"}
    preds = preds[["Res_name","Res_N","Res_type","Res_name_m1",
                    "Res_name_p1","Res_type_m1"]+
                    list(atom_set.intersection(preds.columns))]

    print("Finished reading in %d predicted residues from %s"
                        % (len(preds.index), filename))

    return(preds)

parser = argparse.ArgumentParser(
        description="Analyse error distribution of predicted shifts")
parser.add_argument("SNAPS_path", help="Path to the top-level SNAPS directory.")
parser.add_argument("-p", "--python_cmd", default="python")
parser.add_argument("-N", default=None, help="Limit to first N datasets.")
args = parser.parse_args()

# path = Path("C:/Users/alexh/GitHub/SNAPS/")
path = Path(args.SNAPS_path)

# Import information on the ShiftX2 testset
testset_df = pd.read_table(path/"data/testset/testset.txt", header=None, 
                           names=["ID","PDB","BMRB","Resolution","Length"])
testset_df["obs_file"] = [path/"data/testset/simplified_BMRB"/file 
                      for file in testset_df["BMRB"].astype(str)+".txt"]
testset_df["preds_file"] = [path/"data/testset/shiftx2_results"/file 
                      for file in testset_df["ID"]+"_"+testset_df["PDB"]+".cs"]
testset_df.index = testset_df["ID"]

if args.N is not None:
    testset_df = testset_df.iloc[0:int(args.N),:]

# Import all observed shifts
atom_set = {"H","N","HA","C","CA","CB","C_m1","CA_m1","CB_m1"}

obs_all = None
for i in testset_df.ID:
    obs = import_testset_shifts(testset_df.loc[i, "obs_file"], remove_Pro=False)
    
    #Change Bs to Cs
    obs.loc[obs["Res_type"]=="B", "Res_type"] = "C"

    # Convert wide to long
    obs = obs.melt(id_vars=["SS_name", "Res_N", "Res_type", "Res_type_m1"],
                   value_vars=set(obs.columns).intersection(atom_set), 
                   var_name="Atom_type", value_name="Shift")
    
    obs["ID"] = i

    if obs_all is None:
        obs_all = obs.copy()
    else:
        obs_all = pd.concat([obs_all, obs], ignore_index=True)
    
# Import all predicted shifts
preds_shiftx2 = None
for i in testset_df["ID"]:
    preds = import_pred_shifts(testset_df.loc[i, "preds_file"], 
                                        filetype="shiftx2")
    
    #Change Bs to Cs
    preds.loc[preds["Res_type"]=="B", "Res_name"] = (
            preds.loc[preds["Res_type"]=="B", "Res_name"].str.replace("B","C"))
    preds.loc[preds["Res_type"]=="B", "Res_type"] = "C"
    preds.loc[preds["Res_type_m1"]=="B", "Res_type_m1"] = "C"
    

    #Convert wide to long
    preds = preds.melt(id_vars=["Res_name", "Res_N", "Res_type", "Res_type_m1"],
                   value_vars=set(preds.columns).intersection(atom_set), 
                   var_name="Atom_type", value_name="Shift")
    
    preds["ID"] = i

    if preds_shiftx2 is None:
        preds_shiftx2 = preds.copy()
    else:
        preds_shiftx2 = pd.concat([preds_shiftx2, preds], ignore_index=True)

# Analyse the overall distribution of the real shifts
i_atoms = {"H","N","HA","C","CA","CB"}

# It would be nice to draw the 2.5% and 97.5% quantiles on the data, but I've not found a way to do that and facet nicely.
# low_quantiles = obs_all[obs_all.Atom_type.isin(i_atoms)].groupby("Atom_type").Shift.quantile(0.025)
# high_quantiles = obs_all[obs_all.Atom_type.isin(i_atoms)].groupby("Atom_type").Shift.quantile(0.975)

obs_dist_plot = ggplot(obs_all[obs_all.Atom_type.isin(i_atoms)], aes(x="Shift", fill="Atom_type")) + geom_density()
obs_dist_plot += facet_wrap("Atom_type", scales="free")
obs_dist_plot += scale_x_reverse()
# obs_dist_plot += geom_vline(xintercept=low_quantiles)
# obs_dist_plot += geom_vline(xintercept=high_quantiles)
obs_dist_plot += ggtitle("Distribution of chemical shifts for all atom types in ShiftX2 testset")
obs_dist_plot.save(path/"plots/error_dist/observed shift distribution.pdf", height=200, width=200, units="mm")

# plt = ggplot(obs_all[obs_all.Atom_type.isin(i_atoms)], aes(y="Shift")) + geom_boxplot() 
# plt += facet_wrap("Atom_type", scales="free")
# plt += coord_flip()
# plt.save(path/"plots/error_dist/observed shift boxplot.pdf", height=200, width=200, units="mm")

# Plot correlation between different shifts
obs_all["ID_SS"] = obs_all.ID + "_" + obs_all.SS_name
obs_all_wide = obs_all.pivot(index="ID_SS", columns="Atom_type", values="Shift")
obs_correlation = obs_all_wide.corr()
obs_correlation.to_csv(path/"output"/"error_dist"/"Observed shift correlation.csv")

atoms = list(atom_set)
N = len(atoms)
for i in range(N):
    for j in range(i+1,N):
        plt = ggplot(obs_all_wide, aes(x=atoms[i], y=atoms[j])) + geom_point()
        plt = plt + stat_smooth(method="lm")
        plt = plt + ggtitle("Correlation between observed shifts for atoms "+atoms[i]+" and "+atoms[j]+f". r = {obs_correlation.loc[atoms[i], atoms[j]]:.2f}")
        plt = plt + scale_x_reverse() + scale_y_reverse()
        plt.save(path/"plots/error_dist"/("correlation between "+atoms[i]+" and "+atoms[j]+".pdf"), height=200, width=200, units="mm")

# Analyse the error distribution of each set of predicted shifts
comparison_dict = {"shiftx2":[obs_all, preds_shiftx2]}

for out_dir in comparison_dict:
    obs = comparison_dict[out_dir][0]
    preds = comparison_dict[out_dir][1]

    (path/"plots/error_dist"/out_dir).mkdir(parents=True, exist_ok=True)   # Make output directory

    # Merge the obs and preds dataframes, and clean up
    df = pd.merge(obs, preds, on=["ID","Res_N","Res_type", "Res_type_m1", "Atom_type"], 
                how="outer", suffixes=["_obs","_pred"])
    df_raw = df.copy()
    df = df.dropna(subset=["Shift_obs","Shift_pred"])   # Get rid of lines where either obs or preds is missing

    print(df.groupby("Atom_type").count())      # Print a summary of the imported data

    df["Delta"] = df.Shift_pred - df.Shift_obs

    # Plot all predicted vs observed shifts
    plt = ggplot(df, aes(x="Shift_obs", y="Shift_pred")) + geom_point()
    plt = plt + facet_wrap("Atom_type", scales="free")
    plt = plt + scale_x_reverse() + scale_y_reverse()
    plt.save(path/"plots/error_dist"/out_dir/"predictions vs observations.pdf", height=200, width=200, units="mm")

    # Plot the distribution of observed and predicted shifts
    plt = ggplot(df) 
    plt = plt + geom_density(aes(x="Shift_obs"), fill="red", alpha=0.5) 
    plt = plt + geom_density(aes(x="Shift_pred"), fill="blue", alpha=0.5)
    plt += facet_wrap("Atom_type", scales="free")
    plt += scale_x_reverse()
    plt += ggtitle("Distribution of observed (red) and predicted (blue) chemical shifts")
    plt.save(path/"plots/error_dist"/out_dir/"predicted shift distribution.pdf", height=200, width=200, units="mm")

    # For each atom type, faceted by residue type...
    for a in atom_set:
        limits = [df.loc[df.Atom_type==a,["Shift_obs","Shift_pred"]].max().max(),
                    df.loc[df.Atom_type==a,["Shift_obs","Shift_pred"]].min().min()]
        
        # Plot predicted vs observed shift
        plt = ggplot(df[df.Atom_type==a], aes(x="Shift_obs", y="Shift_pred", color="Delta")) + geom_point()
        plt = plt + geom_abline(intercept=0, slope=1)
        plt = plt + facet_wrap("Res_type")
        plt = plt + xlim(limits) + ylim(limits) # + scale_x_reverse() + scale_y_reverse() 
        plt.save(path/"plots/error_dist"/out_dir/("predictions vs observations by residue - "+a+".pdf"), height=200, width=200, units="mm")
        
        # Plot distributions of observed and predicted shifts
        plt = ggplot(df[df.Atom_type==a]) 
        plt = plt + geom_density(aes(x="Shift_obs"), fill="red", alpha=0.5) 
        plt = plt + geom_density(aes(x="Shift_pred"), fill="blue", alpha=0.5)
        plt += facet_wrap("Res_type")
        plt += scale_x_reverse()
        plt += ggtitle("Distribution of observed (red) and predicted (blue) chemical shifts")
        plt.save(path/"plots/error_dist"/out_dir/("predicted shift distribution by residue - "+a+".pdf"), height=200, width=200, units="mm")
        
        # # Plot the distribution of errors, compared to the overall distribution of observations
        # plt = ggplot(df[df.Atom_type==a], aes(x="Delta")) + geom_density(fill="red", alpha=0.5)
        # plt = plt + geom_density(aes(x="Shift_obs - Shift_obs.median()"), fill="blue", alpha=0.5)
        # plt = plt + facet_wrap("Res_type")
        # plt = plt + scale_x_reverse()
        # plt = plt + ggtitle("Error distribution (red) compared to observed shift distribution (blue)")
        # plt.save(path/"plots/error_dist"/out_dir/("error distribution by residue - "+a+".pdf"), height=200, width=200, units="mm")
        # # I would like to plot the istribution of observed shifts exactly over the error distribution, but haven't found a way

        # Plot overlay of error distributions of each residue type
        plt = ggplot(df[df.Atom_type==a]) + geom_density(aes(x="Delta", color="Res_type"))
        plt = plt + ggtitle("Error distributions for each residue type for atom "+a)
        plt.save(path/"plots/error_dist"/out_dir/("error distribution by residue - "+a+".pdf"), height=200, width=200, units="mm")

        # Plot the prediction error vs observed shift
        plt = ggplot(df[df.Atom_type==a], aes(x="Shift_obs", y="Delta", color="Shift_pred")) + geom_point()
        plt = plt + stat_smooth(method="lm")
        plt = plt + facet_wrap("Res_type")  # , scales="free")
        plt = plt + scale_x_reverse()
        plt.save(path/"plots/error_dist"/out_dir/("delta vs observations by residue - "+a+".pdf"), height=200, width=200, units="mm")
        
        # Plot the prediction error vs predicted shift
        plt = ggplot(df[df.Atom_type==a], aes(x="Shift_pred", y="Delta", color="Shift_obs")) + geom_point()
        plt = plt + stat_smooth(method="lm")
        plt = plt + facet_wrap("Res_type")  # , scales="free")
        plt = plt + scale_x_reverse()
        plt.save(path/"plots/error_dist"/out_dir/("delta vs predictions by residue - "+a+".pdf"), height=200, width=200, units="mm")
        
    ## Calculate correlation between different errors
    df["ID_Res"] = df.ID + "_" + df.Res_name
    delta_wide = df.pivot(index="ID_Res", columns="Atom_type", values="Delta")

    # Output the mean and covariance to .csv files (Can be used with the delta_correlation option in SNAPS config file)
    delta_wide.mean().to_csv(path/"output"/"error_dist"/(out_dir+"_d_mean.csv"))
    delta_wide.cov().to_csv(path/"output"/"error_dist"/(out_dir+"_d_cov.csv"))
    correlation = delta_wide.corr()
    correlation.to_csv(path/"output"/"error_dist"/(out_dir+"_d_corr.csv"))

    # Work out covariance for each residue type, in case it is very different.
    delta_wide_res = {}
    correlation_res = {}
    for r in df.Res_type.unique():
        delta_wide_res[r] = df[df.Res_type==r].pivot(index="ID_Res", columns="Atom_type", values="Delta")
        correlation_res[r] = delta_wide_res[r].corr()
        (correlation_res[r]-correlation).to_csv(path/"output"/"error_dist"/(out_dir+" difference in correlation matrix for residue "+r+".csv"))

    # Plot the correlation between errors for each atom type
    atoms = list(atom_set)
    N = len(atoms)
    for i in range(N):
        for j in range(i+1,N):
            plt = ggplot(delta_wide, aes(x=atoms[i], y=atoms[j])) + geom_point()
            plt = plt + stat_smooth(method="lm")
            plt = plt + ggtitle("Correlation between prediction errors for atoms "+atoms[i]+" and "+atoms[j]+f". r = {correlation.loc[atoms[i], atoms[j]]:.2f}")
            plt = plt + scale_x_reverse() + scale_y_reverse()
            plt.save(path/"plots/error_dist"/out_dir/("correlation between "+atoms[i]+" and "+atoms[j]+".pdf"), height=200, width=200, units="mm")

    # Output the standard deviation of the error for each atom type and residue
    tmp = df.groupby("Atom_type").Delta
    tmp.std().to_csv(path/"output"/"error_dist"/(out_dir+" standard deviation of prediction error (by atom type).csv"))
    (tmp.quantile(0.75) - tmp.quantile(0.25)).to_csv(path/"output"/"error_dist"/(out_dir+" IQR of prediction error (by atom type).csv"))
    tmp2 = df.groupby(["Atom_type","Res_type"]).Delta
    tmp2.std().to_csv(path/"output"/"error_dist"/(out_dir+" standard deviation of prediction error (by atom and residue type).csv"))
    (tmp2.quantile(0.75) - tmp2.quantile(0.25)).to_csv(path/"output"/"error_dist"/(out_dir+" IQR of prediction error (by atom and residue type).csv"))


    # ## For each residue and atom type, calculate a linear fit of predicted vs observed shift
    # df["Res_type_atom"] = df["Res_type"]    # Make a new column with the residue type of the specific atom 
    #                                         # (ie the i-1 residue type for the i-1 atoms)
    # mask = df.Atom_type.isin(["C_m1,CA_m1,CB_m1"])
    # df.loc[mask, "Res_type_atom"] = df.loc[mask, "Res_type_m1"]

    # df["Delta_cor"] = np.nan
    # lm_results = pd.DataFrame(columns=["Atom_type","Res_type","Grad","Offset","Corr_coeff"])

    # for atom in atom_set:
    #     for res in df["Res_type_atom"].unique():
    #         mask = (df["Atom_type"]==atom) & (df["Res_type_atom"]==res)
    #         tmp = df.loc[mask, ["Shift_obs", "Shift_pred"]].dropna(how="any")
    #         try:
    #             lm = linregress(tmp["Shift_obs"], tmp["Shift_pred"])
    #             corr_coeff = tmp.corr().iloc[0,1]
    #             lm_results.loc[atom+"_"+res, :] = [atom, res, lm[0], lm[1], corr_coeff]
    #             # df.loc[mask,"Shift_pred_cor"] = (df.loc[mask,"Shift_pred"] 
    #             #                                 - lm[0]*df.loc[mask,"Shift_pred"]
    #             #                                 - lm[1])
    #             # df.loc[mask,"Delta_cor"] = (df.loc[mask,"Shift_pred_cor"] - 
    #             #                             df.loc[mask,"Shift_pred"])
    #         except:
    #             print("Error: ",res, atom)


    # # Plot the distribution of corrected predictions and original predictions

    

# Compare the errors from different prediction methods