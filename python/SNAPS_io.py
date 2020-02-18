# -*- coding: utf-8 -*-
"""
Functions that convert input text files to pandas tables, or which export the 
calculation results to text files.

@author: Alex
"""

import numpy as np
import pandas as pd
from Bio.SeqUtils import seq1
from Bio import SeqIO

# Set up logging
# This allows logging messages to be registered, but they won't be written 
# anywhere until a handler is defined elsewhere
import logging
logger = logging.getLogger(__name__).addHandler(logging.NullHandler())

#%% Import functions

def import_obs_shifts(filename, filetype, SS_num=False):
    """ Import a chemical shift list
    
    filename: Path to text file containing chemical shifts.
    filetype: Allowed values are "snaps", "ccpn", "sparky", "xeasy", 
        "nmrpipe" or "mars"
        The "ccpn" option is for importing a Resonance table exported from 
        Analysis v2.x. The "snaps" option is for importing an unassigned 
        shift table previously exported from SNAPS
    SS_num: If true, will extract the longest number from the SS_name and 
    treat it as the residue number. Without this, it is not possible to get
    the i-1 shifts for each spin system, in cases where the observed residues 
    have already been assigned.  
    """
    # Import from file
    if filetype=="snaps":
        obs = pd.read_csv(filename, sep="\t")
    elif filetype=="ccpn":
        obs = pd.read_csv(filename, sep="\t")
        obs = obs.loc[:,["Residue", "Assign Name", "Shift"]]
        obs.columns = ["SS_name", "Atom_type", "Shift"]
        obs["Atom_type"] = obs["Atom_type"].str.upper()
    elif filetype=="sparky":
        obs = pd.read_csv(filename, sep="\s+")
        obs = obs.loc[:,["Group", "Atom", "Shift"]]
        obs.columns = ["SS_name", "Atom_type", "Shift"]
        obs.loc[obs["Atom_type"]=="HN", "Atom_type"] = "H"
    elif filetype=="xeasy":
        obs = pd.read_csv(filename, sep="\s+", 
                            header=None, na_values="999.000",
                            names=["i","Shift","SD","Atom_type","SS_name"])
        obs = obs.loc[:, ["SS_name", "Atom_type", "Shift"]]
        obs["SS_name"] = obs["SS_name"].astype(str)
        obs = obs.dropna(subset=["Shift"])
        obs.loc[obs["Atom_type"]=="HN", "Atom_type"] = "H"
    elif filetype=="nmrpipe":
        # Work out where the column names and data are
        with open(filename, 'r') as f:
            for num, line in enumerate(f, 1):
                if line.find("VARS")>-1:
                    colnames_line = num
        
        obs = pd.read_csv(filename, sep="\s+", skiprows=colnames_line+1, 
                            names=["SS_name","Res_type","Atom_type","Shift"])
        obs = obs.loc[:, ["SS_name", "Atom_type", "Shift"]]
        obs["SS_name"] = obs["SS_name"].astype(str)
        obs.loc[obs["Atom_type"]=="HN", "Atom_type"] = "H"
    elif filetype=="mars":
        obs_wide = pd.read_csv(filename, sep="\s+", na_values="-")
        obs_wide = obs_wide.rename(columns={"CO":"C","CO-1":"C_m1",
                                            "CA-1":"CA_m1","CB-1":"CB_m1"})
        obs_wide = obs_wide.drop(columns="HA-1")
        obs_wide["SS_name"] = obs_wide.index.astype(str)
        obs = obs_wide.melt(id_vars="SS_name", 
                            value_vars=["H","HA","N","C","CA",
                                        "CB","C_m1","CA_m1","CB_m1"], 
                            var_name="Atom_type", value_name="Shift")       
#        elif filetype == "nmrstar":
#            tmp = nmrstarlib.read_files(filename)
#            return(tmp)
    else:
        print("import_obs_shifts: invalid filetype '%s'." % (filetype))
        return(None)
    
    # Restrict to backbone atom types
    obs = obs.loc[obs["Atom_type"].isin(["H","HA","N","C","CA","CB",
                                         "C_m1","CA_m1","CB_m1"]),:]
    
    # Convert from long to wide
    obs = obs.pivot(index="SS_name", columns="Atom_type", values="Shift")
    obs.insert(0, "SS_name", obs.index.values)
    obs.index.name = None
    
    if SS_num:
        # Extract residue number from SS_name, and get m1 shifts
        obs["Res_N"] = obs["SS_name"].str.extract(r"(\d+)").astype(int)
        
        obs.index = obs["Res_N"]
        obs_m1 = obs[list({"C","CA","CB"}.intersection(obs.columns))]
        obs_m1.index = obs_m1.index+1
        obs_m1.columns = obs_m1.columns + "_m1"
        obs = pd.merge(obs, obs_m1, how="left", 
                       left_index=True, right_index=True)
        obs = obs.drop(columns="Res_N")
        
        # Set index back to SS_name
        obs.index = obs["SS_name"]
        obs.index.name = None
          
    return(obs)
    
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
    # Import the observed chemical shifts
    obs_long = pd.read_csv(filename, sep="\t")
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
    
def import_sequence(filename):
    """ Imports the protein sequence from a fasta file 
    
    The sequence information is used to fill in any gaps if prediction 
    information is missing. It can also be used if the residue numbering is 
    discontinuous. The sequence file should be in FASTA format, and the 
    sequence id should be either:
        1) the residue number of the first amino acid, or
        2) a comma separated list of residue ranges eg -5:10,15:100 (this is 
        for proteins with discontinuous sequence numbering)
    If the file contains multiple sequences, only the first will be used.
    
    Parameters
    filename: path to file containing sequence information
    """
    
    fasta_records = SeqIO.parse(open(filename),"fasta")
    record1 = next(fasta_records)
    
    # Parse the sequence
    seq_list = list(str(record1.seq))
    
    # Parse the fasta id and make a list of residue numbers
    tmp = record1.id.split(",")
    if len(tmp)==1:     # If id was "123" or "123-456"
        tmp2 = tmp[0].split(":")
        if len(tmp2)==1:        # If id was "123"
            res_N_start = int(tmp2[0])
            res_N_list = list(range(res_N_start, res_N_start+len(seq_list)))
        else:                   # If id was "123-456"
            start, end = tmp2
            res_N_list = list(range(int(start), int(end)+1))
    else:               # If id was "123-200,300-456"
        res_N_list = []
        for x in tmp:
            start, end = x.split(":")
            res_N_list += list(range(int(start), int(end)+1))
    
    # Check res_N_list and seq_list are the same length, and correct if not
    if len(res_N_list) > len(seq_list):
        res_N_list = res_N_list[0:len(seq_list)]
        logger.warning("The residue range provided is longer than the"+
                            " length of the sequence, so has been truncated")
    elif len(res_N_list) < len(seq_list):
        extra_needed = len(seq_list) - len(res_N_list)
        last_N = res_N_list[-1]
        res_N_list += list(range(last_N, last_N + extra_needed))
        logger.warning("The residue range provided is shorter than the"+
                            " length of the sequence, so has been extended")
    
    # Make a dataframe
    seq_df = pd.DataFrame({"Res_N":res_N_list,"Res_type":seq_list})
    seq_df["Res_name"] = seq_df["Res_N"].astype(str) + seq_df["Res_type"]
    seq_df["Res_name"] = seq_df["Res_name"].str.rjust(5)
                                        # Pad Res_name to constant length
    
    seq_df.index = seq_df["Res_name"]
    seq_df.index.name = None
    
    logger.info("Finished importing a sequence of length %d" % 
                     len(seq_df.index))
    
    return(seq_df)
    
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
                logger.warning(
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
                    
        preds_long = pd.read_csv(filename, sep="\s+", names=colnames,
                                   skiprows=colnames_line+1)
        preds_long = preds_long.reindex(columns=["RESID","RESNAME",
                                                 "ATOMNAME","SHIFT"])
        preds_long.columns = ["Res_N","Res_type","Atom_type","Shift"]
        
        # Sparta+ uses HN for backbone amide proton - convert to H
        preds_long.loc[preds_long["Atom_type"]=="HN", "Atom_type"] = "H"
    else:
        logger.error("""Invalid predicted shift type: '%s'. Allowed 
                          options are 'shiftx2' or 'sparta+'""" % (filetype))
        return(None)
    
    logger.info("Imported %d predicted chemical shifts from %s" 
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
    # TODO: This section should probably merge into prepare_obs_preds()
    # If seq_df is missing, create it based on preds
    seq_df = self.seq_df
    
    if seq_df is None:
        seq_df = preds.copy()[["Res_N","Res_type","Res_name"]]
        
        # If there are any missing residue numbers, create them
        min_N = seq_df["Res_N"].min()
        max_N = seq_df["Res_N"].max()
        missing_residue_numbers = (set(range(min_N, max_N+1)).
                                   difference(seq_df["Res_N"]))
        if len(missing_residue_numbers)>0:
            tmp = pd.DataFrame({"Res_N":list(missing_residue_numbers),
                                "Res_type":"X","Res_name":np.NaN})
            tmp["Res_name"] = tmp["Res_N"].astype(str) + tmp["Res_type"]
            tmp["Res_name"] = tmp["Res_name"].str.rjust(5)
            tmp.index = tmp["Res_N"]
            tmp.index.name = None
            seq_df = seq_df.append(tmp).sort_index()
            
    # Add/delete residues from preds so it matches seq_df
    # Only keep predictions that are in seq_df
    tmp = len(preds.index)
    preds = preds[preds["Res_N"].isin(seq_df["Res_N"])]
    tmp2 = tmp - len(preds.index)
    if tmp2>0:
        self.logger.info(("Predictions for %d residues were discarded "
                            "because they were not present in the imported "
                            "sequence file") % tmp2)
        
    # Add predictions for any residue number that is in seq_df but not preds
    tmp = pd.DataFrame(data=seq_df[~seq_df["Res_N"].isin(preds["Res_N"])],columns=preds.columns)
    preds = preds.append(tmp)
    if len(tmp.index)>0:
        self.logger.info(("%d residues from the sequence were missing "
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
        self.logger.warning(
                ("There were inconsistencies between the provided sequence "
                "file and the predicted shifts. The following residues had "
                "inconsistent amino acid types: %s") % 
                ", ".join(ambiguous_res_names))
    
#        preds.index = preds["Res_name"]
#        preds.index.name = None
    seq_df.index = seq_df["Res_name"]
    seq_df.index.name = None
    
    self.seq_df = seq_df.copy()
    
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
    
    logger.info("Finished reading in %d predicted residues from %s"
                     % (len(preds.index), filename))
    
    return(preds)
    
def simulate_pred_shifts(filename, sd, seed=None, 
                         atom_set={"H","N","C","CA","CB","C_m1","CA_m1","CB_m1","HA"}):
    """Generate a 'simulated' predicted shift DataFrame by importing some 
    observed chemical shifts (in 'test' format), and adding Gaussian 
    errors.
    
    sd: dictionary of standard deviation for each atom type 
        eg. {"H":0.1,"N":0.5}
    """
    import numpy.random as rand
    
    if seed is not None:
        rand.seed(seed)
        
    preds = import_testset_shifts(filename, remove_Pro=False)
    
    # Limit columns
    preds = preds[["Res_N","Res_type"]+
                  list({"C","CA","CB","H","N","HA"}.intersection(preds.columns))]
    
    # Add the random shifts
    for atom in atom_set.intersection(preds.columns):
        preds.loc[:,atom] += rand.normal(0, sd[atom], size=len(preds.index))
    
    # Add other columns back in
    preds.insert(1, "Res_name", (preds["Res_N"].astype(str) + 
              preds["Res_type"]))
    preds.loc[:,"Res_name"] = preds["Res_name"].str.rjust(5)
    
    preds.index = preds["Res_N"]
    #preds.index.name=None
    
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
    # I think this is unnecessary, due to duplication at the top of the function?
#    atom_set_all = {"H","N","C","CA","CB","C_m1","CA_m1","CB_m1","HA"}
#    preds = preds[["Res_name","Res_N","Res_type","Res_name_m1",
#                   "Res_name_p1","Res_type_m1"]+
#                  list(atom_set_all.intersection(preds.columns))]
    
    return(preds)
    
#%% Export functions

def output_shiftlist(assign_df, filepath=None, format="sparky", all_preds=None,
                     confidence_list=["High","Medium","Low","Unreliable","Undefined"]):
    """Export a chemical shift list, in a variety of formats
    
    Parameters
    filepath: the file the shifts will be written to
    format: the format of the chemical shift file (sparky, xeasy or nmrpipe)
    all_preds: a DataFrame containing all predictions, from a SNAPS_assigner 
        object. Only required if format='nmrpipe'.
    confidence_list: Only residues with confidence in this list will be output
    """
    
    # Limit confidence levels and discard unneeded columns and rows
    atoms = {"H","N","HA","C","CA","CB"}.intersection(assign_df.columns)
    df_wide = assign_df.loc[assign_df["Confidence"].isin(confidence_list),:]
    logger.info("Chemical shift export: restricted confidence levels to %s" 
                     % ", ".join(confidence_list))
    df_wide.loc[~(df_wide["Dummy_res"] | df_wide["Dummy_SS"]),
                ["Res_N","Res_type"]+list(atoms)]
    #TODO: This discards all i-1 information, even though that may be useful 
    # in some cases (eg for getting proline shifts). Is it possible to make 
    # better use of this?
    
    # Check in case no shifts were selected
    if df_wide.empty:   
        print("No chemical shifts selected for export.")
        logger.warning("Cannot export chemical shifts: No atom types selected")
        # Create an empty file
        df_wide.to_csv(filepath, sep="\t", float_format="%.3f",
                       index=True, header=False)
        return(None)
    
    # Reshape dataframe from wide to long format, sort and remove NAs
    df = df_wide.melt(id_vars=["Res_N","Res_type"], value_vars=atoms, 
                      var_name="Atom_type", value_name="Shift")
    df = df.sort_values(["Res_N","Atom_type"])
    df = df.dropna(subset=["Res_N"])
    df["Res_N"] = df["Res_N"].astype(int)
    
    # Make format-specific modifications and export
    if format=="sparky":
        df["Group"] = df["Res_type"] + df["Res_N"].astype(str)
        
        df = df.rename(columns={"Atom_type":"Atom"})
        df.loc[df["Atom"]=="H","Atom"] = "HN"
        
        nuc_dict = {"HN":"1H", "N":"15N", "HA":"1H",
                    "C":"13C", "CA":"13C", "CB":"13C"}
        df["Nuc"] = [nuc_dict[a] for a in df["Atom"]]
        
        output_df = df[["Group","Atom","Nuc","Shift"]].copy()
        output_df["Sdev"] = 0.0
        output_df["Assignments"] = int(1)
        
        output_df = output_df.dropna()
        
        if filepath is not None:
            output_df.to_csv(filepath, sep="\t", float_format="%.3f",
                             index=False)
    elif format=="xeasy":
        df.loc[df["Atom_type"]=="H","Atom_type"] = "HN"
        
        output_df = df[["Shift","Atom_type","Res_N"]].copy()
        output_df.insert(1, "Sdev", 0)
        
        output_df["Shift"] = output_df["Shift"].fillna(999.0)
        output_df = output_df.dropna()
        output_df = output_df.reset_index(drop=True)
    
        if filepath is not None:
            output_df.to_csv(filepath, sep="\t", float_format="%.3f",
                             index=True, header=False)
    elif format=="nmrpipe":
        df.loc[df["Atom_type"]=="H","Atom_type"] = "HN"
        
        output_df = df[["Res_N", "Res_type", "Atom_type", "Shift"]].copy()
        
        output_df["Shift"] = output_df["Shift"].fillna(9999.0)
        output_df = output_df.dropna()
        
        # NmrPipe requres a sequence, but we don't necessarily have a complete one
        # Piece togther what we can from all_preds, and put X at other positions
        tmp = pd.DataFrame({"Res_N":np.arange(all_preds["Res_N"].min(), 
                                              all_preds["Res_N"].max()+1)})
        tmp = tmp.merge(all_preds[["Res_N","Res_type"]], how="left", on="Res_N")
        tmp["Res_type"] = tmp["Res_type"].fillna("X")
        seq = tmp["Res_type"].sum()     # Convert to a string
        logger.info("Sequence contains %d residues, of which %d have unknown type"
                         % (len(seq), sum(tmp["Res_type"]=="X")))
        
        if filepath is not None:
            f = open(filepath, 'w+')
            f.write("REMARK Chemical shifts (automatically assigned using SNAPS)\n\n")
            f.write("DATA FIRST_RESID %d\n\n" % tmp["Res_N"].min())
            f.write("DATA SEQUENCE %s\n\n" % seq)
            f.write("VARS   RESID RESNAME ATOMNAME SHIFT\n")
            f.write("FORMAT %4d   %1s     %4s      %8.3f\n\n")
            
            for i in output_df.index:
                f.write("%4d %1s %4s %8.3f\n" % (output_df.loc[i, "Res_N"],
                                                 output_df.loc[i, "Res_type"],
                                                 output_df.loc[i, "Atom_type"],
                                                 output_df.loc[i, "Shift"]))
            f.close()
            
    else:
        logger.warning("Cannot export chemical shifts: "+
                            "format string '%s' not recognised." % format)
        return(None)
        
    logger.info("Wrote %d chemical shifts from %d residues to %s" %
                     (len(df.index), 
                      len(df["Res_N"].drop_duplicates()), 
                      filepath))
    return(output_df.to_csv(sep="\t"))