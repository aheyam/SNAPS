#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main SNAPS script for assigning an observed shift list based on predicted shifts

@author: aph516
"""
from tabulate import tabulate

from SNAPS_importer import SNAPS_importer
from SNAPS_assigner import SNAPS_assigner
import logging
from pathlib import Path

import pdb

# For testing
from plotnine import *

def get_arguments(system_args):
    import argparse

    parser = argparse.ArgumentParser(
            description="SNAPS (Simple NMR Assignments from Predicted Shifts)")

    # Mandatory arguments
    parser.add_argument("shift_file",
                        help="A table of observed chemical shifts.")
    parser.add_argument("pred_file",
                        help="A table of predicted chemical shifts.")
    parser.add_argument("output_dir",
                        help="The directory results will be written to.")

    # Information on input files and configuration options
    parser.add_argument("--shift_type",
                        choices=["snaps", "ccpn", "sparky", "mars",
                                 "xeasy", "nmrpipe", "nef", "test"],
                        default="snaps", 
                        help="The format of the observed shift file.")
    parser.add_argument("--pred_type",
                        choices=["shiftx2", "sparta+"],
                        default="shiftx2",
                        help="The format of the predicted shifts")
    parser.add_argument("--pred_seq_offset", type=int, default=0,
                        help="""An offset to apply to the residue numbering in
                        the predicted shifts.""")
    parser.add_argument("-c", "--config_file",
                        default="../config/config.txt",
                        help="A file containing parameters for the analysis.")
    parser.add_argument("--test_aa_classes", default=None,
                        help="""For test data only.
                        A string containing a comma-separated list of the amino acid
                        classes for the i residue, a semicolon, then a list of AA
                        classes for the i-1 residue. No spaces.
                        eg. "ACDEFGHIKLMNPQRSTVWY;G,S,T,AVI,DN,FHYWC,REKPQML" for
                        a sequential HADAMAC """)
    parser.add_argument("--simulate_pred_shifts", action="store_true",
                        help="""If present and if shift_type="test", then use simulated 
                        predicted shifts instead of the ones given in pred_file""")
    parser.add_argument("--sim_pred_multiplier", type=float, default=0.0,
                        help="""If simulate_pred_shifts is True, multiply the atom_95_quartile
                        given in the config file by this value to get the standard deviation 
                        of the simulated shifts for each atom type.""")
    parser.add_argument("--sim_pred_seed", type=int, default=0,
                        help="Random seed used if simulate_pred_shifts is True")
    parser.add_argument("--test", default=None,
                        help="Test SNAPS using example data from the specified test protein")
    #TODO: Need to rethink how SS_class info is imported.

    # Options controlling output files
    parser.add_argument("--shift_output_type", default="sparky",
                        choices=["sparky", "xeasy", "nmrpipe"],
                        help="One or more output formats for chemical shift export")
    parser.add_argument("--shift_output_confidence", nargs="*",
                        choices=["High","Medium","Low","Unreliable","Undefined"],
                        default=["High","Medium","Low","Unreliable","Undefined"],
                        help="""Limits the shiftlist output to assignments with
                        particular confidence levels. More than one level is allowed""")
    parser.add_argument("--strip_plot",
                        action="store_true",
                        help="Output a strip plot to assess assignment quality.")
    parser.add_argument("--hsqc_plot",
                        action="store_true",
                        help="Output an HSQC plot of the assignments.")



    args = parser.parse_args(system_args)
    if args.test is not None:   # For convenience when testing
        # args = parser.parse_args(("data/P3a_L273R/naps_shifts.txt",
        #                           "data/P3a_L273R/shiftx2.cs",
        #                           "output/test.txt",
        #                           "--shift_type","snaps",
        #                           "--pred_type","shiftx2",
        #                           "-c","config/config_yaml_2.txt",
        #                           "-l","output/test.log",
        #                           "--strip_plot_file", "output/strip_plot.htm",
        #                           "--hsqc_plot_file", "output/hsqc_plot.htm",
        #                           "--test"))
        import pandas as pd

        testset_df = pd.read_table("data/testset/testset.txt", header=None,
                                names=["ID","PDB","BMRB","Resolution","Length"])
        testset_df["obs_file"] = [x for x in "data/testset/simplified_BMRB/"+testset_df["BMRB"].astype(str)+".txt"]
        testset_df["preds_file"] = [x for x in "data/testset/shiftx2_results/"+testset_df["ID"]+"_"+testset_df["PDB"]+".cs"]
        testset_df["out_name"] = testset_df["ID"]+"_"+testset_df["BMRB"].astype(str)
        testset_df.index = testset_df["ID"]

        args = parser.parse_args((testset_df.loc[args.test, "obs_file"],
                                  testset_df.loc[args.test, "preds_file"],
                                  "output/test",
                                  "--shift_type","test",
                                  "--pred_type","shiftx2",
                                  "-c","config/config_yaml_2.txt",
                                  "--strip_plot",
                                  "--hsqc_plot",
                                  "--test", args.test))
    return(args)

def runSNAPS(system_args):

    #### Command line arguments
    args = get_arguments(system_args)

    # Create output directory, if it doesn't already exist
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    #### Set up logging
    # Create a logger
    logger = logging.getLogger("SNAPS")
    logger.setLevel(logging.DEBUG)
    # Create a log handler that writes to a specific file.
    # In principle you could have multiple handlers, but here I just have one.
    # Need to explicitly define a handler so it can be explicitly closed
    # once the analysis is complete.
    log_handler = logging.FileHandler(output_dir/"log.txt", mode='w')
    log_handler.setLevel(logging.DEBUG)
    #log_handler.setFormatter(logging.Formatter("%(levelname)s %(asctime)s %(module)s %(funcName)s %(message)s"))
    log_handler.setFormatter(logging.Formatter(
            "%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S"))
    logger.addHandler(log_handler)


    #### Set up the SNAPS_assigner object
    a = SNAPS_assigner()

    # Import config file
    a.read_config_file(args.config_file)

    # Import observed and predicted shifts
    importer = SNAPS_importer()

    if args.shift_type=="test":
        if args.test_aa_classes is None:
            importer.import_testset_shifts(args.shift_file)
        else:
            AA_class, AA_class_m1 = args.test_aa_classes.split(";")
            importer.import_testset_shifts(args.shift_file,
                                           SS_class=AA_class.split(","),
                                           SS_class_m1=AA_class_m1.split(","))
    else:
        importer.import_obs_shifts(args.shift_file, args.shift_type, SS_num=False)

    a.obs = importer.obs
    logger.info("Finished reading in %d spin systems from %s",
                 len(a.obs["SS_name"]), args.shift_file)

    if args.simulate_pred_shifts:
        # Calculate the errors for each atom as the 95% percentile interval multipled by 
        # the sim_pred_multiplier argument
        atom_errors = a.pars["atom_95_percentile"]
        for k in atom_errors.keys():
             atom_errors[k] = atom_errors[k]*args.sim_pred_multiplier

        a.simulate_pred_shifts(args.shift_file, atom_errors, args.sim_pred_seed)
    else:
        a.import_pred_shifts(args.pred_file, args.pred_type, args.pred_seq_offset)

    #### Do the analysis
    a.prepare_obs_preds()
    a.calc_log_prob_matrix()
    a.calc_mismatch_matrix()

    if a.pars["iterate_until_consistent"] == 1:
        a.assign_df = a.find_consistent_assignments(set_assign_df=True)
    else:
        a.assign_from_preds(set_assign_df=True)
        a.add_consistency_info(threshold=a.pars["seq_link_threshold"])
        if (a.pars["alt_assignments"] > 0):
            a.find_alt_assignments(N=a.pars["alt_assignments"])
            
    if a.pars["iterate_until_consistent"] == 2:
        high_conf_assn = a.assign_df.loc[a.assign_df.Confidence=="High", ["Res_name", "SS_name"]]
        node_df = a.find_consistent_assignments_4(threshold=a.pars["seq_link_threshold"], max_iterations=100, verbose=True, init_inc=high_conf_assn)
        best_node = (node_df.N_high + node_df.N_med).idxmax()
        best_matching = node_df.loc[best_node, "Matching"]
        b = a.copy()
        b.make_assign_df(best_matching, set_assign_df=True)
        b.add_consistency_info(threshold=0.2)
        b.assign_df.to_csv(output_dir/"consistent_assign_df.tsv", sep="\t", float_format="%.3f",
                           index=False)
        node_df.to_csv(output_dir/"node_df.tsv", sep="\t", float_format="%.3f", index=False)
        b.plot_strips(output_dir/"strip_plot_consistent.htm", "html")
        
        plt = ggplot(node_df[node_df.Ranked]) + geom_point(aes(x="Iteration",y="ID2", color="N_high+N_med"))
        plt.save(output_dir/"test_history.pdf", verbose=False)
        # breakpoint()


    
    #### Output the results
    
    # Tabulate doesn't account for if some atom types are missing.
    # headings = '''
    #     Res_name Res_N Res_type SS_name Dummy_res Dummy_SS CA CA_pred HA HA_pred H H_pred CB CB_pred
    #      C C_pred N N_pred Log_prob Max_mismatch_m1 Max_mismatch_p1 Num_good_links_m1 
    # '''.split()
    # table = []
    # for df_index, df_row in a.assign_df.iterrows():
    #     table_row = []
    #     table.append(table_row)
    #     for heading in headings:
    #         table_row.append(df_row[heading])

    # with open(args.output_file, 'w') as fp:
    #     print(tabulate(table, tablefmt='plain', headers=headings), file=fp)
    
    a.assign_df.to_csv(output_dir/"assign_df.tsv", sep="\t", float_format="%.3f",
                           index=False)
    logger.info("Finished writing results to assign_df.tsv")

    if (a.pars["alt_assignments"] > 0):
        a.alt_assign_df.to_csv(output_dir/"alt_assign_df.tsv", sep="\t", float_format="%.3f",
                                index=False)
        logger.info("Finished writing alternative assignment results to alt_assign_df.tsv")

    #### Write chemical shift lists
    a.output_shiftlist(output_dir/"assigned_shifts.txt", args.shift_output_type,
                        confidence_list=args.shift_output_confidence)

    #### Make some plots
    plots = []
    if args.hsqc_plot:
        hsqc_plot = a.plot_hsqc(output_dir/"hsqc_plot.htm", "html")
        logger.info("Finished writing HSQC plot to hsqc_plot.htm")
        plots += [hsqc_plot]

    if args.strip_plot:
        strip_plot = a.plot_strips(output_dir/"strip_plot.htm", "html")
        logger.info("Finished writing strip plot to strip_plot.htm")
        plots += [strip_plot]

    # Close the log file
    logger.handlers[0].close()
    logger.removeHandler(logger.handlers[0])


    return(plots)


#%% Run the actual script
if __name__ == '__main__':
    import sys

    runSNAPS(sys.argv[1:])
