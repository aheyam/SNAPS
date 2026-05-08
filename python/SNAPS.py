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

import pdb

def get_arguments(system_args):
    import argparse

    parser = argparse.ArgumentParser(
            description="SNAPS (Simple NMR Assignments from Predicted Shifts)")

    # Mandatory arguments
    parser.add_argument("shift_file",
                        help="A table of observed chemical shifts.")
    parser.add_argument("pred_file",
                        help="A table of predicted chemical shifts.")
    parser.add_argument("output_file",
                        help="The file results will be written to.")

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
    #TODO: Need to rethink how SS_class info is imported.

    # Options controlling output files
    parser.add_argument("-l", "--log_file", default=None,
                        help="A file logging information will be written to.")
    parser.add_argument("--shift_output_file", default=None,
                        help="""The file the assigned shiftlist will be written to.""")
    parser.add_argument("--shift_output_type", default="sparky",
                        choices=["sparky", "xeasy", "nmrpipe"],
                        help="One or more output formats for chemical shift export")
    parser.add_argument("--shift_output_confidence", nargs="*",
                        choices=["High","Medium","Low","Unreliable","Undefined"],
                        default=["High","Medium","Low","Unreliable","Undefined"],
                        help="""Limits the shiftlist output to assignments with
                        particular confidence levels. More than one level is allowed""")
    parser.add_argument("--alt_assignments_output_file", default=None,
                        help="""A file alternative assignments will be written to 
                        (only if alt_assignments is set to >0 in config file)""")
    parser.add_argument("--strip_plot_file",
                        default=None,
                        help="A filename for an output strip plot.")
    parser.add_argument("--hsqc_plot_file",
                        default=None,
                        help="A filename for an output HSQC plot.")



    args = parser.parse_args(system_args)
    if False:   # For convenience when testing
        args = parser.parse_args(("data/P3a_L273R/naps_shifts.txt",
                                  "data/P3a_L273R/shiftx2.cs",
                                  "output/test.txt",
                                  "--shift_type","snaps",
                                  "--pred_type","shiftx2",
                                  "-c","config/config_yaml_2.txt",
                                  "-l","output/test.log",
                                  "--strip_plot_file", "output/strip_plot.htm",
                                  "--hsqc_plot_file", "output/hsqc_plot.htm"))
    return(args)

def runSNAPS(system_args):

    #### Command line arguments
    args = get_arguments(system_args)

    #### Set up logging
    if args.log_file is not None:
        # Create a logger
        logger = logging.getLogger("SNAPS")
        logger.setLevel(logging.DEBUG)
        # Create a log handler that writes to a specific file.
        # In principle you could have multiple handlers, but here I just have one.
        # Need to explicitly define a handler so it can be explicitly closed
        # once the analysis is complete.
        log_handler = logging.FileHandler(args.log_file, mode='w')
        log_handler.setLevel(logging.DEBUG)
        #log_handler.setFormatter(logging.Formatter("%(levelname)s %(asctime)s %(module)s %(funcName)s %(message)s"))
        log_handler.setFormatter(logging.Formatter(
                "%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S"))
        logger.addHandler(log_handler)

    else:
        logging.basicConfig(level=logging.ERROR)
        # logging.basicConfig(level=logging.DEBUG)


    #### Set up the SNAPS_assigner object
    a = SNAPS_assigner()

    # Import config file
    a.read_config_file(args.config_file)
    # breakpoint()


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

    if a.pars["iterate_until_consistent"]:
        a.assign_df = a.find_consistent_assignments(set_assign_df=True)
    else:
        a.assign_from_preds(set_assign_df=True)
        a.add_consistency_info(threshold=a.pars["seq_link_threshold"])
        if (a.pars["alt_assignments"] > 0):
            if args.alt_assignments_output_file is not None:
                a.find_alt_assignments(N=a.pars["alt_assignments"])
            else:
                logger.warning("alt_assignments > 0 but no alt_assignments_output_file defined - skipping alt assignment")

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
    a.assign_df.to_csv(args.output_file, sep="\t", float_format="%.3f",
                           index=False)
    logger.info("Finished writing results to %s", args.output_file)

    if (a.pars["alt_assignments"] > 0) and args.alt_assignments_output_file is not None:
        a.alt_assign_df.to_csv(args.alt_assignments_output_file, sep="\t", float_format="%.3f",
                                index=False)
        logger.info("Finished writing alternative assignment results to %s", args.alt_assignments_output_file)

    #### Write chemical shift lists
    if args.shift_output_file is not None:
        a.output_shiftlist(args.shift_output_file, args.shift_output_type,
                           confidence_list=args.shift_output_confidence)

    #### Make some plots
    plots = []
    if args.hsqc_plot_file is not None:
        hsqc_plot = a.plot_hsqc(args.hsqc_plot_file, "html")
        logger.info("Finished writing HSQC plot to %s", args.hsqc_plot_file)
        plots += [hsqc_plot]

    if args.strip_plot_file is not None:
        strip_plot = a.plot_strips(args.strip_plot_file, "html")
        logger.info("Finished writing strip plot to %s", args.strip_plot_file)
        plots += [strip_plot]

    # Close the log file
    logger.handlers[0].close()
    logger.removeHandler(logger.handlers[0])

    return(plots)


#%% Run the actual script
if __name__ == '__main__':
    import sys

    runSNAPS(sys.argv[1:])
