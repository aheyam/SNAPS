#!/anaconda3/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 14:30:36 2018

@author: aph516
"""

import pandas as pd
from NAPS_assigner import NAPS_assigner
import argparse
from pathlib import Path
from distutils.util import strtobool
import logging

#### User input

# For testing
if False:
    runfile('/Users/aph516/GitHub/NAPS/python/NAPS.py', wdir='/Users/aph516/GitHub/NAPS/python',
            args="../data/testset/simplified_BMRB/4032.txt "+
            "../data/testset/shiftx2_results/A001_1KF3A.cs "+
            "../output/test.txt")
    
    runfile(Path('C:/Users/Alex/GitHub/NAPS/python/NAPS.py'), wdir=Path('C:/Users/Alex/GitHub/NAPS/python'),
            args="../data/testset/simplified_BMRB/4032.txt "+
            "../data/testset/shiftx2_results/A001_1KF3A.cs "+
            "../output/test.txt"+
            " -c config.txt")
    
    runfile(Path('C:/Users/kheyam/Documents/GitHub/NAPS/python/NAPS.py'), wdir=Path('C:/Users/kheyam/Documents/GitHub/NAPS/python'),
            args="../data/testset/simplified_BMRB/4032.txt "+
            "../data/testset/shiftx2_results/A001_1KF3A.cs "+
            "../output/test.txt"+
            " -c config.txt")

parser = argparse.ArgumentParser(description="NMR Assignments from Predicted Shifts")
parser.add_argument("shift_file")
parser.add_argument("pred_file")
parser.add_argument("out_file")
parser.add_argument("-c", "--config_file", default="/Users/aph516/GitHub/NAPS/python/config.txt")
#parser.add_argument("-c", "--config_file", default=Path("C:/kheyam/Documents/GitHub/NAPS/python/config.txt"))
parser.add_argument("-l", "--log_file", default=None)
parser.add_argument("--delta_correlation", action="store_true", 
                    help="If set, account for correlations between prediction errors of different atom types")
parser.add_argument("-a", "--alt_assignments", default=-1, type=int,
                    help="The number of alternative assignments to generate, in addition to the highest ranked.")
parser.add_argument("-p", "--plot_stem", default="/Users/aph516/GitHub/NAPS/plots/plot",
                    help="The base name for plots (do not include a file extension).")

if False:
    args = parser.parse_args(("../data/testset/simplified_BMRB/4032.txt "+
                "../data/testset/shiftx2_results/A001_1KF3A.cs "+
                "../output/test.txt " + 
                "-c ../config/config.txt "+
                "-l ../output/log.txt "+
                "-p ../plots/A001_1KF3A").split())
else:
    args = parser.parse_args()

print(args.config_file)

# Set up logging
if isinstance(args.log_file, str):
    print(args.log_file)
    logging.basicConfig(filename=args.log_file, filemode='w', level=logging.DEBUG,
                        format="%(levelname)s %(asctime)s %(module)s %(funcName)s %(message)s")
else:
    logging.basicConfig(level=logging.ERROR)

#%%

#### Set up the assigner

a = NAPS_assigner()

# Import config file
config = pd.read_table(args.config_file, sep="\s+", comment="#", header=None,
                       index_col=0, names=["Value"])
a.pars["shiftx2_offset"] = int(config.loc["shiftx2_offset"].Value)
a.pars["prob_method"] = config.loc["prob_method"].Value
a.pars["alt_assignments"] = int(config.loc["alt_assignments"].Value)
a.pars["atom_set"] = {s.strip() for s in config.loc["atom_set"].Value.split(",")}
tmp = [s.strip() for s in config.loc["atom_sd"].Value.split(",")]
a.pars["atom_sd"] = dict([(x.split(":")[0], float(x.split(":")[1])) for x in tmp])
a.pars["plot_strips"] = bool(strtobool(config.loc["plot_strips"].Value))
logging.info("Read in configuration from %s.", args.config_file)

# Account for any command line arguments that overide config file
if args.alt_assignments>=0:
    a.pars["alt_assignments"] = args.alt_assignments
if args.delta_correlation:
    a.pars["prob_method"] = "delta_correlation"


a.import_obs_shifts(args.shift_file)
logging.info("Read in %d spin systems from %s.", len(a.obs["SS_name"]), args.shift_file)

a.read_shiftx2(args.pred_file)
logging.info("Read in %d predicted residues from %s.", len(a.preds["Res_name"]), args.pred_file)

#### Do the analysis
a.add_dummy_rows()
a.calc_log_prob_matrix(sf=1, verbose=False)
logging.info("Calculated log probability matrix (%dx%d).", a.log_prob_matrix.shape[0], a.log_prob_matrix.shape[1])
assign_df, best_match_indexes = a.find_best_assignment()
logging.info("Calculated best assignment.")
a.check_assignment_consistency(threshold=0.1)
logging.info("Checked assignment consistency.")

if a.pars["alt_assignments"]>0:
    a.find_alt_assignments2(N=a.pars["alt_assignments"], verbose=False, by_ss=True)
    logging.info("Calculated the %d next best assignments for each spin system", args.alt_assignments)
    a.alt_assign_df.to_csv(args.out_file, sep="\t", float_format="%.3f")
    logging.info("Wrote results to %s", args.out_file)
else:
    a.assign_df.to_csv(args.out_file, sep="\t", float_format="%.3f")
    logging.info("Wrote results to %s", args.out_file)
#tmp = a.find_alt_assignments(best_match_indexes, by_res=False)
#tmp = a.find_alt_assignments2(N=2, verbose=True, by_ss=False)


#%%

#### Make some plots
if a.pars["plot_strips"]:
    plt = a.plot_strips()
    plt.save(args.plot_stem+"_strips.pdf", height=210, width=297, units="mm")
    logging.info("Wrote strip plot to %s", args.plot_stem+"_strips.pdf")