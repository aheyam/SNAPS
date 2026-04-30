# Test analysis of data from Jenny Tomlinson
# exec(open("test.py").read())

from tabulate import tabulate

import pandas as pd
from SNAPS_importer import SNAPS_importer
from SNAPS_assigner import SNAPS_assigner, df_lookup
import logging

import pdb

working_dir = r"C:\Users\chmahey\OneDrive - University of Leeds\Data\Jenny Tomlinson SNAPS"
output_name = r"alphafold_0 withH iterated"

a = SNAPS_assigner()
a.read_config_file("../config/config_yaml_2.txt")
a.pars["iterate_until_consistent"] = True

importer = SNAPS_importer()
importer.import_hsqc_peaks(working_dir+r"\peaklists\hsqc_2.tsv", "ccpn")
importer.import_3d_peaks(working_dir+r"\peaklists\hncaco_2.tsv", "ccpn", "hncaco")
importer.import_3d_peaks(working_dir+r"\peaklists\hnco_2.tsv", "ccpn", "hnco")
importer.import_3d_peaks(working_dir+r"\peaklists\hncacb_2.tsv", "ccpn", "hncacb")
importer.import_3d_peaks(working_dir+r"\peaklists\hnca_2.tsv", "ccpn", "hnca")
importer.import_3d_peaks(working_dir+r"\peaklists\hncocacb_2.tsv", "ccpn", "hncocacb")
importer.import_3d_peaks(working_dir+r"\peaklists\hncoca_2.tsv", "ccpn", "hncoca")
importer.find_shifts_from_peaks()

a.obs = importer.obs
a.import_pred_shifts(working_dir+r"\shift predictions\alphafold_model_0_H.cs", "shiftx2", offset=0)

a.prepare_obs_preds()
a.calc_log_prob_matrix()
a.calc_mismatch_matrix()

if a.pars["iterate_until_consistent"]:
    a.assign_df = a.find_consistent_assignments(set_assign_df=True)
else:
    a.assign_from_preds(set_assign_df=True)
    # breakpoint()
    a.add_consistency_info(threshold=a.pars["seq_link_threshold"])
    # breakpoint()

# headings = '''
#     Res_name Res_N Res_type SS_name Dummy_res Dummy_SS CA CA_pred HA HA_pred H H_pred CB CB_pred
#         C C_pred N N_pred Log_prob Max_mismatch_m1 Max_mismatch_p1 Num_good_links_m1 
# '''.split()
# table = []
# for df_index, df_row in a.assign_df.iterrows():
#     table_row = []
#     table.append(table_row)
#     for heading in headings:
#         table_row.append(df_row[heading])

# with open(working_dir+r"\SNAPS_output.txt", 'w') as fp:
#     print(tabulate(table, tablefmt='plain', headers=headings), file=fp)



#### Write chemical shift lists
a.assign_df.to_csv(working_dir+r"\output\assign_df_"+output_name+".txt")

a.output_shiftlist(working_dir+r"\output\shiftlist_"+output_name+".txt", "sparky",
                    confidence_list=["High","Medium","Low","Unreliable","Undefined"])

#### Make some plots
plots = []

hsqc_plot = a.plot_hsqc(working_dir+r"\output\hsqc_"+output_name+".htm", "html")
plots += [hsqc_plot]


strip_plot = a.plot_strips(working_dir+r"\output\strip_plot_"+output_name+".htm", "html")
plots += [strip_plot]

# Generate a list of residues n the order determined through sequential matching
matching = a.find_seq_assignment()
matching["obs_N"] = matching["i"].str.replace(r"[a-zA-Z]*", "")
matching.index = matching.i_m1
unassigned = list(matching.index)
seq_assignment = []
next_res = unassigned[0]
while len(unassigned)>0:        # Remember that there can be multiple loops in a sequential assignment!
    seq_assignment += [matching.loc[next_res, "i_m1"]]
    unassigned.remove(next_res)
    if matching.loc[next_res, "i"] in unassigned:
        next_res = matching.loc[next_res, "i"]
    elif len(unassigned)>0:
        next_res = unassigned[0]

matching = matching.loc[seq_assignment, :]  # Sort into order
matching.index = matching.i

seq_df = a.assign_df
seq_df.index = seq_df.SS_name
seq_df = seq_df.loc[matching.index, :]

#Calculate mismatches
seq_df["SS_name_m1"] = seq_df.SS_name.shift(1)
seq_df["SS_name_m1"].iloc[0] = seq_df.SS_name.iloc[-1]
seq_df["SS_name_p1"] = seq_df.SS_name.shift(-1)
seq_df["SS_name_p1"].iloc[-1] = seq_df.SS_name.iloc[1]

seq_df["Max_mismatch_m1"] = pd.Series(list(df_lookup(a.mismatch_matrix,
                        seq_df["SS_name_m1"], seq_df["SS_name"])), index=seq_df.index)
seq_df["Max_mismatch_p1"] = pd.Series(list(df_lookup(a.mismatch_matrix,
                        seq_df["SS_name"], seq_df["SS_name_p1"])), index=seq_df.index)
seq_df["Max_mismatch"] = seq_df[["Max_mismatch_m1","Max_mismatch_p1"]].max(axis=1)

# Make a strip plot of the sequential matching

from bokeh.plotting import figure, output_file, save
from bokeh.layouts import gridplot
from bokeh.models.ranges import Range1d
from bokeh.models import WheelZoomTool, LabelSet, ColumnDataSource, Span

tmp_plt = figure(x_range=seq_df["SS_name"].tolist())

plotlist = []


# Show the confidence and assignments
plt = figure(title=("Confidence plot"
                    "(pan/zoom tools can be accessed at top right; "
                    "mouse over the plot to see observed spin system name)"),
            x_range=tmp_plt.x_range,
            y_range = Range1d(0, 1.5),
            tools="xpan, xwheel_zoom,hover,save,reset",
            tooltips=[("Pred", "@Res_name"),("Obs","@SS_name")],
            height=100, width=1000)

# Create a colour map based on confidence
colourmap = {"High":"green",
                "Medium":"yellowgreen",
                "Low":"orange",
                "Unreliable":"red",
                "Undefined":"grey"}

# Plot the peaks
for k in colourmap.keys():
    tmp = ColumnDataSource(seq_df[seq_df["Confidence"]==k])
    plt.vbar(x="SS_name", top=1, width=1,
                color=colourmap[k], legend_label=k, source=tmp)

# Set legend properties
plt.legend.orientation = "horizontal"
plt.legend.location = "top_center"
plt.legend.padding = 0
plt.legend.margin = 0

# Change axis label orientation
plt.xaxis.major_label_orientation = 3.14159/2
plt.axis.visible = False


plotlist = plotlist + [plt]

#Make the mismatch plot
plt = figure(title="Mismatch plot",
                x_range=tmp_plt.x_range,
                y_axis_label="Mismatch (ppm)",
                tools="xpan, xwheel_zoom,save,reset",
                height=200, width=1000)

plt.vbar(x=seq_df["SS_name"],
            top=seq_df["Max_mismatch"],
            width=1)

# Draw a line showing the threshold for mismatches
# threshold_line = Span(location=0.2,
#                         dimension="width", line_dash="dashed",
#                         line_color="red")
# plt.add_layout(threshold_line)

# Change axis label orientation
plt.xaxis.major_label_orientation = 3.14159/2

plotlist += [plt]

# Combine and output
p = gridplot(plotlist, ncols=1)

output_file(working_dir+r"\output\sequential mismatch.htm")
save(p)