from exec import parse_results
import os

#specify the baseline run first, followed by the run to use for each variant
folder=["exeCVXHULL_20240905_122502",
        "clarabel_redundant_20240905_122514",
        "clarabel_reduce_20240905_122520",
        "clarabel_reuse_20240905_142826",
        "clarabel_recycle_20240905_122528",
        "scs_reduce_20240905_122531",
        "scs_reuse_20240905_122536",
        "scs_recycle_20240905_122541",
        "scs_tridiag_reduce_20240905_122545",
        "scs_tridiag_reuse_20240905_122550",
        "scs_tridiag_recycle_20240905_122555"]

#variant names
names=["Clarabel",
       "Clarabel+Reduce",
       "Clarabel+Reduce+Reuse",
       "Clarabel+Reduce+Reuse+Recycle",
       "SCS+Reduce",
       "SCS+Reduce+Reuse",
       "SCS+Reduce+Reuse+Recycle",
       "TSCS+Reduce",
       "TSCS+Reduce+Reuse",
       "TSCS+Reduce+Reuse+Recycle"]
#load data
data=[parse_results.load_folder(os.path.join("Results/medium_2D_Behdani_CETSPs/radius_0.25/run/",f)) for f in folder]

# #convert times to milliseconds

# #all Clarabel times are in seconds
# time_keys=["SOCP Setup Time","SOCP IP Iteration Time","Warm Start Time",'SOCP initialization Time']
# clarabel_indices=list(range(2,5))
# for itemid in clarabel_indices:
#     for size in data[itemid].keys():
#         for case in data[itemid][size]:
#             for key in time_keys:
#                 if key in data[itemid][size][case]:
#                     data[itemid][size][case][key]*=1000

# #SCS Warm Start Time is in seconds, internals times are already milliseconds
# scs_indices=list(range(5,11))
# for itemid in scs_indices:
#     for size in data[itemid].keys():
#         for case in data[itemid][size]:
#             if "Warm Start Time" in data[itemid][size][case]:
#                 data[itemid][size][case]["Warm Start Time"]*=1000


ax=parse_results.compare_bar_plot_stacked_keys_avg_ratio_all_sizes(data[1],"Time total",[data[2],data[3],data[4],
                                                            data[5],data[6],data[7],
                                                            data[8],data[9],data[10]],
                                                           ["CRB+Red","CRB+Red+Reu","CRB+Red+Reu+Rec",
                                                            "SCS+Red","SCS+Red+Reu","SCS+Red+Reu+Rec",
                                                            "TSCS+Red","TSCS+Red+Reu","TSCS+Red+Reu+Rec"],
                                                           [["SOCP Setup Time", 'SOCP initialization Time',"SOCP Iteration Time"],
                                                            ["SOCP Setup Time",'SOCP initialization Time',"SOCP Iteration Time"],
                                                            ["SOCP Setup Time",'SOCP initialization Time',"SOCP Iteration Time","Warm Start Time"]]+
                                                         2*[["SOCP Setup Time", "SOCP Iteration Time"],
                                                            ["SOCP Setup Time","SOCP Iteration Time"],
                                                            ["SOCP Setup Time","SOCP Iteration Time","Warm Start Time"]],"Time total","Other",
                                                            {"SOCP Setup Time":"r","SOCP Iteration Time":"b","Warm Start Time":"g",'SOCP initialization Time':'y',"Other":'m'},{"SOCP initialization Time":"edge"})
ax.set_xlabel("Avg Fraction of Cout+Clarabel Wall Time",fontsize='xx-large')