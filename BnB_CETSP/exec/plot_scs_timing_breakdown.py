from exec import parse_results
import os

#specify the baseline run first, followed by the run to use for each variant
folder=["exeCVXHULL_20240827_164534",
        "clarabel_redundant_20240827_164556",
        "clarabel_dropin_20240827_164609",
        "clarabel_recycling_20240827_164619",
        "clarabel_warmstart_20240827_164633",
        "scs_dropin_20240827_161712",
        "scs_recycling_20240827_191125",
        "scs_warmstart_20240827_161844",
        "scs_tridiag_dropin_20240827_161750",
        "scs_tridiag_recycling_20240827_161813",
        "scs_tridiag_warmstart_20240827_161827"]

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

#convert times to milliseconds

#all Clarabel times are in seconds
time_keys=["SOCP Setup Time","SOCP IP Iteration Time","Warm Start Time",'SOCP initialization Time']
clarabel_indices=list(range(2,5))
for itemid in clarabel_indices:
    for size in data[itemid].keys():
        for case in data[itemid][size]:
            for key in time_keys:
                if key in data[itemid][size][case]:
                    data[itemid][size][case][key]*=1000

#SCS Warm Start Time is in seconds, internals times are already milliseconds
scs_indices=list(range(5,11))
for itemid in scs_indices:
    for size in data[itemid].keys():
        for case in data[itemid][size]:
            if "Warm Start Time" in data[itemid][size][case]:
                data[itemid][size][case]["Warm Start Time"]*=1000


ax=parse_results.compare_bar_plot_stacked_keys_shm_vs_size([data[2],data[3],data[4],
                                                            data[5],data[6],data[7],
                                                            data[8],data[9],data[10]],
                                                           ["Clarabel+Reduce","Clarabel+Reduce+Reuse","Clarabel+Reduce+Reuse+Recycle",
                                                            "SCS+Reduce","SCS+Reduce+Reuse","SCS+Reduce+Reuse+Recycle",
                                                            "TSCS+Reduce","TSCS+Reduce+Reuse","TSCS+Reduce+Reuse+Recycle"],
                                                           [["SOCP Setup Time", 'SOCP initialization Time',"SOCP IP Iteration Time"],
                                                            ["SOCP Setup Time",'SOCP initialization Time',"SOCP IP Iteration Time"],
                                                            ["SOCP Setup Time",'SOCP initialization Time',"SOCP IP Iteration Time","Warm Start Time"]]+
                                                         2*[["SOCP Setup Time", "SOCP Iteration Time"],
                                                            ["SOCP Setup Time","SOCP Iteration Time"],
                                                            ["SOCP Setup Time","SOCP Iteration Time","Warm Start Time"]],
                                                            {"SOCP Setup Time":"r","SOCP Iteration Time":"g","SOCP IP Iteration Time":"g","Warm Start Time":"b",'SOCP initialization Time':'y'},
                                                            log=False)
ax.set_xlabel("Geometric Mean CPU Time (ms)")