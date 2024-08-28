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

ax=parse_results.compare_bar_plot_stacked_keys_shm_vs_size([data[5],data[6],data[7],
                                                            data[8],data[9],data[10]],
                                                           ["SCS+Reduce","SCS+Reduce+Reuse","SCS+Reduce+Reuse+Recycle",
                                                            "TSCS+Reduce","TSCS+Reduce+Reuse","TSCS+Reduce+Reuse+Recycle"],
                                                         2*[["SOCP Setup Time", "SOCP Iteration Time"],
                                                            ["SOCP Setup Time","SOCP Iteration Time"],
                                                            ["SOCP Setup Time","SOCP Iteration Time","Warm Start Time"]],
                                                            {"SOCP Setup Time":"r","SOCP Iteration Time":"g","Warm Start Time":"b"},
                                                            scale_factors_by_key={"Warm Start Time":1000},log=False)
ax.set_xlabel("Geometric Mean CPU Time (ms)")