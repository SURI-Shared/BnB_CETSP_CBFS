from exec import parse_results
import os

#specify the baseline run first, followed by the run to use for each variant
folder=["exeCVXHULL_20240529_125208",
        "clarabel_redundant_20240529_115148",
        "clarabel_dropin_20240528_201056",
        "clarabel_recycling_20240524_192136",
        "clarabel_warmstart_20240530_203854",
        "scs_dropin_20240827_161712",
        "scs_recycling_20240827_161857",
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

ax=parse_results.bar_plot_reductions_in_shm(data[0],data[1:],"Time total",0,names,"Fraction of Cout CPU Time",show_range=True,min_size=10)