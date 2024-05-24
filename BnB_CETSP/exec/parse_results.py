import os
import glob
from collections import defaultdict

import numpy as np

def convert_to(key,type,dict):
    dict[key]=type(dict[key])
def convert_to_float(key,dict):
    try:
        asfloat=float(dict[key])
    except ValueError:
        asfloat=float.fromhex(dict[key])
    dict[key]=asfloat
def convert_to_percent(key,dict):
    percent=float(dict[key][:-1])
    dict[key]=percent
def parse_log(filepath):
    parsed={}
    with open(filepath) as fh:
        for line in fh:
            if ":" in line:
                key,value=line.split(":")
                chopped=value.split()#remove whitespace
                if len(chopped)==1:
                    parsed[key]=chopped[0]
                else:#all lists are of ints
                    parsed[key]=[int(c) for c in chopped]
    to_bool=["Use LAH",'Use Concorde Feasible Solution improvement']
    for key in to_bool:
        convert_to(key,lambda b:bool(int(b)),parsed)
    to_int=["Instance Size",'Number of Resolved Nodes','Iterations to incumbent','SOCPs solved','Solvers made','SOCP Internal Iterations']
    for key in to_int:
        convert_to(key,int,parsed)
    to_float=['Size of the tree','Initial Upper Bound','Average radius size','Value of objective function','Function objective value','Sum of infeasibilities','Total',
              'Time total','Time S.B','Time SOCP','Warm Start Time','    init Time','    construct Time','        solve Time',
              'SOCP Internal Solve Time','SOCP Setup Time','    SOCP Equilibration Time','    SOCP KKT init Time','SOCP initialization Time','SOCP IP Iteration Time','    SOCP KKT Update Time','    SOCP KKT Solve Time']
    for key in to_float:
        convert_to_float(key,parsed)
    to_percent=['Pruned tree percentage','Computed Tree Percentage']
    for key in to_percent:
        convert_to_percent(key,parsed)
    return parsed

def load_folder(folder):
    logfiles=glob.glob(os.path.join(folder,"*.log"))
    parsed=[parse_log(lf) for lf in logfiles]
    by_size=defaultdict(list)
    for p in parsed:
        by_size[p["Instance Size"]].append(p)
    return by_size

def statistics_by_instance_size(size_to_list_of_dicts,key):
    grouped=dict()
    for s in sorted(size_to_list_of_dicts):
        grouped[s]=[d[key] for d in size_to_list_of_dicts[s]]
    stats=dict()
    stats["means"]={s:np.mean(grouped[s]) for s in grouped}
    stats["stds"]={s:np.std(grouped[s]) for s in grouped}
    stats["minimums"]={s:np.min(grouped[s]) for s in grouped}
    stats["maximums"]={s:np.max(grouped[s]) for s in grouped}
    stats["medians"]={s:np.median(grouped[s]) for s in grouped}
    return stats