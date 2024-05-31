import os
import glob
from collections import defaultdict

import numpy as np
from matplotlib import pyplot

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
        if key in parsed:
            convert_to(key,lambda b:bool(int(b)),parsed)
    to_int=["Instance Size",'Number of Resolved Nodes','Iterations to incumbent','SOCPs solved','Solvers made','SOCP Internal Iterations']
    for key in to_int:
        if key in parsed:
            convert_to(key,int,parsed)
    to_float=['Size of the tree','Initial Upper Bound','Average radius size','Value of objective function','Function objective value','Sum of infeasibilities','Total',
              'Time total','Time S.B','Time SOCP','Warm Start Time','    init Time','    construct Time','        solve Time',
              'SOCP Internal Solve Time','SOCP Setup Time','    SOCP Equilibration Time','    SOCP KKT init Time','SOCP initialization Time','SOCP IP Iteration Time','    SOCP KKT Update Time','    SOCP KKT Solve Time']
    for key in to_float:
        if key in parsed:
            convert_to_float(key,parsed)
    to_percent=['Pruned tree percentage','Computed Tree Percentage']
    for key in to_percent:
        if key in parsed:
            convert_to_percent(key,parsed)
    return parsed

def load_folder(folder):
    logfiles=glob.glob(os.path.join(folder,"*.log"))
    parsed=[parse_log(lf) for lf in logfiles]
    by_size=defaultdict(dict)
    for p in parsed:
        by_size[p["Instance Size"]][p["Name"]]=p
    return by_size
def load_folder_by_instance_name(folder):
    logfiles=glob.glob(os.path.join(folder,"*.log"))
    parsed=[parse_log(lf) for lf in logfiles]
    by_name=dict()
    for p in parsed:
        by_name[p["Name"]]=p
    return by_name

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

def shifted_geometric_mean(values,shift=1):
    N=len(values)
    return np.prod(values+shift)**(1/N)-shift

def plot_shm_vs_size(size_data_dict,key,shift,label,ax=None,linespec=None,show_range=False):
    if ax is None:
        fig=pyplot.figure()
        ax=fig.gca()
    y=[]
    ymin=[]
    ymax=[]
    x=list(sorted(size_data_dict.keys()))
    for size in x:
        v=[]
        for data in size_data_dict[size].values():
            v.append(data[key])
        v=np.array(v)
        shm=shifted_geometric_mean(v,shift)
        y.append(shm)
        ymin.append(min(v))
        ymax.append(max(v))
    if linespec is not None:
        if not show_range:
            ax.plot(x,y,linespec,label=label)
        else:
            ax.errorbar(x,y,[ymin,ymax],linespec,label=label)
    else:
        if not show_range:
            ax.plot(x,y,label=label)
        else:
            ax.errorbar(x,y,[ymin,ymax],label=label)
    ax.set_xlabel("Number of Neighborhoods")
    return ax

def bar_plot_shm_vs_size(size_data_dict,key,shift,label,bar_index=0,num_bars=1,ax=None,show_range=False,log=False):
    if ax is None:
        fig=pyplot.figure()
        ax=fig.gca()
    y=[]
    ymin=[]
    ymax=[]
    sizes=list(sorted(size_data_dict.keys()))

    for size in sizes:
        v=[]
        for data in size_data_dict[size].values():
            v.append(data[key])
        v=np.array(v)
        shm=shifted_geometric_mean(v,shift)
        y.append(shm)
        ymin.append(min(v))
        ymax.append(max(v))

    x=np.arange(len(sizes))
    bar_width=1/(num_bars+1)
    if not show_range:
        yerr=None
    else:
        yerr=[ymin,ymax]
    ax.bar(x+bar_width*bar_index,y,bar_width,yerr=yerr,label=label,log=log)
    ax.set_xticks(x+bar_width*(num_bars-1)/2,sizes)
    ax.set_xlabel("Number of Neighborhoods")
    return ax

def plot_reductions_in_shm(nominal_size_data_dict,size_data_dicts,key,shift,labels,ylabel,ax=None):
    if ax is None:
        fig=pyplot.figure()
        ax=fig.gca()
    y=[]
    x=list(sorted(nominal_size_data_dict.keys()))
    for size in x:
        v=[]
        for data in nominal_size_data_dict[size].values():
            v.append(data[key])
        v=np.array(v)
        shm=shifted_geometric_mean(v,shift)
        y.append(shm)
    nominal_shm=np.array(y)

    for i,size_data_dict in enumerate(size_data_dicts):
        y=[]
        for size in x:
            v=[]
            for data in size_data_dict[size].values():
                v.append(data[key])
            v=np.array(v)
            shm=shifted_geometric_mean(v,shift)
            y.append(shm)
        normalized_shm=np.array(y)/nominal_shm
        ax.plot(x,normalized_shm,label=labels[i])
    ax.set_xlabel("Number of Neighborhoods")
    ax.set_ylabel(ylabel)
    ax.legend()
    return ax
def bar_plot_vs_size(size_value_dict,bar_index=0,num_bars=1,ax=None,log=False,ylabel=None,label=None):
    if ax is None:
        fig=pyplot.figure()
        ax=fig.gca()
    x=np.arange(len(size_value_dict))
    bar_width=1/(num_bars+1)
    sizes=list(sorted(size_value_dict.keys()))
    ax.bar(x+bar_width*bar_index,[size_value_dict[s] for s in sizes],bar_width,log=log,label=label)
    ax.set_xticks(x+bar_width*(num_bars-1)/2,sizes)
    ax.set_xlabel("Number of Neighborhoods")
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    return ax

def scatter_plot_percent_delta_vs_size(size_value_nominal,size_value_modified,key,ax=None,ylabel=None,tol=1e-6):
    if ax is None:
        fig=pyplot.figure()
        ax=fig.gca()
    y=[]
    x=[]
    for size in size_value_nominal:
        for instance in size_value_nominal[size]:
            nominal=size_value_nominal[size][instance][key]
            modified=size_value_modified[size][instance][key]
            if abs(nominal)>tol:
                percent_delta=(modified-nominal)/nominal*100
                y.append(percent_delta)
                x.append(size)
    ax.scatter(x,y)
    ax.grid(True,"major","y")
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    ax.set_xlabel("Number of Neighborhoods")
    return ax

def scatter_plot_case1_vs_case2(size_value1,size_value2,key,xlabel,ylabel=None,ax=None):
    if ax is None:
        fig=pyplot.figure()
        ax=fig.gca()
    y=[]
    x=[]
    for size in size_value1:
        for instance in size_value1[size]:
            nominal=size_value1[size][instance][key]
            modified=size_value2[size][instance][key]
            x.append(nominal)
            y.append(modified)
    ax.scatter(x,y)
    ax.plot([min(x+y),max(x+y)],[min(x+y),max(x+y)],":")
    ax.grid(True,"major","both")
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.axis("equal")
    return ax