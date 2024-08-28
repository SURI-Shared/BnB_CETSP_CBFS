import os
import glob
from collections import defaultdict

import numpy as np
from matplotlib import pyplot
import pickle

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
              'SOCP Internal Solve Time','SOCP Setup Time','    SOCP Equilibration Time','    SOCP KKT init Time','SOCP initialization Time','SOCP IP Iteration Time','SOCP Iteration Time','    SOCP KKT Update Time','    SOCP KKT Solve Time']
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

def load_folder_by_instance_name_recursive(folder):
    logfiles=glob.glob(os.path.join(folder,"**/*.log"),recursive=True)
    parsed=[parse_log(lf) for lf in logfiles]
    by_name=dict()
    for i,p in enumerate(parsed):
        folder=os.path.dirname(logfiles[i])
        by_name[os.path.join(folder,p["Name"])]=p
    return by_name

def human_input_lbs(instance_names,save_path,existing_optimal=None,existing_best_alg=None,existing_best_lb=None):
    if existing_optimal is not None:
        oldoptimal=existing_optimal
        oldbestalg=existing_best_alg
        oldbestlb=existing_best_lb
    else:
        oldoptimal=dict()
        oldbestalg=dict()
        oldbestlb=dict()
    try:
        for path in sorted(instance_names):
            if path in oldoptimal:
                continue
            ok=False
            while not ok:
                wasoptimal=input(f"Was {path} optimal: ")
                try:
                    oldoptimal[path]=bool(int(wasoptimal))
                    ok=True
                except:
                    ok=False
                    print("Try again.")
            if bool(int(wasoptimal)):
                continue
            alg=input(f"Best alg for {path}: ")
            oldbestalg[path]=alg
            ok=False
            while not ok:
                try:
                    lb=input(f"Best lb for {path}: ")
                    oldbestlb[path]=float(lb)
                    ok=True
                except:
                    ok=False
                    print("Try again.")
    except:
        pass
    finally:
        with open(save_path,"wb") as fh:
            pickle.dump({"oldoptimal":oldoptimal,"oldbestalg":oldbestalg,"oldbestlb":oldbestlb},fh)
    return {"oldoptimal":oldoptimal,"oldbestalg":oldbestalg,"oldbestlb":oldbestlb}

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

def bar_plot_avg_vs_size(size_data_dict,key,label,bar_index=0,num_bars=1,ax=None,show_range=False,log=False):
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
        y.append(np.mean(v))
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

def bar_plot_multiple_keys_avgs_vs_size(size_data_dict,keys,ax=None,show_range=False,log=False):
    if ax is None:
        fig=pyplot.figure()
        ax=fig.gca()
    for i,key in enumerate(keys):
        ax=bar_plot_avg_vs_size(size_data_dict,key,key,i,len(keys),ax,show_range,log)
    return ax

def bar_plot_multiple_keys_shms_vs_size(size_data_dict,keys,shift=0,ax=None,show_range=False,log=False):
    if ax is None:
        fig=pyplot.figure()
        ax=fig.gca()
    for i,key in enumerate(keys):
        ax=bar_plot_shm_vs_size(size_data_dict,key,shift,key,i,len(keys),ax,show_range,log)
    return ax

def bar_plot_compare_multiple_keys_avgs_vs_size(size_data_dicts,case_names,keys,ax=None,show_range=False,log=False):
    if ax is None:
        fig=pyplot.figure()
        ax=fig.gca()

    ncases=len(size_data_dicts)
    sizes=list(sorted(size_data_dicts[0].keys()))

    x=np.arange(len(sizes))
    num_bars=ncases*len(keys)
    bar_width=1/(num_bars+1)

    bar_index=0
    for key in keys:
        for i,case in enumerate(size_data_dicts):
            y=[]
            ymin=[]
            ymax=[]
            for size in sizes:
                v=[]
                for data in case[size].values():
                    v.append(data[key])
                v=np.array(v)
                y.append(np.mean(v))
                ymin.append(min(v))
                ymax.append(max(v))

            if not show_range:
                yerr=None
            else:
                yerr=[ymin,ymax]
            case_name=case_names[i]
            ax.bar(x+bar_width*bar_index,y,bar_width,yerr=yerr,label=key+" "+case_name,log=log)
            bar_index+=1
    ax.set_xticks(x+bar_width*(num_bars-1)/2,sizes)
    ax.set_xlabel("Number of Neighborhoods")
    ax.legend()
    return ax

def bar_plot_compare_multiple_keys_shms_vs_size(size_data_dicts,case_names,keys,shift=0,ax=None,show_range=False,log=False):
    if ax is None:
        fig=pyplot.figure()
        ax=fig.gca()

    ncases=len(size_data_dicts)
    sizes=list(sorted(size_data_dicts[0].keys()))

    x=np.arange(len(sizes))
    num_bars=ncases*len(keys)
    bar_width=1/(num_bars+1)

    bar_index=0
    for key in keys:
        for i,case in enumerate(size_data_dicts):
            y=[]
            ymin=[]
            ymax=[]
            for size in sizes:
                v=[]
                for data in case[size].values():
                    v.append(data[key])
                v=np.array(v)
                y.append(shifted_geometric_mean(v,shift))
                ymin.append(min(v))
                ymax.append(max(v))

            if not show_range:
                yerr=None
            else:
                yerr=[ymin,ymax]
            case_name=case_names[i]
            ax.bar(x+bar_width*bar_index,y,bar_width,yerr=yerr,label=key+" "+case_name,log=log)
            bar_index+=1
    ax.set_xticks(x+bar_width*(num_bars-1)/2,sizes)
    ax.set_xlabel("Number of Neighborhoods")
    ax.legend()
    return ax

def compare_bar_plot_stacked_keys_shm_vs_size(size_data_dicts,case_names,keys_for_each_case,shift=0,ax=None,show_range=False,log=False,scale_factors_by_key=None):
    if ax is None:
        fig=pyplot.figure()
        ax=fig.gca()
    
    num_bars=len(size_data_dicts)

    for i in range(len(size_data_dicts)):
        ax=bar_plot_stacked_keys_shm_vs_size(size_data_dicts[i],case_names[i],keys_for_each_case[i],num_bars,i,shift,ax,show_range,log,scale_factors_by_key)
    return ax

def bar_plot_stacked_keys_shm_vs_size(size_data_dict,case_name,keys,num_bars,bar_index,shift=0,ax=None,show_range=False,log=False,scale_factors_by_key=None):
    if ax is None:
        fig=pyplot.figure()
        ax=fig.gca()

    sizes=list(sorted(size_data_dict.keys()))

    x=np.arange(len(sizes))
    bar_width=1/(num_bars+1)

    bottom=np.zeros(len(sizes))
    for key in keys:
        y=[]
        ymin=[]
        ymax=[]
        for size in sizes:
            v=[]
            for data in size_data_dict[size].values():
                v.append(data[key])
            v=np.array(v)
            if scale_factors_by_key is not None and key in scale_factors_by_key:
                v*=scale_factors_by_key[key]
            y.append(shifted_geometric_mean(v,shift))
            ymin.append(min(v))
            ymax.append(max(v))

        if not show_range:
            yerr=None
        else:
            yerr=[ymin,ymax]
        ax.bar(x+bar_width*bar_index,y,bar_width,bottom=bottom,yerr=yerr,label=key+" "+case_name,log=log)
        bottom+=y

    ax.set_xticks(x+bar_width*(num_bars-1)/2,sizes)
    ax.set_xlabel("Number of Neighborhoods")
    ax.legend()
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

def bar_plot_reductions_in_shm(nominal_size_data_dict,size_data_dicts,key,shift,labels,ylabel,ax=None,show_range=False,min_size=0):
    if ax is None:
        fig=pyplot.figure()
        ax=fig.gca()
    y=[]
    sizes=[size for size in sorted(nominal_size_data_dict.keys()) if size>=min_size]
    nominal_value_arrays_by_size=dict()
    for size in sizes:
        v=[]
        for data in nominal_size_data_dict[size].values():
            v.append(data[key])
        v=np.array(v)
        shm=shifted_geometric_mean(v,shift)
        y.append(shm)
        nominal_value_arrays_by_size[size]=v
    nominal_shm=np.array(y)

    num_bars=len(size_data_dicts)
    x=np.arange(len(sizes))
    bar_width=1/(num_bars+1)

    for i,size_data_dict in enumerate(size_data_dicts):
        y=[]
        ymin=[]
        ymax=[]
        for size in sizes:
            v=[]
            ratio=[]
            for instance,data in size_data_dict[size].items():
                v.append(data[key])
                nominal=nominal_size_data_dict[size][instance][key]
                ratio.append(data[key]/nominal)
            v=np.array(v)
            shm=shifted_geometric_mean(v,shift)
            y.append(shm)
            ymin.append(min(ratio))
            ymax.append(max(ratio))
        normalized_shm=np.array(y)/nominal_shm
        normalized_min=np.array(ymin)
        normalized_max=np.array(ymax)
        if not show_range:
            yerr=None
        else:
            yerr=[normalized_min,normalized_max]
        ax.bar(x+bar_width*i,normalized_shm,bar_width,yerr=yerr,label=labels[i])
    ax.set_xticks(x+bar_width*(num_bars-1)/2,sizes)
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

def scatter_plot_vs_baseline_value(nominal_size_data_dict,size_data_dicts,key,labels,xlabel,ylabel,ax=None):
    if ax is None:
        fig=pyplot.figure()
        ax=fig.gca()
    sizes=[size for size in sorted(nominal_size_data_dict.keys())]
    nominal_values=np.array([entry[key] for size in sizes for entry in nominal_size_data_dict[size].values()])
    sorted_indices=np.argsort(nominal_values)
    sorted_nominals=nominal_values[sorted_indices]

    for i,size_data_dict in enumerate(size_data_dicts):
        y=np.array([entry[key] for size in sizes for entry in size_data_dict[size].values()])

        ax.scatter(sorted_nominals,y[sorted_indices],label=labels[i])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend()
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