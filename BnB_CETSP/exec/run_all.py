import os
import glob
import git
import time
import subprocess
import multiprocessing
import sys
from itertools import zip_longest,repeat

def write_README(outfolder,**kwargs):
    git_repo_path=os.path.join(os.path.dirname(__file__),"..","..")
    repo=git.Repo(os.path.abspath(git_repo_path))
    commit_name=repo.head.commit.name_rev
    with open(os.path.join(outfolder,"README"),"w") as fh:
        fh.write(f"commit: {commit_name}\n")
        for key,val in kwargs.items():
            fh.write(f"{key}: {val}\n")

def runner(infile,executable,outfolder,args):
    infile_name,ext=os.path.splitext(os.path.basename(infile))
    print(f"Start {infile}")
    log_file_name=os.path.join(outfolder,infile_name+".log")
    with open(log_file_name,"wt") as fh:
        subprocess.run([executable,infile]+args,stdout=fh,stderr=fh)
    print(f"Done {infile}")

def test_planner(processes,executable,in_folders,out_prefixs,args_lists):
    timestr=time.strftime("%Y%m%d_%H%M%S")
    instance_files=[]
    save_paths=[]
    args=[]
    for in_folder,out_prefix,arg_list in zip(in_folders,out_prefixs,args_lists):
        outfolder=os.path.join(out_prefix,executable+"_"+timestr)
        os.makedirs(outfolder,exist_ok=True)
        write_README(outfolder,in_folder=in_folder,executable=executable,cmdline_args=args)

        infiles=glob.glob(os.path.join(in_folder,"*.txt"))
        instance_files.extend(infiles)
        save_paths.extend(repeat(outfolder,len(infiles)))
        args.extend(repeat(arg_list,len(infiles)))
    with multiprocessing.Pool(processes) as pool:
        pool.starmap(runner,zip(instance_files,repeat(executable),save_paths,args))

if __name__=="__main__":
    nprocesses=int(sys.argv[1])
    executable=sys.argv[2]
    save_prefix=sys.argv[3]
    in_folders=["2D","3D","3D","3D","RND","Behdani","Behdani","Behdani"]
    out_folder_names=["2D","3D/OR0.02","3D/OR0.1","3D/OR0.3","RND","Behdani/R0.25","Behdani/R0.5","Behdani/R1"]
    out_prefixes=[os.path.join(save_prefix,ofn) for ofn in out_folder_names]
    make_arg_list=lambda dimension,overlap_ratio: [dimension, overlap_ratio, "14400", "V1", "BeFS", "1", "1"]
    args_lists=[make_arg_list(dim,overlap_ratio) for dim,overlap_ratio in [("2D","1"),("3D","0.02"),("3D","0.1"),("3D","0.3"),("3D","1"),("2D","0.25"),("2D","0.5"),("2D","1")]]
    test_planner(nprocesses,executable,in_folders,out_prefixes,args_lists)