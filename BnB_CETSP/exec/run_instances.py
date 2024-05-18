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

def test_planner(processes,executable,in_folder,out_prefix,args):
    timestr=time.strftime("%Y%m%d_%H%M%S")
    outfolder=os.path.join(out_prefix,executable+"_"+timestr)
    os.makedirs(outfolder,exist_ok=True)
    write_README(outfolder,in_folder=in_folder,executable=executable,cmdline_args=args)

    infiles=glob.glob(os.path.join(in_folder,"*.txt"))
    with multiprocessing.Pool(processes) as pool:
        pool.starmap(runner,zip(infiles,repeat(executable),repeat(outfolder),repeat(args)))

if __name__=="__main__":
    nprocesses=int(sys.argv[1])
    executable=sys.argv[2]
    in_folder=sys.argv[3]
    out_prefix=sys.argv[4]
    args=sys.argv[5:]
    test_planner(nprocesses,executable,in_folder,out_prefix,args)