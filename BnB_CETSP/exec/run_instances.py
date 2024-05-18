import os
import glob
import git
import time
import subprocess
import sys

def write_README(outfolder,**kwargs):
    git_repo_path=os.path.join(os.path.dirname(__file__),"..","..")
    repo=git.Repo(os.path.abspath(git_repo_path))
    commit_name=repo.head.commit.name_rev
    with open(os.path.join(outfolder,"README"),"w") as fh:
        fh.write(f"commit: {commit_name}\n")
        for key,val in kwargs.items():
            fh.write(f"{key}: {val}\n")

def runner(executable, infile_name,log_file_name,args):
    print(f"Start {infile_name}")
    with open(log_file_name,"wt") as fh:
        subprocess.run([executable,infile_name]+args,stdout=fh,stderr=fh)
    print(f"Done {infile_name}")
def test_planner(processes,executable,in_folder,out_prefix,args):
    timestr=time.strftime("%Y%m%d_%H%M%S")
    outfolder=os.path.join(out_prefix,executable+"_"+timestr)
    os.makedirs(outfolder,exist_ok=True)
    write_README(outfolder,in_folder=in_folder,executable=executable,cmdline_args=args)

    infiles=glob.glob(os.path.join(in_folder,"*.txt"))
    outfile_handles=[]
    procs=[]
    for infile in infiles:
        print(infile)
        infile_name,ext=os.path.splitext(os.path.basename(infile))
        log_file_name=os.path.join(outfolder,infile_name+".log")
        fh=open(log_file_name,"wt")
        outfile_handles.append(fh)
        procs.append(subprocess.Popen([executable,infile]+args,stdout=fh,stderr=fh))
    for proc in procs:
        proc.wait()
    for fh in outfile_handles:
        fh.close()

if __name__=="__main__":
    executable=sys.argv[1]
    in_folder=sys.argv[2]
    out_prefix=sys.argv[3]
    args=sys.argv[4:]
    test_planner(executable,in_folder,out_prefix,args)