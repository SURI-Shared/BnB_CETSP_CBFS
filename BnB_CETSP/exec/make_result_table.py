from exec import parse_results
import os
import glob
from collections import deque
import pandas

parent_folder="/home/ggutow/eclipse-workspace/BnB_CETSP_CBFS/BnB_CETSP/Results/16GB"

data=dict()

queue_of_paths=deque([(parent_folder,)])

while len(queue_of_paths)>0:
    current_path=queue_of_paths.popleft()
    contents=glob.glob(os.path.join(*current_path,"*"))
    child_paths=[current_path+(os.path.split(item)[1],) for item in contents if item!=os.path.join(parent_folder,"table.html") and item!=os.path.join(parent_folder,"table.pdf")]
    if os.path.isfile(os.path.join(*child_paths[0])):
        #reached bottom
        print(f"Loading {os.path.join(*current_path)}")
        folder_data=parse_results.load_folder_by_instance_name(os.path.join(*current_path))
        #traverse dict to this place
        current_dict=data
        for level in current_path[:-1]:
            if level not in current_dict:
                current_dict[level]=dict()
            current_dict=current_dict[level]
        current_dict[current_path[-1]]=folder_data
    else:
        queue_of_paths.extend(child_paths)

#remove non-branching levels from the tree of dicts
data=data[parent_folder]
dicts_to_process=deque([data])
while len(dicts_to_process)>0:
    parent=dicts_to_process.popleft()        
    for key in parent:
        while isinstance(parent[key],dict) and len(parent[key].keys())<2:
            parent[key]=list(parent[key].values())[0]
        if isinstance(parent[key],dict):
            dicts_to_process.append(parent[key])

#collapse all but the last two layers (which are instance name and info key respectively)
columns=['Termination Reason','Function objective value','Time total','Time SOCP','Peak Memory Consumption']
def descend_dict(dictionary,iterable_of_keys):
    for key in iterable_of_keys:
        dictionary=dictionary[key]
    return dictionary
tables=dict()
paths=deque((key,) for key in data.keys())
while len(paths)>0:
    current_path=paths.popleft()
    current_item=descend_dict(data,current_path)
    if isinstance(current_item,dict):
        for key in current_item:
            if isinstance(current_item[key],dict) and columns[0] in current_item[key]:
                tables[current_path]=pandas.DataFrame.from_dict(current_item,'index',columns=columns)
                break
            else:
                paths.append(current_path+(key,))

#sort by instance names
for table in tables.values():
    table.sort_index(axis=0,inplace=True)
stylers=[table.style.set_caption(str(os.path.join(*name))) for name, table in tables.items()]

with open(os.path.join(parent_folder,"table.html"),"wt") as fh:
    for styler in stylers:
        fh.write(styler.to_html())
