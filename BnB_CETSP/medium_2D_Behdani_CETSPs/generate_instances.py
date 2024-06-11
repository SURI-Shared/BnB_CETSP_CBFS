import numpy as np

class CETSPFormatError(Exception):
    pass

def parse_cetsp_file(path,is_2D=False,radius=None):
    """
    read a .cetsp file formatted as in https://github.com/CerroneCarmine/CETSP/tree/master
    also handles 2D instances as in https://github.com/UranusR/BnB_CETSP_CBFS/tree/main/BnB_CETSP/2D

    Parameters: path : path-like
                is_2D : bool
                    if True, file is specifying a 2D instance and does not include radii or demand information. radius must be specified
                radius : float
                    if is_2D, every entry of r will be radius
    Returns:    waypoints : np.ndarray
                    n,2 or n,3 positions of waypoints
                r : list[float]
                    radii of neighborhoods
                demand : list[float]
                    demand values for each waypoint (all 0 if is_2D)
                depot : list[float]
                    coordinates of the starting point (x,y,z) or (x,y)
    """
    comment_marker_length=2
    depot_prefixes=["//Depot is ","//Depot: "]
    x=[]
    y=[]
    z=[]
    r=[]
    demand=[]
    depot=[None,None,None]
    with open(path,"r",encoding="UTF-8") as fh:
        for line in fh:
            if line[0:comment_marker_length]=="//":
                #comment
                str_pos=None
                for prefix in depot_prefixes:
                    if line[:len(prefix)]==prefix:
                        #depot line
                        str_pos=line[len(prefix):].split(",")
                        if str_pos[2][-1]!="\n":
                            raise CETSPFormatError(f"Depot specification line beginning '{prefix}' must be followed by three numbers separated by ','.")
                        else:
                            break
                if str_pos is not None:
                    depot[0]=float(str_pos[0])
                    depot[1]=float(str_pos[1])
                    depot[2]=float(str_pos[2][:-1])#drop newline at end
            else:
                str_values=line.split()
                if len(str_values)==5:
                    x.append(float(str_values[0]))
                    y.append(float(str_values[1]))
                    z.append(float(str_values[2]))
                    r.append(float(str_values[3]))
                    demand.append(float(str_values[4]))
                elif is_2D and len(str_values)==2:
                    #the 2D instances only have x,y coordinates, radius is manually chosen
                    x.append(float(str_values[0]))
                    y.append(float(str_values[1]))
                elif len(str_values)!=0:
                    raise CETSPFormatError("Non-comment lines should contain five entries separated by whitespace")
    if is_2D:
        #depot is the first row
        depot=[x[0],y[0],0]
        x=x[1:]
        y=y[1:]
        demand=[0]*len(x)
        r=[radius]*len(x)
    if any((d is None for d in depot)):
        #missing depot specification line
        raise CETSPFormatError("File must contain a line beginning '//Depot is ' specifying the depot location")
    if is_2D:
        waypoints= np.stack([x,y],axis=1)
    else:
        waypoints= np.stack([x,y,z],axis=1)
    return waypoints,np.array(r),np.array(demand),np.array(depot)

class CETSPInstance:
    np_array_names=["depot","waypoints","r","demand"]
    def __init__(self,waypoints,r,demand,depot):
        self.depot=depot
        self.waypoints=waypoints
        self.r=r
        self.demand=demand
    @classmethod
    def from_cetsp_file(cls,path,is_2D=False,radius=None,drop_z=False):
        waypoints,r,demand,depot=parse_cetsp_file(path,is_2D,radius)
        if drop_z:
            waypoints=waypoints[:,:2]
            depot=depot[:2]
        return cls(waypoints,r,demand,depot)
    def to_cetsp_file(self,path,as2D=False):
        with open(path,"wt") as fh:
            self.write_header(fh,as2D)
            if as2D:
                #first line needs to be the depot
                fh.write(f"{self.depot[0]} {self.depot[1]}\n")
            wp_lines=[]
            for i in range(len(self.waypoints)):
                if not as2D:
                    wp_lines.append(f"{self.waypoints[i][0]} {self.waypoints[i][1]} {self.waypoints[i][2]} {self.r[i]} {self.demand[i]}\n")
                else:
                    wp_lines.append(f"{self.waypoints[i][0]} {self.waypoints[i][1]}\n")
            fh.writelines(wp_lines)
            self.write_footer(fh,as2D)
    def write_header(self,fh,as2D=False):
        if not as2D:
            fh.writelines(["//Column order: x, y, z, radius, node demand (for the Close-Enough VRP)\n",
            "//The depot's (x,y,z) is given at file's end\n","// description\n",""])
    def write_footer(self,fh,as2D=False):
        if not as2D:
            fh.writelines(["\n",
            f"//Depot: {self.depot[0]}, {self.depot[1]}, {self.depot[2]}\n",
            f"//Max demand = {np.max(self.demand)}\n",
            f"//Total demand = {np.sum(self.demand)}\n",
            f"//Avg. demand per node = {np.mean(self.demand)}"])
    def to_npz(self,file_path):
        np.savez(file_path,**{arr_name:self.__dict__[arr_name] for arr_name in self.np_array_names})
    @classmethod
    def from_npz(cls,file_path):
        self=cls.__new__(cls)
        #load numpy array attributes from the npz file
        with np.load(file_path) as data:
            for name in cls.np_array_names:
                self.__dict__[name]=data[name]
        return self
    def overlap_ratio(self):
        """
        Computes overlap ratio following the code released for Zhang et al. Results for the CETSP with a Branch and Bound Algorithm.
        They use the average radius of the neighborhoods divided by the largest value of a coordinate of a waypoint or the depot.
        """
        rmean=np.mean(self.r)
        max_waypoint=max(np.max(self.depot),np.max(self.waypoints))
        return rmean/max_waypoint
    @classmethod
    def random_behdani(cls,n_neighborhoods,radius,rng):
        """
        Generate a random 2D CETSP instance with all circles of fixed radius, all waypoints and depot inside [0,16]x[0,10]
        
        as described in:
        Behnam Behdani, J. Cole Smith (2014) An Integer-Programming-Based Approach to the Close-Enough Traveling Salesman
        Problem. INFORMS Journal on Computing 26(3):415-432. https://doi.org/10.1287/ijoc.2013.0574
        """
        waypoints=rng.uniform(np.zeros(2),np.array([16,10]),size=(n_neighborhoods,2))
        depot=rng.uniform(np.zeros(2),np.array([16,10]))
        demand=np.zeros(n_neighborhoods)
        return cls(waypoints,np.full(n_neighborhoods,radius),demand,depot)

rng=np.random.default_rng(0)
folder="/home/ggutow/eclipse-workspace/bayesian-ergodic-search/assets/double_integrator3D_TSPN/medium_2D_Behdani_CETSPs/"
number=30
for N in [5,10,15,20,25,30,35]:
    for i in range(number):
        instance=CETSPInstance.random_behdani(N,0.25,rng)
        filename=folder+f"N{N}_{i}.txt"
        instance.to_cetsp_file(filename,as2D=True)
