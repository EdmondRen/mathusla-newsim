# --------------------------------------------------
# User IO for tracker
# --------------------------------------------------
## The following functions needs to be defined:
## 1. load(filename, *args, **kwargs):
##      """
##      Load the data into a dictionary. 
##
##      INPUT:
##      ---
##      filename: str
##          the filename of the input file
##      *args, **kwargs: 
##          additional optional arguments
##
##      Return:
##      ---
##          data: dict
##              The key of the dictionary is the group of hits.
##              For example, there could be three groups, "inactive", "1" and "2"
##              Group "inactive" will not be processed, and can be used to store floor hits
##              Group "1" and "2" can be used to save hits in top layers and wall tracking layers
##          metadata: dict
##              This dictionary contains metadata. The necessary metadata is the direction of each group
##              Because the tracker takes "y" axis as the default layer direction, it is needed to rotate the coordinates
##              if the layer is actually facing another direction.
##              So, you should provide the coordinate rotation for each group
##              For example: metadata["groups"]["1"] = {"flip_index":[1,2]}
##      """


import math, os, sys, types
import importlib
from collections import namedtuple

# import joblib
import ROOT
import numpy as np

import tracker.datatypes as datatypes

noise_level_simulation = 10 # The simulation noise rate is ten times cosmic
# ------------------------------------------------
# Set the active layers, noise level, detector efficiency
NOISE       = 0.1   # 0.1 times cosmic rate (NOISE*noise_level_simulation=1)
LAYERS      = 6     # Active layers. The code below will automatically decide which layers to keep
EFFICIENCY  = 0.95     # Detector efficiency
# ------------------------------------------------


# Layer ID for each group
# "inactive" group will be ignored in reconstruction
layer_groups = {\
    4:{"inactive": [-1, 0,1,2,3], # All discarded hits will be assigned layer_id of -1
       "1": [6,7,8,9],
       "2": [10,11,12,13]},
    5:{"inactive": [-1, 0,1,2,3],
       "1": [5,6,7,8,9],
       "2": [10,11,12,13,14]},   
    6:{"inactive": [-1,0,1,2,3],
       "1": [4,5,6,7,8,9],
       "2": [10,11,12,13,14,15]},        
    }
layer_group = layer_groups[LAYERS]
rng = np.random.default_rng()
noise_prob = NOISE/noise_level_simulation
assert((noise_prob>=0) & (noise_prob<=1))
efficiency_prob = EFFICIENCY
assert((efficiency_prob>=0) & (efficiency_prob<=1))


# Test:
# import importlib
# io_user = importlib.machinery.SourceFileLoader("*", "io_mu-simulation.py").load_module()
# io_user.load("/project/rrg-mdiamond/tomren/mudata/LLP_geometry_40m///SimOutput/MATHUSLA_LLPfiles_HXX/deweighted_LLP_bb_35/20240410/173948/reconstruction_v2.ROOT")  

# class datatypes:
#     Hit = namedtuple("Hit", ["x", "y", "z", "t", "x_err", "y_err", "z_err", "t_err","layer", "ind"])
#     Track = namedtuple("Track", ["x0", "y0", "z0", "t0", "Ax", "Ay", "Az", "At", "cov", "chi2", "ind", "hits", "hits_filtered"])
#     Vertex = namedtuple("Vertex", ["x0", "y0", "z0", "t0", "cov", "chi2", "tracks"])
#     VertexSeed = namedtuple("VertexSeed",["x0", "y0", "z0", "t0",  "cov", "chi2","dist", "Ntracks", "trackind1","trackind2", "score"])



det_width        = 3.5 # 4.5cm per bar
det_height       = 1 # [cm]
time_resolution  = 1 # [ns], single channel
refraction_index = 1.58

unc_trans = det_width/math.sqrt(12)                  
unc_long = time_resolution*29.979/math.sqrt(2)/refraction_index
unc_time = time_resolution/math.sqrt(2) # ns
unc_thick = det_height/math.sqrt(12) # uncertainty in thickness, cm


class tfile:
    """
    Helper calss for ROOT TFile.
    It wraps the mostly used funcitons to open and read objects from ROOT file.

    
    file = ROOT.TFile.Open(fname_digi)
    tree_name = file.GetListOfKeys()[1].GetName()
    Tree = file.Get(tree_name)
    branches = [Tree.GetListOfBranches()[i].GetName() for i in range(len(Tree.GetListOfBranches()))]
    entries = Tree.GetEntries()    
    
    """
    def __init__(self, filename, verbose=0):
        self.file = ROOT.TFile.Open(filename)
        self._keys = self.file.GetListOfKeys()

    @property
    def keys(self):
        return self._keys        

    def ls(self):
        print(self.file.ls())

    def get(self, key):
        self.object = self.file.Get(key)
        return self.object

    def vector2list(self, vector):
        return [vector.at(i) for i in range(vector.size())]    

    def get_tree(self, key):
        self.tree = self.file.Get(key)
        self.entries = self.tree.GetEntries()
        self.current_entry = -1

        # Read the branch name and type
        # Make a dictionary for each branch
        self.branches = {}
        self.branches_disabled = []
        self.data = {}
        for branch in self.tree.GetListOfBranches():
            # Get the leaf for the branch
            leaf = branch.GetLeaf(branch.GetName())
            # Get the data type of the leaf
            data_type = leaf.GetTypeName()
            
            self.branches[branch.GetName()] = [data_type]
        return self.tree


    def ls_tree(self):
        print("Entries:", self.entries)
        # print("Branches:", self.branches)
        for b in self.branches:
            print(f"Branch: {b:<20}, Type: {self.branches[b]}")

    def disable_branches(self, branches):
        
        pass
            
    def get_entry(self, entry) -> dict:
        self.tree.GetEntry(entry)

        for b in self.branches:
            if "vector" in self.branches[b][0]:
                self.data[b] = self.vector2list(getattr(self.tree, b))
            else:
                self.data[b] = getattr(self.tree, b)
                
        return self.data



def load(filename, printn=2000, start_event=0, end_event=-1, *args, **kwargs):
    # Open the file
    file = ROOT.TFile.Open(filename)
    f1 = tfile(filename)
    tree_name = "digi"
    Tree = file.Get(tree_name)
    branches = [Tree.GetListOfBranches()[i].GetName() for i in range(len(Tree.GetListOfBranches()))]
    entries = Tree.GetEntries()
    print("\nOpening ROOT file", filename)
    print("  Tree:", tree_name)
    # print("  Branches:", branches)
    print("  Entries:", entries)

    if end_event==-1: end_event = entries
    entries_run = [start_event, min(end_event, entries)]
    print("  Entries to read:", entries_run)
    print("  Branches", branches)
    
    # Tell the version of the data
    if 'Digi_detectorID' in branches:
        version = 1  # Simulation as of 2024-05
    elif 'hit_nrg1' in branches:
        version = 2  # UofT test stand data
    else:
        version = 0  # Simulation before 2024-05
        
        
    # Read meatadata
    f1.get_tree("metadata")
    metadata = f1.get_entry(0)
    

    # Make a dictionary
    data_keys=["inactive"]
    if version==0:
        data_keys = ["inactive", "1"]
    elif version==1:
        data_keys = ["inactive", "1", "2"]
    elif version==2:
        data_keys = ["inactive", "1"]                
    data = {key:[] for key in data_keys}

    # Read the file
    for entry in range(*entries_run):
        if (entry+1)%printn==0:  print("  event is ", entry+1)

        for key in data:
            data[key].append([])

        hits = root_hits_extractor(Tree, entry, metadata)
        
        # Group the hits into layer groups
        for hit in hits:
            
            if version==0:
                if  hit.layer<=1 or hit.layer>5:
                    data["inactive"][-1].append(hit)            
                elif hit.layer<=5:
                    data["1"][-1].append(hit)
            
            elif version==1:
                # if  hit.layer<=3: # Two wall and two floor layers
                #     data["inactive"][-1].append(hit)            
                # elif hit.layer<=9: # Top 4 celling tracking layers
                #     data["1"][-1].append(hit)  
                # elif hit.layer<=13: #  4 wall tracking layers
                #     data["2"][-1].append(hit)  
                if  hit.layer in layer_group["inactive"]: # Two wall and two floor layers, plus rejected
                    data["inactive"][-1].append(hit)            
                elif hit.layer in layer_group["1"]: # Top 4 celling tracking layers
                    data["1"][-1].append(hit)  
                elif hit.layer in layer_group["2"]: #  4 wall tracking layers
                    data["2"][-1].append(hit)                 
                    # print(hit.layer) 
                    
            elif version==2:
                if  hit.layer<=4: # Two wall and two floor layers
                    data["1"][-1].append(hit)            

    print("Finished loading file")

    # Make metadata
    metadata={"groups":{}}
    if version==1:
        metadata["groups"]["1"] = {"flip_index":[1,2]}
        metadata["groups"]["2"] = {"flip_index":[1,2]}
    elif version==2:
        metadata["groups"]["1"] = {"flip_index":[1,2]}


    # size = getsize(data)
    # print(f"Memory usage of loaded data {size/1e6:.2f} [MB]")

    return data, metadata


def dump(data, filename):
    
    with open(output_filename,"wb") as f:
        pickle.dump(results, f)    
    # joblib.dump(data,filename+".joblib")

def exists(filename):
    return os.path.exists(filename+".joblib")
    
    


# ----------------------------------------------------------
# Helper functions
# ----------------------------------------------------------

def open(filename):
    file = ROOT.TFile.Open(filename)
    return file

def get_tree(file):
    Tree = file.Get(file.GetListOfKeys()[0].GetName())
    def keys(Tree):
        branch_list = [Tree.GetListOfBranches()[i].GetName() for i in range(len(Tree.GetListOfBranches()))]
        return branch_list
    Tree.keys = types.MethodType(keys,Tree)
    # print(Tree.keys())
    # Tree.keys = branch_list
    
    return Tree    
    
def c2list(x):
    return [x.at(i) for i in range(x.size())]





def make_hits(x, y, z, t, ylayers):
    hits=[]
    for i, layer in enumerate(ylayers):
        if layer%2==1:
            hits.append(datatypes.Hit(x[i], y[i], z[i], t[i], unc_trans, unc_thick, unc_long, unc_time, layer, i, 0, 0))
        else:
            hits.append(datatypes.Hit(x[i], y[i], z[i], t[i], unc_long, unc_thick, unc_trans, unc_time, layer, i, 0, 0))         
            
    return hits  

def make_hits_newsim(x, y, z, t, layer_id, layer_direction, bar_direction, det_id, digi_type):
    hits=[]
    # layer_direction_mod ={0:1, 1:2, 2:0}
    # Digi_type: {0: Geant4 event, 1:cosmic, 2:noise}
    
    for i, layer in enumerate(layer_id):
        unc = [0,0,0]
        other_direction = 0+1+2-layer_direction[i]-bar_direction[i]
        unc[layer_direction[i]] = unc_thick
        unc[bar_direction[i]] = unc_long
        unc[other_direction] = unc_trans
        
        # Sample binomial distribution to decide whether the hit is kept.
        is_noise = (digi_type[i]==2)
        keep_noise = (rng.binomial(1, noise_prob)==1 & is_noise) | (not is_noise)
        keep_hit = (rng.binomial(1, efficiency_prob) and not is_noise) | is_noise
        
        # Unused hits should be assigend a layer number of -1
        if not (keep_noise and keep_hit):
            layer = -1
        hits.append(datatypes.Hit(x[i], y[i], z[i], t[i], unc[0], unc[1], unc[2], unc_time, layer, i,det_id[i], digi_type[i]))
            
    return hits 

def make_hits_teststand(Digi_x, Digi_y, Digi_z, Digi_t, Digi_x_err, Digi_y_err, Digi_z_err, Digi_t_err, Digi_layer, Digi_det_id):
    hits=[]
    # layer_direction_mod ={0:1, 1:2, 2:0}
    
    for i in range(len(Digi_x)):
        hits.append(datatypes.Hit(Digi_x[i], Digi_y[i], Digi_z[i], Digi_t[i], Digi_x_err[i], Digi_y_err[i], Digi_z_err[i], Digi_t_err[i], Digi_layer[i], i, Digi_det_id[i], 1))
            
    return hits 



def root_hits_extractor(Tree, entry, metadata):
    Tree.GetEntry(entry)

    Digi_x = np.array(c2list(Tree.Digi_x))*0.1
    Digi_y = np.array(c2list(Tree.Digi_y))*0.1
    Digi_z = np.array(c2list(Tree.Digi_z))*0.1
    Digi_t = c2list(Tree.Digi_t)
    
    digi_direction = c2list(Tree.Digi_direction)
    Digi_layer_direction = [int(x%10) for x in digi_direction]
    Digi_bar_direction = [int(x/100)%10 for x in digi_direction]    
    
    Digi_det_id = c2list(Tree.Digi_detectorID)
    Digi_layer = [int(x//100000)%1000 for x in Digi_det_id]
    Digi_module = [int(x//100000_000)%1000 for x in Digi_det_id]
    Digi_det = [int(x//100000_000_000)%1000 for x in Digi_det_id]
    for i in range(len(Digi_det)):
        if Digi_det[i]==0:
            Digi_layer[i] +=4
            if Digi_module[i] >=16:
                Digi_layer[i] +=6
        
                

    
    Digi_type = [0]*len(Digi_x)
    return make_hits_newsim(Digi_x, Digi_y, Digi_z, Digi_t, Digi_layer, Digi_layer_direction, Digi_bar_direction, Digi_det_id, Digi_type)

        
        


from types import ModuleType, FunctionType
from gc import get_referents
def getsize(obj):
    """sum size of object & members."""
    # Custom objects know their class.
    # Function objects seem to know way too much, including modules.
    # Exclude modules as well.

    BLACKLIST = type, ModuleType, FunctionType
    if isinstance(obj, BLACKLIST):
        raise TypeError('getsize() does not take argument of type: '+ str(type(obj)))
    seen_ids = set()
    size = 0
    objects = [obj]
    while objects:
        need_referents = []
        for obj in objects:
            if not isinstance(obj, BLACKLIST) and id(obj) not in seen_ids:
                seen_ids.add(id(obj))
                size += sys.getsizeof(obj)
                need_referents.append(obj)
        objects = get_referents(*need_referents)
    return size    
