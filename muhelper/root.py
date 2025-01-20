import ROOT

class tfile:
    """
    Helper calss for ROOT TFile.
    It wraps the mostly used funcitons to open and read objects from ROOT file.

    
    tfile = ROOT.TFile.Open(fname_digi)
    tree_name = tfile.GetListOfKeys()[1].GetName()
    Tree = tfile.Get(tree_name)
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
