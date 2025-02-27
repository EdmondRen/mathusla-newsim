import numpy as np
import scipy as sp
import scipy
from scipy.spatial import ConvexHull

import sys
sys.path.append("../../python")
# from simhelper import root, generator, parser
from simhelper import datatypes, util
# from simhelper import helper_basic as hp





class RRQ:
    CONST_CMS_LOC = np.array([(-70-39/2)*1000, 0, -85.47*1000])
    CONST_CMS_DIRECTION = -CONST_CMS_LOC/np.linalg.norm(CONST_CMS_LOC)
    CONST_VETO_DET = np.array([0])
    CONST_VETO_LAYER = np.array([0,1,2,3])
    CONST_VETO_WALLLAYER = np.array([0,1])
    CONST_VETO_FLOORLAYER = np.array([2,3])
    CONST_NOISE_RATE_SUPRESS = 0.1

    def __init__(self,  event, seed=1):
        ## Random number generator
        RNG = np.random.default_rng(seed=seed)

        ## Get the vertex with most tracks
        ntracks=[]
        for j in range(len(event.vertices)):
            v = event.vertices[j]
            ntracks.append(v.ntracks)
        i_vertex = np.argmax(ntracks)
        self.event=event
        self.v = v =  event.vertices[i_vertex]

        ## Find the CMS to vertex direction
        direction_longi = v.params[:3] - RRQ.CONST_CMS_LOC
        direction_longi = direction_longi/np.linalg.norm(direction_longi)
        direction_horizontal = np.cross(direction_longi, [0,0,1])
        direction_horizontal = direction_horizontal/np.linalg.norm(direction_horizontal)
        direction_other = np.cross(direction_horizontal,direction_longi)    
        self.direction_longi = direction_longi
        self.direction_horizontal = direction_horizontal
        self.direction_other = direction_other

        ## Check if the vertex is mostly with top or side tracks
        tracks_are_top = [self.event.tracks[i].iv_index == 2 for i in self.v.track_ids]
        self.vertex_topfrac  = (sum(tracks_are_top)/len(tracks_are_top))
            
    
    @staticmethod
    def angle(a,b):
        """Find angle between two vectors"""
        opening_angle = np.arccos(np.dot(a, b)/( np.linalg.norm(a) * np.linalg.norm(b) ))
        return opening_angle

    @staticmethod
    def angle_signed(a, b):
        """
        Calculate signed angle between two 2-D vectors
        """
        mag1 = np.linalg.norm(a)
        mag2 = np.linalg.norm(b)
        
        # Compute the dot product and cross product
        dot = np.dot(a, b)/( mag1*mag2 )
        cross = np.cross(a, b)/( mag1*mag2 )
        angle = np.arctan2(cross, dot) # Use atan2 to compute the signed angle
        # angle = angle%(2*np.pi) # Adjust to [0, 2π)
        return angle
    

    @staticmethod
    def min_enclosing_cone(vectors):
        """
        Finds the minimum enclosing cone containing all given 3D vectors using a robust convex hull approach.
        Handles cases where vectors are coplanar, near-coplanar, or when only two vectors are provided.
    
        Parameters:
        vectors (np.ndarray): An (N, 3) array where each row is a 3D vector.
    
        Returns:
        tuple: (axis, angle) where:
            - axis (np.ndarray): The optimal unit vector representing the cone axis.
            - angle (float): The minimum half-angle (in radians) of the enclosing cone.
        """
        num_vectors = len(vectors)
    
        # Normalize the vectors
        unit_vectors = vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    
        if num_vectors == 1:
            # Only one vector: the cone is degenerate (zero angle)
            return unit_vectors[0], 0.0, 0
    
        elif num_vectors == 2:
            # Two vectors: The cone axis is their normalized sum, and angle is half of their separation angle
            axis = (unit_vectors[0] + unit_vectors[1]) / np.linalg.norm(unit_vectors[0] + unit_vectors[1])
            angle = np.arccos(np.dot(unit_vectors[0], unit_vectors[1])) / 2
            return axis, angle, angle
    
        # Compute SVD to check if vectors are coplanar
        U, S, Vt = np.linalg.svd(unit_vectors - np.mean(unit_vectors, axis=0))
    
        if S[2] < 1e-10:  # If the smallest singular value is near zero, vectors are coplanar
            # Coplanar case: Find a normal vector to the plane
            normal = Vt[2]  # The last row of Vt gives the normal direction
    
            # Choose two orthogonal unit vectors spanning the plane
            u = Vt[0]  # First principal direction
            v = Vt[1]  # Second principal direction
    
            # Project unit vectors onto the 2D plane
            projected_2D = np.column_stack((np.dot(unit_vectors, u), np.dot(unit_vectors, v)))
    
            # Compute the convex hull in 2D
            hull = ConvexHull(projected_2D)
    
            # Find the best 2D enclosing cone axis
            centroid_2D = np.mean(projected_2D[hull.vertices], axis=0)
            centroid_2D /= np.linalg.norm(centroid_2D)
    
            # Compute the max angle deviation in 2D
            min_cos_theta = np.min(np.dot(projected_2D[hull.vertices], centroid_2D))
            min_cone_angle = np.arccos(min_cos_theta)
            mean_cone_angle = np.mean(np.arccos(np.dot(projected_2D[hull.vertices], centroid_2D)))
    
            # Convert back to 3D axis
            axis = centroid_2D[0] * u + centroid_2D[1] * v
    
        else:
            # Compute the convex hull of the unit vectors in 3D
            hull = ConvexHull(unit_vectors)
    
            # Find the best cone axis: Compute the centroid of the convex hull
            centroid = np.mean(unit_vectors[hull.vertices], axis=0)
            axis = centroid / np.linalg.norm(centroid)
    
            # Find the maximum angle from axis to any unit vector in the hull
            min_cos_theta = np.min(np.dot(unit_vectors[hull.vertices], axis))
            min_cone_angle = np.arccos(min_cos_theta)
            mean_cone_angle = np.mean(np.arccos(np.dot(unit_vectors[hull.vertices], axis)))
    
        return axis, min_cone_angle, mean_cone_angle
    
    def project_vector(self, vec):
        return [self.direction_longi.dot(vec), self.direction_horizontal.dot(vec), self.direction_other.dot(vec)]


    def get_ndigi_before(self, limit = 100):
        is_before = [(self.v.t0-limit) <digi.t < self.v.t0 for digi in self.event.digis]
        return np.sum(is_before)   

    def get_ndigi_veto_and_active(self, limit = 100):
        is_before = np.array([(self.v.t0-limit) <digi.t < self.v.t0 for digi in self.event.digis])
        is_veto = np.array([(digi.copy_det in RRQ.CONST_VETO_DET) for digi in self.event.digis])
        n_veto = np.sum(is_veto)
        n_active = len(is_veto) - n_veto
        return [n_veto, n_active, sum(is_before&is_veto), sum(is_before&~is_veto)]

    def get_slowest_track(self):
        vs = [self.event.tracks[i].vabs for i in self.v.track_ids]
        vs.sort()
        return np.mean(vs[:2])

    def get_v_downward_track(self):
        vz = [self.event.tracks[i].vdirection[2]<0 for i in self.v.track_ids]
        return np.sum(vz)

    def get_ev_downward_track(self):
        vz = [self.event.tracks[i].vdirection[2]<0 for i in range(len(self.event.tracks))]
        return np.sum(vz)    

    def eval_cone(self):
        """
        Find the properties of the cone formed by the tracks
        * open angle and direction
        * angle difference relative to CMS. This is evaluated on two planes. 
            < 0 means CMS direction is inside the cone.

        Return:  open_angle, axis, angle_diffh, angle_diffv
        open_angle: float
        axis: 3-vec, [x,y,z]
        angle_diffh: float
        angle_diffv: float
        """

        angles_h = []
        angles_v = []
        track_directions = []
        
        for i in self.v.track_ids:
            track_directions.append(self.event.tracks[i].vdirection)
            vdirection_along_cms = self.project_vector(track_directions[-1])
            angles_h.append(self.angle_signed([1,0], [vdirection_along_cms[0], vdirection_along_cms[1]]))
            angles_v.append(self.angle_signed([1,0], [vdirection_along_cms[0], vdirection_along_cms[2]]))

        angles_h_sign = np.sign(angles_h)
        angles_v_sign = np.sign(angles_v)
        angle_diffh = np.min(np.abs(angles_h))
        angle_diffv = np.min(np.abs(angles_v))
        if any(angles_h_sign!=angles_h_sign[0]):
            angle_diffh *= -1
        if any(angles_v_sign!=angles_v_sign[0]):
            angle_diffv *= -1            
            
        # Test the function with a set of 3D vectors
        axis, open_angle_max, open_angle_mean = self.min_enclosing_cone(track_directions)

        angle_diff_abs = self.angle(self.direction_longi, axis)

        return open_angle_max, open_angle_mean, axis, angle_diffh, angle_diffv, angle_diff_abs


    def eval_hits_time(self, ndiv = 4,  z_veto = 8000, z_veto_floor = 3000):
        digi_xyzt = np.array([digi.xyzt for digi in self.event.digis])
        digi_x, digi_y, digi_z, digi_t = digi_xyzt.T
        digi_t_inds = np.argsort(digi_t).astype(int)
        digi_t_sorted = digi_t[digi_t_inds]
        
        mask, tmean, t_std = util.removeoutliers(digi_t, skew_itermax=1, std_itermax=1)
        t_std = max(t_std*2, 300)
        t_std = min(t_std, 300)
        digi_t_inds = digi_t_inds[(digi_t_sorted<(tmean + t_std)) & 
                                  (digi_t_sorted>(tmean - t_std))]
        
        rtn_1 = True
        rtn_2 = 10000
        rtn_3 = True
        if (len(digi_t_inds)>12):
            xs,ys,zs,ts = [],[],[],[]
            for igroup in range(ndiv):
                inds = digi_t_inds[len(digi_t_inds)//ndiv*igroup: len(digi_t_inds)//ndiv*(igroup+1)]
                xs.append(np.mean(np.array(digi_x)[inds]))
                ys.append(np.mean(np.array(digi_y)[inds]))
                zs.append(np.mean(np.array(digi_z)[inds]))
                ts.append(np.mean(np.array(digi_t)[inds]))

            # 1. Check if [vertex is mostly on top tracker] 
            #  and [any duration have verticle position too low]
            if (self.vertex_topfrac>0.5 and any(np.array(zs[:]) < z_veto)) or \
                (sum(np.array(zs) < z_veto_floor)>=2):
                rtn_1 = False

            # 2. Check if the overall trend is going downwards
            dzdt, intercept, r, p, std_err = scipy.stats.linregress(ts[:3], zs[:3])
            dxdz, intercept, r, p, std_err = scipy.stats.linregress(zs[:3], xs[:3])
            dydz, intercept, r, p, std_err = scipy.stats.linregress(zs[:3], ys[:3])
            rtn_2 = dzdt
            
            # 3. Check how well the overall direction align with CMS direction.
            # direction = [dxdz, dydz, 1]
            # angle_to_cms = self.angle(direction, self.direction_longi)
            rtn_3 = np.mean(zs)

        return rtn_1, rtn_2, rtn_3
    

class rq:
    def __init__(self, dictionary):
        self.data = dictionary
        self.cuts = {}  # {name: mask}
        self.cuts_order = {} # {name: index}
        self.cuts_func = {} # {name: func}

    def __getitem__(self, key):
        return self.data[key]

    def add_cut(self, func, name):
        self.cuts[name] = func(self.data)

    def get_cut(self, name):
        return self.cuts[name]


def run_processing(file, entries = -1):

    ## Get metadata
    file.get_tree("metadata_digi")
    metadata_digi = file.get_entry(0)
    
    # Get data
    file.get_tree("data;1")
    print("Entries", file.entries)

    keys = ["Run_number", "Evt_number", "ROOT_entry",
            "event_ntracks","event_nhits","event_nvertices",
            "event_ntrack_reconstructable", "event_ntrack_reconstructable_primary", "event_ntrack_recon", "event_vntrk_max",
            
            # Vertex stats
            "vertex_topfrac",
            "vertex_ntracks", "vertex_ndigi", "vertex_chi2", "vertex_prob", "vertex_residual", "vertex_error", "vertex_residual_longitrans",
            "vertex_ntracklet_0", "vertex_ntracklet_2", "vertex_ntracklet_3+",   
            "vertex_ndownward_track", "event_ndownward_track",
    
            # 
            "event_ndigi_veto", "event_ndigi_active", 
            "vertex_ndigi_veto_before", "vertex_ndigi_active_before", "vertex_slowest_track",
            "vertex_open_angle", "vertex_cms_angle_h", "vertex_cms_angle_v", "vertex_cms_angle",
            "vertex_hits_trend_1", "vertex_hits_trend_2", "vertex_hits_trend_3", 
            
           ] 
    
    res = {key:[] for key in keys}    
    
    if entries==-1:
        entries = file.entries    
    entries = min(entries, file.entries )
    
    for i in range(entries):
        if i%100==0:
            print(i, end="\r")
        
        data = file.get_entry(i)
        event = datatypes.Event(data, metadata_digi, parse_truth=False) 
        rrq = RRQ(event)
    
        # Basic information
        res["ROOT_entry"].append(i)
        res["Run_number"].append(event.Run_number)
        res["Evt_number"].append(event.Evt_number)
    
        # Use the vertex with most tracks
        v =  rrq.v    
        res["vertex_topfrac"].append(rrq.vertex_topfrac)
        res["vertex_ntracks"].append(v.ntracks)
        res["vertex_ndigi"].append(v.nhits)
        res["vertex_chi2"].append(v.chi2)
        res["vertex_prob"].append(scipy.stats.chi2.cdf(res["vertex_chi2"][-1], res["vertex_ntracks"][-1]*3-4))
        res["vertex_residual"].append(v.params - event.genvertices.vertices[0])
        res["vertex_error"].append(np.sqrt(np.diagonal(v.cov)))
        res["vertex_ntracklet_0"].append(v.vertex_ntracklet_0)
        res["vertex_ntracklet_2"].append(v.vertex_ntracklet_2)
        res["vertex_ntracklet_3+"].append(v.vertex_ntracklet_3)
        res["vertex_residual_longitrans"].append(rrq.project_vector(res["vertex_residual"][-1][:3]))
    
        ## Event level information
        res["event_nhits"].append(len(event.tracks))
        res["event_ntracks"].append(len(event.tracks))
        res["event_nvertices"].append(len(event.vertices))
        res["event_ntrack_reconstructable"].append(sum(event.digi_truth_track_nhits>=4))
        res["event_ntrack_reconstructable_primary"].append(sum((event.digi_truth_track_nhits>=4) & (event.digi_truth_track_ids<=len(event.genparticles))))
        res["event_ntrack_recon"].append(len(event.tracks))
    
        ## Calculate some features for cuts
        res["vertex_slowest_track"].append(rrq.get_slowest_track())
        res["vertex_ndownward_track"].append(rrq.get_v_downward_track())
        res["event_ndownward_track"].append(rrq.get_ev_downward_track())
        # Number of hits
        ndigi_veto,ndigi_active, ndigi_veto_before,ndigi_active_before = rrq.get_ndigi_veto_and_active(limit=100)
        res["event_ndigi_veto"].append(ndigi_veto)
        res["event_ndigi_active"].append(ndigi_active)
        res["vertex_ndigi_veto_before"].append(ndigi_veto_before)
        res["vertex_ndigi_active_before"].append(ndigi_active_before)
        
        # Opening angle
        open_angle_max, open_angle_mean, axis, angle_diffh, angle_diffv, angle_diff_abs = rrq.eval_cone()
        res["vertex_open_angle"].append(open_angle_mean)
        res["vertex_cms_angle_h"].append(angle_diffh)
        res["vertex_cms_angle_v"].append(angle_diffv)
        res["vertex_cms_angle"].append(angle_diff_abs)

        r1,r2,r3 = rrq.eval_hits_time()
        res["vertex_hits_trend_1"].append(r1)
        res["vertex_hits_trend_2"].append(r2)
        res["vertex_hits_trend_3"].append(r3)
        
    
    for key in res:
        res[key] = np.array(res[key])

    res["event_ndigi_active_after"] = res["event_ndigi_active"] - res["vertex_ndigi_active_before"]
    res["vertex_ndigi_before"] = res["vertex_ndigi_active_before"] + res["vertex_ndigi_veto_before"]
    print("Finished")

    return rq(res)