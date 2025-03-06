import numpy as np
import scipy as sp
import scipy
from scipy.spatial import ConvexHull
from scipy.optimize import minimize

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
    box_halflen_mm = 10700 * 2; # [mm]
    box_fiducial_shrink = 300;  # [mm] Shrink inside by 30cm
    box_x = [-box_halflen_mm + box_fiducial_shrink, box_halflen_mm - box_fiducial_shrink]
    box_y = [-box_halflen_mm + box_fiducial_shrink, box_halflen_mm - box_fiducial_shrink]
    box_z = [1020 + box_fiducial_shrink, 12000 - box_fiducial_shrink]
    
    def __init__(self,  event, root_data, seed=1):
        ## Random number generator
        RNG = np.random.default_rng(seed=seed)

        ## Get the vertex with most tracks
        i_vertex = -1
        has_valid_vertex = False
        v = None
        if ("Vertex_selected_ind" in root_data):
            # In this case, the skim script already pre-selected the vertex. 
            i_vertex = root_data["Vertex_selected_ind"]
            v =  event.vertices[i_vertex]
            has_valid_vertex = not v.dropped
            
        if not has_valid_vertex:
            ntracks=[v.ntracks for v in event.vertices]

            passed = False;
            current_vertex_ind = -1;
            for iv in np.argsort(ntracks)[::-1]:
                v = event.vertices[iv]
                if v.dropped:
                    continue
                #  0. Vertex is inside the box
                is_inside = (self.box_x[0] < v.params[0] < self.box_x[1] and \
                             self.box_y[0] < v.params[1] < self.box_y[1] and \
                             self.box_z[0] < v.params[2] < self.box_z[1])

                # 1. All tracks are outwards
                # 2. Exists one upward tracks
                is_outward = True;
                is_upward = False;
                for j in v.track_ids:
                    t = event.tracks[j];
                    track_is_outward = (t.iv_index == 2 and t.params_full[7] > 0) or \
                                            (t.iv_index == 0 and t.params_full[7] > 0);
                    track_is_upward = (t.iv_index == 2 and t.params_full[7] > 0) or \
                                            (t.iv_index == 0 and t.params_full[6] > .12);
                    if (not track_is_outward):
                        is_outward = False;
                    if (track_is_upward):
                        is_upward = True;
                passed = (is_inside and is_outward and is_upward);
                current_vertex_ind = iv;
    
                if (passed):
                    has_valid_vertex = True
                    break
                    
        self.has_valid_vertex = has_valid_vertex

        if self.has_valid_vertex:
            self.event=event
            self.v = v
    
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
        # Objective function: maximize the max cosine similarity (minimize max angular deviation)
        def objective(axis):
            axis = axis / np.linalg.norm(axis)  # Ensure unit length
            max_cos_theta = np.max(np.dot(unit_vectors, axis))
            return -max_cos_theta  # We minimize the negative to find the best axis
            
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
    
        # Initial guess for the axis: normalized mean of vectors
        init_axis = np.mean(unit_vectors, axis=0)
        init_axis /= np.linalg.norm(init_axis)
    
        # Objective function: maximize the maximum cosine similarity (minimize max angular deviation)
        def objective(axis):
            axis = axis / np.linalg.norm(axis)  # Ensure unit length
            min_cos_theta = np.min(np.dot(unit_vectors, axis))
            return -min_cos_theta  # We minimize the negative to find the best axis
    
        # Constrain the optimization so that the axis remains a unit vector
        result = minimize(objective, init_axis, method='Powell')
    
        # Extract the optimized axis and compute the enclosing angle
        axis = result.x / np.linalg.norm(result.x)    
        min_cos_theta = np.min(np.dot(unit_vectors, axis))
        min_cone_angle = np.arccos(min_cos_theta)
        mean_cone_angle = np.mean(np.arccos(np.dot(unit_vectors, axis)))            
    
        return init_axis, min_cone_angle, mean_cone_angle
    
    def project_vector(self, vec):
        return [self.direction_longi.dot(vec), self.direction_horizontal.dot(vec), self.direction_other.dot(vec)]


    def get_ndigi_before(self, limit = 100):
        is_before = [(self.v.t0-limit) <digi.t < self.v.t0 for digi in self.event.digis if not digi.dropped]
        return np.sum(is_before)   

    def get_ndigi_veto_and_active_timelimit(self, limit = 100, ):
        """ 
        Return:
        n_veto, 
        n_active, 
        n_veto_pre_veto_timelimit, 
        n_veto_pre_active_timelimit
        """
        
        is_valid = np.array([(not digi.dropped) for digi in self.event.digis])
        is_before = np.array([(self.v.t0-limit) <digi.t < self.v.t0 for digi in self.event.digis])
        is_veto = np.array([(digi.copy_det in RRQ.CONST_VETO_DET) for digi in self.event.digis])
        n_veto = np.sum(is_veto & is_valid)
        n_active = sum(is_valid) - n_veto
        n_veto_pre_veto_timelimit = sum(is_before&is_veto&is_valid)
        n_veto_pre_active_timelimit = sum(is_before&~is_veto&is_valid)
        return [n_veto, n_active, n_veto_pre_veto_timelimit, n_veto_pre_active_timelimit]

    def get_ndigi_veto_and_comp(self, limit = 200, speed_comp_limit = 30):
        """
        Find number of veto hits that are
        * Later than the vertex
        * Compatible with the vertex at speed of light
        """
        is_valid = np.array([(not digi.dropped) for digi in self.event.digis])
        is_after = np.array([self.v.t0 < digi.t < (self.v.t0+limit) for digi in self.event.digis])
        is_above = np.array([self.v.z0 < digi.z for digi in self.event.digis])
        is_veto = np.array([(digi.copy_det in RRQ.CONST_VETO_DET) for digi in self.event.digis])
        is_comp = []
        # for i in np.flatnonzero(is_after&is_veto):
        metrics_dt = []
        metrics_speed = []
        for i in range(len(is_after)):
            digi = self.event.digis[i]
            dr = np.linalg.norm([digi.x - self.v.x0, digi.y - self.v.y0, digi.z - self.v.z0])
            dt = abs(digi.t - self.v.t0)
            metric_dt = abs(dt - dr/295)
            metric_speed = abs(abs(dr/dt) - 295)
            # is_comp.append(np.abs(metric)<5)
            is_comp.append(metric_speed<speed_comp_limit)      
            metrics_dt.append(metric_dt)
            metrics_speed.append(abs(dr/dt)) 

        metrics_dt=np.array(metrics_dt)
        metrics_speed=np.array(metrics_speed)
        metrics_dt=metrics_dt[is_valid&~is_veto&is_after&is_above]
        metrics_speed=metrics_speed[is_valid&~is_veto&is_after&is_above]

        metrics_dt = []
        metrics_speed = []
        

        n_veto_after = sum(is_after&is_veto&is_valid)
        n_veto_after_comp = sum(is_after&is_veto&is_comp&is_valid)
        n_active_after = sum(is_after&~is_veto&is_valid)
        n_active_after_comp = sum(is_after&~is_veto&is_comp&is_above&is_valid)        
        return n_veto_after, n_veto_after_comp, n_active_after, n_active_after_comp, metrics_dt, metrics_speed

    def get_slowest_track(self):
        vs = [self.event.tracks[i].vabs for i in self.v.track_ids]
        vs.sort()
        return np.mean(vs[:2])

    def get_v_downward_track(self):
        vz = [self.event.tracks[i].vdirection[2]<0 for i in self.v.track_ids]
        return np.sum(vz)

    def get_ev_downward_track(self, tlimit=200):
        is_downward = np.array([self.event.tracks[i].vdirection[2]<0 for i in range(len(self.event.tracks)) if not self.event.tracks[i].dropped])
        is_inward = np.array([(self.event.tracks[i].iv_index == 2 and self.event.tracks[i].params_full[7] < 0) for i in range(len(self.event.tracks)) if not self.event.tracks[i].dropped])
        is_sidewall = np.array([(self.event.tracks[i].iv_index == 2) for i in range(len(self.event.tracks)) if not self.event.tracks[i].dropped])
        mask_tlimit = np.array([abs(self.event.tracks[i].params_full[3]-self.v.t0)<tlimit for i in range(len(self.event.tracks)) if not self.event.tracks[i].dropped])

        n_downward_top = sum(is_downward & (~is_sidewall) & mask_tlimit)
        n_downward_side = sum((is_downward & is_sidewall & mask_tlimit) | (is_inward & mask_tlimit))
        track_ids_downward = np.flatnonzero(is_downward | is_inward)
        
        min_downward_dist = []
        for i in track_ids_downward:
            tr = self.event.tracks[i]
            tr_position = tr.at_t(self.v.t0)
            dr = np.linalg.norm(tr_position - self.v.params[:3])
            min_downward_dist.append(dr)
        min_downward_dist = min(min_downward_dist) if len(min_downward_dist)>0 else 40000
        
        return n_downward_top, n_downward_side, min_downward_dist

    def get_vertex_track_dist(self):
        downward_dist = []
        for i in self.v.track_ids:
            tr = self.event.tracks[i]
            tr_position = tr.at_t(self.v.t0)
            dr = np.linalg.norm(tr_position - self.v.params[:3])
            downward_dist.append(dr)     

        return sum(downward_dist)

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
        angle_diffh_mean = np.mean(angles_h)
        angle_diffv_mean = np.mean(angles_v)  
        angle_diffh_span = np.max(angles_h) - np.min(angles_h)
        angle_diffv_span = np.max(angles_v) - np.min(angles_v)
        if any(angles_h_sign!=angles_h_sign[0]):
            angle_diffh *= -1
            # angle_diffh_mean *= -1
        if any(angles_v_sign!=angles_v_sign[0]):
            angle_diffv *= -1        
            # angle_diffv_mean *= -1
            
        # Test the function with a set of 3D vectors
        axis, open_angle_max, open_angle_mean = self.min_enclosing_cone(track_directions)

        angle_diff_abs = self.angle(self.direction_longi, axis)

        return open_angle_max, open_angle_mean, axis, angle_diffh, angle_diffv, angle_diff_abs, angle_diffh_mean, angle_diffv_mean, angle_diffh_span, angle_diffv_span



    def eval_hits_time_exclude(self, ndiv = 2,  z_veto = 8000, z_veto_floor = 3000):
        digi_xyzt = np.array([digi.xyzt for digi in self.event.digis])
        digi_x, digi_y, digi_z, digi_t = digi_xyzt.T
        digi_t_inds = np.argsort(digi_t).astype(int)
        digi_t_sorted = digi_t[digi_t_inds]

       
        mask, tmean, t_std = util.removeoutliers(digi_t, skew_itermax=1, std_itermax=1)
        t_std = max(t_std*2, 300)
        t_std = min(t_std, 300)
        digi_t_inds = digi_t_inds[(digi_t_sorted<(tmean + t_std)) & 
                                  (digi_t_sorted>(tmean - t_std))]

        # Exclude hits that are used in track
        inds_used = []
        for t in self.event.tracks:
            inds_used.extend(list(t.hit_ids)) 
        for i in range(len(digi_t_inds))[::-1]:
            if digi_t_inds[i] in inds_used:
                np.delete( digi_t_inds, [i])
                continue
            if self.event.digis[i].copy_det in RRQ.CONST_VETO_DET:
                np.delete( digi_t_inds, [i])
                continue            
        
        rtn_1 = True
        rtn_2 = 10000
        if (len(digi_t_inds)>6):
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
                (sum(np.array(zs) < z_veto_floor)>=1):
                rtn_1 = False

            # 2. Check if the overall trend is going downwards
            dzdt = (zs[1]-zs[0])/(ts[1]-ts[0])
            rtn_2 = dzdt

        return rtn_1, rtn_2


class CutItem:
    def __init__(self, name, func, cut_order = -1):
        # Descriptive name of the cut
        self.name = name

        # Cut function (to be applied to the dictionary)
        self.func = func
               
        # cut_order: int, the priority of the cut. Lower number will be applied first
        self.cut_order=cut_order
        
        # Event mask
        self.mask = None
        
        # Passage fraction
        self.passage_fraction = None
        
    def get_eff(self):
        return self.passage_fraction

    def get_mask(self):
        return self.mask

    def get_func(self):
        return self.func    
        
    def apply(self, data):
        self.mask = self.func(data)
        self.passage_fraction = sum(self.mask)/len(self.mask)
        return self.mask

class RQ_dict:
    def __init__(self, data):
        self.data = data
        self.cuts_dict = {}  # dict of {name: CutItem}
        self.cuts_name = []
        self.cuts_active = []

        name = "More than 2 tracks"
        func = lambda res: res["vertex_ntracks"]>2
        self.cuts_dict[name] = CutItem(name, func, cut_order = -1)
        self.cuts_dict[name].apply(self.data)
        self.cuts_name.append(name)
        self.mask_2 = self.cuts_dict["More than 2 tracks"].get_mask()
        
        self.add_cut(lambda res: np.ones_like(res["vertex_ntracks"], dtype=bool), "True")   
        self.mask_true = self.cuts_dict["True"].get_mask()

    def __getitem__(self, key):
        return self.data[key]

    def list_cut(self):
        for i in range(len(self.cuts_name)):
            name = self.cuts_name[i]
            mask = self.cuts_dict[name].get_mask()
            print(f"Cut {i:<2}: {name: <20}, npassed {sum(mask):>6}, passage fraction {sum(mask)/len(mask):.7f}")
        pass

    def print_active(self, show=True):
        mask = self.mask_true
        info = dict() # Each entry is a [num passed, frac passed]  pair
        info["Total"] = [len(mask), 1]
        if show:
            print(f"Total events: {len(mask)}")
            
        for i in self.cuts_active:
            name = self.cuts_name[i]
            mask = mask & self.cuts_dict[name].get_mask()
            info[name] = [sum(mask), sum(mask)/len(mask)]
            if show:
                print(f"Cut {i:<2}: {name: <20}, npassed {sum(mask):>6}, passage fraction {sum(mask)/len(mask):.7f}")
        return info        

    def add_cut(self, func, name, cut_order = -1, PRINT=False):
        self.cuts_dict[name] = CutItem(name, func, cut_order = cut_order)
        self.cuts_dict[name].apply(self.data)
        self.cuts_name = list(self.cuts_dict.keys())
        mask = self.cuts_dict[name].get_mask()

        if PRINT:
            print(f"Add cut: {name}, passage fraction {sum(mask)/len(mask)}, (& > 2 tracks): {sum(mask&self.mask_2)/len(self.mask_2)}")

    def get_cut(self, name):
        if type(name) is int:
            name = self.cuts_name[i]
        return self.cuts_dict[name].get_mask()


def run_processing(file, entries = -1, efficiency = 1, min_nhits=4):

    ## Get metadata
    file.get_tree("metadata_digi")
    metadata_digi = file.get_entry(0)
    
    # Get data
    file.get_tree("data;1")
    print("Entries", file.entries)

    keys = ["Run_number", "Evt_number", "ROOT_entry",
            "gen_p3","gen_xyzt", "gen_pdgID",
            "event_ntracks","event_nhits","event_nvertices",
            "event_ntrack_reconstructable", "event_ntrack_reconstructable_primary", "event_vntrk_max",
            "event_track_nhits", "event_track_nhits_upward",
            "event_ndigi_veto", "event_ndigi_active", 
            
            # Vertex stats
            "vertex_topfrac", "vertex_xyzt",
            "vertex_ntracks", "vertex_ndigi", "vertex_chi2", "vertex_prob", "vertex_residual", "vertex_error", "vertex_residual_longitrans",
            "vertex_ntracklet_0", "vertex_ntracklet_2", "vertex_ntracklet_3+",   
            "vertex_ndownward_track", "event_ndownward_track",
    
            # 
            "vertex_ndigi_veto_before_limited", "vertex_ndigi_active_before_limited", "vertex_slowest_track",
            "vertex_ndigi_veto_after", "vertex_ndigi_veto_after_comp",
            "vertex_ndigi_active_after", "vertex_ndigi_active_after_comp",
            "vertex_track_dist",
            
            "vertex_open_angle", 
            "vertex_cms_angle_h", "vertex_cms_angle_v", 
            "vertex_cms_angle_h_mean", "vertex_cms_angle_v_mean", 
            "vertex_cms_angle_h_span", "vertex_cms_angle_v_span", 
            "vertex_cms_angle",
            "vertex_hits_trend_1b", "vertex_hits_trend_2b",
            
            "vertex_comp_metric_dt", "vertex_comp_metric_speed",
           ] 
    
    res = {key:[] for key in keys}    
    
    if entries==-1:
        entries = file.entries    
    entries = min(entries, file.entries )
    
    for i in range(entries):
        if i%100==0:
            print(i, end="\r")
        
        root_data = file.get_entry(i)
        event = datatypes.Event(root_data, metadata_digi, parse_truth=False, detector_efficiency = efficiency, min_nhits=min_nhits)
        rrq = RRQ(event, root_data)

        # Need to have vertex that passed preliminary cuts in RRQ()
        if not rrq.has_valid_vertex:
            continue
    
        # Basic information
        res["ROOT_entry"].append(i)
        res["Run_number"].append(event.Run_number)
        res["Evt_number"].append(event.Evt_number)

        if len(event.genparticles)>0:
            res["gen_p3"].append(event.genparticles[0].momentum)
            res["gen_xyzt"].append(event.genparticles[0].xyzt)
            res["gen_pdgID"].append(event.genparticles[0].pdg)    
    
        # Use the vertex with most tracks
        v =  rrq.v    
        res["vertex_xyzt"].append(v.params)    
        res["vertex_topfrac"].append(rrq.vertex_topfrac)
        res["vertex_ntracks"].append(v.ntracks)
        res["vertex_ndigi"].append(v.nhits)
        res["vertex_chi2"].append(v.chi2)
        res["vertex_prob"].append(scipy.stats.chi2.cdf(res["vertex_chi2"][-1], res["vertex_ntracks"][-1]*3-4))
        if len(event.genvertices.vertices)>0:
            res["vertex_residual"].append(v.params - event.genvertices.vertices[0])
            res["vertex_residual_longitrans"].append(rrq.project_vector(res["vertex_residual"][-1][:3]))
        res["vertex_error"].append(np.sqrt(np.diagonal(v.cov)))
        res["vertex_ntracklet_0"].append(v.vertex_ntracklet_0)
        res["vertex_ntracklet_2"].append(v.vertex_ntracklet_2)
        res["vertex_ntracklet_3+"].append(v.vertex_ntracklet_3)
    
        ## Event level information
        res["event_nhits"].append(event.get_ndigis())
        res["event_ntracks"].append(event.get_ntracks())
        res["event_nvertices"].append(event.get_nvertices())
        res["event_ntrack_reconstructable"].append(sum(event.digi_truth_track_nhits>=4))
        res["event_ntrack_reconstructable_primary"].append(sum((event.digi_truth_track_nhits>=4) & (event.digi_truth_track_ids<=len(event.genparticles))))
        res["event_track_nhits"].append(sum([t.nhits for t in event.tracks if not t.dropped]))
        res["event_track_nhits_upward"].append(sum([t.nhits for t in event.tracks if t.vdirection[2] > 0]))
    
        ## Calculate some features for cuts
        res["vertex_slowest_track"].append(rrq.get_slowest_track())
        res["vertex_ndownward_track"].append(rrq.get_v_downward_track())
        res["event_ndownward_track"].append(rrq.get_ev_downward_track())
        # Number of hits
        ndigi_veto,ndigi_active, ndigi_veto_before,ndigi_active_before = rrq.get_ndigi_veto_and_active_timelimit(limit=100)
        res["event_ndigi_veto"].append(ndigi_veto)
        res["event_ndigi_active"].append(ndigi_active)
        res["vertex_ndigi_veto_before_limited"].append(ndigi_veto_before)
        res["vertex_ndigi_active_before_limited"].append(ndigi_active_before)
        # Number of veto that are later than the vertex and compatible with speed of light
        n_veto_after, n_veto_after_comp, n_active_after, n_active_after_comp, metrics_dt, metrics_speed = rrq.get_ndigi_veto_and_comp(limit=200)    
        res["vertex_ndigi_veto_after"].append(n_veto_after)
        res["vertex_ndigi_veto_after_comp"].append(n_veto_after_comp)
        res["vertex_ndigi_active_after"].append(n_active_after)
        res["vertex_ndigi_active_after_comp"].append(n_active_after_comp)  
        res["vertex_track_dist"].append(rrq.get_vertex_track_dist())

        # Temp:
        res["vertex_comp_metric_dt"].append(metrics_dt)        
        res["vertex_comp_metric_speed"].append(metrics_speed)        
        
        # Opening angle
        open_angle_max, open_angle_mean, axis, angle_diffh, angle_diffv, angle_diff_abs, angle_diffh_mean, angle_diffv_mean, angle_diffh_span, angle_diffv_span = rrq.eval_cone()
        res["vertex_open_angle"].append(open_angle_mean)
        res["vertex_cms_angle_h"].append(angle_diffh)
        res["vertex_cms_angle_v"].append(angle_diffv)
        res["vertex_cms_angle_h_mean"].append(angle_diffh_mean)
        res["vertex_cms_angle_v_mean"].append(angle_diffv_mean)  
        res["vertex_cms_angle_h_span"].append(angle_diffh_span)
        res["vertex_cms_angle_v_span"].append(angle_diffv_span)          
        res["vertex_cms_angle"].append(angle_diff_abs)

        r1,r2 = rrq.eval_hits_time_exclude()
        res["vertex_hits_trend_1b"].append(r1)
        res["vertex_hits_trend_2b"].append(r2)        
        
    
    for key in res:
        try:
            res[key] = np.array(res[key])
        except:
            res[key] = np.array(res[key], dtype="object")
            

    # res["event_ndigi_active_after"] = res["event_ndigi_active"] - res["vertex_ndigi_active_before_limited"]
    res["vertex_ndigi_before_limited"] = res["vertex_ndigi_active_before_limited"] + res["vertex_ndigi_veto_before_limited"]
    print("Finished")

    return RQ_dict(res)