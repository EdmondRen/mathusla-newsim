import copy

import numpy as np

from . import util

class GenParticle:
    def __init__(self, data=None, i=0):
        if data is not None:
            self.xyzt =  np.array([data["Gen_x"][i], data["Gen_y"][i], data["Gen_z"][i], data["Gen_t"][i]])
            self.momentum = np.array([data["Gen_px"][i], data["Gen_py"][i], data["Gen_pz"][i]])
            self.pdg =  data["Gen_pdgID"][i]

class GenVertices:
    def __init__(self, genparticles):
        p_x_unique_inds, unique_indices  = np.unique(np.array([p.xyzt[0] for p in genparticles]), return_index=True)
        self.vertices = [genparticles[i].xyzt for i in unique_indices]

    def plot(self, ind_h, ind_v):
        colors=["cyan", "yellow", "green", "red"]
        for i,vertex in enumerate(self.vertices):
            plt.plot(vertex[ind_h], vertex[ind_v], marker="D", color=colors[i%len(colors)])

class Hit:
    def __init__(self, data=None, i=0):
        if data is not None:
            self.xyzt = np.array([data["Hit_x"][i], data["Hit_y"][i], data["Hit_z"][i], data["Hit_t"][i]])
            self.momentum = np.array([data["Hit_px"][i], data["Hit_py"][i], data["Hit_pz"][i]])
            self.edep = data["Hit_edep"][i]
            self.id = i
            self.pdg = data["Hit_pdgID"][i]
            self.is_primary = False #data["Hit_isprimary"][i]

class Digi:
    def __init__(self, data=None, i=0):
        if data is not None:
            self.xyzt = np.array([data["Digi_x"][i], data["Digi_y"][i], data["Digi_z"][i], data["Digi_t"][i]])
            self.id = i
            self.pdg = data["Digi_pdgID"][i]
            self.is_primary = False
            self.hit_inds = data["Digi_hitInds_unpacked"][i]
            self.direction = data["Digi_direction"][i]
            self.type = data["Digi_type"][i]
            

            try:
                ihit,jhit = self.hit_inds[0],self.hit_inds[1]
            except:
                ihit = self.hit_inds[0]
                jhit=ihit+1 if (ihit<len(data["Hit_y"])-2) else ihit-1
            self.momentum = np.array([data["Hit_px"][ihit], data["Hit_py"][ihit], data["Hit_pz"][ihit]])
            self.xyzt_truth = np.array([data["Hit_x"][ihit], data["Hit_y"][ihit], data["Hit_z"][ihit], data["Hit_t"][ihit]])
            self.velocity =(np.array([data["Hit_x"][ihit], data["Hit_y"][ihit], data["Hit_z"][ihit]]) -\
                            np.array([data["Hit_x"][jhit], data["Hit_y"][jhit], data["Hit_z"][jhit]])) / (data["Hit_t"][ihit] - data["Hit_t"][jhit])
            
            # self.err = np.array([data["Hit_x"][i], data["Hit_y"][i], data["Hit_z"][i], data["Hit_t"][i]])

class TrueTrack:
    def __init__(self, hits, track_id=0, cut_timerange=0.1, cut_nhits = 3, cut_momentum = 10):
        self.hits = hits
        self.pdg = hits[0].pdg
        self.id = track_id
        self.is_primary = hits[0].is_primary

        self.xyzts = np.array([hit.xyzt for hit in hits]).T
        self.ps = np.array([np.linalg.norm(hit.momentum) for hit in hits])

        self.plot_en = True
        if (max(self.xyzts[:,-1]) - min(self.xyzts[:,-1])) < cut_timerange or sum(self.ps>cut_momentum)<cut_nhits:
            self.plot_en=False

class Track:
    def __init__(self, data=None, i=0):
        if data is not None:
            self.params_full = np.array([data["Track_x0"][i], data["Track_y0"][i], data["Track_z0"][i], data["Track_t0"][i], 
                                         data["Track_kx"][i], data["Track_ky"][i], data["Track_kz"][i], data["Track_kt"][i]])
            self.cov = np.array(data["Track_cov_unpacked"][i])
            self.chi2 = data["Track_chi2"][i]
            self.iv_index =  data["Track_iv_ind"][i]
            self.iv_value =  self.params_full[self.iv_index]
            self.iv_error =  data["Track_iv_err"][i]
            self.id =  i
            self.hit_ids =  data["Track_digiInds_unpacked"][i]
            self.params = np.delete(self.params_full, [self.iv_index, self.iv_index+4])
            self.nhits = len(self.hit_ids)


class Vertex:
    def __init__():
        pass


class Event:
    def __init__(self, data_parsed):
        # Get generator particles
        self.genparticles = np.array([GenParticle(data_parsed, i) for i in range(len(data_parsed["Gen_x"]))])
        self.genvertices = GenVertices(self.genparticles)
        
        # Get truth hits
        self.hits = np.array([Hit(data_parsed,i) for i in range(len(data_parsed["Hit_x"]))])
        
        # Get truth track
        track_ids, hits_inds_grouped = util.Utils.groupby(data_parsed["Hit_trackID"], np.arange(len(data_parsed["Hit_x"])))
        self.truetracks = [TrueTrack(self.hits[hit_ids],track_id) for track_id, hit_ids in zip(track_ids, hits_inds_grouped)]
        
        # Get Digits
        self.digis = np.array([Digi(data_parsed,i) for i in range(len(data_parsed["Digi_x"]))])

        # Get recon tracks
        self.tracks = [Track(data_parsed, i) for i in range(len(data_parsed["Track_x0"]))]

        # Get recon vertices
    
    def truetracks_plot(self, ind_h, ind_v):
        for t in self.truetracks:
            if t.plot_en:
                plt.plot(t.xyzts[ind_h], t.xyzts[ind_v], linewidth=1)  