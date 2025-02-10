import copy

import numpy as np
import matplotlib.pyplot as plt

from . import util,parser

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

    def plot(self, ind_h, ind_v, scale = 0.001):
        colors=["cyan"]#, "yellow", "green", "red"]
        for i,vertex in enumerate(self.vertices):
            plt.plot(vertex[ind_h] * scale, vertex[ind_v] * scale, marker="*", color=colors[i%len(colors)])

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
    def __init__(self, data=None, i=0, metadata_digi=None):
        if data is not None:
            self.xyzt = np.array([data["Digi_x"][i], data["Digi_y"][i], data["Digi_z"][i], data["Digi_t"][i]])
            self.id = i
            self.pdg = data["Digi_pdgID"][i]
            self.is_primary = False
            self.hit_inds = data["Digi_hitInds_unpacked"][i]
            self.direction = data["Digi_direction"][i]
            self.type = data["Digi_type"][i]
            self.track_id = data["Digi_trackID"][i]
            
            if (self.type>0 and "Hit_x" in data):
                try:
                    ihit,jhit = self.hit_inds[0],self.hit_inds[1]
                except:
                    ihit = self.hit_inds[0]
                    jhit=ihit+1 if (ihit<len(data["Hit_y"])-2) else ihit-1
                self.momentum = np.array([data["Hit_px"][ihit], data["Hit_py"][ihit], data["Hit_pz"][ihit]])
                self.xyzt_truth = np.array([data["Hit_x"][ihit], data["Hit_y"][ihit], data["Hit_z"][ihit], data["Hit_t"][ihit]])
                self.velocity =(np.array([data["Hit_x"][ihit], data["Hit_y"][ihit], data["Hit_z"][ihit]]) -\
                                np.array([data["Hit_x"][jhit], data["Hit_y"][jhit], data["Hit_z"][jhit]])) / (data["Hit_t"][ihit] - data["Hit_t"][jhit])
            else:
                self.momentum = np.array([0,0,1])
                self.xyzt_truth = np.array([0,0,0,0])
                self.velocity = np.array([0,0,0,0])

            direction_x = round(self.direction//100 % 10)
            direction_y = round(self.direction//10 % 10)
            direction_z = 3-direction_x-direction_y
            # unc_xyzt = [metadata_digi["Uncertainty_x"], metadata_digi["Uncertainty_y"], metadata_digi["Uncertainty_z"], metadata_digi["Uncertainty_t"]]
            self.err4 = np.zeros(4)
            self.err4[direction_x] = metadata_digi["Uncertainty_x"]
            self.err4[direction_y] = metadata_digi["Uncertainty_y"]
            self.err4[direction_z] = metadata_digi["Uncertainty_z"]
            self.err4[3] = metadata_digi["Uncertainty_t"]

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

            self.convert_to_time()
            self.t0=data["Track_t0"][i]

    def convert_to_time(self):
        self.params_time = np.zeros(6);
        k=0
        for j in range(3):
            if (j == self.iv_index):
                self.params_time[j] = self.iv_value;
                self.params_time[j + 3] = 1 / self.params[-1];
            else:
                self.params_time[j] = self.params[k];
                self.params_time[j + 3] = self.params[k + 3] / self.params[-1];
                k+=1

        self.vabs = np.linalg.norm(self.params_time[3:])
    

    def set_digis(self, digis_all):
        ## Find the truth location
        digi = digis_all[self.hit_ids[-1]]
        digi_loc = digi.xyzt_truth
        digi_p = digi.momentum
        digi_v = digi.velocity
        kx, ky = digi_p[0]/digi_p[2], digi_p[1]/digi_p[2]
        kz = 1-np.sqrt(kx**2+ky**2)
        self.params_truth = np.concatenate((np.delete(digi_loc, self.iv_index), [kx,ky, 1/digi_v[2]]))

        ## Find the truth index and pdg 
        digi_trackid = []        
        digi_pdg = []
        digi_times = []
        for i in self.hit_ids:
            digi_trackid.append(digis_all[i].track_id)
            digi_pdg.append(digis_all[i].pdg)
            digi_times.append(digis_all[i].xyzt[-1])

        tid,tid_ind,tid_counts = np.unique(digi_trackid, return_counts=True, return_index=True)
        itrack = np.argmax(tid_counts)

        ## Calculate track purity etc.
        self.track_purity = tid_counts[itrack]/self.nhits
        self.track_id = tid[itrack]
        self.track_pdg = digi_pdg[tid_ind[itrack]]
        self.tmin = min(digi_times)
        self.tmax = max(digi_times)

    def plot(self, ind_h, ind_v, from_time = None, color="k", scale=0.001):
        tmin = self.tmin if from_time is None else from_time
        xs = np.array([self.params_time[ind_h] + self.params_time[ind_h+3]*(tmin-self.t0),
              self.params_time[ind_h] + self.params_time[ind_h+3]*(self.tmax-self.t0)])
        ys = np.array([self.params_time[ind_v] + self.params_time[ind_v+3]*(tmin-self.t0),
              self.params_time[ind_v] + self.params_time[ind_v+3]*(self.tmax-self.t0)] )
        plt.plot(xs * scale,ys * scale, linestyle=":", color=color)


class Vertex:
    def __init__(self, data=None, i=0):
        if data is not None:
            self.params = np.array([data["Vertex_x0"][i], data["Vertex_y0"][i], data["Vertex_z0"][i], data["Vertex_t0"][i]])
            self.cov = np.array(data["Vertex_cov_unpacked"][i])
            self.chi2 = data["Vertex_chi2"][i]
            self.track_ids =  data["Vertex_trackInds_unpacked"][i]
            self.ntracks = len(self.track_ids)
            self.vertex_ntracklet_0 = data["Vertex_tracklet_n0"][i]
            self.vertex_ntracklet_2 = data["Vertex_tracklet_n2"][i]
            self.vertex_ntracklet_3 = data["Vertex_tracklet_n3"][i] + data["Vertex_tracklet_n4p"][i]
    
    def set_tracks(self, tracks_all):
        tracks_purity = []
        for i in self.track_ids:
            t = tracks_all[i]
            tracks_purity.append(t.track_purity)

        self.vertex_purity = np.mean(tracks_purity)

class PDG:    
    name_map = {
         11: ["e-",],
         -11: ["e+",],
         13: [r"mu-",],
         -13: ["mu+",],
         22: [r"gamma",],
         211: ["pi+" ,],
         -211: ["pi-",],
         321: ["K+"  ,],
         -321: ["K-",],
         2112: ["n",],
         -2112: [r"n bar",],
         2212: ["p",],
         -2212: ["p bar",],
        }
    
    @staticmethod
    def str(pdg):
        if pdg in PDG.name_map:
            return PDG.name_map[pdg][0]
        return "pdg: " + str(pdg)          


class Event:
    def __init__(self, data_parsed, metadata_digi):
        ## Add some RRQ
        self.process_recon(data_parsed)
        
        # Get generator particles
        self.genparticles = np.array([GenParticle(data_parsed, i) for i in range(len(data_parsed["Gen_x"]))])
        self.genvertices = GenVertices(self.genparticles)
        
        # Get truth 
        if ("Hit_x" in data_parsed):
            self.HAS_TRUTH = True

            # Truth hits
            self.hits = np.array([Hit(data_parsed,i) for i in range(len(data_parsed["Hit_x"]))])
            # Truth track
            track_ids, hits_inds_grouped = util.Utils.groupby(data_parsed["Hit_trackID"], np.arange(len(data_parsed["Hit_x"])))
            self.truetracks = [TrueTrack(self.hits[hit_ids],track_id) for track_id, hit_ids in zip(track_ids, hits_inds_grouped)]
        else:
            self.HAS_TRUTH = False
        
        # Get Digits
        self.digis = np.array([Digi(data_parsed,i, metadata_digi) for i in range(len(data_parsed["Digi_x"]))])
        self.get_reconstructable()

        # Get recon tracks
        self.tracks = [Track(data_parsed, i) for i in range(len(data_parsed["Track_x0"]))]
        self.digis_used_inds = []
        for t in self.tracks:
            t.set_digis(self.digis)
            self.digis_used_inds.extend(list(t.hit_ids))

        # Get recon vertices
        self.vertices = [Vertex(data_parsed, i) for i in range(len(data_parsed["Vertex_x0"]))]
        self.tracks_used_inds = []
        for t in self.vertices:
            t.set_tracks(self.tracks)
            self.tracks_used_inds.extend(list(t.track_ids))

    def process_recon(self, data):
        data["Digi_hitInds_unpacked"] = parser.unpack_at(data["Digi_hitInds"], divider=-1)
        data["Track_cov_unpacked"] = parser.unpack_cov(data["Track_cov"], dim=6)
        data["Track_digiInds_unpacked"] = parser.unpack_at(data["Track_digiInds"], divider=-1)
        data["Vertex_cov_unpacked"] = parser.unpack_cov(data["Vertex_cov"], dim=4)
        data["Vertex_trackInds_unpacked"] = parser.unpack_at(data["Vertex_trackInds"], divider=-1)

    def get_reconstructable(self):
        ## Find the truth index and pdg 
        digi_trackid = []        
        digi_pdg = []
        for i in range(len(self.digis)):
            digi_trackid.append(self.digis [i].track_id)
            digi_pdg.append(self.digis [i].pdg)

        tid,tid_ind,tid_counts = np.unique(digi_trackid, return_counts=True, return_index=True)

        self.digi_truth_track_ids = tid
        self.digi_truth_track_nhits = tid_counts
    
    def plot_truetracks(self, ind_h, ind_v, scale = 0.001):
        if not self.HAS_TRUTH:
            return
        for t in self.truetracks:
            if t.plot_en:
                plt.plot(t.xyzts[ind_h]*scale, t.xyzts[ind_v]*scale, linewidth=1, alpha=0.9)  

    def plot_digis(self, ind_h, ind_v, scale = 0.001):
        # Plot used hits
        plot_xs = []
        plot_ys = []
        plot_xerrs = []
        plot_yerrs = []
        for d in self.digis:
            if d.type>=0 and d.id in self.digis_used_inds:
                plot_xs.append(d.xyzt[ind_h] * scale)
                plot_ys.append(d.xyzt[ind_v] * scale)
                plot_xerrs.append(d.err4[ind_h] * scale)
                plot_yerrs.append(d.err4[ind_v] * scale)   
        plt.errorbar(plot_xs, plot_ys, xerr=plot_xerrs, yerr=plot_yerrs, fmt="o", color="C0", markersize=3, alpha=0.3, capsize=2)

        # Plot unused hits
        plot_xs = []
        plot_ys = []
        plot_xerrs = []
        plot_yerrs = []
        for d in self.digis:
            if d.type>=0 and d.id not in self.digis_used_inds:
                plot_xs.append(d.xyzt[ind_h] * scale)
                plot_ys.append(d.xyzt[ind_v] * scale)
                plot_xerrs.append(d.err4[ind_h] * scale)
                plot_yerrs.append(d.err4[ind_v] * scale)   
        plt.errorbar(plot_xs, plot_ys, xerr=plot_xerrs, yerr=plot_yerrs, fmt="o", color="salmon", markersize=3, alpha=0.3, capsize=2)

        # Plot noise hits and other hit types
        plot_xs = []
        plot_ys = []
        plot_xerrs = []
        plot_yerrs = []
        for d in self.digis:
            if d.type==-1:
                plot_xs.append(d.xyzt[ind_h] * scale)
                plot_ys.append(d.xyzt[ind_v] * scale)
                plot_xerrs.append(d.err4[ind_h] * scale)
                plot_yerrs.append(d.err4[ind_v] * scale)   
        plt.errorbar(plot_xs, plot_ys, xerr=plot_xerrs, yerr=plot_yerrs, fmt="o", color="grey", markersize=3, alpha=0.3, capsize=2)

    def plot_vertex(self, ind_h, ind_v, plot_tracks=True, lim_x = 100e3, lim_y = 100e3, scale = 0.001):
        icolor=0
        for v in self.vertices:
            x,y = v.params[ind_h], v.params[ind_v]
            if abs(x)>lim_x:
                continue
            if abs(y)>lim_y:
                continue                
                
            plt.scatter(x * scale,y * scale, marker="^", color=f"C{icolor}")
            for track_id in v.track_ids:
                self.tracks[track_id].plot(ind_h,ind_v, from_time = v.params[3], color=f"C{icolor}", scale = scale)
            icolor +=1;

        for t in self.tracks:
            if t.id not in self.tracks_used_inds:
                t.plot(ind_h,ind_v, scale = scale)  

    def plot(self, figsize=(16, 8)):
        fig = plt.figure(figsize=figsize, layout="constrained")
        spec = fig.add_gridspec(4, 7)
        
        axfront = fig.add_subplot(spec[:2, :4])
        axside = fig.add_subplot(spec[2:, :4])
        axtop = fig.add_subplot(spec[:3, 4:])
        
        plt.sca(axfront)
        self.plot_truetracks(0,2)
        self.plot_digis(0,2)
        self.plot_vertex(0,2)
        self.genvertices.plot(0,2)
        plt.xlabel("x (beamline) [m]")
        plt.ylabel("z (verticle) [m]")
        
        plt.sca(axside)
        self.plot_truetracks(1,2)
        self.plot_digis(1,2)
        self.plot_vertex(1,2)
        self.genvertices.plot(1,2)
        plt.xlabel("y (other) [m]")
        plt.ylabel("z (verticle) [m]")
        
        plt.sca(axtop)
        self.plot_truetracks(0,1)
        self.plot_digis(0,1)
        self.plot_vertex(0,1)
        self.genvertices.plot(0,1)
        plt.xlabel("x (beamline) [m]")
        plt.ylabel("y (other) [m]")
        
        plt.show()    

                