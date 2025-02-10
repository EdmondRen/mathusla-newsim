#ifndef tracker_types_HH
#define tracker_types_HH

#include <cstdio>
#include <iostream>
#include <memory>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <TFile.h>
#include <TTree.h>

#include "util_root.hh"

using Eigen::MatrixXd;
using Eigen::Vector3d;
using Eigen::Vector4d;
using Eigen::VectorXd;
using MatrixXdSp = Eigen::SparseMatrix<double>;

namespace Tracker
{

    // Helper functions
    class Helper
    {
    public:
        Helper() {}
        // Remove i^th element from a vector and creat a new one.
        static Eigen::VectorXd removeElement(const Eigen::VectorXd &vec4, int i);
        // Inseart value to a vector at given index i
        static Eigen::VectorXd insertVector(const Eigen::VectorXd &vec4, int i, float val);
        // Function to insert a row at a specific index
        static Eigen::MatrixXd insertRow(const Eigen::MatrixXd &mat, const Eigen::RowVectorXd &row, int insert_pos);
        // Function to insert a column at a specific index
        static Eigen::MatrixXd insertColumn(const Eigen::MatrixXd &mat, const Eigen::VectorXd &col, int insert_pos);
        // Function to insert a row and column of zeros at the i-th index
        static Eigen::MatrixXd insertRowAndColumnOfZeros(const Eigen::MatrixXd &mat, int i);
    };

    class DigiHit
    {
    public:
        DigiHit(double _x, double _y, double _z, double _t,
                double _ex, double _ey, double _ez, double _et,
                int _iv_index, int _layer, int _group, int _id) : iv_index(_iv_index), layer(_layer), group(_group), id(_id)
        {
            vec4 << _x, _y, _z, _t;
            vec4_err << _ex, _ey, _ez, _et;
            update_vec3();
        }

        // Update 3-vector
        // Must be called after changing the vec4 value
        void update_vec3()
        {
            vec3 = Helper::removeElement(vec4, iv_index);
            vec3_err = Helper::removeElement(vec4_err, iv_index);
        }

        // Required data
        Vector4d vec4, vec4_err; // {x,y,z,t} and their uncertainty
        Vector3d vec3, vec3_err; // Independent variable removed
        const int iv_index;      // Index of the independent variable. one of {0,1,2,3}. Use it to decide which one of x,y,z,t is the independent variable.
        int layer;               // Which detector layer this hit belongs to
        int group;               // Which detector grop this hit belongs to

        // Optional
        int id;
        int type;
        int64_t detector_id;

        // Set and Get methods
        inline double x() const { return vec4(0); }
        inline double y() const { return vec4(1); }
        inline double z() const { return vec4(2); }
        inline double t() const { return vec4(3); }
        inline void setx(float x) { vec4(0) = x; }
        inline void sety(float x) { vec4(1) = x; }
        inline void setz(float x) { vec4(2) = x; }
        inline void sett(float x) { vec4(3) = x; }
        inline double ex() const { return vec4_err(0); }
        inline double ey() const { return vec4_err(1); }
        inline double ez() const { return vec4_err(2); }
        inline double et() const { return vec4_err(3); }
        inline void setex(float x) { vec4_err(0) = x; }
        inline void setey(float x) { vec4_err(1) = x; }
        inline void setez(float x) { vec4_err(2) = x; }
        inline void setet(float x) { vec4_err(3) = x; }
        inline double get_step() const { return vec4(iv_index); } // std::cout<<vec4<<" index: "<<(iv_index)<<std::endl;
        inline double get_steperr() const { return vec4_err(iv_index); }
        inline Vector3d get_vec3() const { return vec3; }
        inline Vector3d get_err3() const { return vec3_err; }
    };
    using HitList = std::vector<std::unique_ptr<DigiHit>>;

    class Track
    {
    public:
        Track() {}

        // Required data
        VectorXd params; // {x0, y0, z0, t0, Ax, Ay, Az, At}, except for the ones that is used as independent variable
        MatrixXd cov;    // covariance of {x0, y0, z0, t0, Ax, Ay, Az, At}, except for the ones that is used as independent variable
        float chi2;
        int iv_index; // Index of the independent variable. one of {0,1,2,3}. Use it to decide which one of x,y,z,t is the independent variable.
        int iv_value; // Value of the independent variable
        int iv_error; // Uncertainty of the independent variable

        // Optional data
        int id;
        std::vector<int> hit_ids;

        // Multiple-scattering matrix block
        MatrixXd Q_block;

        // Parameters based on time
        float t0, t0_err;
        VectorXd params_time; // {x0, y0, z0, vx, vy, vz}
        MatrixXd cov_time;    // covariance of {x0, y0, z0, vx, vy, vz}

        // Full 8-D parameters
        VectorXd params_full; // {x0, y0, z0, t0, Ax, Ay, Az, At}
        MatrixXd cov_full;    // covariance of {x0, y0, z0, t0, Ax, Ay, Az, At}

        // Convert to using time as independent variable
        VectorXd convert_to_time();

        // Calculate closest-approach midpoint and distance
        // Return: pair of <closest distance, {x,y,z,t}>
        static std::pair<double, Vector4d> get_closest_midpoint(const Track &track1, const Track &track2);

        // Get the distance at the same time of a given point
        double get_same_time_dist(Vector4d point) const;

        // Get the <dist, chi2> at the same time of a given point
        std::pair<Vector4d, MatrixXd> get_same_time_pos_and_cov(Vector4d point, double speed_constraint = -1.0, bool multiple_scattering = true) const;
        std::pair<double, double> get_same_time_dist_and_chi2(Vector4d point, double speed_constraint = -1.0, bool multiple_scattering = true) const;

        // Get the <dist, chi2> at the same value of independent variable
        std::pair<Vector4d, MatrixXd> get_same_invar_pos_and_cov(Vector4d point, double speed_constraint = -1.0, bool multiple_scattering = true) const;
        std::pair<double, double> get_same_invar_dist_and_chi2(Vector4d point, double speed_constraint = -1.0, bool multiple_scattering = true) const;

        // Get the <dist, cov> at the closest point of approach (CPA) on the track to a given point
        std::pair<Vector4d, MatrixXd> get_cpa_pos_and_cov(Vector4d point, double speed_constraint = -1.0, bool multiple_scattering = true) const;
        std::pair<double, double> get_cpa_dist_and_chi2(Vector4d point, double speed_constraint = -1.0, bool multiple_scattering = true) const;
    };
    using TrackList = std::vector<std::unique_ptr<Track>>;

    class TrackSeed
    {
    public:
        TrackSeed(DigiHit *hit1, DigiHit *hit2) : nhits_found(-1)
        {
            float c = 299.7; // speed of light [mm/ns]
            this->dvec = hit2->vec4.array() - hit1->vec4.array();
            this->dstep = std::abs(hit2->get_step() - hit1->get_step());
            this->dr = dvec.segment(0, 3).norm();
            this->dt = std::abs(dvec(3));
            this->score = std::abs(dr / c - dt);

            // Make a pair and put the earlier one in front
            this->hits = hit1->t() < hit2->t() ? std::make_pair(hit1, hit2) : std::make_pair(hit2, hit1);

            this->chi2prob_found = 1; // Default value of chi2 prob if not specified.
        }

        // Required data
        std::pair<DigiHit *, DigiHit *> hits;
        float score;
        VectorXd dvec;
        float dr, dt, dstep;
        int nhits_found;       // Number of hits found using this seed. Save time when reusing this seed.
        double chi2prob_found; //
        int nhit_occur;

        // float GetScore() { return score; }
    };

    class Vertex
    {
    public:
        Vertex() {}

        // Required data
        VectorXd params; // {x0, y0, z0, t0}
        MatrixXd cov;    // covariance of {x0, y0, z0, t0}
        float chi2;

        // Optional
        int id;
        std::vector<int> track_ids;
        std::vector<int> hit_ids;
    };
    using VertexLilst = std::vector<std::unique_ptr<Vertex>>;

    // Read the file after digitization, and arrange them into <DigiHit> type
    class TreeReaderDigi
    {
    private:
        TFile *File;
        TTree *TreeData;
        TTree *TreeMetadata;

        // Buffer for digitized data to read
        std::vector<float> *Digi_x;
        std::vector<float> *Digi_y;
        std::vector<float> *Digi_z;
        std::vector<float> *Digi_t;
        std::vector<float> *Digi_edep;
        std::vector<int> *Digi_trackID;
        std::vector<int> *Digi_pdgID;
        std::vector<long long> *Digi_detectorID; // Which layer the hit is from. Layer is obtained from the copy number of depth 1 in GENAT4
        std::vector<int> *Digi_type;             // Soure of the event. -1: noise, 0: GUN, 1: PARMA, 2: CRY
        std::vector<int> *Digi_hitInds;          // The index of truth hits of each digitized hit
        std::vector<int> *Digi_direction;        // Indicates the direction of the bar with last three digits . For example, .....012 means x->x, y->y, z->z

        // Buffer for metadata;
        std::string *SimulationName;
        std::string *Geometry;
        std::string *Generator;
        float Uncertainty_t;
        float Uncertainty_x;
        float Uncertainty_y;
        float Uncertainty_z;

        // Container for processed data
        Vector4d digi_unc;
        std::vector<std::unique_ptr<DigiHit>> hits;

        // Counter of current entry
        long long entries_total;
        long long entry_current;

    public:
        TreeReaderDigi(std::string filename);
        ~TreeReaderDigi()
        {
            File->Close();
        }

        std::vector<std::unique_ptr<DigiHit>> GetEntry(long long entry);

        // Group hits by detector copy number
        // detID has format AAABBBCCCDDDD, in which BBB is the detector copy number
        std::unordered_map<int, std::vector<DigiHit *>> ProcessHits(std::vector<std::unique_ptr<DigiHit>> &hits);

        long long GetCurrentEntry() { return entry_current; }
        long long GetEntries() { return entries_total; }
        TTree *GetdataTree() { return TreeData; }
        TTree *GetmetadataTree() { return TreeMetadata; }
    };

    // Write output file
    //  with option to merge the digi and raw data together
    class TreeWriterRecon : public iroot::file::EntryCopy
    {
    private:
        TFile *outputFile;
        TTree *outputTree;
        TTree *outputTreeMetadata;

        TFile *digiFile;
        TTree *digiTree;
        TTree *digiTreeMetadata;

        TFile *simFile;
        TTree *simTree;
        TTree *simTreeMetadata;

        // Buffer for recon data to write
        int SimEntry;
        std::vector<float> Track_x0;
        std::vector<float> Track_y0;
        std::vector<float> Track_z0;
        std::vector<float> Track_t0;
        std::vector<float> Track_kx;
        std::vector<float> Track_ky;
        std::vector<float> Track_kz;
        std::vector<float> Track_kt;
        std::vector<float> Track_cov;
        std::vector<float> Track_chi2;
        std::vector<int> Track_id;
        std::vector<int> Track_iv_ind;
        std::vector<int> Track_iv_err;
        std::vector<int> Track_digiInds;

        std::vector<float> Vertex_x0;
        std::vector<float> Vertex_y0;
        std::vector<float> Vertex_z0;
        std::vector<float> Vertex_t0;
        std::vector<float> Vertex_cov;
        std::vector<float> Vertex_chi2;
        std::vector<int> Vertex_id;
        std::vector<int> Vertex_trackInds;
        std::vector<int> Vertex_tracklet_n0;
        std::vector<int> Vertex_tracklet_n2;
        std::vector<int> Vertex_tracklet_n3;
        std::vector<int> Vertex_tracklet_n4p;

        // Utility for copying raw and digits to the recon file
        // iroot::file::EntryCopy *copier_raw, *copier_digi;

        // Counter of current entry
        long long entries_total;
        long long entry_current;

        // Boolean
        bool EN_COPY_RAW;
        bool EN_COPY_DIGI;

    public:
        TreeWriterRecon(std::string filename_recon, std::string filename_digi = "", std::string filename_sim = "", bool save_raw_reduced = false);
        ~TreeWriterRecon()
        {
            outputFile->Close();
        }

        void SetSimBranches(bool save_raw_reduced);

        int ApplyRecon(TrackList &tracks, VertexLilst &vertices, std::vector<std::unordered_map<int,int>> &track_statsm, int &simulation_entry_number);
        int ApplyCopy(long long entry);
        void Fill();
        void Clear();
        void Write()
        {
            outputFile->cd();
            outputTree->Write();
            outputTreeMetadata->Write();
            if (EN_COPY_DIGI)
                digiTreeMetadata->CloneTree()->Write();
            if (EN_COPY_RAW)
                simTreeMetadata->CloneTree()->Write();
        }
        void Close() { outputFile->Close(); }

        // Buffer for metadata
        //   Put this in public. Too troublesome to write set for each of them
        void FillMetadata() { outputTreeMetadata->Fill(); }
        std::string meta_ReconstructionConfigStr;
    };

} // namespace Tracker

#endif
