import os
import random
from pathlib import Path

import numpy as np
import joblib
import pylorentz

from . import root as iroot
from . import util

def make_event(i, particles, vertex=None):
    """
    INPUT:
        i: event number
        particles: list of particles
            each particle is [PDG, x, y, z, px, py, pz, t]
        vertex: none or vertex
            vertex is [PDG, x, y, z, px, py, pz, t]
            
    NOTE:
        unit is [mm], [MeV]
    RETURN:
        formated string of this event
    """
    
    vertex = [0,0,0,0,0,0,0,0] if vertex is None else vertex # [PID, x,y,z, px,py,pz, t]
    
    evt= f'E {i} {vertex[0]}  {vertex[1]}  {vertex[2]}  {vertex[3]}  {vertex[4]}  {vertex[5]}  {vertex[6]}  {vertex[7]}\n'
    for p in particles:
        evt+=f'P {p[0]} {p[1]} {p[2]} {p[3]} {p[4]} {p[5]} {p[6]} {p[7]}\n'
        
    return evt

def open_file(filename, nevents):
    f = open(filename_filereader, "w")
    f.write(f"# nevents {nevents}\n\n")
    return f



class unit:
    m=100
    cm=1
    

def frame_transform(four_vector_in, boost_direction, beta=None, gamma=None):
    """
    Transform the four-vector (e, p_x, p_y, p_z,) with boost_direction and beta/gamma.
    
    INPUT:
    ---
    e, p_x, p_y, p_z: 
        four-vector to be transformed. Can be momentum or what ever 4-vector
    boost_direction:
        three-vector, the relative velocity of the frame to be transformed to
    beta/gamma:
        at least one needs to be given
        
    RETURN
    ---
    four_vector_out:
        the transformed four-vector
    """
    x, y, z = boost_direction
    e, p_x, p_y, p_z = four_vector_in
    four_vector_out = pylorentz.Momentum4(e, p_x, p_y, p_z).boost(x, y, z, beta=beta, gamma=gamma)
    return four_vector_out
    
    
def twobody_decay(mass, three_momentum, decayproduct_pid=[13,13], rand_seed=None):
    """
    INPUT
    ---
    mass: 
        float, mass of the primary particle [MeV]
    three_momentum: 
        list, three momentum of the primary partile [MeV]
    decayproduct_pid:
        [pid1, pid2], PID of the two decay products. Used to decide the mass of the decay particle.
    rand_seed: 
        If None, then fresh, unpredictable entropy will be pulled from the OS. 
    """
    mass_map={13: 105.7,    # mu-
            11:0.511,       # e-
              211: 139.57,  # pi+/-
              -211: 139.57, # pi+/-
             }
    
    m1=mass_map[decayproduct_pid[0]]
    m2=mass_map[decayproduct_pid[1]]
    M=mass
    P=np.sqrt(np.sum(np.power(three_momentum,2)))
    
    # First, 0->1+2 decay in the rest frame of 0:
    E1 = (M**2+m1**2-m2**2)/(2*M)
    E2 = (M**2+m2**2-m1**2)/(2*M)
    p1 = ((M**2+m1**2-m2**2)**2- 4*M**2*m1**2)**0.5/(2*M)
    p2 = ((M**2+m2**2-m1**2)**2- 4*M**2*m2**2)**0.5/(2*M)
    
    # Generate a randomized direction
    rng = np.random.default_rng(seed=rand_seed)
    p1_x = rng.normal(0,1)
    p1_y = rng.normal(0,1)
    p1_z = rng.normal(0,1)
    mag = np.sqrt(p1_x**2+p1_y**2+p1_z**2)
    p1_x = p1_x/mag*p1
    p1_y = p1_y/mag*p1
    p1_z = p1_z/mag*p1
    
    # Boost the 4-momentum
    boost_direction = -np.array(three_momentum)
    boost_gamma = np.sqrt(M**2+P**2)/M
    vec1 = frame_transform([E1, p1_x, p1_y, p1_z], boost_direction, gamma=boost_gamma)
    vec2 = frame_transform([E2, -p1_x, -p1_y, -p1_z], boost_direction, gamma=boost_gamma)
    
    return vec1, vec2

def twobody_decay_MA(mass, abs_momentum, vertex_xyz, decayproduct_pid=[13,13], rand_seed=None):
    """
    INPUT
    ---
    mass: 
        float, mass of the primary particle [MeV]
    abs_momentum: 
        float, absolute momentum of the primary partile [MeV]
    vertex_xyz: in uni of [cm]
        [x,y,z] in **CMS** coordinates! It's the same coordinates as the HIT_x,y,z in simulation output.
    decayproduct_pid:
        [pid1, pid2], PID of the two decay products. Used to decide the mass of the decay particle.
    rand_seed: 
        If None, then fresh, unpredictable entropy will be pulled from the OS. 
    
    RETURN
    ---
    p4vec_1,p4vec_2:
        Momentum four-vector of the two decay products
    
    """
    r = np.sqrt(np.sum(np.power(vertex_xyz,2)))
    px = abs_momentum/r*vertex_xyz[0]
    py = abs_momentum/r*vertex_xyz[1]
    pz = abs_momentum/r*vertex_xyz[2]
    
    p4vec_1,p4vec_2 = twobody_decay(mass, [px,py,pz], decayproduct_pid=decayproduct_pid, rand_seed=rand_seed) 
    
    return p4vec_1,p4vec_2

def multibody_decay_MA(mass, abs_momentum, vertex_xyz, decayproduct_p4vec_list):
    boost_direction = -np.array(three_momentum)
    boost_gamma = np.sqrt(mass**2+abs_momentum**2)/mass
    
    p4vec_transformed=np.array([frame_transform(particle, boost_direction, gamma=boost_gamma) for particle in decayproduct_p4vec_list])
    return p4vec_transformed 

def generate_sim_script_filereader(events_properties_filename, script_path=None):
    """
    Create simulation script for filereader generator based on the events database file.
    """
    if script_path is None:
        script_path = os.path.splitext(events_properties_filename)[0]+".mac"
        
    with open(events_properties_filename, 'r') as file:
        nlines_read=0
        found_nevents=False
        while (not found_nevents) and nlines_read<=10:
            line_content = file.readline().split()
            nlines_read+=1
            if "nevents" in line_content:
                try:
                    nevents = int(line_content[-1])
                    found_nevents = True
                    break
                except:
                    print("Could not read number of events")
            

    script = "/det/select Box \n"
    script+= "/gen/select file_reader \n"
    script+= f"/gen/file_reader/pathname {events_properties_filename}\n"
    script+= f"/run/beamOn {nevents}"

    with open(script_path, 'w') as file:
        file.write(script)
        
    print("Script saved at", script_path)
    
    return script_path

def generate_twobody_decay_file(filename, mass, abs_momentum, vertex_xyz, decayproduct_pid=[13,13], rand_seed=None, Nevents=10000, OVER_WRITE=False, which_coordinate = "CMS"):
    """
    INPUT
    ---
    filename:
        str, the filename of the event description file to be generated
    mass: 
        float, mass of the primary particle [MeV]
    abs_momentum: 
        float, absolute momentum of the primary partile [MeV]
    vertex_xyz: [x,y,z] or [[x0,x1], [y0,y1], [z0,z1]], in unit of [cm]
        [x,y,z] in **CMS** coordinates! It's the same coordinates as the HIT_x,y,z in simulation output.
        [[x0,x1], [y0,y1], [z0,z1]]: ranges of xyz, will generate uniform distributed xyz. 
    decayproduct_pid:
        [pid1, pid2], PID of the two decay products. Used to decide the mass of the decay particle.
    rand_seed: 
        If None, then fresh, unpredictable entropy will be pulled from the OS. 
    which_coordinate:
        "CMS" or "detector"
    """    
    
    if os.path.exists(filename):
        print("File exists!")
        if not OVER_WRITE:
            print("Please change filename, or set OVER_WRITE=True")
            return
        
    # Transform vertex position to detector coordinate to feed to simulation
    vertex_xyz = np.array(vertex_xyz)
    if vertex_xyz.ndim==1:
        if which_coordinate=="CMS":
            vertex_xyz_det = np.array([vertex_xyz[2],vertex_xyz[0],-vertex_xyz[1]+85.47*unit.m])*10 # turn into mm
            vertex_xyz_cms = vertex_xyz*10
        else:
            vertex_xyz_det = vertex_xyz*10
            vertex_xyz_cms = np.array([vertex_xyz[1],-vertex_xyz[2]+85.47*unit.m,vertex_xyz[0]])*10
        
    
    rng = np.random.default_rng(seed=rand_seed)
    
    
    with open(filename, "w") as file:
        
        
        # first, write the total number of events
        primary_particle_PID = -1

        file.write(f"# nevents {Nevents}\n\n")
        for i in range(Nevents):
            # Generate a random vertex if a range was given
            if vertex_xyz.ndim==2:
                vertex_temp = np.array([rng.uniform(*vertex_xyz[0]), rng.uniform(*vertex_xyz[1]), rng.uniform(*vertex_xyz[2])])
                if which_coordinate=="CMS":
                    # Convert to dtector coordinate and # turn into mm
                    vertex_xyz_det = np.array([vertex_temp[2],vertex_temp[0],-vertex_temp[1]+85.47*unit.m])*10 
                    vertex_xyz_cms = vertex_temp*10
                    
                else:
                    vertex_xyz_det = vertex_temp*10
                    vertex_xyz_cms = np.array([vertex_temp[1],-vertex_temp[2]+85.47*unit.m,vertex_temp[0]])*10
                 
                
            # Calculate the momentum of the parent particle
            r = np.sqrt(np.sum(np.power(vertex_xyz_cms,2)))
            primary_particle_px = abs_momentum/r*vertex_xyz_cms[2]
            primary_particle_py = abs_momentum/r*vertex_xyz_cms[0]
            primary_particle_pz = -abs_momentum/r*vertex_xyz_cms[1]     
            
            # Write parent particle
            file.write(f"n {i}  \t {primary_particle_PID}\t  0.0     0.0    0.0    {primary_particle_px}\t {primary_particle_py}\t {primary_particle_pz}\n")
            rand_seed_i = None if rand_seed is None else rand_seed+i

                    
            # Decay products
            p4vec_1,p4vec_2 = twobody_decay_MA(mass, abs_momentum, vertex_xyz_cms, decayproduct_pid=decayproduct_pid, rand_seed=rand_seed_i)
            # print(vertex_xyz_cms,p4vec_1,p4vec_2)
            # p4vec_1 = [E, px, py,pz] in CMS coordinate. To traslate it to GEANT coordinate: 
            file.write(f"\t {decayproduct_pid[0]}\t  {vertex_xyz_det[0]}\t  {vertex_xyz_det[1]}\t {vertex_xyz_det[2]}\t {p4vec_1[3]}\t {p4vec_1[1]}\t {-p4vec_1[2]}\n")
            file.write(f"\t {decayproduct_pid[1]}\t  {vertex_xyz_det[0]}\t  {vertex_xyz_det[1]}\t {vertex_xyz_det[2]}\t {p4vec_2[3]}\t {p4vec_2[1]}\t {-p4vec_2[2]}\n")
            
            
            
            
def generate_twobody_vertex_range(filename, abs_momentum, vertex_xyz, theta_range, p_unit_pre=[1,1,-1], phi_range=[0,2*np.pi], decayproduct_pid=[13,13], rand_seed=None, Nevents=10000, OVER_WRITE=False):
    """
    INPUT
    ---
    filename:
        str, the filename of the event description file to be generated
    abs_momentum: 
        float, absolute momentum of the two particles [MeV]
    vertex_xyz: [x,y,z] or [[x0,x1], [y0,y1], [z0,z1]], in unit of [cm]
        [x,y,z] in **DETECTOR** coordinates! It's the same coordinates as the x,y,z of the 'Range' generator in simulation output.
        [[x0,x1], [y0,y1], [z0,z1]]: ranges of xyz, will generate uniform distributed xyz. 
    theta_range:
        [theta_low, theta_high]
    decayproduct_pid:
        [pid1, pid2], PID of the two decay products. Used to decide the mass of the decay particle.
    rand_seed: 
        If None, then fresh, unpredictable entropy will be pulled from the OS. 
    which_coordinate:
        "CMS" or "detector"
    """    
    
    if os.path.exists(filename):
        print("File exists!")
        if not OVER_WRITE:
            print("Please change filename, or set OVER_WRITE=True")
            return
        
    # Transform vertex position to detector coordinate to feed to simulation
    vertex_xyz = np.array(vertex_xyz)
    vertex_xyz_det = vertex_xyz*10
    vertex_xyz_cms = np.array([vertex_xyz[1],-vertex_xyz[2]+85.47*unit.m,vertex_xyz[0]])*10
        
    rng = np.random.default_rng(seed=rand_seed)
    with open(filename, "w") as file:
        # first, write the total number of events
        file.write(f"# nevents {Nevents}\n\n")
        
        primary_particle_PID = -1
        r = np.sqrt(np.sum(np.power(vertex_xyz_cms,2)))
        primary_particle_px = abs_momentum/r*vertex_xyz_cms[0]
        primary_particle_py = abs_momentum/r*vertex_xyz_cms[1]
        primary_particle_pz = abs_momentum/r*vertex_xyz_cms[2]     
        
        
        for i in range(Nevents):
            file.write(f"n {i}  \t {primary_particle_PID}\t  {vertex_xyz_det[0]}\t  {vertex_xyz_det[1]}\t {vertex_xyz_det[2]}\t {primary_particle_px}\t {primary_particle_py}\t {primary_particle_pz}\n")
            
            rand_seed_i = None if rand_seed is None else rand_seed+i
            
            pvecs=[]
            for ivec in range(2):
                theta = rng.uniform(*theta_range)
                phi = rng.uniform(*phi_range)
                p_unit=np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), -np.cos(theta)])
                p_unit=p_unit/np.sqrt(np.sum(p_unit**2))*p_unit_pre
                p = abs_momentum*p_unit         
                pvecs.append(p)
            
            # Generate a random vertex if a range was given
            if vertex_xyz.ndim==2:
                vertex_temp = np.array([rng.uniform(*vertex_xyz[0]), rng.uniform(*vertex_xyz[1]), rng.uniform(*vertex_xyz[2])])
                vertex_xyz_det = vertex_temp*10
                    

            file.write(f"\t {decayproduct_pid[0]}\t  {vertex_xyz_det[0]}\t  {vertex_xyz_det[1]}\t {vertex_xyz_det[2]}\t {pvecs[0][0]}\t {pvecs[0][1]}\t {-pvecs[0][2]}\n")
            file.write(f"\t {decayproduct_pid[0]}\t  {vertex_xyz_det[0]}\t  {vertex_xyz_det[1]}\t {vertex_xyz_det[2]}\t {pvecs[1][0]}\t {pvecs[1][1]}\t {-pvecs[1][2]}\n")            
            
            

# -----------------------------------------------------------------------------------
# LLP related

            
def line_cube_intersect(line_P0, line_P1,cube_vertex_xyz, cube_length_xyz, cube_ex, cube_ey, cube_ez):
    """
    Find the intersection of a line and a cube
    
    Line is defined as r(t) = line_P0 + line_P1*t
    Cube is given with 5 variables:
     cube_vertex_xyz: [x0,y0,z0], One vertex of the cube
     cube_length_xyz: [lenth_x,lenth_y,lenth_z]
     cube_ex, cube_ey,cube_ez: [ex1,ey1,ez1], unit vectors of the three sides of the cube
     
    Return:
    ---
    v1,v2: coordinates of the two intersect point.
    If v1==v2==0, there is no intersection.
    
    Test:
    ---
    line_cube_intersect([0,0,0], [2,0,1], [2,2,0],[4,1,2], [-1,0,0],[0,-1,0],[0,0,1] )
    
    """
    
    P0 = np.array(line_P0)
    P1 = np.array(line_P1)
    A = np.array(cube_vertex_xyz)
    
    dot1 = P1@cube_ex
    dot2 = (A-P0)@cube_ex
    if dot1!=0:
        tmin_x = dot2/dot1     if dot1>0 else (cube_length_xyz[0]+dot2)/dot1
        tmax_x = (cube_length_xyz[0]+dot2)/dot1 if dot1>0 else dot2/dot1
    else:
        tmin_x=tmax_x=np.inf
    
    dot1 = P1@cube_ey
    dot2 = (A-P0)@cube_ey
    if dot1!=0:
        tmin_y = dot2/dot1     if dot1>0 else (cube_length_xyz[1]+dot2)/dot1
        tmax_y = (cube_length_xyz[1]+dot2)/dot1 if dot1>0 else dot2/dot1
    else:
        tmin_y=tmax_y=np.inf        
    
    dot1 = P1@cube_ez
    dot2 = (A-P0)@cube_ez
    if dot1!=0:
        tmin_z = dot2/dot1     if dot1>0 else (cube_length_xyz[2]+dot2)/dot1
        tmax_z = (cube_length_xyz[2]+dot2)/dot1 if dot1>0 else dot2/dot1    
    else:
        tmin_z=tmax_z=np.inf        
        
    tmins = np.array([tmin_x,tmin_y,tmin_z])
    tmaxs = np.array([tmax_x,tmax_y,tmax_z])
    tmin  = np.max(tmins[abs(tmins)!=np.inf])
    tmax  = np.min(tmaxs[abs(tmaxs)!=np.inf])
    
    if tmax<tmin:
        return np.array([0,0,0]),np.array([0,0,0])
    else:
        return P0+tmin*P1, P0+tmax*P1
    
# def read_decay_file(filename_products):
#     llp_decay=[]
#     new_decay=[]
#     with open(filename_products ,"r") as f:
#         for line in f:
#             line = line.split(",")
#             # Attach a new event
#             if len(line)==6:
#                 llp_decay.append([])
#             elif len(line)==7:
#                 llp_decay[-1].append([float(line[0]),float(line[1]),float(line[2]),float(line[3]),int(float(line[-2]))])            
                
#     return llp_decay  

def read_decay_file(filename_products):
    """
    Return a list of decay products. Each item is a list of one decay:
    decay 1:
    [particle1:[E, px, py, pz, pid],
    particle1:[E, px, py, pz, pid],
    particle1:[E, px, py, pz, pid],
    ...]
    """
    
    llp_decay=[]
    new_decay=[]
    # A simple state machine:
    # 0: empty line again
    # 1: Make a new particle (header line of new particle)
    # 2: body lines of new particle
    # 3: finish the particle
    pre_state = 0
    current_state = 0
    
    with open(filename_products ,"r") as f:
        for i,line in enumerate(f):
            if i==0:
                continue
                
            line = line.split(",")
            
            # State machine
            if pre_state==0:
                if len(line)<6:
                    current_state=0
                else:
                    current_state=1
            elif pre_state==1:
                if len(line)>0:
                    current_state=2
                else:
                    current_state=0
            elif pre_state==2:
                if len(line)<6:
                    current_state=3
                else:
                    current_state=2
            elif pre_state==3:
                if len(line)<6:
                    current_state=0
                else:
                    current_state=1
            
            
            if current_state==0:
                continue
            elif current_state==1:
                new_decay=[]
            elif current_state==2:
                try:
                    new_decay.append([float(line[0]),float(line[1]),float(line[2]),float(line[3]),int(float(line[5]))])
                except:
                    new_decay=[]
                    current_state=0
            elif current_state==3:
                llp_decay.append(new_decay)
                
                
            pre_state = current_state
    
                
                
    return llp_decay  


def frame_transform(four_vector_in, boost_direction, beta=None, gamma=None):
    """
    Transform the four-vector (e, p_x, p_y, p_z,) with boost_direction and beta/gamma.
    
    INPUT:
    ---
    e, p_x, p_y, p_z: 
        four-vector to be transformed. Can be momentum or what ever 4-vector
    boost_direction:
        three-vector, the relative velocity of the frame to be transformed to
    beta/gamma:
        at least one needs to be given
        
    RETURN
    ---
    four_vector_out:
        the transformed four-vector
    """
    x, y, z = boost_direction
    e, p_x, p_y, p_z = four_vector_in
    four_vector_out = pylorentz.Momentum4(e, p_x, p_y, p_z).boost(x, y, z, beta=beta, gamma=gamma)
    return four_vector_out
                
    
 
def gen_llp(filename_llp4vec, filename_products, filename_output, mX, ctau, N_MC, rand_seed=1, nprint=1000,metadata=None):
    """
    Example:
    
    ```
    N_MC = 10000
    mX=15     # GeV
    ctau=1000 #m
    rand_seed = 1
    nprint=1000
    filename_llp4vec = f"{DATA_DIR}/MATHUSLA_LLPfiles_HXX/All_HXX_LLP4vectors/HXX_LLP4vectors_mX_{mX}_2perevent_unweighted.csv"
    filename_products = f"{DATA_DIR}/H_hadronic_decays_geant/bb_{mX}.txt"
    filename_output = f"scripts/LLP_test_bb_{mX}.txt"  
    
    ```
    """

    rng = np.random.default_rng(seed=rand_seed)
    random.seed(rand_seed)


    ## Decay volume geometry
    ## 100x100:
    # surface_height = 85.47
    # volume_floor_height = surface_height-19
    # volume_top_height = surface_height+6
    # full_module_width=99
    # module_to_CMS = 70
    ## 40x40:
    # surface_height = 85.47
    # volume_floor_height = surface_height+0.8
    # volume_top_height = surface_height+0.8+12.6
    # full_module_width=39
    # module_to_CMS = 70    
    ## 40x40, 6 layers:
    surface_height = 85.47
    volume_floor_height = surface_height+1
    volume_top_height = surface_height+1+11
    full_module_width=39
    module_to_CMS = 70        
    
    
    # The following xyz are all in CMS coordinate!
    cube_vertex_xyz = [full_module_width/2, volume_floor_height, module_to_CMS]
    cube_length_xyz = [full_module_width, volume_top_height-volume_floor_height, full_module_width]
    cube_ex = [-1, 0, 0 ]
    cube_ey = [0, 1, 0]
    cube_ez = [0, 0, 1 ]

    theta_range = np.array([np.arctan(volume_top_height/module_to_CMS),np.arctan(volume_floor_height/(module_to_CMS+full_module_width))])
    eta_range = -np.log(np.tan(theta_range/2))
    phi_range = np.array([-np.arctan(full_module_width*0.5/volume_floor_height),np.arctan(full_module_width*0.5/volume_floor_height)])

    # Load the LLP 4-vector
    print("Reading LLP 4-vector file")
    p4vec = np.loadtxt(filename_llp4vec,delimiter=",", dtype=float)
    # Deal with format difference:
    if len(p4vec[0]) ==5:
        # Exclude the weight column
        weight0 = p4vec[0][0]
        p4vec = p4vec[:,1:]
    else:
        weight0 = 1
    print(f"  {len(p4vec)} 4-vectors")
    print("Reading LLP decay product file")
    llp_decay = read_decay_file(filename_products)
    print(f"  {len(llp_decay)} LLP decay events")

    #-----------------------------
    # setup output file
    if not os.path.exists(os.path.dirname(filename_output)):
        Path(os.path.dirname(filename_output)).mkdir( parents=True, exist_ok=True )
    tree_writer = iroot.tfile_writer("data", filename_output)
    tree_writer.define_branch("Gen_pdgID", 'vector<double>')    
    tree_writer.define_branch("Gen_x", 'vector<float>')    
    tree_writer.define_branch("Gen_y", 'vector<float>')    
    tree_writer.define_branch("Gen_z", 'vector<float>')    
    tree_writer.define_branch("Gen_t", 'vector<float>')    
    tree_writer.define_branch("Gen_px", 'vector<float>')    
    tree_writer.define_branch("Gen_py", 'vector<float>')    
    tree_writer.define_branch("Gen_pz", 'vector<float>')    
    
    # weight 1,2
    N_LLP = len(p4vec)
    N_decay = len(llp_decay)
    weight1 = N_LLP/N_MC
    weight2 = (phi_range[1]-phi_range[0])/2/np.pi

    # Sample events
    # if N_MC>=N_decay:
    #     N_MC=N_decay
    #     print(" Requested event greater than the pool. Limiting to {N_MC} events.")
    # Get random index for decay sample
    i_decay_list=random.sample(range(N_decay), N_decay)
    

    
    # Calculate the maximum probability for the LLP to decay within the volume
    L_min = np.linalg.norm([volume_floor_height, module_to_CMS])
    L_max = np.linalg.norm([volume_top_height, full_module_width*0.5,  module_to_CMS+full_module_width])
    
    print("Calculating kinematics of 4vecs")
    vec_remaining = []
    vec_Boost = []
    vec_Gamma = []
    vec_Pt = []
    vec_Pl = []
    vec_Px = []
    vec_Py = []    
    vec_eta = []
    vec_phi = []
    vec_Pdecay=[]
    vec_L1=[]
    vec_L2=[]
    vec_point1=[]
    vec_point2=[]
    for ivec, vec in enumerate(p4vec[:]):
        P = np.linalg.norm(vec[1:4])
        P_l = vec[3]
        P_t = np.linalg.norm(vec[1:3])
        Boost = P/mX
        Gamma = np.sqrt(mX**2+P**2)/mX        
        eta = np.arctanh(P_l/P)
        
        # 1) eta cut
        if eta<eta_range[0] or eta>eta_range[1]:
            continue            
            
        # 2) phi sampling
        phi = rng.uniform(*phi_range)
        P_x = P_t*np.sin(phi)
        P_y = P_t*np.cos(phi) 
        
        # 3) calculate decay probability within the decay volume
        point1, point2 = line_cube_intersect([0,0,0], [P_x, P_y, P_l], cube_vertex_xyz, cube_length_xyz, cube_ex, cube_ey, cube_ez) # Find the intersect point of the track to the decay volume
        L1, L2 = np.linalg.norm(point1),np.linalg.norm(point2)
        if (L1==L2 and L1==0):
            continue        
        Pdecay = np.exp(-L1/(Boost*ctau)) - np.exp(-L2/(Boost*ctau))
        
        vec_remaining.append(vec)
        vec_Boost.append(Boost)
        vec_Gamma.append(Gamma)
        vec_Pt.append(P_t)
        vec_Pl.append(P_l)
        vec_Px.append(P_x)
        vec_Py.append(P_y)        
        vec_eta.append(eta)
        vec_phi.append(phi)
        vec_Pdecay.append(Pdecay)
        vec_L1.append(L1)
        vec_L2.append(L2)
        vec_point1.append(point1)
        vec_point2.append(point2)
    
    P_decay_max = np.max(vec_Pdecay)
     
    print(f"Selecting {N_MC} events out of {N_decay}")
    n_processed = 0
    vecs_used = 0
    vecs_used_list = []
    all_vertex= []
    while n_processed<N_MC:
        ivec = rng.integers(low=0,high=len(vec_remaining))
        vec = vec_remaining[ivec]
        vecs_used +=1
        
        # Retrieve the already calculated kinematics
        Boost=vec_Boost[ivec]
        Gamma=vec_Gamma[ivec]
        P_l=vec_Pl[ivec]
        P_t=vec_Pt[ivec] 
        P_x=vec_Px[ivec]
        P_y=vec_Py[ivec]         
        eta=vec_eta[ivec]
        phi=vec_phi[ivec]
        L1=vec_L1[ivec]
        L2=vec_L2[ivec]
        P_decay=vec_Pdecay[ivec]
        point1 = vec_point1[ivec]
        point2 = vec_point2[ivec]
        
        
        # 1) eta cut (no longer needed)
        # 2) phi sampling (no longer needed)
        # 3) exponential sampling
        L_sampled = util.rand_exp(Boost*ctau, L1, L2, rng)
        # 4) Drop with probablity of P_decay to deweight the event
        rand_drop = rng.random()
        if rand_drop>(P_decay/P_decay_max):
            continue
            
            
        # This 4-vec passed all cuts now
        point_decay = point1 + (point2-point1)*(L_sampled-L1)/(L2-L1)
        point_decay_G4 = [point_decay[2], point_decay[0], -point_decay[1]+surface_height]
        all_vertex.append(point_decay_G4)


        # 4) Get a LLP decay, and do Lorentz transform
        i_decay = i_decay_list[n_processed%N_decay]
        decays = llp_decay[i_decay]
        decays_transformed = []
        for decay in decays:
            decay_temp = frame_transform(decay[:4],-point_decay, beta=None, gamma=Gamma)
            decays_transformed.append([decay_temp[0], decay_temp[1],decay_temp[2],decay_temp[3], decay[4]]) # E, px, py, pz, PID

        weight_total = weight0*weight1*weight2*P_decay_max


        # # 5) Write event to a file
        # # Remember to transform into Geant coordinate. That's why px, py, pz are swapped, and turned from GeV to MeV, and [m]->[mm]
        # f_out.write(f"### weight {weight_total} {weight0} {weight1} {weight2} {P_decay_max}\n")
        # f_out.write(f"n {n_processed}  -1  0.0 0.0 0.0 {P_l*1000}  {P_x*1000}  {-P_y*1000} \n")
        # for decay in decays_transformed:
        #     f_out.write(f"\t {decay[4]}  {point_decay_G4[0]*1000} {point_decay_G4[1]*1000} {point_decay_G4[2]*1000} {decay[3]*1000}  {decay[1]*1000}  {-decay[2]*1000} \n")

        keys = ["Gen_pdgID","Gen_x","Gen_y","Gen_z","Gen_t","Gen_px","Gen_py","Gen_pz"]
        data = {key: [] for key in keys}
        for decay in decays_transformed:
            data["Gen_pdgID"].append(decay[4])
            data["Gen_x"].append(point_decay_G4[0]*1000 - (module_to_CMS+full_module_width*0.5)*1000)
            data["Gen_y"].append(point_decay_G4[1]*1000)
            data["Gen_z"].append(- point_decay_G4[2]*1000)
            data["Gen_t"].append(0)
            data["Gen_px"].append(decay[3]*1000)
            data["Gen_py"].append(decay[1]*1000)
            data["Gen_pz"].append(decay[2]*1000)
        tree_writer.fill(data)
        
        # update counter and print
        n_processed+=1
        vecs_used_list.append(ivec)
        if n_processed%nprint==0:
            print(f"  Processed {n_processed} events")
            # print("   - number of 4-vector sampled:", vecs_used)
    
    
    tree_writer.write_and_close()            
    
#     vecs_used_weights = vecs_used/N_LLP
#     f_out.write(f"# LLP 4-vec used for {vecs_used_weights}  times on average \n")
#     f_out.write(f"# final weight {vecs_used_weights}\n")
#     f_out.close()
#     # print(f"\n ** LLP 4-vec used for {vecs_used_weights} times on average ** \n")
#     print(f"  Total {n_processed} events processed")
#     print(f"Output saved as {filename_output} for simulation use")    
    
#     # Generate simulation script
#     filename_output_macro = Path(filename_output).with_suffix(".mac")
#     f_out2 = open(filename_output_macro,"w+")
#     f_out2.write(f"""/det/select Box 
# /gen/select file_reader 
# /gen/file_reader/pathname {filename_output}
# /run/beamOn {n_processed}""")
    
#     print(f"Macro file saved at  {filename_output_macro} for simulation use")    
    
    return vecs_used_list
    
    
    
    
    
    
        
    
    
    
# -------------------------------
# Read the filereader file
def read_raw_vertex_weight(filename):
    """
    Return:
    [[x, y, z, px, py, pz],[x, y, z, px, py, pz],...]
    """
    vertices = []
    weights = []
    with open(filename, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break
            if "#" in line:
                if "### weight" in line:
                    line = line.split()
                    weights.append([float(line[2]), float(line[3]), float(line[4])])
                else:
                    continue
            else:
                line = line.split()
                if len(line)>0 and line[0] == "n":
                    while True:
                        line2 = f.readline()
                        if "#" in line2:
                            continue                           
                        line2 = line2.split()
                        if len(line2)>0:
                            break
                    vertex = [float(line2[1]), float(line2[2]), float(line2[3]), float(line[6]), float(line[7]), float(line[8])] # x, y, z, px, py, pz
                    vertices.append(vertex)
                else:
                    continue
    return np.array(vertices), np.array(weights)




def cumsum0(y):
    y=np.array(y)
    return np.cumsum(y) - 0.5*y[0] - 0.5*y

def sample_neutron_energy_zenith(flux_dict, energies, num_samples=1, energy_range=None, box_size_meter = [10,10], seed = 1):
    """
    Samples neutron energy and zenith angle from given differential flux distributions.
    
    Parameters:
        flux_dict (dict):
            - Keys are zenith angles (in degrees).
            - Values are numpy arrays of differential flux corresponding to the energy points.
        energies (numpy array):
            - A numpy array of energy values corresponding to the flux arrays.
        num_samples (int):
            - The number of neutron samples to generate.
        energy_range (tuple, optional):
            - A tuple (E_min, E_max) to specify the energy range to sample from.
    
    Returns:
        events: numpy array of shape (num_samples, 2):
            - Each row contains (sampled xyz, sampled uvw, sampled_ekin [GeV])
        flux: total flux [/cm2/s] within the given energy range
    """
    # Setup random generator
    rng = np.random.default_rng(seed = seed)
    
    # Convert angles to a sorted list
    zenith_angles = np.array(sorted(flux_dict.keys()))
    zenith_radians = zenith_angles * np.pi / 180  # Convert to radians

    # Apply energy range filter if specified
    if energy_range is not None:
        E_min, E_max = energy_range
        energy_mask = (energies >= E_min) & (energies <= E_max)
        filtered_energies = energies[energy_mask]
        filtered_flux_dict = {angle: flux[energy_mask] for angle, flux in flux_dict.items()}
    else:
        filtered_energies = energies
        filtered_flux_dict = flux_dict
    
    # Compute total flux at each zenith angle (integrating over filtered energy)
    total_flux_per_angle = np.array([np.trapz(filtered_flux_dict[angle], filtered_energies) for angle in zenith_angles])
    
    # Compute the solid angle element for each zenith angle
    zenith_radians_diff = np.diff(zenith_radians)
    d_omega = 2 * np.pi * (np.cos(zenith_radians[1:] - 0.5*zenith_radians_diff) - np.cos(zenith_radians[1:] + 0.5*zenith_radians_diff))
    d_omega = np.concatenate(([2 * np.pi * (1-np.cos(0.5*zenith_radians_diff[0])) * 2], d_omega))
    
    # Normalize flux by solid angle to get correct probability distribution
    weighted_flux = total_flux_per_angle * d_omega
    weighted_cdf = cumsum0(weighted_flux)
    weighted_cdf /= weighted_cdf[-1]  # Normalize to [0,1]

    total_flux = sum(weighted_flux)/sum(d_omega)
    
    # Sample continuous zenith angles using inverse transform sampling
    sampled_zenith_angles = np.interp(rng.random(num_samples), weighted_cdf, zenith_radians)
    
    sampled_energies = []

    # print(zenith_radians)
    # plt.plot(zenith_radians, weighted_cdf)
    for angle in sampled_zenith_angles:
        # Interpolate flux for continuous angles
        nearest_indices = np.searchsorted(zenith_radians, angle, side='left')
        if nearest_indices == 0:
            interpolated_flux = filtered_flux_dict[zenith_angles[0]]
        elif nearest_indices == len(zenith_angles):
            interpolated_flux = filtered_flux_dict[zenith_angles[-1]]
        else:
            lower_angle, upper_angle = zenith_angles[nearest_indices - 1], zenith_angles[nearest_indices]
            lower_flux, upper_flux = filtered_flux_dict[lower_angle], filtered_flux_dict[upper_angle]
            interpolated_flux = lower_flux + (upper_flux - lower_flux) * (angle - lower_angle * np.pi / 180) / ((upper_angle - lower_angle) * np.pi / 180)
        
        # Normalize the flux to create a probability distribution for energy sampling
        cumulative_distribution = cumsum0(interpolated_flux)
        cumulative_distribution /= cumulative_distribution[-1]  # Normalize to [0,1]
        
        # Sample energy using inverse transform sampling
        random_value = rng.random()
        sampled_energy = np.interp(random_value, cumulative_distribution, filtered_energies)
        sampled_energies.append(sampled_energy)

    sampled_xs = (rng.random(num_samples)-0.5) * box_size_meter[0]
    sampled_ys = (rng.random(num_samples)-0.5) * box_size_meter[1]
    sampled_zs = np.zeros_like(sampled_xs)

    phis = rng.random(num_samples) * 2 * np.pi
    sampled_us = np.sin(sampled_zenith_angles) * np.cos(phis)
    sampled_vs = np.sin(sampled_zenith_angles) * np.sin(phis)
    sampled_ws = -np.cos(sampled_zenith_angles)
    
    
    # Return samples as a NumPy array
    return np.column_stack((sampled_xs, sampled_ys ,sampled_zs, sampled_us, sampled_vs, sampled_ws, sampled_energies)), total_flux

def gen_cosmic(filename_output, particle = "n" ,num_samples=1, energy_range=[1e5,1e12], rand_seed=1,  nprint=1000,metadata=None, box_size_meter = [10,10], box_z = 0):
    """
    
    Parameters:
        filename_output: str
        particle: str, one of {"n", "p"}

    Return:
        total_flux: float
            total flux in this energy range
    
    Example:
    """

    # setup output file
    if not os.path.exists(os.path.dirname(filename_output)):
        Path(os.path.dirname(filename_output)).mkdir( parents=True, exist_ok=True )
    tree_writer = iroot.tfile_writer("data", filename_output)
    tree_writer.define_branch("Gen_pdgID", 'vector<double>')    
    tree_writer.define_branch("Gen_x", 'vector<float>')    
    tree_writer.define_branch("Gen_y", 'vector<float>')    
    tree_writer.define_branch("Gen_z", 'vector<float>')    
    tree_writer.define_branch("Gen_t", 'vector<float>')    
    tree_writer.define_branch("Gen_px", 'vector<float>')    
    tree_writer.define_branch("Gen_py", 'vector<float>')    
    tree_writer.define_branch("Gen_pz", 'vector<float>')   

    # Load flux table
    if particle=="n":
        data = joblib.load(os.path.dirname(os.path.realpath(__file__)) + "/data/flux_table_neutron.joblib")
        egrid = data["egrid"]
        flux_table = data["flux_table"]
        pdgid = 2112
        mass = 0.939
    elif particle=="p":
        data = joblib.load(os.path.dirname(os.path.realpath(__file__)) + "/data/flux_table_proton.joblib")        
        egrid = data["egrid"]
        flux_table = data["flux_table"]
        pdgid = 2212
        mass = 0.938
    else:
        raise ValueError('Set particle to one of {"n", "p"} ')

    samples, total_flux = sample_neutron_energy_zenith(flux_table, egrid, num_samples=num_samples, energy_range=energy_range, box_size_meter = box_size_meter)

    print(f"Total flux in selected energy range : {total_flux} /cm2/s")

    for i in range(num_samples):
        keys = ["Gen_pdgID","Gen_x","Gen_y","Gen_z","Gen_t","Gen_px","Gen_py","Gen_pz"]
        data = {key: [] for key in keys}
        sample = samples[i]
        ek = sample[-1]
        p = np.sqrt(ek**2 + 2 * ek * mass)
        
        data["Gen_pdgID"].append(pdgid)
        data["Gen_x"].append(sample[0] * 1000) # m->mm
        data["Gen_y"].append(sample[1] * 1000) # m->mm
        data["Gen_z"].append(box_z * 1000) # m->mm
        data["Gen_t"].append(0)
        data["Gen_px"].append(p*sample[0+3] * 1000) # GeV-> MeV
        data["Gen_py"].append(p*sample[1+3] * 1000) # GeV-> MeV
        data["Gen_pz"].append(p*sample[2+3] * 1000) # GeV-> MeV
        tree_writer.fill(data)
    
    tree_writer.write_and_close()    

    return total_flux

