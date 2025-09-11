import sys, os
import importlib
from importlib import reload
# importlib.import_module(module_name)
sys.path.append("../../python")

import matplotlib.pyplot as plt
import numpy as np

from simhelper import root
from simhelper import generator
reload(generator)

DATA_DIR = "/project/6075887/MATHUSLA/simulation/llp"

N_MC = 10000
rand_seed = 1
nprint=2000

for ctau in [1,30,10,30,100,300,1000,3000,10000]:  # in mm 
    print(f"---------------ctau={ctau} m-----------------")
    for mX in [15,25,35,45,55]:
        print(f"---------------mX={mX} GeV-----------------")
        filename_llp4vec = f"{DATA_DIR}/0_decay_files_raw/MATHUSLA_LLPfiles_HXX/All_HXX_LLP4vectors/HXX_LLP4vectors_mX_{mX}_2perevent_unweighted.csv"
        filename_products = f"{DATA_DIR}/0_decay_files_raw/MATHUSLA_LLPfiles_HXX/H_hadronic_decays_geant/bb_{mX}.txt"
        filename_output = f"{DATA_DIR}/1_decay_files_transformed/MATHUSLA_LLPfiles_HXX/LLP_bb_{mX}_ctau_{ctau}.root"        
        vec_used = generator.gen_llp(filename_llp4vec, filename_products, filename_output, mX, ctau, N_MC, rand_seed=1, nprint=nprint,metadata=None)
