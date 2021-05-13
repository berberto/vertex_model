import itertools
# get_ipython().magic(u'matplotlib inline')
import numpy as np
import matplotlib
# %matplotlib qt5
import matplotlib.pyplot as plt
from matplotlib import animation, rc
# from IPython.display import HTML
import vertex_model as model
from vertex_model.run_select import run_simulation_INM, definecolors
import vertex_model.initialisation as init

import time
timestart = time.time()

type_=4
# L_point = [-.4, -0.2,-0.3, -0.05, 0.075, 0.15]
# G_point = [0.14, 0.12, 0.1, 0.065, 0.04, 0.02]
# i=4
G=0.075 # G_point[i]
L=0.04 # L_point[i]
K=1.0

rand = np.random.RandomState(1999) #random number to choose Lambda
params = [K,G,L]  # K=x[0],G=x[1],L=x[2]

no_simulations = 1
#Cell_cycle_lengths_progenitors = []
#Cell_cycle_lengths_pois = []
#final_history = []
for i in range(no_simulations):
    history = run_simulation_INM(params,1,rand,type_)  # timend = 1 for test
    # Rows 32-49 -> copy in ANALYSIS
    """
    last = history[-1] # get state of the cells object from the previous step -> this + following rows -> analysis
    gone = last.empty()
    divided = gone
    diffntd = gone
    if 'poisoned' in last.properties:
	    pois = last.properties['poisoned']==1
	    healty = last.properties['poisoned']==0
	    divided = gone & healty
	    diffntd = gone & pois
    progenitor_deadcells_last = np.where(divided)[0]
    differentiated_last = np.where(diffntd)[0]
    correct_cell_cycle_lengths = last.properties['age'][progenitor_deadcells_last]
    pois_cell_cycle_lengths = last.properties['age'][differentiated_last]
    
    #new_cell_cycle_lengths = last.properties['age'][deadcells_last] # cell cycle lengths = age of deadcells
    Cell_cycle_lengths_progenitors = np.append(Cell_cycle_lengths_progenitors, correct_cell_cycle_lengths)
    Cell_cycle_lengths_pois = np.append(Cell_cycle_lengths_pois, pois_cell_cycle_lengths)
    final_history = np.append(final_history, history) # get a final history array

    print(last.properties['nucl_pos'])
    print(min(last.properties['nucl_pos']))
    print(max(last.properties['nucl_pos']))
    print(last.properties['A0_final'])
    print(len(last.properties['final_nucl_pos']))
    """

timeend = time.time()
print(timeend - timestart)

import dill
with open ("history_model_1.pkl", "wb") as file:
    dill.dump(history, file)



