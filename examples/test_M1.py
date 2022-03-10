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
import sys

import time
import dill


type_=4

# parameters of the vertex model
G=0.075
L=0.04
K=1.0

# parameters of the nucleus A-B stochastic dynamics
# sys.argv to pass from command line
# call as, e.g.: python test_M1.py 100 0.15
try:
	k=float(sys.argv[1]) # 100
	D=float(sys.argv[2]) # 0.15
	r=float(sys.argv[3])
except:
	print("Missing command line parameters. Execute as\n")
	print("   python test_M1.py  <k>  <D>  <r = a/k>\n")
	exit()

# parameters of the crowding force
s=0.2
a=r*k

# adding the nuclear movement parameters to 'params'
params = [K,G,L,k,D,s,a]

def run_sim (run_id):
    timestart = time.time()
    history = run_simulation_INM(params,1,rand,type_)  # timend = 1 for test
    timeend = time.time()

    with open (f"model_1_k{k}_D{D}_run_{run_id}.pkl", "wb") as file:
        dill.dump(history, file)

    return timeend - timestart

np.random.seed(1999)
rand = np.random.RandomState(1999)
no_simulations = 5

for run in range(1,no_simulations+1):
	run_sim(run)




