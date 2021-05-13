import itertools
import numpy as np
import matplotlib.pyplot as plt
import vertex_model as model
import vertex_model.initialisation as init
from vertex_model.forces import TargetArea, Tension, Perimeter, Pressure

import dill
with open("history_model_1.pkl", "rb") as file:
    history = dill.load(file)

last = history[-1]

# Cell cycle length
"""
deadcells = np.where(last.empty())[0]
cell_cycle_lengths = last.properties['age'][deadcells]
plt.hist(cell_cycle_lengths, bins = 50, density=True)
plt.title('Cell cycle length distribution Model 1')
plt.show()
"""


"""
plt.hist(last.properties['nucl_pos'], bins = 50)
plt.xlim(0,1)
plt.title("Nucleus position Model 1")
plt.show()
"""

alivecells = np.where(~last.empty())[0]


plt.hist(last.mesh.area[alivecells], bins = 50)
plt.title('Last Area Distribution Model 1')
plt.show()



print(len(alivecells))
print(len(last.mesh.area))
print(len(last.mesh.area[alivecells]))


