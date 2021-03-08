
###### Global constants of the system

global dt, viscosity, t_G1, t_G2, A_c, t_S, T1_eps, P, microns, time_hours, expansion_constant, pos_d
dt=0.001            #time step
viscosity= 0.02  #viscosity*dv/dt = F
T1_eps = 0.01
A_c=1.3 #critical area

"""60 hours Dorsal"""
pos_d="Dorsal"
J=60
t_G1=0.4            #Time proportion in G1 phase
t_S = 0.5*(2.0/3.0)            #Time proportion in S phase
t_G2 = 0.5*(1.0/3.0)           #Time proportion in G2 phase 
t_M=1.0-(t_G1+t_S+t_G2)
P=1


experimental_perimeter = 9.93 #average, Dorsal e10.5
simulation_perimeter = 2.17  #average, Dorsal 60 hours
microns = experimental_perimeter/simulation_perimeter 
experimental_cell_cycle = 13.0 #average, hours 60 hours
simulation_cell_cycle = 100.0 #average step units 
time_hours = experimental_cell_cycle/simulation_cell_cycle 

expansion_constant = 50
