
# # All the function to run a simulation

#Import libraries

#########################################################################PARALELL################################################

#########################################################################PARALELL################################################
# %matplotlib tk 
#get_ipython().magic(u'matplotlib') #to use model.animate and see video alive
import itertools
import numpy as np
import matplotlib.pyplot as plt
import vertex_model as model
import vertex_model.initialisation as init
from vertex_model.forces import TargetArea, Tension, Perimeter, Pressure
import os
import seaborn as sns
import warnings
warnings.filterwarnings('ignore') #Don't show warnings
from .Gobal_Constant import dt, viscosity, t_G1, t_G2, t_S, A_c, J, pos_d, T1_eps, P, microns, time_hours, expansion_constant #file with necessary constants

diff_rate_hours=0.05 #differentiation rate (1/h) 



def check_estage():
    print("Running a %s hours"%J + " %s"%pos_d_v)
    print("dt=%s"%dt)            #time step
    print("viscosity=%s" %viscosity)  #viscosity*dv/dt = F
    print("A_c=%s"%A_c) #critical area
    print("T1_eps =%s"%T1_eps)
    

# run simulation - input is a generator = simulation => returns array
def run(simulation,N_step,skip):
    return [cells.copy() for cells in itertools.islice(simulation,0,int(N_step),int(skip))] # itertools.islice looks through the generator's elements. => list of cell objects = the history



def division_axis(mesh,face_id,rand):
    """Choose a random division axis (given as a pair of boundary edges to bisect) for the given cell.
    
    The first edge is chosen randomly from the bounding edges of the cell with probability proportional 
    to edge length. The second edge is then fixed to be n_edge/2 from the first. 
    """
    edges = mesh.boundary(face_id)
    if edges==[-1]:
        print('here')
        os._exit(1)
    p = np.cumsum(mesh.length[edges])
    e0 = p.searchsorted(rand.rand()*p[-1])
    return edges[e0],edges[e0-len(edges)//2]  

def bin_by_xpos(cells,percentiles):
    vx = cells.mesh.vertices[0]
    #simple 'midpoint' as mean of vertex positions
    # np.bincount = Count number of occurrences of each value in array of non-negative ints.
    # 1.0, 0.5; mid_x / width?
    mid_x = np.bincount(cells.mesh.face_id_by_edge,weights=vx)
    counts = np.maximum(np.bincount(cells.mesh.face_id_by_edge),1.0)
    mid_x = mid_x / counts 
    width = cells.mesh.geometry.width
    return np.searchsorted(percentiles,(mid_x/width + 0.5) % 1.0)   

#simulation without division - thermalisation without ages of cells and INM - evolves a mesh
def basic_simulation(cells,force,dt=dt,T1_eps=0.04): # a generator
    while True:
        cells.mesh , number_T1 = cells.mesh.transition(T1_eps)
        F = force(cells)/viscosity
        expansion = 0.05*np.average(F*cells.mesh.vertices,1)*dt
        dv = dt*model.sum_vertices(cells.mesh.edges,F) 
        cells.mesh = cells.mesh.moved(dv).scaled(1.0+ expansion)
        yield cells # outputs the cells object but continues to run (unlike return)!

# simulation with division and INM (no differentiation rate domain) # type 0
def simulation_with_division(cells,force,dt=dt,T1_eps=T1_eps,lifespan=100.0,rand=None): #(cells,force,dt=0.001,T1_eps=0.04,lifespan=100.0,rand=Non
    properties = cells.properties
    properties['parent'] = cells.mesh.face_ids #save the ids to control division parents-daugthers 
    properties['ageingrate'] = np.random.normal(1.0/lifespan,0.2/lifespan,len(cells)) #degradation rate per each cell
    properties['ids_division'] = [] #save ids of the cell os the division when it's ready per each time step
    properties['force_x'] = []
    properties['force_y'] = []
    properties['T1_angle_pD'] = []
    properties['Division_angle_pD'] = []
    expansion = np.array([0.0,0.0])
    while True:
        #cells id where is true the division conditions: living cells & area greater than 2 & age cell in mitosis 
        ready = np.where(~cells.empty() & (cells.mesh.area>=A_c) & (cells.properties['age']>=(t_G1+t_S+t_G2)))[0]  
        if len(ready): #these are the cells ready to undergo division at the current timestep
            properties['ageingrate'] =np.append(properties['ageingrate'], np.abs(np.random.normal(1.0/lifespan,0.2/lifespan,2*len(ready))))
            properties['age'] = np.append(properties['age'],np.zeros(2*len(ready)))
            properties['parent'] = np.append(properties['parent'],np.repeat(properties['parent'][ready],2))  # Daugthers and parent have the same ids
            properties['ids_division'] = ready
            edge_pairs = [division_axis(cells.mesh,cell_id,rand) for cell_id in ready] #New edges after division 
            cells.mesh = cells.mesh.add_edges(edge_pairs) #Add new edges in the mesh
            for i in range(len(ready)):
                commun_edges = np.intersect1d(cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-2*(i+1)))[0]],cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-1-2*i))[0]])
                division_new_edge=np.where(cells.mesh.length==np.max(commun_edges))[0]
                properties['Division_angle_pD']= np.append(properties['Division_angle_pD'],cells.mesh.edge_angle[division_new_edge][0])
        properties['age'] = properties['age']+dt*properties['ageingrate'] #add time step depending of the degradation rate 
        
        #Calculate z nuclei position (Apical-Basal movement), depending of the cell cycle phase time and age of the cell
        N_G1=1-1.0/t_G1*properties['age'] 
        N_S=0 
        N_G2=1.0/(t_G2)*(properties['age']-(t_G1+t_S))
        properties['zposn'] = np.minimum(1.0,np.maximum(N_G1,np.maximum(N_S,N_G2)))

        
        """Target area function depending age and z nuclei position"""
        properties['A0'] = (properties['age']+1.0)*0.5*(1.0+properties['zposn']**2)
        
        #############BAETTI VID
        
        #ath BAETTI D VID
        cells.mesh , number_T1, d= cells.mesh.transition(T1_eps)  #check edges verifing T1 transition
        if len(number_T1)>0:
            for ii in number_T1:
                index = cells.mesh.face_id_by_edge[ii]
                properties['T1_angle_pD'] =np.append(properties['T1_angle_pD'],cells.mesh.edge_angle[ii])
        F = force(cells)/viscosity  #force per each cell force= targetarea+Tension+perimeter+pressure_boundary 

        #EG BAETTI VID
        len_modified=np.matrix.copy(cells.mesh.length)
        #len_modified[np.where((properties['parent_group'][cells.mesh.face_id_by_edge]==1) & np.logical_not(properties['parent_group'][cells.mesh.face_id_by_edge[cells.mesh.edges.reverse]]==0))]*=0.12
        
        tens = (0.5*cells.by_edge('Lambda', 'Lambda_boundary')/len_modified)*cells.mesh.edge_vect
        tens= tens - tens.take(cells.mesh.edges.prev, 1)
        F+=tens

        dv = dt*model.sum_vertices(cells.mesh.edges,F) #movement of the vertices using eq: viscosity*dv/dt = F
        properties['force_x'] = F[0]*viscosity
        properties['force_y'] = F[1]*viscosity   
        
        if hasattr(cells.mesh.geometry,'width'):
            expansion[0] = expansion_constant*np.average(F[0]*cells.mesh.vertices[0])*dt/(cells.mesh.geometry.width**2)
        if hasattr(cells.mesh.geometry,'height'): #Cylinder mesh doesn't have 'height' argument
            expansion[1] = np.average(F[1]*cells.mesh.vertices[1])*dt/(cells.mesh.geometry.height**2)
        cells.mesh = cells.mesh.moved(dv).scaled(1.0+expansion)
        yield cells 

def pairwise_crowding_force (delta_z, a=1.0, s=0.1):
    '''
    Calculates the pairwise contributions for the crowding force as a 
    function of the difference in AB nuclear position between neighbouring
    cells.
    It accepts as input an array (and the function is applied element-wise)
    '''
    y = delta_z/s
    return a * np.exp(-y**2/2.)*y

def crowding_force (cells, a=1.0, s=0.1):
    '''
    For every cell calculate the total crowding force
    as the sum over its neighbours of the pairwise crowding force

    '''

    # Get the positions of all the nuclei.
    nucl_pos = cells.properties['nucl_pos']

    '''
    Get the ids of cells by their edges. The arrays obtained have
    a length equal to the number of edges, so many their entries
    are repeated as many times as the neighbours of each cell
    '''
    cell_ids = cells.mesh.face_id_by_edge   # ids of faces/cells
    neig_ids = cell_ids[cells.mesh.reverse] # ids of their neighbours

    '''
    checking if the cells alive correspond to the ids. This should **always**
    be the case when there is NO differentiation.
    '''
    # print(np.where(~cells.empty())[0] == np.sort(np.unique(cell_ids)))

    '''
    position of nuclei in cells and their neighbours.
    Same as above, many entries in "z" and "zn" are repeated, because we are
    going through cells via the edges
    '''
    z  = np.take(nucl_pos, cell_ids) # cells of interest
    zn = np.take(nucl_pos, neig_ids) # their neighbours

    '''
    calculate the difference in nucl_pos between neighbouring cells.
    Here there are NOT going to be repeated entries, in general (if not due to
    coincidences): "delta_z" is a property of the **pair** of neighbouring
    cells.
    '''
    delta_z = z - zn   # the difference in z for all pairs
    f_pairs = pairwise_crowding_force(delta_z, a=a, s=s)  # the corresponding pairwise forces

    # define an array as long as nucl_pos, where we calculate the force
    force = np.zeros_like(nucl_pos)

    '''
    Calculate the total crowding force by looping over **cells**, instead of edges.
    It may be more efficient, because the outer loop is over n_cells indices,
    instead of ~n_cells^2 (number of pairs, twice), and the "np.where" should be much
    faster than a simple elemeng-by-element search.
    '''
    for cell in np.unique(cell_ids):
        neigh_cells = np.where(neig_ids == cell)[0] # check the ids of the neighbours of "cell"
        force[cell] = np.sum(f_pairs[neigh_cells])  # sum the corresponding contributions

    # '''
    # Alternatively, loop over edges.
    # The following loop looks like it's over cells, but it actually is over
    # **pairs** of cells. This is because "cell_ids" is taken from "face_id_by_edge"
    # as discussed above.
    # '''
    # for i, cell in enumerate(cell_ids):
    #     '''
    #     "i" is the index of the **pair** of cells
    #     "cell" is the index of the cell of interest: it is used to
    #         access the corresponding element of the "force" array.
    #         The same index will appear as many times as its neighbours
    #     '''
    #     force[cell] += f_pairs[i]

    return force



        
 # simulation with division and INM (no differentiation rate domain) # type 0 Rebeca model 1 -> type 4
def simulation_with_division_model_1(cells,force,dt=dt,T1_eps=T1_eps,lifespan=100.0,rand=None, k=40, D=0.1): #(cells,force,dt=0.001,T1_eps=0.04,lifespan=100.0,rand=Non
    properties = cells.properties
    properties['parent'] = cells.mesh.face_ids #save the ids to control division parents-daugthers 
    properties['ageingrate'] = np.random.normal(1.0/lifespan,0.2/lifespan,len(cells)) #degradation rate per each cell
    properties['ids_division'] = [] #save ids of the cell os the division when it's ready per each time step
    properties['force_x'] = []
    properties['force_y'] = []
    properties['T1_angle_pD'] = []
    properties['Division_angle_pD'] = []
    properties['nucl_pos'] = properties['zposn'].copy()
    expansion = np.array([0.0,0.0])
    while True:
        #cells id where is true the division conditions: living cells & area greater than 2 & age cell in mitosis 
        ready = np.where(~cells.empty() & (cells.mesh.area>=A_c) & (cells.properties['age']>=(t_G1+t_S+t_G2)))[0]  
        if len(ready): #these are the cells ready to undergo division at the current timestep
            properties['ageingrate'] = np.append(properties['ageingrate'], np.abs(np.random.normal(1.0/lifespan,0.2/lifespan,2*len(ready))))
            properties['age'] = np.append(properties['age'],np.zeros(2*len(ready)))
            properties['parent'] = np.append(properties['parent'],np.repeat(properties['parent'][ready],2))  # Daugthers and parent have the same ids
            properties['nucl_pos'] = np.append(properties['nucl_pos'], np.ones(2*len(ready)))
            properties['ids_division'] = ready
            edge_pairs = [division_axis(cells.mesh,cell_id,rand) for cell_id in ready] #New edges after division 
            cells.mesh = cells.mesh.add_edges(edge_pairs) #Add new edges in the mesh
            for i in range(len(ready)):
                commun_edges = np.intersect1d(cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-2*(i+1)))[0]],cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-1-2*i))[0]])
                division_new_edge=np.where(cells.mesh.length==np.max(commun_edges))[0]
                properties['Division_angle_pD']= np.append(properties['Division_angle_pD'],cells.mesh.edge_angle[division_new_edge][0])
        properties['age'] = properties['age']+dt*properties['ageingrate'] #add time step depending of the degradation rate 
        
        #Calculate z nuclei position (Apical-Basal movement), depending of the cell cycle phase time and age of the cell
        N_G1=1-1.0/t_G1*properties['age'] 
        N_S=0 
        N_G2=1.0/(t_G2)*(properties['age']-(t_G1+t_S))
        properties['zposn'] = np.minimum(1.0,np.maximum(N_G1,np.maximum(N_S,N_G2)))

        properties['nucl_pos'] += k*(properties['zposn'] - properties['nucl_pos'])*dt + np.sqrt(2*D*dt)*np.random.randn(len(properties['zposn']))
        
        
        """Target area function depending age and z nuclei position"""
        properties['A0'] = (properties['age']+1.0)*0.5*(1.0+properties['nucl_pos']**2)
        
        #############BAETTI VID
        
        #ath BAETTI D VID
        cells.mesh , number_T1, d= cells.mesh.transition(T1_eps)  #check edges verifing T1 transition
        if len(number_T1)>0:
            for ii in number_T1:
                index = cells.mesh.face_id_by_edge[ii]
                properties['T1_angle_pD'] =np.append(properties['T1_angle_pD'],cells.mesh.edge_angle[ii])
        F = force(cells)/viscosity  #force per each cell force= targetarea+Tension+perimeter+pressure_boundary 

        #EG BAETTI VID
        len_modified=np.matrix.copy(cells.mesh.length)
        #len_modified[np.where((properties['parent_group'][cells.mesh.face_id_by_edge]==1) & np.logical_not(properties['parent_group'][cells.mesh.face_id_by_edge[cells.mesh.edges.reverse]]==0))]*=0.12
        
        tens = (0.5*cells.by_edge('Lambda', 'Lambda_boundary')/len_modified)*cells.mesh.edge_vect
        tens= tens - tens.take(cells.mesh.edges.prev, 1)
        F+=tens

        dv = dt*model.sum_vertices(cells.mesh.edges,F) #movement of the vertices using eq: viscosity*dv/dt = F
        properties['force_x'] = F[0]*viscosity
        properties['force_y'] = F[1]*viscosity   
        
        if hasattr(cells.mesh.geometry,'width'):
            expansion[0] = expansion_constant*np.average(F[0]*cells.mesh.vertices[0])*dt/(cells.mesh.geometry.width**2)
        if hasattr(cells.mesh.geometry,'height'): #Cylinder mesh doesn't have 'height' argument
            expansion[1] = np.average(F[1]*cells.mesh.vertices[1])*dt/(cells.mesh.geometry.height**2)
        cells.mesh = cells.mesh.moved(dv).scaled(1.0+expansion)
        yield cells 
       
        
        
        

# simulation with division, INM and no differentiation rate - type 0 with both pD and pMN
def simulation_with_division_clone(cells,force,dt=dt,T1_eps=T1_eps,lifespan=100.0,rand=None): #(cells,force,dt=0.001,T1_eps=0.04,lifespan=100.0,rand=Non
    properties = cells.properties
    properties['parent'] = cells.mesh.face_ids #save the ids to control division parents-daugthers 
    properties['ageingrate'] = np.random.normal(1.0/lifespan,0.2/lifespan,len(cells)) #degradation rate per each cell
    properties['ids_division'] = [] #save ids of the cell os the division when it's ready per each time step
    properties['force_x'] = []
    properties['force_y'] = []
    properties['T1_angle_pMN'] = []
    properties['T1_angle_pD'] = []
    properties['Division_angle_pMN'] = []
    properties['Division_angle_pD'] = []
    properties['ed']=[]
    properties['T1_edge']=[]
    properties['ids_removed']=[]

    expansion = np.array([0.0,0.0])
    while True:
        #cells id where is true the division conditions: living cells & area greater than 2 & age cell in mitosis 
        ready = np.where(~cells.empty() & (cells.mesh.area>=A_c) & (cells.properties['age']>=(t_G1+t_S+t_G2)))[0]  
        if len(ready): #these are the cells ready to undergo division at the current timestep
            properties['ageingrate'] =np.append(properties['ageingrate'], np.abs(np.random.normal(1.0/lifespan,0.2/lifespan,2.0*len(ready))))
            properties['age'] = np.append(properties['age'],np.zeros(2*len(ready)))
            properties['parent'] = np.append(properties['parent'],np.repeat(properties['parent'][ready],2))  # Daugthers and parent have the same ids
            properties['ids_division'] = ready
            properties['parent_group'] = np.append(properties['parent_group'],np.repeat(properties['parent_group'][ready],2)) #use to draw clones
            edge_pairs = [division_axis(cells.mesh,cell_id,rand) for cell_id in ready] #New edges after division
            properties['ed'].append(edge_pairs)
            cells.mesh = cells.mesh.add_edges(edge_pairs) #Add new edges in the mesh
            for i in range(len(ready)):
                commun_edges = np.intersect1d(cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-2*(i+1)))[0]],cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-1-2*i))[0]])
                division_new_edge=np.where(cells.mesh.length==np.max(commun_edges))[0]
                properties['la']=[commun_edges, division_new_edge]
                if properties['parent_group'][cells.mesh.n_face-2*(i+1)]==1:
                    properties['Division_angle_pMN']= np.append(properties['Division_angle_pMN'],cells.mesh.edge_angle[division_new_edge][0])
                else:
                    properties['Division_angle_pD']= np.append(properties['Division_angle_pD'],cells.mesh.edge_angle[division_new_edge][0])
        properties['age'] = properties['age']+dt*properties['ageingrate'] #add time step depending of the degradation rate 
        
        """Calculate z nuclei position (Apical-Basal movement), depending of the cell cycle phase time and age of the cell"""
        N_G1=1-1.0/t_G1*properties['age'] #nuclei position in G1 phase
        N_S=0
        N_G2=1.0/(t_G2)*(properties['age']-(t_G1+t_S))  #nuclei position in G2 and S phase
        properties['zposn'] = np.minimum(1.0,np.maximum(N_G1,np.maximum(N_S,N_G2)))
        
        
        """Target area function depending age and z nuclei position"""
        properties['A0'] = (properties['age']+1.0)*0.5*(1.0+properties['zposn']**2)
        

        cells.mesh , number_T1, ids_removed= cells.mesh.transition(T1_eps)
        if len(ids_removed)>0:
            properties['ids_removed'].append(ids_removed)

        if len(number_T1)>0:
            for ii in number_T1:
                properties['T1_edge']=np.append(properties['T1_edge'], ii)
                index = cells.mesh.face_id_by_edge[ii]
                if properties['parent_group'][index]==1:
                    properties['T1_angle_pMN'] =np.append(properties['T1_angle_pMN'],cells.mesh.edge_angle[ii])
                else:
                    properties['T1_angle_pD'] =np.append(properties['T1_angle_pD'],cells.mesh.edge_angle[ii])
        F = force(cells)/viscosity  #force per each cell force= targetarea+Tension+perimeter+pressure_boundary 
        dv = dt*model.sum_vertices(cells.mesh.edges,F) #movement of the vertices using eq: viscosity*dv/dt = F
        properties['force_x'] = F[0]*viscosity
        properties['force_y'] = F[1]*viscosity   
        if hasattr(cells.mesh.geometry,'width'):
            expansion[0] = expansion_constant*np.average(F[0]*cells.mesh.vertices[0])*dt/(cells.mesh.geometry.width**2)
        if hasattr(cells.mesh.geometry,'height'): #Cylinder mesh doesn't have 'height' argument
            expansion[1] = np.average(F[1]*cells.mesh.vertices[1])*dt/(cells.mesh.geometry.height**2)
        cells.mesh = cells.mesh.moved(dv).scaled(1.0+expansion)
        yield cells 


# simulation with division and all the cells are differentated as pNM diff_rate_hours (1/h) - type 1
def simulation_with_division_clone_differentiation(cells,force,dt=dt,T1_eps=T1_eps,lifespan=100.0,rand=None): #(cells,force,dt=0.001,T1_eps=0.04,lifespan=100.0,rand=Non
    properties = cells.properties
    properties['parent'] = cells.mesh.face_ids #save the ids to control division parents-daugthers 
    properties['ageingrate'] = np.random.normal(1.0/lifespan,0.2/lifespan,len(cells)) #degradation rate per each cell
    properties['ids_division'] = [] #save ids of the cell os the division when it's ready per each time step
    properties['poisoned'] = np.zeros(len(cells)) ### to add diferenciation rate in PMN

    expansion = np.array([0.0,0.0])
    while True:
        #cells id where is true the division conditions: living cells & area greater than 2 & age cell in mitosis 
        ready = np.where(~cells.empty() & (cells.mesh.area>=A_c) & (cells.properties['age']>=(t_G1+t_S+t_G2)))[0]  
        if len(ready): #these are the cells ready to undergo division at the current timestep
            properties['ageingrate'] =np.append(properties['ageingrate'], np.abs(np.random.normal(1.0/lifespan,0.2/lifespan,2.0*len(ready))))
            properties['age'] = np.append(properties['age'],np.zeros(2*len(ready)))
            properties['parent'] = np.append(properties['parent'],np.repeat(properties['parent'][ready],2))  # Daugthers and parent have the same ids
            properties['ids_division'] = ready
            properties['parent_group'] = np.append(properties['parent_group'],np.repeat(properties['parent_group'][ready],2)) #use to draw clones
            properties['poisoned'] = np.append(properties['poisoned'], np.zeros(2*len(ready))) ### to add diferenciation rate in PMN
            edge_pairs = [division_axis(cells.mesh,cell_id,rand) for cell_id in ready] #New edges after division 
            cells.mesh = cells.mesh.add_edges(edge_pairs) #Add new edges in the mesh
        ###### Defferenciation rate
        properties['poisoned'] = properties['poisoned'] - (properties['poisoned']-1) * (~(cells.empty()) & (rand.rand(len(cells)) < (dt*diff_rate_hours*time_hours)))
        properties['age'] = properties['age']+dt*properties['ageingrate'] #add time step depending of the degradation rate 
        
        """Calculate z nuclei position (Apical-Basal movement), depending of the cell cycle phase time and age of the cell"""
        N_G1=1-1.0/t_G1*properties['age'] #nuclei position in G1 phase
        N_S=0
        N_G2=1.0/(t_G2)*(properties['age']-(t_G1+t_S))  #nuclei position in G2 and S phase
        properties['zposn'] = np.minimum(1.0,np.maximum(N_G1,np.maximum(N_S,N_G2)))
        
        
        """Target area function depending age and z nuclei position"""
        properties['A0'] = (properties['age']+1.0)*0.5*(1.0+properties['zposn']**2)*(1.0-cells.properties['poisoned'])
        
        cells.mesh , number_T1= cells.mesh.transition(T1_eps)  #check edges verifing T1 transition
        F = force(cells)/viscosity  #force per each cell force= targetarea+Tension+perimeter+pressure_boundary 
        dv = dt*model.sum_vertices(cells.mesh.edges,F) #movement of the vertices using eq: viscosity*dv/dt = F
        
        if hasattr(cells.mesh.geometry,'width'):
            expansion[0] = expansion_constant*np.average(F[0]*cells.mesh.vertices[0])*dt/(cells.mesh.geometry.width**2)
        if hasattr(cells.mesh.geometry,'height'): #Cylinder mesh doesn't have 'height' argument
            expansion[1] = np.average(F[1]*cells.mesh.vertices[1])*dt/(cells.mesh.geometry.height**2)
        cells.mesh = cells.mesh.moved(dv).scaled(1.0+expansion)
        yield cells 


# simulation with division with INM and 2 diferent populations (with and without differentiation rate) - type 2
def simulation_with_division_clone_differenciation_3stripes(cells,force,dt=dt,T1_eps=T1_eps,lifespan=100.0,rand=None): #(cells,force,dt=0.001,T1_eps=0.04,lifespan=100.0,rand=Non

    
    properties = cells.properties
    properties['parent'] = cells.mesh.face_ids #save the ids to control division parents-daugthers 
    properties['ageingrate'] = np.random.normal(1.0/lifespan,0.2/lifespan,len(cells)) #degradation rate per each cell
    #floor plate cells
    no=len(properties['ageingrate'][np.where(properties['parent_group']==3)])
    properties['ageingrate'][np.where(properties['parent_group']==3)]=np.random.normal(0.5/lifespan,0.1/lifespan,no)

    properties['ids_division'] = [] #save ids of the cell os the division when it's ready per each time step
    properties['ids_division_1'] = [] #save ids of the cell os the division when it's ready per each time step and region
    properties['ids_division_02'] = [] #save ids of the cell os the division when it's ready per each time step and region
    properties['poisoned'] = np.zeros(len(cells)) ### to add diferenciation rate in PMN
    # properties['differentiation_rate']= np.zeros(len(cells),dtype=int)
    properties['force_x'] = []
    properties['force_y'] = []
    properties['T1_angle_pMN'] = []
    properties['T1_angle_pD'] = []
    properties['Division_angle_pMN'] = []
    properties['Division_angle_pD'] = []
    properties['deleted_edges']=[]
    properties['edges_division']=[]
    properties['T1_edge']=[]
    properties['nucl_pos'] = []
    expansion = np.array([0.0,0.0])
    
    # added by Rebeca:
    D = 0.15
    k = 100 # a constant
    a=0.2 # a constant controlling the size of force
    s=0.2
    
    
    while True:
        #cells id where is true the division conditions: living cells & area greater than 2 & age cell in mitosis 
        
        ready = np.where((~cells.empty()  &(cells.mesh.area>=A_c) & (cells.properties['age']>=(t_G1+t_S+t_G2)) & (cells.properties['zposn']>=0.75)) | (~cells.empty() & (cells.properties['parent_group']==3) &(cells.mesh.area>=0.4) & (cells.properties['age']>=(t_G1+t_S+t_G2)) & (cells.properties['zposn']>=0.75)))[0] 

        if len(ready): #these are the cells ready to undergo division at the current timestep
            #properties['ageingrate'] =np.append(properties['ageingrate'], np.repeat(properties['ageingrate'][ready],2))
            properties['age'] = np.append(properties['age'],np.zeros(2*len(ready)))
            properties['parent'] = np.append(properties['parent'],np.repeat(properties['parent'][ready],2))  # Daugthers and parent have the same ids
            properties['ids_division'] = ready
            properties['parent_group'] = np.append(properties['parent_group'],np.repeat(properties['parent_group'][ready],2)) #use to draw clones
            
            properties['poisoned'] = np.append(properties['poisoned'], np.zeros(2*len(ready))) ### to add diferenciation rate in PMN
            edge_pairs = [division_axis(cells.mesh,cell_id,rand) for cell_id in ready]  #New edges after division 
            properties['edges_division'].append(edge_pairs)
            cells.mesh = cells.mesh.add_edges(edge_pairs) #Add new edges in the mesh  
            

           
            
            for i in range(len(ready)):
                commun_edges = np.intersect1d(cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-2*(i+1)))[0]],cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-1-2*i))[0]])
                division_new_edge=np.where(cells.mesh.length==np.max(commun_edges))[0]
                
                
                if properties['parent_group'][ready[i]]==3:
                    properties['ageingrate'] = np.append(properties['ageingrate'], np.abs(np.random.normal(1.0/lifespan,0.2/lifespan,2)))
                else:
                    properties['ageingrate'] = np.append(properties['ageingrate'], np.abs(np.random.normal(1.0/lifespan,0.2/lifespan,2)))
                   # print len(properties['ageingrate'])
                #if properties['parent_group'][cells.mesh.n_face-2*(i+1)]==1:
                 #   properties['Division_angle_pMN']= np.append(properties['Division_angle_pMN'],cells.mesh.edge_angle[division_new_edge][0])
                #else:
                 #   properties['Division_angle_pD']= np.append(properties['Division_angle_pD'],cells.mesh.edge_angle[division_new_edge][0])
            #for ids in ready:
             #   if properties['parent_group'][ids]==1:
              #      properties['ids_division_1'] = np.append(properties['ids_division_1'], ids)
               # else:
                #    properties['ids_division_02'] = np.append(properties['ids_division_02'], ids)
        ###### Differentiation rate
        properties['differentiation_rate'] = time_hours*dt*(np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))[properties['parent_group']] #Used 0.02, 0.0002 & 1/13
        properties['poisoned'] = properties['poisoned'] - (properties['poisoned']-1) * (~(cells.empty()) & (rand.rand(len(cells)) < properties['differentiation_rate']))
        
        # add age ingrate only for alive cells
        alivecells = np.where(~cells.empty())[0]
        properties['age'][alivecells] += dt*properties['ageingrate'][alivecells] #add time step depending of the degradation rate 
        
        # define arrays for each phase:
        N_G1_new = np.ones_like(properties['age']) # vector with ones for as many cells (ages) as they are at the current timestep
        N_S_new = np.zeros_like(properties['age'])
        N_G2_new = np.zeros_like(properties['age'])
        N_M_new = np.ones_like(properties['age'])
        
        N_G1=1-1.0/t_G1*properties['age'] #nuclei position in G1 phase
        N_G1_new = N_G1_new + k*(N_G1 - N_G1_new)*dt + np.sqrt(2*D*dt)*np.random.randn(len(properties['age'])) #nuclei position in G1 phase
        N_S=0
        N_S_new = N_S_new + k*(N_S - N_S_new)*dt + np.sqrt(2*D*dt)*np.random.randn(len(properties['age']))
        
        N_G2=1.0/(t_G2)*(properties['age']-(t_G1+t_S))  #nuclei position in G2 and S phase
        N_G2_new = N_G2_new + k*(N_G2 - N_G2_new)*dt + np.sqrt(2*D*dt)*np.random.randn(len(properties['age']))
        
        N_M = 1
        N_M_new = N_M_new + k*(N_M - N_M_new)*dt + np.sqrt(2*D*dt)*np.random.randn(len(properties['age']))
        
        properties['zposn'] = np.minimum(1.0,np.maximum(N_G1,np.maximum(N_S,N_G2)))
        properties['nucl_pos'] = np.minimum(1.0,np.maximum(N_G1_new,np.maximum(N_S_new,N_G2_new))) # new position
        
        
          
        properties['zposn'][np.where(properties['parent_group']==3)]=0
       

        
        """Target area function depending age and z nuclei position"""
        properties['A0'] = (properties['age']+1.0)*0.5*(1.0+properties['zposn']**2)*(1.0-cells.properties['poisoned'])

        
        cells.mesh , number_T1, del_edges= cells.mesh.transition(T1_eps)
        
        properties['deleted_edges'].append(del_edges)
        #if len(number_T1)>0:
         #   for ii in number_T1:
          #      properties['T1_edge']=np.append(properties['T1_edge'], ii)
               # index = cells.mesh.face_id_by_edge[ii]
                #if properties['parent_group'][index]==1:
                 #   properties['T1_angle_pMN'] =np.append(properties['T1_angle_pMN'],cells.mesh.edge_angle[ii])
                #else:
                 #   properties['T1_angle_pD'] =np.append(properties['T1_angle_pD'],cells.mesh.edge_angle[ii])
        F = force(cells)/viscosity  #force per each cell force= targetarea+Tension+perimeter+pressure_boundary 

        #ADDING THE FLOOR PLATE
        len_modified=np.matrix.copy(cells.mesh.length)
        #this next line can be use for modifying all edges at the boundary
        #len_modified[np.where( ((properties['parent_group'][cells.mesh.face_id_by_edge]==3) & np.logical_not(properties['parent_group'][cells.mesh.face_id_by_edge[cells.mesh.edges.reverse]]==3)) | ((properties['parent_group'][cells.mesh.face_id_by_edge]==2) & (properties['parent_group'][cells.mesh.face_id_by_edge[cells.mesh.edges.reverse]]==3)) | ((properties['parent_group'][cells.mesh.face_id_by_edge]==4) & (properties['parent_group'][cells.mesh.face_id_by_edge[cells.mesh.edges.reverse]]==3)) )]*=0.004
        
        #modifying tension w.r.t. area
        for n in np.where( ((properties['parent_group'][cells.mesh.face_id_by_edge]==3) & np.logical_not(properties['parent_group'][cells.mesh.face_id_by_edge[cells.mesh.edges.reverse]]==3)))[0]:
            if abs(cells.mesh.area[cells.mesh.face_id_by_edge[n]]-cells.mesh.area[cells.mesh.face_id_by_edge[cells.mesh.edges.reverse[n]]])>0.4:
                len_modified[n]*=0.002
   
        tens = (0.5*cells.by_edge('Lambda', 'Lambda_boundary')/len_modified)*cells.mesh.edge_vect
        tens= tens - tens.take(cells.mesh.edges.prev, 1)
        F+=tens

        dv = dt*model.sum_vertices(cells.mesh.edges,F) #movement of the vertices using eq: viscosity*dv/dt = F
        #properties['force_x'] = F[0]*viscosity
        #properties['force_y'] = F[1]*viscosity
        if hasattr(cells.mesh.geometry,'width'):
            expansion[0] = expansion_constant*np.average(F[0]*cells.mesh.vertices[0])*dt/(cells.mesh.geometry.width**2)
        if hasattr(cells.mesh.geometry,'height'): #Cylinder mesh doesn't have 'height' argument
            expansion[1] = np.average(F[1]*cells.mesh.vertices[1])*dt/(cells.mesh.geometry.height**2)
        cells.mesh = cells.mesh.moved(dv).scaled(1.0+expansion)
        
        yield cells 

# simulation with division with INM and all the cells high without differentiation rate - type 3 ??
def simulation_with_division_clone_whole_tissue_differenciation(cells,force,dt=dt,T1_eps=T1_eps,lifespan=100.0,rand=None): #(cells,force,dt=0.001,T1_eps=0.04,lifespan=100.0,rand=None
    properties = cells.properties
    properties['parent'] = cells.mesh.face_ids #save the ids to control division parents-daugthers 
    properties['ageingrate'] = np.random.normal(1.0/lifespan,0.2/lifespan,len(cells)) #degradation rate per each cell
    properties['ids_division'] = [] #save ids of the cell os the division when it's ready per each time step
    properties['ids_division_1'] = [] #save ids of the cell os the division when it's ready per each time step and region
    properties['ids_division_02'] = [] #save ids of the cell os the division when it's ready per each time step and region
    properties['poisoned'] = np.zeros(len(cells)) ### to add diferenciation rate in PMN
    # properties['differentiation_rate']= np.zeros(len(cells),dtype=int)
    properties['force_x'] = []
    properties['force_y'] = []
    properties['T1_angle_pMN'] = []
    properties['T1_angle_pD'] = []
    properties['Division_angle_pMN'] = []
    properties['Division_angle_pD'] = []
    properties['ids_removed']=[]
    expansion = np.array([0.0,0.0])
    while True:
        #cells id where is true the division conditions: living cells & area greater than 2 & age cell in mitosis 
        ready = np.where(~cells.empty() & (cells.mesh.area>=A_c) & (cells.properties['age']>=(t_G1+t_S+t_G2)))[0]  
        if len(ready): #these are the cells ready to undergo division at the current timestep
            properties['ageingrate'] =np.append(properties['ageingrate'], np.abs(np.random.normal(1.0/lifespan,0.2/lifespan,2.0*len(ready)))) # use only for alive cells at current timestep
            properties['age'] = np.append(properties['age'],np.zeros(2*len(ready)))
            properties['parent'] = np.append(properties['parent'],np.repeat(properties['parent'][ready],2))  # Daugthers and parent have the same ids
            properties['ids_division'] = ready
            properties['parent_group'] = np.append(properties['parent_group'],np.repeat(properties['parent_group'][ready],2)) #use to draw clones
            properties['poisoned'] = np.append(properties['poisoned'], np.zeros(2*len(ready))) ### to add diferenciation rate in PMN
            # properties['differentiation_rate'] = np.append(properties['differentiation_rate'], diff_rate_hours*np.ones(2*len(ready)))          
            edge_pairs = [division_axis(cells.mesh,cell_id,rand) for cell_id in ready] #New edges after division 
            cells.mesh = cells.mesh.add_edges(edge_pairs) #Add new edges in the mesh
            for i in range(len(ready)):
                commun_edges = np.intersect1d(cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-2*(i+1)))[0]],cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-1-2*i))[0]])
                division_new_edge=np.where(cells.mesh.length==np.max(commun_edges))[0]
                if properties['parent_group'][cells.mesh.n_face-2*(i+1)]==1:
                    properties['Division_angle_pMN']= np.append(properties['Division_angle_pMN'],cells.mesh.edge_angle[division_new_edge][0])
                else: #antes estaba mal y ponia pD!!!! las simulaciones guardadas estan mal
                    properties['Division_angle_pMN']= np.append(properties['Division_angle_pMN'],cells.mesh.edge_angle[division_new_edge][0])
            for ids in ready:
                if properties['parent_group'][ids]==1:
                    properties['ids_division_1'] = np.append(properties['ids_division_1'], ids)
                else:#antes estaba mal y ponia 0_2!!!! las simulaciones guardadas estan mal
                    properties['ids_division_1'] = np.append(properties['ids_division_1'], ids)
        ###### Defferentiation rate
        properties['differentiation_rate'] = time_hours*dt*(np.array([0.0,diff_rate_hours,0.0,0.0,0.0,diff_rate_hours, 0.0]))[properties['parent_group']] #Used 0.02, 0.0002 & 1/13
        properties['poisoned'] = properties['poisoned'] - (properties['poisoned']-1) * (~(cells.empty()) & (rand.rand(len(cells)) < properties['differentiation_rate']))
        
        # add age ingrate only for alive cells
        alivecells = np.where(~cells.empty())[0]
        properties['age'][alivecells] += dt*properties['ageingrate'][alivecells] #add time step depending of the degradation rate 
        
        N_G1=1-1.0/t_G1*properties['age'] #nuclei position in G1 phase
        N_S=0
        N_G2=1.0/(t_G2)*(properties['age']-(t_G1+t_S))  #nuclei position in G2 and S phase
        properties['zposn'] = np.minimum(1.0,np.maximum(N_G1,np.maximum(N_S,N_G2)))
        
        
        """Target area function depending age and z nuclei position"""
        properties['A0'] = (properties['age']+1.0)*0.5*(1.0+properties['zposn']**2)*(1.0-cells.properties['poisoned'])
        
        cells.mesh , number_T1, ids_removed= cells.mesh.transition(T1_eps)

        properties['ids_removed'].append(ids_removed) #####BAETTI VID
        if len(number_T1)>0:
            for ii in number_T1:
                index = cells.mesh.face_id_by_edge[ii]
                if properties['parent_group'][index]==1:
                    properties['T1_angle_pMN'] =np.append(properties['T1_angle_pMN'],cells.mesh.edge_angle[ii])
                else:#antes estaba mal y ponia pD!!!! las simulaciones guardadas estan mal
                    properties['T1_angle_pMN'] =np.append(properties['T1_angle_pMN'],cells.mesh.edge_angle[ii])
        F = force(cells)/viscosity  #force per each cell force= targetarea+Tension+perimeter+pressure_boundary 
        dv = dt*model.sum_vertices(cells.mesh.edges,F) #movement of the vertices using eq: viscosity*dv/dt = F
        properties['force_x'] = F[0]*viscosity
        properties['force_y'] = F[1]*viscosity
    
        
        if hasattr(cells.mesh.geometry,'width'):
            expansion[0] = expansion_constant*np.average(F[0]*cells.mesh.vertices[0])*dt/(cells.mesh.geometry.width**2)
        if hasattr(cells.mesh.geometry,'height'): #Cylinder mesh doesn't have 'height' argument
            expansion[1] = np.average(F[1]*cells.mesh.vertices[1])*dt/(cells.mesh.geometry.height**2)
        cells.mesh = cells.mesh.moved(dv).scaled(1.0+expansion)
        yield cells 
        
    
# simulation with division with INM and 2 diferent populations (with and without differentiation rate) - type 4 (type 2 modified by Rebeca)
def simulation_with_division_clone_differenciation_3stripes_model_1(cells,force,dt=dt,T1_eps=T1_eps,lifespan=100.0,rand=None): #(cells,force,dt=0.001,T1_eps=0.04,lifespan=100.0,rand=Non

    
    properties = cells.properties
    properties['parent'] = cells.mesh.face_ids #save the ids to control division parents-daugthers 
    properties['ageingrate'] = np.random.normal(1.0/lifespan,0.2/lifespan,len(cells)) #degradation rate per each cell
    #floor plate cells
    no=len(properties['ageingrate'][np.where(properties['parent_group']==3)])
    properties['ageingrate'][np.where(properties['parent_group']==3)]=np.random.normal(0.5/lifespan,0.1/lifespan,no)

    properties['ids_division'] = [] #save ids of the cell os the division when it's ready per each time step
    properties['ids_division_1'] = [] #save ids of the cell os the division when it's ready per each time step and region
    properties['ids_division_02'] = [] #save ids of the cell os the division when it's ready per each time step and region
    properties['poisoned'] = np.zeros(len(cells)) ### to add diferenciation rate in PMN
    # properties['differentiation_rate']= np.zeros(len(cells),dtype=int)
    properties['force_x'] = []
    properties['force_y'] = []
    properties['T1_angle_pMN'] = []
    properties['T1_angle_pD'] = []
    properties['Division_angle_pMN'] = []
    properties['Division_angle_pD'] = []
    properties['deleted_edges']=[]
    properties['edges_division']=[]
    properties['T1_edge']=[]
    properties['nucl_pos']=[]
    expansion = np.array([0.0,0.0])
    
    # added by Rebeca:
    D = 0.15
    k = 100 # a constant
    a=0.2 # a constant controlling the size of force
    s=0.2
    

    
    while True:
        #cells id where is true the division conditions: living cells & area greater than 2 & age cell in mitosis 
        
        # added by Rebeca: properties['zposn'] >= 0.75
        ready = np.where((~cells.empty()  &(cells.mesh.area>=A_c) & (cells.properties['age']>=(t_G1+t_S+t_G2)) & (cells.properties['nucl_pos']>=0.75)) | (~cells.empty() & (cells.properties['parent_group']==3) &(cells.mesh.area>=0.4) & (cells.properties['age']>=(t_G1+t_S+t_G2)) & (cells.properties['nucl_pos']>=0.75)))[0]  

        if len(ready): #these are the cells ready to undergo division at the current timestep
            #properties['ageingrate'] =np.append(properties['ageingrate'], np.repeat(properties['ageingrate'][ready],2))
            properties['age'] = np.append(properties['age'],np.zeros(2*len(ready)))
            properties['parent'] = np.append(properties['parent'],np.repeat(properties['parent'][ready],2))  # Daugthers and parent have the same ids
            properties['ids_division'] = ready
            properties['parent_group'] = np.append(properties['parent_group'],np.repeat(properties['parent_group'][ready],2)) #use to draw clones
            
            properties['poisoned'] = np.append(properties['poisoned'], np.zeros(2*len(ready))) ### to add diferenciation rate in PMN
            edge_pairs = [division_axis(cells.mesh,cell_id,rand) for cell_id in ready]  #New edges after division 
            properties['edges_division'].append(edge_pairs)
            cells.mesh = cells.mesh.add_edges(edge_pairs) #Add new edges in the mesh
            
       
            
            for i in range(len(ready)):
                commun_edges = np.intersect1d(cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-2*(i+1)))[0]],cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-1-2*i))[0]])
                division_new_edge=np.where(cells.mesh.length==np.max(commun_edges))[0]
                
                

              
                if properties['parent_group'][ready[i]]==3:
                    properties['ageingrate'] = np.append(properties['ageingrate'], np.abs(np.random.normal(1.0/lifespan,0.2/lifespan,2)))
                else:
                    properties['ageingrate'] = np.append(properties['ageingrate'], np.abs(np.random.normal(1.0/lifespan,0.2/lifespan,2)))
                   # print len(properties['ageingrate'])
                #if properties['parent_group'][cells.mesh.n_face-2*(i+1)]==1:
                 #   properties['Division_angle_pMN']= np.append(properties['Division_angle_pMN'],cells.mesh.edge_angle[division_new_edge][0])
                #else:
                 #   properties['Division_angle_pD']= np.append(properties['Division_angle_pD'],cells.mesh.edge_angle[division_new_edge][0])
            #for ids in ready:
             #   if properties['parent_group'][ids]==1:
              #      properties['ids_division_1'] = np.append(properties['ids_division_1'], ids)
               # else:
                #    properties['ids_division_02'] = np.append(properties['ids_division_02'], ids)
        ###### Differentiation rate
        properties['differentiation_rate'] = time_hours*dt*(np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))[properties['parent_group']] #Used 0.02, 0.0002 & 1/13
        properties['poisoned'] = properties['poisoned'] - (properties['poisoned']-1) * (~(cells.empty()) & (rand.rand(len(cells)) < properties['differentiation_rate']))
        
        # add age ingrate only for alive cells -> not necessarily needed as all cells are alive in type 2
        alivecells = np.where(~cells.empty())[0]
        properties['age'][alivecells] += dt*properties['ageingrate'][alivecells] #add time step depending of the degradation rate 
        
        
        N_G1_new = np.ones_like(properties['age']) # vector with ones for as many cells (ages) as they are at the current timestep
        N_S_new = np.zeros_like(properties['age'])
        N_G2_new = np.zeros_like(properties['age'])
        N_M_new = np.ones_like(properties['age'])
        
        N_G1=1-1.0/t_G1*properties['age']
        N_G1_new = N_G1_new + k*(N_G1 - N_G1_new)*dt + np.sqrt(2*D*dt)*np.random.randn(len(properties['age'])) #nuclei position in G1 phase
       
        N_S=0 
        N_S_new = N_S_new + k*(N_S - N_S_new)*dt + np.sqrt(2*D*dt)*np.random.randn(len(properties['age']))
        
        N_G2=1.0/(t_G2)*(properties['age']-(t_G1+t_S))  #nuclei position in G2 and S phase
        N_G2_new = N_G2_new + k*(N_G2 - N_G2_new)*dt + np.sqrt(2*D*dt)*np.random.randn(len(properties['age']))
        
        N_M = 1
        N_M_new = N_M_new + k*(N_M - N_M_new)*dt + np.sqrt(2*D*dt)*np.random.randn(len(properties['age']))
        properties['zposn'] = np.minimum(1.0,np.maximum(N_G1,np.maximum(N_S,N_G2)))
        properties['nucl_pos'] = np.minimum(1.0,np.maximum(N_G1_new,np.maximum(N_S_new,N_G2_new)))
        
        
        properties['zposn'][np.where(properties['parent_group']==3)]=0

        
        """Target area function depending age and z nuclei position"""
        properties['A0'] = (properties['age']+1.0)*0.5*(1.0+properties['zposn']**2)*(1.0-cells.properties['poisoned'])

        
        cells.mesh , number_T1, del_edges= cells.mesh.transition(T1_eps)
        
        properties['deleted_edges'].append(del_edges)
        #if len(number_T1)>0:
         #   for ii in number_T1:
          #      properties['T1_edge']=np.append(properties['T1_edge'], ii)
               # index = cells.mesh.face_id_by_edge[ii]
                #if properties['parent_group'][index]==1:
                 #   properties['T1_angle_pMN'] =np.append(properties['T1_angle_pMN'],cells.mesh.edge_angle[ii])
                #else:
                 #   properties['T1_angle_pD'] =np.append(properties['T1_angle_pD'],cells.mesh.edge_angle[ii])
        F = force(cells)/viscosity  #force per each cell force= targetarea+Tension+perimeter+pressure_boundary 

        #ADDING THE FLOOR PLATE
        len_modified=np.matrix.copy(cells.mesh.length)
        #this next line can be use for modifying all edges at the boundary
        #len_modified[np.where( ((properties['parent_group'][cells.mesh.face_id_by_edge]==3) & np.logical_not(properties['parent_group'][cells.mesh.face_id_by_edge[cells.mesh.edges.reverse]]==3)) | ((properties['parent_group'][cells.mesh.face_id_by_edge]==2) & (properties['parent_group'][cells.mesh.face_id_by_edge[cells.mesh.edges.reverse]]==3)) | ((properties['parent_group'][cells.mesh.face_id_by_edge]==4) & (properties['parent_group'][cells.mesh.face_id_by_edge[cells.mesh.edges.reverse]]==3)) )]*=0.004
        
        #modifying tension w.r.t. area
        for n in np.where( ((properties['parent_group'][cells.mesh.face_id_by_edge]==3) & np.logical_not(properties['parent_group'][cells.mesh.face_id_by_edge[cells.mesh.edges.reverse]]==3)))[0]:
            if abs(cells.mesh.area[cells.mesh.face_id_by_edge[n]]-cells.mesh.area[cells.mesh.face_id_by_edge[cells.mesh.edges.reverse[n]]])>0.4:
                len_modified[n]*=0.002
   
        tens = (0.5*cells.by_edge('Lambda', 'Lambda_boundary')/len_modified)*cells.mesh.edge_vect
        tens= tens - tens.take(cells.mesh.edges.prev, 1)
        F+=tens

        dv = dt*model.sum_vertices(cells.mesh.edges,F) #movement of the vertices using eq: viscosity*dv/dt = F
        #properties['force_x'] = F[0]*viscosity
        #properties['force_y'] = F[1]*viscosity
        if hasattr(cells.mesh.geometry,'width'):
            expansion[0] = expansion_constant*np.average(F[0]*cells.mesh.vertices[0])*dt/(cells.mesh.geometry.width**2)
        if hasattr(cells.mesh.geometry,'height'): #Cylinder mesh doesn't have 'height' argument
            expansion[1] = np.average(F[1]*cells.mesh.vertices[1])*dt/(cells.mesh.geometry.height**2)
        cells.mesh = cells.mesh.moved(dv).scaled(1.0+expansion)
        
        yield cells         
        


#_clones
def definecolors(cells):
    peach = '#eed5b7'
    light_blue ='#87cefa'
    pink = '#ffc0cb'
    light_green = '#98fb98'
    import matplotlib.colors as colors
    vv=sns.color_palette("hls", 10)
    v=[colors.rgb2hex(colorrgb) for colorrgb in vv]
    palette = np.array([light_green, pink,light_green,'g','r','g','m','c','',peach])
    palette = np.array([v[1],v[0],v[1], v[1],v[4],v[5],v[6],v[7],v[8],v[9],peach])
    # palette = np.array([peach, 'b', 'r','g','r','b', peach,'w',light_blue,pink,light_green,pink,light_blue,'w'])
    #colors = np.array([1 if x==71 else 2 if x==80 else 3 if x==49 else 0 for x in cells.properties['parent']])
    #colors = cells.properties['parent_group']*np.array([1 if x in sampleset else 0 for x in cells.properties['parent']])
    colors = cells.properties['parent_group']
    return palette[colors]

"""Run simulation and save data functions"""
"""
def run_simulation(x):
    K=x[0]
    G=x[1]
    L=x[2]
    rand = np.random #np.random.RandomState(123456) #I have modified the random function because RamdomState takes always the same numbers
    mesh = init.cylindrical_hex_mesh(10,10,noise=0.2,rand=rand)
    cells = model.Cells(mesh,properties={'K':K,'Gamma':G,'P':0.0,'boundary_P':P,'Lambda':L, 'Lambda_boundary':0.5})
    cells.properties['age'] = np.random.rand(len(cells))
    force = TargetArea() + Tension() + Perimeter() + Pressure()
    history = run(simulation_with_division(cells,force,rand=rand),500.0/dt,1.0/dt)
    # model.animate_video_mpg(history)
    return history
"""
def run_simulation_INM(x, timend,rand, sim_type):
    global dt
    #sim_type 0 simulation_with_division_clone (no differentiation rate)
    #sim_type 1 simulation_with_division_clone_differentiation (all differentiation rate)
    #sim_type 2 simulation_with_division_clone_differenciation_3stripes (2 population with and without diffentiation rate)
    #sim_type 3 simulation_with_division_clone_whole_tissue_differenciation (differentiation rate everywhere)
    # print(dt)
    K=x[0]
    G=x[1]
    L=x[2]
    rand1 = np.random.RandomState(123456) #I have modified the random function because RamdomState takes always the same numbers
    #mesh = init.cylindrical_hex_mesh(2,2,noise=0.2,rand=rand1)
    mesh = init.toroidal_hex_mesh(20,20,noise=0.2,rand=rand1)
    cells = model.Cells(mesh,properties={'K':K,'Gamma':G,'P':0.0,'boundary_P':P,'Lambda':L, 'Lambda_boundary':0.5})
    cells.properties['age'] = np.random.rand(len(cells))

    force = TargetArea()  + Perimeter() + Pressure()
    
    

####eg er ad baeta thessu vid
    cells.properties['parent_group'] = np.zeros(len(cells),dtype=int) #use to draw clone
    cells.properties['parent_group'] = bin_by_xpos(cells,np.cumsum([0.475,0.5,0.475]))

    print("Start thermalization...")
    history1 = run(simulation_with_division(cells,force,rand=rand1),10.0/dt,0.5/dt)
    cells = history1[-1].copy() #last timestep in the 'thermalization' phase -> use this as the initial condition for the actual run below
    # sampleset = np.random.choice(cells.mesh.face_ids,20, replace=False) #take 7 different ramdon cells to follow clones

 
    cells.properties['parent_group'] = np.zeros(len(cells),dtype=int) #use to draw clone
    cells.properties['parent'] = cells.mesh.face_ids #save the ids to control division parents-daugthers 
    cells.properties['parent_group'] = bin_by_xpos(cells,np.cumsum([0.225, 0.15, 0.075, 0.1, 0.075, 0.15, 0.225]))

    print(f"Start simulation of type {sim_type}...")
    # sampleset = np.random.choice(cells.mesh.face_ids[cells.properties['parent_group']==1],10, replace=False) #take 10 different ramdon cells to follow clones
    if sim_type == 0:
        history = run(simulation_with_division_clone(cells,force,rand=rand),(timend)/dt,1.0/dt)
        history[-1].properties['parent_group'] = np.zeros(len(history[-1].properties['parent_group']),dtype=int)
    if sim_type == 1:
        history = run(simulation_with_division_clone_differentiation(cells,force,rand=rand),(timend)/dt,1.0/dt)
        history[-1].properties['parent_group'] = np.zeros(len(history[-1].properties['parent_group']),dtype=int)+1
    if sim_type == 2:
        #we take ventral and dorsal time per phase cell cycle if we are in the 2 pop part, because pNM are ventral and pD are dorsal
        history = run(simulation_with_division_clone_differenciation_3stripes(cells,force,rand=rand),(timend)/dt,1.0/dt)
        cells.properties['parent_group'] = cells.properties['parent_group'] #+np.array([3 if x in sampleset else 0 for x in cells.properties['parent']])
    if sim_type == 3:
        history = run(simulation_with_division_clone_whole_tissue_differenciation(cells,force,rand=rand),(timend)/dt,1.0/dt)
        history[-1].properties['parent_group'] = np.zeros(len(history[-1].properties['parent_group']),dtype=int)
    if sim_type == 4:
        history = run(simulation_with_division_model_1(cells,force,rand=rand),(timend)/dt,1.0/dt)
   
 
        
    return history
   




def save_data(I,history,outputdir):

    # """Save create the folder to save the data"""
    # outputdir="output_Dorsal%s"%J+"_A_c%s_viscosity003"%int(10*A_c)
    # if not os.path.exists(outputdir): # if the folder doesn't exist create it
    #     os.makedirs(outputdir)
    """Information of Area, Perimeter, Neigbours: time_mean and end. And force in time and final age distribution"""
    force = TargetArea()  + Perimeter() + Pressure()
    generation_a=[cells.mesh.area for cells in history] # collects a list of arrays each of which contains area of each cell at any given time 
    generation_p=[cells.mesh.perimeter for cells in history]
    generation_n=[np.bincount(cells.mesh.face_id_by_edge) for cells in history]
    generation_f=[force(cells) for cells in history] 
    death=[cells.empty() for cells in history]
    #mean variables in time
    area_mean=[]
    perimeter_mean=[]
    neigh_mean=[]
    force_mean=[]
    force_units=[]
    area_total=[]
    number_cells = []
    ids_face_by_edge = []
    for i in range(0,len(history)-1):
        valid=np.where(~death[i] & (generation_a[i]>0))[0] # -> valid checks which cells haven't died
        number_cells = np.append(number_cells, len(valid))
        area_mean=np.append(area_mean,np.mean(generation_a[i][valid]))
        area_total=np.append(area_total,np.sum(generation_a[i][valid]))
        perimeter_mean=np.append(perimeter_mean,np.mean(generation_p[i][valid]))
        neigh_mean=np.append(neigh_mean,np.mean(generation_n[i][valid]))
        # force_mean=np.append(force_mean,np.mean(np.sqrt(generation_f[i][0][valid]**2+generation_f[i][1][valid]**2)))
        # force_units=np.append(force_units, np.sum(np.sqrt(generation_f[i][0][valid]**2+generation_f[i][1][valid]**2))/9.0)


    #end value of Area, Perimeter and Neighbour of the tissue
    """Define files per different values of area division condition"""
    outputdirname=outputdir+'/'
    outfile_a=outputdirname+"area_mean_%0.3f"%I
    outfile_a_total=outputdirname+"area_total_%0.3f"%I
    outfile_p=outputdirname+"perimeter_mean_%0.3f"%I
    outfile_n=outputdirname+"neigh_mean_%0.3f"%I
    # outfile_f=outputdirname+"force_mean_%0.3f"%I
    # outfile_f_units=outputdirname+"force_units_%0.3f"%I

    with open(outfile_a,"w") as tfile:
        np.savetxt(tfile,area_mean)
    with open(outfile_a_total,"w") as tfile:
        np.savetxt(tfile,area_total)
    with open(outfile_p,"w") as tfile:
        np.savetxt(tfile,perimeter_mean)
    with open(outfile_n,"w") as tfile:
        np.savetxt(tfile,neigh_mean)
    # with open(outfile_f,"w") as tfile:
    #     np.savetxt(tfile,force_mean)
    # with open(outfile_f_units,"w") as tfile:
    #     np.savetxt(tfile,force_units)
        
    vert_x = cells.mesh.vertices[0]
    vert_y = cells.mesh.vertices[1]
    ids_face_by_edge = cells.mesh.face_id_by_edge  
    outfile_a_end=outputdirname+"area_end_%0.3f"%I
    outfile_p_end=outputdirname+"perimeter_end_%0.3f"%I
    outfile_n_end=outputdirname+"neigh_end_%0.3f"%I
    outfile_age_end=outputdirname+"age_end_%0.3f"%I
    outfile_nc=outputdirname+"number_cells_%0.3f"%I
    outfile_x_vert = outputdirname + "vertices_x_end_%0.3f"%I
    outfile_y_vert = outputdirname + "vertices_y_end_%0.3f"%I
    outfile_ids_face_edge = outputdirname + "ids_face_edge_%0.3f"%I
    
    valid=np.where(~death[-1] & (generation_a[-1]>0))[0]    
    with open(outfile_a_end,"w") as tfile:
        np.savetxt(tfile,np.array(generation_a[-1][valid]))
    with open(outfile_p_end,"w") as tfile:
        np.savetxt(tfile,np.array(generation_p[-1][valid]))
    with open(outfile_n_end,"w") as tfile:
        np.savetxt(tfile,np.array(generation_n[-1][valid]))
    with open(outfile_age_end,"w") as tfile:
        np.savetxt(tfile,np.array(cells.properties['age'][valid]))
    with open(outfile_nc,"w") as tfile:
        np.savetxt(tfile,number_cells)
    with open(outfile_x_vert,"w") as tfile:
        np.savetxt(tfile, vert_x)
    with open(outfile_y_vert,"w") as tfile:
        np.savetxt(tfile, vert_y)
    with open(outfile_ids_face_edge,"w") as tfile:
        np.savetxt(tfile, ids_face_by_edge)
