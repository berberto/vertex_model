
# coding: utf-8

# In[2]:

import itertools
import numpy as np
from matplotlib.collections import PolyCollection
import matplotlib.pyplot as plt
import os
from .permutations import cycles
#from matplotlib import animation


def _draw_edges(mesh, ax):
    w = mesh.vertices - mesh.vertices.take(mesh.edges.rotate, 1)  # winding
    to_draw = mesh.edges.ids[(mesh.edges.ids < mesh.edges.reverse) | (np.abs(w[0])+np.abs(w[1]) > 0.1)]
    start, end = mesh.vertices.take(to_draw, 1), mesh.vertices.take(mesh.edges.next[to_draw], 1)

    n = np.empty(len(start[0]))
    n.fill(np.nan)
    x = np.dstack([start[0], end[0], n]).ravel()
    y = np.dstack([start[1], end[1], n]).ravel()

    ax.plot(x, y, 'k-', linewidth=1.0)

def _draw_edges_non(mesh, ax):
    w = mesh.vertices - mesh.vertices.take(mesh.edges.rotate, 1)  # winding
    to_draw = mesh.edges.ids[(mesh.edges.ids < mesh.edges.reverse) | (np.abs(w[0])+np.abs(w[1]) > 0.1)]
    start, end = mesh.vertices.take(to_draw, 1), mesh.vertices.take(mesh.edges.next[to_draw], 1)

    n = np.empty(len(start[0]))
    n.fill(np.nan)
    x = np.dstack([start[0], end[0], n]).ravel()
    y = np.dstack([start[1], end[1], n]).ravel()

def _draw_midpoints(cells, ax):
    s = cells.vertices/(np.maximum(np.bincount(cells.edges.cell), 1)[cells.edges. cell])
    sx, sy = np.bincount(cells.edges.cell, weights=s[0]), np.bincount(cells.edges. cell, weights=s[1])
    mx, my = sx[cells.edges.cell], sy[cells.edges.cell]

    ax.scatter(mx, my, s=3.0, c='b', edgecolors='none', marker='o')


_PALETTES = {  # .<=3      4       5       6       7       8     >=9    edges
    'Default': '#000000 #33cc33 #ffff19 #ccccb2 #005cb8 #cc2900 #4a0093 #b0fc3e',
    'CB':      '#edf8fb #edf8fb #bfd3e6 #9ebcda #8c96c6 #8856a7 #810f7c #000000',
}
_PALETTES = {name: np.array([clr.split()[0]]*4+clr.split()[1:])
             for name, clr in _PALETTES.items()}


def _draw_faces(mesh, ax, facecolors, edgecolor='k'):
    order, labels = cycles(mesh.edges.next)
    counts = np.bincount(labels)

    vs = mesh.vertices.T.take(order, 0)

    faces, face_ids = [], []
    cell_ids = mesh.face_id_by_edge[order]
    boundary = mesh.boundary_faces if mesh.has_boundary() else []

    for (i, c) in zip(counts, np.cumsum(counts)):
        cell_id = cell_ids[c-i]
        if cell_id in boundary:
            continue
        faces.append(vs[c-i:c])
        face_ids.append(cell_ids[c-i])

    coll = PolyCollection(faces, facecolors=facecolors[face_ids], edgecolors=edgecolor,
                          linewidths=2.0)
    ax.add_collection(coll)

def _draw_faces_no_edge(mesh, ax, facecolors):
    order, labels = cycles(mesh.edges.next)
    counts = np.bincount(labels)

    vs = mesh.vertices.T.take(order, 0)

    faces, face_ids = [], []
    cell_ids = mesh.face_id_by_edge[order]
    boundary = mesh.boundary_faces if mesh.has_boundary() else []

    for (i, c) in zip(counts, np.cumsum(counts)):
        cell_id = cell_ids[c-i]
        if cell_id in boundary:
            continue
        faces.append(vs[c-i:c])
        face_ids.append(cell_ids[c-i])

    coll = PolyCollection(faces, facecolors=facecolors[face_ids])
    ax.add_collection(coll)

def _draw_geometry(geometry, ax=None):
    # Torus
    if hasattr(geometry, 'width') and hasattr(geometry, 'height'):
        w, h = geometry.width, geometry.height
        # ax.add_patch(plt.Rectangle((-0.5*w, -0.5*h), w, h, fill=False, linewidth=2.0))


def draw(cells, ax=None, size=None):
    if not ax:
        fig = plt.figure()
        ax = fig.gca()  
    ax.cla()
    facecolors = cells.properties.get('color', None)

    mesh = cells.mesh.recentre()

    if facecolors is None:
        _draw_edges(mesh, ax)
    else:
        _draw_faces(mesh, ax, facecolors)
# _draw_midpoints(cells,ax)
    _draw_geometry(mesh.geometry, ax)

    ax.set_xticks([])
    ax.set_yticks([])

    size = size or 2.0*np.max(mesh.vertices[0])
    lim = [-0.55*size, 0.55*size]
    ax.set_xlim(lim)
    ax.set_ylim(lim)

    plt.draw()

def animate(cells_array, facecolours='Default'):
    v_max = np.max((np.max(cells_array[0].mesh.vertices), np.max(cells_array[-1].mesh.vertices)))
    size = 2.0*v_max
    #to don't draw 
    fig = plt.figure()#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAa
    ax = fig.gca()#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAa
    fig.set_size_inches(6, 6)#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAa
    for cells in cells_array:#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAa
        draw(cells, ax, size)#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAa
  #hasta aqui to don't draw        

def draw1(cells, ax=None, size=None):
    if not ax:
        fig = plt.figure()
        ax = fig.gca()  
    ax.cla()
    facecolors = cells.properties.get('color', None)

    mesh = cells.mesh.recentre()

    if facecolors is None:
        _draw_edges_non(mesh, ax)
    else:
        _draw_faces_no_edge(mesh, ax, facecolors)
    _draw_midpoints(cells,ax)
    _draw_geometry(mesh.geometry, ax)

    ax.set_xticks([])
    ax.set_yticks([])

    size = size or 2.0*np.max(mesh.vertices[0])
    lim = [-0.55*size, 0.55*size]
    ax.set_xlim(lim)
    ax.set_ylim(lim)

    # plt.draw()

#Doesn't work
# def animate_video1(cells_array): 
#     fig = plt.figure(); 
#     ax = fig.gca(); 
#     fig.set_size_inches(6,6); 
#     i=0
#     # frames=[]
#     # for cells in cells_array:
#     #     draw(cells,ax,size)
#     #     i=i+1
#     #     frame="images/image%03i.png" % i
#     #     fig.savefig(frame,dpi=300)
#     #     frames.append(frame)

#     # call the animator.  blit=True means only re-draw the parts that have changed.
#     anim = animation.FuncAnimation(fig, lambda i: cells_array[i], #init_func=init,
#                                frames=200, interval=20, blit=True)
#     anim.save('animation.mp4',fps=30)
        
def animate_video_mpg(cells_array,name_file,facecolours='Default'):    
    v_max = np.max((np.max(cells_array[0].mesh.vertices), np.max(cells_array[-1].mesh.vertices)))
    size = 2.0*v_max
    outputdir="images"
    if not os.path.exists(outputdir): # if the folder doesn't exist create it
        os.makedirs(outputdir)
    fig = plt.figure(); 
    ax = fig.gca(); 
    fig.set_size_inches(6,6); 
    i=0
    frames=[]
    for cells in cells_array:
        draw(cells,ax,size)
        i=i+1
        frame="images/image%03i.png" % i
        fig.savefig(frame,dpi=500)
        frames.append(frame)  
    os.system("mencoder 'mf://images/image*.png' -mf type=png:fps=20 -ovc lavc -lavcopts vcodec=wmv2 -oac copy  -o " + name_file)  
    # os.system("ffmpeg -framerate 5/1 -i images/image%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p test.mp4") #for Mac computer
    for frame in frames: os.remove(frame)

def animate_video_mpg_zoom(cells_array,name_file,facecolours='Default'):    
    v_max = np.max((np.max(cells_array[0].mesh.vertices[1]), np.max(cells_array[-1].mesh.vertices[1])))
    size = 2.0*v_max
    outputdir="images"
    if not os.path.exists(outputdir): # if the folder doesn't exist create it
        os.makedirs(outputdir)
    fig = plt.figure(); 
    ax = fig.gca(); 
    fig.set_size_inches(6,6); 
    i=0
    frames=[]
    for cells in cells_array:
        draw(cells,ax,size)
        i=i+1
        frame="images/zoom_image%03i.png" % i
        fig.savefig(frame,dpi=500)
        frames.append(frame)  
    os.system("mencoder 'mf://images/zoom_image*.png' -mf type=png:fps=20 -ovc lavc -lavcopts vcodec=wmv2 -oac copy  -o " + name_file)  
    # os.system("ffmpeg -framerate 5/1 -i images/image%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p test.mp4") #for Mac computer
    for frame in frames: os.remove(frame)

# In[2]:

#os.system("mencoder images/image*.png -mf type=png:fps=20 -ovc lavc -lavcopts vcodec=wmv2 -oac copy  -o  clones.mpg")

