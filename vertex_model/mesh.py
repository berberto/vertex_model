# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import functools
import copy
import os
import numpy as np

from .permutations import cycle


# ====Cells data structure====
def cached_property(func):
    """A decorator for caching class properties."""
    @functools.wraps(func)
    def getx(self):
        try:
            return self._cache[func]
        except AttributeError:
            self._cache = {}
        except KeyError:
            pass
        self._cache[func] = res = func(self)
        return res
    return property(getx)


_MAX_CELLS = 5*10**4


class Edges(object):
    """
    Attributes:
        ids: The array of integers [0,1,...,n_edge-1] giving the ID of each (directed) edge.
        rotate: An array of n_edge integers giving the ID of the next edge around a vertex.
        rotate2: An array of n_edge integers giving the ID of the previous edge around a vertex.
        reverse: An array of n_edge integers giving the ID of the reversed edge.
    """
    # share fixed arrays between objects for efficiency
    IDS = np.arange(6*_MAX_CELLS)
    ROTATE = np.roll(IDS.reshape(-1, 3), 2, axis=1).ravel()
    ROTATE2 = ROTATE[ROTATE]

    def __init__(self, reverse):
        self.reverse = reverse
        n = len(reverse)
        self.rotate = self.ROTATE[:n]
        self.rotate2 = self.ROTATE2[:n]
        self.ids = self.IDS[:n]

        # make 'immutable' so that various computed properties can be cached
        reverse.setflags(write=False)

    def __len__(self):
        return len(self.reverse)

    @cached_property
    def next(self):
        """An array of n_edge integers giving the ID of the next edge around the (left) face."""
        return self.rotate[self.reverse]

    @cached_property
    def prev(self):
        """An array of n_edge integers giving the ID of the previous edge around the (left) face."""
        return self.reverse[self.rotate2]


class MeshGeometry(object):
    def scaled(self, scale):
        return self

    def recentre(self, mesh):
        return mesh


class Plane(MeshGeometry):
    pass


class Torus(MeshGeometry):
    def __init__(self, width, height):
        self.width = width
        self.height = height

    def scaled(self, scale):
        return Torus(self.width*scale[0], self.height*scale[1])

    def recentre(self, mesh):
        counts = mesh.edge_counts()[mesh.face_id_by_edge]
        for i, L in enumerate((self.width, self.height)):
            mid = np.bincount(mesh.face_id_by_edge, mesh.vertices[i])[mesh.face_id_by_edge]
            mesh.vertices[i, mid > 0.5*L*counts] -= L
            mesh.vertices[i, mid < -0.5*L*counts] += L
        return mesh


class Cylinder(MeshGeometry):
    def __init__(self, width):
        self.width = width

    def scaled(self, scale):
        return Cylinder(self.width*scale[0])

    def recentre(self, mesh):
        counts = mesh.edge_counts()[mesh.face_id_by_edge]
        L = self.width
        mid = np.bincount(mesh.face_id_by_edge, mesh.vertices[0])[mesh.face_id_by_edge]
        mesh.vertices[0, mid > 0.5*L*counts] -= L
        mesh.vertices[0, mid < -0.5*L*counts] += L
        return mesh


class Mesh(object):
    """
    An mesh in the 'half-edge' representation.

    meshes are intended to be immutable so that we can cache various computed properties (such as
    lengths of edges) for efficiency. Methods to move vertices, or change the topology all return new
    Meshs object, leaving the original intact.

    Args:
        edges: An Edges object holding topological information about the network of cells.
        vertices: A (2,n_edge) array of floats giving the x,y position of the start vertex of each edge.
        face_ids: An (n_edge) array of integers holding the face_id for each edge
        geometry: An instance of MeshGeometry, holding information about the background geometry
        n_face: total number of cells including 'dead' cells which are no longer in edges network.
    """

    def __init__(self, edges, vertices, face_id_by_edge, geometry, n_face=None, boundary_faces=None):
        self.edges = edges
        self.vertices = vertices
        self.face_id_by_edge = face_id_by_edge
        face_id_by_edge.setflags(write=False)
        self.n_face = n_face or np.max(face_id_by_edge) + 1
        self.geometry = geometry
        self.boundary_faces = boundary_faces
        
    def copy(self):
        mesh = copy.copy(self)
        mesh._cache = {}
        return mesh

    def scaled(self, scale):
        mesh = self.copy()
        mesh.vertices = mesh.vertices * scale[:, None]
        mesh.geometry = mesh.geometry.scaled(scale)
        return mesh

    def recentre(self):
        mesh = self.copy()
        return mesh.geometry.recentre(mesh)

    def has_boundary(self):
        return self.boundary_faces is not None

    @cached_property
    def boundary_edges(self):
        edges = np.hstack([self.boundary(face) for face in self.boundary_faces])
        return edges

    @cached_property
    def face_ids(self):
        return np.arange(self.n_face)

    @cached_property
    def edge_vect(self):
        v = self.vertices
        dv = v.take(self.edges.next, 1) - v
        if self.has_boundary():
            if np.any(self.boundary_edges==-1):
                # print self.boundary_edges
                os._exit(1)
            dv[:, self.boundary_edges] = -dv.take(self.edges.reverse[self.boundary_edges], 1)
        return dv

    @cached_property
    def length(self):
        vect = self.edge_vect
        return np.sqrt(vect[0]**2 + vect[1]**2)

    @cached_property
    def edge_angle(self):
        vect = self.edge_vect
        angle = np.arctan(vect[1]/vect[0])
        return angle

    @cached_property
    def area(self):
        return np.bincount(self.face_id_by_edge, 0.5*np.cross(self.vertices.T, self.edge_vect.T), self.n_face)

    @cached_property
    def perimeter(self):
        return np.bincount(self.face_id_by_edge, self.length, self.n_face)

    @cached_property
    def d_area(self):
        dv = self.edge_vect
        dA = np.empty_like(dv)
        dA[0], dA[1] = dv[1], -dv[0]
        return 0.5*(dA + dA.take(self.edges.prev, 1))

    @cached_property
    def d_perimeter(self):
        dL = self.edge_vect/self.length
        return dL - dL.take(self.edges.prev, 1)

    def edge_counts(self):
        return np.bincount(self.face_id_by_edge)

    @cached_property
    def _edge_lookup(self):
        lookup = np.empty(self.n_face, int)
        lookup.fill(-1)
        lookup[self.face_id_by_edge] = self.edges.ids
        return lookup

    def boundary(self, cell_id):
        edge = self._edge_lookup[cell_id]
        if edge == -1:
            return [-1]
            os._exit(1) #force end program because edge is out of the edges ids #Pilar have added it!!!
        if edge != -1: 
            return cycle(self.edges.next, edge)

    def moved(self, dv):
        """Caller is responsible for making dv consistent for different copies of a vertex."""
        mesh = self.copy()
        mesh.vertices = mesh.vertices + dv
        return mesh

    def transition(self, eps):
        return _transition(self, eps)

    def add_edges(self, edge_pairs):
        return _add_edges(self, edge_pairs)


def _remove(edges, reverse, vertices, face_id_by_edge):
    """Removes the given edges and their orbit under rotation (ie removes complete vertices.)

    Args:
        edges: An array of edge IDs to remove
        reverse: see mesh.edges.reverse
        vertices: see mesh.vertices
        face+ids: see mesh.face_ids

    Returns:
        A tuple of reverse,vertex,cell arrays with the edges removed.
    """
    es = np.unique(edges//3*3)
    to_del = np.dstack([es, es+1, es+2]).ravel()  # orbit of the given edges under rotation
    reverse = np.delete(reverse, to_del)
    reverse = (np.cumsum(np.bincount(reverse))-1)[reverse]  # relabel to get a perm of [0,...,N-1]
    vertices = np.delete(vertices, to_del, 1)
    face_id_by_edge = np.delete(face_id_by_edge, to_del)
    return reverse, vertices, face_id_by_edge, to_del


# def sum_vertices(edges,dv):
#    return np.repeat(dv[:,::3]+dv[:,1::3]+dv[:,2::3],3,1)
def sum_vertices(edges, dv):
    return dv + np.take(dv, edges.rotate, 1) + np.take(dv, edges.rotate2, 1)


def _T1(edge, eps, rotate, reverse, vertices, face_id_by_edge):
    e0 = edge
    e1 = rotate[edge]
    e2 = rotate[e1]
    e3 = reverse[edge]
    e4 = rotate[e3]
    e5 = rotate[e4]

    before = np.array([e1, e2, e4, e5])
    after = np.array([e2, e4, e5, e1])

    after_r = reverse.take(after)
    reverse[before] = after_r
    reverse[after_r] = before

    dv = vertices[:, e4]-vertices[:, e0]
    l = 1.01*eps/np.sqrt(dv[0]*dv[0]+dv[1]*dv[1])
    dw = [dv[1]*l, -dv[0]*l]

    # for i in [0, 1]:
    #     dp = 0.5*(dv[i]+dw[i])
    #     dq = 0.5*(dv[i]-dw[i])
    #     v = vertices[i]
    #     v[before] = v.take(after) + np.array([dp, -dq, -dp, dq])
    #     v[e0] = v[e4] + dw[i]
    #     v[e3] = v[e1] - dw[i]

    # change the coordinates of the vertices
    for i in [0, 1]:
        dp = 0.5*(dv[i]+dw[i])
        dq = 0.5*(dv[i]-dw[i])
        vertices[i,before] = vertices[i,after] + np.array([dq, -dp, -dq, dp])
        vertices[i,e0] = vertices[i,e4] - dw[i] # e4 = next[e0]
        vertices[i,e3] = vertices[i,e1] + dw[i] # e1 = next[e3]

    face_id_by_edge[before] = face_id_by_edge.take(after)
    face_id_by_edge[e0] = face_id_by_edge[e4]
    face_id_by_edge[e3] = face_id_by_edge[e1]

    return np.hstack([before, after_r])


def _transition(mesh, eps):
    edges = mesh.edges
    half_edges = edges.ids[edges.ids < edges.reverse]  # half the edges have id < reverse_id
    dv = mesh.edge_vect.take(half_edges, 1)
    short_edges = set(half_edges[np.sum(dv*dv, 0) < eps*eps])
    ids_t1=half_edges[np.sum(dv*dv, 0) < eps*eps]
    
    if not short_edges:
        return mesh, ids_t1, [] #####BAETTI THESSU VID
    reverse, vertices, face_id_by_edge = edges.reverse.copy(), mesh.vertices.copy(), mesh.face_id_by_edge.copy()
    rotate = edges.rotate
    # Do T1 transitions
    # to avoid nasty edge cases, we don't allow T1's to happen on adjacent edges
    # and delay to the next timestep if necessary.
    # A better approach would be to take multiple partial timesteps.
    boundary_edges = mesh.boundary_edges if mesh.has_boundary() else []
    while short_edges:
        edge = short_edges.pop()
        if edge in boundary_edges:
            edge = reverse[edge]
        neighbours = _T1(edge, eps, rotate, reverse, vertices, face_id_by_edge)
        for x in neighbours:
            short_edges.discard(x) 

    nxt = rotate[reverse]
    #c=[]
    # Remove collapsed (ie two-sided) faces.
    edg_rem=[]
    while True:
        nxt = rotate[reverse]
        #c.append()
        two_sided = np.where(nxt[nxt] == edges.ids[:len(nxt)])[0]
        if not len(two_sided):
            break
        while np.any(reverse[reverse[rotate[two_sided]]] != reverse[rotate[nxt[two_sided]]]):
            reverse[reverse[rotate[two_sided]]] = reverse[rotate[nxt[two_sided]]]
        prev_face_id_by_edge = face_id_by_edge
        reverse, vertices, face_id_by_edge, to_del = _remove(two_sided, reverse, vertices, face_id_by_edge)
        ids_removed = np.setdiff1d(prev_face_id_by_edge,face_id_by_edge)
        edg_rem.append(to_del)
        #print to_del, ids_removed
        # if ~(ids_t1==np.delete(ids_t1,ids_removed)):
        #     print 'Ids T1 to remove:', ids_t1, ids_removed, np.delete(ids_t1,ids_removed)
        # print(ids_removed, to_del)
    mesh = mesh.copy()
    mesh.edges = Edges(reverse)
    mesh.vertices = vertices
    mesh.face_id_by_edge = face_id_by_edge
    return mesh, ids_t1, edg_rem  #######BAETTI THESSU VID


def _add_edges(mesh, edge_pairs):
    n_edge_old = len(mesh.edges)
    n_edge_new = n_edge_old + 6*len(edge_pairs)

    n_face_old = mesh.n_face

    reverse = np.resize(mesh.edges.reverse, n_edge_new)
    vertices = np.empty((2, n_edge_new))
    vertices[:, :n_edge_old] = mesh.vertices
    face_id_by_edge = np.resize(mesh.face_id_by_edge, n_edge_new)
    rotate = Edges.ROTATE

    v = vertices.T  # easier to work with transpose here
    n = n_edge_old
    for i, (e1, e2) in enumerate(edge_pairs):
        # find midpoints: rotate[reverse] = next
        v1, v2 = 0.5*(v.take(rotate[reverse[[e1, e2]]], 0) + v.take([e1, e2], 0))
        v[[n, n+1, n+2]] = v1
        v[[n+3, n+4, n+5]] = v2

        a = [n, n+1, n+2, n+3, n+5]
        b = [e1, n+4, reverse[e1], e2, reverse[e2]]
        reverse[a], reverse[b] = b, a

        for j, edge in enumerate((e1, e2)):
            face_id = n_face_old + 2*i + j
            while face_id_by_edge[edge] != face_id:
                face_id_by_edge[edge] = face_id
                edge = rotate[reverse[edge]]  # rotate[reverse] = next

        re1, re2 = rotate[[e1, e2]]
        # winding
        face_id_by_edge[n] = face_id_by_edge[re1]
        v[n] += v[re1]-v[e1]
        face_id_by_edge[n+3] = face_id_by_edge[re2]
        v[n+3] += v[re2]-v[e2]
        n += 6

    mesh = mesh.copy()
    mesh.edges = Edges(reverse)
    mesh.face_id_by_edge = face_id_by_edge
    mesh.vertices = vertices
    mesh.n_face = n_face_old + 2*len(edge_pairs)
    return mesh

# <codecell>


# <codecell>


