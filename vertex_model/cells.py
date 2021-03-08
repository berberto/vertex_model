# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import numpy as np


class Cells(object):
    def __init__(self, mesh, properties=None):
        self.mesh = mesh
        self.properties = properties or {}

    def copy(self):
        # mesh is 'immutable' (so we can cache shared computations) => no need to copy it
        return Cells(self.mesh, self.properties.copy())

    def __len__(self):
        return self.mesh.n_face

    def empty(self):
        # self.mesh._edge_lookup[self.mesh.area<0]=-1 ###########anyadido por mi!!!!!!!!! 
        return self.mesh._edge_lookup == -1

    def by_face(self, property_name, boundary_property_name=None):
        value = self.properties[property_name]
        if self.mesh.has_boundary():
            value = make_array(value, len(self))
            boundary_value = self.properties[boundary_property_name] if boundary_property_name else 0.0
            value[self.mesh.boundary_faces] = boundary_value
        return value

    def by_edge(self, property_name, boundary_property_name=None):
        value_by_face = self.by_face(property_name, boundary_property_name)
        if not hasattr(value_by_face, 'shape'):  # scalar
            return value_by_face
        return value_by_face.take(self.mesh.face_id_by_edge)


def make_array(value, n):
    if hasattr(value, 'shape'):
        return value
    expanded = np.empty(n, type(value))
    expanded.fill(value)
    return expanded

