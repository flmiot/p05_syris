#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 14:19:40 2017

@author: ubuntu
"""
import re
import numpy as np
import quantities as q

from syris.bodies.mesh import Mesh
from syris.geometry import Trajectory
from util import get_material


def read_blender_obj(filename, objects=None):
    """Read blender wavefront *filename*, extract only *objects* which are object indices."""
    remainder = open(filename, 'r').read()
    triangles = None
    face_start = 0
    i = 0

    while remainder:
        remainder, v, f, m = _extract_object2(remainder)
        if objects is None or i in objects:
            if triangles is None:
                triangles = v[f - face_start].transpose()
            else:
                triangles = np.concatenate((triangles, v[f - face_start].transpose()), axis=1)
        face_start += len(v)
        i += 1

    return triangles

def _extract_object2(txt):
    """Extract an object from string *txt*."""
    face_start = txt.index('usemtl ')
    if 'v ' not in txt[face_start:]:
        obj_end = None
    else:
        obj_end = face_start + txt[face_start:].index('v ')
    subtxt = txt[:obj_end]

    pattern = r'{} (?P<x>.*) (?P<y>.*) (?P<z>.*)'
    v_pattern = re.compile(pattern.format('v'))
    f_pattern = re.compile(pattern.format('f'))
    vertices = np.array(re.findall(v_pattern, subtxt)).astype(np.float32)
    faces = np.array(re.findall(f_pattern, subtxt)).astype(np.int32).flatten() - 1
    material = re.findall('usemtl .*', subtxt)[0][7:]

    remainder = txt[obj_end:] if obj_end else None

    return remainder, vertices, faces, material

def read_wavefront_obj( filename, scale, center ):
    """ Returns list with objects of type syris.bodies.mesh"""
    remainder = open(filename, 'r').read()
    face_start = 0
    i = 0
    Meshes = []

    while remainder:
        triangles = None
        remainder, v, f, mat = _extract_object2( remainder )
        if triangles is None:
            triangles = v[f - face_start].transpose()
        else:
            triangles = np.concatenate((triangles, v[f - face_start].transpose()), axis=1)
        tr = Trajectory( center )
        m = get_material( mat )
        c = [(0,0,0)] * q.m
        Meshes.extend([ Mesh(triangles * scale, tr, material = m, center = c) ])
        face_start += len(v)
        i += 1

    return Meshes
    