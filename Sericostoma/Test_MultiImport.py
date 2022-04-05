# -*- coding: utf-8 -*-
"""
Created on Sun Apr 02 00:08:17 2017

@author: hambu
"""

import os, sys
import time
import matplotlib.pyplot as plt
import numpy as np
import quantities as q
import scipy.misc as sc
import syris


import logging

from syris.physics import propagate
from syris.bodies.base import CompositeBody
from syris.bodies.mesh import Mesh, read_blender_obj
from syris.devices.cameras import Camera
from syris.devices.detectors import Detector
from syris.devices.lenses import Lens
from syris.devices.filters import Scintillator
from syris.devices.sources import make_topotomo
from syris.geometry import Trajectory
from syris.gpu.util import get_host
from syris.math import fwnm_to_sigma
from util import get_material, show

#ext_folder = '/home/ubuntu/Dropbox/MA/TEIL02/extensions/extensions'
ext_folder = 'C:\Users\hambu\Dropbox\MA\TEIL02\extensions\extensions'

if ext_folder not in sys.path:
    sys.path.insert(0, ext_folder)
#from extensions import read_wavefront_obj, TomoExperiment
from importmulti import read_collada

OUTPUT = 'raw/phantom'
OBJ_PATH = 'dummy.obj'
NO_OF_IMAGES = 5
THETA_MIN = 0 
THETA_MAX = 180 
n = 256


LOG = logging.getLogger(__name__)
syris.init(platform_name = "Intel(R) OpenCL", device_index = 0) 

shape = (n, n)
energies = np.arange(6, 30, 1) * q.keV
dE = energies[1] - energies[0]                    

# === CAMERA ===
vis_wavelengths = np.arange(500, 700) * q.nm
camera = Camera(12 * q.um, .1, 500, 23, 32, shape, exp_time=500 * q.ms, fps=2 / q.s,
                quantum_efficiencies=0.5 * np.ones(len(vis_wavelengths)),
                wavelengths=vis_wavelengths, dtype=np.float32)

# === LENS ===
lens = Lens(1, f_number=1.4, focal_length=50 * q.mm, transmission_eff=0.7, sigma=None)

# === SCINTILLATOR ===
x = vis_wavelengths.rescale(q.nm).magnitude
dx = x[1] - x[0]
sigma = fwnm_to_sigma(50)
emission = np.exp(-(x - 450) ** 2 / (2 * sigma ** 2)) / (sigma * np.sqrt(2 * np.pi)) * dx
luag = get_material('luag.mat')
scintillator = Scintillator(50 * q.um,
                            luag,
                            14 * np.ones(len(energies)) / q.keV,
                            energies,
                            emission / q.nm,
                            vis_wavelengths,
                            1.84)

# === DETECTOR ===
detector = Detector(scintillator, lens, camera)

# === BENDING MAGNET SOURCE ===
source_trajectory = Trajectory([(n / 2, n / 2, 0)] * detector.pixel_size)
bm = make_topotomo(dE=dE, trajectory=source_trajectory, pixel_size=detector.pixel_size)
bm.trajectory.bind(detector.pixel_size)

# == SAMPLE ===
tris = read_collada('C:\Users\hambu\Desktop\untitled.dae', 
                    units = 1 * q.mm)
tris2 = read_blender_obj('C:/Users/hambu/Desktop/box.obj')

m = get_material("glass.mat")
mesh_col = Mesh(tris * q.mm, source_trajectory, material = m)
mesh_obj = Mesh(tris2 * q.mm, source_trajectory, material = m)
#mesh = meshes[0]


ps = 12*q.um
meshes = [mesh_col, mesh_obj]

for mesh in meshes:
    #mesh.bind_trajectory(ps)
    mesh.translate([n / 2, n/2, 0]* ps)


    
shape = (n,n)
#u = mesh.transfer((256, 256), 12 * q.um, 10 * q.keV).get()
pc = propagate([mesh_col], shape, [10*q.keV], 0 *q.m, ps, t = 0*q.s, detector = detector).get()
#po = propagate([mesh_obj], shape, [10*q.keV], 0 *q.m, ps, t = None, detector = detector).get()
#show(u.real)
#show(u.imag)
show(pc)
#show(po)
plt.show()