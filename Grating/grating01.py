# -*- coding: utf-8 -*-
"""


@author: otteflor
"""

import os, sys
import time
import matplotlib.pyplot as plt
import numpy as np
import quantities as q
import scipy.misc as sc
import syris


import logging

from syris.physics import propagate, transfer_many
from syris.devices.sources import make_topotomo
from syris.bodies.base import CompositeBody
from syris.bodies.mesh import Mesh, make_cube, read_wavefront_obj
from syris.devices.filters import Grating
from syris.geometry import Trajectory
from syris.gpu.util import get_host
from syris.physics import energy_to_wavelength
from syris.math import fwnm_to_sigma
from util import get_material, show
from syris.materials import make_henke
import tifffile as tf

OUTPUT = None
ss = 1
n = 1024 * ss
ps = .1*q.um / ss

LOG = logging.getLogger(__name__)
syris.init()
shape = (n,n)

energies = np.arange(6, 30, 1) * q.keV
dE = energies[1] - energies[0]

# === SAMPLE ===
meshes = read_wavefront_obj("grating_mesh.obj",
                            1*q.um,
                            [(n / 2, n / 2, 0)] * ps,
                            iterations = 2)
tr = Trajectory([(n / 2, n / 2, 0)] * ps)
cb = CompositeBody(tr, bodies = meshes)
cb.bind_trajectory(ps)
cb.translate((n / 2, n / 2, 0) * ps)


cube_tris = make_cube() / q.m * 300 * q.um
#mat_nickel = make_henke('Nickel', energies, density = 7.81 *q.g / (q.cm)**3, formula = 'Ni')
#mat_nickel.save('data/nickel.mat')
mat = get_material("nickel.mat")
cube = Mesh(cube_tris, trajectory = tr, material = mat)
cube.trajectory.bind(ps)
cube.translate((n / 2, n / 2, 0) * ps)

grating_period = 2 * q.um
grating = Grating(100 * q.um, mat, 0.5, 32 * ss)

# === BENDING MAGNET SOURCE ===

bm = make_topotomo(dE=dE, trajectory=tr, pixel_size=ps)
bm.trajectory.bind(ps)

energy = 8 * q.keV
lam = energy_to_wavelength(energy)
talbot_distance = lam / (1-np.sqrt(1-lam**2 / grating_period**2))
print(talbot_distance.simplified)
distances = [0, 0.25, 0.5, 0.75, 1]* talbot_distance
for d in distances:
  intensity = propagate([bm, cb], shape, [energy], d*q.m, ps).get()
  show(intensity)
  plt.show()
#amplitude = transfer_many([bm, grating], shape, ps, 20* q.keV, exponent = True).get()
#amplitude = bm.transfer(n, ps, 20 * q.keV, exponent = True).get()



if OUTPUT:
  path = os.path.join( OUTPUT, 'proj_{:>05}.tif').format(i)
  #sc.imsave(path, image)
  tf.imsave( path, image.astype(np.float32) )
