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
import tifffile as tf


ext_folder = '/home/ubuntu/Dropbox/MA/TEIL02/extensions/extensions'

if ext_folder not in sys.path:
    sys.path.insert(0, ext_folder)
from extensions import read_wavefront_obj, TomoExperiment

OUTPUT = 'raw/phantom'
OBJ_PATH = 'sericostoma.obj'
NO_OF_IMAGES = 180
THETA_MIN = 0
THETA_MAX = 180
START_I = 71
n = 2048


LOG = logging.getLogger(__name__)
syris.init( )

shape = (n, n)
energies = np.arange(6, 30, 1) * q.keV
dE = energies[1] - energies[0]

# === CAMERA ===
vis_wavelengths = np.arange(500, 700) * q.nm
camera = Camera(11 * q.um, .1, 500, 23, 32, shape, exp_time=500 * q.ms, fps=2 / q.s,
                quantum_efficiencies=0.5 * np.ones(len(vis_wavelengths)),
                wavelengths=vis_wavelengths, dtype=np.float32)

# === LENS ===
lens = Lens(9, f_number=1.4, focal_length=50 * q.mm, transmission_eff=0.7, sigma=None)

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
meshes = read_wavefront_obj(OBJ_PATH,
                            1 * q.mm,
                            [(n / 2, n / 2, 0)] * detector.pixel_size)

tr = Trajectory([(n / 2, n / 2, 0)] * detector.pixel_size)
cb = CompositeBody(tr, bodies = meshes)
cb.bind_trajectory(detector.pixel_size)

# === MAKE EXPERIMENT ===
ex = TomoExperiment( cb, detector, 0 * q.cm, energies)


# == CONDUCT EXPERIMENT
if OUTPUT is not None and not os.path.exists( OUTPUT ):
    os.makedirs( OUTPUT, mode=0o755)

t_0 = 0 * q.s
t_1 = NO_OF_IMAGES / detector.camera.fps
st = time.time()
mpl_im = None

if START_I > 0:
    lin = np.linspace(THETA_MIN, THETA_MAX, NO_OF_IMAGES)
    THETA_MIN = lin[START_I]
    NO_OF_IMAGES = NO_OF_IMAGES - START_I + 1

for i, proj in enumerate(ex.make_tomography( NO_OF_IMAGES, THETA_MIN, THETA_MAX, 1*q.s)):

    if START_I > 0:
        i = i + START_I

    image = get_host(proj)
    msg = '===== COMPUTED PROJECTION # {}'
    LOG.error(msg.format(i))
    show(image)
    plt.show()

    if OUTPUT:
        path = os.path.join( OUTPUT, 'proj_{:>05}.tif').format(i)
        #sc.imsave(path, image)
	tf.imsave( path, image.astype(np.float32) )
