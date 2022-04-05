# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 11:49:11 2017

@author: otteflor
"""

import syris
import matplotlib.pyplot as plt
import numpy as np
import quantities as q
from syris.geometry import Trajectory
from syris.devices.sources import BendingMagnet
from syris.bodies.mesh import make_cube

syris.init(platform_name = 'Intel(R) OpenCL', device_index = 1)

n = 512
shape = (n, n)
ps = 1 * q.um
dE = 1 * q.keV
energies = np.arange(5, 30, dE.magnitude) * q.keV
tr = Trajectory([(n/2, n/2, 0)] * ps, velocity=10 * q.um / q.s, pixel_size=ps)
bm = BendingMagnet(1.5 * q.GeV, 100 * q.mA, 7 * q.T, 30 * q.m, dE, (200, 800) * q.um, ps, tr)

mesh = make_cube()

