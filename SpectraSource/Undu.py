# -*- coding: utf-8 -*-
"""


@author: otteflor
"""

import syris
import matplotlib.pyplot as plt
import numpy as np
import quantities as q
from syris.geometry import Trajectory
from syris.devices.sources import SpectraSource
from syris.physics import propagate
from util import show

syris.init()
filename = '/home/ubuntu/Dropbox/MA/TEIL02/Spectra/syris_undulator/en7859p04_60m_256a11um_k1p0.dta'

n = 256
shape = (n, n)
ps = 11 * q.um
energy = [7859.04] * q.keV
tr = Trajectory([(n/2, n/2, 0)] * ps, velocity=10 * q.um / q.s, pixel_size=ps)
source = SpectraSource(filename, 60 * q.m, (500,500)* q.um, ps, tr)
intensity = propagate([source], shape, energy, 0*q.m, ps)
show(intensity.get())
plt.show()
