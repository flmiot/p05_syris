# -*- coding: utf-8 -*-
"""


@author: otteflor
"""

OUTPUT = 'raw/phantom'
OBJ_PATH = 'phantom_alex.obj'
NO_OF_IMAGES = 360
n = 1024

import os
import time
import matplotlib.pyplot as plt
import numpy as np
import quantities as q
import scipy.misc as sc
import syris

import pyopencl.array as cl_array
import logging
import syris.config as cfg
import syris.math as smath
import syris.imageprocessing as ip
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
from syris.experiments import Experiment
from syris.math import fwnm_to_sigma
from util import get_material, show
import syris.geometry as geom
from extensions import read_wavefront_obj
import tifffile as tf


import logging
LOG = logging.getLogger(__name__)
syris.init( )

class PetraExperiment(object):

    """A virtual synchrotron experiment base class."""

    def __init__(self, samples, source, detector, propagation_distance, energies):
        self.source = source
        self.samples = samples
        self.detector = detector
        self.propagation_distance = propagation_distance
        self.energies = energies
        self._time = None

    @property
    def time(self):
        """Total time of all samples."""
        if self._time is None:
            self._time = max([obj.trajectory.time for obj in self.samples
                              if obj.trajectory is not None])

        return self._time

    def get_next_time(self, t, pixel_size):
        """Get next time from *t* for all the samples."""
        return min([obj.get_next_time(t, pixel_size) for obj in self.samples])

    def make_source_blur(self, shape, pixel_size, queue=None, block=False):
        """Make geometrical source blurring kernel with *shape* (y, x) size and *pixel_size*. Use
        OpenCL command *queue* and *block* if True.
        """
        l = self.source.sample_distance
        size = self.source.size
        width = (self.propagation_distance * size[1] / l).simplified.magnitude
        height = (self.propagation_distance * size[0] / l).simplified.magnitude
        sigma = (smath.fwnm_to_sigma(height, n=2), smath.fwnm_to_sigma(width, n=2)) * q.m

        return ip.get_gauss_2d(shape, sigma, pixel_size=pixel_size, fourier=True,
                               queue=queue, block=block)

    def compute_intensity(self, t_0, t_1, shape, pixel_size, queue=None, block=False):
        """Compute intensity between times *t_0* and *t_1*."""
        exp_time = (t_1 - t_0).simplified.magnitude
        image = propagate(self.samples, shape, self.energies, self.propagation_distance,
                          pixel_size, detector=self.detector, t=None) * exp_time

        return image

    def setup_sample(self, centerpoint = None):
        
        psm = self.detector.pixel_size.simplified.magnitude
        point = (n * psm / 2, n * psm / 2, 0) * q.m
        if centerpoint:
            point = centerpoint  
        for ind, obj in enumerate(self.samples):
            if obj != bm:
                self.samples[ind].translate(point) 
                    
        

        
    def make_sequence(self, t_start, t_end, shape=None, shot_noise=True, amplifier_noise=True,
                      source_blur=True, queue=None):
        """Make images between times *t_start* and *t_end*."""
        if queue is None:
            queue = cfg.OPENCL.queue
        shape_0 = self.detector.camera.shape
        if shape is None:
            shape = shape_0
        ps_0 = self.detector.pixel_size
        ps = shape_0[0] / float(shape[0]) * ps_0
        fps = self.detector.camera.fps
        frame_time = 1 / fps
        times = np.arange(t_start.simplified.magnitude, t_end.simplified.magnitude,
                          frame_time.simplified.magnitude) * q.s
        image = cl_array.Array(queue, shape, dtype=cfg.PRECISION.np_float)
        source_blur_kernel = None
        if source_blur:
            source_blur_kernel = self.make_source_blur(shape, ps, queue=queue, block=False)

        fmt = 'Making sequence with shape {} and pixel size {} from {} to {}'
        LOG.debug(fmt.format(shape, ps, t_start, t_end))

        for i, t_0 in enumerate(times):
            image.fill(0)
            t = t_0
            t_next = self.get_next_time(t, ps)
            
            # Turn sample
            for ind, obj in enumerate(self.samples):
                if isinstance(obj, CompositeBody):
                    self.samples[ind].rotate(0.5*q.deg , geom.Y_AX)    
            
            
            while t_next < t_0 + frame_time:
                LOG.debug('Motion blur: {} -> {}'.format(t, t_next))
                image += self.compute_intensity(t, t_next, shape, ps)
                t = t_next
                t_next = self.get_next_time(t, ps)
            image += self.compute_intensity(t, t_0 + frame_time, shape, ps)
            if source_blur:
                image = ip.ifft_2(ip.fft_2(image) * source_blur_kernel).real
            camera_image = self.detector.camera.get_image(image, shot_noise=shot_noise,
                                                          amplifier_noise=amplifier_noise)
            LOG.debug('Image: {} -> {}'.format(t_0, t_0 + frame_time))
            yield camera_image




shape = (n, n)
energies = np.arange(6, 30, 1) * q.keV
dE = energies[1] - energies[0]                    

# === CAMERA ===
vis_wavelengths = np.arange(500, 700) * q.nm
camera = Camera(11 * q.um, .1, 500, 23, 32, shape, exp_time=25 * q.ms, fps=2000 / q.s,
                quantum_efficiencies=0.5 * np.ones(len(vis_wavelengths)),
                wavelengths=vis_wavelengths, dtype=np.float32)

# === LENS ===
lens = Lens(10, f_number=1.4, focal_length=50 * q.mm, transmission_eff=0.7, sigma=None)

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

# == SAMPLE ===
meshes = read_wavefront_obj(OBJ_PATH, 
                            0.5 * q.mm, 
                            [(n / 2, n / 2, 0)] * detector.pixel_size) 

tr = Trajectory([(n / 2, n / 2, 0)] * detector.pixel_size)
cb = CompositeBody(tr, bodies = meshes)
cb.translate([n / 2, n / 2, 0] * detector.pixel_size)
cb.bind_trajectory(detector.pixel_size)

# === MAKE EXPERIMENT ===
ex = PetraExperiment( [bm], bm, detector, 0 * q.cm, energies)   


            
# == CONDUCT EXPERIMENT
if OUTPUT is not None and not os.path.exists( OUTPUT ):
    os.makedirs( OUTPUT, mode=0o755)
    
t_0 = 0 * q.s
t_1 = NO_OF_IMAGES / detector.camera.fps
st = time.time()
mpl_im = None

for i, proj in enumerate(ex.make_sequence(t_0, t_1)):
    image = get_host(proj)
    msg = '===== COMPUTED PROJECTION # {}'
    LOG.error(msg.format(i))      
    #show(image)
    #plt.show()

    if OUTPUT:
        path = os.path.join( OUTPUT, 'flat_{:>05}.tif').format(i)
        tf.imsave( path, image.astype(np.float32) )
