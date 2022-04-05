"""Show different propagators."""
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import quantities as q
from numpy.fft import fftfreq
import syris
from syris.physics import energy_to_wavelength, compute_propagator
from syris.bodies.mesh import make_cube, Mesh
from syris.bodies.isosurfaces import MetaBall 
from util import get_material
from syris.geometry import Trajectory
from syris.imageprocessing import fft_2, ifft_2
from util import show


def compute_fourier_propagator(n, lam, z, ps, fresnel=True):
    try:
        ps[0]
    except:
        ps = (ps.magnitude, ps.magnitude) * ps.units

    lam = lam.rescale(q.m).magnitude
    z = z.rescale(q.m).magnitude
    ps = ps.rescale(q.m).magnitude

    freqs = fftfreq(n)
    f = np.tile(freqs, [n, 1])
    g = np.copy(f.transpose())
    f /= ps[1]
    g /= ps[0]

    if fresnel:
        result = np.exp(1j * 2 * np.pi / lam * z) * \
            np.exp(-1j * np.pi * lam * z * (f ** 2 + g ** 2))
    else:
        result = np.exp(1j * 2 * np.pi / lam * z * np.sqrt(1 - (f * lam) ** 2 - (g * lam) ** 2))

    return result


def main():
    """main"""
    syris.init(platform_name = 'Intel(R) OpenCL', device_index = 1, double_precision=False)
    n = 4096
    ps = 512. / 4096 * 0.5 * q.um
    energy = 10 * q.keV
    lam = energy_to_wavelength(energy)
    # Compute the sampling limit for given n, ps and lam
    ca = (lam / 2 / ps).simplified.magnitude
    tga = np.tan(np.arccos(ca))
    distance = (tga * n * ps / 2).simplified
    print 'Propagation distance:', distance
    
    cub = make_cube() / q.m * n *0.1 * ps
    cmat = get_material( 'glass.mat' )
    ctra = Trajectory([(n/2, n/2, 0)] * ps, ps)
    cube = Mesh(cub, ctra, cmat)
    cube.translate((n/2, n/2, 0) * ps)
    sphere = MetaBall(ctra, n *0.1 * ps, cmat)
    u = cube.transfer(shape = (n,n), pixel_size = ps, energy = energy)
    distances = np.linspace(0, 100, 12) * q.cm
    fig = plt.figure() 
    for i, dis in enumerate(distances):
        propagator = compute_propagator(n, dis, lam, ps, fresnel = True, mollified = False)
        fft_2(u)
        u *= propagator
        ifft_2(u)
        intensity = u.real
        fig.add_subplot( 4,3,i+1 )
        plt.imshow(intensity.get(), cmap = "gray")
        plt.title(distances[i])
        
    plt.colorbar()
    plt.show()

    


if __name__ == '__main__':
    main()
