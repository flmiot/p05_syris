import numpy as np
import syris
import quantities as q
from syris.materials import make_henke


syris.init(platform_name = 'Intel(R) OpenCL', device_index = 1)
energies = np.arange(6, 30, 1) * q.keV
mat_polye = make_henke('Polyethylene', energies, density = 0.94 *q.g / (q.cm)**3, formula = 'C2H4')
mat_polys = make_henke('Polystyrene', energies, density = 1.00 *q.g / (q.cm)**3, formula = 'C8H8')
mat_pa6 = make_henke('Polyamid-6', energies, density = 1.084 *q.g / (q.cm)**3, formula = 'C6H11NO')
mat_pmma = make_henke('PMMA', energies, density = 1.18 *q.g / (q.cm)**3, formula = 'C5O2H8')

mat_polye.save('data/polyethylene.mat')
mat_polys.save('data/polystyrene.mat')
mat_pa6.save('data/polycaprolactam.mat')
mat_pmma.save('data/pmma.mat')