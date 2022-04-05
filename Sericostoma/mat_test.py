import numpy as np
import syris
import quantities as q
from syris.materials import make_henke


syris.init()
energies = np.arange(6, 30, 1) * q.keV
mat_chitin = make_henke('Chitin', energies, density = 1.37 *q.g / (q.cm)**3, formula = 'C8H13NO5')
mat_water = make_henke('Water', energies, density = 1.00 *q.g / (q.cm)**3, formula = 'H2O')

mat_chitin.save('data/chitin.mat')
mat_water.save('data/water.mat')