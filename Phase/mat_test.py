import numpy as np
import syris
import quantities as q
from syris.materials import make_henke, Material


syris.init()
energies = np.arange(6, 30, 1) * q.keV
#mat_chitin = make_henke('Chitin', energies, density = 1.37 *q.g / (q.cm)**3, formula = 'C8H13NO5')
#mat_water = make_henke('Water', energies, density = 1.00 *q.g / (q.cm)**3, formula = 'H2O')
#mat_au = make_henke('Aluminium', energies, density = 19.3 *q.g / (q.cm)**3, formula = 'Au')
#mat_polye = make_henke('Polyethylene', energies, density = 0.94 *q.g / (q.cm)**3, formula = 'C2H4')
#mat_polys = make_henke('Polystyrene', energies, density = 1.00 *q.g / (q.cm)**3, formula = 'C8H8')
#mat_pa6 = make_henke('Polyamid-6', energies, density = 1.084 *q.g / (q.cm)**3, formula = 'C6H11NO')
#mat_pmma = make_henke('PMMA', energies, density = 1.18 *q.g / (q.cm)**3, formula = 'C5O2H8')
#mat_al = make_henke('Aluminium', energies, density = 2.7 *q.g / (q.cm)**3, formula = 'Al')
mat_urethan = make_henke('Urethan', energies, density = 0.98*q.g / (q.cm)**3, formula = 'C3H7NO2')
mat_iodine = make_henke('Iodine', energies, density = 4.94*q.g / (q.cm)**3, formula = 'I')
#mat_polye.save('data/polyethylene.mat')
#mat_polys.save('data/polystyrene.mat')
#mat_pa6.save('data/polycaprolactam.mat')
#mat_pmma.save('data/pmma.mat')
#mat_chitin.save('data/chitin.mat')
#mat_water.save('data/water.mat')
#mat_au.save('data/au.mat')
#mat_al.save('data/al.mat')
mat_urethan.save('data/urethan.mat')
mat_iodine.save('data/iodine.mat')

# create 12 mixtures of urethan and iodine
for i in np.arange(1,13):
    temp_mat = Material("Mixture", mat_urethan.refractive_indices * (1 - i * 0.005) + mat_iodine.refractive_indices * i * 0.005, energies) 
    filename = 'data/urethane_iodine_{:>05}.mat'.format(i * 5)
    temp_mat.save(filename)

