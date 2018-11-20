from __future__ import division

import numpy as np

from matscipy.neighbours import neighbour_list
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.calculator import PropertyNotImplementedError

# default parameters for the Morse potential come from the following article:
# De Wette et al. Physics Letters vol 23 (1966) http://www.sciencedirect.com/science/article/pii/0031916366900138

class MorsePotential(Calculator):
    implemented_properties = ['energy', 'forces', 'stress', 'local_energy']
    default_parameters = {'D': 0.16156,
                          'alpha': 2.0926,
                          'r0': 2.6163,
                          'rc': 8.0}
    nolabel = True

    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)

    def calculate(self, atoms=None,
                  properties=['energy', 'local_energy'],
                  system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        if atoms is None:
            atoms = self.atoms
        natoms = len(atoms)

        D = self.parameters.D
        alpha = self.parameters.alpha
        r0 = self.parameters.r0
        rc = self.parameters.rc

        i, j, dr, r = neighbour_list('ijDd', atoms, rc)
        tmp = np.exp(-alpha*(r-r0))
        de =  0.5 * D * ( tmp*tmp - 2.0*tmp )
        local_energy = np.bincount(i, de, minlength=natoms)
        energy = local_energy.sum()

        tmpd = -alpha * tmp
        df =  (D * ( 2.0*tmp*tmpd - 2.0 * tmpd) / r)[:, np.newaxis] * dr
        forces = np.zeros((natoms, 3))
        for kk in range(3):
            forces[:, kk] = np.bincount(i, weights=df[:, kk],
                                        minlength=natoms)

        if 'stress' in properties:
            if self.atoms.number_of_lattice_vectors == 3:
                virial = np.zeros(6)
                if len(i) > 0:
                    virial = 0.5*np.array([dr[:,0]*df[:,0],               # xx
                                           dr[:,1]*df[:,1],               # yy
                                           dr[:,2]*df[:,2],               # zz
                                           dr[:,1]*df[:,2],               # yz
                                           dr[:,0]*df[:,2],               # xz
                                           dr[:,0]*df[:,1]]).sum(axis=1)  # xy
                self.results['stress'] = virial/atoms.get_volume()
            else:
                raise PropertyNotImplementedError

        self.results['energy'] = energy
        self.results['free_energy'] = energy
        self.results['forces'] = forces
        self.results['local_energy'] = local_energy
