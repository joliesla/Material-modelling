from __future__ import division

import numpy as np

from ase.neighborlist import NeighborList
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

    def morse_pair_energy(self,r):
        if r > self.parameters.rc:
            return 0.0
        tmp = np.exp(-self.parameters.alpha*(r-self.parameters.r0))
        return  self.parameters.D * ( tmp*tmp - 2.0*tmp )

    def morse_pair_energy_deriv(self,r):
        if r > self.parameters.rc:
            return 0.0
        tmp = np.exp(-self.parameters.alpha*(r-self.parameters.r0))
        tmpd = -self.parameters.alpha * tmp
        return self.parameters.D * ( 2.0*tmp*tmpd - 2.0 * tmpd)
    
    def calculate(self, atoms=None,
                  properties=['energy', 'local_energy'],
                  system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
            
        natoms = len(self.atoms)
            
        D = self.parameters.D
        alpha = self.parameters.alpha
        r0 = self.parameters.r0
        rc = self.parameters.rc
                
        if 'numbers' in system_changes:
            self.nl = NeighborList([rc / 2.0] * natoms, skin = 2.0, self_interaction=False)
                    
        self.nl.update(self.atoms)
                    
        positions = self.atoms.positions
        cell = self.atoms.cell
                
        e0 = self.morse_pair_energy(rc)
                    
        energy = 0.0
        local_energy = np.zeros((natoms,1))
        forces = np.zeros((natoms, 3))
        stress = np.zeros((3, 3))

        m = lambda x: self.morse_pair_energy(x)
        mvec = np.vectorize(m)
        md = lambda x: self.morse_pair_energy_deriv(x)
        mdvec = np.vectorize(md)
        
        for a1 in range(natoms):
            neighbors, offsets = self.nl.get_neighbors(a1)
            cells = np.dot(offsets, cell)
            if len(neighbors) > 0:
                d = positions[neighbors] + cells - positions[a1]
                r2 = (d**2).sum(1)
                r = np.sqrt(r2)
                evec = mvec(r)
                tmp = evec.sum()
                local_energy[a1] += tmp/2.0
                energy += tmp
                f = (mdvec(r) / r)[:, np.newaxis] * d
                forces[a1] += f.sum(axis=0)
                for a2, f2, le in zip(neighbors, f, evec):
                    forces[a2] -= f2
                    local_energy[a2] += le/2.0
                stress -= np.dot(f.T, d)
                            
        if 'stress' in properties:
            if self.atoms.number_of_lattice_vectors == 3:
                stress += stress.T.copy()
                stress *= -0.5 / self.atoms.get_volume()
                self.results['stress'] = stress.flat[[0, 4, 8, 5, 2, 1]]
            else:
                raise PropertyNotImplementedError
                                
        self.results['energy'] = energy
        self.results['free_energy'] = energy
        self.results['forces'] = forces
        self.results['local_energy'] = local_energy
                                                                                                                                                                                                                                                                                                                                                                            
