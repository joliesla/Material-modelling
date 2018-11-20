import numpy as np
import json
from ase import Atoms
from ase.data import atomic_numbers
from ase.units import m
from ase.calculators.singlepoint import SinglePointCalculator

sysnames = {'configuration_periodic_dimensions',
            'simulation_cell',
            'atom_labels',
            'atom_positions'}


def toarray(ev):
    if ev is None:
        return None

    return np.array(ev['flatValues']).reshape(ev['valuesShape'])

def read_nomad_parse_events(fd, index):  # grrr index!  WTF
    data = json.load(fd)

    images = []
    systems = {}

    for ev in data['events']:
        event = ev['event']
        if event == 'openSection':
            if ev['metaName'] == 'section_system':
                sys_gindex = ev['gIndex']
                assert sys_gindex not in systems
                systems[sys_gindex] = {}

            if ev['metaName'] == '':
                pass

        elif event == 'addArrayValues':
            if ev['metaName'] in sysnames:
                systems[sys_gindex][ev['metaName']] = ev

        elif event == 'closeSection':
            if ev['metaName'] == 'section_system':
                sys_gindex = None

    for gindex in systems:
        system = systems[gindex]
        labels = toarray(system['atom_labels'])
        positions = toarray(system['atom_positions']) * m
        cell = toarray(system.get('simulation_cell'))
        pbc = toarray(system.get('configuration_periodic_dimensions'))

        numbers = [atomic_numbers.get(label, 0) for label in labels]
        atoms = Atoms(numbers, positions=positions, pbc=pbc)
        if cell is not None:
            atoms.cell = cell * m
        images.append(atoms)
    return images
